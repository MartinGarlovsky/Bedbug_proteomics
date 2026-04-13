#!/usr/bin/env python3
"""
Region-Constrained Motif-Based Active Site Analysis

1. Uses alignment to identify the conserved catalytic domain region
2. Searches for motif only within that region
3. Scores based on strict chemistry rules

This prevents false positives from finding motif patterns in wrong parts of the sequence.
"""

import os
import sys
import json
import csv
from pathlib import Path
import subprocess
import tempfile

try:
    from Bio import SeqIO, AlignIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    print("ERROR: BioPython required")
    sys.exit(1)


class RegionConstrainedMotifFinder:
    """
    Find active sites with motif search constrained to catalytic domain
    """
    
    # Known reference positions to define search region
    REFERENCE_INFO = {
        'Q9V3D8': {'active_site_center': 360, 'search_window': 150},
        'P00727': {'active_site_center': 310, 'search_window': 150}
    }
    
    # Motif pattern — residue sets are grounded in M17 LAP crystal structures and
    # site-directed mutagenesis data, NOT inferred from general chemistry.
    #
    # Primary structural references:
    #   blLAP  = bovine lens LAP (Burley et al. 1992 J Mol Biol 224:113; PDB 1LAM)
    #   PfLAP  = P. falciparum M17 LAP (McGowan et al. 2010 PNAS 107:2449; PDB 3KBK)
    #   tLAP-A = tomato LAP-A mutagenesis (Gu & Walling 2002 Eur J Biochem 269:1630)
    #   Dorus  = Dorus et al. 2011 BMC Genomics 12:177 (S-LAP family classification)
    #
    # blLAP active site (Burley et al. 1992):
    #   Site 2 (tight Zn):  Asp255, Asp332, Glu334 coordinate; Lys262 and Arg336 catalytic
    #   Site 1 (loose Zn):  Asp273, Glu334 (bridging), Lys250 coordinate
    #
    # Numbering below uses the Dorus et al. S-LAP / M17 consensus positions.
    # Each entry has:
    #   offset     — residues from metal_1 anchor (Lys327 equivalent)
    #   tolerance  — positional search window ±n residues
    #   canonical  — exact M17 LAP consensus residue(s)           → CONSERVED
    #   compatible — residues with structural/biochemical support
    #                for metal coordination in M17 LAPs            → COMPATIBLE
    #   (anything else)                                            → NON_FUNCTIONAL
    #
    # Catalytic sites have NO compatible class: mutagenesis shows all tested
    # substitutions at Lys339 reduce activity (Gu & Walling 2002), and Arg413
    # is essential for transition-state stabilisation (Sträter & Lipscomb 1995).
    MOTIF_PATTERN = {

        # Lys327 (blLAP Lys250): coordinates site 1 Zn via side-chain amino group.
        # Gln: observed in functionally classified active M17 LAPs (Dorus 2011,
        #   Cluster II S-LAPs). The amide nitrogen of Gln can substitute for the
        #   amino group of Lys in metal coordination.
        # Arg: basic, same charge as Lys, but not directly observed coordinating
        #   metal in any M17 LAP structure — excluded.
        'metal_1': {
            'offset': 0,  'tolerance': 1,
            'canonical':   ['K'],
            'compatible':  ['Q'],
        },

        # Asp332 (blLAP Asp255): coordinates site 2 tight Zn via carboxylate oxygen.
        # Asn: amide oxygen can coordinate metals; conservative substitution of Asp;
        #   observed in S-LAP Cluster I (Dorus 2011).
        # Glu: acidic, conservative; carboxylate oxygen retains coordinating ability.
        # His: coordinates metals via imidazole nitrogen; observed in loopin-1
        #   (Dorus 2011); well-established metal ligand in metalloenzymes generally.
        # Cys: thiol coordinates metals; observed in S-LAP5 (Dorus 2011).
        'metal_2': {
            'offset': 5,  'tolerance': 1,
            'canonical':   ['D'],
            'compatible':  ['N', 'E', 'H', 'C'],
        },

        # Asp350 (blLAP Asp273): coordinates site 1 loose Zn via carboxylate oxygen.
        # Asn: amide oxygen can coordinate; conservative substitution (Dorus 2011).
        # Glu: acidic conservative substitution; retains carboxylate.
        # Ser/Thr: hydroxyl oxygen can coordinate divalent metals; observed in
        #   S-LAP Cluster I variants (Dorus 2011); partial support from tLAP-A
        #   E429S retaining trace activity (Gu & Walling 2002, different position
        #   but indicates hydroxyl has some coordinating capacity).
        # Cys: thiol coordinates metals (Dorus 2011, S-LAP5).
        # Arg: positively charged, no carboxylate/hydroxyl/thiol — cannot coordinate
        #   metal at a site requiring an oxygen/sulphur donor. Removed.
        'metal_3': {
            'offset': 23, 'tolerance': 2,
            'canonical':   ['D'],
            'compatible':  ['N', 'E', 'S', 'T', 'C'],
        },

        # Lys409 / Glu411 — NOTE on numbering:
        # Dorus et al. use S-LAP positions 409 and 411, but their functional roles
        # map to blLAP as follows:
        #   offset 82 (Dorus "Lys409") = blLAP Glu334 — bridges BOTH Zn sites
        #   offset 84 (Dorus "Glu411") = blLAP Lys409 — site 1 loose coordination
        # The mutagenesis data (Gu & Walling 2002) applies to the blLAP Glu334
        # equivalent, which is at offset 82.

        # offset 82 = blLAP Lys409: site 1 loose Zn coordination.
        # Less conserved and less mutagenesis-constrained than Glu334.
        # Lys→Asp and Lys→Ser substitutions observed in active M17 LAPs (Dorus 2011).
        'metal_4': {
            'offset': 82, 'tolerance': 2,
            'canonical':   ['K'],           # blLAP Lys409 equivalent
            'compatible':  ['D', 'S'],      # observed in active M17 LAPs (Dorus 2011)
        },

        # offset 84 = blLAP Glu334: bridges both Zn sites via two carboxylate oxygens.
        # The most mutagenesis-constrained metal-binding position.
        # tLAP-A Glu429 (Gu & Walling 2002):
        #   Asp (E→D): retained >95% activity — strongly validated compatible.
        #   Val/Ser: <3% activity — not compatible.
        'metal_5': {
            'offset': 84, 'tolerance': 2,
            'canonical':   ['E'],           # blLAP Glu334 equivalent
            'compatible':  ['D'],           # E→D: >95% activity retained (Gu & Walling 2002)
        },

        # Lys339 (blLAP Lys262): catalytic base; activates the nucleophilic water.
        # Mutagenesis (tLAP-A Lys354, Gu & Walling 2002): most substitutions reduce
        # activity substantially. No substitution shown to retain full activity.
        # NO compatible class — any divergence scored NON_FUNCTIONAL.
        'cat_1': {
            'offset': 12, 'tolerance': 1,
            'canonical':   ['K'],
            'compatible':  [],
        },

        # Arg413 (blLAP Arg336): stabilises the gem-diolate transition state
        # (Sträter & Lipscomb 1995). Essential for catalysis; no mutagenesis data
        # supports any substitution retaining activity.
        # NO compatible class — any divergence scored NON_FUNCTIONAL.
        'cat_2': {
            'offset': 86, 'tolerance': 2,
            'canonical':   ['R'],
            'compatible':  [],
        },
    }
    
    def __init__(self, fasta_file):
        """Initialize"""
        self.fasta_file = fasta_file
        self.sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}
        self.alignment = None
        self.catalytic_regions = {}
        
        print(f"Loaded {len(self.sequences)} sequences")
    
    def select_reference(self):
        """Select best available reference"""
        for ref_id in ['Q9V3D8', 'P00727']:
            if ref_id in self.sequences:
                print(f"Using reference: {ref_id}")
                return ref_id
        
        print("ERROR: No reference sequence found!")
        return None
    
    def perform_alignment(self):
        """Align all sequences"""
        print("\nPerforming multiple sequence alignment...")
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_in:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_out:
                tmp_in_path = tmp_in.name
                tmp_out_path = tmp_out.name
        
        try:
            SeqIO.write([SeqIO.SeqRecord(seq=seq, id=sid, description='') 
                        for sid, seq in self.sequences.items()], 
                       tmp_in_path, "fasta")
            
            try:
                subprocess.run(['muscle', '-in', tmp_in_path, '-out', tmp_out_path],
                             check=True, capture_output=True, text=True)
            except FileNotFoundError:
                print("ERROR: MUSCLE not found. Install with: conda install -c bioconda muscle")
                return False
            
            self.alignment = AlignIO.read(tmp_out_path, "fasta")
            print(f"✓ Alignment complete")
            return True
            
        finally:
            if os.path.exists(tmp_in_path):
                os.remove(tmp_in_path)
            if os.path.exists(tmp_out_path):
                os.remove(tmp_out_path)
    
    def identify_catalytic_regions(self, ref_id):
        """
        Map catalytic domain region from reference to all sequences, and record
        the exact alignment-derived position of the metal_1 anchor (ref Lys327)
        in each sequence.

        Using the alignment to pin the anchor prevents the motif scanner from
        latching onto a different Lys elsewhere in the catalytic window.
        """
        print(f"\nMapping catalytic domain regions from {ref_id}...")

        ref_info = self.REFERENCE_INFO[ref_id]
        center = ref_info['active_site_center']
        window = ref_info['search_window']

        # Find reference record in alignment
        ref_record = None
        for record in self.alignment:
            if record.id == ref_id:
                ref_record = record
                break
        if not ref_record:
            return

        # Build reference seq-pos → alignment-col map
        ref_seq_to_aln = {}
        seq_pos = 0
        for aln_pos, residue in enumerate(str(ref_record.seq)):
            if residue != '-':
                seq_pos += 1
                ref_seq_to_aln[seq_pos] = aln_pos

        # Catalytic window boundaries in alignment columns
        start_seq = max(1, center - window)
        end_seq   = min(seq_pos, center + window)
        start_aln = ref_seq_to_aln.get(start_seq, 0)
        end_aln   = ref_seq_to_aln.get(end_seq, len(ref_record.seq))

        # Alignment column of metal_1 anchor in the reference
        # metal_1 offset=0, so anchor = the Lys at 'center' adjusted to the
        # actual motif start.  We find it by searching for the highest-scoring
        # motif in the REFERENCE first to pin the exact anchor column.
        # Simpler: use the reference's own catalytic region start as a proxy,
        # then find the real anchor col by scoring the reference motif.
        ref_anchor_aln_col = None  # will be set below after mapping

        # Map window to each sequence; also build aln-col → seq-pos map per seq
        for record in self.alignment:
            seq_id  = record.id
            aln_seq = str(record.seq)

            aln_to_seq = {}
            s = 0
            for aln_pos, residue in enumerate(aln_seq):
                if residue != '-':
                    s += 1
                    aln_to_seq[aln_pos] = s

            seq_positions = [aln_to_seq[c] for c in range(start_aln, end_aln + 1)
                             if c in aln_to_seq]
            if not seq_positions:
                continue

            self.catalytic_regions[seq_id] = {
                'start':       min(seq_positions),
                'end':         max(seq_positions),
                'aln_to_seq':  aln_to_seq,   # store for anchor lookup
            }
            print(f"  {seq_id}: positions {min(seq_positions)}-{max(seq_positions)}")

        # Now find the reference anchor column by scoring the reference's own motif
        if ref_id in self.catalytic_regions:
            ref_region = self.catalytic_regions[ref_id]
            ref_seq    = self.sequences[ref_id]
            best_score = -1
            best_anchor_seq_pos = None
            for pos in range(ref_region['start'], ref_region['end'] - 90):
                if ref_seq[pos - 1] not in ['K', 'R']:
                    continue
                motif = self.score_motif_at_position(ref_seq, pos)
                if motif and motif['overall_score'] > best_score:
                    best_score = motif['overall_score']
                    best_anchor_seq_pos = pos
            if best_anchor_seq_pos is not None:
                # Map that seq position back to an alignment column
                ref_aln_to_seq = {v: k for k, v in ref_region['aln_to_seq'].items()}
                ref_anchor_aln_col = ref_aln_to_seq.get(best_anchor_seq_pos)
                print(f"\n  Reference anchor: seq pos {best_anchor_seq_pos} "
                      f"→ alignment col {ref_anchor_aln_col}")

        # Store the alignment-mapped anchor seq-position for every sequence
        if ref_anchor_aln_col is not None:
            for seq_id, region in self.catalytic_regions.items():
                anchor_seq_pos = region['aln_to_seq'].get(ref_anchor_aln_col)
                region['anchor_pos'] = anchor_seq_pos
                if anchor_seq_pos:
                    print(f"  {seq_id}: anchor at seq pos {anchor_seq_pos} "
                          f"(residue: {self.sequences[seq_id][anchor_seq_pos - 1]})")
                else:
                    print(f"  {seq_id}: anchor position gapped in alignment")
    
    # Classification tiers (mirrors the published figure colour scheme)
    CONSERVED     = 'CONSERVED'     # red   — matches M17 LAP consensus exactly
    COMPATIBLE    = 'COMPATIBLE'    # pink  — known metal ligand in other M17 LAPs
    NON_FUNCTIONAL = 'NON_FUNCTIONAL'  # white — neither

    # Numeric scores per tier used for overall scoring
    TIER_SCORE = {CONSERVED: 1.0, COMPATIBLE: 0.5, NON_FUNCTIONAL: 0.0}

    def _classify_residue(self, residue, site_name):
        """Return (tier, score) for a residue at a given site."""
        config = self.MOTIF_PATTERN[site_name]
        if residue in config['canonical']:
            return self.CONSERVED, self.TIER_SCORE[self.CONSERVED]
        if residue in config['compatible']:
            return self.COMPATIBLE, self.TIER_SCORE[self.COMPATIBLE]
        return self.NON_FUNCTIONAL, self.TIER_SCORE[self.NON_FUNCTIONAL]

    def score_motif_at_position(self, sequence, start_pos):
        """Score a motif candidate using the three-tier literature-curated scheme.

        Classification per site:
          CONSERVED      — residue matches the M17 LAP consensus (score 1.0)
          COMPATIBLE     — known metal ligand in other M17 LAPs   (score 0.5)
          NON_FUNCTIONAL — neither                                 (score 0.0)

        Catalytic sites (cat_1, cat_2) have no COMPATIBLE class.
        A small positional penalty (0.05 per offset step) is applied when the
        best match is not at the exact expected position.

        Unique position constraint: sites are evaluated in ascending offset order.
        Once a sequence position is claimed by one site, it cannot be claimed by
        another. This prevents metal_4 and metal_5 (offset 82 and 84) from
        colliding onto the same residue when the tolerance window is wide relative
        to their spacing.
        """
        sites_found = {}
        claimed_positions = set()

        # Process sites in ascending offset order so earlier offsets have priority
        sorted_sites = sorted(self.MOTIF_PATTERN.items(), key=lambda x: x[1]['offset'])

        for site_name, config in sorted_sites:
            expected_pos = start_pos + config['offset']
            tolerance    = config['tolerance']

            best_match = None
            best_score = -1

            for offset in range(-tolerance, tolerance + 1):
                check_pos = expected_pos + offset
                if not (1 <= check_pos <= len(sequence)):
                    continue
                if check_pos in claimed_positions:
                    continue  # position already claimed by an earlier site

                residue = sequence[check_pos - 1]
                tier, score = self._classify_residue(residue, site_name)

                position_penalty = abs(offset) * 0.05
                adjusted = score - position_penalty

                if adjusted > best_score:
                    best_score = adjusted
                    best_match = {
                        'position':             check_pos,
                        'residue':              residue,
                        'offset_from_expected': offset,
                        'tier':                 tier,
                        'score':                score,
                    }

            if best_match:
                sites_found[site_name] = best_match
                claimed_positions.add(best_match['position'])

        if not sites_found:
            return None

        metal_sites = {k: v for k, v in sites_found.items() if k.startswith('metal')}
        cat_sites   = {k: v for k, v in sites_found.items() if k.startswith('cat')}

        metal_score = sum(s['score'] for s in metal_sites.values()) / 5 if metal_sites else 0
        cat_score   = sum(s['score'] for s in cat_sites.values())   / 2 if cat_sites   else 0

        return {
            'start_pos':      start_pos,
            'sites':          sites_found,
            'n_sites_found':  len(sites_found),
            'metal_score':    metal_score,
            'catalytic_score': cat_score,
            'overall_score':  (metal_score * 0.6) + (cat_score * 0.4),
        }
    
    def find_best_motif_in_region(self, seq_id, sequence):
        """
        Score the motif using the alignment-pinned anchor position.

        The anchor (metal_1, Lys327 equivalent) is identified by mapping the
        reference anchor column through the multiple sequence alignment, so we
        never accidentally latch onto a different Lys elsewhere in the window.

        Falls back to scanning if no alignment anchor is available.
        """
        if seq_id not in self.catalytic_regions:
            return None

        region = self.catalytic_regions[seq_id]
        anchor_pos = region.get('anchor_pos')

        if anchor_pos is not None:
            # Use the alignment-derived anchor directly
            motif = self.score_motif_at_position(sequence, anchor_pos)
            return motif
        else:
            # Fallback: scan the region (used if anchor column was gapped)
            print(f"  Warning: {seq_id} has no alignment anchor — falling back to scan")
            best_motif  = None
            best_score  = -1
            for pos in range(region['start'], region['end'] - 86):
                if sequence[pos - 1] not in ['K', 'R', 'Q']:  # Q for Gln-type anchor
                    continue
                motif = self.score_motif_at_position(sequence, pos)
                if motif and motif['overall_score'] > best_score:
                    best_score = motif['overall_score']
                    best_motif = motif
            return best_motif
    
    def assess_activity(self, motif):
        """Classify activity using tier counts, mirroring the published figure logic.

        LIKELY_ACTIVE:   all 5 metal sites CONSERVED or COMPATIBLE,
                         both catalytic sites CONSERVED
        POSSIBLY_ACTIVE: ≥3 metal sites CONSERVED/COMPATIBLE,
                         ≥1 catalytic site CONSERVED
        UNLIKELY_ACTIVE: does not meet the above
        NOT_FOUND:       no motif detected in the search region
        """
        if not motif:
            return 'NOT_FOUND'

        sites = motif['sites']

        metal_sites = {k: v for k, v in sites.items() if k.startswith('metal')}
        cat_sites   = {k: v for k, v in sites.items() if k.startswith('cat')}

        metal_functional = sum(
            1 for s in metal_sites.values()
            if s['tier'] in (self.CONSERVED, self.COMPATIBLE)
        )
        cat_conserved = sum(
            1 for s in cat_sites.values()
            if s['tier'] == self.CONSERVED
        )

        if metal_functional >= 5 and cat_conserved == 2:
            return 'LIKELY_ACTIVE'
        elif metal_functional >= 3 and cat_conserved >= 1:
            return 'POSSIBLY_ACTIVE'
        else:
            return 'UNLIKELY_ACTIVE'
    
    def analyze_all(self, ref_id):
        """Analyze all sequences"""
        results = []
        
        print("\n" + "="*80)
        print("ANALYZING SEQUENCES")
        print("="*80 + "\n")
        
        for seq_id, sequence in self.sequences.items():
            print(f"Analyzing {seq_id}...", end=" ")
            
            motif = self.find_best_motif_in_region(seq_id, sequence)
            assessment = self.assess_activity(motif)
            
            result = {
                'protein_id': seq_id,
                'reference_id': ref_id,
                'motif_found': motif is not None,
                'assessment': assessment
            }
            
            if motif:
                result.update({
                    'motif_start': motif['start_pos'],
                    'n_sites_found': motif['n_sites_found'],
                    'overall_score': motif['overall_score'],
                    'metal_score': motif['metal_score'],
                    'catalytic_score': motif['catalytic_score'],
                    'sites': motif['sites'],
                    'search_region': self.catalytic_regions.get(seq_id, {})
                })
                print(f"{assessment} (score: {motif['overall_score']:.2f})")
            else:
                result.update({
                    'overall_score': 0,
                    'metal_score': 0,
                    'catalytic_score': 0,
                    'n_sites_found': 0
                })
                print("NO MOTIF")
            
            results.append(result)
        
        return results


def write_alignment_outputs(finder):
    """Write alignment files for manual inspection in three formats.

    1. catalytic_domain_alignment.fasta — full MSA, aligned FASTA
       Open in: Jalview (jalview.org) — recommended; AliView; SeaView
    2. catalytic_domain_alignment.aln   — full MSA, ClustalW format
       Open in: BioEdit; Clustal Omega viewer
    3. catalytic_domain_region.fasta    — catalytic window only, aligned FASTA
       Open in same viewers, or paste into ESPript 3 (espript.ibcp.fr) for figures

    Recommended workflow:
      Open catalytic_domain_alignment.fasta in Jalview
      → Colour → Zappo (physicochemical groups, charge inversions stand out)
      → add annotation row marking the 7 motif positions
    """
    if finder.alignment is None:
        print("  ! No alignment available — skipping alignment outputs")
        return

    AlignIO.write(finder.alignment, 'catalytic_domain_alignment.fasta', 'fasta')
    print("  ✓ catalytic_domain_alignment.fasta  (open in Jalview / AliView / SeaView)")

    AlignIO.write(finder.alignment, 'catalytic_domain_alignment.aln', 'clustal')
    print("  ✓ catalytic_domain_alignment.aln    (open in BioEdit / Clustal viewer)")

    # Trimmed catalytic-window-only alignment
    ref_id = None
    ref_record = None
    for sid in ['Q9V3D8', 'A1Z9G3', 'P00727']:
        for record in finder.alignment:
            if record.id == sid:
                ref_id = sid
                ref_record = record
                break
        if ref_id:
            break

    if ref_id and ref_id in finder.catalytic_regions:
        region = finder.catalytic_regions[ref_id]
        seq_pos = 0
        aln_start = aln_end = None
        for aln_pos, residue in enumerate(str(ref_record.seq)):
            if residue != '-':
                seq_pos += 1
                if seq_pos == region['start'] and aln_start is None:
                    aln_start = aln_pos
                if seq_pos == region['end']:
                    aln_end = aln_pos
        if aln_start is not None and aln_end is not None:
            trimmed = finder.alignment[:, aln_start:aln_end + 1]
            AlignIO.write(trimmed, 'catalytic_domain_region.fasta', 'fasta')
            print("  ✓ catalytic_domain_region.fasta    (catalytic window only — paste into ESPript 3)")
        else:
            print("  ! Could not determine catalytic window columns — skipping trimmed alignment")
    else:
        print("  ! Reference not found in alignment — skipping trimmed alignment")


def write_outputs(results, finder=None):
    """Write outputs"""
    print("\n" + "="*80)
    print("Writing outputs...")
    print("="*80 + "\n")
    
    with open('region_constrained_results.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'protein_id', 'assessment', 'overall_score',
            'metal_score', 'catalytic_score', 'n_sites_found',
            'motif_start', 'search_region_start', 'search_region_end', 'anchor_pos'
        ])
        
        for r in results:
            region = r.get('search_region', {})
            writer.writerow([
                r['protein_id'],
                r['assessment'],
                f"{r['overall_score']:.4f}",
                f"{r['metal_score']:.4f}",
                f"{r['catalytic_score']:.4f}",
                r['n_sites_found'],
                r.get('motif_start', 'NA'),
                region.get('start', 'NA'),
                region.get('end', 'NA'),
                region.get('anchor_pos', 'NA'),
            ])
    
    print("  ✓ region_constrained_results.tsv")
    
    with open('region_constrained_details.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'protein_id', 'site_name', 'position', 'residue',
            'canonical_residues', 'compatible_residues', 'offset', 'tier', 'score'
        ])
        
        for r in results:
            if 'sites' in r:
                for site_name, site_info in r['sites'].items():
                    canonical   = RegionConstrainedMotifFinder.MOTIF_PATTERN[site_name]['canonical']
                    compatible  = RegionConstrainedMotifFinder.MOTIF_PATTERN[site_name]['compatible']
                    writer.writerow([
                        r['protein_id'],
                        site_name,
                        site_info['position'],
                        site_info['residue'],
                        '/'.join(canonical),
                        '/'.join(compatible) if compatible else 'none',
                        site_info['offset_from_expected'],
                        site_info.get('tier', 'NA'),
                        f"{site_info['score']:.4f}",
                    ])
    
    print("  ✓ region_constrained_details.tsv")
    
    # Strip internal alignment data before serialising
    for r in results:
        if 'search_region' in r:
            r['search_region'] = {k: v for k, v in r['search_region'].items()
                                  if k != 'aln_to_seq'}

    with open('region_constrained_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print("  ✓ region_constrained_results.json")

    if finder is not None:
        write_alignment_outputs(finder)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nUsage: python region_constrained_analysis.py sequences.fasta")
        return
    
    finder = RegionConstrainedMotifFinder(sys.argv[1])
    
    ref_id = finder.select_reference()
    if not ref_id:
        return
    
    if not finder.perform_alignment():
        return
    
    finder.identify_catalytic_regions(ref_id)
    
    results = finder.analyze_all(ref_id)
    
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80 + "\n")
    
    likely = [r for r in results if r['assessment'] == 'LIKELY_ACTIVE']
    possibly = [r for r in results if r['assessment'] == 'POSSIBLY_ACTIVE']
    unlikely = [r for r in results if r['assessment'] == 'UNLIKELY_ACTIVE']
    
    print(f"LIKELY ACTIVE:    {len(likely)}")
    print(f"POSSIBLY ACTIVE:  {len(possibly)}")
    print(f"UNLIKELY ACTIVE:  {len(unlikely)}")
    
    write_outputs(results, finder=finder)
    
    print("\n" + "="*80)
    print("DONE!")
    print("="*80)


if __name__ == "__main__":
    main()
