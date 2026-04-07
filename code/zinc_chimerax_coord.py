from chimerax.core.commands import run
from chimerax.atomic import all_atoms
import numpy as np
import os
import csv

# ── USER SETTINGS ────────────────────────────────────────────────────────────
CIF_DIR    = "/Users/martingarlovsky/Dropbox/TUDAZ/Cimex/Bedbug_proteomics/SLAP_activity/af3_comparisons/zn_cifs"
OUTPUT_CSV = "/Users/martingarlovsky/Dropbox/TUDAZ/Cimex/Bedbug_proteomics/SLAP_activity/af3_comparisons/zinc_summary_chimerax.csv"
CUTOFF     = 3.5    # Angstrom cutoff for initial contact search
# ─────────────────────────────────────────────────────────────────────────────

# ── M17 blLAP COORDINATION CHEMISTRY (from PDB 1LAP, Burley et al. 1992) ─────
#
# Two zinc sites with DIFFERENT coordinating residues:
#
#   Zn1 (exchangeable / "loose" site):
#       Lys(NZ) + Asp(OD) + Glu(OE) — 3 sidechain contacts
#       KEY FEATURE: Lys is ONLY at this site
#
#   Zn2 (tight site):
#       Asp(OD) x2 + Glu(OE) + Asp(O backbone) — 3 sidechain + 1 backbone
#       NO Lys at this site
#
# AF3 note: predicted Zn–ligand distances are typically 0.3–0.6 Å longer
# than crystal structures because AF3 does not energy-minimise metal geometry.
# Thresholds are relaxed accordingly.

M17_SIDECHAIN_ATOMS = {
    "LYS": {"NZ"},
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}

# AF3-adjusted distance thresholds (Angstrom)
DIST_IDEAL_MAX  = 2.8   # ideal for AF3 predicted structures
DIST_ACCEPT_MAX = 3.2   # acceptable upper limit for AF3
PLDDT_MIN       = 70.0

# Score thresholds for final call
SCORE_ACTIVE    = 6
SCORE_UNCERTAIN = 3
# ─────────────────────────────────────────────────────────────────────────────


def is_canonical_sidechain(atom):
    res  = atom.residue.name.upper()
    name = atom.name.upper()
    return res in M17_SIDECHAIN_ATOMS and name in M17_SIDECHAIN_ATOMS[res]


def deduplicate_bidentate(coordinating):
    """
    For bidentate residues (Asp/Glu with both OD1/OD2 or OE1/OE2 in range),
    keep only the closer oxygen. This prevents double-counting a single
    coordinating residue as two contacts and inflating coordination number.
    Returns deduplicated list and a note about any bidentate residues found.
    """
    # Group canonical sidechain contacts by (chain, residue name, res number)
    groups = {}
    non_grouped = []
    for c in coordinating:
        if c["canonical_sidechain"] and c["residue"] in ("ASP", "GLU"):
            key = (c["chain"], c["residue"], c["res_num"])
            groups.setdefault(key, []).append(c)
        else:
            non_grouped.append(c)

    deduped = list(non_grouped)
    bidentate_notes = []
    for key, contacts in groups.items():
        if len(contacts) > 1:
            # Keep only the closest oxygen
            closest = min(contacts, key=lambda x: x["distance"])
            deduped.append(closest)
            bidentate_notes.append(
                f"{key[1]}{key[2]} bidentate ({[round(c['distance'],3) for c in contacts]}A, kept {closest['distance']}A)"
            )
        else:
            deduped.append(contacts[0])

    return deduped, bidentate_notes


def classify_site(coordinating_deduped):
    """
    Classify which zinc site this is based on coordinating residues.
    Zn1 (exchangeable): has Lys
    Zn2 (tight):        no Lys, has multiple Asp
    Returns 'Zn1_exchangeable', 'Zn2_tight', or 'unknown'
    """
    canon = [c for c in coordinating_deduped if c["canonical_sidechain"]]
    res_types = {c["residue"] for c in canon}
    if "LYS" in res_types:
        return "Zn1_exchangeable"
    elif "ASP" in res_types:
        return "Zn2_tight"
    return "unknown"


def score_zinc_site(coordinating_raw, mean_plddt):
    """
    Score a zinc site against M17 LAP expectations.
    Applies per-site rules based on whether Lys is present.

    Scoring rubric (max 9):
      Zn1 (Lys present):
        +3  LYS coordinates zinc               (M17 Zn1 hallmark)
        +2  ASP present
        +1  GLU present
        +2  deduplicated coord number = 3      (expected for Zn1)
        +1  all contacts within DIST_IDEAL_MAX
        +1  pLDDT >= 70

      Zn2 (no Lys):
        +3  2+ ASP contacts                    (Zn2 has two Asp)
        +2  GLU present                        (bridging Glu334)
        +1  backbone Asp contact present       (Asp332 carbonyl in crystal)
        +2  deduplicated coord number = 3-4
        +1  all contacts within DIST_IDEAL_MAX
        +1  pLDDT >= 70 (replaces backbone check if backbone contact absent)

      Penalties:
        -1  any contact > DIST_ACCEPT_MAX
    """
    score = 0
    flags = []

    deduped, bidentate_notes = deduplicate_bidentate(coordinating_raw)
    if bidentate_notes:
        flags.append(f"bidentate_deduplicated: {'; '.join(bidentate_notes)}")

    canon   = [c for c in deduped if c["canonical_sidechain"]]
    nonsc   = [c for c in deduped if not c["canonical_sidechain"]]
    res_types = {c["residue"] for c in canon}
    site    = classify_site(deduped)
    n       = len(canon)

    if site == "Zn1_exchangeable":
        # Zn1 scoring
        if "LYS" in res_types:
            score += 3
        else:
            flags.append("LYS_MISSING")

        asp_contacts = [c for c in canon if c["residue"] == "ASP"]
        if asp_contacts:
            score += 2
        else:
            flags.append("ASP_MISSING")

        if "GLU" in res_types:
            score += 1
        else:
            flags.append("GLU_MISSING")

        if n == 3:
            score += 2
        elif n == 2:
            score += 1
            flags.append(f"low_coord_number_{n}")
        else:
            flags.append(f"unexpected_coord_number_{n}_for_Zn1")

    else:
        # Zn2 scoring (tight site — no Lys expected)
        asp_contacts = [c for c in canon if c["residue"] == "ASP"]
        if len(asp_contacts) >= 2:
            score += 3
        elif len(asp_contacts) == 1:
            score += 1
            flags.append("only_1_ASP_contact_Zn2_expects_2")
        else:
            flags.append("ASP_MISSING")

        if "GLU" in res_types:
            score += 2
        else:
            flags.append("GLU_MISSING_Zn2")

        # Backbone Asp carbonyl contact (expected for Zn2)
        backbone_asp = [c for c in nonsc if c["residue"] == "ASP" and c["atom_name"] == "O"]
        if backbone_asp:
            score += 1
            flags.append("backbone_Asp_O_present (expected for Zn2)")
        else:
            # Give pLDDT point here instead if pLDDT is good
            pass

        if 3 <= n <= 4:
            score += 2
        elif n == 2:
            score += 1
            flags.append(f"low_coord_number_{n}")
        else:
            flags.append(f"unexpected_coord_number_{n}_for_Zn2")

    # Distance quality (applies to both sites)
    bad_contacts = [c for c in canon if c["distance"] > DIST_ACCEPT_MAX]
    good_contacts = [c for c in canon if c["distance"] <= DIST_IDEAL_MAX]
    if not bad_contacts and len(good_contacts) == len(canon):
        score += 1
    if bad_contacts:
        flags.append(f"contacts_>{DIST_ACCEPT_MAX}A: {[c['distance'] for c in bad_contacts]}")
        score -= 1

    # pLDDT
    if mean_plddt is not None:
        if mean_plddt >= PLDDT_MIN:
            score += 1
        else:
            flags.append(f"low_plddt_{round(mean_plddt, 1)}")

    if score >= SCORE_ACTIVE:
        call = "LIKELY_ACTIVE"
    elif score >= SCORE_UNCERTAIN:
        call = "UNCERTAIN"
    else:
        call = "LIKELY_INACTIVE"

    return score, site, call, flags, deduped


# ── MAIN LOOP ─────────────────────────────────────────────────────────────────
results = []
cif_files = sorted([f for f in os.listdir(CIF_DIR) if f.endswith(".cif")])
print(f"Found {len(cif_files)} CIF files to process.\n")

for cif_file in cif_files:
    protein_name = cif_file.replace(".cif", "")
    filepath     = os.path.join(CIF_DIR, cif_file)
    print(f"Processing: {protein_name}")

    run(session, f"open {filepath}")
    atoms = all_atoms(session)

    zn_atoms     = [a for a in atoms if a.element.name == "Zn"]
    num_zn_total = len(zn_atoms)

    if not zn_atoms:
        print(f"  No zinc found")
        results.append({
            "protein": protein_name, "zinc_found": False,
            "zn_index": None, "num_zn_total": 0, "site_type": None,
            "num_canonical_contacts": None, "coordinating_residues_raw": None,
            "coordinating_residues_deduped": None,
            "canonical_residue_types": None, "lys_present": False,
            "mean_distance": None, "mean_plddt_at_site": None,
            "score": 0, "call": "LIKELY_INACTIVE", "flags": "no_zinc_found",
        })
        run(session, "close all")
        continue

    for i, zn in enumerate(zn_atoms):
        zn_pos = np.array(zn.scene_coord)
        coordinating_raw = []

        for atom in atoms:
            if atom.element.name == "Zn":
                continue
            if atom.element.name not in ("N", "O", "S"):
                continue
            atom_pos = np.array(atom.scene_coord)
            dist = float(np.linalg.norm(zn_pos - atom_pos))
            if dist <= CUTOFF:
                coordinating_raw.append({
                    "residue":             atom.residue.name.upper(),
                    "res_num":             atom.residue.number,
                    "chain":               atom.residue.chain_id,
                    "atom_name":           atom.name.upper(),
                    "element":             atom.element.name,
                    "distance":            round(dist, 3),
                    "canonical_sidechain": is_canonical_sidechain(atom),
                    "bfactor":             atom.bfactor,
                })

        plddt_vals = [c["bfactor"] for c in coordinating_raw if c["bfactor"] is not None]
        mean_plddt = round(sum(plddt_vals) / len(plddt_vals), 1) if plddt_vals else None

        score, site, call, flags, deduped = score_zinc_site(coordinating_raw, mean_plddt)

        canon_deduped   = [c for c in deduped if c["canonical_sidechain"]]
        canon_res_types = sorted({c["residue"] for c in canon_deduped})
        lys_present     = "LYS" in canon_res_types

        def fmt_contacts(contacts):
            return "; ".join(
                f"{c['chain']}:{c['residue']}{c['res_num']}({c['atom_name']})={c['distance']}A"
                + ("*" if c["canonical_sidechain"] else "")
                for c in contacts
            )

        raw_summary    = fmt_contacts(coordinating_raw)
        deduped_summary = fmt_contacts(deduped)

        distances = [c["distance"] for c in deduped]
        mean_dist = round(sum(distances) / len(distances), 3) if distances else None
        flag_str  = " | ".join(flags) if flags else "none"

        print(f"  Zn[{i+1}] ({site}): score={score}/9 ({call}) | lys={lys_present} | plddt={mean_plddt}")
        if flags:
            print(f"    flags: {flag_str}")

        results.append({
            "protein":                      protein_name,
            "zinc_found":                   True,
            "zn_index":                     i + 1,
            "num_zn_total":                 num_zn_total,
            "site_type":                    site,
            "num_canonical_contacts":       len(canon_deduped),
            "coordinating_residues_raw":    raw_summary,
            "coordinating_residues_deduped": deduped_summary,
            "canonical_residue_types":      ", ".join(canon_res_types),
            "lys_present":                  lys_present,
            "mean_distance":                mean_dist,
            "mean_plddt_at_site":           mean_plddt,
            "score":                        score,
            "call":                         call,
            "flags":                        flag_str,
        })

    run(session, "close all")

# ── WRITE CSV ─────────────────────────────────────────────────────────────────
fieldnames = [
    "protein", "zinc_found", "zn_index", "num_zn_total", "site_type",
    "num_canonical_contacts", "coordinating_residues_raw",
    "coordinating_residues_deduped", "canonical_residue_types", "lys_present",
    "mean_distance", "mean_plddt_at_site", "score", "call", "flags",
]

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

# ── SUMMARY ───────────────────────────────────────────────────────────────────
calls = [r["call"] for r in results if r["zinc_found"]]
print(f"\n{'='*60}")
print(f"Done. {len(cif_files)} structures processed.")
print(f"  LIKELY_ACTIVE:   {calls.count('LIKELY_ACTIVE')}")
print(f"  UNCERTAIN:       {calls.count('UNCERTAIN')}")
print(f"  LIKELY_INACTIVE: {calls.count('LIKELY_INACTIVE')}")
print(f"Results written to: {OUTPUT_CSV}")
