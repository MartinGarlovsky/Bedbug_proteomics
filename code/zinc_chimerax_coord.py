from chimerax.core.commands import run
from chimerax.atomic import all_atoms
import numpy as np
import os
import csv

# ── USER SETTINGS ────────────────────────────────────────────────────────────
CIF_DIR    = "/Users/martingarlovsky/Dropbox/TUDAZ/Cimex/Bedbug_proteomics/SLAP_activity/af3_comparisons/zn_cifs"
OUTPUT_CSV = "/Users/martingarlovsky/Dropbox/TUDAZ/Cimex/Bedbug_proteomics/SLAP_activity/af3_comparisons/zinc_summary.csv"
CUTOFF     = 3.5
# ─────────────────────────────────────────────────────────────────────────────

# ── M17 blLAP COORDINATION CHEMISTRY (PDB 1LAP, Burley et al. 1992) ──────────
M17_SIDECHAIN_ATOMS = {
    "LYS": {"NZ"},
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}

DIST_IDEAL_MAX      = 2.8
DIST_ACCEPT_MAX     = 3.2
PLDDT_MIN           = 70.0
ZN_ZN_BINUCLEAR_MAX = 4.5
ZN_ZN_BORDERLINE    = 7.0
SCORE_CANONICAL  = 6
SCORE_UNCERTAIN     = 3
# ─────────────────────────────────────────────────────────────────────────────


def is_canonical_sidechain(atom):
    res  = atom.residue.name.upper()
    name = atom.name.upper()
    return res in M17_SIDECHAIN_ATOMS and name in M17_SIDECHAIN_ATOMS[res]


def deduplicate_bidentate(coordinating):
    groups      = {}
    non_grouped = []
    for c in coordinating:
        if c["canonical_sidechain"] and c["residue"] in ("ASP", "GLU"):
            key = (c["chain"], c["residue"], c["res_num"])
            groups.setdefault(key, []).append(c)
        else:
            non_grouped.append(c)
    deduped         = list(non_grouped)
    bidentate_notes = []
    for key, contacts in groups.items():
        if len(contacts) > 1:
            closest = min(contacts, key=lambda x: x["distance"])
            deduped.append(closest)
            bidentate_notes.append(
                f"{key[1]}{key[2]} bidentate "
                f"({[round(c['distance'],3) for c in contacts]}A, kept {closest['distance']}A)"
            )
        else:
            deduped.append(contacts[0])
    return deduped, bidentate_notes


def classify_site(coordinating_deduped):
    canon     = [c for c in coordinating_deduped if c["canonical_sidechain"]]
    res_types = {c["residue"] for c in canon}
    if "LYS" in res_types:
        return "Zn1_exchangeable"
    elif "ASP" in res_types:
        return "Zn2_tight"
    return "unknown"


def score_zinc_site(coordinating_raw, mean_plddt):
    score = 0
    flags = []

    deduped, bidentate_notes = deduplicate_bidentate(coordinating_raw)
    if bidentate_notes:
        flags.append(f"bidentate_deduplicated: {'; '.join(bidentate_notes)}")

    canon     = [c for c in deduped if c["canonical_sidechain"]]
    non_canon = [c for c in deduped if not c["canonical_sidechain"]]
    res_types = {c["residue"] for c in canon}
    site      = classify_site(deduped)
    n         = len(canon)

    if site == "Zn1_exchangeable":
        if "LYS" in res_types:
            score += 3
        else:
            flags.append("LYS_MISSING")
        if "ASP" in res_types:
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
        backbone_asp = [c for c in non_canon if c["residue"] == "ASP" and c["atom_name"] == "O"]
        if backbone_asp:
            score += 1
            flags.append("backbone_Asp_O_present (expected for Zn2)")
        if 3 <= n <= 4:
            score += 2
        elif n == 2:
            score += 1
            flags.append(f"low_coord_number_{n}")
        else:
            flags.append(f"unexpected_coord_number_{n}_for_Zn2")

    bad_contacts  = [c for c in canon if c["distance"] > DIST_ACCEPT_MAX]
    good_contacts = [c for c in canon if c["distance"] <= DIST_IDEAL_MAX]
    if not bad_contacts and len(good_contacts) == len(canon):
        score += 1
    if bad_contacts:
        flags.append(f"contacts_>{DIST_ACCEPT_MAX}A: {[c['distance'] for c in bad_contacts]}")
        score -= 1

    if non_canon:
        non_desc = ", ".join(
            f"{c['residue']}{c['res_num']}({c['atom_name']})" for c in non_canon
        )
        flags.append(f"non_canonical_contacts: {non_desc}")

    if mean_plddt is not None:
        if mean_plddt >= PLDDT_MIN:
            score += 1
        else:
            flags.append(f"low_plddt_{round(mean_plddt, 1)}")

    if score >= SCORE_CANONICAL:
        call = "CANONICAL"
    elif score >= SCORE_UNCERTAIN:
        call = "UNCERTAIN"
    else:
        call = "NON_CANONICAL"

    return score, site, call, flags, deduped


def assess_binuclear_site(zn_data):
    if len(zn_data) != 2:
        return {
            "zn_zn_distance":        None,
            "binuclear_distance_ok": None,
            "bridging_ligand":       None,
            "binuclear_call":        "single_zinc_only" if len(zn_data) == 1 else "no_zinc",
        }

    dist = round(float(np.linalg.norm(zn_data[0]["zn_pos"] - zn_data[1]["zn_pos"])), 3)

    if dist <= ZN_ZN_BINUCLEAR_MAX:
        dist_ok        = True
        binuclear_call = "BINUCLEAR_OK"
    elif dist <= ZN_ZN_BORDERLINE:
        dist_ok        = False
        binuclear_call = "BORDERLINE_SEPARATION"
    else:
        dist_ok        = False
        binuclear_call = "SEPARATE_POCKETS"

    def contact_keys(zn_entry):
        return {
            (c["chain"], c["residue"], c["res_num"])
            for c in zn_entry["deduped"] if c["canonical_sidechain"]
        }

    shared = contact_keys(zn_data[0]) & contact_keys(zn_data[1])

    if shared:
        bridging = "; ".join(f"{k[1]}{k[2]}(chain {k[0]})" for k in shared)
    else:
        bridging = "none"
        if binuclear_call == "BINUCLEAR_OK":
            binuclear_call = "BINUCLEAR_OK_no_shared_ligand"

    return {
        "zn_zn_distance":        dist,
        "binuclear_distance_ok": dist_ok,
        "bridging_ligand":       bridging,
        "binuclear_call":        binuclear_call,
    }


def get_protein_overall_call(per_zn, binuclear_call):
    """
    Single protein-level verdict integrating both per-site calls and
    the binuclear assessment.

    A functional M17 LAP requires:
      - Both zinc sites properly coordinated
      - Both zincs in the same pocket (binuclear)

    Decision logic:
      SEPARATE_POCKETS                             → LIKELY_INACTIVE
      BINUCLEAR_OK + both CANONICAL                → LIKELY_ACTIVE
      BINUCLEAR_OK + one NON_CANONICAL             → LIKELY_INACTIVE
      BINUCLEAR_OK + otherwise                     → UNCERTAIN
      BORDERLINE or no_shared_ligand + both CANONICAL → UNCERTAIN
      BORDERLINE or no_shared_ligand + otherwise   → LIKELY_INACTIVE
      single_zinc / no_zinc                        → LIKELY_INACTIVE
    """
    calls = [z["call"] for z in per_zn]

    if binuclear_call in ("no_zinc", "single_zinc_only"):
        return "LIKELY_INACTIVE"

    if binuclear_call == "SEPARATE_POCKETS":
        return "LIKELY_INACTIVE"

    if binuclear_call == "BINUCLEAR_OK":
        if all(c == "CANONICAL" for c in calls):
            return "LIKELY_ACTIVE"
        elif any(c == "NON_CANONICAL" for c in calls):
            return "LIKELY_INACTIVE"
        else:
            return "UNCERTAIN"

    if binuclear_call in ("BORDERLINE_SEPARATION", "BINUCLEAR_OK_no_shared_ligand"):
        if all(c == "CANONICAL" for c in calls):
            return "UNCERTAIN"
        else:
            return "LIKELY_INACTIVE"

    return "UNCERTAIN"


def fmt_contacts(contacts):
    return "; ".join(
        f"{c['chain']}:{c['residue']}{c['res_num']}({c['atom_name']})={c['distance']}A"
        + ("*" if c["canonical_sidechain"] else "")
        for c in contacts
    )


# ── MAIN LOOP ─────────────────────────────────────────────────────────────────
results   = []
cif_files = sorted([f for f in os.listdir(CIF_DIR) if f.endswith(".cif")])
print(f"Found {len(cif_files)} CIF files to process.\n")

for cif_file in cif_files:
    protein_name = cif_file.replace(".cif", "")
    filepath     = os.path.join(CIF_DIR, cif_file)
    print(f"Processing: {protein_name}")

    run(session, f"open {filepath}")
    atoms        = all_atoms(session)
    zn_atoms     = [a for a in atoms if a.element.name == "Zn"]
    num_zn_total = len(zn_atoms)

    if not zn_atoms:
        print(f"  No zinc found")
        results.append({
            "protein": protein_name, "zinc_found": False,
            "zn_index": None, "num_zn_total": 0, "site_type": None,
            "num_canonical_contacts": None,
            "coordinating_residues_raw": None,
            "coordinating_residues_deduped": None,
            "canonical_residue_types": None, "lys_present": False,
            "mean_distance": None, "mean_plddt_at_site": None,
            "score": 0, "call": "NON_CANONICAL", "flags": "no_zinc_found",
            "zn_zn_distance": None, "binuclear_distance_ok": None,
            "bridging_ligand": None, "binuclear_call": "no_zinc",
            "overall_call": "LIKELY_INACTIVE",
        })
        run(session, "close all")
        continue

    # ── Per-zinc scoring ──────────────────────────────────────────────────────
    per_zn = []

    for i, zn in enumerate(zn_atoms):
        zn_pos           = np.array(zn.scene_coord)
        coordinating_raw = []

        for atom in atoms:
            if atom.element.name == "Zn":
                continue
            if atom.element.name not in ("N", "O", "S"):
                continue
            atom_pos = np.array(atom.scene_coord)
            dist     = float(np.linalg.norm(zn_pos - atom_pos))
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

        per_zn.append({
            "zn_index":   i,
            "zn_pos":     zn_pos,
            "deduped":    deduped,
            "raw":        coordinating_raw,
            "score":      score,
            "site":       site,
            "call":       call,
            "flags":      flags,
            "mean_plddt": mean_plddt,
        })

    # ── Binuclear + protein-level overall call ────────────────────────────────
    binuclear    = assess_binuclear_site(per_zn)
    overall_call = get_protein_overall_call(per_zn, binuclear["binuclear_call"])

    print(f"  Zn-Zn: {binuclear['zn_zn_distance']} A | "
          f"binuclear: {binuclear['binuclear_call']} | "
          f"bridging: {binuclear['bridging_ligand']} | "
          f"OVERALL: {overall_call}")

    # ── Write one row per zinc (overall_call identical on both rows) ──────────
    for entry in per_zn:
        deduped         = entry["deduped"]
        canon_deduped   = [c for c in deduped if c["canonical_sidechain"]]
        canon_res_types = sorted({c["residue"] for c in canon_deduped})
        lys_present     = "LYS" in canon_res_types
        distances       = [c["distance"] for c in deduped]
        mean_dist       = round(sum(distances) / len(distances), 3) if distances else None
        flag_str        = " | ".join(entry["flags"]) if entry["flags"] else "none"

        print(f"  Zn[{entry['zn_index']+1}] ({entry['site']}): "
              f"score={entry['score']}/10 ({entry['call']}) | "
              f"lys={lys_present} | plddt={entry['mean_plddt']}")

        results.append({
            "protein":                       protein_name,
            "zinc_found":                    True,
            "zn_index":                      entry["zn_index"] + 1,
            "num_zn_total":                  num_zn_total,
            "site_type":                     entry["site"],
            "num_canonical_contacts":        len(canon_deduped),
            "coordinating_residues_raw":     fmt_contacts(entry["raw"]),
            "coordinating_residues_deduped": fmt_contacts(deduped),
            "canonical_residue_types":       ", ".join(canon_res_types),
            "lys_present":                   lys_present,
            "mean_distance":                 mean_dist,
            "mean_plddt_at_site":            entry["mean_plddt"],
            "score":                         entry["score"],
            "call":                          entry["call"],
            "flags":                         flag_str,
            "zn_zn_distance":                binuclear["zn_zn_distance"],
            "binuclear_distance_ok":         binuclear["binuclear_distance_ok"],
            "bridging_ligand":               binuclear["bridging_ligand"],
            "binuclear_call":                binuclear["binuclear_call"],
            "overall_call":                  overall_call,  # same for both zinc rows
        })

    run(session, "close all")

# ── WRITE CSV ─────────────────────────────────────────────────────────────────
fieldnames = [
    "protein", "zinc_found", "zn_index", "num_zn_total", "site_type",
    "num_canonical_contacts", "coordinating_residues_raw",
    "coordinating_residues_deduped", "canonical_residue_types", "lys_present",
    "mean_distance", "mean_plddt_at_site", "score", "call", "flags",
    "zn_zn_distance", "binuclear_distance_ok", "bridging_ligand",
    "binuclear_call", "overall_call",
]

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

# ── SUMMARY ───────────────────────────────────────────────────────────────────
# One row per protein for summary (use zn_index == 1 to avoid double-counting)
protein_rows  = [r for r in results if r["zinc_found"] and r["zn_index"] == 1]
overall_calls = [r["overall_call"] for r in protein_rows]
bi_calls      = [r["binuclear_call"] for r in protein_rows]
per_calls     = [r["call"] for r in results if r["zinc_found"]]

print(f"\n{'='*60}")
print(f"Done. {len(cif_files)} structures processed.")
print(f"\nPer-site calls (all zinc ions):")
print(f"  CANONICAL:       {per_calls.count('CANONICAL')}")
print(f"  UNCERTAIN:       {per_calls.count('UNCERTAIN')}")
print(f"  NON_CANONICAL:   {per_calls.count('NON_CANONICAL')}")
print(f"\nBinuclear assessment (per protein):")
print(f"  BINUCLEAR_OK:                  {bi_calls.count('BINUCLEAR_OK')}")
print(f"  BINUCLEAR_OK_no_shared_ligand: {bi_calls.count('BINUCLEAR_OK_no_shared_ligand')}")
print(f"  BORDERLINE_SEPARATION:         {bi_calls.count('BORDERLINE_SEPARATION')}")
print(f"  SEPARATE_POCKETS:              {bi_calls.count('SEPARATE_POCKETS')}")
print(f"\nOverall call (protein level):")
print(f"  LIKELY_ACTIVE:   {overall_calls.count('LIKELY_ACTIVE')}")
print(f"  UNCERTAIN:       {overall_calls.count('UNCERTAIN')}")
print(f"  LIKELY_INACTIVE: {overall_calls.count('LIKELY_INACTIVE')}")
print(f"\nResults written to: {OUTPUT_CSV}")
