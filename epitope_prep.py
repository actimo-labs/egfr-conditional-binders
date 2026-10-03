#!/usr/bin/env python3
"""
epitope_prep.py — EGFR domain III epitope prep and titratable-anchor mapper.

Stage 0 of the Anthropic x Adaptyv 2026 EGFR conditional-binder pipeline.

Takes the competition reference structure (PDB 6ARU, chain A) and emits the
hotspot spec that BindCraft / RFdiffusion / ProteinMPNN all consume, annotated
with everything the pH gate needs:

  - per-residue relative SASA (Shrake-Rupley, Tien et al. 2013 max-ASA)
  - epitope face, derived from contacts with a reference partner chain
    (cetuximab Fab H/L in 6ARU) or from an explicit residue list
  - human/mouse conservation, for the cross-reactivity objective
  - PROPKA3 pKa for every titratable group
  - Class A anchors (target Asp/Glu -> binder His)
  - Class B anchors (target His    -> binder Asp/Glu)   <- preferred, switch on target
  - anchor PATCHES: mutually-spannable sets of >=3 anchors, since one binder
    must engage them all (spec v2 section 5.1)

Using the cetuximab chains only to DEFINE the epitope is permitted under the
competition's de novo rule. No coordinates or sequence from the Fab are carried
into any design.

Usage
-----
    # get the structure and the mouse ortholog first (network-restricted
    # containers cannot reach RCSB/UniProt; run these on your own machine):
    wget https://files.rcsb.org/download/6ARU.pdb
    wget -O mouse_egfr.fasta "https://rest.uniprot.org/uniprotkb/Q01279.fasta"

    python epitope_prep.py 6ARU.pdb \
        --target-chain A \
        --partner-chains H L \
        --mouse-fasta mouse_egfr.fasta \
        --outdir out/

Outputs (in --outdir)
---------------------
    target_clean.pdb      target chain, standard residues only
    residues.csv          every domain III residue, fully annotated
    anchors.csv           titratable anchors, ranked
    patches.json          spannable anchor sets, ranked  <- the design spec
    hotspots.json         BindCraft / RFdiffusion target settings
    mpnn_bias.json        ProteinMPNN per-position bias template
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, asdict, field
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------

# Tien et al. 2013 (PLoS ONE 8:e80635) theoretical max ASA, angstrom^2
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}

# EGFR precursor numbering (UniProt P00533-1): signal peptide is 1-24,
# so mature residue n == precursor residue n + 24.
PRECURSOR_OFFSET = 24

# side-chain tip atoms used to place the complementary binder residue
TIP_ATOMS = {
    "ASP": ("OD1", "OD2", "CG"),
    "GLU": ("OE1", "OE2", "CD"),
    "HIS": ("ND1", "NE2", "CE1"),
}

CLASS_A_RES = {"ASP", "GLU"}   # target anion  -> binder His
CLASS_B_RES = {"HIS"}          # target His    -> binder Asp/Glu


# --------------------------------------------------------------------------
# data model
# --------------------------------------------------------------------------

@dataclass
class Residue:
    resnum: int                 # as numbered in the input file
    precursor_num: int          # UniProt P00533-1 numbering
    resname: str
    one_letter: str
    chain: str
    rel_sasa: float = 0.0
    abs_sasa: float = 0.0
    in_domain3: bool = False
    is_epitope: bool = False
    min_partner_dist: float = float("inf")
    conserved_mouse: str = "unknown"   # "yes" | "no" | "gap" | "unknown"
    mouse_resname: str = ""
    pka: float | None = None
    model_pka: float | None = None
    anchor_class: str = ""             # "A" | "B" | ""
    anchor_score: float = 0.0
    tip_xyz: list[float] = field(default_factory=list)
    out_vec: list[float] = field(default_factory=list)


# --------------------------------------------------------------------------
# structure handling
# --------------------------------------------------------------------------

class StandardAA(Select):
    """Keep only standard amino acids of one chain. Drops waters, glycans, ions."""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return is_aa(residue, standard=True) and residue.id[0] == " "

    def accept_atom(self, atom):
        # drop altloc B+ and hydrogens (PROPKA re-protonates anyway)
        return atom.element != "H" and (
            not atom.is_disordered() or atom.get_altloc() in (" ", "A")
        )


def load_structure(pdb_path: Path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure("target", str(pdb_path))


def write_clean_chain(structure, chain_id: str, out_path: Path) -> Path:
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), select=StandardAA(chain_id))
    return out_path


def chain_residues(structure, chain_id: str) -> list:
    model = next(structure.get_models())
    if chain_id not in model:
        raise SystemExit(f"chain {chain_id!r} not in structure; found {[c.id for c in model]}")
    return [r for r in model[chain_id] if is_aa(r, standard=True) and r.id[0] == " "]


# --------------------------------------------------------------------------
# SASA
# --------------------------------------------------------------------------

def compute_sasa(clean_pdb: Path) -> dict[int, float]:
    """Shrake-Rupley SASA on the ISOLATED target chain."""
    struct = load_structure(clean_pdb)
    ShrakeRupley().compute(struct[0], level="R")
    return {
        r.id[1]: float(r.sasa)
        for r in struct[0].get_residues()
        if is_aa(r, standard=True)
    }


# --------------------------------------------------------------------------
# epitope definition
# --------------------------------------------------------------------------

def partner_contacts(
    structure, target_chain: str, partner_chains: list[str], cutoff: float
) -> dict[int, float]:
    """Min heavy-atom distance from each target residue to any partner chain atom."""
    model = next(structure.get_models())
    partner_atoms = []
    for cid in partner_chains:
        if cid not in model:
            print(f"  ! partner chain {cid!r} absent, skipping", file=sys.stderr)
            continue
        partner_atoms += [
            a for a in model[cid].get_atoms()
            if a.element != "H" and a.get_parent().id[0] == " "
        ]
    if not partner_atoms:
        return {}

    pcoords = np.array([a.coord for a in partner_atoms])
    out = {}
    for res in chain_residues(structure, target_chain):
        rc = np.array([a.coord for a in res if a.element != "H"])
        if rc.size == 0:
            continue
        d = np.linalg.norm(rc[:, None, :] - pcoords[None, :, :], axis=-1).min()
        out[res.id[1]] = float(d)
    return out


# --------------------------------------------------------------------------
# conservation
# --------------------------------------------------------------------------

def read_fasta(path: Path) -> str:
    seq = []
    for line in path.read_text().splitlines():
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq)


def map_conservation(residues: list[Residue], mouse_seq: str) -> None:
    """Global pairwise align human chain sequence to the mouse ortholog."""
    human_seq = "".join(r.one_letter for r in residues)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = None
    aligner.match_score = 2
    aligner.mismatch_score = -1

    aln = aligner.align(human_seq, mouse_seq)[0]
    h_aln, m_aln = aln[0], aln[1]

    hi = 0
    for hc, mc in zip(h_aln, m_aln):
        if hc == "-":
            continue
        res = residues[hi]
        if mc == "-":
            res.conserved_mouse = "gap"
        else:
            res.mouse_resname = mc
            res.conserved_mouse = "yes" if mc == hc else "no"
        hi += 1


# --------------------------------------------------------------------------
# PROPKA
# --------------------------------------------------------------------------

def run_propka(clean_pdb: Path) -> dict[int, tuple[float, float]]:
    """Run PROPKA3, parse the SUMMARY block. Returns {resnum: (pKa, model_pKa)}."""
    exe = shutil.which("propka3") or shutil.which("propka")
    if exe is None:
        print("  ! propka3 not on PATH, skipping pKa annotation", file=sys.stderr)
        return {}

    with tempfile.TemporaryDirectory() as td:
        work = Path(td) / clean_pdb.name
        shutil.copy(clean_pdb, work)
        try:
            subprocess.run(
                [exe, work.name], cwd=td, check=True,
                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=900,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            err = getattr(e, "stderr", b"") or b""
            if isinstance(err, bytes):
                err = err.decode("utf-8", "replace")
            print("  ! propka failed, skipping pKa annotation", file=sys.stderr)
            for line in err.strip().splitlines()[-8:]:
                print(f"      {line}", file=sys.stderr)
            return {}

        pka_file = Path(td) / (work.stem + ".pka")
        if not pka_file.exists():
            return {}
        return parse_pka(pka_file.read_text())


def parse_pka(text: str) -> dict[int, tuple[float, float]]:
    out: dict[int, tuple[float, float]] = {}
    in_summary = False
    for line in text.splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if in_summary:
            if line.strip().startswith("-") or "Free energy" in line:
                break
            parts = line.split()
            # e.g.  ASP  15 A     3.61       3.80
            if len(parts) >= 5 and parts[0].isalpha() and parts[1].lstrip("-").isdigit():
                try:
                    out[int(parts[1])] = (float(parts[3]), float(parts[4]))
                except ValueError:
                    continue
    return out


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------

def tip_geometry(res) -> tuple[list[float], list[float]]:
    """Side-chain tip centroid and outward unit vector (tip - CA)."""
    names = TIP_ATOMS.get(res.get_resname(), ())
    pts = [res[n].coord for n in names if n in res]
    if not pts or "CA" not in res:
        return [], []
    tip = np.mean(pts, axis=0)
    vec = tip - res["CA"].coord
    n = np.linalg.norm(vec)
    if n < 1e-6:
        return [round(float(x), 3) for x in tip], []
    return (
        [round(float(x), 3) for x in tip],
        [round(float(x), 3) for x in (vec / n)],
    )


# --------------------------------------------------------------------------
# anchor scoring
# --------------------------------------------------------------------------

def score_anchor(r: Residue, args) -> float:
    """0-100. Exposure + epitope proximity + conservation + pKa window."""
    score = 0.0

    # exposure, 0-30
    score += 30.0 * min(r.rel_sasa / 0.50, 1.0)

    # epitope proximity, 0-30
    if r.min_partner_dist < float("inf"):
        score += 30.0 * max(0.0, 1.0 - (r.min_partner_dist / (args.epitope_cutoff * 2)))
    elif r.is_epitope:
        score += 30.0

    # mouse conservation, 0-30. Hard requirement for objective 2.
    if r.conserved_mouse == "yes":
        score += 30.0
    elif r.conserved_mouse == "conservative":   # D<->E, charge preserved
        score += 26.0
    elif r.conserved_mouse == "unknown":
        score += 12.0

    # pKa window, 0-10. We want groups that actually titrate near the assay window.
    if r.pka is not None:
        if r.anchor_class == "B":          # target His, want it protonatable at 6.5
            score += 10.0 * max(0.0, 1.0 - abs(r.pka - 6.5) / 2.0)
        else:                              # target carboxylate, want it ionised at 6.5
            score += 10.0 if r.pka < 5.5 else 10.0 * max(0.0, 1.0 - (r.pka - 5.5) / 2.0)

    return round(score, 1)


def find_patches(anchors: list[Residue], span: float, min_size: int,
                 epitope_reach: float = 10.0) -> list[dict]:
    """Greedy: anchor sets whose tips all sit within `span` A of each other."""
    usable = [a for a in anchors if a.tip_xyz]
    if len(usable) < min_size:
        return []

    coords = np.array([a.tip_xyz for a in usable])
    dmat = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)

    patches, seen = [], set()
    for i, seed in enumerate(usable):
        # members within span of the seed, best-scoring first
        cand = sorted(
            [j for j in range(len(usable)) if dmat[i, j] <= span],
            key=lambda j: -usable[j].anchor_score,
        )
        # enforce mutual span, not just span-to-seed
        members: list[int] = []
        for j in cand:
            if all(dmat[j, k] <= span for k in members):
                members.append(j)
        if len(members) < min_size:
            continue

        key = tuple(sorted(members))
        if key in seen:
            continue
        seen.add(key)

        mem = [usable[j] for j in members]
        n_b = sum(1 for m in mem if m.anchor_class == "B")
        n_cons = sum(1 for m in mem if m.conserved_mouse in ("yes", "conservative"))
        n_epi = sum(1 for m in mem if m.min_partner_dist <= epitope_reach)
        mean_score = float(np.mean([m.anchor_score for m in mem]))

        patches.append({
            "n_anchors": len(mem),
            "n_class_B": n_b,
            "n_conserved": n_cons,
            "n_epitope_anchors": n_epi,
            "max_span_A": round(float(dmat[np.ix_(members, members)].max()), 2),
            "mean_anchor_score": round(mean_score, 1),
            # spec v2: >=3 pairs, >=1 class B, all conserved
            "meets_spec_v2": len(mem) >= 3 and n_b >= 1 and n_cons == len(mem),
            "centroid": [round(float(x), 3) for x in coords[members].mean(axis=0)],
            "anchors": [
                {
                    "precursor_num": m.precursor_num,
                    "file_num": m.resnum,
                    "resname": m.resname,
                    "class": m.anchor_class,
                    "binder_partner": "HIS" if m.anchor_class == "A" else "ASP/GLU",
                    "rel_sasa": round(m.rel_sasa, 3),
                    "pka": m.pka,
                    "conserved_mouse": m.conserved_mouse,
                    "score": m.anchor_score,
                    "tip_xyz": m.tip_xyz,
                    "out_vec": m.out_vec,
                }
                for m in sorted(mem, key=lambda x: -x.anchor_score)
            ],
        })

    # Epitope coverage outranks anchor count. A seven-anchor patch on the back
    # face of domain III is worth less than a three-anchor patch sitting in the
    # functional epitope: more anchors cannot rescue the wrong surface.
    patches.sort(
        key=lambda p: (p["meets_spec_v2"], p["n_epitope_anchors"] > 0,
                       p["n_class_B"], p["mean_anchor_score"],
                       p["n_epitope_anchors"], p["n_anchors"]),
        reverse=True,
    )
    return patches


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb", type=Path, help="input PDB (e.g. 6ARU.pdb)")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--partner-chains", nargs="*", default=["H", "L"],
                    help="chains used ONLY to define the epitope face (cetuximab Fab in 6ARU)")
    ap.add_argument("--epitope-cutoff", type=float, default=5.0,
                    help="heavy-atom distance for epitope contact, A")
    ap.add_argument("--epitope-residues", nargs="*", type=int, default=None,
                    help="explicit epitope residues (precursor numbering); overrides partner contacts")
    ap.add_argument("--domain3", nargs=2, type=int, default=[310, 480],
                    metavar=("START", "END"), help="domain III window, precursor numbering")
    ap.add_argument("--numbering", choices=["mature", "precursor"], default="mature",
                    help="numbering used in the input PDB (6ARU is mature)")
    ap.add_argument("--mouse-fasta", type=Path, default=None,
                    help="mouse EGFR FASTA (UniProt Q01279)")
    ap.add_argument("--min-rel-sasa", type=float, default=0.20,
                    help="surface-exposure cutoff for anchor eligibility")
    ap.add_argument("--patch-span", type=float, default=25.0,
                    help="max tip-to-tip distance within one anchor patch, A")
    ap.add_argument("--min-patch", type=int, default=3,
                    help="minimum anchors per patch (spec v2: 3)")
    ap.add_argument("--outdir", type=Path, default=Path("out"))
    args = ap.parse_args()

    if not args.pdb.exists():
        raise SystemExit(f"not found: {args.pdb}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    offset = PRECURSOR_OFFSET if args.numbering == "mature" else 0
    d3_lo, d3_hi = args.domain3

    print(f"[1] loading {args.pdb}")
    structure = load_structure(args.pdb)

    print(f"[2] cleaning chain {args.target_chain}")
    clean = write_clean_chain(structure, args.target_chain, args.outdir / "target_clean.pdb")

    print("[3] SASA")
    sasa = compute_sasa(clean)

    print("[4] epitope face")
    if args.epitope_residues:
        contacts = {}
        epitope_set = set(args.epitope_residues)
    else:
        contacts = partner_contacts(structure, args.target_chain,
                                    args.partner_chains, args.epitope_cutoff)
        epitope_set = set()

    residues: list[Residue] = []
    for res in chain_residues(structure, args.target_chain):
        rn = res.id[1]
        three = res.get_resname()
        one = THREE_TO_ONE.get(three.capitalize(), THREE_TO_ONE.get(three, "X"))
        pnum = rn + offset
        abs_s = sasa.get(rn, 0.0)
        rel_s = abs_s / MAX_ASA.get(one, 200.0)
        dist = contacts.get(rn, float("inf"))
        tip, vec = tip_geometry(res)

        residues.append(Residue(
            resnum=rn, precursor_num=pnum, resname=three, one_letter=one,
            chain=args.target_chain,
            abs_sasa=round(abs_s, 2), rel_sasa=round(rel_s, 3),
            in_domain3=(d3_lo <= pnum <= d3_hi),
            is_epitope=(pnum in epitope_set) or (dist <= args.epitope_cutoff),
            min_partner_dist=round(dist, 2) if dist < float("inf") else float("inf"),
            tip_xyz=tip, out_vec=vec,
        ))

    n_epi = sum(1 for r in residues if r.is_epitope)
    print(f"    {len(residues)} residues, {n_epi} on the epitope face")

    print("[5] mouse conservation")
    if args.mouse_fasta and args.mouse_fasta.exists():
        map_conservation(residues, read_fasta(args.mouse_fasta))
        cons = sum(1 for r in residues if r.conserved_mouse == "yes")
        print(f"    {cons}/{len(residues)} identical to mouse "
              f"({100*cons/max(len(residues),1):.1f}%)")
    else:
        print("    ! no --mouse-fasta; objective 2 cannot be gated. "
              "Fetch: https://rest.uniprot.org/uniprotkb/Q01279.fasta")

    print("[6] PROPKA3")
    pkas = run_propka(clean)
    for r in residues:
        if r.resnum in pkas:
            r.pka, r.model_pka = pkas[r.resnum]
    print(f"    {len(pkas)} titratable groups")

    print("[7] anchors")
    anchors: list[Residue] = []
    for r in residues:
        if not r.in_domain3 or r.rel_sasa < args.min_rel_sasa:
            continue
        if r.resname in CLASS_A_RES:
            r.anchor_class = "A"
        elif r.resname in CLASS_B_RES:
            r.anchor_class = "B"
        else:
            continue

        # D<->E across species preserves the anion, so a class A anchor is
        # functionally conserved even though the residues differ. Class B has
        # no equivalent: His->anything loses titratability outright.
        if (r.anchor_class == "A" and r.conserved_mouse == "no"
                and r.mouse_resname in ("D", "E")):
            r.conserved_mouse = "conservative"

        r.anchor_score = score_anchor(r, args)
        anchors.append(r)

    anchors.sort(key=lambda x: -x.anchor_score)
    n_a = sum(1 for a in anchors if a.anchor_class == "A")
    n_b = sum(1 for a in anchors if a.anchor_class == "B")
    print(f"    {len(anchors)} anchors  (class A {n_a} -> binder His, "
          f"class B {n_b} -> binder Asp/Glu)")

    print("[8] patches")
    patches = find_patches(anchors, args.patch_span, args.min_patch,
                           epitope_reach=args.epitope_cutoff * 2)
    viable = [p for p in patches if p["meets_spec_v2"]]
    print(f"    {len(patches)} patches, {len(viable)} meet spec v2 "
          f"(>=3 anchors, >=1 class B, all mouse-conserved)")

    # ---------------- outputs ----------------
    with open(args.outdir / "residues.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(residues[0]).keys()))
        w.writeheader()
        for r in residues:
            row = asdict(r)
            row["tip_xyz"] = ";".join(map(str, r.tip_xyz))
            row["out_vec"] = ";".join(map(str, r.out_vec))
            if row["min_partner_dist"] == float("inf"):
                row["min_partner_dist"] = ""
            w.writerow(row)

    with open(args.outdir / "anchors.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["precursor_num", "file_num", "resname", "class", "binder_partner",
                    "rel_sasa", "pka", "conserved_mouse", "mouse_res",
                    "min_partner_dist", "score"])
        for a in anchors:
            w.writerow([a.precursor_num, a.resnum, a.resname, a.anchor_class,
                        "HIS" if a.anchor_class == "A" else "ASP/GLU",
                        a.rel_sasa, a.pka, a.conserved_mouse, a.mouse_resname,
                        "" if a.min_partner_dist == float("inf") else a.min_partner_dist,
                        a.anchor_score])

    (args.outdir / "patches.json").write_text(json.dumps({
        "spec": {
            "min_anchors": args.min_patch,
            "min_class_B": 1,
            "require_all_mouse_conserved": True,
            "patch_span_A": args.patch_span,
        },
        "n_patches": len(patches),
        "n_meeting_spec": len(viable),
        "patches": patches[:25],
    }, indent=2))

    # BindCraft / RFdiffusion target settings. Hotspots in FILE numbering, which is
    # what those tools read off the PDB.
    top = patches[0] if patches else None
    hot = [a["file_num"] for a in top["anchors"]] if top else \
          [a.resnum for a in anchors[:6]]
    epi = sorted(r.resnum for r in residues if r.is_epitope and r.in_domain3)

    (args.outdir / "hotspots.json").write_text(json.dumps({
        "starting_pdb": str((args.outdir / "target_clean.pdb").resolve()),
        "chains": args.target_chain,
        "target_hotspot_residues": ",".join(f"{args.target_chain}{n}" for n in hot),
        "epitope_residues_file_numbering": epi,
        "numbering_note": f"file numbering; precursor = file + {offset}",
        "anchor_patch": top,
    }, indent=2))

    # ProteinMPNN bias template: which binder positions to push toward the
    # complementary residue. Binder indices are filled in after backbone generation.
    (args.outdir / "mpnn_bias.json").write_text(json.dumps({
        "_comment": "Per-position bias for the BINDER chain. Map each target anchor to "
                    "the nearest binder position after backbone generation, then apply "
                    "the bias below at that position. Feed via --bias_by_res_jsonl.",
        "bias_recipe": [
            {
                "target_anchor": f"{a['resname']}{a['precursor_num']}",
                "class": a["class"],
                "target_tip_xyz": a["tip_xyz"],
                "target_out_vec": a["out_vec"],
                "binder_bias": ({"H": 3.0} if a["class"] == "A"
                                else {"D": 2.5, "E": 2.5}),
                "rationale": ("protonated binder His salt-bridges target carboxylate at pH 6.5"
                              if a["class"] == "A" else
                              "binder carboxylate salt-bridges protonated target His at pH 6.5"),
            }
            for a in (top["anchors"] if top else [])
        ],
    }, indent=2))

    # ---------------- report ----------------
    print("\n  top anchors")
    print(f"  {'res':>10}  {'cls':>3}  {'relSASA':>7}  {'pKa':>5}  {'mouse':>7}  {'score':>5}")
    for a in anchors[:12]:
        pk = f"{a.pka:.2f}" if a.pka is not None else "  -  "
        print(f"  {a.resname}{a.precursor_num:<6}  {a.anchor_class:>3}  "
              f"{a.rel_sasa:>7.3f}  {pk:>5}  {a.conserved_mouse:>7}  {a.anchor_score:>5.1f}")

    if viable:
        p = viable[0]
        names = ", ".join(f"{a['resname']}{a['precursor_num']}({a['class']})"
                          for a in p["anchors"])
        print(f"\n  best spec-compliant patch: {p['n_anchors']} anchors, "
              f"{p['n_class_B']} class B, {p['n_epitope_anchors']} on epitope, "
              f"span {p['max_span_A']} A")
        print(f"    {names}")
    else:
        print("\n  ! no patch meets spec v2. Loosen --patch-span or --min-rel-sasa, "
              "or cut the pH allocation per spec v2 section 10.")

    print(f"\n  written to {args.outdir.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())