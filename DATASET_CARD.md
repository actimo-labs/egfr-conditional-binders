---
license: odc-by
task_categories:
  - tabular-regression
tags:
  - biology
  - protein
  - structural-biology
  - protein-design
  - epitope-mapping
  - egfr
pretty_name: EGFR Domain III Titratable Anchor Map
size_categories:
  - n<1K
configs:
  - config_name: residues
    data_files: residues.csv
  - config_name: anchors
    data_files: anchors.csv
---

# EGFR Domain III Titratable Anchor Map

Residue-level annotation of the human EGFR domain III epitope face, built for designing **pH-conditional** binders — ones that engage at pH 6.5 and release at pH 7.4.

Derived from PDB **6ARU** chain A (EGFR ectodomain, cetuximab-Fab-bound, 3.2 Å). Numbering is UniProt **P00533-1** precursor numbering throughout.

## Why this exists

Designing an acidic-ON binder means engineering titratable contacts that gain favourable interaction when protonated. The linkage ceiling is `10^ΔpH` ≈ 8× per site over 6.5 → 7.4, so clearing a "no detectable binding" bar needs **three titratable pairs minimum**. Finding three that a single binder can span simultaneously turns out to be the binding constraint, and nothing published maps them.

This dataset is that map.

## Contents

**`residues.csv`** — every domain III residue

| column | meaning |
|---|---|
| `precursor_num` | UniProt P00533-1 numbering |
| `resname` | three-letter code |
| `rel_sasa` | relative solvent accessibility (Shrake-Rupley / Tien et al. 2013 max-ASA) |
| `is_epitope` | within 5 Å of the bound cetuximab Fab |
| `min_partner_dist` | min heavy-atom distance to Fab, Å |
| `conserved_mouse` | identical in mouse EGFR (UniProt Q01279) |
| `pka` | PROPKA3 predicted pKa |

**`anchors.csv`** — surface titratable groups eligible as design anchors

| class | target residue | complementary binder residue | switch sits on |
|---|---|---|---|
| **A** | Asp / Glu | His | binder |
| **B** | His | Asp / Glu | target |

Class B is the better anchor where available: the switching residue is on EGFR, so it is present in both assay conditions by construction, and it spares the binder from carrying four histidines (desolvation cost, aggregation risk, pI drift). Domain III carries H370 and H433, the residues the cross-reactive pH-dependent antibody 14C07/G5V2 exploits.

**`patches.json`** — anchor sets within mutual spanning distance, ranked. Flags `meets_spec_v2` for patches with ≥3 anchors, ≥1 Class B, all mouse-conserved. Includes side-chain tip coordinates and outward vectors for placing the complementary binder residue.

## Headline numbers

| | |
|---|---|
| Residues annotated | 609 |
| On the cetuximab contact face | 26 |
| Identity to mouse EGFR | 88.8% (published: 88%) |
| Surface titratable anchors in domain III | 19 (16 Class A, 3 Class B) |
| Patches meeting ≥3 pairs, ≥1 Class B, all conserved | 4 |

**HIS433** — pKa 6.22, relative SASA 0.648, identical in mouse, inside the cetuximab contact footprint. ~37% protonated at pH 6.5 against ~6% at 7.4. The best available anchor on the protein, and the only one servicing all three design objectives at once.

**HIS358** — pKa 6.26, conserved, but outside the Fab footprint.

**HIS383** — excluded. Arginine in mouse, so titratability is lost entirely.

## How it was built

[`epitope_prep.py`](https://github.com/) — Biopython for structure handling and SASA, PROPKA3 for pKa, global pairwise alignment against the mouse ortholog for conservation. Epitope face derived from the bound Fab rather than a literature residue list.

Using cetuximab to *define* the epitope is permitted under the Adaptyv competition's de novo rule. No Fab coordinates or sequence enter any design.

## Caveats

- pKa values are PROPKA3 predictions on a single crystal conformation. Treat as a tier-1 screen, not ground truth — a static structure cannot see a salt bridge that dissociates in 5 ns of MD.
- Conservation comes from sequence alignment, not a structural superposition of the mouse ortholog.
- `rel_sasa` is computed on the isolated chain, so residues buried by domains II/IV in the tethered conformation will read as more exposed than they are in the full ectodomain.

## Citation

If this is useful, cite the structure (Sickmier et al., PDB 6ARU) and PROPKA3 (Olsson et al., JCTC 2011).