# egfr-conditional-binders

Design pipeline for the [Anthropic × Adaptyv 2026 EGFR challenge](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr) — de novo binders that bind human EGFR at pH 6.5, release at pH 7.4, and cross-react with mouse.

Work in progress, built in the open. Results from the wet lab get published under ODC-BY either way, so there's no reason not to show the working.

## The constraint that shapes everything

For one titratable site, linkage gives:

```
Ka(pH) ∝ (1 + 10^(pKa_bound − pH)) / (1 + 10^(pKa_free − pH))
```

Push `pKa_bound` high and `pKa_free` low and the ratio across the assay window collapses to **10^ΔpH**. Over 6.5 → 7.4 that is **~8× maximum, per site**. No tuning beats it.

| Titratable pairs | Max shift | ΔΔG |
|---|---|---|
| 1 | ~8× | 1.2 kcal/mol |
| 2 | ~63× | 2.5 kcal/mol |
| 3 | ~500× | 3.7 kcal/mol |
| 4 | ~4000× | 4.9 kcal/mol |

"No detectable binding at pH 7.4" needs ≥100×. So the spec is **three titratable pairs minimum, four preferred**. One engineered histidine is thermodynamically incapable of clearing the bar, which is the predictable way to lose this challenge.

This is why FcRn uses three His–carboxylate pairs and not one.

## Two consequences

**The switch is ΔpKa, not the presence of His.** A histidine at pKa 6.5 free and 6.5 bound gives zero pH dependence no matter how good the salt bridge looks in the model. Target `pKa_free ≤ 6.0`, `pKa_bound ≥ 8.0`, net proton uptake on binding ≥ 2.

**Generative models are pH-blind.** ProteinMPNN, RFdiffusion, AF2/AF3 and Boltz all use a 20-token alphabet, output heavy atoms only, carry no HID/HIE/HIP distinction and take no pH input. Whatever they encode about His–carboxylate geometry is a pH-averaged smear over PDB entries crystallised anywhere from pH 4.5 to 8.5. So don't ask the generator to discover pH selectivity — force the placement, and let physics (PROPKA → Rosetta → MD) do all the discrimination.

## What's here now

`epitope_prep.py` — epitope prep and titratable-anchor mapper. Takes PDB 6ARU and emits the hotspot spec that BindCraft, RFdiffusion and ProteinMPNN consume, annotated with everything the pH gate needs.

- Per-residue relative SASA (Shrake-Rupley, Tien et al. max-ASA)
- Epitope face derived from the bound cetuximab Fab chains rather than a guessed residue list
- Human/mouse conservation, for the cross-reactivity objective
- PROPKA3 pKa on every titratable group
- **Class A** anchors: target Asp/Glu → binder His
- **Class B** anchors: target His → binder Asp/Glu — preferred, because the switching residue sits on the target and is therefore present in both assay conditions by construction
- **Patch finding**: anchor sets that one binder can actually span simultaneously, with mutual distance enforced rather than distance-to-seed

That last part is the bit that bites. Individual anchors are easy. Finding three that a single binder can engage at once is the real constraint, and if the structure doesn't offer them, the honest move is to cut the pH allocation rather than ship designs that can't work.

Using the cetuximab chains to *define* the epitope is permitted under the competition's de novo rule. No Fab coordinates or sequence propagate into any design.

## What it found on 6ARU

609 residues in the EGFR ectodomain chain, 26 on the cetuximab contact face, **88.8% identity to mouse** — matching the published 88% for the ectodomain, which is a decent check that the alignment is doing what it should.

19 surface titratable anchors in domain III: 16 Class A, 3 Class B.

**HIS433 is the result.** pKa 6.22, relative SASA 0.648, identical in mouse, and inside the cetuximab contact footprint. That pKa is close to ideal for a 6.5 / 7.4 switch — roughly 37% protonated at pH 6.5 against 6% at 7.4, a ~6× occupancy swing from the residue alone before any pKa shift contributed by a binder carboxylate. One residue servicing affinity, cross-reactivity and pH selectivity at once.

**HIS358** is a second Class B anchor at pKa 6.26, similarly conserved, but outside the Fab footprint.

**HIS383 is excluded** — it is arginine in mouse, which loses titratability entirely rather than merely changing identity. The kind of thing strict conservation filtering is for.

Four anchor patches meet the ≥3-pair spec at a 25 Å span. Widening to 32 Å, which is nearer a VHH paratope's reach than a miniprotein's, brings HIS433 and GLU412 into one region together with both Class B histidines.

Caveat worth stating plainly: these are PROPKA3 predictions on one crystal conformation. They are a tier-one screen, not measurements. A static structure cannot see a salt bridge that dissociates in 5 ns of MD.

## Not done yet

- BindCraft / RFantibody launcher with GCS checkpointing on Spot A100s
- Rosetta ddG at explicit HIS / HIS_D / HIS_E / HIS_P states
- MM/GBSA salt-bridge occupancy, protonated vs neutral
- Scoring funnel (ipSAE ≥ 0.61, interface energetics, liability scrub)

## Run it

Needs **Python 3.12 or 3.13**. On 3.14, `propka` 3.5.1 crashes in its parameter parser — PEP 649 changed `__annotations__` resolution and propka reads it directly, so it fails before touching your structure.

```bash
python3.12 -m venv .venv && source .venv/bin/activate
python -m pip install biopython propka

curl -O https://files.rcsb.org/download/6ARU.pdb
curl -L -o mouse_egfr.fasta "https://rest.uniprot.org/uniprotkb/Q01279.fasta"

python epitope_prep.py 6ARU.pdb \
    --target-chain A --partner-chains B C \
    --mouse-fasta mouse_egfr.fasta \
    --outdir out/
```

In 6ARU, chain A is EGFR and B/C are the cetuximab Fab. Sanity check on any other structure: `grep ^ATOM file.pdb | cut -c22 | sort -u`.

Outputs `residues.csv`, `anchors.csv`, `patches.json`, `hotspots.json`, `mpnn_bias.json`, `target_clean.pdb`.

6ARU uses mature numbering (residue 1 = Leu25), so precursor = file + 24. If your structure differs, pass `--numbering precursor`.

The annotated anchor map is published separately as a dataset: **[hf.co/datasets/<user>/egfr-domain3-anchor-map](https://huggingface.co/datasets)**

## License

Code MIT. Competition results ODC-BY, per Adaptyv's terms.