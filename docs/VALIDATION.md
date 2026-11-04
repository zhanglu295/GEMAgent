# Validation record

Date: 2026-10-18

Completed in `gemagent-core`:

- editable install and CLI smoke tests passed;
- official `gapseq medium` and official `gapseq fill` adapter checks passed;
- candidate filtering and gapfill-added evidence extraction passed;
- `python -m pytest -q`: **61 tests passed**;
- no runtime decision code or prompt contains an `observed` truth label.

The test suite covers:

- expand without a BacDive/AGORA truth panel;
- direct positive claims producing both `ADD` and `CREATE`;
- legacy `us-review` alias remaining expand-only;
- prune features carrying draft, gapfilled, GPR-checked, and expanded lineage;
- absence of an `observed` column;
- pathway completeness and transporter homology fields;
- direct GTDB lineage and single-row MIMEdb parsing;
- full feature/prompt context propagation;
- independent block rollback for growth and positive-claim guards;
- a successful block being retained when other blocks fail.

R. albus 7 pilot (real LLM prune):

- feature table: 178 rows;
- triage: 88 `PROTECTED`, 90 `BLOCK_ELIGIBLE`, 0 `UNCERTAIN`;
- GTDB lineage and MIMEdb fields present;
- decisions: 178 `KEEP`, 0 `BLOCK`, 0 `UNCERTAIN`; 88 were deterministic and 90 were LLM adjudications;
- outputs include `features.tsv`, `decisions.json/tsv`, `run_manifest.json`,
  pruned SBML, and a 215 KB pruned RDS;
- base growth: 0.2754108295; the LLM returned `KEEP` for all 90 eligible calls, so no blocks were applied.

## 12-strain split-review rerun

All 12 BacDive strains were rerun through the corrected expand -> real-LLM
prune pipeline. The final evaluation contains 248 panel calls.

| method | F1 | accuracy | TP | FP | TN | FN |
|---|---:|---:|---:|---:|---:|---:|
| AGORA2 as-is | 0.694 | 0.726 | 77 | 24 | 103 | 44 |
| gapseq as-is | 0.534 | 0.556 | 63 | 52 | 75 | 58 |
| gapseq + us-review | 0.689 | 0.669 | 91 | 52 | 75 | 30 |
| GEMAgent from genome after expand | 0.659 | 0.633 | 88 | 58 | 69 | 33 |
| GEMAgent from genome after LLM prune | 0.659 | 0.633 | 88 | 58 | 69 | 33 |

Across the complete models, prune made 2,346 decisions and applied 462 BLOCK
calls after the two guards. None of those blocks changed any of the 248 panel
calls; consequently expand and final are identical on this evaluation panel.
The result gains recall over raw gapseq but over-predicts relative to AGORA2.

Implementation fixes found during the rerun:

- Cobra removes `M_`/`R_` prefixes when reading SBML; expand now looks up and
  creates model IDs in Cobra's namespace, preventing duplicate species and
  ensuring ADD/CREATE reactions are truly written back.
- FVA in prune is single-process to avoid file-descriptor exhaustion on larger
  models.
- Positive-claim guard FVA is batched per block instead of one solver call per
  protected exchange.
- Regression tests cover the exported SBML readback, ADD/CREATE presence,
  and RDS exporter availability.

R. albus 7 rerun: expand solves cellulose, xylan, lichenan, and glucomannan,
reaching TP 6 / FP 2 / TN 7 / FN 3 (accuracy 0.722, F1 0.706). The real LLM
prune applies 44 blocks outside the panel and leaves all 18 panel calls unchanged.

## Still not run

- any NJC19/NJC20 rerun. NJC remains paused pending the revised evaluation.

## Matched-medium gapseq control

The 12 gapseq models were rebuilt with each strain's predicted medium and then
passed through the same expand and real-LLM prune pipeline, claims, and 248
BacDive observations as the from-genome arm.

| method | F1 | accuracy | TP | FP | TN | FN |
|---|---:|---:|---:|---:|---:|---:|
| AGORA2 as-is | 0.694 | 0.726 | 77 | 24 | 103 | 44 |
| gapseq as-is, matched medium | 0.517 | 0.532 | 62 | 57 | 70 | 59 |
| gapseq + split review | 0.657 | 0.629 | 88 | 59 | 68 | 33 |
| GEMAgent from genome + split review | 0.659 | 0.633 | 88 | 58 | 69 | 33 |

Under the matched-medium, same-review comparison, the from-genome arm is
marginally higher: +0.002 F1 and +0.004 accuracy, with one fewer FP. Gapseq
prune applied 101 blocks and from-genome prune 462 blocks; neither changed any
of the 248 panel calls.

The 12 gapseq and from-genome draft reaction sets are identical. The earlier
0.015 from-genome deficit came from using FT.csv for every gapseq fill instead
of the strain-specific predicted medium; it was not a difference in the review
implementation. Gapseq models do not carry GEMAgent GPR/KOfam/gapfill-review
provenance, so the prune skill treats those inputs as neutral.
