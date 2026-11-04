# Architecture

The base reconstruction is delegated to the installed official `gapseq` CLI:

`find -> find-transport -> draft -> medium -> fill`

The package owns only post-draft review orchestration and review skills:

- `biomass`: KOFam + medium evidence -> biomass decision -> adjusted draft.
- `gapfill-review`: review reactions added by official fill -> exclusion list ->
  filtered candidate RDS -> official `gapseq fill` -> evidence diff.
- `gpr`: pure-OR branch review with KOfam/EC/reaction context -> LLM
  KEEP/REMOVE/UNCERTAIN for every selected branch -> REMOVE-only application,
  never empty a GPR -> reviewed SBML/RDS and model-invariant QC.
- `us-review-expand`: positive literature claims -> `ADD` / `CREATE` / `NO_ACTION`.
- `us-review-prune`: existing calls + four-generation model evidence -> `KEEP` / `BLOCK` / `UNCERTAIN`.

Canonical implementations: `impl/agent/agent_us_expand.py` (expand) and `impl/agent/agent_us_prune.py` (prune). The retired single-pass `agent_us_review/gapseq_us_review.py` was removed in 2f4b81b; older checkouts can still recover it from git history.

The gapfill filter and evidence exporter are small data-format helpers. They do
not contain gapfilling logic and never modify the installed gapseq source.

## Scope

This repository contains the reconstruction pipeline and the review skills it
runs. Dataset-specific experiment drivers, prompt-replay utilities, and large
evaluation outputs are distributed separately from the core package.

## Model namespaces

Reconstruction families disagree on species ids, compartment encodings and
exchange names. `gemagent.namespace` owns those differences (`seed` for
gapseq/ModelSEED, `bigg` for CarveMe, `agora` for AGORA2) and can detect the
family from an SBML file. Every review skill receives an adapter instead of
hard-coding prefixes; `gemagent namespace inspect MODEL` reports what was
detected and which exchanges cannot be resolved.

## Truth boundary

Truth panels, BacDive observations, and evaluation labels are excluded from
runtime prompts and feature tables. Expansion uses the claim direction supplied
in the claim table. Pruning uses model lineage, gapfill provenance, GPR/KOfam,
pathway completeness, network consistency, explicit literature statements,
GTDB-Tk lineage, and MIMEdb context.

Positive claims are a protection set. Explicit negative claims are evidence;
absence of literature is neutral. A missing GTDB or MIMEdb field is `unknown`
and remains neutral.

## Prune triage and guards

`PROTECTED` calls have a positive claim or strong GPR/pathway support.
`BLOCK_ELIGIBLE` calls have an explicit negative claim or at least two
independent risk categories. `UNCERTAIN` calls are retained by default.

Only `BLOCK_ELIGIBLE` rows are sent to the LLM. Every LLM `BLOCK` must cite
non-empty evidence references. Blocks are applied one at a time and checked
against two guards:

1. growth must remain at least 90% of the pre-prune baseline;
2. every initially capable positive-claim call must remain capable.

Failed blocks are rolled back individually; successful blocks are retained.
