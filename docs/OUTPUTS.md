# Output contract

Each build writes a self-contained workspace. The manifest is the primary
entry point for downstream automation.

## Manifest

`run_manifest.json` records:

- workflow status and recipe;
- input genome and inferred genome identifier;
- resolved gapseq, Rscript, and evidence paths;
- errors raised by a failed stage;
- generated top-level output files.

Individual review stages write their own manifests and decision records.

## Review decisions

Decision tables use stage-specific JSON and TSV representations. They preserve
the candidate identifier, review outcome, applied state, evidence summary, and
model validation result. Evaluation truth is not included in runtime feature
tables or prompts.

The uptake/secretion pruning stage writes:

```text
us_review_prune/
├── features.tsv
├── decisions.json
├── decisions.tsv
└── run_manifest.json
```

## Model lineage

The full workflow can produce the following model sequence:

```text
<genome>-draft.RDS
<genome>-gapfilled.RDS
<genome>-gapfilled.xml
<genome>-gapfilled-gprchecked.xml
<genome>-gapfilled-gprchecked-usreviewed.xml
<genome>-gapfilled-gprchecked-usreviewed-pruned.xml
<genome>-gapfilled-gprchecked-usreviewed-pruned.RDS
```

A stage may be omitted by its corresponding `--no-*` flag. Downstream stages
use the newest compatible model available in the workspace.

## Review logs and audits

Review stages write logs and machine-readable audit records under their local
subdirectories. Use these files, not console output, as the source of truth for
reproducibility.
