# OR2-LLM-v1 as the core GPR strategy

`gemagent-core-v2` replaces the former default GPR reviewer with the frozen
OR2-LLM-v1 protocol used by the paper evaluation.

## Protocol

- Select pure-OR rules with at least `min_or` branches (default: 2).
- Use KOfam KO definitions, gene EC groups, reaction names and organism context
  as evidence. Deterministic EC grouping does not decide the final edit.
- Ask the LLM to adjudicate every selected branch as KEEP, REMOVE or UNCERTAIN.
- Apply REMOVE only. KEEP and UNCERTAIN leave the branch unchanged.
- Never remove every gene from a rule. Do not edit reactions, metabolites,
  stoichiometry or flux bounds.
- Do not supply essentiality labels, curated GPRs or other benchmark truth.

The original EC-consensus prompt/reviewer and REMOVE/ADD entrypoints were removed
from this version. `gpr_evidence.py` remains only as an evidence-grouping helper.

## Execution contract

`step_gpr_check` exports a neutral pre-review model, builds the OR2 audit,
freezes LLM decisions, applies REMOVE with the SBML-appropriate applier, writes
SBML/RDS/sybil outputs, and requires `gpr_check/invariant.tsv` to report
`ALL_PASS=YES`. Missing candidates produce an empty review and leave the model
unchanged. Missing evidence, missing LLM configuration, a failed model export, or a
failed QC check is a hard failure rather than a silent success. Individual
review calls that still fail after retries become `UNCERTAIN`; if every
selected reaction fails, the stage exits non-zero.

LLM credentials are read from `LLM_ENV_FILE` or the repository-root `.env`.
No credentials are stored in the evaluation archive.
