# Configuration

GEMAgent resolves runtime dependencies from environment variables. Shell values
take precedence over `.env`. Copy `.env.example` to `.env` and set the paths for
the target installation.

## Required for reviewed builds

| Variable | Description |
|---|---|
| `GAPSEQ_BIN` | Path to the official `gapseq` executable. If omitted, GEMAgent searches `PATH` and the active Python environment. |
| `GEMAGENT_RSCRIPT` | Path to `Rscript`. If omitted, GEMAgent searches `PATH` and the active Python environment. |
| `GEMAGENT_DB_DIR` | Root directory of the GEMAgent evidence bundle. |

## LLM configuration

| Variable | Description |
|---|---|
| `LLM_API_KEY` | API key for an OpenAI-compatible chat-completions service. |
| `LLM_BASE_URL` | Base URL for the service. |
| `LLM_MODEL` | Model or endpoint identifier used by review agents. |
| `LLM_MODEL_THINK` | Optional model identifier used by reasoning-enabled GPR review. |
| `LLM_ENV_FILE` | Optional explicit dotenv file for the standalone OR2 reviewer. |
| `GEMAGENT_LLM_USAGE_LOG` | Optional file that receives token-usage lines. |

The default code path enables reasoning for review agents. A caller can disable
it for a specific bound model with:

```python
llm.bind(thinking="disabled")
```

## Optional data paths

| Variable | Description |
|---|---|
| `GEMAGENT_EVAL_ROOT` | Root of optional evaluation data. |
| `GEMAGENT_US_DATA_DIR` | Directory containing `agora2_ref/` and `us_review/data/` for the standalone scan utility. |
| `GEMAGENT_US_REVIEW_DIR` | Directory containing an exchange-name map. |
| `GEMAGENT_SEED_METABOLITES` | Path to `seed_metabolites_edited.tsv` for formal-charge repair. |
| `GEMAGENT_GAPSEQ_DAT_DIR` | Explicit gapseq data directory when automatic discovery is insufficient. |

## Runtime tuning

| Variable | Default | Description |
|---|---:|---|
| `GEMAGENT_BIOMASS_WORKERS` | `8` | Worker count for biomass review. |
| `GEMAGENT_GPR_WORKERS` | `16` | Worker count for GPR review. |
| `GEMAGENT_REVIEW_RETRIES` | `2` | Retry count for bounded review failures. |
| `GEMAGENT_DRAFT_HIGH_EVI_BS` | `200` | High-evidence bit-score threshold passed to gapseq draft. |
| `GEMAGENT_DRAFT_MIN_BS_CORE` | `100` | Core bit-score threshold passed to gapseq draft. |
| `GEMAGENT_FIND_THREADS` | `1` | Thread count for the gapseq find stage. |

## Validation

Check the resolved executable paths without running a reconstruction:

```bash
gemagent build /path/to/genome.fna.gz --out runs/check --recipe base-only --dry-run
```

The command prints the gapseq commands and the resolved configuration without
executing them.
