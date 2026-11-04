# GEMAgent

GEMAgent is an evidence-guided curation framework for genome-scale metabolic
models (GEMs). It uses the official [gapseq](https://github.com/jotech/gapseq)
pipeline to reconstruct a draft model and then applies specialized review
agents to biomass precursor requirements, gapfill additions, gene-protein-
reaction (GPR) associations, and uptake/secretion capabilities.

The package provides a command-line interface and a Python API. It is designed
for reproducible model construction and for research on evidence-grounded
metabolic model curation.

## What GEMAgent does

GEMAgent separates reconstruction from review. The standard full workflow is:

```text
genome
  -> official gapseq annotation, draft construction, medium prediction, and fill
  -> biomass reaction review
  -> gapfill review
  -> GPR review
  -> literature-supported uptake/secretion expansion
  -> uptake/secretion pruning with model-consistency guards
  -> reviewed SBML, RDS, decision records, and provenance manifests
```

The review stages are modular. Each stage can be run independently, records its
evidence and decisions, and checks the resulting model before accepting an edit.

## Design principles

- **Official gapseq remains the reconstruction engine.** GEMAgent does not
  vendor or patch gapseq.
- **Review is bounded.** Agents can modify only the artifacts assigned to their
  review scope.
- **Evidence is explicit.** Decisions reference model evidence, annotations, or
  literature statements instead of hidden benchmark labels.
- **Evaluation labels stay outside runtime prompts.** Accuracy and truth metrics
  are computed only after model edits have been recorded.
- **Every run is auditable.** Workspaces contain decisions, model lineage,
  input hashes, executable paths, and a run manifest.

## Requirements

- Python 3.10 or newer
- `gapseq` available on `PATH` or configured with `GAPSEQ_BIN`
- R and the R packages required by gapseq, available through `Rscript`
- A GEMAgent-curated evidence bundle containing `biomass/`, `media/`,
  `mimedb/`, `metabolite/`, and `reaction/`
- An OpenAI-compatible chat-completions endpoint for LLM-backed review stages

For the complete Python and R environment used by the project, use the provided
[`Dockerfile`](Dockerfile).

## Installation

### Local Python environment

```bash
git clone <repository-url>
cd GEMAgent

python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate

python -m pip install --upgrade pip
python -m pip install -e .
gemagent --help
```

Install development dependencies and run the test suite:

```bash
python -m pip install -e ".[dev]"
pytest
```

### Docker

```bash
docker build -t gemagent .
docker run --rm gemagent --help
```

For real model construction, mount the evidence bundle and any external data
directories into the container, then provide the same environment variables
listed below.

## Configuration

Copy the example configuration and replace every placeholder:

```bash
cp .env.example .env
```

| Variable | Required | Purpose |
|---|---:|---|
| `LLM_API_KEY` | LLM reviews | API key for the OpenAI-compatible endpoint |
| `LLM_BASE_URL` | LLM reviews | Chat-completions base URL |
| `LLM_MODEL` | LLM reviews | Model or endpoint identifier |
| `LLM_MODEL_THINK` | Optional | Model identifier used for reasoning-enabled GPR review |
| `GAPSEQ_BIN` | Recommended | Absolute path to the official `gapseq` executable |
| `GEMAGENT_RSCRIPT` | Recommended | Absolute path to `Rscript` |
| `GEMAGENT_DB_DIR` | Reviews | Root of the GEMAgent evidence bundle |
| `GEMAGENT_EVAL_ROOT` | Optional | Evaluation data and legacy exchange-map discovery |
| `GEMAGENT_US_REVIEW_DIR` | Optional | Directory containing `exchange_map.tsv` |
| `GEMAGENT_SEED_METABOLITES` | Optional | SEED metabolite table used by the charge repair helper |

Additional performance and diagnostic variables are documented in
[`docs/CONFIGURATION.md`](docs/CONFIGURATION.md).

Environment variables passed directly to the shell take precedence over values
in `.env`.

## Evidence bundle

The public package does not redistribute the GEMAgent-curated evidence bundle.
Set `GEMAGENT_DB_DIR` to a directory with this layout:

```text
gemagent-data/
├── biomass/
├── media/
├── mimedb/
├── metabolite/
└── reaction/
```

See [`docs/DATA.md`](docs/DATA.md) for the expected inputs and exchange-map
resolution rules.

## Quick start

### Build with the full review workflow

```bash
gemagent build /path/to/genome.fna.gz \
  --out runs/my-strain \
  --recipe full \
  --org-name "Example organism" \
  --kofam-tsv /path/to/kofam.tsv \
  --claims-dir /path/to/claims \
  --exchange-map /path/to/us_review/data/exchange_map.tsv
```

The `full` recipe runs:

```text
official find
-> find-transport
-> draft
-> KOfam/medium
-> biomass review
-> official fill
-> gapfill review
-> GPR review
-> uptake/secretion expansion
-> uptake/secretion pruning
```

### Build the unreviewed gapseq baseline

```bash
gemagent build /path/to/genome.fna.gz \
  --out runs/my-strain-base \
  --recipe base-only
```

The `base-only` recipe uses only the official gapseq path and no review agents.
It is the appropriate baseline when comparing reviewed and unreviewed models.

### Inspect execution without running tools

```bash
gemagent build /path/to/genome.fna.gz \
  --out runs/my-strain \
  --recipe base-only \
  --dry-run
```

## Review skills

List the registered review agents:

```bash
gemagent skill list
```

Run an individual review:

```bash
gemagent skill run biomass \
  --workspace runs/my-strain \
  --genome-id my-strain \
  --kofam-tsv /path/to/kofam.tsv

gemagent skill run gapfill-review \
  --workspace runs/my-strain \
  --genome-id my-strain

gemagent skill run gpr \
  --workspace runs/my-strain \
  --genome-id my-strain \
  --kofam-tsv /path/to/kofam.tsv

gemagent skill run us-review-expand \
  --workspace runs/my-strain \
  --genome-id my-strain \
  --claims-dir /path/to/claims

gemagent skill run us-review-prune \
  --workspace runs/my-strain \
  --genome-id my-strain
```

| Skill | Input scope | Output |
|---|---|---|
| `biomass` | Draft model, medium, KOfam annotations | Reviewed biomass precursor requirements |
| `gapfill-review` | Draft, filled model, medium, candidate tables | Exclusions followed by an official gapseq refill |
| `gpr` | Filled model, KOfam/EC evidence, reaction context | REMOVE-only branch review and invariant checks |
| `us-review-expand` | Positive literature claims and exchange map | `ADD`, `CREATE`, or `NO_ACTION` |
| `us-review-prune` | Expanded model plus model-lineage evidence | `KEEP`, `BLOCK`, or `UNCERTAIN` |

`us-review` is retained as a compatibility alias for `us-review-expand`; it does
not perform pruning.

### Optional Deep Research readers

The standalone phenotype-reader scripts in `agent_us_review/` use helper modules
for text windowing and evidence validation. Point GEMAgent at the helper package
before running those scripts:

```bash
export GEMAGENT_DEEP_RESEARCH_DIR=/path/to/agent_deepresearch
```

The core `gemagent build` workflow consumes claim tables directly and does not
require these helper modules.

## Python API

```python
import gemagent

gemagent.build(
    "/path/to/genome.fna.gz",
    "runs/my-strain",
    recipe="full",
    org_name="Example organism",
    kofam_tsv="/path/to/kofam.tsv",
    claims_dir="/path/to/claims",
    exchange_map="/path/to/us_review/data/exchange_map.tsv",
)

gemagent.run_skill(
    "us-review-prune",
    "runs/my-strain",
    genome_id="my-strain",
)
```

## Evaluation

Evaluation code consumes model predictions only after decisions have been
written. For an uptake/secretion table:

```bash
gemagent eval-us \
  --truth evaluation/truth.tsv \
  --calls runs/my-strain/us_review_prune/decisions.tsv \
  --derived-from-decisions \
  --out evaluation/results
```

The command reports the confusion matrix, accuracy, sensitivity, specificity,
precision, recall, F1, F2, false-positive rate, balanced accuracy, and MCC.

## Outputs

A completed workspace contains model artifacts and provenance records:

```text
runs/my-strain/
├── _input/
├── logs/
├── biomass_check/
├── gapfill_review/
├── gpr_check/
├── us_review/
├── us_review_prune/
├── my-strain-draft.RDS
├── my-strain-gapfilled.xml
├── my-strain-gapfilled-gprchecked.xml
├── my-strain-gapfilled-gprchecked-usreviewed.xml
├── my-strain-gapfilled-gprchecked-usreviewed-pruned.xml
└── run_manifest.json
```

Exact intermediate filenames depend on which review stages are enabled. See
[`docs/OUTPUTS.md`](docs/OUTPUTS.md) for the decision and provenance contracts.

## Repository layout

```text
.
├── docs/                 Architecture, data, configuration, and validation notes
├── src/gemagent/         Python package and model-review implementations
├── tests/                Unit and protocol tests
├── Dockerfile            Reproducible gapseq + Python+R environment
├── pyproject.toml        Package metadata and dependencies
└── .env.example          Public configuration template
```

## Documentation

- [Architecture](docs/ARCHITECTURE.md)
- [Configuration](docs/CONFIGURATION.md)
- [Evidence data](docs/DATA.md)
- [GPR review protocol](docs/OR2_GPR.md)
- [Output contract](docs/OUTPUTS.md)
- [Validation summary](docs/VALIDATION.md)

## Reproducibility

- The official gapseq executable and R environment are recorded in
  `run_manifest.json`.
- Input and output files are hashed where the workflow requires provenance.
- Decision tables retain the evidence and model context used by each review.
- Outputs are written to stable workspace paths so downstream evaluation does
  not depend on local absolute paths.
- The test suite can run without LLM credentials where the reviewed code path is
  not exercised.

## Reporting issues and contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for development setup and pull-request
guidance. Report security-sensitive issues according to
[`SECURITY.md`](SECURITY.md).

## Citation

If you use GEMAgent in research, cite the accompanying manuscript. Citation
metadata is provided in [`CITATION.cff`](CITATION.cff).

## License

GEMAgent is released under the [MIT License](LICENSE).
