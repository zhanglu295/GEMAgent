# Contributing

Thank you for contributing to GEMAgent.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
pytest
```

## Pull requests

- Keep changes scoped to one behavior or documentation area.
- Add or update tests for executable behavior.
- Preserve the boundary between runtime review logic and evaluation labels.
- Do not commit credentials, `.env` files, licensed solver files, large model
  outputs, or private evaluation data.
- Update documentation when CLI arguments, environment variables, or output
  contracts change.
- Run `pytest` before opening a pull request.

## Code style

- Support Python 3.10 and newer.
- Use type hints on public Python APIs.
- Keep paths relative to configuration values or the package root.
- Prefer deterministic transformations and record any stochastic behavior.
- Keep model decisions auditable through explicit decision and provenance files.

## Data and model contributions

Large evidence bundles, reconstructed models, and benchmark outputs should be
published through an archival repository and linked from the documentation.
Do not add them directly to the Git repository unless they are small,
redistributable test fixtures.

## Reporting bugs

Include the command, GEMAgent version, Python version, gapseq version, operating
system, and a minimal traceback. Remove API keys, credentials, private paths,
and unpublished data before posting logs.
