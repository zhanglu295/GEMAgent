# Evidence data

The package does not vendor the GEMAgent-curated evidence bundle. Set:

```bash
export GEMAGENT_DB_DIR=/path/to/gemagent-data
```

The directory must contain the gapseq-compatible evidence files used by the
reviews:

```text
biomass/       biomass template JSON files
media/         medium prediction rules / curated media
mimedb/        organism context profiles
metabolite/    metabolite tables
reaction/      seed reaction tables
```

`us-review-expand` also needs an exchange-name map. It can be supplied with
`--exchange-map`, through `GEMAGENT_US_REVIEW_DIR`, or under
`$GEMAGENT_EVAL_ROOT/us_review/data/exchange_map.tsv`. Claims may instead carry
an explicit `exchange_id` / `gapseq_exchange` column.

`us-review-prune` accepts an explicit GTDB-Tk summary (`--gtdb-summary`) and a
MIMEdb TSV (`--mimedb-tsv`). Supported GTDB files include a normal
`user_genome`/`gtdb_taxonomy` table, a two-column genome/taxonomy table, or a
single lineage line. A single-row MIMEdb profile is accepted directly. Missing
fields become `unknown` and are neutral in triage.

The taxonomy tables are optional for workflows where `--org-name` and
`--taxonomy` are supplied explicitly. Evaluation truth is never placed in the
prune feature table or prompt.
