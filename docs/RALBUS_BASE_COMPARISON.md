# R. albus 7 base comparison

Date: 2026-10-18

Both arms used the same existing draft, reaction weights, reaction-to-gene table,
and pathway table. The comparison isolates the changed post-draft path.

- Official arm: `gapseq medium`, then official `gapseq fill`.
- Core arm: `gemagent build --recipe base-only --from-draft`, which calls the
  same official medium and fill through the adapter.

Results:

| Check | Official | gemagent-core |
|---|---:|---:|
| Medium SHA256 | `c2094e08...49cd9` | `c2094e08...49cd9` |
| Reactions | 1,500 | 1,500 |
| Exchange reactions | 152 | 152 |
| Reaction-set difference | none | none |
| Final XML SHA256 | `53c17796...0ab8` | `53c17796...0ab8` |

The final XML files are byte-identical. Therefore the earlier disagreement was
not caused by the standard gapseq annotation/draft step or by the official
medium/fill logic; it came from the legacy custom medium/gapfill path. The
four-review full model was intentionally not rerun in this comparison.