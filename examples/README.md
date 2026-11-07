# Examples

Real output from the server, captured by `run_examples.py`. Nothing here is
written by hand: the script builds the test files from `make_test_data.py`,
starts `server.py` over stdio, and calls each tool with an MCP client.

| File | Contents |
|---|---|
| [`mcp_session.md`](mcp_session.md) | Full transcript: the tools the server offers, then every request and the response the client received |
| [`outputs/`](outputs/) | The JSON result of each successful call |
| [`server.log`](server.log) | The server's log for the same session |

## Results at a glance

| Call | Key values |
|---|---|
| `vcf_stats(test.vcf)` | 5 records: 3 SNPs, 2 indels, Ts/Tv 2.0 |
| `bam_stats(test.bam)` | 25 reads, 23 mapped (92%), 3 duplicates |
| `bam_stats(test.bam, region="chr1:1001-1100")` | 20 reads (duplicates excluded), 50% of bases covered, mean depth 10 |
| `h5ad_summary(small.h5ad, column="cell_type")` | 200 cells x 100 genes, raw counts; 100 T, 60 B, 40 NK cells |
| `fastq_stats(sample_R1.fastq)` | 1,000 reads of 240-251 bp, 20 N bases, 40 lowercase bases |
| `fastq_stats(sample_R1.fastq, max_reads=100)` | First 100 reads only, marked `"sampled": true` |

## Error handling

Invalid requests return a message that says what went wrong, so the client
can correct the call:

| Request | Message returned |
|---|---|
| A path outside the allowed folder | `/etc/passwd is outside the allowed folder test_data` |
| A region with the wrong contig name | `Contig '1' is not in the BAM header. Contigs in this file: ['chr1']` |
| An `obs` column that does not exist | `Column 'louvain' not found. Available: ['cell_type']` |

The log records every call with its duration. Refused paths are logged as
warnings and failures with a full traceback.

## Regenerating

```bash
python examples/run_examples.py
```

Paths in the saved files are shortened (`test_data/...`, `./...`,
`<python-env>/...`) so they do not depend on the machine that ran the script.
Timestamps and durations will differ between runs; the summary values will not.
