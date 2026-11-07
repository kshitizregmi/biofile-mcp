# Examples

Real output from the server, captured by `run_examples.py`. Nothing here is
written by hand: the script builds the test files from `make_test_data.py`,
starts `server.py` over stdio, and calls each tool with an MCP client, the
same way Claude or another assistant would.

The test files are tiny and made up, with values chosen in advance (for
example, a quality drop placed at read position 171), so every number below
can be checked by eye. For what each tool returns on real sequencing data,
and what the numbers mean, see the [format guides](../docs/README.md).

| File | Contents |
|---|---|
| [`mcp_session.md`](mcp_session.md) | What the AI client sees: the list of tools the server offers, then every request and the response it received |
| [`outputs/`](outputs/) | The JSON result of each successful call, one file per call |
| [`server.log`](server.log) | What the server writes on its side for the same session |

## Results and what they show

| Call | Result | What it shows |
|---|---|---|
| `vcf_stats(test.vcf)` | 5 records: 3 SNPs, 2 indels, Ts/Tv 2.0; samples TUMOR and NORMAL | Variant counts by type and per-sample genotypes. NORMAL has one missing genotype. [VCF guide](../docs/vcf.md) |
| `bam_stats(test.bam)` | 25 reads, 23 mapped (92%), 3 duplicates | The whole-file read counts. [BAM guide](../docs/bam.md) |
| `bam_stats(test.bam, region="chr1:1001-1100")` | 20 reads, 50% of bases covered, mean depth 10 | Coverage of one region. Duplicates are left out, so 20 of the 23 mapped reads count |
| `bam_qc(test.bam)` | Single-end, 0% error rate, every read at MAPQ 60, coverage 23 over 50 bases | The full alignment report. Duplicates count toward coverage here, as in `samtools stats` |
| `bam_depth_profile(test.bam, region="chr1:1001-1100", max_points=10)` | 20x over 1001–1050 and 0x after, in 10 bins of 10 bp | Where in the region the reads are; half of it has none |
| `fastq_stats(sample_R1.fastq)` | 1,000 reads of 240–251 bp, 71.6% Q30, 20 N bases, 40 lowercase bases | Run-level read quality. Lowercase bases mean the reads were edited. [FASTQ guide](../docs/fastq.md) |
| `fastq_stats(sample_R1.fastq, max_reads=100)` | First 100 reads only, marked `"sampled": true` | How to get a quick estimate on a large file |
| `fastq_quality_profile(sample_R1.fastq, max_points=26)` | Quality in 26 bins of 10 bp; the median falls below Q30 in the bin starting at 171 | Where to trim: after about 170 bases |
| `h5ad_summary(small.h5ad, column="cell_type")` | 200 cells × 100 genes, raw counts; 100 T, 60 B and 40 NK cells | What the file contains and how many cells are in each group. [h5ad guide](../docs/h5ad.md) |
| `h5ad_column_stats(small.h5ad, column="cell_type")` | Categorical: 100 / 60 / 40 cells | One `obs` column described on its own |

## Error handling

The last three calls in the session are wrong on purpose, to show that the
server refuses them safely and says what to fix, so the assistant can correct
its request:

| Request | Message returned | Why |
|---|---|---|
| `vcf_stats(path="/etc/passwd")` | `/etc/passwd is outside the allowed folder test_data` | The tools can only read inside `BIOFILE_ALLOWED_DIR` |
| `bam_stats(path="test.bam", region="1:1001-1100")` | `Contig '1' is not in the BAM header. Contigs in this file: ['chr1']` | The BAM names its chromosome `chr1`, not `1` |
| `h5ad_summary(path="small.h5ad", column="louvain")` | `Column 'louvain' not found. Available: ['cell_type']` | The file has no `louvain` column |

## Reading the server log

Every call writes a `call` line with its arguments, then one line with the
outcome:

```
INFO     call bam_stats path=test.bam args={'region': None}     ← request received
INFO     done bam_stats in 0.00s                                ← succeeded, with run time
WARNING  refused vcf_stats: /etc/passwd is outside ...          ← blocked by the folder check
WARNING  rejected bam_stats after 0.00s: Contig '1' is not ...  ← invalid argument
ERROR    failed <tool> after ...s  + traceback                  ← an unexpected error (none here)
```

Invalid arguments are a normal part of an assistant's work, so they get a
one-line warning. A full traceback is written only for unexpected errors,
where it is needed to find the bug. The log goes to stderr and to the file
named by `BIOFILE_LOG_FILE`, never to stdout, which carries the MCP messages.

## Regenerating

```bash
python examples/run_examples.py
```

Paths in the saved files are shortened (`test_data/...`, `./...`,
`<python-env>/...`) so they do not depend on the machine that ran the script.
Timestamps and durations will differ between runs; the summary values will not.
