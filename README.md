# biofile-mcp

[![tests](https://github.com/kshitizregmi/biofile-mcp/actions/workflows/tests.yml/badge.svg)](https://github.com/kshitizregmi/biofile-mcp/actions/workflows/tests.yml)

An MCP server that lets AI assistants (Claude Code, Claude Desktop, Cursor, local models)
read and summarize genomics files: VCF, BAM, FASTQ and single-cell h5ad.
It returns small JSON summaries, so large files never go into the model's context.
Files stay on your machine; the tools are read-only.

## Tools

| Tool | Input | Returns |
|---|---|---|
| `vcf_stats` | VCF / VCF.gz | samples, records, SNPs, indels, multiallelic sites, Ts/Tv (via bcftools stats) |
| `bam_stats` | BAM (+ optional region) | total / mapped / duplicate reads; region read count, % covered, mean depth (via samtools) |
| `h5ad_summary` | h5ad (+ optional column) | cells, genes, metadata, embeddings, raw-counts check, cells per category |
| `fastq_stats` | FASTQ / FASTQ.gz | reads, lengths, Q20/Q30, GC, N, lowercase and per-base counts (via pyfastx) |

See [`examples/`](examples/) for real output: the JSON each tool returns, a full
MCP session transcript and the matching server log.

## Project layout

```
biofile-mcp/
├── biofile/              # plain Python functions: all the logic lives here
│   ├── vcf.py
│   ├── bam.py
│   ├── h5ad.py
│   └── fastq.py
├── server.py             # MCP layer: exposes the functions as tools, with logging
├── main.py               # run the functions from the terminal, without MCP
├── make_test_data.py     # tiny test files with known answers
├── tests/                # pytest suite, including a real MCP client test
├── examples/             # captured server output: results, session transcript, log
├── environment.yml       # conda environment
├── setup.sh              # one-step setup: env + tests + Claude Code
└── .github/workflows/    # CI: runs the tests on Linux and macOS
```

## Quick setup with conda (recommended)

```bash
cd biofile-mcp
bash setup.sh                       # Claude may read ./test_data
bash setup.sh ~/path/to/your/data   # or choose your data folder
```

This creates the `biofile` conda env, makes test files, runs the test suite and
registers the server with Claude Code. Claude Code then starts the server with
the env's own Python, so you never need to `conda activate` it yourself.
Rerun it any time: it updates the env and re-registers the server.

## Manual setup (Linux or macOS; on Windows use WSL, because pysam doesn't support Windows)

```bash
cd biofile-mcp
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt -r requirements-dev.txt
```

Requires Python 3.12+.

## Run the tests

```bash
python -m pytest -v
```

The suite checks every tool against files with known answers (e.g. Ts/Tv 2.0,
92% mapped reads, 100/60/40 cells per type), error cases (empty or invalid
files, wrong contig names, unknown columns), the allowed-folder checks
(including `../` escapes), logging, and an end-to-end test where a real MCP
client starts `server.py` over stdio and calls the tools, the same way Claude does.

## Try the functions from the terminal (no MCP)

```bash
python make_test_data.py
python main.py test_data/test.vcf
python main.py test_data/test.bam --region chr1:1001-1100
python main.py test_data/small.h5ad --column cell_type
python main.py test_data/sample_R1.fastq --max-reads 100
```

## Test the MCP server by hand (no AI)

```bash
BIOFILE_ALLOWED_DIR=$PWD/test_data mcp dev server.py
```
This opens the MCP Inspector in your browser, where you can call each tool yourself.

## Connect to Claude Code

```bash
claude mcp add biofile \
  -e BIOFILE_ALLOWED_DIR=/full/path/to/your/data \
  -- /full/path/to/python /full/path/to/biofile-mcp/server.py
```

Check it's connected with `claude mcp list`, then ask Claude Code things like:
- "What's the Ts/Tv ratio of test.vcf?"
- "How many cell types are in small.h5ad and how many cells each?"
- "Is the quality of sample_R1.fastq good enough, or should I trim?"

## Configuration and logs

| Variable | Default | Meaning |
|---|---|---|
| `BIOFILE_ALLOWED_DIR` | current folder | only files inside this folder can be read |
| `BIOFILE_LOG_LEVEL` | `INFO` | `DEBUG`, `INFO`, `WARNING` or `ERROR` |
| `BIOFILE_LOG_FILE` | none | also write logs to this file |

Every tool call is logged with its arguments, duration and outcome:

```
2026-11-06 21:40:12,517 INFO biofile: call vcf_stats path=test.vcf args={}
2026-11-06 21:40:12,548 INFO biofile: done vcf_stats in 0.03s
2026-11-06 21:40:15,102 WARNING biofile: refused vcf_stats: /etc/passwd is outside the allowed folder ...
```

Logs go to stderr (and the log file, if set), never to stdout, because over stdio
stdout carries the MCP protocol messages. To keep a log file, add
`-e BIOFILE_LOG_FILE=/full/path/to/biofile.log` to the `claude mcp add` command.

## Safety

- Read-only: no tool writes, moves or deletes files.
- Only files inside `BIOFILE_ALLOWED_DIR` can be read; other paths, including `../` escapes, are refused.
- Errors are returned to the model as clear messages so it can correct its call.
