# biofile-mcp

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

## Project layout

```
biofile-mcp/
├── biofile/            # plain Python functions: all the logic lives here
│   ├── vcf.py
│   ├── bam.py
│   ├── h5ad.py
│   └── fastq.py
├── server.py           # MCP layer: exposes the functions as tools
├── main.py             # run the functions from the terminal, without MCP
├── make_test_data.py   # tiny test files with known answers
├── environment.yml     # conda environment
└── setup.sh            # one-step setup: env + tests + Claude Code
```

## Quick setup with conda (recommended)

```bash
cd biofile-mcp
bash setup.sh                       # Claude may read ./test_data
bash setup.sh ~/path/to/your/data   # or choose your data folder
```

This creates the `biofile` conda env, makes test files, checks all 4 tools and
registers the server with Claude Code. Claude Code then starts the server with
the env's own Python, so you never need to `conda activate` it yourself.
Rerun it any time: it updates the env and re-registers the server.

## Manual setup (Linux or macOS; on Windows use WSL, because pysam doesn't support Windows)

```bash
cd biofile-mcp
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Requires Python 3.10+.

## 1. Test the functions (no MCP)

```bash
python make_test_data.py
python main.py test_data/test.vcf
python main.py test_data/test.bam --region chr1:1001-1100
python main.py test_data/small.h5ad --column cell_type
python main.py test_data/sample_R1.fastq --max-reads 100
```

## 2. Test the MCP server by hand (no AI)

```bash
BIOFILE_ALLOWED_DIR=$PWD/test_data mcp dev server.py
```
This opens the MCP Inspector in your browser, where you can call each tool yourself.

## 3. Connect to Claude Code

```bash
claude mcp add biofile \
  -e BIOFILE_ALLOWED_DIR=/full/path/to/your/data \
  -- /full/path/to/biofile-mcp/.venv/bin/python /full/path/to/biofile-mcp/server.py
```

Check it's connected with `claude mcp list`, then ask Claude Code things like:
- "What's the Ts/Tv ratio of test.vcf?"
- "How many cell types are in small.h5ad and how many cells each?"
- "Is the quality of sample_R1.fastq good enough, or should I trim?"

## Safety

- Read-only: no tool writes, moves or deletes files.
- Only files inside `BIOFILE_ALLOWED_DIR` can be read; other paths are refused.
- Errors are returned to the model as clear messages so it can correct its call.
