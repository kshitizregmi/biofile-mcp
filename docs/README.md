# biofile guide

This guide explains each file type biofile reads: what the file is, what is
inside it, which files you need, what every tool does with it, and how to read
the result. Every example on these pages comes from real public data, and
every picture is the web interface opened on that data.

| If you have | Read | Tools |
|---|---|---|
| Sequencing reads (`.fastq`, `.fq`, `.fastq.gz`) | [FASTQ guide](fastq.md) | `fastq_stats`, `fastq_quality_profile` |
| Aligned reads (`.bam` and its `.bai`) | [BAM guide](bam.md) | `bam_stats`, `bam_qc`, `bam_depth_profile` |
| Variant calls (`.vcf`, `.vcf.gz`, `.bcf`) | [VCF guide](vcf.md) | `vcf_stats` |
| Single-cell data (`.h5ad`) | [h5ad guide](h5ad.md) | `h5ad_summary`, `h5ad_column_stats` |

## Where these files come from

A typical sequencing project moves through the four formats in this order:

```
sequencer ──▶ FASTQ ──align to a reference──▶ BAM ──call variants──▶ VCF
                 │
                 └──single-cell pipeline (Cell Ranger, STARsolo, ...) ──▶ h5ad
```

Checking each file before the next step catches problems early: a bad run is
cheaper to spot in the FASTQ than after alignment and variant calling.

## Which files each tool needs

Most tools take one file. The exceptions are BAM region queries, which also
need the BAM's index.

| Tool | Required | Also needed | Not needed |
|---|---|---|---|
| `fastq_stats`, `fastq_quality_profile` | one FASTQ file | nothing | For paired-end data, run each of `_R1` and `_R2` separately |
| `bam_stats`, `bam_qc` on the whole file | the BAM | nothing | the index |
| `bam_stats`, `bam_qc` with `region`, and `bam_depth_profile` | the BAM | its index, `sample.bam.bai` (or `sample.bai`, or `.csi`) in the same folder | the reference genome |
| `vcf_stats` | the VCF or BCF | nothing | the `.tbi`/`.csi` index, the reference |
| `h5ad_summary`, `h5ad_column_stats` | one `.h5ad` file | nothing | |

No tool writes anything next to your files. To make a missing BAM index, run
`samtools index sample.bam` (the BAM must be sorted by position). The web
interface makes one for you when you upload a BAM without its index.

## The example data

All pages use the same five public files, about 25 MB in total:

```bash
mkdir -p real_data && cd real_data
curl -L -O https://raw.githubusercontent.com/nf-core/test-datasets/ampliseq/testdata/1_S103_L001_R1_001.fastq.gz -O https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam -O https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam.bai -O https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/vcf/NA12878_GIAB.chr22.vcf.gz -O https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad
cd ..
```

| File | What it is |
|---|---|
| `1_S103_L001_R1_001.fastq.gz` | 2,500 reads from a 16S rRNA gene (bacterial) amplicon run, from the nf-core/ampliseq test data |
| `test.paired_end.sorted.bam` + `.bai` | 5,644 paired-end human reads aligned to a 40 kb piece of chromosome 22, from the nf-core test data |
| `NA12878_GIAB.chr22.vcf.gz` | 81 variants on chromosome 22 from the Genome in a Bottle benchmark for sample NA12878 |
| `pbmc3k.h5ad` | 2,638 human blood cells (PBMC 3k), already clustered and labelled with cell types |

## Three ways to run a tool

Every page shows the same three options for each tool. All three run the same
code and return the same JSON.

**1. Ask an AI assistant** that has biofile connected (see the main
[README](../README.md#connect-to-an-ai-assistant)), in plain words:

> What's the Ts/Tv ratio of NA12878_GIAB.chr22.vcf.gz?

**2. Call the tool from the terminal** with the MCP Inspector (needs Node.js).
`BIOFILE_ALLOWED_DIR` is the folder the tool may read from, and `path` is the
file inside it:

```bash
npx @modelcontextprotocol/inspector --cli -e BIOFILE_ALLOWED_DIR=$PWD/real_data python server.py --method tools/call --tool-name vcf_stats --tool-arg path=NA12878_GIAB.chr22.vcf.gz
```

**3. Open the file in the web interface** (`streamlit run app.py`), choose
**Path on this computer** in the sidebar, paste the full path and press
**Summarize**. The pictures on each page show what you get.

## Words used on these pages

| Word | Meaning |
|---|---|
| Read | One short piece of DNA or RNA sequence produced by the sequencer, usually 50–300 bases |
| Phred quality (Q) | The sequencer's confidence in a base. Q20 means a 1 in 100 chance the base is wrong, Q30 means 1 in 1,000, Q40 means 1 in 10,000 |
| Reference genome | The standard genome sequence that reads are compared against, such as GRCh38 for human |
| Contig | One sequence in the reference: a chromosome (`chr22`) or a smaller piece |
| Depth / coverage | How many reads cover a base. "30x" means 30 reads on average |
| Paired-end | Each DNA fragment is read from both ends, giving two reads (`R1` and `R2`) that belong together |
| Region | A stretch of the reference written `contig:start-end`, such as `chr22:1500-5000`. Positions start at 1 |

## Rebuilding the pictures

The pictures in `docs/screenshots/real/` are made by
[`make_screenshots.py`](make_screenshots.py), which opens each example file in
the web interface and saves one picture per tab. Run it after changing the
charts so the guides stay accurate.
