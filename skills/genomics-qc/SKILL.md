---
name: genomics-qc
description: Quality control and summaries of sequencing files with the biofile tools. Use when the user asks to check, summarize, QC, compare or explain a VCF, BCF, BAM, FASTQ or h5ad (AnnData) file, asks whether a run, library, alignment, call set or single-cell dataset looks good, or asks about Ts/Tv, Q30, mapping rate, duplicates, insert size, coverage, genotypes or cell QC metrics.
---

# Genomics file QC with biofile

The biofile tools compute summaries locally with bcftools, samtools, pyfastx
and anndata. They are read-only and can only open files inside the user's
data folder: `${user_config.data_dir}`. Paths may be given relative to that
folder.

## 1. Find the file

If the user names a file, use it. Otherwise list the data folder with your
normal file tools and pick by extension:

| Extension | Format | Notes |
| --- | --- | --- |
| `.vcf`, `.vcf.gz`, `.bcf` | Variant calls | `.tbi`/`.csi` index files are not inputs |
| `.bam` | Alignments | `bam_depth_profile` and region queries need a `.bai` or `.csi` next to it |
| `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz` | Reads | Paired files usually end in `_R1`/`_R2` |
| `.h5ad` | AnnData (single cell) | Opened in backed mode, so large files are fine |

If a file sits outside the data folder, tell the user and suggest moving it
or changing the data folder with `/config`. Do not try to work around the
restriction.

## 2. Pick the tool for the question

Start with the overview tool, then go deeper only where the question needs it.

| Question | Tool | Key arguments |
| --- | --- | --- |
| How many variants, what kinds, Ts/Tv, per-sample genotypes | `vcf_stats` | `path` |
| Read counts, mapping rate, duplicates | `bam_stats` | `path`, optional `region` |
| Full alignment QC: MAPQ, insert size, error rate, coverage, per contig | `bam_qc` | `path`, optional `region` |
| Depth along one region | `bam_depth_profile` | `path`, `region` (at most 10 Mb), `min_depth` |
| Read count, length, Q20/Q30, GC, N content | `fastq_stats` | `path`, optional `max_reads` |
| Per-position quality, base composition, per-read quality and GC | `fastq_quality_profile` | `path`, `max_reads` (default 100,000) |
| Cells, genes, metadata columns, embeddings, cell-type counts | `h5ad_summary` | `path`, optional `column` |
| Distribution of one metadata column (e.g. `pct_counts_mt`) | `h5ad_column_stats` | `path`, `column` |

Large files:

- `fastq_stats` reads every record unless `max_reads` is set. For files over
  a few GB, run it with `max_reads=1000000` first and say the numbers come
  from a sample. The result's `sampled` field tells you whether it was.
- `bam_qc` on a whole genome reads the entire file. For a quick look, pass a
  `region` such as a single chromosome, then offer the full run.
- Regions use samtools syntax: `chr1`, `chr1:1000000-2000000`. Contig names
  must match the BAM header exactly (`chr1` and `1` are different). The tool
  error lists the valid names if one is wrong.

## 3. Interpret the numbers

Judge each metric against what is normal for that assay, not against a single
universal cutoff. Read `${CLAUDE_PLUGIN_ROOT}/skills/genomics-qc/references/qc-thresholds.md`
for typical ranges and the usual causes when a value is off. Ask the user for
the assay (WGS, exome, RNA-seq, amplicon, single-cell) and organism when the
interpretation depends on it and the file does not make it clear.

Rules that hold across formats:

- A value outside the typical range is a flag to explain, not a verdict. Say
  what would cause it and how the user could confirm.
- Never report a metric the tools did not return. If a question needs data
  the tools do not provide (for example adapter content or variant
  annotation), say so and name the tool that would answer it (FastQC or
  fastp, VEP or snpEff, and so on).
- Distinguish sampled from full-file results.
- Tiny test files give unstable ratios. Mention it when a file has fewer than
  about 1,000 reads or 100 variants.

## 4. Report

Lead with a one-line verdict, then the evidence:

1. **Verdict**: pass, pass with notes, or needs attention, in one sentence.
2. **Key metrics**: a short table of the numbers that drove the verdict,
   each with its typical range.
3. **Flags**: each issue with the metric, the likely cause and a concrete
   next step.
4. **What was not checked**, when it matters for the user's question.

Round percentages to one decimal place and large counts with thousands
separators. Keep the report short; offer the full tool output on request.

When comparing files (R1 vs R2, tumour vs normal, before vs after
filtering), run the same tool on each and present the metrics side by side.
