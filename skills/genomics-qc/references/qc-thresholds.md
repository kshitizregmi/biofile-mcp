# Typical QC ranges

These are rules of thumb for short-read Illumina data unless stated
otherwise. Instruments, kits, read lengths and pipelines shift them, so use
them to decide what to question, not to pass or fail a sample on their own.
Field names in `code` are the keys the biofile tools return.

## FASTQ (`fastq_stats`, `fastq_quality_profile`)

| Metric | Typical | Worth a closer look | Common causes |
| --- | --- | --- | --- |
| `q30_pct` | Above 80% on recent instruments; above 75% is often the vendor spec for 2x150 | Below 70% | Over-clustered flow cell, ageing reagents, long reads on older chemistry |
| `median_below_q30_from` | Late in the read or never | Within the first half of the read | Phasing, bubbles, a bad cycle. Quality trimming helps |
| Per-position median drop at the 3' end | Gradual | Sharp drop to Q20 or below | Normal for R2 and long reads; trim if downstream tools are sensitive |
| `gc_pct` | Close to the organism: human about 41%, mouse about 42%, E. coli about 51% | More than 5 points off, or a second peak in `read_gc` | Contamination, adapter dimers, strong PCR bias |
| `n_bases` / `total_bases` | Below 0.1% | Above 1% | Failed cycles; check `base_pct` for the position |
| `base_pct` by position | Flat lines after the first ~10 bases | Uneven lines across the whole read | Expected for amplicon (16S, ITS), low-diversity and bisulfite libraries; otherwise adapters or overrepresented sequences |
| Read lengths | One length for raw reads | A spread of lengths | Already trimmed, or long-read data |

Amplicon data (16S, ITS) normally shows skewed per-position composition and
often a steeper quality drop in R2. Do not flag those as failures.

## BAM (`bam_stats`, `bam_qc`, `bam_depth_profile`)

| Metric | Typical | Worth a closer look | Common causes |
| --- | --- | --- | --- |
| `mapped_pct` | Above 95% for human WGS against a matching reference; RNA-seq aligners usually 85–95% | Below 90% for WGS | Wrong reference or build, contamination, adapter read-through, unmerged short inserts |
| `properly_paired_pct` | Above 90% | Below 80% | Wrong insert size settings, structural rearrangement, poor library |
| Duplicates / reads | WGS PCR-free below 10%; PCR WGS 10–20%; exome 10–30% | Above 30% for DNA | Low input, over-amplification, a small library sequenced deeply. High duplicates are normal in RNA-seq for highly expressed genes and in amplicon data |
| `mq0` / mapped | Below 5% for WGS | Above 10% | Repeats, segmental duplications, a reference with alt contigs and a non-alt-aware aligner |
| `error_rate` (mismatches / bases mapped) | 0.1–0.6% | Above 1% | Poor run quality, wrong reference, contamination, very divergent sample |
| `insert_size_mean`, `insert_size_sd` | Depends on the library prep; SD usually below 30% of the mean | Mean below read length, or a very wide or two-peaked distribution | Over-fragmentation, mixed libraries, adapter read-through |
| `outward_pairs` / pairs | Below 1% for standard paired-end | Above 5% | Mate-pair library, structural variants, chimeric reads |
| `pairs_on_different_contigs` / pairs | Below 2% | Above 5% | Chimeras, translocations, contamination, or reads on decoy contigs |
| `coverage.pct_bases_at_least` | WGS 30x: over 90% of bases at 10x or more | Large share of bases at zero | Low depth, GC bias, capture failure (exome), the region simply not targeted |

`bam_qc` counts duplicate reads in coverage, so its depth is slightly higher
than tools that drop duplicates first. Say so when comparing against other
reports.

For targeted panels and exomes, genome-wide coverage percentages are not
meaningful. Ask for the target BED and use `bam_depth_profile` on key genes.

## VCF (`vcf_stats`)

| Metric | Typical | Worth a closer look | Common causes |
| --- | --- | --- | --- |
| `ts_tv` | Human WGS 2.0–2.1; exome 2.8–3.3 | WGS below 1.8 or exome below 2.5 | False positive SNPs; random errors push Ts/Tv toward 0.5 |
| Het / hom-alt per sample | Human WGS about 1.5–2.0, lower in less diverse populations | Far above 2.5 | Contamination, mixed samples, reference bias |
| `singletons` per sample in a cohort | Similar across samples of similar ancestry | One sample far above the rest | Poor DNA, contamination, or genuinely different ancestry |
| `missing` per sample | Near zero for joint-called cohorts | Above 5% | Low coverage or a failed sample |
| `indel_lengths` | Peak at 1 bp, falling quickly | In coding regions, no enrichment of multiples of 3 | Indel false positives; frameshift excess in exomes |
| `qual` distribution | Most sites high | A large low-quality tail | Unfiltered call set; check whether filters were applied |
| `genotype_depth` | Most genotypes at 10x or more for WGS | Many genotypes below 10x | Low coverage; genotype calls there are less reliable |

Small subsets (single genes, a few hundred variants) give noisy ratios.

## h5ad (`h5ad_summary`, `h5ad_column_stats`)

Single-cell QC thresholds depend on tissue, chemistry and the pipeline, so
look at the distributions rather than fixed cutoffs.

| Column (common names) | What to check | Common causes when off |
| --- | --- | --- |
| `pct_counts_mt`, `percent_mito` | Usually below 5–10% for PBMC; up to 20% can be normal for heart, kidney or muscle | Dying or broken cells |
| `n_genes_by_counts`, `n_genes`, `nFeature_RNA` | A single main peak; a low tail of empty droplets, a high tail of doublets | Ambient RNA, doublets |
| `total_counts`, `n_counts`, `nCount_RNA` | Tracks `n_genes` | Very high values suggest doublets |
| Cell-type counts | Proportions plausible for the tissue | A missing expected type, or one cluster dominating, can mean a sample or annotation problem |

`X_looks_like` reports whether the main matrix holds raw counts or
normalized values. Many methods (differential expression with count models,
doublet detection, scVI) need raw counts, often stored in `layers["counts"]`
or `.raw`.
