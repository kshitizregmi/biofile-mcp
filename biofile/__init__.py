"""Read-only summaries of common genomics file formats.

Each function takes a file path and returns a small, JSON-serializable dict,
so the same code can back a command-line tool or an MCP server.
"""
from biofile.bam import bam_depth_profile, bam_qc, bam_stats
from biofile.fastq import fastq_quality_profile, fastq_stats
from biofile.formats import detect_format
from biofile.h5ad import h5ad_column_stats, h5ad_embedding, h5ad_summary
from biofile.vcf import vcf_stats

__all__ = [
    "vcf_stats",
    "bam_stats",
    "bam_depth_profile",
    "bam_qc",
    "h5ad_summary",
    "h5ad_column_stats",
    "h5ad_embedding",
    "fastq_stats",
    "fastq_quality_profile",
    "detect_format",
]
