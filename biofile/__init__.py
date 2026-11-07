"""Read-only summaries of common genomics file formats.

Each function takes a file path and returns a small, JSON-serializable dict,
so the same code can back a command-line tool or an MCP server.
"""
from biofile.bam import bam_stats
from biofile.fastq import fastq_stats
from biofile.h5ad import h5ad_summary
from biofile.vcf import vcf_stats

__all__ = ["vcf_stats", "bam_stats", "h5ad_summary", "fastq_stats"]
