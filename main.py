"""Run the biofile functions from the terminal, without MCP or Claude.

Examples:
    python main.py test_data/test.vcf
    python main.py test_data/test.bam --region chr1:1001-1100
    python main.py test_data/small.h5ad --column cell_type
    python main.py test_data/sample_R1.fastq --max-reads 100
"""
import argparse
import json

from biofile import bam_stats, fastq_stats, h5ad_summary, vcf_stats

parser = argparse.ArgumentParser(description="Summarize a genomics file")
parser.add_argument("path")
parser.add_argument("--region", help="BAM only, e.g. chr1:1001-1100")
parser.add_argument("--column", help="h5ad only, e.g. cell_type")
parser.add_argument("--max-reads", type=int, help="FASTQ only")
args = parser.parse_args()

name = args.path.lower()

if name.endswith((".vcf", ".vcf.gz", ".bcf")):
    result = vcf_stats(args.path)
elif name.endswith(".bam"):
    result = bam_stats(args.path, region=args.region)
elif name.endswith(".h5ad"):
    result = h5ad_summary(args.path, column=args.column)
elif name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
    result = fastq_stats(args.path, max_reads=args.max_reads)
else:
    raise SystemExit(f"Unsupported file type: {args.path}")

print(json.dumps(result, indent=2))
