"""Command-line interface to the biofile summaries.

Picks the summary function from the file extension and prints the result as
JSON. Useful for checking results directly, without an MCP client.

Examples::

    python main.py test_data/test.vcf
    python main.py test_data/test.bam --region chr1:1001-1100
    python main.py test_data/small.h5ad --column cell_type
    python main.py test_data/sample_R1.fastq --max-reads 100
"""
import argparse
import json

from biofile import bam_stats, fastq_stats, h5ad_summary, vcf_stats

VCF_SUFFIXES = (".vcf", ".vcf.gz", ".bcf")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize a VCF, BAM, FASTQ or h5ad file.")
    parser.add_argument("path", help="file to summarize")
    parser.add_argument("--region", help="BAM only: region such as chr1:1001-1100")
    parser.add_argument("--column", help="h5ad only: obs column to count, such as cell_type")
    parser.add_argument("--max-reads", type=int, help="FASTQ only: stop after this many reads")
    args = parser.parse_args()

    name = args.path.lower()

    if name.endswith(VCF_SUFFIXES):
        result = vcf_stats(args.path)
    elif name.endswith(".bam"):
        result = bam_stats(args.path, region=args.region)
    elif name.endswith(".h5ad"):
        result = h5ad_summary(args.path, column=args.column)
    elif name.endswith(FASTQ_SUFFIXES):
        result = fastq_stats(args.path, max_reads=args.max_reads)
    else:
        parser.error(f"unsupported file type: {args.path}")

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
