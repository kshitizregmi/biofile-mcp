"""Command-line interface to the biofile summaries.

Detects the file format (from the extension, or from the file's contents if
the extension is unfamiliar) and prints the summary as JSON. Useful for
checking results directly, without an MCP client.

Examples::

    python main.py test_data/test.vcf
    python main.py test_data/test.bam --region chr1:1001-1100
    python main.py test_data/small.h5ad --column cell_type
    python main.py test_data/sample_R1.fastq --max-reads 100
"""
import argparse
import json

from biofile import bam_stats, detect_format, fastq_stats, h5ad_summary, vcf_stats


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize a VCF, BAM, FASTQ or h5ad file.")
    parser.add_argument("path", help="file to summarize")
    parser.add_argument("--region", help="BAM only: region such as chr1:1001-1100")
    parser.add_argument("--column", help="h5ad only: obs column to count, such as cell_type")
    parser.add_argument("--max-reads", type=int, help="FASTQ only: stop after this many reads")
    args = parser.parse_args()

    fmt = detect_format(args.path)

    if fmt == "VCF":
        result = vcf_stats(args.path)
    elif fmt == "BAM":
        result = bam_stats(args.path, region=args.region)
    elif fmt == "h5ad":
        result = h5ad_summary(args.path, column=args.column)
    elif fmt == "FASTQ":
        result = fastq_stats(args.path, max_reads=args.max_reads)
    else:
        parser.error(f"not a recognized VCF, BAM, FASTQ or h5ad file: {args.path}")

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
