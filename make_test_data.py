"""Create small test files with known answers in ./test_data."""
import random
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pysam import samtools

out = Path("test_data")
out.mkdir(exist_ok=True)

# VCF: 3 SNPs (2 transitions, 1 transversion), 2 indels
(out / "test.vcf").write_text(
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=248956422>\n"
    '##FILTER=<ID=LowQual,Description="Low quality">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n"
    "chr1\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n"
    "chr1\t2000\t.\tC\tT\t60\tPASS\t.\tGT\t1/1\t0/1\n"
    "chr1\t3000\t.\tG\tC\t20\tLowQual\t.\tGT\t0/1\t./.\n"
    "chr1\t4000\t.\tAT\tA\t45\tPASS\t.\tGT\t0/1\t0/0\n"
    "chr1\t5000\t.\tC\tCAG\t70\tPASS\t.\tGT\t1/1\t0/0\n"
)

# BAM: 20 normal reads, 3 duplicates, 2 unmapped
seq, qual = "A" * 50, "I" * 50
sam = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chr1\tLN:248956422", "@RG\tID:rg1\tSM:TUMOR"]
sam += [f"read{i}\t0\tchr1\t1001\t60\t50M\t*\t0\t0\t{seq}\t{qual}\tRG:Z:rg1" for i in range(20)]
sam += [f"dup{i}\t1024\tchr1\t1001\t60\t50M\t*\t0\t0\t{seq}\t{qual}\tRG:Z:rg1" for i in range(3)]
sam += [f"unmapped{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\tRG:Z:rg1" for i in range(2)]
(out / "test.sam").write_text("\n".join(sam) + "\n")
samtools.sort("-o", str(out / "test.bam"), str(out / "test.sam"))
samtools.index(str(out / "test.bam"))
(out / "test.sam").unlink()

# FASTQ: 1000 reads, N in every 50th read, 4 lowercase bases in every 100th read
random.seed(1)
with open(out / "sample_R1.fastq", "w") as f:
    for i in range(1000):
        L = random.randint(240, 251)
        s = "".join(random.choice("ACGT") for _ in range(L))
        if i % 50 == 0:
            s = s[:10] + "N" + s[11:]
        if i % 100 == 0:
            s = s[:20] + "acgt" + s[24:]
        q = "".join(chr(33 + (37 if j < L * 0.7 else random.randint(10, 30))) for j in range(L))
        f.write(f"@read{i}\n{s}\n+\n{q}\n")

# h5ad: 200 cells x 100 genes, raw counts, 3 cell types (100/60/40)
rng = np.random.default_rng(0)
counts = sp.csr_matrix(rng.poisson(1.0, size=(200, 100)).astype(np.float32))
obs = pd.DataFrame(
    {"cell_type": pd.Categorical(["T cell"] * 100 + ["B cell"] * 60 + ["NK cell"] * 40)},
    index=[f"cell_{i}" for i in range(200)],
)
adata = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=[f"GENE{i}" for i in range(100)]))
adata.obsm["X_umap"] = rng.normal(size=(200, 2))
adata.write_h5ad(out / "small.h5ad")

print("Created:", sorted(p.name for p in out.iterdir()))
