# scanpy 1.11.0
# metacells 0.9.5
# refering to 'https://tanaylab.github.io/metacells-vignettes/iterative.html' for building this pipline

import scanpy as sc
import metacells as mc
import argparse

parser = argparse.ArgumentParser(description="A script that performing metacells with .h5ad file")

parser.add_argument("-i", "--input", type=str, help=".h5ad file")
parser.add_argument("-o", "--output", type=str, help=".h5ad file")

args = parser.parse_args()

# Load data
adata = sc.read_h5ad(args.input)

# The next decision we need to make is which genes to exclude from the data, by their names. The poster children for this are mytochondrial genes and strong sex-chromosome genes.
EXCLUDED_GENE_NAMES = ["XIST", "MALAT1"]  # Sex-specific genes.
EXCLUDED_GENE_PATTERNS = ["MT-.*"]        # Mitochondrial.

mc.pl.exclude_genes(adata, excluded_gene_names=EXCLUDED_GENE_NAMES, excluded_gene_patterns=EXCLUDED_GENE_PATTERNS)

# The next decision we need to make is which cells to exclude due to containing too many UMIs in the excluded genes. If a cell contains "too many" excluded (mainly mytochondrial) gene UMIs, this may indicate a badly sampled cell, leading to very skewed results. Again, the exact threshold depends on both the technology and the dataset. Here we resort to looking at the distribution of the fraction of excluded genes in each cell, and manually picking the threshold.
mc.tl.compute_excluded_gene_umis(adata)