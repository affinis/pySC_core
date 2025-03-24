# anndata 0.11.3
# scanpy 1.11.0
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description="A script that convert .loom file into .h5ad file.")

parser.add_argument("-i", "--input", type=str, help="cellranger filtered feature bc matrix dir name with features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz")
parser.add_argument("-o", "--output", type=str, help="h5ad filename, with .h5ad")

args = parser.parse_args()

# Path to the Cell Ranger output directory
cellranger_dir = args.input

# Read the data
adata = sc.read_10x_mtx(cellranger_dir, var_names='gene_symbols')

adata.write(args.output)
