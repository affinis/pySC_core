# loompy 3.0.6
# scanpy 1.11.0
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description="A script that convert .loom file into .h5ad file.")

parser.add_argument("-i", "--input", type=str, help="loom file")
parser.add_argument("-o", "--output", type=str, help="h5ad filename")

args = parser.parse_args()

# Load the Loom file
adata = sc.read_loom(args.input)

# Save as H5AD
adata.write(args.output)