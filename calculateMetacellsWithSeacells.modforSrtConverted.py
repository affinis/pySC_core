# numpy                     1.24.4
# pandas                    2.0.3
# scanpy                    1.9.8
# seaborn                   0.13.2
# SEACells                  0.3.3
# matplotlib                3.7.5

import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 256

#parse arguments
parser = argparse.ArgumentParser(description="calculate metacell with seurat object converted file")

parser.add_argument("-i", "--input", type=str, help="seurat object converted filtered feature bc matrix dir name with features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz")
parser.add_argument("-r", "--ratio", type=float, help="num metacells/num total cells",default=0.01)
parser.add_argument("-a", "--cell_annotations", type=str, help="an annotation file")
parser.add_argument("-o", "--output", type=str, help="output directory")
parser.add_argument("-p", "--cell_prefix", type=str, help="sample prefix")

args = parser.parse_args()

# Path to the Cell Ranger output directory
cellranger_dir = args.input

# Read the data
print("reading sc data...")
adata=sc.read_10x_mtx(cellranger_dir, var_names='gene_symbols')

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Saving count data
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
print("normalizing data...")
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

# finding highly variable genes
print("finding highly variable genes...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# performing PCA
print("performing PCA...")
sc.tl.pca(adata)

# compute neighborhood graph
sc.pp.neighbors(adata)

# run UMAP
print("performing UMAP...")
sc.tl.umap(adata)

# read cell annotations
print("reading cell annotation...")
cell_annotations = pd.read_csv(args.cell_annotations, sep='\t')

# add prefix to cell barcodes
adata.obs.index=adata.obs.index.astype(str)

# Check if all cells in the TSV exist in adata
missing_cells = set(cell_annotations['index']) - set(adata.obs_names)
if missing_cells:
    print(f"Warning: {len(missing_cells)} cells in TSV not found in AnnData.")

# Filter annotations to only include cells present in adata
cell_annotations = cell_annotations[cell_annotations['index'].isin(adata.obs_names)]

# Set cell_id as index (if not already)
cell_annotations = cell_annotations.set_index('index')

# Ensure the order matches adata.obs
cell_annotations = cell_annotations.reindex(adata.obs_names)

# Add cell_type to adata.obs
adata.obs['manual.level1'] = cell_annotations['manual.level1']
adata.obs['manual.level2'] = cell_annotations['manual.level2']
adata.obs['predicted.celltype.l1.5'] = cell_annotations['predicted.celltype.l1.5']
adata.obs['manual_NI'] = cell_annotations['manual_NI']

# Copy the counts to ".raw" attribute of the anndata since it is necessary for downstream analysis
# This step should be performed after filtering 
raw_ad = sc.AnnData(adata.X)
raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
adata.raw = raw_ad

## Core parameters 
print("start to run seacells...")
n_SEACells = round(args.ratio*adata.n_obs)
build_kernel_on = 'X_pca' # key in ad.obsm to use for computing metacells
                          # This would be replaced by 'X_svd' for ATAC data

## Additional parameters
n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells

model = SEACells.core.SEACells(adata, 
                  build_kernel_on=build_kernel_on, 
                  n_SEACells=n_SEACells, 
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5)

model.construct_kernel_matrix()
M = model.kernel_matrix

# Initialize archetypes
model.initialize_archetypes()

# Plot the initilization to ensure they are spread across phenotypic space
# SEACells.plot.plot_initialization(adata, model)

model.fit(min_iter=10, max_iter=50)

SEACell_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')
SEACell_soft_ad = SEACells.core.summarize_by_soft_SEACell(adata, model.A_, celltype_label='predicted.celltype.l1.5',
                                                          summarize_layer='raw', minimum_weight=0.05)

# Save counts as CSV (genes × cells)
pd.DataFrame(
    SEACell_ad.X.toarray().T,
    index=SEACell_ad.var_names,
    columns=SEACell_ad.obs_names
).to_csv(args.output+"/"+args.cell_prefix+"_"+"metacells.csv")

# Save metadata
adata.obs[["SEACell"]]=args.cell_prefix+"_"+adata.obs[["SEACell"]]
pd.DataFrame(adata.obs).to_csv(args.output+"/"+args.cell_prefix+"_"+"metadata.csv")

# Save metacellplot
plot2d_name=args.output+"/"+args.cell_prefix+"_"+"metacells2Dplot.pdf"
SEACells.plot.plot_2D(adata, key='X_umap', colour_metacells=True,show=False,save_as=plot2d_name)

SEACell_size_name=args.output+"/"+args.cell_prefix+"_"+"metacells_size_plot.pdf"
SEACells.plot.plot_SEACell_sizes(adata, bins=5,show=False,save_as=SEACell_size_name)
