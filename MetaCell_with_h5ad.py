# scanpy 1.11.0
# anndata 0.11.3
# metacells 0.9.5
# refering to 'https://tanaylab.github.io/metacells-vignettes/iterative.html' for building this pipline

import scanpy as sc
import argparse
import anndata as ad             # For reading/writing AnnData files
#import matplotlib.pyplot as plt  # For plotting
import metacells as mc           # The Metacells package
import numpy as np               # For array/matrix operations
import pandas as pd              # For data frames
import os                        # For filesystem operations
import seaborn as sb             # For plotting
import scipy.sparse as sp        # For sparse matrices
import shutil                    # for filesystem operations
from math import hypot           # For plotting
from typing import *             # For type annotations

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

# pick a maximal fraction of excluded UMIs in each cell.
PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION = 0.25

# exclude the cells we have chosen
mc.pl.exclude_cells(adata, properly_sampled_min_cell_total=PROPERLY_SAMPLED_MIN_CELL_TOTAL,
                    properly_sampled_max_cell_total=PROPERLY_SAMPLED_MAX_CELL_TOTAL,
                    properly_sampled_max_excluded_genes_fraction=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION,
)

# extract clean data
clean = mc.pl.extract_clean_data(full, name="clean.data")
mc.ut.top_level(clean)
print(f"Clean: {clean.n_obs} cells, {clean.n_vars} genes")

# compute the 1st iteration metacells
cells = clean
clean = None  # Allow it to be gc-ed
mc.ut.set_name(cells, "data.iteration-1.cells")
print(f"Iteration 1: {cells.n_obs} cells, {cells.n_vars} genes")

# A crucial decision when running metacells is the list of genes are lateral, that is, should not be used to group cells together. The poster child for this are cell-cycle genes. These genes are strong and any clustering algorithm will therefore prefer to group together cells in the same cell-cycle state, at the expense of mixing up (reasonably close) other cell states, which are what we are actually interested in. Note that lateral genes are still used in deviant cell detection, that is, each lateral gene should still have a consistent level in all the cells in each metacell.
BASE_LATERAL_GENE_NAMES = [
    "AURKA", "MCM3", "MCM4", "MCM7", "MKI67", "PCNA", "RRM2", "SMC4", "TPX2",  # Cell-cycle
    "FOS", "HSP90AB1", "TXN",                                                  # Stress
]
BASE_LATERAL_GENE_PATTERNS = ["RP[LS].*"]  # Ribosomal
# We'll reuse this through the iterations.
# It is just a thin wrapper for mark_lateral_genes,
# and optionally also shows the results.
def update_lateral_genes(
    *,
    names: List[str] = [],
    patterns: List[str] = [],
    op: str = "set",
    show: bool = True
) -> None:
    mc.pl.mark_lateral_genes(
        cells,
        lateral_gene_names=names,
        lateral_gene_patterns=patterns,
        op=op
    )

    if not show:
        return
    
    lateral_genes_mask = mc.ut.get_v_numpy(cells, "lateral_gene")
    lateral_gene_names = set(cells.var_names[lateral_genes_mask])
    
    print(sorted([
        name for name in lateral_gene_names
        if not name.startswith("RPL") and not name.startswith("RPS")
    ]))

    print(f"""and {len([
        name for name in lateral_gene_names if name.startswith("RPL") or name.startswith("RPS")
    ])} RP[LS].* genes""")

update_lateral_genes(names=BASE_LATERAL_GENE_NAMES, patterns=BASE_LATERAL_GENE_PATTERNS)

# Since our initial list was very partial, we would like to extend it to include any genes highly-correlated with them (as long as this makes sense, biologically speaking).
mc.pl.relate_to_lateral_genes(cells)

# Finally, we need to decide on how much parallelization to use. This is a purely technical decision - it does not affect the results, just the performance.
# Either use the guesstimator:
max_parallel_piles = mc.pl.guess_max_parallel_piles(cells)
# Or, if running out of memory manually override:
# max_paralle_piles = ...
print(max_parallel_piles)
mc.pl.set_max_parallel_piles(max_parallel_piles)

# Assigning cells to metacells
with mc.ut.progress_bar():
    mc.pl.divide_and_conquer_pipeline(cells)

# Collecting the metacells
metacells = mc.pl.collect_metacells(cells, name="data.iteration-1.metacells")
print(f"Iteration 1: {metacells.n_obs} metacells, {metacells.n_vars} genes")

#Saving the data
metacells.write_h5ad("../output/iterative/iteration-1/test.metacells.h5ad")
    






