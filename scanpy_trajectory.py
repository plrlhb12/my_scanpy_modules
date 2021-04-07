from google.colab import drive
drive.mount("/content/drive")

! pip3 install scanpy
! pip3 install leidenalg

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as pl
from matplotlib import rcParams
import os
os.chdir("/content/drive/Shared drives/CARD/projects/iNDI/line_prioritization/projects_lirong/Florian_data/")

sc.logging.print_versions()
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80,frameon=False, figsize=(3, 3))

# in the original tutorial, there are two comparisons on the graph visualization: "draw_graph on pca"   vs  "draw_graph on diffusion map"
# 1.draw_graph on pca space
# 2.draw_graph on diffusion map space
# 3. 

# 1. Use forced-directed graph drawing (draw_graph) on PCA space to visualize
# The regular flowchart is pca (19182, 50)----computing neigborhoods----tl.draw_graph(19182, 2)----pl.draw_graph
# It is a class of long-established algorithms for visualizing graphs, only two dimention FA1 and FA2 
# An alternative to tSNE and preserves the topology of the data better, and suggested for visualizing scRNA-seq data
adata = sc.read_h5ad("hypothalamic_exclude_heg.h5ad")
# if the adata has already calculated pca and neigborhoods, skip the firs two steps
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='leiden_0.6_peng', legend_loc='on data', title= "graph_peng", size=30, palette=sc.pl.palettes.vega_20_scanpy)


# 2. Denoising the graph using draw_graph on diffusion map space instead of PCA space to visualize
# The regular flowchart is dm (19182, 15)----computing neigborhoods----tl.draw_graph(19182, 2)----pl.draw_graph
# Computing neigbor distances on a few diffusion components amounts to denoising the graph - we just take a few of the first spectral components. 
# This is not a necessary step, neither for PAGA, nor clustering, nor pseudotime estimation.
# drawbacks: a lot of the branches are overplotted.
sc.tl.diffmap(adata)
# Use the parameter of key_added="diff_map_neighbor" if don't want to replace previous neighbors, connectivities, and distances genreated by calculating neighbors on PCA space
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')
sc.tl.draw_graph(adata)
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(10, 10))
sc.pl.draw_graph(adata, color='leiden_0.6_panglao', legend_loc='on data', title= "graph_panglao", size=30, palette=sc.pl.palettes.vega_20_scanpy)

# Using PAGA graph (partition-associated graph abstract), a coarse-grained and simplified (abstracted) graph. 
# Non-significant edges in the coarse- grained graph are thresholded away
# Use clusters computered by leiden or louvain
# it is better run these two commons together to avoid a potential error of index out of range
sc.tl.paga(adata, groups='leiden_0.6_panglao')
sc.pl.paga(adata, color=['leiden_0.6_panglao'] )


# Now is trajectory analysis
# Recomputing the embedding using PAGA-initialization

# the drwa_graph was init_pos on paga, which has been calculated using the groups="leiden_0.6_panglao"
# if we want to init use another categorical, e.g., "leiden_0.6", we have to recalculate and replot the paga using
#sc.tl.paga(adata, groups='leiden_0.6_peng')
#sc.pl.paga(adata, color='leiden_0.6_peng' )
sc.tl.draw_graph(adata, init_pos='paga')


