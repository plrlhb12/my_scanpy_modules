#! /usr/bin/env python

# cell annotation using the reference marker list and a rank_genes_groups at a specific resolution

import h5py
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import json

def ingest_indi(adata_brain, indi):
    
    adata_indi = sc.read_h5ad("../all/"+indi+"_concat_raw.h5ad")
    adata_brain.obs['leiden_0.6'] = adata_brain.obs['leiden_0.6_annotation']
    var_names = adata_brain.var_names.intersection(adata_indi.var_names)
    adata_brain = adata_brain[:, var_names]
    adata_indi = adata_indi[:, var_names]
    
    sc.pp.pca(adata_brain)
    sc.pp.neighbors(adata_brain)
    sc.tl.umap(adata_brain)
    
    sc.tl.ingest(adata_indi, adata_brain, obs='leiden_0.6')
    sc.pl.umap(adata_indi, color=['leiden_0.6'], save="_ingest_"+indi)
    print(adata_indi.obs["leiden_0.6"].value_counts())

def main():
    
    parser = argparse.ArgumentParser(description="Arguments for annotate cell types for clusters")
    # # optional argument
    # parser.add_argument("-i", "--input_file", type=str, help="path of the input of 'after_ranking_gene.h5ad'", default="after_ranking_gene.h5ad")
    parser.add_argument("-i", "--indi", type=str, help="path of the input of indi_concat_raw.h5ad")

    # parser.add_argument("-m", "--marker_ref_path", type=str, help="path of panglao reference markers", \
    #     default="/Users/pengl7/Desktop/scanpy_modules/reference_markers/marker_panglao_dic.p")
    # parser.add_argument("-o", "--out", type=str, help="path of the anndata object to be saved", default="after_annotated.h5ad")
    # parser.add_argument("-d", "--dpi", type=int, help="resolution of the output figure", default=80)
    # parser.add_argument("-s", "--figsize", type=float, nargs=2, help="size of output figure, use 2 numbers, e.g., 2 2")
    # parser.add_argument("-f", "--figure_type", type=str, help="define the export type of plot_type, e.g., png, pdf, or svg", default="pdf")
    # parser.add_argument("-p", "--project", type=str, help="give the project name", default="")
    # parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="default is show=True; provide no, false, or 0 to block print to screen")
    # parser.add_argument("-k", "--key", type=str, help="Choose the key of leiden clustering", default="leiden_0.4")
    # parser.add_argument("-r", "--ranking_key", type=str, help="Choose the key of ranking_genes_groups to be compared to marker_ref, \
    #     e.g., ranking_genes_groups", default="ranking_genes_groups_r0.6")
    # parser.add_argument("-n", "--new_cluster_names", type=str, nargs="+", help="provide the cell type name corresponding to each cluster")

    # # parser.set_defaults(func=annotate) #??
    args = parser.parse_args()
    indi = args.indi
    # # args.func(args)
    # input_file = args.input_file
    # marker_ref_path = args.marker_ref_path
    # dpi = args.dpi
    # out = args.out
    # figsize = args.figsize
    # figure_type = args.figure_type
    # show = args.show
    # project = args.project if (args.project == "") else ("_" + args.project)
    # ranking_key = args.ranking_key
    # key = args.key
    # new_cluster_names = args.new_cluster_names

    # print("\nThe arguments are: ", args)

    # set scanpy parameters
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_version_and_date()
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80, facecolor='white')

    # get the h5ad after ranking genes

    adata_brain = sc.read_h5ad("after_annotated.h5ad")

    ingest_indi(adata_brain, indi)


if __name__ == "__main__":
    main()