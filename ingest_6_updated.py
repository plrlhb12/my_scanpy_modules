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

def ingest_indi(adata_ref, adata_test, key=None, out=None):
    """
    Project reference annoations to the test dataset
    """
    # temporaray change the leiden_r key using the annotated cell types.
    # however, if the reference adata hasn't been annotated, comment this line
    adata_ref.obs[key] = adata_ref.obs[key+'_annotation']

    var_names = adata_ref.var_names.intersection(adata_test.var_names)
    adata_ref = adata_ref[:, var_names]
    adata_test = adata_test[:, var_names]
    
    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)
    
    sc.tl.ingest(adata_test, adata_ref, obs=key)
    sc.pl.umap(adata_test, color=[key], save="_ingest_"+test)
    print(adata_test.obs[key].value_counts())
    adata_test.write_h5ad(out)

def main():
    
    parser = argparse.ArgumentParser(description="Arguments for annotate cell types for clusters")
    # # optional argument
    # parser.add_argument("-i", "--input_file", type=str, help="path of the input of 'after_ranking_gene.h5ad'", default="after_ranking_gene.h5ad")
    parser.add_argument("-r", "--reference", type=str, help="path of the reference h5ad", default="after_annotated.h5ad")
    parser.add_argument("-t", "--test", type=str, help="path of the surfix of the test h5ad")

    # parser.add_argument("-m", "--marker_ref_path", type=str, help="path of panglao reference markers", \
    #     default="/Users/pengl7/Desktop/scanpy_modules/reference_markers/marker_panglao_dic.p")
    parser.add_argument("-o", "--out", type=str, help="path of the test anndata object to be saved", default="after_annotated_ingested.h5ad")
    # parser.add_argument("-d", "--dpi", type=int, help="resolution of the output figure", default=80)
    # parser.add_argument("-s", "--figsize", type=float, nargs=2, help="size of output figure, use 2 numbers, e.g., 2 2")
    # parser.add_argument("-f", "--figure_type", type=str, help="define the export type of plot_type, e.g., png, pdf, or svg", default="pdf")
    # parser.add_argument("-p", "--project", type=str, help="give the project name", default="")
    # parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="default is show=True; provide no, false, or 0 to block print to screen")
    parser.add_argument("-k", "--key", type=str, help="Choose the key of leiden clustering", default="leiden_0.4")
    # parser.add_argument("-r", "--ranking_key", type=str, help="Choose the key of ranking_genes_groups to be compared to marker_ref, \
    #     e.g., ranking_genes_groups", default="ranking_genes_groups_r0.6")
    # parser.add_argument("-n", "--new_cluster_names", type=str, nargs="+", help="provide the cell type name corresponding to each cluster")

    # # parser.set_defaults(func=annotate) #??
    args = parser.parse_args()
    reference = args.reference
    test = args.test
    key = args.key
    out = args.out
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
    adata_ref = sc.read_h5ad(reference)
    adata_test = sc.read_h5ad(test)
    ingest_test(adata_ref, adata_test, key=key, out=out)


if __name__ == "__main__":
    main()

# example:
# python ingest_6_updated.py -r after_annotated_brain.h5ad -t NGN2_concat_raw.h5ad -o NGN2_after_annotated_ingested.h5ad -k leiden_0.6
