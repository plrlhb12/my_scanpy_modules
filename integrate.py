#!/usr/bin/env python3

import argparse
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import inspect

def read(*files):
    """
    Read multiple files into a adata list.
    """
    print("read these files", files)
    adatas=[]
    for i in files:
        if i.endswith(".h5"):
            adata = sc.read_10x_h5(i)
        elif i.endswith(".h5ad"):
            adata = sc.read_h5ad(i)
        else:
            # There is a potential bug in saving file later if use cache=True
            # here input is the directory name which contains matrix file
            adata = sc.read_10x_mtx(i, var_names='gene_symbols', cache=True)
        adata.var_names_make_unique()
        adatas.append(adata)
    return adatas

def concat(adatas, batch_categories=None, project=None):
    """
    Concatenate several adata into 1.
    """
    print("Concatente", len(adatas), "files")
    print()
    print(adatas[0])
    adata_concat = adatas[0].concatenate(*adatas[1:], join="inner", batch_categories=batch_categories)
      
    if project is None:
        adata_concat.write_h5ad("concat_raw.h5ad")
    else:
        adata_concat.write_h5ad(project+"_concat_raw.h5ad")
        
    return adata_concat

def ingest(adatas):
    """
    Integrate adata using the method of ingest. Please define the reference data.
    """
    pass

def bbknn(adatas):
    """
    Integrate adata using the method of bbknn.
    """
    pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", help="provide the files to be combined")
    parser.add_argument("-b", "--batch_categories", nargs="*",help="provide the batch labels to save, default is numbers if not providing any")
    parser.add_argument("-p", "--project", help="provide the name prefix to save")

    args = parser.parse_args()
    print("Theses are the arguments used", *vars(args))
    print()
    files = args.files
    batch_categories = args.batch_categories
    project = args.project

    adatas = read(*files)
    concat(adatas, batch_categories=batch_categories, project=project)


if __name__ == "__main__":
    main()