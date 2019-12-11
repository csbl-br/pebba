import os
import sys

import pandas as pd

from convertido_final.pebba import _cutoff_path
from convertido_final.pebba import _get_pathway
from convertido_final.utils import read_gmt_hier
from convertido_final.visualization import create_interactive_plot


def pebba(deg_file,
          gmt_file,
          gene_col="Gene.symbol",
          logFC_col="logFC",
          pvalue_col="P.Value",
          min_genes=100,
          max_genes=1500,
          p_cut=0.2,
          verbose=True,
          analysis_name=None,
          results_dir="Results",
          force=False):
    validate_range_of_inputs(min_genes, max_genes, p_cut)

    deg = get_deg(deg_file)
    deg = rename_deg_columns(deg, gene_col, logFC_col, pvalue_col)

    # Get information from all unique terms
    gene_to_pathway_relationship_df = read_gmt_hier(gmt_file)  # hier is useless,
    # TODO: pass the utils gmt reader into the run_enrich function to improve modularity
    # TODO: gene_to_pathway_relationship_df and dict_genes_by_pathway contain the same information, transform into one thing only
    f = lambda df: [gene for gene in df["gene"]]
    dict_genes_by_pathway = dict(gene_to_pathway_relationship_df.groupby(["term"]).apply(f))

    directions = ["up", "down", "any"]

    dfs = dict()
    paths = dict()
    cut_paths = dict()
    for direction in directions:
        if verbose:
            print(direction + "\nGetting Pathways")
        dfs[direction], paths[direction] = _get_pathway(deg, direction, min_genes, max_genes, p_cut)
        if verbose:
            print("Getting Pathway Cutoff")

        cut_paths[direction] = _cutoff_path(paths[direction], p_cut, direction)

    if verbose:
        print("Saving heatmaps")

    create_results_directory(results_dir, force)
    if analysis_name is None:
        analysis_name = set_analysis_name(deg_file)

    for direction in directions:
        create_interactive_plot(paths[direction], dict_genes_by_pathway, direction, analysis_name, cut_paths[direction],
                                results_dir)

    if verbose:
        print("Done")


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


def validate_range_of_inputs(min_genes, max_genes, p_cut):
    if min_genes < 50 or min_genes > 2900:
        sys.exit("Variable min_genes must be between 50 and 2900 genes")

    if max_genes < 100 or max_genes > 3000:
        sys.exit("Variable max_genes must be between 100 and 3000 genes")

    if p_cut < 0.00001 or p_cut > 1:
        sys.exit("Variable p_cut must be between 0.00001 and 1")


def create_results_directory(results_dir, force):
    results_dir = os.path.abspath(results_dir)
    if not os.path.exists(results_dir):
        os.makedirs("Results/Tables")
        os.makedirs("Results/Heatmaps")
    else:
        if not force:
            sys.exit("Stopping analysis: " + results_dir + " already exists! Use force=True to overwrite.")


def get_deg(file_in):
    """
    TODO: do a better job at documentation in this function,
    TODO: take the name of the columns to rename everything and stop having to pass around column names as parameters
    DEG = Differentially Expressed Genes

    :param file_in:
    :return:
    """
    if isinstance(file_in, str):
        deg_list = pd.read_csv(file_in, sep="\t")

    elif isinstance(file_in, pd.DataFrame):
        deg_list = file_in

    else:
        raise TypeError("Path to file or Pandas DataFrame expected.")

    # Remove rows that do not have a valid gene symbol
    return deg_list.dropna()


def rename_deg_columns(deg, gene_col, logFC_col, pvalue_col):
    deg = deg.rename(columns={gene_col: "Gene.symbol",
                              logFC_col: "logFC",
                              pvalue_col: "P.Value"})

    return deg


def set_analysis_name(file_in):
    """
    If a name is not provided, generate one from the file name or give the generic name "PEBBA_analysis"
    :param file_in: string or Pandas DataFrame
    :return: string
    """
    if isinstance(file_in, str):
        return os.path.splitext(os.path.basename(file_in))[0]  # pega o basename e tira a extensao
    else:
        return "PEBBA_analysis"
