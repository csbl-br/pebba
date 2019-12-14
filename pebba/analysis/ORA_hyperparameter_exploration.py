import numpy as np
import pandas as pd
import sys

from pebba.analysis.ORA import run_ORA


def generate_ORA_dataframe(deg, dict_genes_by_pathway, direction, min_genes, max_genes):
    """
    """
    all_genes = deg["Gene.symbol"]

    top_genes = get_top_genes(deg, direction, max_genes)

    ORA_dataframe = run_ORA_analysis_for_different_top_genes_ranges(
        min_genes, max_genes, top_genes, all_genes, dict_genes_by_pathway
    )
    # Aqui tinha toda a parte de fazer sumario,
    # dar merge com df original e ordenar ele por "first top cut significant" e dropar tudo depois.
    # Como essa ordenacao e sumario n fazem o menor sentido eu deletei tudo.

    return ORA_dataframe


def get_top_genes(deg, direction, max_genes):
    top_genes = sort_deg_by_direction(direction, deg, max_genes)
    top_genes["pi_value"] = calculate_pi_value(top_genes)
    top_genes = sort_top_genes_by_pi_value(top_genes)
    return top_genes


def sort_deg_by_direction(direction, deg_list, max_genes):
    if direction == "up":
        top_genes = deg_list.sort_values(by="logFC", ascending=False).head(n=max_genes)
    elif direction == "down":
        top_genes = deg_list.sort_values(by="logFC", ascending=True).head(n=max_genes)
    elif direction == "any":
        deg_list["logFC"] = deg_list["logFC"].astype(np.float64)
        deg_list["logFC"] = deg_list["logFC"].abs()
        top_genes = deg_list.sort_values(by="logFC", ascending=True).head(n=max_genes)
    else:
        sys.exit("Invalid direction argument")
    top_genes["Gene.symbol"] = top_genes["Gene.symbol"].astype(str)
    return top_genes


def calculate_pi_value(top_genes):
    return top_genes["logFC"].apply(abs) * (-top_genes["P.Value"].apply(np.log10))


def sort_top_genes_by_pi_value(top_genes):
    return top_genes.sort_values(by="pi_value", ascending=False).reset_index(drop=True)


def run_ORA_analysis_for_different_top_genes_ranges(
    min_genes, max_genes, top_genes, all_genes, gmt_file
):
    ORA_results = []
    for i in range(
        min_genes, max_genes + 1, 50
    ):  # max_genes +1 to include the number max_genes in the range
        top_i_genes = top_genes.loc[0:i, "Gene.symbol"]
        ORA_result_i = run_ORA(top_i_genes, all_genes, gmt_file)
        ORA_result_i.columns = [str(i)]  # ["term", str(i)]
        # ORA_result_i = ORA_result_i.set_index("term", drop=True)

        ORA_results.append(ORA_result_i)

    ORA_dataframe = pd.concat(ORA_results, axis=1, join="outer")

    ORA_dataframe.fillna(1.0)
    ORA_dataframe = (ORA_dataframe.apply(np.log10)) * (-1)

    return ORA_dataframe
