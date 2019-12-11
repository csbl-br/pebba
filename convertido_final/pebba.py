import numpy as np
import pandas as pd
from pypathway import ORA, GMTUtils
import sys


def _cutoff_path(path_table, p_cut, direction):
    df_index = path_table.columns
    df = pd.DataFrame()
    df["MaxR"] = path_table.max()
    df["SumR"] = path_table.sum()
    path_cut_p = np.log10(p_cut) * (-1)

    # How many pathways above path_cut_p (freq)
    how_many_pathways_above_cut = calculate_how_many_pathways_above_cut(path_table, path_cut_p, axis=0)
    n_rows = len(path_table.index)  ######
    df["times"] = how_many_pathways_above_cut / n_rows

    df.columns = ["maximum_MinuslogP_" + direction,
                  "sum_MinuslogP_" + direction,
                  "times_significant_" + direction]
    return df


###############################################################################################################################################

def _get_pathway(deg, direction, min_genes, max_genes, p_cut):
    """
    TODO: gmt information, who passes to who and when and how
    """
    all_genes = deg["Gene.symbol"]  ## Empty values (non-annotated genes) will be removed

    top_genes = get_top_genes(direction, deg, max_genes)
    top_genes["pi_value"] = calculate_pi_value(top_genes)
    top_genes = sort_top_genes_by_pi_value(top_genes)

    ORA_dataframe = get_ORA_dataframe(min_genes, max_genes, top_genes, all_genes, gmt_file)
    ORA_dataframe.fillna(1.0)  # acho q n eh mais necessario, ORA faz sozinho
    ORA_dataframe = (ORA_dataframe.apply(np.log10)) * (-1)

    # TODO: refatorar daqui em diante
    ORA_summary = summarizes_ORA_information(ORA_dataframe, p_cut, direction)

    ORA_dataframe = pd.concat([ORA_summary, ORA_dataframe], axis=1)

    ORA_dataframe = ORA_dataframe.sort_values(by="FirstTopCut_significant_" + direction, ascending=False)
    ORA_dataframe = ORA_dataframe.drop(labels=["TopCut_highestMinuslogP_" + direction,
                                               "maximum_MinuslogP_" + direction,
                                               "sum_MinuslogP_" + direction,
                                               "times_significant_" + direction,
                                               "FirstTopCut_significant_" + direction,
                                               "PEBBA_score_" + direction], axis=1)
    return ORA_summary, ORA_dataframe


def get_top_genes(direction, deg_list, max_genes):
    if (direction == "up"):
        top_genes = deg_list.sort_values(by="logFC", ascending=False).head(n=max_genes)
    elif (direction == "down"):
        top_genes = deg_list.sort_values(by="logFC", ascending=True).head(n=max_genes)
    elif (direction == "any"):
        deg_list["logFC"] = deg_list["logFC"].astype(np.float64)
        deg_list["logFC"] = deg_list["logFC"].abs()
        top_genes = deg_list.sort_values(by="logFC", ascending=True).head(n=max_genes)
    else:
        sys.exit("Invalid direction argument")
    top_genes["Gene.Symbol"] = top_genes["Gene.Symbol"].astype(str)
    return top_genes


def calculate_pi_value(top_genes):
    return top_genes["logFC"].apply(abs) * (- top_genes["P.Value"].apply(np.log10))


def sort_top_genes_by_pi_value(top_genes):
    return top_genes.sort_values(by="pi_value", ascending=False).reset_index(drop=True)


def get_ORA_dataframe(min_genes, max_genes, top_genes, all_genes, gmt_file):
    ORA_results = []
    for i in range(min_genes, max_genes + 1, 50):  # max_genes +1 to include the number max_genes in the range
        top_i_genes = top_genes.loc[0:i, "Gene.symbol"]
        ORA_result_i = _run_enrich(top_i_genes, all_genes, gmt_file)
        ORA_result_i.columns = ["term", str(i)]
        ORA_result_i = ORA_result_i.set_index("term", drop=True)

        ORA_results.append(ORA_result_i)

    return pd.concat(ORA_results, axis=1, join="outer")


def summarizes_ORA_information(merge_p2, p_cut, direction):
    path_cut_p = np.log10(p_cut) * (-1)

    NG = merge_p2.idxmax(axis=1)  # O recorte de genes q apresentou o maior p valor possui NG genes
    NG = NG.astype(np.int64)
    p_max = merge_p2.max(axis=1)
    p_sum = merge_p2.sum(axis=1)

    num_columns_merge_p2 = merge_p2.shape[1]
    how_many_pathways_above_cut = calculate_how_many_pathways_above_cut(merge_p2, path_cut_p, axis=1)

    times = how_many_pathways_above_cut / num_columns_merge_p2

    ES3 = (1 - np.exp(- p_max) / (1 + (0.1 * np.sqrt(NG))))

    first = merge_p2.apply(first_column_above_path_cut_p, axis=1, path_cut_p=path_cut_p)
    first = first.apply(lambda x: merge_p2.columns[x] if x != 0 else 0)

    dicionario = {"TopCut_highestMinuslogP_" + direction: NG,
                  "maximum_MinuslogP_" + direction: p_max,
                  "sum_MinuslogP_" + direction: p_sum,
                  "times_significant_" + direction: times,
                  "FirstTopCut_significant_" + direction: first,
                  "PEBBA_score_" + direction: ES3}

    df = pd.DataFrame(dicionario)
    df["FirstTopCut_significant_" + direction] = df["FirstTopCut_significant_" + direction].astype(np.int64)

    return df

def _run_enrich(top_genes, all_genes, gmt_file):
    term2gene = GMTUtils.parse_gmt_file(gmt_file)
    df = ORA.run(top_genes, all_genes, term2gene).df
    df = df[["name", "fdr"]]
    return df


##############################################################

def first_column_above_path_cut_p(row, path_cut_p):
    for cont, element in enumerate(row):
        if element > path_cut_p:
            return cont

    return 0


def calculate_how_many_pathways_above_cut(df, path_cut_p, axis):
    f = lambda x: x > path_cut_p
    how_many_pathways_above_cut = df.apply(f, axis=1).sum(axis=axis)
    return how_many_pathways_above_cut
