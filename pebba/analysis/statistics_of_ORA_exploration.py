import pandas as pd
import numpy as np


def generate_ORA_statistics(path_table, p_cut, direction):
    # TODO: the only thing used in here is times_significant, talk to helder and delete everything we won't use
    df = pd.DataFrame()
    df["MaxR"] = path_table.max()
    df["SumR"] = path_table.sum()
    path_cut_p = np.log10(p_cut) * (-1)

    # How many pathways above path_cut_p (freq)
    how_many_pathways_above_cut = calculate_how_many_pathways_above_cut(
        path_table, path_cut_p, axis=0
    )
    n_rows = len(path_table.index)  ######
    df["times"] = how_many_pathways_above_cut / n_rows

    df.columns = [
        "maximum_MinuslogP_" + direction,
        "sum_MinuslogP_" + direction,
        "times_significant_" + direction,
    ]
    return df


def first_column_above_path_cut_p(row, path_cut_p):
    for cont, element in enumerate(row):
        if element > path_cut_p:
            return cont

    return 0


def calculate_how_many_pathways_above_cut(df, path_cut_p, axis):
    f = lambda x: x > path_cut_p
    how_many_pathways_above_cut = df.apply(f, axis=1).sum(axis=axis)
    return how_many_pathways_above_cut
