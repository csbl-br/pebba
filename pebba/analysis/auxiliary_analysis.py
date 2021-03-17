def calculate_how_many_above_cut(df, path_cut_p, axis_sum=0):
    n_above_cut = (df > path_cut_p).sum(axis=axis_sum)
    return n_above_cut
