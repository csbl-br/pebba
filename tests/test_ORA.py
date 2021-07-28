import pandas as pd
from pebba.analysis.auxiliary_analysis import calculate_how_many_above_cut
from pebba.analysis.ORA_hyperparameter_exploration import _reorder_df_by_mean_expression


def test_n_above_cut():
    fake_data = {
        "100": [1, 2, 3],
        "200": [5, 2, 0],
        "300": [1, 2, 2],
        "400": [1, 4, 3],
        "index": ["a", "b", "c"],
    }
    fake_df = pd.DataFrame(fake_data).set_index("index")
    # by cuts
    assert all(
        [1, 1, 2] == calculate_how_many_above_cut(fake_df, path_cut_p=2.5, axis_sum=1)
    )
    # by pathways
    assert all(
        [1, 1, 0, 2]
        == calculate_how_many_above_cut(fake_df, path_cut_p=2.5, axis_sum=0)
    )


def test_reorder_df():
    fake_data = {
        "100": [1, 7, 3],
        "200": [1, 7, 3],
        "300": [1, 7, 3],
        "400": [1, 7, 3],
        "index": ["a", "b", "c"],
    }
    fake_df = pd.DataFrame(fake_data).set_index("index")

    assert all(["b", "c", "a"] == _reorder_df_by_mean_expression(fake_df).index)
