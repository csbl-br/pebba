import pandas as pd
from pebba.analysis.auxiliary_analysis import calculate_how_many_above_cut


def test_n_above_cut():
    fake_data = {
        "a": [1, 2, 3],
        "b": [5, 2, 0],
        "c": [1, 2, 2],
        "d": [1, 4, 3],
        "index": ["100", "200", "300"],
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
