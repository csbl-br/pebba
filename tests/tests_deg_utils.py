# import pytest
import pandas as pd

from pebba.utils.deg_utils import preprocess_deg


def test_preprocess_deg():
    fake_data = {"a": [1, 2, 3], "b": [1, 2, 3], "c": [1, 2, 3], "d": [1, 2, 3]}
    fake_deg = pd.DataFrame(fake_data)

    preprocessed_fake_deg = preprocess_deg(
        fake_deg, gene_col="c", logFC_col="a", pvalue_col="d"
    )
    processed_names = list(preprocessed_fake_deg.columns)
    intended_names = ["logFC", "b", "Gene.symbol", "P.Value"]
    for processed_name, intended_name in zip(processed_names, intended_names):
        assert processed_name == intended_name
