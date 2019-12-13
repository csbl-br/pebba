import pytest

from pebba.main import validate_range_of_inputs


@pytest.mark.parametrize(
    "min_genes, max_genes, p_cut",
    [
        (3, 2000, 0.1),
        (3000000, 2000, 0.1),
        (200, 7, 0.1),
        (200, 123451235, 0.1),
        (200, 700, 0),
        (200, 700, 2),
        (200, -700, 0.01),
    ],
)
def test_validate_range_of_inputs(min_genes, max_genes, p_cut):
    with pytest.raises(ValueError):
        validate_range_of_inputs(min_genes, max_genes, p_cut)
