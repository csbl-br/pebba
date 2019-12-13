import pandas as pd


def get_deg(deg_file, gene_col, logFC_col, pvalue_col):
    deg = read_deg(deg_file)
    deg = preprocess_deg(deg, gene_col, logFC_col, pvalue_col)
    return deg


def read_deg(file_in):
    """
    DEG = Differentially Expressed Genes

    :param file_in:
    :return:
    """
    if isinstance(file_in, str):
        deg = pd.read_csv(file_in, sep="\t")

    elif isinstance(file_in, pd.DataFrame):
        deg = file_in

    else:
        raise TypeError("Path to file or Pandas DataFrame expected.")

    # Remove rows that do not have a valid gene symbol
    return deg


def preprocess_deg(deg, gene_col, logFC_col, pvalue_col):
    deg = deg.dropna()
    deg = deg.rename(
        columns={gene_col: "Gene.symbol", logFC_col: "logFC", pvalue_col: "P.Value"}
    )
    return deg
