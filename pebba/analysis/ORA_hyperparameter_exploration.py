from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from functools import partial
from pebba.analysis.ORA import run_ORA
from pebba.analysis.top_genes import get_top_genes


def generate_ORA_dataframe(deg, genes_by_pathway, direction, min_genes, max_genes):
    """
    TODO: deixar usuario escolher nunmero de max workers
    """
    number_of_genes_in_deg = deg["Gene.symbol"].size
    top_genes = get_top_genes(deg, direction, max_genes)

    with ProcessPoolExecutor(max_workers=2) as executor:
        teste = partial(
            run_ORA_analysis_for_different_top_genes_ranges,
            top_genes=top_genes,
            genes_by_pathway=genes_by_pathway,
            number_of_genes_in_deg=number_of_genes_in_deg,
        )
        # TODO: eh possivel q as iteracoes tenham duracoes muito diferentes? comecar de tras pra frente parece uma ideia bem melhor na real
        # DONE: n pareceu ter uma melhora significativa
        # range(min_genes, max_genes + 1, 50)
        ORA_results = executor.map(teste, range(max_genes, min_genes - 1, -50))

    ORA_dataframe = pd.concat(ORA_results, axis=1, join="outer")
    ORA_dataframe.fillna(1.0)
    ORA_dataframe = (ORA_dataframe.apply(np.log10)) * (-1)

    return ORA_dataframe


def run_ORA_analysis_for_different_top_genes_ranges(
    i, top_genes, genes_by_pathway, number_of_genes_in_deg
):
    top_i_genes = top_genes.loc[0:i, "Gene.symbol"]
    ORA_result_i = run_ORA(top_i_genes, number_of_genes_in_deg, genes_by_pathway)
    ORA_result_i.columns = [str(i)]
    return ORA_result_i
