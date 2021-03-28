from scipy import stats
from statsmodels.sandbox.stats import multicomp
import pandas as pd
import numpy as np


def run_ORA(top_genes, number_of_genes_in_deg, genes_by_pathway, pathway_sizes):
    size_of_top_genes = len(top_genes)
    degs_by_pathway = get_number_of_degs_by_pathway(top_genes, genes_by_pathway)
    result = run_hypergeom_tests(
        number_of_genes_in_deg,
        size_of_top_genes,
        degs_by_pathway,
        pathway_sizes,
    )
    fdr_bh = apply_fdr_correction(result, list(genes_by_pathway.keys()))

    return fdr_bh


def get_number_of_degs_by_pathway(top_genes, gmt_file):
    number_of_degs_by_pathway = np.array(
        [len(set(genes) & set(top_genes)) for pathway, genes in gmt_file.items()]
    )
    return number_of_degs_by_pathway


def run_hypergeom_tests(
    number_of_genes_in_deg,
    size_of_top_genes,
    number_of_degs_by_pathway,
    pathway_sizes,
):
    # TODO: pq esse menos 1? pypathway usa isso e eu n sei pq
    result = stats.hypergeom.sf(
        number_of_degs_by_pathway - 1,
        number_of_genes_in_deg,
        pathway_sizes,
        size_of_top_genes,
    )
    return result


def apply_fdr_correction(result, pathways):
    _, o, _, _ = multicomp.multipletests(result, method="fdr_bh")
    fdr = {pathways[i]: o[i] for i in range(len(result))}
    return pd.DataFrame.from_dict(data=fdr, orient="index", columns=["to_be_changed"])
