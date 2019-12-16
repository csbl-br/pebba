from scipy import stats
from statsmodels.sandbox.stats import multicomp
import pandas as pd


def run_ORA(top_genes, number_of_genes_in_deg, genes_by_pathway):
    size_of_top_genes = len(top_genes)
    differentially_expressed_genes_by_pathway = get_number_of_differentially_expressed_genes_by_pathway(
        top_genes, genes_by_pathway
    )

    result = run_hypergeom_tests(
        number_of_genes_in_deg,
        size_of_top_genes,
        genes_by_pathway,
        differentially_expressed_genes_by_pathway,
    )
    fdr_bh = apply_fdr_correction(result)

    return fdr_bh


def get_number_of_differentially_expressed_genes_by_pathway(top_genes, gmt_file):
    """
    key = pathway
    value =  list of genes members of that pathway and reported as differentially expressed

    """
    number_of_differentially_expressed_genes_by_pathway = {
        pathway: len(set(genes) & set(top_genes)) for pathway, genes in gmt_file.items()
    }
    return number_of_differentially_expressed_genes_by_pathway


def run_hypergeom_tests(
    number_of_genes_in_deg,
    size_of_top_genes,
    genes_by_pathway,
    number_of_differentially_expressed_genes_by_pathway,
):
    # TODO: ver se stats.hypergeom.sf ja n aceita os negocios todos em paralelo
    # DONE: aceita, possivel refatoracao util no futuro
    # TODO: pq esse menos 1? pypathway usa isso e eu n sei pq
    # TODO: pass len(genes_by_pathway) already calculated, do it once at the beggining
    result = {
        pathway: stats.hypergeom.sf(
            number_of_differentially_expressed_genes_by_pathway[pathway] - 1,
            number_of_genes_in_deg,
            len(genes_by_pathway[pathway]),
            size_of_top_genes,
        )
        for pathway, genes in genes_by_pathway.items()
    }
    return result


def apply_fdr_correction(result):
    _, o, _, _ = multicomp.multipletests(list(result.values()), method="fdr_bh")
    fdr = {list(result.keys())[i]: o[i] for i in range(len(list(result.keys())))}
    return pd.DataFrame.from_dict(data=fdr, orient="index", columns=["to_be_changed"])
