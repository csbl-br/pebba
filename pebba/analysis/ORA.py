from scipy import stats
from statsmodels.sandbox.stats import multicomp
import pandas as pd


def run_ORA(top_genes, all_genes, gmt_file):
    genes_by_pathway = generate_genes_by_pathway(all_genes, gmt_file)

    differenttially_expressed_genes_by_pathway = generate_differenttially_expressed_genes_by_pathway(
        top_genes, gmt_file
    )
    result = run_hipergeometric_tests(
        all_genes,
        top_genes,
        genes_by_pathway,
        differenttially_expressed_genes_by_pathway,
    )
    fdr_bh = apply_fdr_correction(result)

    return fdr_bh


def generate_genes_by_pathway(all_genes, gmt_file):
    """
    key = pathways
    value = list of genes members of that pathway
    """

    # TODO: make the intersection between genes in the deg and the gmt before, outside and only once
    genes_by_pathway = {
        pathway: list(set(genes) & set(all_genes))
        for pathway, genes in gmt_file.items()
    }
    return genes_by_pathway


def generate_differenttially_expressed_genes_by_pathway(top_genes, gmt_file):
    """
    key = pathway
    value =  list of genes members of that pathway and reported as differentially expressed

    """
    differenttially_expressed_genes_by_pathway = {
        pathway: list(set(genes) & set(top_genes))
        for pathway, genes in gmt_file.items()
    }
    return differenttially_expressed_genes_by_pathway


def run_hipergeometric_tests(
    all_genes, top_genes, genes_by_pathway, differenttially_expressed_genes_by_pathway
):
    len_all_genes = len(all_genes)
    len_top_genes = len(top_genes)
    # TODO: ver se stats.hypergeom.sf ja n aceita os negocios todos em paralelo
    result = {
        pathway: stats.hypergeom.sf(
            len(differenttially_expressed_genes_by_pathway[pathway])
            - 1,  # TODO: pq esse menos 1? pypathway usa isso e eu n sei pq
            len_all_genes,
            len(genes_by_pathway[pathway]),
            len_top_genes,
        )
        for pathway, genes in genes_by_pathway.items()
    }
    return result


def apply_fdr_correction(result):
    _, o, _, _ = multicomp.multipletests(list(result.values()), method="fdr_bh")
    fdr = {list(result.keys())[i]: o[i] for i in range(len(list(result.keys())))}
    return pd.DataFrame.from_dict(data=fdr, orient="index", columns=["to_be_changed"])
