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
    # by pathway, genes member of the pathway
    genes_by_pathway = {
        pathway: list(set(genes) & set([str(x) for x in all_genes]))
        for pathway, genes in gmt_file.items()
    }
    return genes_by_pathway


def generate_differenttially_expressed_genes_by_pathway(top_genes, gmt_file):
    # by pathway, genes that are differentially expressed
    differenttially_expressed_genes_by_pathway = {
        pathway: list(set(genes) & set([str(x) for x in top_genes]))
        for pathway, genes in gmt_file.items()
    }
    return differenttially_expressed_genes_by_pathway


def run_hipergeometric_tests(
    all_genes, top_genes, genes_by_pathway, differenttially_expressed_genes_by_pathway
):
    # TODO: calcular esses lens antes
    # TODO: ver se stats.hypergeom.sf ja n aceita os negocios todos em paralelo
    result = {
        pathway: stats.hypergeom.sf(
            len(differenttially_expressed_genes_by_pathway[pathway]) - 1,
            len(all_genes),
            len(genes_by_pathway[pathway]),
            len(top_genes),
        )
        for pathway, genes in genes_by_pathway.items()
    }
    return result


def apply_fdr_correction(result):
    _, o, _, _ = multicomp.multipletests(list(result.values()), method="fdr_bh")
    fdr = {list(result.keys())[i]: o[i] for i in range(len(list(result.keys())))}
    return pd.DataFrame(fdr)
