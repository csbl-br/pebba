import numpy as np


def get_top_genes(deg, direction, max_genes):
    top_genes = sort_deg_by_direction(direction, deg, max_genes)
    top_genes["pi_value"] = calculate_pi_value(top_genes)
    top_genes = sort_top_genes_by_pi_value(top_genes)
    return top_genes


def sort_deg_by_direction(direction, deg_list, max_genes):
    if direction == "up":
        top_genes = deg_list.sort_values(by="logFC", ascending=False).head(n=max_genes)
    elif direction == "down":
        top_genes = deg_list.sort_values(by="logFC", ascending=True).head(n=max_genes)
    elif direction == "any":
        deg_list["logFC"] = deg_list["logFC"].astype(np.float64)
        deg_list["logFC"] = deg_list["logFC"].abs()
        top_genes = deg_list.sort_values(by="logFC", ascending=True).head(n=max_genes)
    else:
        raise ValueError("Invalid direction argument")
    top_genes["Gene.symbol"] = top_genes["Gene.symbol"].astype(str)
    return top_genes


def calculate_pi_value(top_genes):
    return top_genes["logFC"].apply(abs) * (-top_genes["P.Value"].apply(np.log10))


def sort_top_genes_by_pi_value(top_genes):
    return top_genes.sort_values(by="pi_value", ascending=False).reset_index(drop=True)
