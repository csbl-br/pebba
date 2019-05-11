import numpy as np
import pandas as pd


def _get_directional_cutoff(direction,deg_list, logFC_col, pvalue_col, min_genes, max_genes):
    if (direction == "down"):
        ascending = True # no original, decreasing = False
    else:
         ascending = False
            
    #pega o deg_list e ordena de maneira decrescente ou crescente usando o logFC_col como chave
    # ai ele pega só as n=max genes primeiras linhas e retorna os valores de logFC
    top = deg_list.sort(logFC_col, ascending = ascending)[[logFC_col , pvalue_col]]
    top = top[:max_genes]
    
    top["pi_value"] = top[logFC_col].apply(abs) * (- top[pvalue_col].apply(np.log10))
    top = top.sort(columns = "pi_value" , ascending = False)
    df1 = pd.DataFrame(columns = ["minFC","minP","minPi", "topcut"])
    for i in range(min_genes, max_genes,50):
        top_genes = top.ix[0:i]
        minFC = min(abs(top_genes[logFc]))
        maxP = max(top_genes[pvalue_col])
        minP = - np.log10(maxP)
        minPi = min(top_genes.ix[i,3])##### essa passagem n faz sentido
        #topcut=i
        row_X = {"minFC":minFC, "minP": minP , "minPi":minPi,"topcut":i}
        df1.concat(row_X)
    df1.set_index("topcut")
    return df1



def _get_cutoff(deg_list, logFC_col, pvalue_col, min_genes, max_genes):
    dirs = ["down", "up"]
    res_up = get_directional_cutoff("up",deg_list, logFC_col, pvalue_col, min_genes, max_genes)
    res_down = get_directional_cutoff("down",deg_list, logFC_col, pvalue_col, min_genes, max_genes)
    res = res_down.merge(resp_up,on ="topcut", how="outer" , suffixes = ("_down","_up"))
    
    res["fc"] = min(res["minFC_down"],res["minFC_up"])
    res["p"] = min(res["minP_down"],res["minP_up"])
    res["pi"] = min(res["minPi_down"],res["minPi_up"])
    
    #ver se isso vai funcionar, sepa vale a pena ja passar tudo com os nomes certos pra ter menos trabalho e mais limpo
    res.columns= ["TopCut", "minimum_log2fc_down", "minimum_MinuslogP_down",
                    "minimum_Pi_down", "minimum_log2fc_up", "minimum_MinuslogP_up",
                    "minimum_Pi_up", "minimum_log2fc_combined",
                    "minimum_MinuslogP_combined", "minimum_Pi_combined"]
    return res





def _cutoff_path(path_table, p_cut, direction):
    df = path_table
    n_rows =len(df.index)
    df["MaxR"] = df.max(axis=1)
    df["SumR"] = df.sum(axis=1)
    path_cut_p = np.log10(p_cut) * (-1)
    #How many pathways above path_cut_p (freq)
    
    f = lambda x: x > path_cut_p
    how_many_pathways_above_cut =  df.apply(f,axis=1).count()
    df["times"] = how_many_pathways_above_cut / n_rows
    
    df.columns = ["maximum_MinuslogP_"+ direction,
                  "sum_MinuslogP_"+ direction,
                  "times_significant_"+ direction]
    return df
