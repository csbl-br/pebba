import numpy as np
import pandas as pd
from pypathway import ORA, GMTUtils





def _get_cutoff(deg_list, logFC_col, pvalue_col, min_genes, max_genes):
    
    dirs = ["down", "up"]
    
    res_up = _get_directional_cutoff("up",deg_list, logFC_col, pvalue_col, min_genes, max_genes)
    res_down = _get_directional_cutoff("down",deg_list, logFC_col, pvalue_col, min_genes, max_genes)
    
    res = res_down.merge(res_up,on ="TopCut", how="outer" , suffixes = ("_down","_up"))
    
    
    res["minimum_log2fc_combined"] = res[['minimum_log2fc_down','minimum_log2fc_up']].min(axis=1)
    res["minimum_MinuslogP_combined"] = res[['minimum_MinuslogP_down','minimum_MinuslogP_up']].min(axis=1) 
    res["minimum_Pi_combined"] = res[['minimum_Pi_down','minimum_Pi_up']].min(axis=1)
    
    
    return res


def _get_directional_cutoff(direction,deg_list, logFC_col, pvalue_col, min_genes, max_genes):
    if (direction == "down"):
        ascending = True # no original, decreasing = False
    else:
        ascending = False
            
    #pega o deg_list e ordena de maneira decrescente ou crescente usando o logFC_col como chave
    # ai ele pega só as n=max genes primeiras linhas e retorna os valores de logFC
    top = deg_list.sort_values(by = logFC_col, ascending = ascending)[[logFC_col , pvalue_col]]
    top = top[:max_genes]
    
    top["pi_value"] = top[logFC_col].apply(abs) * (- top[pvalue_col].apply(np.log10))
    top = top.sort_values(by = "pi_value" , ascending = False)
    df = pd.DataFrame(columns = ["minimum_log2fc","minimum_MinuslogP","minimum_Pi", "TopCut"])
    rows = []
    for i in range(min_genes, max_genes,50):
        top_genes = top.iloc[0:i]
        minFC = min(abs(top_genes[logFC_col]))
        maxP = max(top_genes[pvalue_col])
        minP = - np.log10(maxP)
        
        minPi = min(top_genes["pi_value"])
        
        ##### ja foi ordenado, entao posso so pegar o elemento especifico ao invez de procurar o elemento de novo
 #       minPi = min(top_genes.iloc[i,3])##### essa passagem n faz sentido

        row = {"minimum_log2fc":minFC, "minimum_MinuslogP": minP , "minimum_Pi":minPi,"TopCut":i}
        rows.append(row)
    df = pd.DataFrame(rows)
    df.set_index("TopCut")
    return df



###############################################################################################################################################


def _cutoff_path(path_table, p_cut, direction):
    
    df_index = path_table.columns
    df = pd.DataFrame()
    df["MaxR"] = path_table.max()
    df["SumR"] = path_table.sum()
    path_cut_p = np.log10(p_cut) * (-1)
    
    
    #How many pathways above path_cut_p (freq)    
    how_many_pathways_above_cut = calculate_how_many_pathways_above_cut(path_table,path_cut_p, axis =0)
    n_rows =len(path_table.index) ######
    df["times"] = how_many_pathways_above_cut / n_rows
    
    df.columns = ["maximum_MinuslogP_"+ direction,
                  "sum_MinuslogP_"+ direction,
                  "times_significant_"+ direction]
    return df


###############################################################################################################################################


def _get_pathway(merge_p, term2gene, all_genes, deg_list,gene_col, logFC_col, pvalue_col, direction, min_genes, max_genes, p_cut):
   
    top = get_top(direction, deg_list, max_genes, logFC_col)
    
    top["pi_value"] = top[logFC_col].apply(abs) * (- top[pvalue_col].apply(np.log10))
    top = top.sort_values(by = "pi_value", ascending =False).reset_index(drop = True)
    
      
    pathGs = []  #melhorar isso
    for i in range(min_genes , max_genes , 50):
        top_genes = top.loc[0:i,gene_col].astype(str)
        pathG = _run_enrich(top_genes, all_genes, gmt_file)
        pathG.columns = [ "term" , str(i)]
        pathG = pathG.set_index("term",drop=True)
        pathGs.append(pathG)
    merge_p = pd.concat(pathGs, axis=1, join = "outer")
    merge_p.fillna(1.0) # acho q n eh mais necessario, ORA faz sozinho
    
    merge_p2 = ( merge_p.apply(np.log10) )*(-1)
    
    path_cut_p = np.log10(p_cut)*(-1)
    
    df = summarizes_ORA_information(merge_p2, path_cut_p,direction)
    merge_p2 = pd.concat([df, merge_p2], axis=1)    

    
    merge_p2 = merge_p2.sort_values(by = "FirstTopCut_significant_" + direction , ascending = False)
    merge_p2 = merge_p2.drop(labels = ["TopCut_highestMinuslogP_" + direction ,
                  "maximum_MinuslogP_" + direction ,
                  "sum_MinuslogP_" + direction ,
                  "times_significant_" + direction  ,
                  "FirstTopCut_significant_" + direction , 
                  "PEBBA_score_" + direction], axis =1 )
    
    #refatorar toda essa nojeira legada
    
    return  df , merge_p2
    

def get_top(direction, deg_list, max_genes, logFC_col):
    
    if(direction == "up"):
        top = deg_list.sort_values(by = logFC_col, ascending = False).head(n=max_genes)
    elif(direction =="down"):
        top = deg_list.sort_values(by= logFC_col, ascending = True).head(n=max_genes)
    elif(direction =="any"):
        deg_list[logFC_col] = deg_list[logFC_col].astype(np.float64) 
        deg_list[logFC_col] = deg_list[logFC_col].abs()
        top = deg_list.sort_values(by = logFC_col, ascending = True).head(n=max_genes)
    else:
        sys.exit("Invalid direction argument")
    return top


def summarizes_ORA_information(merge_p2, path_cut_p, direction) :


    NG = merge_p2.idxmax(axis = 1) # O recorte de genes q apresentou o maior p valor possui NG genes
    NG = NG.astype(np.int64)
    p_max = merge_p2.max(axis = 1)
    p_sum = merge_p2.sum(axis = 1)
    
    num_columns_merge_p2 = merge_p2.shape[1]
    how_many_pathways_above_cut = calculate_how_many_pathways_above_cut(merge_p2,path_cut_p,axis =1)
     
    times = how_many_pathways_above_cut / num_columns_merge_p2
    
    ES3 = (1 - np.exp(- p_max) / (1 + (0.1 * np.sqrt(NG)) ) )
   
    first = merge_p2.apply(first_column_above_path_cut_p , axis = 1, path_cut_p=path_cut_p )
    first = first.apply(lambda x: merge_p2.columns[x] if x !=0 else 0 )
 
    dicionario = {"TopCut_highestMinuslogP_" + direction : NG,
                  "maximum_MinuslogP_" + direction : p_max ,
                  "sum_MinuslogP_" + direction : p_sum,
                  "times_significant_" + direction : times ,
                  "FirstTopCut_significant_" + direction : first, 
                  "PEBBA_score_" + direction : ES3}
    
    df = pd.DataFrame(dicionario)
    df["FirstTopCut_significant_" + direction] = df["FirstTopCut_significant_" + direction].astype(np.int64)
    
    return df


def _run_enrich(top_genes, all_genes, gmt_file):
    term2gene = GMTUtils.parse_gmt_file(gmt_file)
    df = ORA.run(top_genes, all_genes, term2gene).df
    df = df[["name", "fdr"]]
    return df   


def first_column_above_path_cut_p(row, path_cut_p):
    for cont , element in enumerate(row):
        if element > path_cut_p:
            return cont
    
    return 0


def calculate_how_many_pathways_above_cut(df, path_cut_p,axis):
    f = lambda x: x > path_cut_p
    how_many_pathways_above_cut =  df.apply(f,axis=1).sum(axis=axis)
    return how_many_pathways_above_cut   

###############################################################################################################################################
# ## backup

# In[ ]:


#
# def calculate_how_many_pathways_above_cut(df, path_cut_p):
#     f = lambda x: x > path_cut_p
#     how_many_pathways_above_cut =  df.apply(f,axis=1).sum(axis=1)
#     return how_many_pathways_above_cut   

# # # \_cutoff_path
# # In[59]:

# def _cutoff_path(path_table, p_cut, direction):
    
#     df_index = path_table.columns
#     df = pd.DataFrame()
#     df["MaxR"] = path_table.max()
#     df["SumR"] = path_table.sum()
#     path_cut_p = np.log10(p_cut) * (-1)
    
    
#     #How many pathways above path_cut_p (freq)    
#     how_many_pathways_above_cut = calculate_how_many_pathways_above_cut(path_table,path_cut_p)
#     n_rows =len(path_table.index) ######
#     print(how_many_pathways_above_cut)
#     df["times"] = how_many_pathways_above_cut / n_rows
    
#     df.columns = ["maximum_MinuslogP_"+ direction,
#                   "sum_MinuslogP_"+ direction,
#                   "times_significant_"+ direction]
    
#     df.set_index(df_index) #drop True inutil (?)
# #     print(df)
#     return df




# def summarizes_ORA_information(merge_p2, path_cut_p, direction) :


#     NG = merge_p2.idxmax(axis = 1) # O recorte de genes q apresentou o maior p valor possui NG genes
#     NG = NG.astype(np.int64)
#     p_max = merge_p2.max(axis = 1)
#     p_sum = merge_p2.sum(axis = 1)
    
#     num_columns_merge_p2 = merge_p2.shape[1]
#     how_many_pathways_above_cut = calculate_how_many_pathways_above_cut(merge_p2,path_cut_p)
     
#     times = how_many_pathways_above_cut / num_columns_merge_p2
    
#     ES3 = (1 - np.exp(- p_max) / (1 + (0.1 * np.sqrt(NG)) ) )
   
#     first = merge_p2.apply(first_column_above_path_cut_p , axis = 1, path_cut_p=path_cut_p )
#     first = first.apply(lambda x: merge_p2.columns[x] if x !=0 else 0 )
 
#     dicionario = {"TopCut_highestMinuslogP_" + direction : NG,
#                   "maximum_MinuslogP_" + direction : p_max ,
#                   "sum_MinuslogP_" + direction : p_sum,
#                   "times_significant_" + direction : times ,
#                   "FirstTopCut_significant_" + direction : first, 
#                   "PEBBA_score_" + direction : ES3}
    
#     df = pd.DataFrame(dicionario)
#     df["FirstTopCut_significant_" + direction] = df["FirstTopCut_significant_" + direction].astype(np.int64)
    
#     return df

