import os
import sys
import pandas as pd
import numpy as np # checar se uso, n lembro mais


#from pypathway import GMTUtils



#from pypathway import ORA




def pebba(file_in, 
          gmt_file, 
          gene_col="Gene.symbol",
          logFC_col="logFC",
          pvalue_col="P.Value",
          min_genes=100,
          max_genes=1500,
          p_cut=0.2,
          verbose=True,
          analysis_name= None, 
          results_dir="Results",
          force=False):

    validates_inputs(min_genes,max_genes,p_cut)

    create_results_directory(results_dir,force)

    

    
    ## Get information from all unique terms
    
    term2gene , path_desc , merge_p  = read_gmt_hier(gmt_file) #utils.read_gmt_hier(gmt_file)    
   
    ############################################
    f = lambda df: [gene for gene in df["gene"] ]
    dict_genes_por_via = dict(term2gene.groupby(["term"]).apply(f) )
    ################################################
    
# alterei funcao pra retornar o merge_p, 
# ta retornando no formato de vetor mas talvez tenha q mudar para dataframe ou serie
    
    
    
    if( isinstance(file_in , str)):
        deg_list = pd.read_csv(file_in,  sep = "\t")#header true no original, so tirei o header pra o python inferir
        if(analysis_name is None):
            analysis_name =   os.path.splitext(os.path.basename(file_in) )[0] # pega o basename e tira a extensao
     
    elif(isinstance(file_in , pd.DataFrame)):
        deg_list = file_in
    

    if(analysis_name is None):
        analysis_name = "PEBBA_analysis"
    
    ## Remove rows that do not have a valid gene symbol
    deg_list = deg_list.dropna()
    ## Get background genes as a character vector
    ## Empty values (non-annotated genes) will be removed
    all_genes = deg_list["Gene.symbol"] # passar pra lista?
    
       # Get cutoff values -------------------------------------------------------
    if(verbose):
        print("Getting cutoff")
        
    table_cut = _get_cutoff(deg_list, logFC_col, pvalue_col, min_genes, max_genes)
    
    directions = ["up", "down" , "any"]
    
    cut_paths =dict()
    dfs = dict()
    paths = dict()
    for direction in directions:
        if (verbose):
            print(direction + "\nGetting Pathways")
        dfs[direction] , paths[direction] = _get_pathway(merge_p, term2gene, all_genes,
                            deg_list, gene_col, logFC_col,
                            pvalue_col, direction,
                            min_genes, max_genes, p_cut)
        if (verbose):
            print("Getting Pathway Cutoff")
        
        cut_paths[direction] = _cutoff_path(paths[direction], p_cut, direction)
        
        
    
    
    ## Save heatmaps
    if(verbose):
        print("Saving heatmaps")    
        
    for direction in directions:
        create_interactive_plot(paths[direction], dict_genes_por_via , direction, analysis_name, cut_paths[direction], results_dir) 
    
    if(verbose): 
        print("Done")
        
        
        
 



def validates_inputs(min_genes,max_genes,p_cut):

    if(min_genes < 50 or min_genes > 2900):
        sys.exit("Variable min_genes must be between 50 and 2900 genes")
          
    if(max_genes < 100 or max_genes > 3000):
        sys.exit("Variable max_genes must be between 100 and 3000 genes")  
    
    if(p_cut < 0.00001 or p_cut > 1):
        sys.exit("Variable p_cut must be between 0.00001 and 1")       







def create_results_directory(results_dir,force):
    results_dir = os.path.abspath(results_dir)
    if not os.path.exists(results_dir): 
        os.makedirs("Results/Tables") 
        os.makedirs("Results/Heatmaps")      
    else:    
        if( not force): 
            sys.exit("Stopping analysis: ", results_dir, " already exists! Use force=True to overwrite.")
