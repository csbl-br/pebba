import pandas as pd

def read_gmt_hier(file_name):
    '''
    file_name string -> dataframe and dictionary
    
    Reads .gmt file and returns a Pandas DataFrame and a dictionary with the information of the gmt file in a format ready to be used.
    
    '''
    gmt_names = []
    gmt_desc  = []
    gmt_genes = []
    res =pd.DataFrame()

    with open(file_name, 'r') as f:
        # separar cada elemento separado por tab e guardar eles
        for line in f:
            gmt_names.append(line.split("\t")[0])
            gmt_desc.append(line.split("\t")[1])
            gmt_genes.append(line.split("\t")[2:])

    for i in range(len(gmt_genes)) :

        # apagar \n presente no ultimo gene de cada lista (artefato da leitura do arquivo)
        gmt_genes[i][-1] = gmt_genes[i][-1].replace("\n", "")

        # Poem na forma de um dataframe, cada linha um gene e suas informações relativas 
        temp= pd.DataFrame({'term': [gmt_names[i]]*len(gmt_genes[i]), 'hier':  [gmt_desc[i]]*len(gmt_genes[i]), 'gene' : gmt_genes[i] })
        res = pd.concat([res,temp])

    # reseta o indice
    res = res.reset_index(drop=True)    

    #relação entre nomes e descricões  (é pra isso q essa variavel serve? no original é um datafra esquisito, achei q assim ia ser mais otimizado)   
    path_desc = dict(zip(gmt_names,gmt_desc))

    return res, path_desc , gmt_names
