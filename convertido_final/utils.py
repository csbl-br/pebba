import pandas as pd

def read_gmt_hier(file_name):
    '''
    file_name string -> dataframe
    Reads .gmt file and returns a Pandas DataFrame
    '''
    pathway_names = []
    gene_names = []

    with open(file_name, 'r') as file:
        # separar cada elemento separado por tab e guardar eles
        for line in file:
            pathway_names.append(line.split("\t")[0])
            gene_names.append(line.split("\t")[2:])
    temp = []
    for i in range(len(gene_names)):
        # apagar \n presente no ultimo gene de cada lista (artefato da leitura do arquivo)
        gene_names[i][-1] = gene_names[i][-1].replace("\n", "")
        temp.append(pd.DataFrame({'term': [pathway_names[i]] * len(gene_names[i]),
                                  'gene': gene_names[i]}))
        # TODO: checar se posso só tirar essa declaracao de pd.DataFrame e só dps no final dar um concatzao em tudo direto nos dicionarios

    res = pd.concat(temp)
    res = res.reset_index(drop=True)
    return res
