import pandas as pd


def read_gmt(file_name):
    """
    file_name string -> dataframe
    Reads .gmt file and returns a Pandas DataFrame
    """
    pathway_names = []
    gene_names = []

    with open(file_name, "r") as file:
        # separar cada elemento separado por tab e guardar eles
        for line in file:
            pathway_names.append(line.split("\t")[0])
            gene_names.append(line.split("\t")[2:])
    temp = []
    for i in range(len(gene_names)):
        # apagar \n presente no ultimo gene de cada lista (artefato da leitura do arquivo)
        gene_names[i][-1] = gene_names[i][-1].replace("\n", "")
        temp.append(
            pd.DataFrame(
                {"term": [pathway_names[i]] * len(gene_names[i]), "gene": gene_names[i]}
            )
        )
        # TODO: checar se posso só tirar essa declaracao de pd.DataFrame e só dps no final dar um concatzao em tudo direto nos dicionarios

    gene_to_pathway_relationship_df = pd.concat(temp)
    gene_to_pathway_relationship_df = gene_to_pathway_relationship_df.reset_index(
        drop=True
    )

    # TODO: pass the utils gmt reader into the run_ORA function to improve modularity
    # TODO: gene_to_pathway_relationship_df and dict_genes_by_pathway contain the same information, transform into one thing only
    f = lambda df: [gene for gene in df["gene"]]
    dict_genes_by_pathway = dict(
        gene_to_pathway_relationship_df.groupby(["term"]).apply(f)
    )

    return dict_genes_by_pathway


# TODO: check if this implementation by pypathway is really the same as I did, if so it will spare me time refactoring. See license and how to cite them.
# def parse_gmt_file(file):
#     '''
#     parse a local gmt file,
#     the file should be presented like:
#     setName\tsource[optional]\tgenes....
#
#     :param file: the file path
#     :return: the parsed dict
#     '''
#     with open(file) as fp:
#         con = fp.read()
#     return {x.split('\t')[0]: [t for t in x.split('\t')[2:] if t] for x in con.split('\n') if x}
