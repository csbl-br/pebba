def get_gmt(gmt_file, all_genes_in_deg):
    dict_of_genes_by_pathway = read_gmt(gmt_file)
    dict_of_genes_by_pathway = preprocess_gmt(
        dict_of_genes_by_pathway, all_genes_in_deg
    )
    return dict_of_genes_by_pathway


def read_gmt(file):
    """
    parse a local gmt file,
    the file should be presented like:
    setName\tsource[optional]\tgenes....

    :param file: the file path
    :return: the parsed dict
    """
    with open(file) as fp:
        file = fp.read()
    return {
        line.split("\t")[0]: [t for t in line.split("\t")[2:] if t]
        for line in file.split("\n")
        if line
    }


def preprocess_gmt(dict_of_genes_by_pathway, all_genes_in_deg):
    """
    Genes not present in the DEG will be excluded of the analysis,
    in essence shortening the size of a given pathway.
    This may skew the results, users be advised to check your inputs as we assume they are consistent.
    TODO: como dar credito pro pypathway e questão das licensas
    TODO: decidir se essa eh a politica q eu vou adotar mesmo, talvez adicionar mais um disclaimer em algum lugar
    TODO: pelo menos ter uma maneira de recuperar o q foi jogado fora pra deixar registrado
    """
    genes_by_pathway = {
        pathway: list(set(genes) & set(all_genes_in_deg))
        for pathway, genes in dict_of_genes_by_pathway.items()
    }
    return genes_by_pathway
