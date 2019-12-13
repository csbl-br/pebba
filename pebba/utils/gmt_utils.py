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
