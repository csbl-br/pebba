from pebba.main import pebba


# TODO: colocar profiler_output.txt no gitignore e gitar tudo q tem aqui


# From the terminal:
# python -m cProfile -o profiler_output.txt profiling.py
# snakeviz profiler_output.txt
if __name__ == "__main__":
    pebba(
        "data/GSE49757_Septic_vs_Healthy.txt",
        "data/Reactome_2016_15and100Genes.gmt",
        force=True,
    )
