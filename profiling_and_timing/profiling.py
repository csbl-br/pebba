import sys, os

dir_path = os.path.dirname(os.path.realpath(__file__))
package_path = os.path.split(dir_path)[0]
sys.path.append(package_path)


from pebba.main import pebba

print("ok")
# From the terminal:
# python -m cProfile -o profiler_output.txt profiling.py
# snakeviz profiler_output.txt
if __name__ == "__main__":
    pebba(
        "../data/GSE49757_Septic_vs_Healthy.txt",
        "../data/Reactome_2016_15and100Genes.gmt",
        force=True,
    )
