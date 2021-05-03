import argparse
from pebba.main import pebba


parser = argparse.ArgumentParser(description="A tool for ORA meta-analysis")

parser.add_argument("deg", type=str, help="path to deg file")
parser.add_argument("gmt", type=str, help="path to gmt file")

parser.add_argument(
    "--gene",
    dest="gene_col",
    type=str,
    default="Gene.symbol",
    help="genes column name on the deg file",
)
parser.add_argument(
    "--logfc",
    dest="logFC_col",
    type=str,
    default="logFC",
    help="logFC column name on the deg file",
)
parser.add_argument(
    "--pval",
    dest="pvalue_col",
    type=str,
    default="P.Value",
    help="p_value column name on the deg file",
)
parser.add_argument(
    "--min",
    dest="min_genes",
    type=int,
    default=100,
    help="minimum number of genes to be considered in the analysis",
)
parser.add_argument(
    "--max",
    dest="max_genes",
    type=int,
    default=1500,
    help="maximum number of genes to be considered in the analysis",
)
parser.add_argument(
    "--p_cut",
    type=float,
    default=0.2,
    help="TODO",
)
parser.add_argument(
    "--drop_cut",
    type=float,
    default=0.05,
    help="minimum value of enrichment for a pathway to be included on the heatmap",
)
parser.add_argument(
    "--name",
    dest="analysis_name",
    type=str,
    help="name under which the analysis will be saved",
)
parser.add_argument(
    "-o",
    "--output_dir",
    dest="results_dir",
    type=str,
    default="Results",
    help="name of the directory in which to save the analysis",
)
parser.add_argument(
    "--no_plot",
    help="do not plot the heatmaps",
    dest="make_plots",
    action="store_false",
)
parser.add_argument(
    "--no_csv",
    help="do not return the dataframes",
    dest="save_csv",
    action="store_false",
)
parser.add_argument(
    "-f",
    "--force",
    help="if an analysis with the same name already exists, force it to be overwritten",
    action="store_true",
)


args = vars(parser.parse_args())
pebba(**args)
