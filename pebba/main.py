import os
import sys

from pebba.utils.gmt_utils import get_gmt
from pebba.utils.deg_utils import get_deg
from pebba.analysis.statistics_of_ORA_exploration import generate_ORA_statistics
from pebba.analysis.ORA_hyperparameter_exploration import generate_ORA_dataframe
from pebba.visualization import create_interactive_plot


def pebba(
    deg,
    gmt,
    gene_col="Gene.symbol",
    logFC_col="logFC",
    pvalue_col="P.Value",
    min_genes=100,
    max_genes=1500,
    p_cut=0.2,
    analysis_name=None,
    results_dir="Results",
    force=False,
):
    validate_range_of_inputs(min_genes, max_genes, p_cut)

    create_results_directory(results_dir, force)
    if analysis_name is None:
        analysis_name = set_analysis_name(deg)

    deg = get_deg(deg, gene_col, logFC_col, pvalue_col)
    dict_of_genes_by_pathway = get_gmt(gmt, all_genes_in_deg=deg["Gene.symbol"])

    ORA_dataframe_all_directions = get_ORA_dataframes(
        deg, dict_of_genes_by_pathway, min_genes, max_genes
    )
    statistics_for_plot = get_statistics_for_plots(ORA_dataframe_all_directions, p_cut)

    create_interactive_plots(
        ORA_dataframe_all_directions,
        dict_of_genes_by_pathway,
        statistics_for_plot,
        analysis_name,
        results_dir,
    )

    return


def validate_range_of_inputs(min_genes, max_genes, p_cut):
    if min_genes < 50 or min_genes > 2900:
        raise ValueError("Variable min_genes must be between 50 and 2900 genes")

    if max_genes < 100 or max_genes > 3000:
        raise ValueError("Variable max_genes must be between 100 and 3000 genes")

    if p_cut < 0.00001 or p_cut > 1:
        raise ValueError("Variable p_cut must be between 0.00001 and 1")


def create_results_directory(results_dir, force):
    results_dir = os.path.abspath(results_dir)
    if not os.path.exists(results_dir):
        os.makedirs("Results/Heatmaps")
    else:
        if not force:
            sys.exit(
                "Stopping analysis: "
                + results_dir
                + " already exists! Use force=True to overwrite."
            )


def set_analysis_name(file_in):
    """
    If a name is not provided, generate one from the file name or give the generic name "PEBBA_analysis"
    :param file_in: string or Pandas DataFrame
    :return: string
    """
    if isinstance(file_in, str):
        return os.path.splitext(os.path.basename(file_in))[0]
    else:
        return "PEBBA_analysis"


####################################################################################
# TODO: On a second round of refactoring, change the program structure to get rid of this functions.
# The idea is to generate the ORA dataframe, then the statistics and then save the plot uninterupted for one direction at a time,
# instead of doing 3 times each step

# TODO: Substitute sys.exit calls by exceptions


def get_ORA_dataframes(deg, dict_genes_by_pathway, min_genes, max_genes):
    directions = ["up", "down", "any"]
    ORA_dataframe_all_directions = {
        direction: generate_ORA_dataframe(
            deg, dict_genes_by_pathway, direction, min_genes, max_genes
        )
        for direction in directions
    }

    return ORA_dataframe_all_directions


def get_statistics_for_plots(ORA_dataframe_all_directions, p_cut):
    directions = ["up", "down", "any"]
    statistics_for_plot = {
        direction: generate_ORA_statistics(
            ORA_dataframe_all_directions[direction], p_cut, direction
        )
        for direction in directions
    }

    return statistics_for_plot


def create_interactive_plots(
    ORA_dataframe_all_directions,
    dict_genes_by_pathway,
    statistics_for_plot,
    analysis_name,
    results_dir,
):
    directions = ["up", "down", "any"]
    for direction in directions:
        create_interactive_plot(
            ORA_dataframe_all_directions[direction],
            dict_genes_by_pathway,
            direction,
            analysis_name,
            statistics_for_plot[direction],
            results_dir,
        )
