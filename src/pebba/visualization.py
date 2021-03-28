import numpy as np
import plotly.graph_objs as go
from plotly.offline import plot

from pebba.analysis.auxiliary_analysis import calculate_how_many_above_cut


def create_interactive_plot(
    df,
    dict_genes_por_via,
    direction,
    analysis_name,
    results_dir,
    p_cut,
    output_type="file",  # or div
):
    heatmap = create_heatmap(df)

    path_cut_p = np.log10(p_cut) * (-1)
    barplot1 = create_barplot_pathway_counts(df, path_cut_p)
    barplot2 = create_barplot_genescut_count(df, path_cut_p)

    data = [heatmap, barplot1, barplot2]

    layout = generate_layout()
    figure = go.Figure(data=data, layout=layout)

    plot(
        figure,
        filename=results_dir + "/Heatmaps/" + analysis_name + "_" + direction + ".html",
        output_type=output_type,
        config={
            "displaylogo": False,
            "modeBarButtonsToRemove": ["pan2d", "toggleSpikelines"],
        },
    )


def generate_layout():

    # this dict manipulation will only work on python>=3.5
    # https://stackoverflow.com/questions/38987/how-do-i-merge-two-dictionaries-in-a-single-expression-taking-union-of-dictiona
    base_layout = dict(
        mirror=True,
        showline=True,
        linewidth=1,
        linecolor="rgb(33, 27, 22)",
    )

    layout = go.Layout(
        # paper_bgcolor="rgba(255,255,255,0)",  # background color for the paper of the plot
        # plot_bgcolor="rgb(255,255,255)",  # background color for the plot
        # hoverlabel=dict(
        #     # bgcolor="black", #too black for my taste #TODO find a good color profile
        # ),
        template="plotly_white",
        showlegend=False,
        xaxis={
            **base_layout,
            "domain": [0.18, 1],
            "matches": "x2",
            "showticklabels": False,
        },
        xaxis2={
            **base_layout,
            "domain": [0.18, 1],
            "showticklabels": False,
            "anchor": "y2",
        },
        xaxis3={
            **base_layout,
            "domain": [0, 0.15],
            "anchor": "y3",
            "autorange": "reversed",
        },
        yaxis={
            **base_layout,
            "domain": [0.30, 1],
            "matches": "y3",
        },
        yaxis2={
            **base_layout,
            "domain": [0, 0.25],
            "autorange": "reversed",
            "anchor": "x2",
        },
        yaxis3={
            **base_layout,
            "domain": [0.30, 1],
            "anchor": "x3",
        },
    )
    return layout


def create_heatmap(df):
    NGs = df.columns.tolist()
    pathways = df.index.tolist()
    values = [df[column].tolist() for column in df]

    trace = go.Heatmap(
        z=values,
        y=NGs,
        x=pathways,
        colorscale=["rgb(255,255,255)", "rgb(229, 45, 39)", "rgb(179, 18, 23)"],
        colorbar={
            "len": 0.7,
            "y": 1,
            "yanchor": "top",
        },
        hovertemplate="<b>Pathway: </b>%{x} <br>"
        + "<b>Gene cut: </b>%{y} <br>"
        + "<b>Enrichment Score: </b>%{z}",
        name="",
    )
    return trace


def create_barplot_pathway_counts(df, score_cut):
    barplot = go.Bar(
        x=df.index.tolist(),
        y=calculate_how_many_above_cut(df, path_cut_p=score_cut, axis_sum=1).tolist(),
        orientation="v",
        xaxis="x2",
        yaxis="y2",
        marker={
            "color": "rgb(135, 57, 57)",  # red
        },
        hovertemplate="<b>Pathway: </b>%{x} <br>"
        + "<b>Nº of times enrichment was detected: </b>%{y}",
        name="",
    )
    return barplot


def create_barplot_genescut_count(df, score_cut):
    barplot = go.Bar(
        x=calculate_how_many_above_cut(df, path_cut_p=score_cut, axis_sum=0).tolist(),
        y=df.columns.tolist(),
        orientation="h",
        xaxis="x3",
        yaxis="y3",
        marker={
            "color": "rgb(135, 57, 57)",  # red
        },
        # TODO pick a better name than gene cut and use it across the code
        hovertemplate="<b>Gene cut: </b>%{y} <br>"
        + "<b>Nº of times enrichment was detected: </b>%{x}",
        name="",
    )
    return barplot
