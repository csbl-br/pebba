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
    drop_cut,
    output_type="file",  # or div
):

    # drop every pathway that is practically not enriched at all
    df = df[df.apply(lambda row: not all(row < drop_cut), axis=1)]

    heatmap_color, bar_color = pick_colors(direction)
    heatmap = create_heatmap(df, heatmap_color)
    path_cut_p = np.log10(p_cut) * (-1)
    barplot1 = create_barplot_pathway_counts(df, path_cut_p, bar_color)
    barplot2 = create_barplot_genescut_count(df, path_cut_p, bar_color)

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


def pick_colors(direction):
    colorscales = {
        "up": ["rgb(255,255,255)", "rgb(229, 45, 39)", "rgb(179, 18, 23)"],  # red
        "down": ["rgb(255,255,255)", "rgb(47, 187, 237)", "rgb(41, 128, 185)"],  # blue
        "any": ["rgb(255,255,255)", "rgb(91, 91, 102)", "rgb(33, 33, 36)"],  # grey
        # "any": ["rgb(255,255,255)", "rgb(158, 39, 227)", "rgb(103, 17, 173)"],  # purple
    }
    colors = {
        "up": "rgb(135, 57, 57)",  # red
        "down": "rgb(126, 139, 158)",  # blue
        "any": "rgb(91, 91, 102)",  # grey
        # "any": "rgb(109, 10, 166)",  # purple
    }
    return colorscales[direction], colors[direction]


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


def create_heatmap(df, colorscale):
    NGs = df.columns.tolist()
    pathways = df.index.tolist()
    values = [df[column].tolist() for column in df]

    trace = go.Heatmap(
        z=values,
        y=NGs,
        x=pathways,
        colorscale=colorscale,
        colorbar={
            "len": 0.7,
            "y": 1,
            "yanchor": "top",
        },
        hovertemplate="<b>Pathway: </b>%{x} <br>"
        + "<b>Nº of genes considered: </b>%{y} <br>"
        + "<b>Enrichment Confidence Score: </b>%{z}",
        name="",
    )
    return trace


def create_barplot_pathway_counts(df, score_cut, color):
    barplot = go.Bar(
        x=df.index.tolist(),
        y=calculate_how_many_above_cut(df, path_cut_p=score_cut, axis_sum=1).tolist(),
        orientation="v",
        xaxis="x2",
        yaxis="y2",
        marker={"color": color},
        hovertemplate="<b>Pathway: </b>%{x} <br>"
        + "<b>Nº of times enrichment was detected: </b>%{y}",
        name="",
    )
    return barplot


def create_barplot_genescut_count(df, score_cut, color):
    barplot = go.Bar(
        x=calculate_how_many_above_cut(df, path_cut_p=score_cut, axis_sum=0).tolist(),
        y=df.columns.tolist(),
        orientation="h",
        xaxis="x3",
        yaxis="y3",
        marker={"color": color},
        hovertemplate="<b>Nº of genes considered: </b>%{y} <br>"
        + "<b>Nº of times enrichment was detected: </b>%{x}",
        name="",
    )
    return barplot
