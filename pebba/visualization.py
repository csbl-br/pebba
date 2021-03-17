import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

from pebba.analysis.auxiliary_analysis import calculate_how_many_above_cut


def create_interactive_plot(
    df,
    dict_genes_por_via,
    direction,
    analysis_name,
    results_dir,
    output_type="file",  # or div
    score_cut=1.0,  # TODO choose a sensible default for the score_cut and allow the user to input it
):
    heatmap = create_heatmap(df)
    barplot1 = create_barplot_pathway_counts(df, score_cut)
    barplot2 = create_barplot_genescut_count(df, score_cut)

    data = [heatmap, barplot1, barplot2]

    layout = go.Layout(
        showlegend=False,
        # paper_bgcolor="rgba(0,0,0,0)",  # background color for the paper of the plot
        plot_bgcolor="rgb(255,255,255)",  # background color for the plot
        xaxis=dict(
            domain=[0.18, 1],
            showticklabels=False,
            mirror=True,
            showline=True,
            linewidth=1,
            linecolor="rgb(33, 27, 22)",
            matches="x2",
        ),
        xaxis2=dict(
            domain=[0.18, 1],
            showticklabels=False,
            mirror=True,
            showline=True,
            linewidth=1,
            linecolor="rgb(33, 27, 22)",
        ),
        xaxis3=dict(
            domain=[0, 0.15],
            anchor="y3",
            mirror=True,
            showline=True,
            linewidth=1,
            linecolor="rgb(33, 27, 22)",
            gridcolor="rgb(200, 200, 200)",
            autorange="reversed",
        ),
        yaxis=dict(
            domain=[0.30, 1],
            mirror=True,
            showline=True,
            linewidth=1,
            linecolor="rgb(33, 27, 22)",
            matches="y3",
        ),
        yaxis2=dict(
            domain=[0, 0.25],
            autorange="reversed",
            anchor="x2",
            mirror=True,
            showline=True,
            linewidth=1,
            linecolor="rgb(33, 27, 22)",
            gridcolor="rgb(200, 200, 200)",
        ),
        yaxis3=dict(
            domain=[0.30, 1],
            anchor="x3",
            mirror=True,
            showline=True,
            linewidth=1,
            linecolor="rgb(33, 27, 22)",
        ),
    )
    figure = go.Figure(data=data, layout=layout)

    plot(
        figure,
        filename=results_dir + "/Heatmaps/" + analysis_name + "_" + direction + ".html",
        output_type=output_type,
    )


def create_heatmap(df):
    NGs = df.columns.tolist()
    pathways = df.index.tolist()
    values = [df[column].tolist() for column in df]

    trace = go.Heatmap(
        z=values,
        y=NGs,
        x=pathways,
        colorscale=
        # ["rgb(255,255,255)", "rgb(229, 45, 39)", "rgb(179, 18, 23)"], # red
        ["rgb(255,255,255)", "rgb(47, 187, 237)", "rgb(41, 128, 185)"],  # blue
        colorbar={
            "len": 0.7,
            "y": 1,
            "yanchor": "top",
            # "ticklabelposition": "inside top",
            # "tickfont": {"color": "rgb(0,0,0)"},
        },
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
            "color":
            # "rgb(135, 57, 57)", #red
            "rgb(126, 139, 158)",  # blue
        },
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
            "color":
            # "rgb(135, 57, 57)", #red
            "rgb(126, 139, 158)",  # blue
        },
    )
    return barplot
