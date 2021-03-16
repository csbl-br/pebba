import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot


# init_notebook_mode(connected=True)


def create_interactive_plot(
    df, dict_genes_por_via, direction, analysis_name, statistics_for_plot, results_dir
):
    heatmap = create_heatmap(df)
    barplot = create_barplot(dict_genes_por_via, df)
    line_plot = plot_number_of_enriched_pathways(statistics_for_plot, direction)

    data = [heatmap, barplot, line_plot]

    layout = go.Layout(
        showlegend=False,
        # paper_bgcolor="rgba(0,0,0,0)",  # background color for the paper of the plot
        plot_bgcolor="rgb(255,255,255)",  # background color for the plot
        xaxis=dict(
            domain=[0, 0.80],
            showticklabels=False,
            mirror=True,
            showline=True,
            linewidth=1,
            linecolor="rgb(33, 27, 22)",
            matches="x2",
        ),
        xaxis2=dict(
            domain=[0, 0.80],
            showticklabels=False,
            mirror=True,
            showline=True,
            linewidth=1,
            linecolor="rgb(33, 27, 22)",
        ),
        xaxis3=dict(
            domain=[0.85, 1],
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
            "x": 0.8,
            "len": 0.7,
            "y": 1,
            "yanchor": "top",
            # "ticklabelposition": "inside top",
            # "tickfont": {"color": "rgb(0,0,0)"},
        },
    )
    return trace


def create_barplot(dict_genes_por_via, df):
    barplot = go.Bar(
        x=df.index.tolist(),
        y=[len(dict_genes_por_via[via]) for via in df.index.tolist()],
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


def plot_number_of_enriched_pathways(statistics_for_plot, direction):
    barplot = go.Bar(
        x=statistics_for_plot["times_significant_" + direction].tolist(),
        y=statistics_for_plot.index.tolist(),
        orientation="h",
        xaxis="x3",
        yaxis="y3",
        marker={
            "color":
            # "rgb(135, 57, 57)",
            "rgb(126, 139, 158)",
        },
    )
    return barplot
