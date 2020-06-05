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
        xaxis=dict(domain=[0, 0.80], showticklabels=False),
        xaxis2=dict(domain=[0, 0.80], showticklabels=False),
        xaxis3=dict(domain=[0.85, 1], anchor="y3"),
        yaxis=dict(domain=[0.30, 1]),
        yaxis2=dict(domain=[0, 0.25], autorange="reversed", anchor="x2"),
        yaxis3=dict(domain=[0.30, 1], anchor="x3"),
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

    trace = go.Heatmap(z=values, y=NGs, x=pathways)
    return trace


def create_barplot(dict_genes_por_via, df):
    barplot = go.Bar(
        x=df.index.tolist(),
        y=[len(dict_genes_por_via[via]) for via in df.index.tolist()],
        orientation="v",
        xaxis="x2",
        yaxis="y2",
    )
    return barplot


def plot_number_of_enriched_pathways(statistics_for_plot, direction):
    plot = go.Scatter(
        x=statistics_for_plot["times_significant_" + direction].tolist(),
        y=statistics_for_plot.index.tolist(),
        orientation="h",
        mode="lines+markers",
        name="lines+markers",
        xaxis="x3",
        yaxis="y3",
    )
    return plot