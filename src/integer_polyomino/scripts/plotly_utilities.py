import pandas as pd
import plotly.graph_objs as go
import seaborn as sns
import plotly
import itertools
import cufflinks as cf

def metric_to_piechart(metric):
    labels = ["irobustness", "evolvability", "unbound"]
    values = [metric[label] for label in labels]

    return go.Pie(labels = labels, values = values)

def row_to_pie(df, row):
    labels = ['irobustness', 'evolvability', 'robust_evolvability', 'unbound']
    values = [df.loc[row, label] for label in labels]

    values[0] -= values[2]/2
    values[1] -= values[2]/2

    return go.Pie(labels=labels, values=values)

def row_to_hbar(df, row):
    labels = ['irobustness', 'evolvability', 'robust_evolvability', 'unbound']
    values = [df.loc[row, label] for label in labels]

    values[0] -= values[2]/2
    values[1] -= values[2]/2

    return go.Bar(y=labels, x=values, orientation='h')

def pretty_genome_pie(df, row):

    trace = row_to_pie(df, row)
    layout = go.Layout({'title': 'Metric Distribution for ' + str(df.loc[row, "genome"]) + ' with phenotype ' + str(df.loc[row, 'pIDs'])})

    return go.Figure([trace], layout=layout)

def pretty_genome_hbar(df, row):

    trace = row_to_hbar(df, row)
    layout = go.Layout({'title': 'Metric Distribution for ' + str(df.loc[row, "genome"]) + ' with phenotype ' + str(df.loc[row, 'pIDs'])})

    return go.Figure([trace], layout=layout)


def pretty_set_pie(df, row):

    trace = row_to_pie(df, row)
    layout = go.Layout({'title': 'Metric Distribution for ' + str(df.loc[row, 'pIDs'])})

    return go.Figure([trace], layout=layout)


def metric_subplots(dff, selected_row_indices, parameters, barmode):
    metrics = ['srobustness', 'irobustness', 'evolvability',  'complex_diversity', 'diversity', 'robust_evolvability', 'complex_evolvability', 'rare', 'unbound']
    subplot_coord = list(itertools.product(range(1, 5), range(1, 3)))
    # genome_file = "GenomeMetrics_N{n_genes}_C{gen_colour}_T{threshold}_B{phenotype_builds}_Cx{high_colour}_J{n_jiggle}_D{determinism}_S{low_colour}".format(**parameters)
    titles = []
    for row_index in (selected_row_indices or []):
        for metric in metrics:
            titles.append('Isomorphic ' + metric + ' Distribution of pID set : ' + str(eval(dff['pIDs'][row_index])))
    fig = plotly.tools.make_subplots(
        rows=max(4, 4 * len(selected_row_indices)), cols=2,
        # subplot_titles=titles,
        shared_xaxes=False)
    fig_index = -4
    for row_index in (selected_row_indices or []):
        fig_index += 4
        pID = dff.loc[row_index, 'pIDs']
        df_genome = dff[dff['pIDs'] == pID]
        # new_colors = ['rgb' + str(rgb) for rgb in sns.color_palette("hls", df_genome.Iso_index.max() + 1)]
        # for iso_index in df_genome.Iso_index.unique():
        for metric, coord in zip(metrics, subplot_coord):
            fig.append_trace({
                'x': df_genome[metric],
                'type': 'histogram',
                # 'name': str(df_genome['original'].unique()[0]),
                # 'marker': dict(color=new_colors[iso_index]),
                'showlegend': True if metric == 'srobustness' else False
                }, coord[0] + fig_index, coord[1])
    for index, metric in zip(range(1, 9), metrics):
        fig['layout']['xaxis' + str(index)].update(title=metric)
        fig['layout']['yaxis' + str(index)].update(title='Number of genomes')
    fig['layout']['barmode'] = barmode
    fig['layout']['showlegend'] = True
    fig['layout']['height'] = max(fig_index, 1) * 400 * 4
    return fig

def distribution_metrics_phenotype_trace(pID, file_name, metrics, hdf):
    traces = []
    with pd.HDFStore(hdf, mode='r') as store:
        df = store.select(file_name, where='pIDs == pID')
    new_colors = ['rgb' + str(rgb) for rgb in sns.color_palette("hls", df.Iso_index.max() + 1)]
    for iso_index in df.Iso_index.unique():
        for metric in metrics:
            traces.append({
                'x': df[df['Iso_index'] == iso_index][metric],
                'type': 'histogram',
                'name': str(df[df['Iso_index'] == iso_index]['original'].unique()[0]),
                'marker': dict(color=new_colors[iso_index]),
                'legendgroup': str(iso_index),
                'showlegend': True if metric == 'srobustness' else False})
    return traces


def distribution_metric_phenotype_trace(pID, file_name, metric, hdf):
    traces = []
    with pd.HDFStore(hdf, mode='r') as store:
        df = store.select(file_name, where='pIDs == pID')
    new_colors = ['rgb' + str(rgb) for rgb in sns.color_palette("hls", df.Iso_index.max() + 1)]
    for iso_index in df.Iso_index.unique():
        traces.append({
            'x': df[df['Iso_index'] == iso_index][metric],
            'type': 'histogram',
            'name': str(df[df['Iso_index'] == iso_index]['original'].unique()[0]),
            'marker': dict(color=new_colors[iso_index]),
            'legendgroup': str(iso_index),
            'showlegend': True})
    return traces

def distribution_metric_phenotype_trace_dup(pID, file_names, metric, hdf):
    traces = []
    for file in file_names:
        with pd.HDFStore(hdf, mode='r') as store:
            df = store.select(file_name, where='pIDs == pID')
        new_colors = ['rgb' + str(rgb) for rgb in sns.color_palette("hls", df.Iso_index.max() + 1)]
        for iso_index in df.Iso_index.unique():
            traces.append({
                'x': df[df['Iso_index'] == iso_index][metric],
                'type': 'histogram',
                'name': str(df[df['Iso_index'] == iso_index]['original'].unique()[0]),
                'marker': dict(color=new_colors[iso_index]),
                'legendgroup': str(iso_index),
                'showlegend': True})
    return traces


def distribution_metric_phenotype_dup_trace(pID, genome_file, duplicate_file, metric, hdf):
    traces = []
    with pd.HDFStore(hdf, mode='r') as store:
        df = store.select(genome_file, where='pIDs == pID')
        dfd = store.select(duplicate_file, where='pIDs == pID')
    new_colors = ['rgb' + str(rgb) for rgb in sns.color_palette("hls", df.Iso_index.max() + dfd.Iso_index.max() + 2)]
    for iso_index in df.Iso_index.unique():
        traces.append({
            'x': df[df['Iso_index'] == iso_index][metric],
            'type': 'histogram',
            'name': str(df[df['Iso_index'] == iso_index]['original'].unique()[0]),
            'marker': dict(color=new_colors[iso_index]),
            'legendgroup': str(iso_index),
            'showlegend': True})
    for iso_index in dfd.Iso_index.unique():
        traces.append({
            'x': dfd[dfd['Iso_index'] == iso_index][metric],
            'type': 'histogram',
            'name': str(dfd[dfd['Iso_index'] == iso_index]['original'].unique()[0]),
            'marker': dict(color=new_colors[int(iso_index + df.Iso_index.max() + 1)]),
            'legendgroup': str(df.Iso_index.max() + iso_index + 1),
            'showlegend': True})
    return traces

def distribution_metrics_phenotype(pID, file_name, hdf, barmode='group'):
    fig = plotly.tools.make_subplots(rows=5, cols=2, shared_xaxes=False)
    subplot_coords = list(itertools.product(range(1, 6), range(1, 3)))
    metrics = ['srobustness', 'irobustness', 'evolvability',  'complex_diversity', 'diversity', 'robust_evolvability', 'complex_evolvability', 'rare', 'unbound']
    traces = distribution_metrics_phenotype_trace(pID=pID, file_name=file_name, metrics=metrics, hdf=hdf)


    for index in range(0, len(traces) // len(metrics)):
        for trace, coords in zip(traces[index * len(metrics):(index+1) * len(metrics)], subplot_coords):
            fig.append_trace(trace, coords[0], coords[1])

    # Layout settings
    fig['layout']['title'] = 'Metric Distributions for Phenotype : ' + pID
    for index, metric in zip(range(1, 10), metrics):
        fig['layout']['xaxis' + str(index)].update(title=metric)
        fig['layout']['yaxis' + str(index)].update(title='Number of genomes')
    fig['layout']['barmode'] = barmode
    fig['layout']['showlegend'] = True
    fig['layout']['height'] = 300 * 5
    return fig

def scatter_metric_size(df, metric, max_size = 10, deterministic=False, multi=False, **kwargs):
    mini_df = df[df['max_size'] <= max_size]

    if (deterministic or multi):
        mini_df = mini_df[mini_df['pIDs'].apply(lambda x: len(eval(x))) == 1]

    fig = mini_df.iplot(kind='scatter', mode='markers', x='max_size', y=metric, text='pIDs', asFigure=True, xTitle='max_size', yTitle=metric, colors='#1f77b4', **kwargs)

    if multi:
        mini_df = df[df['max_size'] <= max_size]
        mini_df = mini_df[mini_df['pIDs'].apply(lambda x: len(eval(x))) > 1]
        fig.add_trace(mini_df.iplot(kind='scatter', mode='markers', x='max_size', y=metric, text='pIDs', asFigure=True, colors='#ff7f0e', **kwargs)['data'][0])
        fig.data[0].name = 'deterministic'
        fig.data[1].name = 'multi'

    fig['layout']['hovermode'] = 'closest'

    return fig

def dual_metric_size(df, metrics, max_size = 10, **kwargs):
    mini_df = df[df['max_size'] <= max_size]

    fig = mini_df.iplot(kind='scatter', mode='markers', x='max_size', y=metrics, text='pIDs', asFigure=True,
                        xTitle='max_size', yTitle='Neighbourhood %', **kwargs)

    fig['layout']['hovermode'] = 'closest'

    return fig

def scatter_metrics(df, xMetric, yMetric, max_size = 10, deterministic=False, multi=False, **kwargs):
    mini_df = df[df['max_size'] <= max_size]

    if (deterministic or multi):
        mini_df = mini_df[mini_df['pIDs'].apply(lambda x: len(eval(x))) == 1]

    fig = mini_df.iplot(kind='scatter', mode='markers', x=xMetric, y=yMetric, text='pIDs', asFigure=True,
                        xTitle=xMetric, yTitle=yMetric, **kwargs)

    if multi:
        mini_df = df[df['max_size'] <= max_size]
        mini_df = mini_df[mini_df['pIDs'].apply(lambda x: len(eval(x))) > 1]
        fig.add_trace(mini_df.iplot(kind='scatter', mode='markers', x=xMetric, y=yMetric, text='pIDs', asFigure=True,
                            xTitle=xMetric, yTitle=yMetric, **kwargs)['data'][0])
        fig.data[0].name = 'deterministic'
        fig.data[1].name = 'multi'

    fig['layout']['hovermode'] = 'closest'

    return fig

def all_metrics_size(df, metrics, suffixe='_norm', max_size = 10, **kwargs):

    ytraces = [metric + suffixe for metric in metrics]
    mini_df = df[df['max_size'] < max_size]
    fig = mini_df.iplot(kind='scatter', subplots=True, shape=(5, 2),  x='max_size', y=ytraces, mode='markers', text='pIDs',
                        dimensions=(1600, 1600), asFigure=True, **kwargs)

    fig['layout']['hovermode'] = 'closest'

    return fig
