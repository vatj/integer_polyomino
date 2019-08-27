#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import itertools
from plotly.colors import DEFAULT_PLOTLY_COLORS


# In[ ]:


import plotly.graph_objs as go


# In[ ]:


def reshape_phenotypes(array):

    return array[0], array[1], str(tuple([array[0], array[1]])), array[2], array[3], array[4:].reshape(array[3], array[2]).astype(np.int8)


def import_phenotypes(file_str=str()):

    init = pd.read_csv(file_str, header=None, squeeze=False, converters={0 : lambda x: np.fromstring(x, dtype=np.int8, sep=' ')}).squeeze()
    data = pd.DataFrame()
    data['size'], data['subindex'], data['pIDs'], data['height'], data['width'], data['visual'] = list(zip(*init.apply(reshape_phenotypes)))
    data.set_index('pIDs', inplace=True)

    return data

def pID_to_shape(height=2, width=2, visual=np.array([[0,1],[2,1]]), **kwargs):
    color_index = 0
    layout = go.Layout(height=height * 200, width=width * 200,
                       xaxis=dict(range=[0, width], showgrid=False, zeroline=False, showline=False, ticks='', showticklabels=False),
                       yaxis=dict(range=[0, height], showgrid=False, zeroline=False, showline=False, ticks='', showticklabels=False),
                       hovermode='closest')
    fig = go.Figure(layout=layout)
    offset = 0.005
    for x, y in itertools.product(range(0, width), range(0, height)):
        if visual[x, y]:
            color_index = visual[x, y] - 1
            fig.add_trace(dict(x =[x + offset, x + 1 - offset, x + 1 - offset, x + offset, x + offset],
                     y=[y + offset, y + offset, y + 1 - offset, y + 1 - offset, y + offset],
                     mode='lines', showlegend=False, fill='toself', fillcolor=DEFAULT_PLOTLY_COLORS[color_index],
                     line=dict(width=3, color='black')))
    return fig

def pID_to_traces(height=2, width=2, visual=np.array([[0,1],[2,1]]), **kwargs):
    color_index = 0
    offset = 0.005
    traces = []
    for x, y in itertools.product(range(0, width), range(0, height)):
        if visual[x, y]:
            color_index = visual[x, y] - 1
            traces.append(dict(x =[x + offset, x + 1 - offset, x + 1 - offset, x + offset, x + offset],
                     y=[y + offset, y + offset, y + 1 - offset, y + 1 - offset, y + offset],
                     mode='lines', showlegend=False, fill='toself', fillcolor=DEFAULT_PLOTLY_COLORS[color_index],
                     line=dict(width=3, color='black')))
    return traces



# In[ ]:





# In[ ]:
