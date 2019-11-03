import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

import plotly.graph_objects as go
import scipy.io

from torusmesh import TorusMesh
from funcpath import FuncPath

# read data
matdata = scipy.io.loadmat('dash_demo_data.mat')
x_t, y_t, z_t =  matdata['torus_data']
data = matdata['path_data']
end_points = matdata['end_points']

# generates base figure
fig = go.Figure(
    data=[go.Surface(x=x_t,
    y=y_t,
    z=z_t,
    opacity=0.50,
    name = 'base',
    showscale = False
    )],
    layout=go.Layout(
        uirevision=1,
        height = 500
		)
    )

# draw end points
fig.add_trace(go.Scatter3d(
    x=end_points[0],
    y=end_points[1],
    z=end_points[2],
    mode='markers',
    name='Vertices',
    marker = dict(symbol = 'x',  size = 4)
    ))

# draw initial path
fig.add_trace(go.Scatter3d(
    x=data[0,0],
    y=data[1,0],
    z=data[2,0],
    mode='lines',
    name='Path',
    ))


# external CSS stylesheets
external_stylesheets = [
    'https://codepen.io/chriddyp/pen/bWLwgP.css',
    {
        'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO',
        'crossorigin': 'anonymous'
    }
]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


server = app.server

app.layout = html.Div(children=[
    html.H1(children='Functional Iteration on Paths on the Torus'),

    # html.Div(children='''
    #     Particle dynamics on the Torus
    # '''),

    # html.Button(id='iterate_button', n_clicks=0, children='Iterate'),

    dcc.Graph(
        id='torus',
		figure = fig
	),


    html.H3(id='show-iteration'),

    dcc.Slider(
        id='iteration-slider',
        min=0,
        max=60,
        step=1,
        value=0,
        marks={i * 5: '{}'.format(i * 5) for i in range(13)},
    )
],
style={"max-width": "800px", "margin": "auto"},)


@app.callback(
    [Output('torus', 'figure'),
    Output('show-iteration','children')],
    [Input('iteration-slider', 'value')],
    [State('torus', 'figure')])
def update_figure(value, fig):
    # note: graph_objects become dict when passed as an argument into callback function
    new = {
        'mode': 'lines',
        'name': 'lines',
        'x': data[0,value],
        'y': data[1,value],
        'z': data[2,value],
        'type': 'scatter3d'}
    fig['data'][2] = new
    return fig, '{}-th iteration'.format(value)

if __name__ == '__main__':
    app.run_server(debug=True)

# TODO: plot end points
# TODO: host page
# TODO: hide colorbar
# TODO: user-input end points
# TODO: continuously iterate (interval?)
# TODO: make lines more visible