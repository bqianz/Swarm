import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

import plotly.graph_objects as go

from torusmesh import TorusMesh
from funcpath import FuncPath

n = 50
c, a = 2, 1

torus = TorusMesh(n,c,a)

start = (0,0)
end = (5,5)

paths = torus.initial_path(start,end)
num_paths = len(paths)

# generates base figure
fig = torus.plotly_draw_go()

# draw initial path

for i in range(num_paths):
    paths[i].plotly_draw_path_go(fig)

app = dash.Dash(__name__)

app.layout = html.Div(children=[
    html.H1(children='Functional Iteration on Paths on the Torus'),

    # html.Div(children='''
    #     Particle dynamics on the Torus
    # '''),

    html.Button(id='iterate_button', n_clicks=0, children='Iterate'),

    dcc.Graph(
        id='torus',
		figure = fig
	)

])


@app.callback(
    Output('torus', 'figure'),
    [Input('iterate_button', 'n_clicks')],
    [State('torus', 'figure')]
    )
def update_figure(n_clicks, fig):
    # note: graph_objects become dict when passed as an argument into callback function

    # in this order because callback is fired at n_clicks = 0
    for i in range(num_paths):
        paths[i].plotly_update_path(i,fig)
        paths[i].functional_iteration()
    return fig


if __name__ == '__main__':
    app.run_server(debug=True)

# TODO: plot end points
# TODO: host page
# TODO: hide colorbar
# TODO: user-input end points
# TODO: continuously iterate (interval?)
# TODO: make lines more visible