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

p1, p2, p3 = torus.initial_path(start,end)

# generates base figure
fig = torus.plotly_figure_realtime()

app = dash.Dash(__name__)

app.layout = html.Div(children=[
    html.H1(children='Torus'),

    html.Div(children='''
        Particle dynamics on the Torus
    '''),

    # html.Div(id = 'check_type', children='''
    #     fig is of type "{}"
    # '''.format(type(fig))),

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
    p1.plotly_path_dict(fig)
    p1.functional_iteration()
    return fig


# TODO: make it so angle stays in user-input in layout
# TODO: show multiple lines

# @app.callback(
#     Output('check_type', 'children'),
#     [Input('iterate_button', 'n_clicks')],
#     [State('torus', 'figure')]
#     )
# def check_figure_type(n_clicks, fig):
#     return'''
#         fig is of type "{}"
#     '''.format(type(fig))

# @app.callback(
#     Output('check_type', 'children'),
#     [Input('iterate_button', 'n_clicks')],
#     [State('torus', 'figure')]
#     )
# def check_figure_type(n_clicks, fig):
#     print(fig['data'][1])
#     return 0


if __name__ == '__main__':
    app.run_server(debug=True)