# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
# import plotly.figure_factory as FF

from torusmesh import TorusMesh
from funcpath import FuncPath

n = 30
c, a = 2, 1

torus = TorusMesh(n,c,a)

start = (0,0)
end = (5,5)

p1, p2, p3 = torus.initial_path(start,end)

fig = torus.plotly_figure_realtime()

p1.plotly_path(fig)


app = dash.Dash(__name__)

app.layout = html.Div(children=[
    html.H1(children='Torus'),

    html.Div(children='''
        Particle dynamics on the Torus
    '''),

    dcc.Graph(
        id='example-graph',

		# figure = FF.create_trisurf(x=x, y=y, z=z,simplices=simplices, title="Torus", aspectratio=dict(x=1, y=1, z=0.3))
		figure = fig
	)

])

if __name__ == '__main__':
    app.run_server(debug=True)