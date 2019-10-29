import plotly.graph_objects as go

from torusmesh import TorusMesh
from funcpath import FuncPath

n = 30
c, a = 2, 1

torus = TorusMesh(n,c,a)

start = (0,0)
end = (5,5)

p1, p2, p3 = torus.initial_path(start,end)

# generates base figure
fig = torus.plotly_figure_realtime()

p1.plotly_path(fig)

def check_figure_type(fig):
    return'''
        fig is of type "{}"
    '''.format(type(fig))

print(check_figure_type(fig))