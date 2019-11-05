from torusmesh import TorusMesh
from funcpath import FuncPath

# 1 outer equator
# 2 inner equator
# 3 meridian
# 4 how close does v have to be to the inner equator?
# 5 general

n = 50
c, a = 2, 1

torus = TorusMesh(n,c,a)

ex1 = ((0,0),(torus.half, 0))

paths = torus.initial_path(ex1[0],ex1[1])
num_paths = len(paths)