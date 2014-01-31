import scipy.io as sio
import numpy as np
from bokeh.plotting import *
from math import floor

try:
    import matplotlib.pyplot as plt
except:
    raise

import networkx as nx

output_file("read.html")

name = 'test.mtx'
info = sio.mminfo(name)
n = info[0]
nnz = info[2]

print("matrix:    {}".format(name))
print("size:      {} x {}".format(n,n))
print("nonzeros:  {} (density: {})".format(nnz, round(float(nnz) / (n * n),6)))
print("nnz type:  {} {}".format(info[5],info[4]))
assert(n==info[1])

print("reading matrix ..."),
A = sio.mmread(name)
print("done\nconverting matrix to csr ..."),
A = A.tocsr()
print("done")

G = nx.to_networkx_graph(A)
step = 0
numsamples = 100
p = 10

procbin = float(n)/p
sampbin = math.floor(float(n)/numsamples)

front = np.zeros(p)
history = np.zeros((p, numsamples))


for node in nx.bfs_edges(G, 1):
  step += 1
  targnode = floor(node[1] / procbin)
  front[targnode] += 1
  if (step % numsamples == 0) :
    for proc in range(p):
      history[proc, step/sampbin] = front[proc]

x = np.linspace(0, 100, 100) 

hold()

for proc in range(p):
  line(x, history[proc,:], color="#0000FF", x_axis_type = "datetime", 
    tools="pan,zoom,resize", width=1200,height=300, 
    title = 'Load')
 
xaxis()[0].axis_label = "Chunk"
yaxis()[0].axis_label = "Load"  
show()