import scipy.io as sio
import numpy as np
import math
from bokeh.plotting import *

try:
    import matplotlib.pyplot as plt
except:
    raise

import networkx as nx

# This is a test to explore the memory access patterns of a parallel graph traversal

def getProc(vert, bindiv):
  return math.floor(vert/bindiv)

output_file("proc.html", title="bfs proc example")

name = 'test.mtx'
info = sio.mminfo(name)
n = info[0]
nnz = info[2]

p = 10
numsamples = 100

print("matrix:    {}".format(name))
print("size:      {} x {}".format(n,n))
print("nonzeros:  {} (density: %{})".format(nnz, 100*round(float(nnz) / (n * n),6)))
print("nnz type:  {} {}".format(info[5],info[4]))
assert(n==info[1])

A = sio.mmread(name)
print(" - > converting to csr ... "),
A = A.tocsr()
print("done\n - > converting to networkx graph ..."),
G = nx.to_networkx_graph(A)
print("done")

front = np.zeros(p)
loads = np.zeros((p,numsamples))
bindiv = float(n)/p
samplediv = math.floor(float(n)/numsamples)
step = 0
x = np.linspace(0,n,numsamples)

for val in nx.bfs_edges(G,20):
  front[getProc(val[1],bindiv)] += 1
  step += 1
  if (step % samplediv == 0):
    for proc in range (p):
      loads[proc][step/samplediv] = front[proc]
      front[proc] = 0

hold()
for proc in range (p) :         
  line(x, loads[proc], color="#0000FF", x_axis_type = "datetime", 
    tools="pan,zoom,resize", width=1200,height=300, 
    title = 'Proc Activity') 
xaxis()[0].axis_label = "Step"
yaxis()[0].axis_label = "Num Nodes Active"
show()
