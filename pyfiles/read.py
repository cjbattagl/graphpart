#goal: communication model for 1D and 2D BFS,
#and finally a hybrid 1D/2D

import scipy.io as sio
import numpy as np
import networkx as nx
from bokeh.plotting import *
from math import floor

try:
    import matplotlib.pyplot as plt
except:
    raise
    
def my_bfs_edges(G,source):
    visited = set([source])
    stack = [(source,iter(G[source]))]
    while stack:
        parent,children = stack[0]
        try:
            child = next(children)
            if child not in visited:
                yield parent,child
                visited.add(child)
                stack.append((child,iter(G[child])))
        except StopIteration:
            stack.pop(0)

output_file("read.html")
name = 'ca-AstroPh.mtx'
info = sio.mminfo(name)
n = info[0]
nnz = info[2]

print("matrix: {}".format(name))
print("size: {} x {}".format(n,n))
print("nonzeros: {} (density: {})".format(nnz, round(float(nnz) / (n * n),6)))
print("nnz type: {} {}".format(info[5],info[4]))
assert(n==info[1])
print("reading matrix ..."),
A = sio.mmread(name)
print("done\nconverting matrix to csr ..."),
A = A.tocsr()
print("done")

G = nx.to_networkx_graph(A)
step = 0
numsamples = n
p = 4
s = 500
procbin = float(n)/p
sampbin = floor(float(n)/numsamples)
front = np.zeros(p)
history = np.zeros((p, numsamples))

def share_same_proc(node1, node2):
  return (floor(node1 / procbin) == (floor(node2 / procbin)

for node in my_bfs_edges(G, s):
  step += 1
  targnode = floor(node[1] / procbin)
  if (targnode != floor(node[0] / procbin)) :
    front[targnode] += 1
  if (step % sampbin == 0) :
    for proc in range(p):
      history[proc, step/sampbin] = front[proc]
      
T = nx.bfs_tree(G,s)
succ = T.successors(s)
nextlevel = []
while (len(succ) > 0):
  for node in succ:
    
    nextlevel = nextlevel + T.successors(node)
    

print(front)

x = np.linspace(0, numsamples, numsamples)
hold()

for proc in range(p):
  line(x, history[proc,:], color="#0000FF",
    title = 'Load')
 
xaxis()[0].axis_label = "Chunk"
yaxis()[0].axis_label = "Load"
#show()

print("pct. communicated edges: %{}".format(round(100*np.sum(front)/step),5))
print("exp. communicated edges: %{}".format(round(100*((float(p)-1)/p),5)))
