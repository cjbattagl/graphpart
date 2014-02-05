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

output_file("read.html")
mat_name = 'ca-HepTh.mtx'
mat_info = sio.mminfo(mat_name)
mat_n = mat_info[0]
mat_nnz = mat_info[2]

print("matrix: {}".format(mat_name))
print("size: {} x {}".format(mat_n, mat_n))
print("nonzeros: {} (density: {})".format(mat_nnz, round(float(mat_nnz) / (mat_n * mat_n),6)))
print("nnz type: {} {}".format(mat_info[5], mat_info[4]))
assert(mat_n == mat_info[1])
print("reading matrix ..."),
A = sio.mmread(mat_name)
print("done\nconverting matrix to csr ..."),
A = A.tocsr()
print("done")

G = nx.to_networkx_graph(A)
numsamples = mat_n  #number of samples to take
numprocs = 4        #number of processors to simulate
s = 400             #starting vertex
procbin = float(mat_n)/numprocs
sampbin = floor(float(mat_n)/numsamples)
front = np.zeros(numprocs)
history = np.zeros((numprocs, numsamples))

step = 0

def share_same_proc(node1, node2):
  return (floor(node1 / procbin) == floor(node2 / procbin))
  
def do_shit():
  for node in bfs_edges(G, s):
    step += 1
    targnode = floor(node[1] / procbin)
    if (targnode != floor(node[0] / procbin)) :
      front[targnode] += 1
    if (step % sampbin == 0) :
      for proc in range(numprocs):
        history[proc, step/sampbin] = front[proc]
        
def do_1d(G, s):
  iter = 0
  numnodes = 0
  T = nx.bfs_tree(G,s)
  nextlevel = set()
  succ = T.successors(s)
  while (len(succ) > 0):
    for node in succ:
      nextlevel = nextlevel.union(T.successors(node))
      numnodes += 1
    succ = list(nextlevel)
    nextlevel.clear()
  print(numnodes)

do_1d(G,s)

#print(front)

x = np.linspace(0, numsamples, numsamples)
hold()

for proc in range(numprocs):
  line(x, history[proc,:], color="#0000FF",
    title = 'Load')
 
xaxis()[0].axis_label = "Chunk"
yaxis()[0].axis_label = "Load"
#show()

print("pct. communicated edges: %{}".format(round(100*np.sum(front)/step),5))
print("exp. communicated edges: %{}".format(round(100*((float(numprocs)-1)/numprocs),5)))
