#goal: communication model for 1D and 2D BFS,
#and finally a hybrid 1D/2D
import scipy.io as sio
import numpy as np
from collections import OrderedDict
import pandas as pd
import networkx as nx
from bokeh.plotting import *
from math import floor

try:
    import matplotlib.pyplot as plt
except:
    raise

############ Read in matrix, find basic stats ###########
output_file("read.html")
mat_name = 'ep.mtx'
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
gnumprocs = 4        #number of processors to simulate
s = 400             #starting vertex
gprocbin = float(mat_n)/gnumprocs
gsampbin = floor(float(mat_n)/numsamples)
front = np.zeros(gnumprocs)
history = np.zeros((gnumprocs, numsamples))

def proc_of_node(node):
  return int(floor(node / gprocbin))

def share_same_proc(node1, node2):
  return (proc_of_node(node1) == proc_of_node(node2))
  
def do_basic_1d():
  step = 0
  for node in bfs_edges(G, s):
    step += 1
    targnode = floor(node[1] / gprocbin)
    if (~share_same_proc(node[0], node[1])) :
      front[targnode] += 1
    if (step % gsampbin == 0) :
      for proc in range(gnumprocs):
        history[proc, step/gsampbin] = front[proc]       

# return nodes at each level
def BFS_levels(G, s):
  if G:
    T = nx.bfs_tree(G,s)
    nextlevel = set()
    succ = T.successors(s)
    while (len(succ) > 0):
      for node in succ:
        nextlevel = nextlevel.union(T.successors(node))
      succ = list(nextlevel)
      nextlevel.clear()
      yield succ

# print number of nodes at each level
def compute_frontier(G, s):
  for level in BFS_levels(G, s):
    print(len(level))
    
gmaxlevels = 14
vec_children = np.zeros(gmaxlevels)
vec_failedchildren = np.zeros(gmaxlevels)
vec_peers = np.zeros(gmaxlevels)
vec_parents = np.zeros(gmaxlevels)
    
# return no. of communicated edges at each level
def BFS_levels_comm(G, s):
  if G:
    T = nx.bfs_tree(G,s)
    nextlevel = set()
    succ = T.successors(s)
    prevlevel = set()
    prevlevel.add(s)
    level = 0
    
    children = [set() for proc in range(gnumprocs)] # ALL children
    failedchildren = [0 for proc in range(gnumprocs)]
    peers = [set() for proc in range(gnumprocs)]
    failed_parents = [set() for proc in range(gnumprocs)]

    while (len(succ) > 0):
      level += 1
      comm = 0
      for node in succ:
        node_succ = set(T.successors(node))
        nodeproc = proc_of_node(node)
        nodeneighbs = set(G.neighbors(node))
        # all children
        children[nodeproc] = children[nodeproc].union(nodeneighbs.intersection(node_succ))      
        # number of failed children
        failedchildren[nodeproc] += len(children[nodeproc].intersection(nodeneighbs.intersection(node_succ))) 
        # peer
        peers[nodeproc] = peers[nodeproc].union(nodeneighbs.intersection(succ))
        # failed parent
        failed_parents[nodeproc] = failed_parents[nodeproc].union(nodeneighbs.intersection(prevlevel))   
        for targnode in node_succ:
          if (share_same_proc(node, targnode) != 1) :
            comm += 1
        nextlevel = nextlevel.union(node_succ)
      prevlevel = succ
      succ = list(nextlevel)
      nextlevel.clear()
      if (level < gmaxlevels):
        vec_children[level] += len(succ)
        for i in range(gnumprocs):
          vec_failedchildren[level] += failedchildren[i]
          vec_peers[level] += len(peers[i])
          vec_parents[level] += len(failed_parents[i])
          children[i].clear()
          peers[i].clear()
          failed_parents[i].clear()
          failedchildren[i] = 0
      yield comm
      
def compute_comm(G, s):
  for comm in BFS_levels_comm(G, s):
    print(comm)

    
compute_comm(G,s)
#compute_frontier(G, s)

print("pct. communicated edges: %{}".format(round(100*np.sum(front)/step),5))
print("exp. communicated edges: %{}".format(round(100*((float(gnumprocs)-1)/gnumprocs),5)))

import colorsys
def get_N_HexCol(N=5):
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    RGB_tuples = map(lambda x: tuple(map(lambda y: int(y * 255),x)),RGB_tuples)
    HEX_tuples = map(lambda x: tuple(map(lambda y: chr(y).encode('hex'),x)), RGB_tuples)
    HEX_tuples = map(lambda x: "".join(x), HEX_tuples)
    return HEX_tuples

def procplot():
  N = gnumprocs
  categories = ['children', 'failedchildren', 'peers', 'parents']
  data = {}
  data['x'] = np.arange(gmaxlevels)
  data['parents'] = vec_parents.copy()
  data['children'] = vec_children.copy()
  data['failedchildren'] = vec_failedchildren.copy()
  data['peers'] = vec_peers.copy() 
  df = pd.DataFrame(data)
  df = df.set_index(['x'])
  colors = get_N_HexCol(N)
  def stacked(df, categories):
      areas = OrderedDict()
      last = np.zeros(len(df[categories[0]]))
      for cat in categories:
        next = last + df[cat]
        areas[cat] = np.hstack((last[::-1], next))
        last = next
      return areas
  output_file("edges.html", title="bfs edges")
  areas = stacked(df, categories)
  #colors = brewer["Spectral"][len(areas)]
  x2 = np.hstack((data['x'][::-1], data['x']))
  patches([x2 for a in areas], list(areas.values()), color=colors, alpha=0.8, line_color=None)
  show()

procplot()
