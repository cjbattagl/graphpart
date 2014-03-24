#goal: communication model for 1D and 2D BFS,
#and finally a hybrid 1D/2D

#1D 2D hybrid: choose tuning degree. sort by degree
#streaming partition the low degree nodes to the p processors
#2D partition the hi degree nodes
#each node locally stores the visited and distance vectors for their own low-degree nodes
#each node replicates these vectors for high-degree nodes (thus nodes can avoid recommunicating them)
#show that this reduces communication enormously because hi-degree nodes dominate

import scipy.io as sio
import numpy as np
from collections import OrderedDict
import pandas as pd
import networkx as nx
from bokeh.plotting import *
from math import floor
from stacked_graph import *

try:
    import matplotlib.pyplot as plt
except:
    raise

############ Read in matrix, find basic stats ###########
output_file("read.html")
mat_name = 'ca-AstroPh.mtx'
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
s = 500             #starting vertex
ghideg = 60
gprocbin = float(mat_n)/gnumprocs
gsampbin = floor(float(mat_n)/numsamples)
front = np.zeros(gnumprocs)
history = np.zeros((gnumprocs, numsamples))

def proc_of_node(node):
  return int(floor(node / gprocbin))

def degree_stats(G):
  max = 0
  for (node, deg) in G.degree_iter():
    if (deg > max):
      max = deg
      print(max)
  print max
  
#degree_stats(G)

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
vec_hideg = np.zeros(gmaxlevels)
    
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
      hidegedges = 0
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
        for targnode in nodeneighbs:
          if (G.degree(targnode) > ghideg):
            hidegedges += 1
        for targnode in node_succ:
          if (share_same_proc(node, targnode) != 1) :
            comm += 1
        nextlevel = nextlevel.union(node_succ)
      prevlevel = succ
      succ = list(nextlevel)
      nextlevel.clear()
      if (level < gmaxlevels):
        #vec_children[level] += len(succ)
        vec_hideg[level] += hidegedges
        for i in range(gnumprocs):
          vec_failedchildren[level] += failedchildren[i]
          vec_children[level] += len(children[i])
          vec_peers[level] += len(peers[i])
          vec_parents[level] += len(failed_parents[i])

          children[i].clear()
          peers[i].clear()
          failed_parents[i].clear()
          failedchildren[i] = 0
          hidegedges = 0
      yield comm
      
def compute_comm(G, s):
  for comm in BFS_levels_comm(G, s):
    print(comm)

    
#compute_comm(G,s)
#compute_frontier(G, s)

#vec_children = np.zeros(gmaxlevels)
#vec_failedchildren = np.zeros(gmaxlevels)
#vec_peers = np.zeros(gmaxlevels)
#vec_parents = np.zeros(gmaxlevels)
#vec_hideg = np.zeros(gmaxlevels)

if 1 == 1: #__name__ == '__main__':
    pl.clf()
    N_dsets = 4
    T = gmaxlevels
    amp = 1
    fade = .15
    dsets = [vec_children, vec_failedchildren, vec_peers, vec_parents]
    print(dsets)
    #for i in xrange(N_dsets):
    #    this_dset = np.zeros(T)
    #    t_onset = np.random.randint(.9*T)-T/3
    #    if t_onset >= 0:   
    #        remaining_t = np.arange(T-t_onset)     
    #    else:
    #        remaining_t = np.arange(T)-t_onset
    #   this_dset[max(t_onset,0):]=np.exp(-.15*np.random.gamma(10,.1)*remaining_t)\
    #                       * remaining_t * np.random.gamma(6,.2)# * np.cos(-fade*remaining_t*np.random.gamma(10,.1))**2
    #   dsets.append(this_dset)
    stacked_graph(dsets, baseline_fn = min_weighted_wiggles, color_seq='random')

pl.show()
