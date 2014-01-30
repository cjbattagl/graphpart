import scipy.io as sio

try:
    import matplotlib.pyplot as plt
except:
    raise

import networkx as nx

name = 'test.mtx'
info = sio.mminfo(name)
n = info[0]
nnz = info[2]

print("matrix:    {}".format(name))
print("size:      {} x {}.".format(n,n))
print("nonzeros:  {} (density: {}).".format(nnz, round(float(nnz) / (n * n),6)))
print("nnz type:  {} {}.".format(info[5],info[4]))
assert(n==info[1])

A = sio.mmread(name)
A = A.tocsr()

#G = nx.from_numpy_matrix(A)
