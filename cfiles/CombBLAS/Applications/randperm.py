import numpy as np
import scipy.sparse as sparse
import scipy.io

print scipy.io.mminfo('./email-Enron.mtx')
A = scipy.io.mmread('./email-Enron.mtx')


#rr = range(10)
#np.random.shuffle(rr)

#for i in range(N):
#    M[:,i] = M[p,i]
#for i in range(N):
#    M[i,:] = M[i,p]
