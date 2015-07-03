
sparse=[]
M0={0:16250., 3:15625., 4:625.}
sparse.append(M0)
M1={1:45509., 5:32400., 6:2500.,  7:10609.}
sparse.append(M1)
M2={2:58634., 3:15625., 5:32400., 7:10609.}
sparse.append(M2)
M3={3:31250.}
sparse.append(M3)
M4={4:625.}
sparse.append(M4)
M5={5:64800.}
sparse.append(M5)
M6={6:2500.}
sparse.append(M6)
M7={7:21218.}
sparse.append(M7)
print "**************************Sparse Normal Matrix****************************"
for i in sparse: print i

import time
from elimination_graph import graph
import numpy as np
c= graph(sparse)
print "**************************Sparse Decomposition****************************"
t1=time.time()
L=c.cholesky_UTU_sparse()
print "sparse time", time.time()-t1
for i in L: print i


m1=[[ 16250.,      0.,      0.,  15625.,    625.,      0.,      0.,      0.],
    [     0.,  45509.,      0.,      0.,      0.,  32400.,   2500.,  10609.],
    [     0.,      0.,  58634.,  15625.,      0.,  32400.,      0.,  10609.],
    [ 15625.,      0.,  15625.,  31250.,      0.,      0.,      0.,      0.],
    [   625.,      0.,      0.,      0.,    625.,      0.,      0.,      0.],
    [     0.,  32400.,  32400.,      0.,      0.,  64800.,      0.,      0.],
    [     0.,   2500.,      0.,      0.,      0.,      0.,   2500.,      0.],
    [     0.,  10609.,  10609.,      0.,      0.,      0.,      0.,  21218.]]


print "**************************Dense Normal Matrix****************************"
print np.array(m1)
def cholesky_UTU(A):
    import math
    L = [[0.0] * len(A) for _ in xrange(len(A))]
    c=0
    for i in xrange(len(A)):
      s=0
      for k in xrange(i):
        s = s+L[k][i]**2
        c+=1
#      s=round(s,4)
      L[i][i] = math.sqrt(A[i][i] - s)
      for j in xrange(i+1, len(A)):
        s=0
        for k in xrange(i):
          s = s+L[k][i] * L[k][j]
          c+=1
        L[i][j]=(1.0 / L[i][i] * (A[i][j] - s))
    print "Total n of operations", c
    return L
print "**************************Dense Decomposition****************************"
t1=time.time()
d=cholesky_UTU(m1)
print "dense time", time.time()-t1
d=np.array(d)
np.set_printoptions(precision=2,suppress=True, linewidth=200, threshold=10000)
print d

