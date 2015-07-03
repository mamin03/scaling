from __future__ import division
import math

class graph:               # the graph data structure that perform sparse Cholesky decomposition
   def __init__(self, A):  # the constructor will get the normal matrix as a list of dictionary objects (look at the test file)
     self.g_of_A=[]        # the graph is an array of vertices, each vertex is an array of edges
     self.n=len(A)         # the dimention of the array
     self.L_column=[]      # stores the nonzero in each column
     self.L=A              # Initiate L with A because all nonzero elements in A are nonzero in L 
     self.dependencies={}  # This is a dictionary will carry the elements we need to loop over for calculating the Lij terms 
     for i in A:           # constructing graph data structure. Each row is a vertex 
       self.g_of_A.append(i.keys()) # there is an edge between vertex i and j if the j element in the ith row is nonzero

     for i in range(self.n):
       for j in A[i].keys()[1:]:
         self.g_of_A[j].append(i)  # if there is an edge between i and j then there is and edge between j and i (undirected graph)
     self.set_non_zero_columns()   # this function will set the non zero elements in L without calculating their values

   def remove_vertex(self, ID):    # this function remove a vertex with ID=ID. 
      edges=self.g_of_A[ID][1:]    # since each vertex is connected to itself start from the second edge
      for i in range(len(edges)):  # loop over each element 
        for j in range(i+1, len(edges)): # and connect it with the following elements
          if edges[j] not in self.g_of_A[edges[i]]: # make sure this connection hasn't been added before
            self.L[edges[i]][edges[j]]=0            # new connection mean a fill in in the L matrix that is initialized to zero value
            self.g_of_A[edges[i]].append(edges[j])  # the graph has to be updated with the new connections
      for i in self.g_of_A: # this should be improved. You may loop only on the vertices in edges
        if ID in i:         # I will do that in the C++ code
         i.remove(ID)


   def set_non_zero_columns(self): # The performance is much faster if we have another array that stores
     self.set_nonzero_L()          # the nonzero elements in each column as well
     for i in range(self.n):       # initialize an array of columns
       col=[]
       self.L_column.append(col)
     for i in range(self.n):       # loop over the nonzero terms in each row and store them in the corresponding column 
       for j in self.L[i].keys()[1:]:
         self.L_column[j].append(i)
     for i in range(self.n):       # we could save more time here by storing the elements needed to calculate the Lij terms
       for j in self.L[i].keys()[1:]: # we are voiding any calculations that mutiply by zero
         list=[]
         for k in self.L_column[i]:
           Lij=self.L[k].get(j)
           if Lij is not None:
               list.append(k)
         self.dependencies[(i,j)]=list # each (i,j) term has a list contains the elements you need to loop over

   def set_nonzero_L(self):   # this function is just calling remove_vertex for all vertices in the graph
     for i in range(self.n):
       self.remove_vertex(i)

   def cholesky_UTU_sparse(self): # this is the function that calculate the decomposed sparse matrix
     n=self.n
     for i in xrange(n):   # loop over all diagonal elements (diagonal elements are nonzero)
       s=0
       for k in self.L_column[i]: # to calculate the (i,i) element you need to sum the squares of all elements above it in the same column i
         s+=self.L[k][i]**2       # self.L_column[i] stores the indices of the nonzero elements in the column i
       self.L[i][i]=math.sqrt(self.L[i][i] - s) # the initial value of self.L[i][i] is A[i][i] which is needed for calculating L[i][i]
       for j in self.L[i].keys()[1:]:  # after calculating the (i,i) element start calculating (i,j) elements
         s=0
         for k in self.dependencies[(i,j)]: # for each (i,j) elements loop only of the terms that the product 
             s+=self.L[k][i]*self.L[k][j]   # self.L[k][i]*self.L[k][j] is nonzero
         self.L[i][j]=(1.0 / self.L[i][i] * (self.L[i][j] - s))  # finally evaluate the (i,j) term 
     return self.L
