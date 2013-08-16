from networkx import *
import scipy.sparse.linalg as lin
import scipy.sparse as sp
import numpy as np
import matplotlib.pylab as plt

def cavity(G):
  # the index of the edge in the dictionary
  index = 0
  DG = DiGraph(G)
  #edges = G.edges()
  B = Graph()
  for e1 in DG.edges_iter():
    for e2 in DG.edges_iter():
      k, l = e1
      i, j = e2
      if j == k and i != l:
        B.add_edge(e1, e2)
  return B

c = 100
k = 200
n = 10000
print "making graph"
G = fast_gnp_random_graph(n, c/float(n))
print "done making graph"

print "making hub"
G.add_node(n)
#for i in range(k):
for i in range(c):
  G.add_edge(n, i)

G.add_node(n+1)
#for i in range(k+2):
for i in range(2*c):
  G.add_edge(n+1, n-(i+1))

print "done making hub"

'''
print "making cavity"
BG = cavity(G)
print "done with cavity"
'''


'''
for i in range(k):
  G.add_edge(1, i)
'''

'''
eigvals = adjacency_spectrum(G)
print "done finding spectrum"
r = [e.real for e in eigvals]
i = [e.imag for e in eigvals]
plt.hist(r, bins=100)
plt.show()
exit()
'''

A = to_scipy_sparse_matrix(G, dtype='d')
I = sp.identity(A.shape[0])
I_minus_D = sp.lil_matrix(A.shape)
for node,deg in G.degree_iter():
  I_minus_D[node,node] = 1.0-deg
print "done making submatrices"

crazy = sp.bmat([[None,I],[I_minus_D,A]])
print "done making crazy"


eig = lin.eigsh(A, k=3, which="LA")
#eig = lin.eigh(adjacency_matrix(G), eigvals=(n-1, n+1))
print eig[0]
print eig[1][n]
print eig[1][n+1]
print sum(eig[1][:,2])

#eig2 = lin.eigs(to_scipy_sparse_matrix(BG, dtype='d'), k=1, which = "LR")
eig2 = lin.eigs(crazy, k=1, which="LR")
print eig2[1][n].real
print eig2[1][n+1].real
print sum(eig2[1][0:n+2]).real

#for e in BG.nodes_iter():
  #print e
exit()
f = open("flow.txt", 'w+')
for v in eig2[1]:
  f.write(str(v)+'\n')
f.close()

exit()
f = open("adj.txt", "w+")
for v in eig[1]:
  f.write(str(v.real)+"\n")

#print eigenvector_centrality(G)
'''
for v in range(len(eig[1][0])):
  f.write(str(eig[1][:,v]))
'''
