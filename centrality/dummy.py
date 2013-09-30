from networkx import *
import scipy.sparse.linalg as lin
import scipy.sparse as sp
import matplotlib.pylab as plt
import numpy as np

def hashi_eig(A, G):
  I = sp.identity(A.shape[0])
  I_minus_D = sp.lil_matrix(A.shape)
  n = len(G.nodes())
  for node,deg in G.degree_iter():
    I_minus_D[node,node] = 1.0-deg
  crazy = sp.bmat([[None,I],[I_minus_D,A]])
  eig = lin.eigs(crazy, k=1, which="LR")[1][:n,0]
  root_total = np.sqrt(sum(x*x for x in eig))
  return [x/root_total for x in eig]

def dummy(adj, graph, position, colors=None):
    if colors == None:
        colors = []
        for n in graph.nodes():
            if n == 0:
                colors.append('b')
            elif n in graph.neighbors(0):
                colors.append('g')
            else:
                colors.append('r')
    eig = lin.eigsh(adj, k=1, which='LA')[1][:,0]
    print sum(x**2 for x in eig)
    print eig[:100]
    sizes = map(lambda x:abs(x*x*1000000),eig)

    eigh = hashi_eig(adj, graph)
    print sum(x.real**2 for x in eigh)
    print eigh[:100]
    sizesh = map(lambda x:abs(x.real*x.real*1000000),eigh)
    f1 = plt.figure(figsize=(7,7))
    print 'done finding eig'
    networkx.draw(graph, linewidths=.1,pos=position, node_color = colors, width=.05, node_size = sizes, with_labels=False)

    f2 = plt.figure(figsize=(7,7))
    networkx.draw(graph, linewidths=.1,pos=position, node_color = colors, width=.05, node_size = sizesh, with_labels=False)
    
    return eig, eigh, f1, f2
    
def dummy2(adj, graph):
    eig = lin.eigsh(adj, k=1, which='LA')
    sizes = map(lambda x:abs(x*x*1000000),eig[1][:,0])

    I = sp.identity(adj.shape[0])
    I_minus_D = sp.lil_matrix(adj.shape)
    for node,deg in graph.degree_iter():
      I_minus_D[node,node] = 1.0-deg
    crazy = sp.bmat([[None,I],[I_minus_D,adj]])
    
    eigh = lin.eigs(crazy, k=1, which="LR")
    sizesh = map(lambda x:abs(x.real*x.real*1000000*200),eigh[1][:,0:len(sizes)])
    return sizes[0], sizesh[0]
