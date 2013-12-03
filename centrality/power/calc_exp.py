import powerlaw
from networkx import *

import scipy.sparse.linalg as lin
import scipy.sparse as sp
import numpy as np
import networkx as nx
import matplotlib.pylab as plt
from math import log
import os


def hashi_eig(A, degrees):
  I = sp.identity(A.shape[0])
  I_minus_D = sp.lil_matrix(A.shape)
  deg_list = list(degrees)
  for node,(_,deg) in enumerate(deg_list):
    I_minus_D[node,node] = 1.0-deg
  crazy = sp.bmat([[None,I],[I_minus_D,A]])
  eig = lin.eigs(crazy, k=1, which="LR")[1][:len(deg_list),0]
  root_total = np.sqrt(sum(x*x for x in eig))
  return [x/root_total for x in eig]

#for min_degree in range(50, 60, 1):
for f in os.listdir('networks4'):
  if f == 'citing_cited_ids.txt' or 'most' in f:
    continue
  #print(f, min_degree)
  #G = nx.read_gml('networks/%s'%f)
  G = nx.read_edgelist('networks4/%s'%f)
  #nodes_to_delete = []
  #for node in G.nodes():
  #  if len(G[node]) < min_degree:
  #    nodes_to_delete.append(node)
  #for node in nodes_to_delete:
  #  G.remove_node(node)
  #print(len(G.nodes()))
  data = [d for x,d in G.degree_iter()]
  '''
  hist = filter(lambda (a,b):b!=0 and a!=0,enumerate(nx.degree_histogram(G)))
  x = [log(degree) for degree, count in hist]
  y = [log(count) for degree, count in hist]
  '''
  results = powerlaw.Fit(data, discrete=True)
  #print(results.find_xmin())
  #results.power_law.plot_pdf()
  '''
  results.plot_cdf()
  plt.show()
  results.power_law.plot_cdf()
  plt.show()
  '''
  print(f, results.power_law.alpha)
  '''
  plt.figure()
  plt.title(f)
  plt.scatter(x, y)
  plt.show()
  '''
  #R, p = results.distribution_compare('power_law', 'lognormal')

  A = to_scipy_sparse_matrix(G, dtype='d', format='csr')
  print 'sparsified'
  eig = lin.eigsh(A, k=1, which='LA')
  print 'found  first eig'
  op = 0
  op = sum([x**4 for x in eig[1][:,0]])**(.5)
  eigh = hashi_eig(A, G.degree_iter())
  oph = sum([x**4 for x in eigh])**(.5)
  print(op, oph)
#plt.show()
