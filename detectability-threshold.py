import decellebeliefprop as bp
import networkx as nx
import scipy.sparse as sp
import scipy.sparse.linalg as lin
import numpy as np
import matplotlib.pyplot as plt
import cProfile
import pdb

def do():
  n_c = 40
  n_p = 40
  n = float(n_c+n_p)
  runs = 1
  samples = 3
  min_diff = 10
  max_diff = 20
  #x = np.linspace(.6, .6, samples)
  x = np.linspace(min_diff, max_diff, samples)
  #samples = 100
  #x = np.linspace(.5, 1, samples)
  y = [0.] * samples
  for dummy in range(runs):
    #for index, p_cc in enumerate(x):
    for index, diff in enumerate(x):
      c_cp = diff / 2 + 3
      c_cc = c_cp + diff/2
      c_pp = c_cp - diff/2
      print dummy, p_cc
      p_cp = diff
      p_pp = 1 - p_cc
      #p_cp = 1 - p_cc
      #p_pp = p_cc

      G = nx.Graph(bp.make_graph(n_c, n_p, p_cc, p_cp, p_pp))

      #pdb.set_trace()

      gamma = [n_c/n,n_p/n]
      omega = [[p_cc*n, p_cp*n],[p_cp*n, p_pp*n]]

      A = np.array(nx.adjacency_matrix(G))
      print 'starting bp'
      ass = bp.bp_fixed_params(gamma,omega,G,5)
      nx.draw(G, node_color=ass)
      print ass
      plt.show()
      print 'ending bp'

      score = 0
      for index2, val in enumerate(ass):
        if index2 < n_c and val > .5:
          score += 1
        elif val < .5:
          score += 1
      y[index] += score
    print y

  y = [val / ((n_c+n_p)*runs) for val in y]
  #plt.plot(x, y)
  #plt.show()

do()

#print perturbation
#print phi
#print "messages"
#print phi
#print omega
#print messages
#nx.draw(G)
