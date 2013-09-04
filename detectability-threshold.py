import beliefprop as bp
import networkx as nx
import scipy.sparse as sp
import scipy.sparse.linalg as lin
import numpy as np
import matplotlib.pyplot as plt
import cProfile

def do():
  n_c = 20
  n_p = 20
  runs = 1
  samples = 3
  x = np.linspace(.5, .8, samples)
  #samples = 100
  #x = np.linspace(.5, 1, samples)
  y = [0.] * samples
  for dummy in range(runs):
    for index, p_cc in enumerate(x):
      print dummy, p_cc
      #p_cp = .5
      #p_pp = 1 - p_cc
      p_cp = 1 - p_cc
      p_pp = p_cc

      G = nx.Graph(bp.make_graph(n_c, n_p, p_cc, p_cp, p_pp))

      gamma = [n_c/float(n_c+n_p),n_p/float(n_c+n_p)]
      omega = [[p_cc, p_cp],[p_cp, p_pp]]

      A = np.array(nx.adjacency_matrix(G))
      print 'starting bp'
      ass = bp.bp_fixed_params(gamma,omega,A,10)
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

cProfile.run('do()')

#print perturbation
#print phi
#print "messages"
#print phi
#print omega
#print messages
#nx.draw(G)
