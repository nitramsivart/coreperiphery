import networkx as nx
import numpy as np
#from oldbeliefprop import bp
from decellebeliefprop import bp
import matplotlib.pylab as plt

G = nx.read_gml('test/karate.gml')
A = np.array(nx.adjacency_matrix(G))

q = 2
n = len(G.nodes())
gamma =[0.5,0.5]
#omega =[(0.9,0.8),(0.8,.1)]
a = .5
omega =[(0.22532584,0.09263349),(0.09263349,0.22532584)]
omega =[(0.1912654,0.07863089),(0.07863089,0.1912653)]
omega= [(0.500001,0.5),(0.5,0.5)]
omega = [(0.13494815, 0.13494805),(0.13494805, 0.13494815)] # cp
omega= [(0.21366782,0.05622837),(0.05622837,0.21366782)] # com
omega =[(.7601,0.2),(0.2,.76)] # community
omega =[[.6,0.3],[0.3,.5]] # core periph
tmax = 100
print "edges:", len(G.edges())
print 'nodes:', len(G.nodes())
# the following loop turns p's into c's
for u,v in [(0,0),(0,1),(1,0),(1,1)]:
  omega[u][v] = omega[u][v] * len(G.nodes())
ass,phi,gamma,omega,messages = bp(gamma,omega,A,tmax)

'''
for u in range(15):
  s = ""
  for v in range(5):
    s += "%d -> %d: " % (u,v)
    for r in range(q):
      val = (messages[v][u][r] - .50001) / (messages[0][1][0] - .50001)
      s += "%s, " % val
    s += '\n'
  print s,
'''
print phi

nx.draw(G,node_color =ass)
plt.show()
