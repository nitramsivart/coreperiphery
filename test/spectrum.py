import sage.all
import networkx as nx
import numpy as np
import matplotlib.pylab as plt
from scipy import linalg

def cavity(A, flow=False):
  # the index of the edge in the dictionary
  index = 0
  edges = {}
  for i in range(len(A)):
    for j in range(len(A)):
      if A[i][j] == 1:
        edges[index] = (i, j)
        index += 1
  print len(edges)
  B = np.zeros([len(edges), len(edges)])
  for e1 in range(len(B)):
    for e2 in range(len(B)):
      k, l = edges[e1]
      i, j = edges[e2]
      if j == k and i != l:
        if flow:
          B[e1][e2] = 0 if sum(A[j])==1 else 1./float(sum(A[j])-1)
        else:
          B[e1][e2] = 1
  return B, edges

def generate_tree():
  size = 256 
  A = np.zeros([size, size])
  A[0][1] = A[0][2] = A[1][0] = A[2][0] = 1
  for i in range(1, size):
    for j in range(i, size):
      if sum(A[i]) < 3 and sum(A[j]) == 0:
        A[i][j] = A[j][i] = 1
    
  return A

def make_2core(A):
  while True:
    for i in range(len(A)):
      if sum(A[i]) == 1:
        A = np.delete(A, i, axis=1)
        A = np.delete(A, i, axis=0)
        break
      if i == len(A)-1:
        return A

G = nx.read_gml('karate.gml')
A = np.array(nx.adjacency_matrix(G))
#A = make_2core(A)
#G = nx.Graph(A)
#A = generate_tree()
#G = nx.Graph(A)
B, edges = cavity(A, True)


'''
q = 2
n = len(G.nodes())
gamma =[0.5,0.5]
omega =[(0.25,0.04),(0.04,.25)]
tmax = 80
ass,phi,gamma,omega,messages = bp(gamma,omega,A,tmax)

for u in range(5):
  s = ""
  for v in range(5):
    s += "%d -> %d: " % (u,v)
    for r in range(q):
      s += "%f, " % messages[u][v][r]
    s += '\n'
  print s,

nx.draw(G,node_color =ass)
plt.show()
'''
one = np.ones([len(B),len(B)])/len(B)
#print one
m = sage.all.matrix(B-one)
values = np.array(m.eigenvalues())
#print values
vectors = m.eigenvectors_right()
'''
#print vectors
first = values.argmax()
f1 = open('vectors.txt', 'w+')
for v in vectors:
  f1.write(str(v))
f1.close()
'''
'''
print values[first], vectors[first]
print sum(B.dot(vectors[first]))
print sum(vectors[first])
'''
#print B.dot(np.ones([len(B)]))


#values[first] = 0
second = values.argmax()
#print values[second], vectors[second]
colors = np.zeros([len(A)])
#counts = np.zeros([len(A)])
for i, e in enumerate(vectors[second][1][0]):
  colors[edges[i][1]] += e

print colors
colors = [1 if c>0 else -1 for c in colors]

nx.draw(G, node_color=colors)
plt.show()
plt.savefig('test.png')
#plt.scatter([k.real for k in values], [k.imag for k in values])
#plt.plot([k.real for k in evs])
#plt.plot([k.imag for k in evs])
#plt.show()
