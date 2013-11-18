from networkx import *
from custom_networkx import *
import scipy.sparse.linalg as lin
import scipy.sparse as sp
import numpy as np
import matplotlib.pylab as plt
import time
import glob
import cProfile
import random

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

def hashi_eigval(A, degrees):
  I = sp.identity(A.shape[0])
  I_minus_D = sp.lil_matrix(A.shape)
  deg_list = list(degrees)
  for node,(_,deg) in enumerate(deg_list):
    I_minus_D[node,node] = 1.0-deg
  crazy = sp.bmat([[None,I],[I_minus_D,A]])
  print(lin.eigs(crazy, k=1, which="LR",return_eigenvectors=False)[0])

c = 6.
d = 100
'''
ij = np.random.randint(n, size=(2, n*c))
ji = np.array([ij[1],ij[0]])
data = np.random.randint(2, size=n*c)
#print(ij[0], ji[0])
ij = np.array([np.append(ij[0], ji[0]),
              np.append(ij[1], ji[1])])
print(len(data))
data = np.append(data, data)
print(len(data))
X = sp.coo.coo_matrix((data,ij), (n, n))
#print(ij)
#print(data)

#make a
#ij = sp.bmat([[np.random.randint(n, size=d)],[None]])
#ij = sp.bmat([[None],[np.random.randint(n, size=d)]])
#data = np.random.randint(2, size=d)
#print(ij)
#print(data)
#a = sp.coo.coo_matrix((data,ij), (n,1))
#a_T = sp.coo.coo_matrix((data,ji), (1,n))
#a = np.random.binomial(1,d/n,size=(n,1))
#np.append(a, 0)
#print(a)
a_T = np.transpose(a)
'''

'''
X = sp.lil_matrix((n,n))
a = sp.lil_matrix((n,1))
a_T = sp.lil_matrix((1,n))
for i in range(n):
  print(i)
  if random.random() < d/n:
    a[i,0] = 1
    a_T[0,i] = 1
  for j in range(n):
    if random.random() < c/n:
      X[i,j] = 1
'''
#print(X)
#X = sp.coo_matrix((n,n))
#for index in range(n*c):
#dot_arr = []
#diff_arr = []
mean0 = []
mean1 = []
mean2 = []
mean3 = []
std1 = []
std2 = []
std3 = []
n_arr = range(20000, 21000, 5)
for n in n_arr:
  trials = 1
  #center = np.empty((n,n))
  #center.fill(c/n)
  #u = np.empty(n)
  #u.fill(1./(n ** .5))
  #u = np.append(u, 0)
  dot = 0
  diff = 0
  for dummy in range(trials):
    A = to_scipy_sparse_matrix(erdos_renyi_graph(n, c/n))
    a = np.random.binomial(1,d/n,size=(n,1))
    a_T = np.transpose(a)
    #X = A - center
    #B = sp.bmat([[X,a],[a_T, None]], dtype='d')
    A = sp.bmat([[A,a],[a_T, None]], dtype='d', format='csr')
    #eigB = lin.eigsh(B, k=1, which='LA')
    eigA = lin.eigsh(A, k=1, which='LA')[1][:,0]
    #dot += np.dot(u, eigB[1][:,0]) ** 2
    #diff += (eigB[0] - eigA[0]) ** 2
    
    #compiling node lists
    #print('start')
    print(n)
    d0 = set([n])
    d1 = set([index for index, item in enumerate(a) if item == 1])
    d2 = set([])
    d3 = set([])
    for node in d1:
      d2.update(A[node].indices)
    d2 = d2.difference(d1)
    d2 = d2.difference(d0)
    for node in d2:
      d3.update(A[node].indices)
    d3 = d3.difference(d2)
    d3 = d3.difference(d1)
    d0_evals = np.array(abs(eigA[n]))
    d1_evals = np.array([])
    d2_evals = np.array([])
    d3_evals = np.array([])
    for node in d1:
      d1_evals = np.append(d1_evals, abs(eigA[node]))
    for node in d2:
      d2_evals = np.append(d2_evals, abs(eigA[node]))
    for node in d3:
      d3_evals = np.append(d3_evals, abs(eigA[node]))

    mean0.append(d0_evals)
    mean1.append(np.mean(d1_evals))
    mean2.append(np.mean(d2_evals))
    mean3.append(np.mean(d3_evals))
    std1.append(np.std(d1_evals))
    std2.append(np.std(d2_evals))
    std3.append(np.std(d3_evals))
    print(len(d1))
    print(len(d2))
    print(len(d3))
  #print(n)
  #print('%f %f' % (dot, diff))
  #dot_arr.append(n * dot/trials)
  #diff_arr.append(n * diff/trials)

plt.figure()
plt.plot(n_arr, mean0)
plt.figure()
plt.plot(n_arr, mean1)
plt.figure()
plt.plot(n_arr, mean2)
plt.figure()
plt.plot(n_arr, mean3)
plt.figure()
plt.plot(n_arr, std1)
plt.figure()
plt.plot(n_arr, std2)
plt.figure()
plt.plot(n_arr, std3)
plt.show()
