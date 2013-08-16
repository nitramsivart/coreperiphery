from networkx import *
import scipy.sparse.linalg as lin
import scipy.sparse as sp
import numpy as np
import matplotlib.pylab as plt
import time
import glob
from datetime import datetime

def check_crossover():
  runs = 100
  samples = 200
  eig1 = [0] * samples
  eig2 = [0] * samples
  eig3 = [0] * samples
  f = open('eig-data-%d.txt' % (time.time() % 10000), 'w+')
  for i in range(runs):
    for k in range(1,samples+1):
      c = 10
      n = 10000
      G = fast_gnp_random_graph(n, c/float(n))
      print "done making graph"

      for i in range(1,k+1):
        G.add_edge(0, i)
      print "done making hub"

      A = to_scipy_sparse_matrix(G, dtype='d')
      print "done with sparse matrix"
      ''' Below here was how we did the e-vector centrality
      eig = lin.eigs(A, k=2, which='LR')
      '''
      #this is hashimoto style
      I = sp.identity(A.shape[0])
      I_minus_D = sp.lil_matrix(A.shape)
      for node,deg in G.degree_iter():
        I_minus_D[node,node] = 1.0-deg
      print "done making submatrices"

      crazy = sp.bmat([[None,I],[I_minus_D,A]])
      print "done making crazy"
      eig = lin.eigsh(crazy, k=3, which="LA")
      print eig[0]
      eig1[k-1] += abs(eig[0][0].real)
      eig2[k-1] += abs(eig[0][1].real)
      eig3[k-1] += abs(eig[0][2].real)
  for i in range(samples):
    eig1[i] /= float(runs)
    eig2[i] /= float(runs)
    eig3[i] /= float(runs)
    f.write(str(eig1[i]) + ' ' + str(eig2[i]) + '\n')

  f.close()
  plt.plot(eig1)
  plt.plot(eig2)
  plt.plot(eig3)
  plt.show()

def plot_localization():
  runs = 50
  samples = 200
  heights = [0] * samples
  f = open('localization-data.txt', 'w+')
  for i in range(runs):
    for k in range(1, samples+1):
      print i, k
      c = 10
      n = 10000

      G = fast_gnp_random_graph(n, c/float(n))
      for j in range(1,k+1):
        G.add_edge(0, j)
      A = to_scipy_sparse_matrix(G, dtype='d')

      #eig = lin.eigs(A, k=1, which='LR')

      #this is hashimoto style
      I = sp.identity(A.shape[0])
      I_minus_D = sp.lil_matrix(A.shape)
      for node,deg in G.degree_iter():
        I_minus_D[node,node] = 1.0-deg
      crazy = sp.bmat([[None,I],[I_minus_D,A]])

      eig = lin.eigsh(crazy, k=1, which="LA")
      height = eig[1][n]
      base = sum(map(lambda x : x*x, eig[1][n+1:2*n]))
      heights[k-1] += height**2 / base

  for i in range(samples):
    heights[i] /= float(runs)
    f.write(str(heights[i]) + '\n')
  f.close()
  #plt.ylim(0,1)
  plt.plot(heights)
  plt.show()

def mega_graph():
  dirlist = glob.glob('./eig-data-*')
  print dirlist
  eig1 = [0] * 199
  eig2 = [0] * 199
  for d in dirlist:
    for i, line in enumerate(open(d)):
      eig1[i] += float(line.split()[0])
      eig2[i] += float(line.split()[1])
  eig1 = [x / float(len(dirlist)) for x in eig1]
  eig2 = [x / float(len(dirlist)) for x in eig2]
  plt.plot(eig1)
  plt.plot(eig2)
  plt.show()

def power_law():
  runs = 16000
  samples = 50
  x = [0]*samples
  y = [0]*samples
  y2 = [0]*samples
  print runs, samples
  for k in range(runs):
    for i, gamma in enumerate(np.linspace(2.0, 3.0, samples)):
      print k, i
      powerlaw = lambda n : utils.powerlaw_sequence(n, gamma)
      seq = utils.create_degree_sequence(100000,powerlaw)
      G = configuration_model(seq)
      A = to_scipy_sparse_matrix(G, dtype='d')
      eig = lin.eigs(A, k=2, which='LR')
      #print eig[0]
      x[i] += gamma
      y[i] += eig[0][0]
      y2[i] += eig[0][1]
  x = [j / float(runs) for j in x]
  y = [j / float(runs) for j in y]
  y2 = [j / float(runs) for j in y2]
  print y
  print y2
  plt.scatter(x, y)
  plt.scatter(x, y2)
  plt.show()

def process(str_list):
  samples = 20
  data = [0]*samples
  for s in str_list:
    for i, s2 in enumerate(s.split(',')):
      data[i] += float(''.join([c for c in s2 if c in '1234567890.']))
  data = [x / float(len(data)) for x in data]
  x = np.linspace(2.3, 2.7, samples)
  print x
  print data
  plt.scatter(x, data)
  plt.show()

power_law()
#plot_localization()
#mega_graph()
#check_crossover()
exit()
#check_crossover()
#check_crossover()
#check_crossover()

process(['[(294.94739625359296+0j), (263.44188958649829+0j), (145.98087487508664+0j), (223.56558588822088+0j), (146.78448141055148+0j), (137.68663415093155+0j), (136.6380596409843+0j), (150.38724278748569+0j), (124.73681797906652+0j), (106.6800547054679+0j), (94.961111792606175+0j), (70.891874362811379+0j), (63.948317806234677+0j), (72.90887358977497+0j), (68.889269250434467+0j), (43.171358140744452+0j), (46.309960108955842+0j), (63.737041231094118+0j), (41.070838356851844+0j), (29.906044029924907+0j)]','[(210.20539286257824+0j), (270.76472519800404+0j), (252.62434460665091+0j), (152.29045508571639+0j), (188.02607555742088+0j), (142.0181001533075+0j), (122.82352550371969+0j), (100.15248341399872+0j), (121.75488998217483+0j), (97.574557859891414+0j), (119.53499950609329+0j), (55.090715631996815+0j), (59.148050483331971+0j), (52.6371939748175+0j), (62.584812012787673+0j), (36.839871425798364+0j), (78.397897222027879+0j), (47.832674954494152+0j), (75.952684927421487+0j), (30.716184018217412+0j)]','[(254.27699405318739+0j), (247.61267893664422+0j), (164.00164249891699+0j), (241.31807693162119+0j), (194.14505842441002+0j), (171.56418556948321+0j), (147.0539728027768+0j), (151.13425608432993+0j), (102.58571944481749+0j), (112.59051015939637+0j), (91.045586710108779+0j), (66.029926465684454+0j), (63.264919992517363+0j), (85.605991370541645+0j), (45.354382349422657+0j), (75.803396019657328+0j), (49.648262475823799+0j), (45.632517448933974+0j), (31.897356244448886+0j), (42.008866726839287+0j)]','[(289.41491999157017+0j), (221.35884189781058+0j), (286.55029728683036+0j), (158.22094190122849+0j), (227.91794440876433+0j), (148.67005002164862+0j), (134.66715649238691+0j), (112.55427977017713+0j), (159.4572176052414+0j), (101.2723518920769+0j), (85.387318904628458+0j), (98.706066802263408+0j), (80.203713388789339+0j), (64.664601831830325+0j), (68.743487578090424+0j), (62.126450753508813+0j), (41.101597970462478+0j), (31.269552133897868+0j), (46.292788464365877+0j), (37.031618886613515+0j)]'])


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
