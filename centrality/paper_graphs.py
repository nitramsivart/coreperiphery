from networkx import *
from custom_networkx import *
import scipy.sparse.linalg as lin
import scipy.sparse as sp
import numpy as np
import matplotlib.pylab as plt
import time
import glob
import random
from datetime import datetime
import cProfile


def check_crossover(hashi=False):
  runs = 1
  samples = 100
  eig1 = [0] * samples
  eig2 = [0] * samples
  eig3 = [0] * samples
  f = open('eig-data-%d.txt' % (time.time() % 10000), 'w+')
  hubrange = range(50,150)
  for i in range(runs):
    for index, hubsize in enumerate(hubrange):
      print i, hubsize
      c = 10
      n = 50000
      G = fast_gnp_random_graph(n, c/float(n))

      for j in range(1,hubsize+1):
        G.add_edge(0, j)

      A = to_scipy_sparse_matrix(G, dtype='d')

      if not hashi:
        eig = lin.eigsh(A, k=2, which='LA', return_eigenvectors=False)
      else:
        #this is hashimoto style
        I = sp.identity(A.shape[0])
        I_minus_D = sp.lil_matrix(A.shape)
        for node,deg in G.degree_iter():
          I_minus_D[node,node] = 1.0-deg

        crazy = sp.bmat([[None,I],[I_minus_D,A]])
        print "done making crazy"
        eig = lin.eigs(crazy, k=1, which="LR", return_eigenvectors=False)
      eig1[index] += abs(eig[0].real)
      eig2[index] += abs(eig[1].real)
      #eig3[k-1] += abs(eig[2].real)
  for i in range(samples):
    eig1[i] /= float(runs)
    eig2[i] /= float(runs)
    #eig3[i] /= float(runs)
    f.write(str(eig1[i]) + ' ' + str(eig2[i]) + '\n')

  f.close()
  #plt.ylim(6,22)
  plt.plot(hubrange, eig1)
  plt.plot(hubrange, eig2)
  #plt.plot(eig3)
  plt.show()

def plot_localization(hashi=False):
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

      if not hashi:
        eig = lin.eigsh(A, k=1, which='LA')
        height = eig[1][0]
        base = sum(map(lambda x : x*x, eig[1][1:]))
        heights[k-1] += height**2 / base
      else:
        I = sp.identity(A.shape[0])
        I_minus_D = sp.lil_matrix(A.shape)
        for node,deg in G.degree_iter():
          I_minus_D[node,node] = 1.0-deg
        crazy = sp.bmat([[None,I],[I_minus_D,A]])

        eig = lin.eigs(crazy, k=1, which="LA")
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

def localization_correlation():
  f = open('localization-correlation.txt', 'w+')

  vals1 = [{},{}]
  counts1 = [{},{}]
  vals2 = [{},{}]
  counts2 = [{},{}]
  var10 = [[[],[]],[[],[]]]
  c = 10
  n = 2000000
  G = fast_gnp_random_graph(n, c/float(n))

  eig = [None,None]
  for index, hubsize in enumerate([50, 160]):

    for j in range(1,hubsize+1):
      G.add_edge(0, j)
    A = to_scipy_sparse_matrix(G, dtype='d')
    print 'done building graph'

    eig[index] = lin.eigsh(A, k=2, which='LA')
    print eig[index][0]
    print 'done getting eigs'
    for pos,elem in enumerate(eig[index][1][:,0][hubsize+1:]):
      deg = G.degree(hubsize+1+pos)
      if deg is 10:
        var10[index][0].append(elem)
      vals1[index][deg] = vals1[index].get(deg, 0) + elem
      counts1[index][deg] = counts1[index].get(deg, 0) + 1
    for pos,elem in enumerate(eig[index][1][:,1][hubsize+1:]):
      deg = G.degree(hubsize+1+pos)
      if deg is 10:
        var10[index][1].append(elem)
      vals2[index][deg] = vals2[index].get(deg, 0) + elem
      counts2[index][deg] = counts2[index].get(deg, 0) + 1
    print 'done counting'
    for key, val in counts1[index].iteritems():
      vals1[index][key] /= float(val)
    for key, val in counts2[index].iteritems():
      vals2[index][key] /= float(val)

  for i, (e1, e2, e3, e4) in enumerate(zip(eig[0][1][:,0], eig[1][1][:,0],
                                           eig[0][1][:,1], eig[1][1][:,1])):
    deg = G.degree(i)
    f.write("%d, %f, %f, %f, %f\n" % (deg,e1,e2,e3,e4))
  f.close()
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.spines['top'].set_color('none')
  ax.spines['bottom'].set_color('none')
  ax.spines['left'].set_color('none')
  ax.spines['right'].set_color('none')
  ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
  
  
  ax1 = fig.add_subplot(221)
  ax1.set_title('2nd evect, no localization, std = %f'%np.std(var10[0][0]))
  plt.plot(vals1[0].keys(), vals1[0].values())
  ax2 = fig.add_subplot(222)
  ax2.set_title('2nd evect, with localization, std = %f'%np.std(var10[1][0]))
  plt.plot(vals1[1].keys(), vals1[1].values())
  ax3 = fig.add_subplot(223)
  ax3.set_title('1st evect, no localization, std = %f'%np.std(var10[0][1]))
  plt.plot(vals2[0].keys(), vals2[0].values())
  ax4 = fig.add_subplot(224)
  ax4.set_title('1st evect, with localization, std = %f'%np.std(var10[1][1]))
  ax.set_xlabel('node degree')
  plt.plot(vals2[1].keys(), vals2[1].values())
  plt.show()
  print counts2[0]

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

def power_law(evalues=True,localization=False):
  runs = 1000
  samples = 40
  nodes = 50000
  print runs, samples, nodes
  y = [0]*samples
  y2 = [0]*samples
  heights = [0] * samples
  gammarange = np.linspace(2.0, 3.0, samples)
  for runnum in range(runs):
    for i, gamma in enumerate(gammarange):
      print runnum, i
      #powerlaw = lambda n : utils.powerlaw_sequence(n, gamma)
      #seq = utils.create_degree_sequence(nodes,powerlaw)
      seq = custom_create_degree_sequence(nodes, gamma)
      print 'done w degree seq'
      G = custom_configuration_model(seq)
      A = to_scipy_sparse_matrix(G, dtype='d')
      print 'done making graph'
      eig = lin.eigsh(A, k=1, which='LA')
      print 'done computing eigs'
      if localization:
        height = max(eig[1])
        base = sum(eig[1])
        #base = reduce(lambda x,y: x+y*y, eig[1])
        heights[i] += height/base
        #heights[i] += height**2 / base
      if evalues:
        y[i] += eig[0][0]
        y2[i] += eig[0][1]
  if localization:
    heights = [j / float(runs) for j in heights]
    print heights
    plt.scatter(gammarange, heights)
    plt.show()
  if evalues:
    y = [j / float(runs) for j in y]
    y2 = [j / float(runs) for j in y2]
    print y
    print y2
    plt.scatter(gammarange, y)
    plt.scatter(gammarange, y2)
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


#cProfile.run('power_law(True,True)')
#power_law(False,True)
#plot_localization()
#mega_graph()
#check_crossover()
localization_correlation()
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
