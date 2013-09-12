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
  runs = 1
  samples = 200
  heights = [0] * samples
  f = open('localization-data.txt', 'w+')
  for i in range(runs):
    c = 100
    n = 100000
    G = fast_gnp_random_graph(n, c/float(n))
    for k in range(140, samples+1):
      print i, k

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

        eig = lin.eigs(crazy, k=1, which="LR")
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

def mega():
  runs = 1
  min_hub = 1
  max_hub = 200
  samples = max_hub - min_hub + 1
  eval1 = [0] * samples
  eval2 = [0] * samples
  evalh = [0] * samples
  op = [0] * samples
  oph = [0] * samples
  nb = [0] * samples
  nbh = [0] * samples
  c = 7
  n = 1000
  f = open('mega-data-%d-%d.txt' % (c, n), 'w+')
  print c, n
  for dummy_counter in range(runs):
    G = fast_gnp_random_graph(n, c/float(n))
    for k in range(1, min_hub):
      G.add_edge(0, k)
    A = to_scipy_sparse_matrix(G, dtype='d',format='csr')
    for i, k in enumerate(range(min_hub, max_hub+1)):
      print dummy_counter, k
      G.add_edge(0, k)
      A[0, k] = 1
      A[k, 0] = 1

      #non hashi
      eig = lin.eigsh(A, k=2, which='LA')
      eval1[i] += eig[0][1]
      eval2[i] += eig[0][0]
      op[i] += abs(eig[1][0,1])
      total = 0.
      for j in range(1,k+1):
        total += abs(eig[1][j,1])
      nb[i] += total/k

      #hashi
      I = sp.identity(A.shape[0])
      I_minus_D = sp.lil_matrix(A.shape)
      for node,deg in G.degree_iter():
        I_minus_D[node,node] = 1.0-deg
      crazy = sp.bmat([[None,I],[I_minus_D,A]])

      eig = lin.eigs(crazy, k=1, which="LR")
      evalh[i] += eig[0][0]
      oph += abs(eig[1][0])
      total = 0.
      for j in range(1,k+1):
        total += abs(eig[1][0+j])
      nbh[i] += total/k

  x = range(min_hub, max_hub+1)
  plt.figure()
  plt.plot(x, eval1)
  plt.plot(x, eval2)
  plt.plot(x, evalh)
  plt.figure()
  plt.plot(x, evalh)
  plt.figure()
  plt.plot(x, op)
  plt.figure()
  plt.plot(x, oph)
  plt.figure()
  plt.plot(x, nb)
  plt.figure()
  plt.plot(x, nbh)
  plt.show()

  for i in range(samples):
    f.write("%d,%r,%r,%r,%r,%r,%r,%r\n" % 
            (i,eval1[i],eval2[i],evalh[i],op[i],oph[i],nb[i],nbh[i]))
  f.close()
                                      
def interact():
  def dummy(adj, graph, position, colors):
    eig = lin.eigsh(adj, k=1, which='LA')
    l = map(lambda x:abs(x*x*40000),eig[1][:,0])
    plt.figure(figsize=(13,13))
    networkx.draw(graph, linewidths=.1,pos=position, node_color = colors, width=.05, node_size = l, with_labels=False)

  min_hub = 25
  max_hub = 100
  c = 7
  n = 1000
  f = open('interact-data-%d-%d.txt' % (c, n), 'w+')
  print c, n
  G = fast_gnp_random_graph(n, c/float(n))
  colors = ['r'] * n
  colors2 = ['r'] * n
  colors[0] = colors2[0] = 'b'
  #pos = graphviz_layout(G, prog='sfdp', root=0)
  for k in range(1, min_hub+1):
    G.add_edge(0, k)
    colors[k] = 'g'
    colors2[k] = 'g'

  G2 = G.copy()
  for k in range(min_hub+1, max_hub+1):
    G2.add_edge(0, k)
    colors2[k] = 'g'
  pos = graphviz_layout(G2, prog='neato', root=0)
  A = to_scipy_sparse_matrix(G, dtype='d')
  A2 = to_scipy_sparse_matrix(G2, dtype='d')
  dummy(A,G,pos,colors)
  dummy(A2,G2,pos,colors2)
  plt.show()

def movie():

  min_hub = 1
  max_hub = 200
  c = 7
  n = 400000
  print c, n
  G = fast_gnp_random_graph(n, c/float(n))
  print 'done w graph'
  colors = ['r'] * n
  colors[0] = 'b'
  #pos = graphviz_layout(G, prog='sfdp', root=0)
  for k in range(1, max_hub+1):
    G.add_edge(0, k)
    colors[k] = 'g'
  #pos = graphviz_layout(G, prog='neato', root=0)
  xmax = None
  ymax = None
  for hubsize in reversed(range(min_hub, max_hub+1)):
    print hubsize
    if hubsize != max_hub:
      G.remove_edge(0, hubsize+1)
      colors[hubsize+1] = 'r'
    A = to_scipy_sparse_matrix(G, dtype='d')
    eig = lin.eigsh(A, k=1, which='LA')
    print 'done with eig'
    #l = map(lambda x:abs(x*x*80000),eig[1][:,0])
    l = map(lambda x:abs(x), eig[1][:,0])
    deg = [d for x,d in G.degree_iter()]
    plt.figure(figsize=(13,13))
    plt.scatter(deg, l, c=colors)
    #networkx.draw(G, linewidths=.02,pos=pos, node_color = colors, width=.01, node_size = l, with_labels=False)
    plt.title(str(hubsize))
    if hubsize == max_hub:
      xmax = G.degree(0)
      ymax = abs(eig[1][0,0])
    plt.axis([0,xmax,0,ymax])
    plt.savefig('img-400000-7/%d'%hubsize)

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

#mega()
#cProfile.run('interact()')
#interact()
movie()
#cProfile.run('power_law(True,True)')
#power_law(False,True)
#plot_localization()
#mega_graph()
#check_crossover()
#localization_correlation()
exit()
#check_crossover()
#check_crossover()
#check_crossover()



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
