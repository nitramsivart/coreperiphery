from networkx import *
from custom_networkx import *
import scipy.sparse.linalg as lin
import scipy.sparse as sp
import numpy as np
import matplotlib.pylab as plt
import time
import glob
import cProfile

def hashi_eig(A, degrees):
  I = sp.identity(A.shape[0])
  I_minus_D = sp.lil_matrix(A.shape)
  deg_list = list(degrees)
  for node,deg in deg_list:
    I_minus_D[node,node] = 1.0-deg
  crazy = sp.bmat([[None,I],[I_minus_D,A]])
  eig = lin.eigs(crazy, k=1, which="LR")[1][:len(deg_list),0]
  root_total = np.sqrt(sum(x*x for x in eig))
  return [x/root_total for x in eig]

def check_crossover(hashi=False):
  runs = 1
  samples = 100
  eig1 = [0] * samples
  eig2 = [0] * samples
  f = open('eig-data-%d.txt' % (time.time() % 10000), 'w+')
  hubrange = range(50,150)
  for i in range(runs):
    for index, hubsize in enumerate(hubrange):
      print i, hubsize
      c = 10
      n = 50000
      G = fast_gnp_random_graph(n, c/float(n))
      for nbr in list(G[0]):
        G.remove_edge(0,nbr)

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
    for nbr in list(G[0]):
      G.remove_edge(0,nbr)
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
        eig = hashi_eig(A, G.degree_iter())
        height = eig[0]
        base = sum(map(lambda x : x*x, eig))
        heights[k-1] += height**2 / base

  for i in range(samples):
    heights[i] /= float(runs)
    f.write(str(heights[i]) + '\n')
  f.close()
  #plt.ylim(0,1)
  plt.plot(heights)
  plt.show()

def mega(num=2000000):
  runs = 1
  min_hub = 1
  max_hub = 160
  samples = max_hub - min_hub + 1
  eval1 = [0] * samples
  eval2 = [0] * samples
  evalh = [0] * samples
  op = [0] * samples
  oph = [0] * samples
  nb = [0] * samples
  nbh = [0] * samples
  c = 10
  n = num
  f = open('for_paper/mega-data-%d-%d.txt' % (c, n), 'w+')
  print c, n
  for dummy_counter in range(runs):
    G = fast_gnp_random_graph(n, c/float(n))
    for nbr in list(G[0]):
      G.remove_edge(0,nbr)
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
      op[i] += sum(map(lambda x : x**4, eig[1][:,1]))**(1./4.)
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
      root_total = np.sqrt(sum(x*x for x in eig[1][:n,0]))
      eig = [x/root_total for x in eig[1][:n,0]]
      oph[i] += sum(map(lambda x : x**4, eig))**(1./4.)
      total = 0.
      for j in range(1,k+1):
        total += abs(eig[j])
      nbh[i] += total/k

  for i in range(samples):
    f.write("%d %r %r %r %r %r %r %r\n" % 
            (i,eval1[i],eval2[i],evalh[i],op[i],oph[i],nb[i],nbh[i]))
  f.close()
  return

  x = range(min_hub, max_hub+1)
  plt.figure()
  plt.plot(x, eval1)
  plt.plot(x, eval2)
  plt.plot(x, evalh)
  plt.savefig("for_paper/eig.png")
  plt.figure()
  plt.plot(x, evalh)
  plt.savefig("for_paper/eigh.png")
  plt.figure()
  plt.plot(x, op)
  plt.savefig("for_paper/op.png")
  plt.figure()
  plt.plot(x, oph)
  plt.savefig("for_paper/oph.png")
  plt.figure()
  plt.plot(x, nb)
  plt.savefig("for_paper/nb.png")
  plt.figure()
  plt.plot(x, nbh)
  plt.savefig("for_paper/nbh.png")
  #plt.show()
  return eval1, eval2, evalh, op, oph, nb, nbh

def write_to_file(num=2000000):
  min_hub = 20
  middle_hub = 70
  max_hub = 120
  c = 10
  n = num
  print c, n
  G = fast_gnp_random_graph(n, c/float(n))
  for nbr in list(G[0]):
    G.remove_edge(0,nbr)
  print 'done making graph'

  for k in range(1, min_hub+1):
    G.add_edge(0, k)
  G1 = G.copy()
  A1 = to_scipy_sparse_matrix(G, dtype='d')

  
  eig1 = lin.eigsh(A1, k=1, which='LA')[1][:,0]
  eig1h = hashi_eig(A1, G.degree_iter())
  print "done edges 1"

  for k in range(min_hub+1, middle_hub+1):
    G.add_edge(0, k)
  G2 = G.copy()
  A2 = to_scipy_sparse_matrix(G, dtype='d')
  
  eig2 = lin.eigsh(A2, k=1, which='LA')[1][:,0]
  eig2h = hashi_eig(A2, G.degree_iter())
  print "done edges 2"

  for k in range(middle_hub+1, max_hub+1):
    G.add_edge(0, k)
  A3 = to_scipy_sparse_matrix(G, dtype='d')

  eig3 = lin.eigsh(A3, k=1, which='LA')[1][:,0]
  eig3h = hashi_eig(A3, G.degree_iter())
  print "done edges 3"

  f = open('for_paper/all-eigs-%d-%d.txt' % (c, n), 'w+')
  for i in range(len(eig3)):
    d = tuple(map(lambda x:abs(x.real),(G.degree(i),eig1[i], eig2[i], eig3[i], eig1h[i], eig2h[i], eig3h[i])))
    f.write("%d\t%f\t%f\t%f\t%f\t%f\t%f\n" % d)
  f.close()

  '''
  write_dot(G1, "G1-%d-%d" % (c,n))
  write_dot(G2, "G2-%d-%d" % (c,n))
  write_dot(G, "G3-%d-%d" % (c,n))
  '''
                                      
def interact(num=2000000):
  def dummy(adj, graph, degrees, position, colors, nodes, name):
    print len(degrees)
    write_adjlist(graph,"for_paper/%s.adjlist"%name)
    eig = lin.eigsh(adj, k=1, which='LA')[1][:,0]
    #sizes = [x**2*1000000 for x in eig[1][:,0]]

    eigh = hashi_eig(adj, degrees)
    f = open("for_paper/%s-eigs.txt"%name, "w+")
    for e, eh in zip(eig, eigh):
      f.write("%r %r\n" % (abs(e), abs(eh.real)))
    f.close()
    #sizesh = [x.real**2*1000000 for x in eigh]
    #plt.figure(figsize=(7,7))
    print 'done finding eig'
    #networkx.draw(graph, linewidths=.1,pos=position, node_color = colors, width=.05, node_size = sizes, with_labels=False)
    #plt.savefig("%s.png"%name)

    #plt.figure(figsize=(7,7))
    #networkx.draw(graph, linewidths=.1,pos=position, node_color = colors, width=.05, node_size = sizesh, with_labels=False)
    #plt.savefig("%s-h.png"%name)

  min_hub = 20
  middle_hub = 70
  max_hub = 120
  c = 10
  n = num
  #f = open('interact-data-%d-%d.txt' % (c, n), 'w+')
  print c, n
  G = fast_gnp_random_graph(n, c/float(n))
  for nbr in list(G[0]):
    G.remove_edge(0,nbr)
  print 'done making graph'
  colors = ['r'] * n
  colors2 = ['r'] * n
  colors3 = ['r'] * n
  colors[0] = colors2[0] = colors3[0] = 'b'
  for k in range(1, min_hub+1):
    G.add_edge(0, k)
  for k in G.neighbors(0):
    colors[k] = 'g'
    colors2[k] = 'g'
    colors3[k] = 'g'
    for j in G.neighbors(k):
      if colors[j] == 'r':
        colors[j] = 'red'
        colors2[j] = 'red'
        colors3[j] = 'red'
  degrees1 = list(G.degree_iter())

  G2 = G.copy()
  for k in range(min_hub+1, middle_hub+1):
    G2.add_edge(0, k)
    colors2[k] = 'g'
    colors3[k] = 'g'
    for j in G2.neighbors(k):
      if colors2[j] == 'r':
        colors2[j] = 'red'
        colors3[j] = 'red'
  degrees2 = list(G2.degree_iter())
  G3 = G2.copy()
  for k in range(middle_hub+1, max_hub+1):
    G3.add_edge(0, k)
    colors3[k] = 'g'
    for j in G3.neighbors(k):
      if colors3[j] == 'r':
        colors3[j] = 'red'
  degrees3 = list(G3.degree_iter())

  print 'done adding edges'
  #nodes = filter(lambda n : colors3[n] != 'r', G.nodes())

  A = to_scipy_sparse_matrix(G, dtype='d')
  A2 = to_scipy_sparse_matrix(G2, dtype='d')
  A3 = to_scipy_sparse_matrix(G3, dtype='d')
  print 'done sending to matrix'

  #return G, G2, G3, A, A2, A3

  count = 0
  final_nodes = {}
  final_colors1 = []
  final_colors2 = []
  final_colors3 = []
  for v1 in G3.neighbors(0):
    for v2 in G3.neighbors(v1):
      if colors3[v2] == 'r':
        colors3[v2] = 'red'
      #for v3 in G3.neighbors(v2):
        #if colors3[v3] == 'r':
          #colors3[v3] = 'red'
  for n in range(len(colors3)):
    if colors3[n] != 'r':
      final_nodes[count] = n
      final_colors1.append(colors[n])
      final_colors2.append(colors2[n])
      final_colors3.append(colors3[n])
      count += 1
    else:
      G.remove_node(n)
      G2.remove_node(n)
      G3.remove_node(n)
  print len(G3.nodes())
  pos = graphviz_layout(G3)
  f = open("for_paper/positions.txt", "w+")
  f.write(str(pos))
  f.close()

  print 'done positioning'
  dummy(A,G,degrees1,pos,final_colors1,nodes,"1")
  dummy(A2,G2,degrees2,pos,final_colors2,nodes,"2")
  dummy(A3,G3,degrees3,pos,final_colors3,nodes,"3")
  plt.show()

def movie():
  min_hub = 10
  max_hub = 100
  c = 6
  n = 500000
  print c, n
  G = fast_gnp_random_graph(n, c/float(n))
  for nbr in list(G[0]):
    G.remove_edge(0,nbr)
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
    
    #hashi
    '''
    I = sp.identity(A.shape[0])
    I_minus_D = sp.lil_matrix(A.shape)
    for node,deg in G.degree_iter():
      I_minus_D[node,node] = 1.0-deg
    crazy = sp.bmat([[None,I],[I_minus_D,A]])
    eig = lin.eigs(crazy, k=1, which="LR")
    root_total = sqrt(sum(x*x for x in eigh)) # needs some editing
    eigh = [x/root_total for x in eigh]
    '''
    eig = lin.eigsh(A, k=1, which='LA')
    print 'done with eig'
    #l = map(lambda x:abs(x*x*80000),eig[1][:,0])
    l = map(lambda x:abs(x), eig[1][:n,0])
    print len(l)
    deg = [d for x,d in G.degree_iter()]
    plt.figure(figsize=(7,7))
    plt.scatter(deg, l, c=colors)
    #networkx.draw(G, linewidths=.02,pos=pos, node_color = colors, width=.01, node_size = l, with_labels=False)
    plt.title(str(hubsize))
    if hubsize == max_hub:
      xmax = G.degree(0)
      ymax = abs(eig[1][0,0])
    plt.axis([0,xmax,0,ymax])
    plt.savefig('img-500000-6/%d'%hubsize)

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
  for nbr in list(G[0]):
    G.remove_edge(0,nbr)

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

def simplest_power_law(num=100000):
  runs = 10
  samples = 100
  nodes = num
  print runs, samples, nodes
  heights = [0] * samples
  heightsh = [0] * samples
  gammarange = np.linspace(2.2, 2.8, samples)
  for runnum in range(runs):
    print heights
    for i, gamma in enumerate(gammarange):
      print runnum, i
      seq = custom_create_degree_sequence(nodes, gamma)
      root_total = sum(x**2 for x in seq)**(.5)
      #G = custom_configuration_model(seq)
      heights[i] += sum((x/root_total)**4 for x in seq)
  heights = [(j / float(runs))**(.25) for j in heights]
  print heights
  plt.plot(gammarange, heights)
  plt.show()

def simpler_power_law(num=2000000):
  runs = 1000
  samples = 100
  nodes = num
  print runs, samples, nodes
  heights = [0] * samples
  heightsh = [0] * samples
  gammarange = np.linspace(2.2, 2.8, samples)
  for runnum in range(runs):
    if True or runnum % 10 == 0:
      print heights
      print heightsh
    for i, gamma in enumerate(gammarange):
      print runnum, i
      seq = custom_create_degree_sequence(nodes, gamma)
      G = custom_configuration_model(seq)
      A = to_scipy_sparse_matrix(G, dtype='d')
      eig = lin.eigsh(A, k=1, which='LA')[1][:,0]
      # hashi
      eigh = hashi_eig(A, G.degree_iter())
      heights[i] += sum(x**4 for x in eig)
      heightsh[i] += sum(x.real**4 for x in eigh)
  f = open("for_paper/power-law-%d.txt" % num, 'w+')
  heights = [j / float(runs) for j in heights]
  heightsh = [j / float(runs) for j in heightsh]
  print heights
  print heightsh
  f.write("%d %d %d\n" %(runs, samples, nodes))
  for i in range(samples):
    f.write("%r %r" % (heights[i], heightsh[i]))
  f.close()

def power_law(evalues=True,localization=False, num = 10000):
  runs = 1000
  samples = 100
  nodes = num
  print runs, samples, nodes
  y = [0]*samples
  y2 = [0]*samples
  heights = [0] * samples
  gammarange = np.linspace(2.2, 2.8, samples)
  for runnum in range(runs):
    print heights
    for i, gamma in enumerate(gammarange):
      print runnum, i
      #powerlaw = lambda n : utils.powerlaw_sequence(n, gamma)
      #seq = utils.create_degree_sequence(nodes,powerlaw)
      seq = custom_create_degree_sequence(nodes, gamma)
      G = custom_configuration_model(seq)
      A = to_scipy_sparse_matrix(G, dtype='d')
      eig = lin.eigsh(A, k=1, which='LA')
      print 'done computing eigs'
      if localization:
        heights[i] += sum(map(lambda x : x**4, eig[1][:,0]))
        '''
        height = max(eig[1])
        base = sum(eig[1])
        #base = reduce(lambda x,y: x+y*y, eig[1])
        heights[i] += height/base
        #heights[i] += height**2 / base
        '''
      if evalues:
        y[i] += eig[0][0]
        y2[i] += eig[0][1]
  if localization:
    heights = [j / float(runs) for j in heights]
    print heights
    #plt.scatter(gammarange, heights)
    #plt.show()
  if evalues:
    y = [j / float(runs) for j in y]
    y2 = [j / float(runs) for j in y2]
    print y
    print y2
    plt.scatter(gammarange, y)
    plt.scatter(gammarange, y2)
    plt.show()
  return heights

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


#interact(1000000)
#mega(1000000)
#write_to_file(1000000)
simplest_power_law(1000000)
#mega(2000000)
#cProfile.run('interact()')
#interact()
#write_to_file()
#movie()
#cProfile.run('power_law(True,True)')
#power_law(False,True)
#plot_localization()
#mega_graph()
#check_crossover()
#localization_correlation()
