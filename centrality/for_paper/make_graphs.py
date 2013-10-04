import matplotlib.pylab as plt
import numpy as np
import cProfile
import networkx as nx

def eig_plots():
  string = "data/mega-data-10-1000000.txt"
  hubsize = []
  eig1 = []
  eig2 = []
  eigh = []
  op = []
  oph = []
  nb = []
  nbh = []
  for l in open(string):
    vals = [abs(complex(s).real) for s in l.split()]
    hubsize.append(vals[0])
    eig1.append(vals[1])
    eig2.append(vals[2])
    eigh.append(vals[3])
    op.append(vals[4])
    oph.append(vals[5])
    nb.append(vals[6])
    nbh.append(vals[7])

  plt.plot(hubsize, eig1, color='k', ls='-',linewidth=5.0)
  plt.plot(hubsize, eig2, color='g', ls=':',linewidth=5.0)
  plt.plot(hubsize, eigh, color='r', ls='--',linewidth=5.0)
  plt.savefig("eig.eps", bbox_inches="tight")
  plt.figure()
  plt.plot(hubsize, op, color='k', ls='-',linewidth=5.0)
  plt.plot(hubsize, oph, color='r', ls='--',linewidth=5.0)
  plt.savefig("op.eps", bbox_inches="tight")
  plt.figure()
  plt.plot(hubsize, nb, color='k', ls='-',linewidth=5.0)
  plt.plot(hubsize, nbh, color='r', ls='--',linewidth=5.0)
  plt.savefig("nb.eps", bbox_inches="tight")
  plt.show()

def power_plots():
  filename = 'data/power_law.txt'
  gamma_range = list(np.linspace(2.2, 2.8, 100))
  total_weight = 0.
  power_adj = [0.] * 100
  power_hashi = [0.] * 100
  for index, line in enumerate(open(filename)):
    if index % 3 == 0:
      total_weight += int(line)
    elif index % 3 == 1:
      power_adj = [x+y for x, y in zip(eval(line), power_adj)]
    else:
      power_hashi = [x+y for x, y in zip(eval(line), power_hashi)]
  power_adj = [(x/total_weight)**(1./4.) for x in power_adj]
  power_hashi = [(x/total_weight)**(1./4.) for x in power_hashi]
  plt.plot(gamma_range, power_adj)
  plt.figure()
  plt.plot(gamma_range, power_hashi)
  plt.show()

def localization(name):
  position = eval(open('data/positions.txt').readline())
  #print pos
  G = nx.read_adjlist('data/%s.adjlist'%name)
  mapping = {str(i):i for i in range(1000000)}                                    
  print len(G.nodes())
  H = nx.relabel_nodes(G, mapping)
  nodes = H.nodes()
  #print nodes
  '''
  for n in nodes:
    if n % 20 != 0:
      H.remove_node(n)
  '''
  counts = []
  eig = []
  eigh = []
  #f = open('%s-trimmed.txt'%name, 'w+')
  for l in open('data/%s-trimmed.txt'%name):
    strings = l.split()
    counts.append(int(strings[0]))
    eig.append(27000*float(strings[1]))
    eigh.append(27000*float(strings[2]))

  max_val = max(eig)
  max_valh = max(eigh)
  difference = max_val - min(eig)
  differenceh = max_valh - min(eigh)
  colors = []
  shapes = []
  transparency = []
  transparencyh = []
  plt.figure(figsize=(5,5))
  '''
  max_cut = 1.05
  min_cut = 0
  xmax= max_cut*max(xx for xx,yy in position.values())
  xmin= min_cut*min(xx for xx,yy in position.values())
  ymax= max_cut*max(yy for xx,yy in position.values())
  ymin= min_cut*min(yy for xx,yy in position.values())
  plt.xlim(xmin,xmax)
  plt.ylim(ymin,ymax)
  '''

  nx.draw(H,edgelist=H.edges(), linewidths=.01, width=.1, pos=position,
          nodelist=[0], node_shape='o', node_size=eig[0],
          node_color='k', with_labels=False)

  for index, count in enumerate(counts):
    if count in H.neighbors(0):
      nx.draw(H,edgelist=[], linewidths=.01, width=.1, pos=position,
              nodelist=[count], node_shape='o', node_size=eig[index],
              node_color='r', alpha=.65, with_labels=False)

  for index, count in enumerate(counts):
    if count not in H.neighbors(0) and count != 0:
      nx.draw(H,edgelist=[], linewidths=.01, width=.1, pos=position,
              nodelist=[count], node_shape='o', node_size=eig[index],
              node_color='b', alpha=.65, with_labels=False)

  print "done drawing"
  plt.savefig("%s.png"%name, dpi=100, bbox_inches='tight')
#power_plots()
#exit()
#cProfile.run('localization("1")')
localization("1")
localization("2")
localization("3")
#plt.show()
