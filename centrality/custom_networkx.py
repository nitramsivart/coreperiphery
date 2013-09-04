from networkx import *
import random
import scipy.sparse as sp

def custom_create_degree_sequence(nodes, gamma):
  max_tries = 50
  tries=0
  max_deg=nodes
  while tries < max_tries:
    trialseq=utils.powerlaw_sequence(nodes, gamma)
    # round to integer values in the range [0,max_deg]
    seq=[min(max_deg, max( int(round(s)),0 )) for s in trialseq]
    # if graphical return, else throw away and try again
    #if is_valid_degree_sequence(seq, 'eg'):
    if sum(seq)%2 == 0:
      return seq
    tries+=1
  raise NetworkXError("Exceeded max (%d) attempts at a valid sequence."%max_tries)

def custom_configuration_model(deg_sequence):
  create_using = Graph()

  # start with empty N-node graph
  N=len(deg_sequence)

  # allow multiedges and selfloops
  #G=sp.dok_matrix((N,N),dtype='d')
  G=empty_graph(N,create_using)

  # build stublist, a list of available degree-repeated stubs
  # e.g. for deg_sequence=[3,2,1,1,1]
  # initially, stublist=[1,1,1,2,2,3,4,5]
  # i.e., node 1 has degree=3 and is repeated 3 times, etc.
  stublist=[]
  for n in range(N):
    for i in range(deg_sequence[n]):
      stublist.append(n)

  # shuffle stublist and assign pairs by removing 2 elements at a time
  random.shuffle(stublist)
  while stublist:
    n1 = stublist.pop()
    n2 = stublist.pop()
    #G[n1,n2] = 1
    #G[n2,n1] = 1
    G.add_edge(n1,n2)
  return G#sp.coo_matrix(G)
