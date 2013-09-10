# by Travis Martin

import numpy as np
import itertools as it
import random
import operator
import pdb
from math import exp

# gives the probability of having the edge (or lack of edge) adj
# and the block probability connection omega
def edge_prob(adj, omega, n):
  #assert 0 <= omega <= 1, omega
  assert adj == 1 or adj == 0, adj
  #return (omega**adj)*np.exp(-omega)
  return omega if adj == 1. else 1.-omega


# adj is the adjacency matrix of the graph
# messages is an N x N x [#types] array, giving marginal type probs
# omega is the block affinity matrix
# gamma is the prior on probability of being in a certain group
def update_messages(adj, messages, omega, gamma):
  new_messages = np.copy(messages)
  n_range = xrange(len(messages)) # number of nodes
  k_range = xrange(len(messages[0][0])) # number of types
  l = list(it.product(n_range, n_range))
  random.shuffle(l)
  for u, v in l:
    # here we calculate the message passed from u->v, for each type r
    # gives probability that u is type r, without the presence of v
    if u == v:
      continue

    norm_total = 0.
    prod_total = 0.
    sum_total = 0.
    for r in k_range:
      prod_total = 0.
      for w in n_range:
        if w == u or w == v:
          continue

        sum_total = 0.
        for s in k_range:
          assert 0 <= messages[w][u][s] <= 1, messages[w][u][s]
          sum_total += messages[w][u][s] * edge_prob(adj[w][u], omega[r][s],len(n_range))
        #assert 0 <= sum_total <= 1
        prod_total += np.log(sum_total)
      new_messages[u][v][r] = gamma[r]*np.exp(prod_total)
      norm_total += new_messages[u][v][r]
    for r in k_range:
      new_messages[u][v][r] = new_messages[u][v][r] / norm_total
      #messages[u][v][r] = new_messages[u][v][r]

    assert 0 <= new_messages[u][v][0] <= 1, pdb.set_trace()
  return new_messages

def get_message_field(adj, messages, omega, gamma):
  n = len(messages)
  k = len(messages[0][0])
  message_field = np.zeros([n, k])
  for u, r in it.product(range(n), range(k)):
    prod_total = 0.
    for w in range(n):
      if w == u:
        continue

      sum_total = (messages[w][u][0] * edge_prob(adj[w][u], omega[r][0], n) + 
                   messages[w][u][1] * edge_prob(adj[w][u], omega[r][1], n))

      prod_total += np.log(sum_total)
    message_field[u][r] = gamma[r]*np.exp(prod_total)

  # normalize
  for u in range(n):
    sum_total = sum(message_field[u])
    for r in range(k):
      message_field[u][r] /= sum_total
      assert 0 <= message_field[u][r] <= 1
  return message_field


def get_joint_dist(adj, messages, omega):
  n = len(messages)
  k = len(messages[0][0])
  joint_dist = np.zeros([n, n, k, k])
  for u,v,r,s in it.product(range(n), range(n), range(k), range(k)):
    if u == v:
      continue
    joint_dist[u][v][r][s] = (edge_prob(adj[u][v], omega[r][s], n) *
                              messages[u][v][r] * messages[v][u][s])

  # normalize
  for u, v in it.product(range(n), range(n)):
    if u==v:
      continue
    sum_total = sum(sum(joint_dist[u][v]))
    for r, s in it.product(range(k), range(k)):
      joint_dist[u][v][r][s] /= sum_total
  return joint_dist


def update_parameters(adj, messages, message_field, joint_dist, omega_old):
  n = len(messages)
  k = len(messages[0][0])
  omega = np.zeros([k,k])
  gamma = np.zeros([k])
  # update gamma
  for r in range(k):
    for u in range(n):
      gamma[r] += message_field[u][r]
    # normalize gamma
    gamma[r] /= n

    assert gamma[r] <= 1, gamma[r]

  #calculate normalizing constants according to decelle
  z = np.zeros([n,n])
  for u, v in it.product(range(n), range(n)):
    sum_total = 0.
    for b in range(k):
      for a in range(b):
        sum_total += omega_old[a][b] * (messages[u][v][a] * messages[v][u][b] + messages[u][v][b]*messages[v][u][a])
    for a in range(k):
      sum_total += omega_old[a][a] * messages[u][v][a] * messages[v][u][a]
    z[u][v] += sum_total

  pdb.set_trace()

  # update omega
  for r in range(k):
    for s in range(k):
      for u, v in it.product(range(n), range(n)):
        if u == v or adj[u][v] == 0:
          continue
        omega[r][s] += omega_old[r][s] * (messages[u][v][r] * messages[v][u][s] + messages[u][v][s]*messages[v][u][r]) / z[u][v]
      omega[r][s] /= (gamma[r] * gamma[s] * n)
  print sum(sum(omega))

  return omega, gamma

def update_parameters2(adj, messages, message_field, joint_dist, omega_old):
  print "start update_parameters2"
  n = len(messages)
  k = len(messages[0][0])
  omega = np.zeros([k,k])
  gamma = np.zeros([k])

  # update gamma
  for r in range(k):
    for u in range(n):
      gamma[r] += message_field[u][r]
    # normalize gamma
    gamma[r] /= n

    assert gamma[r] <= 1, gamma[r]

  degrees = np.zeros([k])
  #calculate degrees by type:
  for r in range(k):
    sum_total = 0.
    for u in range(n):
      degrees[r] += sum(adj[u]) * message_field[u][r]
      sum_total += message_field[u][r]
    degrees[r] /= sum_total


  # update omega
  for r in range(k):
    for s in range(k):
      for u, v in it.product(range(n), range(n)):
        if u == v or adj[u][v] == 0:
          continue
        omega[r][s] += adj[u][v] * joint_dist[u][v][r][s]
      print r, s, omega[r][s], degrees[r], degrees[s]
      omega[r][s] /= (degrees[r] * degrees[s])
      omega[r][s] *= n
  print sum(sum(omega))
  print "end update_parameters2"

  return omega, gamma

def make_graph(n1, n2, p11, p12, p22):
  size = n1 + n2
  mat = np.zeros((size, size))
  for i,j in it.product(range(size), range(size)):
    if i == j:
      continue
    if i < n1 and j < n1:
      mat[i][j] = mat[j][i] = 1 if random.random() < p11 else 0
    elif i < n1 or j < n1:
      mat[i][j] = mat[j][i] = 1 if random.random() < p12 else 0
    else:
      mat[i][j] = mat[j][i] = 1 if random.random() < p22 else 0

  return mat

# returns a graph where all n1 and n2 nodes
# are indistinguishable in terms of degree
def fixed_deg_graph(n1, n2, p11, p12, p22):
  size = n1+n2
  mat = np.zeros((size,size))
  for i in range(size):
    for j in range(size):
      return None
      


def make_complete_graph(nodes):
  mat = np.zeros((nodes, nodes))
  for i, j in it.product(range(nodes), range(nodes)):
    mat[i][j] = 1 if i != j else 1
  return mat

def bp(gamma, omega, adj, tmax):
  nodes = len(adj)
  types = len(gamma)
  messages = np.zeros([nodes,nodes,types])
  print omega
  for u,v in it.product(range(nodes), range(nodes)):
    
    a = .9 if u < 15 or u in [17,18,20,22,31] else .1
    #a = .1 if u in [0,1,2,32,33] else .5
    #a = .51
    messages[u][v][0] = a
    messages[u][v][1] = 1-a
  for i in range(tmax):
    messages = update_messages(adj, messages, omega, gamma)
    joint_dist = get_joint_dist(adj, messages, omega)
    message_field = get_message_field(adj, messages, omega, gamma)
    omega, gamma = update_parameters(adj, messages, message_field, joint_dist, omega)
    omega2, gamma2 = update_parameters2(adj, messages, message_field, joint_dist, omega)
    #print message_field
    print gamma, gamma2
    print omega
    print omega2
  ass = [max(enumerate(l), key=operator.itemgetter(1))[0] 
             for l in message_field]
  ass = [l[0] for l in message_field] 
  return ass, message_field, gamma, omega, messages

def bp_fixed_params(gamma, omega, adj, tmax):
  nodes = len(adj)
  types = len(gamma)
  messages = np.zeros([nodes,nodes,types])
  for u,v in it.product(range(nodes), range(nodes)):
    a = .5 + (.01 if u % 2 == 0 else -0.01)
    messages[u][v][0] = a
    messages[u][v][1] = 1-a
  for i in range(tmax):
    messages = update_messages(adj, messages, omega, gamma)
  message_field = get_message_field(adj, messages, omega, gamma)
  ass = [l[0] for l in message_field] 
  return ass

# attempt to follow the progress of the bp and see convergence
def graphbp(gamma, omega, adj, tmax):
  nodes = len(adj)
  types = len(gamma)
  messages = np.zeros([nodes,nodes,types])
  for u,v in it.product(range(nodes), range(nodes)):
    a = random.random()
    messages[u][v][0] = a
    messages[u][v][1] = 1-a
  xcoords = []
  ycoords = []
  message_field = get_message_field(adj, messages, omega, gamma)
  for i in range(tmax):
    print i
    newx = sum(message_field[:(nodes/2)])[0]
    newy = sum(message_field[(nodes/2):])[1]
    print newx
    print newy
    xcoords.append(newx)
    ycoords.append(newy)
    messages = update_messages(adj, messages, omega, gamma)
    joint_dist = get_joint_dist(adj, messages, omega)
    message_field = get_message_field(adj, messages, omega, gamma)
  return xcoords, ycoords

def linear_expansion(gamma, omega, adj, perturb):
  nodes = len(adj)
  types = len(gamma)
  messages = np.zeros([nodes,nodes,types])
  for u,v in it.product(range(nodes), range(nodes)):
    a = .49# + perturb[u] * .001
    messages[u][v][0] = a
    messages[u][v][1] = 1-a
  message_field = get_message_field(adj, messages, omega, gamma)

  for i in range(0):
    messages = update_messages(adj, messages, omega, gamma)
    joint_dist = get_joint_dist(adj, messages, omega)
    message_field = get_message_field(adj, messages, omega, gamma)
    omega, gamma = update_parameters(adj, messages, message_field, joint_dist)
  return message_field, gamma, omega, messages

def back_bp(gamma, omega, adj, tmax):
  nodes = len(adj)
  types = len(gamma)
  messages = np.zeros([nodes,nodes,types])
  difference = None
  for u,v in it.product(range(nodes), range(nodes)):
    
    #a = .9 if u < 15 or u in [17,18,20,22,31] else .1
    a = .5
    messages[u][v][0] = a
    messages[u][v][1] = 1-a
  for i in range(tmax):
    #print 'mes', messages
    try:
      new_messages = update_messages(adj, messages, omega, gamma)
      #print 'new_mes', new_messages
      difference = (new_messages-messages)/10.
      #print 'diff', difference
      messages = messages - difference
      for u,v in it.product(range(nodes), range(nodes)):
        if messages[u][v][0] < 0:
          messages[u][v][0] = 0
          messages[u][v][1] = 1
        elif messages[u][v][0] > 1:
          messages[u][v][0] = 1
          messages[u][v][1] = 0

      #print 'mes', messages
      #joint_dist = get_joint_dist(adj, messages, omega)
      message_field = get_message_field(adj, messages, omega, gamma)
      #omega, gamma = update_parameters(adj, messages, message_field, joint_dist)
      print 'mes_field\n', message_field
    except AssertionError:
      print "assertion error"
      #print messages
      exit(1)
  #print difference
  ass = [max(enumerate(l), key=operator.itemgetter(1))[0] 
             for l in message_field]
  ass = [l[0] for l in message_field] 
  return ass, message_field, messages

def bp_find_zeros(gamma, omega, adj, tmax):
  nodes = len(adj)
  types = len(gamma)
  messages = np.zeros([nodes,nodes,types])
  difference = None
  min_mess = None
  min_sumdiff = nodes * nodes
  for i in range(tmax):
    print i
    for u,v in it.product(range(nodes), range(nodes)):
      a = random.random()
      messages[u][v][0] = a
      messages[u][v][1] = 1-a
    new_messages = update_messages(adj, messages, omega, gamma)
    difference = (new_messages-messages)
    sumdiff = sum(sum(sum(difference)))
    if sumdiff < min_sumdiff:
      min_mess = messages
      min_sumdiff = sumdiff
  message_field = get_message_field(adj, min_mess, omega, gamma)
  ass = [max(enumerate(l), key=operator.itemgetter(1))[0] 
             for l in message_field]
  ass = [l[0] for l in message_field] 
  return ass, message_field, messages
