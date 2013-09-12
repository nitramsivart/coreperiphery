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
  assert adj == True or adj == False, adj
  #return (omega**adj)*np.exp(-omega)
  return omega if adj else 1.-omega/float(n)


# adj is the adjacency matrix of the graph
# messages is an N x N x [#types] array, giving marginal type probs
# omega is the block affinity matrix
# gamma is the prior on probability of being in a certain group
def update_messages(adj, messages, omega, gamma, tmax):
  #new_messages = np.copy(messages)
  n = float(len(messages))
  n_range = xrange(len(messages)) # number of nodes
  k_range = xrange(len(messages[0][0])) # number of types
  l = list(it.product(n_range, n_range))
  #random.shuffle(l)

  h = [0,0]
  message_field = get_message_field(adj, messages, omega, gamma)
  for u, r, s in it.product(n_range, k_range, k_range):
    h[s] += omega[r][s] * message_field[u][r] / float(n)

  for t in range(tmax):
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
        for w in adj.neighbors(u):
          if w == v:
            continue

          sum_total = 0.
          for s in k_range:
            assert 0 <= messages[w][u][s] <= 1, messages[w][u][s]
            #sum_total += messages[w][u][s] * edge_prob(adj[w][u], omega[r][s],len(n_range))
            sum_total += messages[w][u][s] * omega[r][s]
          #assert 0 <= sum_total <= 1
          prod_total += np.log(sum_total)
        messages[u][v][r] = gamma[r]*np.exp(prod_total - h[r])
        norm_total += messages[u][v][r]

      norm_total = 0.
      old_message_field = np.copy(message_field[u])
      for r in k_range: #normalizing messages
        messages[u][v][r] = messages[u][v][r] / norm_total
        if messages[u][v][r] > 1:
          messages[u][v][r] = 1
        #messages[u][v][r] = new_messages[u][v][r]

        #update message_field
        prod_total = 1.
        for w in adj.neighbors(u):
          sum_total = 0.
          for s in k_range:
            sum_total += messages[w][u][r] * omega[r][s]
          prod_total *= sum_total
        message_field[u][r] = gamma[r] * prod_total * np.exp(-h[r])
        norm_total += message_field[u][r]

      #normalize message_field
      for r in k_range:
        message_field[u][r] /= norm_total
      for r, s in it.product(k_range, k_range):
        h[r] += (message_field[u][s] - old_message_field[s]) * omega[r][s] / n

      #assert message_field[u][0] != old_message_field[0]
      assert 0 <= messages[u][v][0] <= 1, pdb.set_trace()
    return messages

def get_message_field(adj, messages, omega, gamma):
  n = len(messages)
  k = len(messages[0][0])
  message_field = np.zeros([n, k])
  for u, r in it.product(range(n), range(k)):
    prod_total = 0.
    for w in range(n):
      if w == u:
        continue

      sum_total = (messages[w][u][0] * edge_prob(u in adj[w], omega[r][0], n) + 
                   messages[w][u][1] * edge_prob(u in adj[w], omega[r][1], n))

      prod_total += np.log(sum_total)
    message_field[u][r] = gamma[r]*np.exp(prod_total)

  # normalize
  for u in range(n):
    sum_total = sum(message_field[u])
    for r in range(k):
      message_field[u][r] /= sum_total
      assert 0 <= message_field[u][r] <= 1
  return message_field


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

def bp_fixed_params(gamma, omega, adj, tmax):
  nodes = len(adj)
  types = len(gamma)
  messages = np.zeros([nodes,nodes,types])
  for u,v in it.product(range(nodes), range(nodes)):
    a = .5 + (.01 if u % 2 == 0 else -0.01)
    messages[u][v][0] = a
    messages[u][v][1] = 1-a
  messages = update_messages(adj, messages, omega, gamma, tmax)
  message_field = get_message_field(adj, messages, omega, gamma)
  ass = [l[0] for l in message_field] 
  return ass
