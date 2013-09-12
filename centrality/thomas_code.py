from networkx import *
import scipy.sparse.linalg as lin
import scipy.sparse as sp
import matplotlib.pylab as plt


# computes and graphs several values as hubsize increases past the
# transition point
def mega():
  min_hub = 1
  max_hub = 150
  samples = max_hub - min_hub + 1

  # highest eigenvalue
  eval1 = [0] * samples
  # second highest eigenvalue
  eval2 = [0] * samples
  # highest eigenvalue, hashimoto
  evalh = [0] * samples

  # order parameter: the leading eigenvector value for the hub vertex. This is used
  # in the Newman/Rao paper.
  # I've used other things for this order parameter, such as hub value squared,
  # or hub value squared / sum of other values squared
  op = [0] * samples
  # order parameter for the hashimoto matrix
  oph = [0] * samples
  # the average leading eigenvector value of the hub neighbors
  nb = [0] * samples
  # same, but for hashimoto matrix
  nbh = [0] * samples

  c = 10
  n = 20000
  f = open('mega-data-%d-%d.txt' % (c, n), 'w+')
  print c, n
  G = fast_gnp_random_graph(n, c/float(n))
  for k in range(1, min_hub):
    G.add_edge(0, k)
  A = to_scipy_sparse_matrix(G, dtype='d',format='csr')
  for i, k in enumerate(range(min_hub, max_hub+1)):
    print k
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
  plt.title("eigenvalues (adj & hashimoto)")
  plt.xlabel("hub size")
  plt.plot(x, eval1)
  plt.plot(x, eval2)
  plt.plot(x, evalh)
  plt.figure()
  plt.title("eigenvalues (hashimoto)")
  plt.xlabel("hub size")
  plt.plot(x, evalh)
  plt.figure()
  plt.title("order parameter")
  plt.xlabel("hub size")
  plt.plot(x, op)
  plt.figure()
  plt.title("order parameter (hashimoto)")
  plt.xlabel("hub size")
  plt.plot(x, oph)
  plt.figure()
  plt.title("average hub neighbor value")
  plt.xlabel("hub size")
  plt.plot(x, nb)
  plt.figure()
  plt.title("average hub neighbor value (hashimoto)")
  plt.xlabel("hub size")
  plt.plot(x, nbh)
  plt.show()

  for i in range(samples):
    f.write("%d,%r,%r,%r,%r,%r,%r,%r\n" % 
            (i,eval1[i],eval2[i],evalh[i],op[i],oph[i],nb[i],nbh[i]))
  f.close()

mega()
