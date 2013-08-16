import beliefprop as bp
import networkx as nx
import scipy.sparse as sp
import scipy.sparse.linalg as lin
import numpy as np
import matplotlib.pyplot as plt
import turtle

n_c = 40
n_p = 40
p_cc = .6
p_cp = .3
p_pp = .25

G = nx.Graph(bp.make_graph(n_c, n_p, p_cc, p_cp, p_pp))

gamma = [n_c/float(n_c+n_p),n_p/float(n_c+n_p)]
omega = [[p_cc, p_cp],[p_cp, p_pp]]

A = np.array(nx.adjacency_matrix(G))
xcoords, ycoords = bp.graphbp(gamma,omega,A,5)

for x, y in zip(xcoords, ycoords):
  turtle.goto(x*5, y*5)
turtle.exitonclick()
print xcoords
print ycoords
