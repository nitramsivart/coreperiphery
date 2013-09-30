from networkx import *
import scipy.sparse.linalg as lin

# Generate a power law graph
seq = None
while seq == None or not is_valid_degree_sequence(seq):
  seq=utils.create_degree_sequence(10000,utils.powerlaw_sequence,exponent=2.9)

# Convert graph to adjacency matrix
G = Graph(configuration_model(seq))

A = to_scipy_sparse_matrix(G, dtype='d')

# select the eigenvector corresponding to the highest eigenvalue
eig = lin.eigsh(A, k=1, which='LA')[1][:,0]
print eig
# there are negative and positive values!?!?!?
