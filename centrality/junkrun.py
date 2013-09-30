from dummy import *
from paper_graphs import *
import matplotlib.pylab as plt

interact()
'''
results = interact()
G1, G2, G3, A1, A2, A3 = results

safe = set([])
for n in G3.neighbors(0):
    safe.add(n)
    for j in G3.neighbors(n):
        safe.add(j)
for n in range(len(G3.nodes())):
    if n not in safe:
        G1.remove_node(n)
        G2.remove_node(n)
        G3.remove_node(n)
        
pos = graphviz_layout(G3)
results1 = dummy(A1, G1, pos)
results2 = dummy(A2, G2, pos)
results3 = dummy(A3, G3, pos)
plt.show()
'''
