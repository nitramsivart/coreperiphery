from pydot import *

graph = Dot()
for l in open("G1"):
  (n1, n2, x) = l.split(' ')
  n1 = Node(n1)
  n2 = Node(n2)
  graph.add_node(n1)
  graph.add_node(n2)
  graph.add_edge(Edge(n1, n2))
print "done reading"
graph.write('picture.png', prog='sfdp')
