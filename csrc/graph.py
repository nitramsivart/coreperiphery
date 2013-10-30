import matplotlib.pyplot as plt

x = []
y = []
for l in open('data'):
  x.append(2*(20-int(l.split()[0])))
  y.append(float(l.split()[1])/200000.)

plt.ylabel('fraction correct')
plt.xlabel('c_cc - c_pp')
plt.axis([0,40,.5,1])
plt.scatter(x, y)
plt.show()
