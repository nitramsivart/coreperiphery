

sums = [0] * 6
for l in open('all-eigs-10-1000000.txt'):
  s = l.split()
  for index, val in enumerate(s):
    sums[index] += float(val)**2

print sums
