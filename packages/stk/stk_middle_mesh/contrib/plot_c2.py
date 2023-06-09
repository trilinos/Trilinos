#!/usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.tri
import matplotlib.pyplot as plt
import sys

# plots the output of printVertEdges for two files on the same pplot

if len(sys.argv) == 3:
  fname = sys.argv[1];
  fname2 = sys.argv[2];
else:
  print("Usage: {} fname1 fname2".format(sys.argv[0]))
  exit(1)


fig, ax = plt.subplots()

# plot first file
verts = np.loadtxt(fname + "_verts.txt")
edges = np.loadtxt(fname + "_edges.txt", int)

print("verts = \n", verts)
print("edges = \n", edges)


ax.plot(verts[:, 0], verts[:, 1], 'ro')

for i in range(edges.shape[0]):
  print("edge between verts ", edges[i, 0], " and ", edges[i, 1])
  print(verts[edges[i, :], 1])
  ax.plot(verts[edges[i, :], 0], verts[edges[i, :], 1], 'r-')

#for i in range(verts.shape[0]):
#  ax.annotate("%3.2f, %3.2f" % (verts[i, 0], verts[i, 1]), (verts[i, 0], verts[i, 1]))


# plot second file
verts = np.loadtxt(fname2 + "_verts.txt")
edges = np.loadtxt(fname2 + "_edges.txt", int)

print("verts = \n", verts)
print("edges = \n", edges)


ax.plot(verts[:, 0], verts[:, 1], 'bo')

for i in range(edges.shape[0]):
  print("edge between verts ", edges[i, 0], " and ", edges[i, 1])
  print(verts[edges[i, :], 1])
  ax.plot(verts[edges[i, :], 0], verts[edges[i, :], 1], 'b-')

fig.savefig(fname + "_both.png", dpi=300)
plt.show()
