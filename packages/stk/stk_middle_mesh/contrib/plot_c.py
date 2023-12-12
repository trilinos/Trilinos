#!/usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.tri
import matplotlib.pyplot as plt
import sys

# plots the output of printVertEdges

if len(sys.argv) == 2:
  fname = sys.argv[1];
else:
  print("Usage: {} fname".format(sys.argv[0]))
  exit(1)


verts = np.loadtxt(fname + "_verts.txt")
edges = np.loadtxt(fname + "_edges.txt", int)

print("verts = \n", verts)
print("edges = \n", edges)

#tri = matplotlib.tri.Triangulation(verts[:, 0], verts[:, 1], tris)

fig, ax = plt.subplots()
ax.plot(verts[:, 0], verts[:, 1], 'ro')

for i in range(verts.shape[0]):
  ax.annotate(i, (verts[i, 0], verts[i, 1]))

for i in range(edges.shape[0]):
  print("edge between verts ", edges[i, 0], " and ", edges[i, 1])
  print("v1 = ", verts[edges[i, 0], :], ", v2 = ", verts[edges[i, 1], :])
  ax.plot(verts[edges[i, :], 0], verts[edges[i, :], 1], 'b-')

fig.savefig(fname + ".png", dpi=300)

plt.show()
