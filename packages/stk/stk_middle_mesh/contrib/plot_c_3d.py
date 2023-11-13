#!/usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.tri
from mpl_toolkits import mplot3d
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

print("vert radii")
for i in range(verts.shape[0]):
  print("vert ", i, ", r = ", np.sqrt(verts[i, 0]**2 + verts[i, 1]**2 + verts[i, 2]**2))


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(verts[:, 0], verts[:, 1], verts[:, 2], 'ro')

#for i in range(verts.shape[0]):
#  ax.text(verts[i, 0], verts[i, 1], verts[i,2], '%s' % (str(i)))


for i in range(edges.shape[0]):
  print("edge between verts ", edges[i, 0], " and ", edges[i, 1])
  print(verts[edges[i, :], 1])
  ax.plot3D(verts[edges[i, :], 0], verts[edges[i, :], 1], verts[edges[i, :], 2], 'b-')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

fig.savefig(fname + ".png", dpi=300)
plt.show()
