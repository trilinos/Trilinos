#!/usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.tri
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import sys

# plots the output of printVertEdges

if len(sys.argv) == 3:
  fname1 = sys.argv[1];
  fname2 = sys.argv[2];
else:
  print("Usage: {} fname".format(sys.argv[0]))
  exit(1)


verts1 = np.loadtxt(fname1 + "_verts.txt")
edges1 = np.loadtxt(fname1 + "_edges.txt", int)

verts2 = np.loadtxt(fname2 + "_verts.txt")
edges2 = np.loadtxt(fname2 + "_edges.txt", int)


#print("verts = \n", verts)
#print("edges = \n", edges)

#print("vert radii")
#for i in range(verts.shape[0]):
#  print("vert ", i, ", r = ", np.sqrt(verts[i, 0]**2 + verts[i, 1]**2 + verts[i, 2]**2))


fig = plt.figure()
ax = plt.axes(projection='3d')


ax.plot3D(verts2[:, 0], verts2[:, 1], verts2[:, 2], 'go')

for i in range(verts2.shape[0]):
  ax.text(verts2[i, 0], verts2[i, 1], verts2[i,2], '%s' % (str(i)))

for i in range(edges2.shape[0]):
  print("edge between verts ", edges2[i, 0], " and ", edges2[i, 1])
  #print(verts[edges[i, :], 1])
  ax.plot3D(verts2[edges2[i, :], 0], verts2[edges2[i, :], 1], verts2[edges2[i, :], 2], 'y-')

ax.plot3D(verts1[:, 0], verts1[:, 1], verts1[:, 2], 'ro')

for i in range(edges1.shape[0]):
  print("edge between verts ", edges1[i, 0], " and ", edges1[i, 1])
  #print(verts[edges[i, :], 1])
  ax.plot3D(verts1[edges1[i, :], 0], verts1[edges1[i, :], 1], verts1[edges1[i, :], 2], 'b-')





ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

fig.savefig(fname1 + "_both.png", dpi=300)
plt.show()
