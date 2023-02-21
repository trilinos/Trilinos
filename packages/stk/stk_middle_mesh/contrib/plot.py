#!/usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.tri
import matplotlib.pyplot as plt
import sys

# plots the output of printTIN

if len(sys.argv) == 2:
  fname = sys.argv[1];
else:
  print("Usage: {} fname".format(sys.argv[0]))
  exit(1)


verts = np.loadtxt(fname + "_verts.txt")
tris = np.loadtxt(fname + "_triangles.txt", int)

print("verts = \n", verts)
print("triangles = \n", tris)

tri = matplotlib.tri.Triangulation(verts[:, 0], verts[:, 1], tris)

fig, ax = plt.subplots()
ax.triplot(tri)

fig.savefig(fname + ".png", dpi=300)
