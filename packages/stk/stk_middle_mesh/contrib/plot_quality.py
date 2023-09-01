#!/usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.tri
import matplotlib.pyplot as plt
import sys

# plots the output of printVertEdges

if len(sys.argv) == 3:
  fname = sys.argv[1];
  fname2 = sys.argv[2];
else:
  print("Usage: {} fname".format(sys.argv[0]))
  exit(1)


qualities = np.loadtxt(fname)
qualities = np.sort(qualities)

qualities2 = np.loadtxt(fname2)
qualities2 = np.sort(qualities2)


fig, ax = plt.subplots()
ax.plot(qualities, 'ro')
ax.plot(qualities2, 'bo')
ax.set_xlabel("vert")
ax.set_ylabel("Mesh Quality")
ax.legend(["before", "after"])

fig.savefig(fname + ".png", dpi=300)
