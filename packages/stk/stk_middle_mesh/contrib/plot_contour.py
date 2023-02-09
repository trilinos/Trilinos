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

xi1_vals = np.loadtxt(fname + "_xi1_vals.txt")
xi2_vals = np.loadtxt(fname + "_xi2_vals.txt")
f_vals  = np.loadtxt(fname + "_func_vals.txt")

print("min val = ", np.amin(f_vals))

X, Y = np.meshgrid(xi1_vals, xi2_vals);


fig, ax = plt.subplots()
c = ax.contourf(X, Y, f_vals, 20)

ax.set_xlabel("xi1")
ax.set_ylabel("xi2")
fig.colorbar(c, ax=ax)

fig.savefig(fname + ".png", dpi=300)
