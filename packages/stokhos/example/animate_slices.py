#!/usr/bin/env python

# Script to generate an animation of the Cijk triple product tensor
# moving through k.  It uses an executable in stokhos/example to generate
# the sparsity patterns, and so needs to be run in that directory.
#
# The script can generate an mp4 movie of the animation by setting
# the save_figs variable below to True.  Creating the moving requires
# ffmpeg.
# 
# Note, if a case is run that has more than 999 terms, then the %03d in
# the filenames below needs to be changed (to e.g., %04d).

from matplotlib.pylab import *
from os import system, getcwd
from scipy.io import mmread
from scipy.sparse import csr_matrix
from subprocess import call

p = 3   # polynomial order
d = 3   # number of random variables (stocahstic dimension)
save_figs = True # save figures into a movie

# Generate slices
system("rm -f A_*.mm frame_*.png")
exe = getcwd() + "/Stokhos_sparsity_slice_example.exe"
args = "--drop=1e-12 --full --order=" + str(p) + " --dimension=" + str(d)
n = call([exe] + args.split())

# Turn on interactive mode
ion()

# Make animation of slices
B = csr_matrix((n,n))
for i in range(0,n):
    A = mmread("A_"+str(i)+".mm").tocsr()
    B = B + A
    clf()
    spy(A,markersize=2)
    title("order = "+str(p)+", dimension = "+str(d)+", k = "+str(i))
    draw()
    if save_figs:
        savefig("frame_%03d.png" % i, dpi=100)
clf()
spy(B,markersize=2)
title("order = "+str(p)+", dimension = "+str(d)+", final matrix")
draw()
if save_figs:
    savefig("frame_%03d.png" % n, dpi=100)

# Make movie using ffmpeg
if save_figs:
    movie_name = "Cijk_"+str(p)+"_"+str(d)+".mp4";
    system("ffmpeg -y -i frame_%03d.png "+movie_name)
    print "\nCreated movie in file",movie_name

# clean up
system("rm -f A_*.mm frame_*.png")
