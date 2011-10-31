#!/usr/bin/env python

# This script requires scipy (including numpy) and matplotlib

# You must run Stokhos_Linear2D_Diffusion_PCE_Example.exe and 
# Stokhos_Linear2D_Diffusion_PCE_Interlaced_Example.exe first!

from scipy.io import mmread
from numpy import transpose, reshape
from numpy.linalg import norm
from matplotlib.pyplot import *

# number of grid points in each spatial dimension
n = 32

# number of stochastic terms
P = 10

# read in block, interlaced solutions
x = mmread('stochastic_solution.mm')
x2 = mmread('stochastic_solution_interlaced.mm')

# permute interlaced solution to block
xx2 = reshape(transpose(reshape(x2,[P,n*n],'F')),[P*n*n,1],'F')

# compute error
e = norm(x-xx2)

print 'Error between block and interlace solutions is',e

# read in interlaced operator
A = mmread('stochastic_operator_interlaced.mm').tocsr()

# display interlaced operator
figure(1)
clf()
spy(A,markersize=1)

figure(2)
clf()
spy(A[400:500,400:500],markersize=1)
show()
