#!/usr/bin/env python
#-------------------------------------------------------------------------------

import sys, os
import string
import math

t0 = 0.0
t1 = 1.0
Nt = 1001
dt = (t1-t0)/(Nt-1.0)

a = 0.0
f = 1.0
L = 1.0

gamma0 = 0
gamma1 = 1

phi = math.atan((f/L)*(gamma0-a)/gamma1) - f/L*t0
b   = gamma1/((f/L)*math.cos((f/L)*t0+phi))

print 'phi = %g' % (phi)
print 'b   = %g' % (b)

fout = open('sincos.dat', 'w')
for i in range(Nt):
  t  = dt*i
  x0 = a + b*math.sin((f/L)*t+phi)
  x1 = b*(f/L)*math.cos((f/L)*t+phi)
  line = '%g %g %g\n' % (t, x0, x1)
  fout.write(line)

fout.close()
