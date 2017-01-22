#! /usr/bin/env python

import numpy
import pylab

# This should set variables nx, ny, nt, sf, ff
execfile("params.py")

# Build the mesh
x = numpy.arange(nx) / (nx-1.0)
y = numpy.arange(ny) / (ny-1.0)
xm, ym = numpy.meshgrid(x, y)

for n in range(0, nt+1, ff):
    try:
        fname = "u%d.bin" % n
        u = numpy.fromfile(fname)
        u.shape = (nx,ny)
        fname = "v%d.bin" % n
        v = numpy.fromfile(fname)
        v.shape = (nx,ny)
        print "Time step", n
        print "    Plotting u"
        pylab.contourf(xm, ym, u)
        pylab.show()
        print "    Plotting v"
        pylab.contourf(xm,ym, v)
        pylab.show()
    except IOError:
        print "File '%s' not found" % fname
        break
