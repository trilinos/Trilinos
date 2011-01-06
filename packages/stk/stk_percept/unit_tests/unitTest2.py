""" This code requires SymPy, an open source package for symbolic computation
It also depends on NumPy.  See scipy.org and sympy.org

Compute the l2-norm of x*y*z over the cube[-1/2:1/2, -1/2:1/2, -1/2:1/2], after rotating
it around the z-axis by 30 degrees
"""

from sympy import symbols, Symbol, factorial, Rational, zeros, div, eye, \
        integrate, diff, pprint

import numpy
x, y, z = symbols('x y z')

def main():
   # same as unitTest1.py, but rotate the mesh around the z axis by 30 degrees first, then integrate
   deg = 30.0
   theta = numpy.pi*deg/180.0
   rm = numpy.array([[numpy.cos(theta),-numpy.sin(theta),0.0],[numpy.sin(theta),numpy.cos(theta),0.0],[0.0,0.0,1.0]])
   print rm
   print numpy.dot(rm,[x,y,z])[0]
   intf = integrate( ((numpy.dot(rm,[x,y,z])[0])*(numpy.dot(rm,[x,y,z])[1])*(numpy.dot(rm,[x,y,z])[2]))**2, (x, -0.5, 0.5))
   intf = integrate(intf, (y, -0.5, 0.5))
   intf = integrate(intf, (z, -0.5, 0.5))
   #  the answers should be as shown:
   print intf, " ", intf**(0.5), "0.000318287037037037   0.0178406008037016"

if __name__ == "__main__":
   main()
