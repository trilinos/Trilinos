""" This code requires SymPy, an open source package for symbolic computation
It also depends on NumPy.  See scipy.org and sympy.org

Compute the l2-norm of x*y*z over the cube[-1/2:1/2, -1/2:1/2, -1/2:1/2]

"""

from sympy import symbols, Symbol, factorial, Rational, zeros, div, eye, \
        integrate, diff, pprint
import numpy
x, y, z, t, tg = symbols('x y z t tg')

def main():

   # same as unitTest_h1Norm.py, but rotate the mesh around the z axis by 30 degrees first, then integrate
   deg = 30.0
   theta = numpy.pi*deg/180.0
   rm = numpy.array([[numpy.cos(theta),-numpy.sin(theta),0.0],[numpy.sin(theta),numpy.cos(theta),0.0],[0.0,0.0,1.0]])
   print rm
   print numpy.dot(rm,[x,y,z])[0]

   t = ((numpy.dot(rm,[x,y,z])[0]) + 2*(numpy.dot(rm,[x,y,z])[1]) + 3*(numpy.dot(rm,[x,y,z])[2]) )**2
   #tg = ((numpy.dot(rm,[1,0,0])[0])**2 + (numpy.dot(rm,[0,2,0])[1])**2 + (numpy.dot(rm,[0,0,3])[2])**2)
   tg = (1**2+2**2+3**2)
   print "t= " , t, " tg= " , tg

   intf = integrate(t + tg , (x, -0.5, 0.5))

   intf = integrate(intf, (y, -0.5, 0.5))
   intf = integrate(intf, (z, -0.5, 0.5))
   #  the answers should be as shown:
   print "expected: 15.1666666666667   3.89444048184931"
   print "  actual: ", intf, " ", intf**(0.5)

if __name__ == "__main__":
   main()
