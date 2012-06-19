""" This code requires SymPy, an open source package for symbolic computation
It also depends on NumPy.  See scipy.org and sympy.org

Compute the l2-norm of x*y*z over the cube[-1/2:1/2, -1/2:1/2, -1/2:1/2]

"""

from sympy import symbols, Symbol, factorial, Rational, zeros, div, eye, \
        integrate, diff, pprint
import numpy
x, y, z = symbols('x y z')

def main():
   intf = integrate( (x+2*y+3*z)**2 + (1**2+2**2+3**2), (x, -0.5, 0.5))

   # same as unitTest_h1Norm.py, but rotate the mesh around the z axis by 30 degrees first, then integrate
   deg = 30.0
   theta = numpy.pi*deg/180.0
   rm = numpy.array([[numpy.cos(theta),-numpy.sin(theta),0.0],[numpy.sin(theta),numpy.cos(theta),0.0],[0.0,0.0,1.0]])
   print rm
   print numpy.dot(rm,[x,y,z])[0]

   intf = integrate(
      ((numpy.dot(rm,[x,y,z])[0]) + 2*(numpy.dot(rm,[x,y,z])[1]) + 3*(numpy.dot(rm,[x,y,z])[2]) )**2 +
      ((numpy.dot(rm,[1,0,0])[0])**2 + (numpy.dot(rm,[0,2,0])[1])**2 + (numpy.dot(rm,[0,0,3])[2])**2),
      (x, -0.5, 0.5))
   intf = integrate(intf, (y, -0.5, 0.5))
   intf = integrate(intf, (z, -0.5, 0.5))
   #  the answers should be as shown:
   print intf, " ", intf**(0.5), "13.9166666666667   3.73050488093323"

if __name__ == "__main__":
   main()
