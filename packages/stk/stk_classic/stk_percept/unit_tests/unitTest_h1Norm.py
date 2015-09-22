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
   intf = integrate(intf, (y, -0.5, 0.5))
   intf = integrate(intf, (z, -0.5, 0.5))
   #  the answers should be as shown:
   print "expected: 15.1666666666667   3.89444048184931"
   print "  actual: ", intf, " ", intf**(0.5)

if __name__ == "__main__":
   main()
