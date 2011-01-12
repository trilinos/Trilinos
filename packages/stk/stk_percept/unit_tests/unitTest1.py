""" This code requires SymPy, an open source package for symbolic computation
It also depends on NumPy.  See scipy.org and sympy.org

Compute the l2-norm of x*y*z over the cube[-1/2:1/2, -1/2:1/2, -1/2:1/2]

"""

from sympy import symbols, Symbol, factorial, Rational, zeros, div, eye, \
        integrate, diff, pprint
import numpy
x, y, z = symbols('x y z')

def main():
   intf = integrate( (x*y*z)**2, (x, -0.5, 0.5))
   intf = integrate(intf, (y, -0.5, 0.5))
   intf = integrate(intf, (z, -0.5, 0.5))
   #  the answers should be as shown:
   print intf, " ", intf**(0.5), "0.000578703703703704   0.0240562612162344"

if __name__ == "__main__":
   main()
