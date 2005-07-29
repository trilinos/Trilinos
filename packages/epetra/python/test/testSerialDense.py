#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#              PyTrilinos.Epetra: Python Interface to Epetra
#                   Copyright (2005) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

# Imports.  Users importing an installed version of PyTrilinos should use the
# "from PyTrilinos import ..." syntax.  Here, the setpath module adds the build
# directory, including "PyTrilinos", to the front of the search path.  We thus
# use "import ..." for Trilinos modules.  This prevents us from accidentally
# picking up a system-installed version and ensures that we are testing the
# build module.
import sys

try:
    import setpath
    import Epetra
except ImportError:
    from PyTrilinos import Epetra
    print >>sys.stderr, "Using system-installed Epetra"

import unittest
from   Numeric    import *

##########################################################################

class EpetraSerialDenseTestCase(unittest.TestCase):
    "TestCase class for Epetra SerialDense Vectors, Matrices and Solvers"

    def setUp(self):
        self.size = 4
        self.rows = 2
        self.cols = 4

    def testDefaultVectorConstructor(self):
        "Test Epetra.SerialDenseVector default constructor"
        sdv = Epetra.SerialDenseVector()
        self.assertEqual(sdv.CV(), Epetra.Copy)
        self.assertEqual(sdv.Length(), 0)

    def testVectorSizeResize(self):
        "Test Epetra.SerialDenseVector Size and Resize methods"
        sdv = Epetra.SerialDenseVector()
        sdv.Size(3*self.size)
        self.assertEqual(sdv.Length(), 3*self.size)
        sdv.Resize(self.size)
        self.assertEqual(sdv.Length(), self.size)

    def testVectorSizedConstructor(self):
        "Test Epetra.SerialDenseVector sized constructor"
        sdv = Epetra.SerialDenseVector(self.size)
        self.assertEqual(sdv.CV(), Epetra.Copy)
        self.assertEqual(sdv.Length(), self.size)

    def testVectorCopyConstructor(self):
        "Test Epetra.SerialDenseVector copy constructor"
        sdv1 = Epetra.SerialDenseVector(self.size)
        sdv2 = Epetra.SerialDenseVector(sdv1)
        self.assertEqual(sdv1.Length(), sdv2.Length())

    def testVectorIndexErrors(self):
        "Test Epetra.SerialDenseVector index errors "
        sdv = Epetra.SerialDenseVector(self.size)
        self.assertRaises(TypeError, sdv.__getitem__, 0,1)
        self.assertRaises(TypeError, sdv.__setitem__, 0,1,3.14)

    def testMatrixDefaultConstructor(self):
        "Test Epetra.SerialDenseMatrix default constructor"
        sdm = Epetra.SerialDenseMatrix()
        self.assertEqual(sdm.CV(), Epetra.Copy)
        self.assertEqual(sdm.M(), 0)
        self.assertEqual(sdm.N(), 0)

    def testMatrixShapeReshape(self):
        "Test Epetra.SerialDenseMatrix Shape and Reshape methods"
        sdm = Epetra.SerialDenseMatrix()
        sdm.Shape(self.rows,self.cols)
        self.assertEqual(sdm.M(), self.rows)
        self.assertEqual(sdm.N(), self.cols)
        sdm.Reshape(2*self.rows,2*self.cols)
        self.assertEqual(sdm.M(), 2*self.rows)
        self.assertEqual(sdm.N(), 2*self.cols)

    def testMatrixSizedConstructor(self):
        "Test Epetra.SerialDenseMatrix sized constructor"
        sdm = Epetra.SerialDenseMatrix(self.rows,self.cols)
        self.assertEqual(sdm.M(), self.rows)
        self.assertEqual(sdm.N(), self.cols)

    def testMatrixCopyConstructor(self):
        "Test Epetra.SerialDenseMatrix copy constructor"
        sdm1 = Epetra.SerialDenseMatrix(self.rows,self.cols)
        sdm2 = Epetra.SerialDenseMatrix(sdm1)
        self.assertEqual(sdm1.M(), sdm2.M())
        self.assertEqual(sdm1.N(), sdm2.N())

    def testMatrixNormOneNormInf(self):
        "Test Epetra.SerialDenseMatrix NormOne and NormInf methods"
        scalar = 2.0
        size   = self.size
        sdm    = Epetra.SerialDenseMatrix(size,size)
        for i in range(size):
            for j in range(size):
                sdm[i,j] = scalar
        self.assertEqual(sdm.NormOne(), scalar*size)
        self.assertEqual(sdm.NormInf(), scalar*size)

    def testMatrixIndexErrors(self):
        "Test Epetra.SerialDenseMatrix index errors "
        sdm = Epetra.SerialDenseMatrix(self.size,self.size)
        #self.assertRaises(TypeError, sdm.__getitem__, 0,1,2)
        self.assertRaises(TypeError, sdm.__setitem__, 0,1,2,3.14)

    def testSolver(self):
        "Test Epetra.SerialDenseSolver"
        size = self.size
        sdm  = Epetra.SerialDenseMatrix(size,size)
        for i in range(size):
            if (i>0): sdm[i,i-1] = 1
            sdm[i,i] = -2
            if (i<size-1): sdm[i,i+1] = 1
        inv  = Epetra.SerialDenseMatrix(sdm)
        sys  = Epetra.SerialDenseSolver()
        sys.SetMatrix(inv)
        sys.Invert()
        idty = Epetra.SerialDenseMatrix(size,size)
        idty.Multiply("N","N",1,sdm,inv,0)
        if "assertAlmostEqual" in dir(unittest.TestCase):
            for i in range(size):
                for j in range(size):
                    if i==j: self.assertAlmostEqual(idty[i,j],1.0,10)
                    else:    self.assertAlmostEqual(idty[i,j],0.0,10)

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraSerialDenseTestCase))

    # Run the test suite
    print >>sys.stderr, \
          "\n**************************\nTesting Epetra.SerialDense\n**************************\n"
    unittest.TextTestRunner(verbosity=2).run(suite)
