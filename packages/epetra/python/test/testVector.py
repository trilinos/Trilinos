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

# Imports
import setpath
import unittest
from   Numeric    import *
from   PyTrilinos import Epetra

##########################################################################

class EpetraVector1DTestCase(unittest.TestCase):
    "TestCase class for 1D Vector objects"

    def setUp(self):
        self.length     = 9
        self.scale      = 1.0 / (self.length-1)
        self.comm       = Epetra.SerialComm()
        self.map        = Epetra.Map(self.length,0,self.comm)
        self.numPyArray = arange(self.length) * self.scale
        
    def testBlockMap(self):
        "Test Epetra.Vector 1D BlockMap constructor"
        epetraVector = Epetra.Vector(self.map)
        self.assertEqual(len(epetraVector), self.length)
        for i in range(self.length):
            self.assertEqual(epetraVector[i], 0.0)
        epetraVector[:] = self.numPyArray
        for i in range(self.length):
            self.assertEqual(epetraVector[i], self.numPyArray[i])

    def testBlockMap_Numeric(self):
        "Test Epetra.Vector 1D BlockMap, Numeric array constructor"
        epetraVector = Epetra.Vector(self.map, self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,) = epetraVector.shape
        for i in range(ni):
            self.assertEqual(epetraVector[i], self.numPyArray[i])

    def testNumeric(self):
        "Test Epetra.Vector 1D Numeric array constructor"
        epetraVector = Epetra.Vector(self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,) = epetraVector.shape
        for i in range(ni):
            self.assertEqual(epetraVector[i], self.numPyArray[i])

    def testPyObject(self):
        "Test Epetra.Vector 1D PyObject constructor"
        epetraVector = Epetra.Vector([0,1,2,3])
        self.assertEqual(epetraVector.MyLength(), 4)
        (ni,) = epetraVector.shape
        for i in range(ni):
            self.assertEqual(epetraVector[i], float(i))

#     def testNorm1(self):
#         "Test Epetra.Vector Norm1 method"
#         epetraVector  = Epetra.Vector([-1,2,-3,4])
#         (status,norm) = epetraVector.Norm1()
#         self.assertEqual(status,0)
#         self.assertEqual(norm, array([10.0]))

##########################################################################

class EpetraVector2DTestCase(unittest.TestCase):
    "TestCase class for 2D Vector objects"

    def setUp(self):
        self.length     = 3 * 3
        self.scale      = 1.0 / (self.length-1)
        self.comm       = Epetra.SerialComm()
        self.map        = Epetra.Map(self.length,0,self.comm)
        self.numPyArray = arange(self.length) * self.scale
        self.numPyArray.shape = (3,self.length/3)
        
    def testBlockMap(self):
        "Test Epetra.Vector 2D BlockMap constructor"
        epetraVector = Epetra.Vector(self.map)
        self.assertEqual(len(epetraVector), self.length)
        for i in range(self.length):
            self.assertEqual(epetraVector[i], 0.0)
        epetraVector[:] = self.numPyArray.flat
        for i in range(self.length):
            self.assertEqual(epetraVector[i], i*self.scale)

    def testBlockMap_Numeric(self):
        "Test Epetra.Vector 2D BlockMap, Numeric array constructor"
        epetraVector = Epetra.Vector(self.map, self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                self.assertEqual(epetraVector[i,j], self.numPyArray[i,j])

    def testNumeric(self):
        "Test Epetra.Vector 2D Numeric array constructor"
        epetraVector = Epetra.Vector(self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                self.assertEqual(epetraVector[i,j], self.numPyArray[i,j])

    def testPyObject(self):
        "Test Epetra.Vector 2D PyObject constructor"
        epetraVector = Epetra.Vector([[0,1],[2,3]])
        self.assertEqual(epetraVector.MyLength(), 4)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                k = nj*i + j
                self.assertEqual(epetraVector[i,j], float(k))

    def testNonContiguous(self):
        "Test Epetra.Vector 2D noncontiguous Numeric array constructor"
        nonContigArray = swapaxes(self.numPyArray,0,1)
        self.assertEqual(nonContigArray.iscontiguous(), False)
        epetraVector = Epetra.Vector(nonContigArray)
        self.assertEqual(epetraVector.iscontiguous(), True)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                self.assertEqual(epetraVector[i,j], self.numPyArray[j,i])

##########################################################################

class EpetraVector3DTestCase(unittest.TestCase):
    "TestCase class for 3D Vector objects"

    def setUp(self):
        (ni,nj,nk)      = (2,3,4)
        self.length     = ni * nj * nk
        self.scale      = 1.0 / (self.length-1)
        self.comm       = Epetra.SerialComm()
        self.map        = Epetra.Map(self.length,0,self.comm)
        self.numPyArray = arange(self.length) * self.scale
        self.numPyArray.shape = (ni,nj,nk)

    def testBlockMap(self):
        "Test Epetra.Vector 3D BlockMap constructor"
        epetraVector = Epetra.Vector(self.map)
        self.assertEqual(len(epetraVector), self.length)
        for i in range(self.length):
            self.assertEqual(epetraVector[i], 0.0)
        epetraVector[:] = self.numPyArray.flat
        for i in range(self.length):
            self.assertEqual(epetraVector[i], i*self.scale)

    def testBlockMap_Numeric(self):
        "Test Epetra.Vector 3D BlockMap, Numeric array constructor"
        epetraVector = Epetra.Vector(self.map, self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    self.assertEqual(epetraVector[i,j,k], self.numPyArray[i,j,k])

    def testNumeric(self):
        "Test Epetra.Vector 3D Numeric array constructor"
        epetraVector = Epetra.Vector(self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    self.assertEqual(epetraVector[i,j,k], self.numPyArray[i,j,k])

    def testPyObject(self):
        "Test Epetra.Vector 3D PyObject constructor"
        epetraVector = Epetra.Vector([[[0,1],[2,3]],[[4,5],[6,7]]])
        self.assertEqual(epetraVector.MyLength(), 8)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    l = nk*nj*i + nk*j + k
                    self.assertEqual(epetraVector[i,j,k], float(l))

    def testNonContiguous01(self):
        "Test Epetra.Vector 3D noncontig(0,1) Numeric array constructor"
        nonContigArray = swapaxes(self.numPyArray,0,1)
        self.assertEqual(nonContigArray.iscontiguous(), False)
        epetraVector = Epetra.Vector(nonContigArray)
        self.assertEqual(epetraVector.iscontiguous(), True)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    self.assertEqual(epetraVector[i,j,k], self.numPyArray[j,i,k])

    def testNonContiguous02(self):
        "Test Epetra.Vector 3D noncontig(0,2) Numeric array constructor"
        nonContigArray = swapaxes(self.numPyArray,0,2)
        self.assertEqual(nonContigArray.iscontiguous(), False)
        epetraVector = Epetra.Vector(nonContigArray)
        self.assertEqual(epetraVector.iscontiguous(), True)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    self.assertEqual(epetraVector[i,j,k], self.numPyArray[k,j,i])

    def testNonContiguous12(self):
        "Test Epetra.Vector 3D noncontig(1,2) Numeric array constructor"
        nonContigArray = swapaxes(self.numPyArray,1,2)
        self.assertEqual(nonContigArray.iscontiguous(), False)
        epetraVector = Epetra.Vector(nonContigArray)
        self.assertEqual(epetraVector.iscontiguous(), True)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    self.assertEqual(epetraVector[i,j,k], self.numPyArray[i,k,j])

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraVector1DTestCase))
    suite.addTest(unittest.makeSuite(EpetraVector2DTestCase))
    suite.addTest(unittest.makeSuite(EpetraVector3DTestCase))

    # Run the test suite
    print "\n*********************\nTesting Epetra.Vector\n*********************\n"
    unittest.TextTestRunner(verbosity=2).run(suite)
