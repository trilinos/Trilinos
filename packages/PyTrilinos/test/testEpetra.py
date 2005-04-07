#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                  PyTrilinos: Rapid Prototyping Package
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

class EpetraObjectTestCase(unittest.TestCase):
    "TestCase for Epetra_Objects"

    def setUp(self):
        self.object = Epetra.Object()

    def testLabel(self):
        "Test Epetra.Object Label method"
        self.assertEqual(self.object.Label(), 'Epetra::Object')

    def testSetLabel(self):
        "Test Epetra.Object SetLabel method"
        label = 'New Label'
        self.object.SetLabel(label)
        self.assertEqual(self.object.Label(), label)

    def testGetTracebackMode(self):
        "Test Epetra.Object GetTracebackMode method"
        self.assertEqual(self.object.GetTracebackMode(), 1)

    def testSetTracebackMode(self):
        "Test Epetra.Object SetTracebackMode method"
        tracebackMode = 2
        self.object.SetTracebackMode(tracebackMode)
        self.assertEqual(self.object.GetTracebackMode(), tracebackMode)

##########################################################################

class EpetraSerialCommTestCase(unittest.TestCase):
    "TestCase class for SerialComm communicator objects"

    def setUp(self):
        self.comm = Epetra.SerialComm()

    def testMyPID(self):
        "Test Epetra.SerialComm MyPID method"
        self.assertEqual(self.comm.MyPID()  , 0)

    def testNumProc(self):
        "Test Epetra.SerialComm NumProc method"
        self.assertEqual(self.comm.NumProc(), 1)

##########################################################################

class EpetraBlockMapTestCase(unittest.TestCase):
    "TestCase class for BlockMap objects"

    def setUp(self):
        self.comm = Epetra.SerialComm()
        self.map  = Epetra.BlockMap(2,2,0,self.comm)

    def testNumGlobalElements(self):
        "Test Epetra.BlockMap NumGlobalElements method"
        self.assertEqual(self.map.NumGlobalElements(), 2)

    def testNumMyElements(self):
        "Test Epetra.BlockMap NumMyElements method"
        self.assertEqual(self.map.NumMyElements(), 2)

    def testElementSize(self):
        "Test Epetra.BlockMap ElementSize method"
        self.assertEqual(self.map.ElementSize(), 2)

    def testIndexBase(self):
        "Test Epetra.BlockMap IndexBase method"
        self.assertEqual(self.map.IndexBase(), 0)

    def testNumGlobalPoints(self):
        "Test Epetra.BlockMap NumGlobalPoints method"
        self.assertEqual(self.map.NumGlobalPoints(), 4)

    def testNumMyPoints(self):
        "Test Epetra.BlockMap NumMyPoints method"
        self.assertEqual(self.map.NumMyPoints(), 4)

    def testMinMyElementSize(self):
        "Test Epetra.BlockMap MinMyElementSize method"
        self.assertEqual(self.map.MinMyElementSize(), self.map.ElementSize())

    def testMaxMyElementSize(self):
        "Test Epetra.BlockMap MaxMyElementSize method"
        self.assertEqual(self.map.MaxMyElementSize(), self.map.ElementSize())

    def testMinElementSize(self):
        "Test Epetra.BlockMap MinElementSize method"
        self.assertEqual(self.map.MinElementSize(), self.map.ElementSize())

    def testMaxElementSize(self):
        "Test Epetra.BlockMap MaxElementSize method"
        self.assertEqual(self.map.MaxElementSize(), self.map.ElementSize())

    def testConstantElementSize(self):
        "Test Epetra.BlockMap ConstantElementSize method"
        self.assertEqual(self.map.ConstantElementSize(), True)

    def testDistributedGlobal(self):
        "Test Epetra.BlockMap DistributedGlobal method"
        self.assertEqual(self.map.DistributedGlobal(), False)

    def testMinAllGID(self):
        "Test Epetra.BlockMap MinAllGID method"
        self.assertEqual(self.map.MinAllGID(), 0)

    def testMaxAllGID(self):
        "Test Epetra.BlockMap MaxAllGID method"
        self.assertEqual(self.map.MaxAllGID(), 1)

    def testMinMyGID(self):
        "Test Epetra.BlockMap MinMyGID method"
        self.assertEqual(self.map.MinMyGID(), 0)

    def testMaxMyGID(self):
        "Test Epetra.BlockMap MaxMyGID method"
        self.assertEqual(self.map.MaxMyGID(), 1)

    def testMinLID(self):
        "Test Epetra.BlockMap MinLID method"
        self.assertEqual(self.map.MinLID(), 0)

    def testMaxLID(self):
        "Test Epetra.BlockMap MaxLID method"
        self.assertEqual(self.map.MaxLID(), 1)

    def testIDs(self):
        "Test Epetra.BlockMap local and global IDs"
        for i in range(self.map.NumMyElements()):
            self.assertEqual(self.map.LID(i)  , self.map.GID(i))
            self.assertEqual(self.map.MyGID(i), True           )
            self.assertEqual(self.map.MyLID(i), True           )

##########################################################################

class EpetraMapTestCase(unittest.TestCase):
    "TestCase class for Map objects"

    def setUp(self):
        self.comm = Epetra.SerialComm()
        self.map  = Epetra.Map(4,0,self.comm)

    def testNumGlobalElements(self):
        "Test Epetra.Map NumGlobalElements method"
        self.assertEqual(self.map.NumGlobalElements(), 4)

    def testNumMyElements(self):
        "Test Epetra.Map NumMyElements method"
        self.assertEqual(self.map.NumMyElements(), 4)

    def testElementSize(self):
        "Test Epetra.Map ElementSize method"
        self.assertEqual(self.map.ElementSize(), 1)

    def testIndexBase(self):
        "Test Epetra.Map IndexBase method"
        self.assertEqual(self.map.IndexBase(), 0)

    def testNumGlobalPoints(self):
        "Test Epetra.Map NumGlobalPoints method"
        self.assertEqual(self.map.NumGlobalPoints(), 4)

    def testNumMyPoints(self):
        "Test Epetra.Map NumMyPoints method"
        self.assertEqual(self.map.NumMyPoints(), 4)

    def testMinMyElementSize(self):
        "Test Epetra.Map MinMyElementSize method"
        self.assertEqual(self.map.MinMyElementSize(), self.map.ElementSize())

    def testMaxMyElementSize(self):
        "Test Epetra.Map MaxMyElementSize method"
        self.assertEqual(self.map.MaxMyElementSize(), self.map.ElementSize())

    def testMinElementSize(self):
        "Test Epetra.Map MinElementSize method"
        self.assertEqual(self.map.MinElementSize(), self.map.ElementSize())

    def testMaxElementSize(self):
        "Test Epetra.Map MaxElementSize method"
        self.assertEqual(self.map.MaxElementSize(), self.map.ElementSize())

    def testConstantElementSize(self):
        "Test Epetra.Map ConstantElementSize method"
        self.assertEqual(self.map.ConstantElementSize(), True)

    def testSameAs(self):
        "Test Epetra.Map SameAs method"
        self.assertEqual(self.map.SameAs(self.map), True)

    def testDistributedGlobal(self):
        "Test Epetra.Map DistributedGlobal method"
        self.assertEqual(self.map.DistributedGlobal(), False)

    def testMinAllGID(self):
        "Test Epetra.Map MinAllGID method"
        self.assertEqual(self.map.MinAllGID(), 0)

    def testMaxAllGID(self):
        "Test Epetra.Map MaxAllGID method"
        self.assertEqual(self.map.MaxAllGID(), 3)

    def testMinMyGID(self):
        "Test Epetra.Map MinMyGID method"
        self.assertEqual(self.map.MinMyGID(), 0)

    def testMaxMyGID(self):
        "Test Epetra.Map MaxMyGID method"
        self.assertEqual(self.map.MaxMyGID(), 3)

    def testMinLID(self):
        "Test Epetra.Map MinLID method"
        self.assertEqual(self.map.MinLID(), 0)

    def testMaxLID(self):
        "Test Epetra.Map MaxLID method"
        self.assertEqual(self.map.MaxLID(), 3)

    def testIDs(self):
        "Test Epetra.Map global and local IDs"
        for i in range(self.map.NumMyElements()):
            self.assertEqual(self.map.LID(i)  , self.map.GID(i))
            self.assertEqual(self.map.MyGID(i), True           )
            self.assertEqual(self.map.MyLID(i), True           )

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

class EpetraCrsGraphTestCase(unittest.TestCase):
    "TestCase class for Epetra CrsGraphs"

    def setUp(self):
        self.size = 11
        self.comm = Epetra.SerialComm()
        self.map  = Epetra.Map(self.size, 0, self.comm)

    def testConstructor1(self):
        "Test Epetra.CrsGraph constructor with fixed number of indices per row"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.map, 3)

    def testInsertGlobalIndices(self):
        "Test Epetra.CrsGraph InsertGlobalIndices method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.map, 3)
        crsg.InsertGlobalIndices(0,2,array([0,1]))
        for i in range(1,self.size-1):
            crsg.InsertGlobalIndices(i,3,array([i-1,i,i+1]))
        crsg.InsertGlobalIndices(self.size-1,2,array([self.size-2,self.size-1]))

    def testInsertMyIndices(self):
        "Test Epetra.CrsGraph InsertMyIndices method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.map, 3)
        crsg.InsertMyIndices(0,2,array([0,1]))
        for i in range(1,self.size-1):
            crsg.InsertMyIndices(i,3,array([i-1,i,i+1]))
        crsg.InsertMyIndices(self.size-1,2,array([self.size-2,self.size-1]))

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
    suite.addTest(unittest.makeSuite(EpetraObjectTestCase     ))
    suite.addTest(unittest.makeSuite(EpetraSerialCommTestCase ))
    suite.addTest(unittest.makeSuite(EpetraBlockMapTestCase   ))
    suite.addTest(unittest.makeSuite(EpetraMapTestCase        ))
    suite.addTest(unittest.makeSuite(EpetraVector1DTestCase   ))
    suite.addTest(unittest.makeSuite(EpetraVector2DTestCase   ))
    suite.addTest(unittest.makeSuite(EpetraVector3DTestCase   ))
    suite.addTest(unittest.makeSuite(EpetraCrsGraphTestCase   ))
    suite.addTest(unittest.makeSuite(EpetraSerialDenseTestCase))

    # Run the test suite
    unittest.TextTestRunner(verbosity=2).run(suite)
