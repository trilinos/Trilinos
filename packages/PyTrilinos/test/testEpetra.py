#! /usr/bin/env python

# Imports
import setpath
import unittest
from   Numeric    import *
from   PyTrilinos import Epetra

class EpetraObjectTestCase(unittest.TestCase):
    "TestCase for Epetra_Objects"

    def setUp(self):
        self.object = Epetra.Object()

    def testLabel(self):
        "Test the Object Label method"
        self.assertEqual(self.object.Label(), 'Epetra::Object')

    def testSetLabel(self):
        "Test the Object SetLabel method"
        label = 'New Label'
        self.object.SetLabel(label)
        self.assertEqual(self.object.Label(), label)

    def testGetTracebackMode(self):
        "Test the Object GetTracebackMode method"
        self.assertEqual(self.object.GetTracebackMode(), 1)

    def testSetTracebackMode(self):
        "Test the Object SetTracebackMode method"
        tracebackMode = 2
        self.object.SetTracebackMode(tracebackMode)
        self.assertEqual(self.object.GetTracebackMode(), tracebackMode)

class EpetraSerialCommTestCase(unittest.TestCase):
    "TestCase class for SerialComm communicator objects"

    def setUp(self):
        self.comm = Epetra.SerialComm()

    def testMyPID(self):
        "Test the SerialComm MyPID method"
        self.assertEqual(self.comm.MyPID()  , 0)

    def testNumProc(self):
        "Test the SerialComm NumProc method"
        self.assertEqual(self.comm.NumProc(), 1)

class EpetraBlockMapTestCase(unittest.TestCase):
    "TestCase class for BlockMap objects"

    def setUp(self):
        self.comm = Epetra.SerialComm()
        self.map  = Epetra.BlockMap(2,2,0,self.comm)

    def testNumGlobalElements(self):
        "Test the BlockMap NumGlobalElements method"
        self.assertEqual(self.map.NumGlobalElements(), 2)

    def testNumMyElements(self):
        "Test the BlockMap NumMyElements method"
        self.assertEqual(self.map.NumMyElements(), 2)

    def testElementSize(self):
        "Test the BlockMap ElementSize method"
        self.assertEqual(self.map.ElementSize(), 2)

    def testIndexBase(self):
        "Test the BlockMap IndexBase method"
        self.assertEqual(self.map.IndexBase(), 0)

    def testNumGlobalPoints(self):
        "Test the BlockMap NumGlobalPoints method"
        self.assertEqual(self.map.NumGlobalPoints(), 4)

    def testNumMyPoints(self):
        "Test the BlockMap NumMyPoints method"
        self.assertEqual(self.map.NumMyPoints(), 4)

    def testMinMyElementSize(self):
        "Test the BlockMap MinMyElementSize method"
        self.assertEqual(self.map.MinMyElementSize(), self.map.ElementSize())

    def testMaxMyElementSize(self):
        "Test the BlockMap MaxMyElementSize method"
        self.assertEqual(self.map.MaxMyElementSize(), self.map.ElementSize())

    def testMinElementSize(self):
        "Test the BlockMap MinElementSize method"
        self.assertEqual(self.map.MinElementSize(), self.map.ElementSize())

    def testMaxElementSize(self):
        "Test the BlockMap MaxElementSize method"
        self.assertEqual(self.map.MaxElementSize(), self.map.ElementSize())

    def testConstantElementSize(self):
        "Test the BlockMap ConstantElementSize method"
        self.assertEqual(self.map.ConstantElementSize(), True)

    def testDistributedGlobal(self):
        "Test the BlockMap DistributedGlobal method"
        self.assertEqual(self.map.DistributedGlobal(), False)

    def testMinAllGID(self):
        "Test the BlockMap MinAllGID method"
        self.assertEqual(self.map.MinAllGID(), 0)

    def testMaxAllGID(self):
        "Test the BlockMap MaxAllGID method"
        self.assertEqual(self.map.MaxAllGID(), 1)

    def testMinMyGID(self):
        "Test the BlockMap MinMyGID method"
        self.assertEqual(self.map.MinMyGID(), 0)

    def testMaxMyGID(self):
        "Test the BlockMap MaxMyGID method"
        self.assertEqual(self.map.MaxMyGID(), 1)

    def testMinLID(self):
        "Test the BlockMap MinLID method"
        self.assertEqual(self.map.MinLID(), 0)

    def testMaxLID(self):
        "Test the BlockMap MaxLID method"
        self.assertEqual(self.map.MaxLID(), 1)

    def testIDs(self):
        "Test the BlockMap local and global IDs"
        for i in range(self.map.NumMyElements()):
            self.assertEqual(self.map.LID(i)  , self.map.GID(i))
            self.assertEqual(self.map.MyGID(i), True           )
            self.assertEqual(self.map.MyLID(i), True           )

class EpetraMapTestCase(unittest.TestCase):
    "TestCase class for Map objects"

    def setUp(self):
        self.comm = Epetra.SerialComm()
        self.map  = Epetra.Map(4,0,self.comm)

    def testNumGlobalElements(self):
        "Test the Map NumGlobalElements method"
        self.assertEqual(self.map.NumGlobalElements(), 4)

    def testNumMyElements(self):
        "Test the Map NumMyElements method"
        self.assertEqual(self.map.NumMyElements(), 4)

    def testElementSize(self):
        "Test the Map ElementSize method"
        self.assertEqual(self.map.ElementSize(), 1)

    def testIndexBase(self):
        "Test the Map IndexBase method"
        self.assertEqual(self.map.IndexBase(), 0)

    def testNumGlobalPoints(self):
        "Test the Map NumGlobalPoints method"
        self.assertEqual(self.map.NumGlobalPoints(), 4)

    def testNumMyPoints(self):
        "Test the Map NumMyPoints method"
        self.assertEqual(self.map.NumMyPoints(), 4)

    def testMinMyElementSize(self):
        "Test the Map MinMyElementSize method"
        self.assertEqual(self.map.MinMyElementSize(), self.map.ElementSize())

    def testMaxMyElementSize(self):
        "Test the Map MaxMyElementSize method"
        self.assertEqual(self.map.MaxMyElementSize(), self.map.ElementSize())

    def testMinElementSize(self):
        "Test the Map MinElementSize method"
        self.assertEqual(self.map.MinElementSize(), self.map.ElementSize())

    def testMaxElementSize(self):
        "Test the Map MaxElementSize method"
        self.assertEqual(self.map.MaxElementSize(), self.map.ElementSize())

    def testConstantElementSize(self):
        "Test the Map ConstantElementSize method"
        self.assertEqual(self.map.ConstantElementSize(), True)

    def testSameAs(self):
        "Test the Map SameAs method"
        self.assertEqual(self.map.SameAs(self.map), True)

    def testDistributedGlobal(self):
        "Test the Map DistributedGlobal method"
        self.assertEqual(self.map.DistributedGlobal(), False)

    def testMinAllGID(self):
        "Test the Map MinAllGID method"
        self.assertEqual(self.map.MinAllGID(), 0)

    def testMaxAllGID(self):
        "Test the Map MaxAllGID method"
        self.assertEqual(self.map.MaxAllGID(), 3)

    def testMinMyGID(self):
        "Test the Map MinMyGID method"
        self.assertEqual(self.map.MinMyGID(), 0)

    def testMaxMyGID(self):
        "Test the Map MaxMyGID method"
        self.assertEqual(self.map.MaxMyGID(), 3)

    def testMinLID(self):
        "Test the Map MinLID method"
        self.assertEqual(self.map.MinLID(), 0)

    def testMaxLID(self):
        "Test the Map MaxLID method"
        self.assertEqual(self.map.MaxLID(), 3)

    def testIDs(self):
        "Test the Map global and local IDs"
        for i in range(self.map.NumMyElements()):
            self.assertEqual(self.map.LID(i)  , self.map.GID(i))
            self.assertEqual(self.map.MyGID(i), True           )
            self.assertEqual(self.map.MyLID(i), True           )


class EpetraVector1DTestCase(unittest.TestCase):
    "TestCase class for 1D Vector objects"

    def setUp(self):
        self.length     = 9
        self.scale      = 1.0 / (self.length-1)
        self.comm       = Epetra.SerialComm()
        self.map        = Epetra.Map(self.length,0,self.comm)
        self.numPyArray = arange(self.length) * self.scale
        
    def testBlockMap(self):
        "Test the 1D Vector BlockMap constructor"
        epetraVector = Epetra.Vector(self.map)
        self.assertEqual(len(epetraVector), self.length)
        for i in range(self.length):
            self.assertEqual(epetraVector[i], 0.0)
        epetraVector[:] = self.numPyArray
        for i in range(self.length):
            self.assertEqual(epetraVector[i], self.numPyArray[i])

    def testBlockMap_Numeric(self):
        "Test the 1D Vector BlockMap, Numeric array constructor"
        epetraVector = Epetra.Vector(self.map, self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,) = epetraVector.shape
        for i in range(ni):
            self.assertEqual(epetraVector[i], self.numPyArray[i])

    def testNumeric(self):
        "Test the 1D Vector Numeric array constructor"
        epetraVector = Epetra.Vector(self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,) = epetraVector.shape
        for i in range(ni):
            self.assertEqual(epetraVector[i], self.numPyArray[i])

    def testPyObject(self):
        "Test the 1D Vector PyObject constructor"
        epetraVector = Epetra.Vector([0,1,2,3])
        self.assertEqual(epetraVector.MyLength(), 4)
        (ni,) = epetraVector.shape
        for i in range(ni):
            self.assertEqual(epetraVector[i], float(i))

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
        "Test the 2D Vector BlockMap constructor"
        epetraVector = Epetra.Vector(self.map)
        self.assertEqual(len(epetraVector), self.length)
        for i in range(self.length):
            self.assertEqual(epetraVector[i], 0.0)
        epetraVector[:] = self.numPyArray.flat
        for i in range(self.length):
            self.assertEqual(epetraVector[i], i*self.scale)

    def testBlockMap_Numeric(self):
        "Test the 2D Vector BlockMap, Numeric array constructor"
        epetraVector = Epetra.Vector(self.map, self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                self.assertEqual(epetraVector[i,j], self.numPyArray[i,j])

    def testNumeric(self):
        "Test the 2D Vector Numeric array constructor"
        epetraVector = Epetra.Vector(self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                self.assertEqual(epetraVector[i,j], self.numPyArray[i,j])

    def testPyObject(self):
        "Test the 2D Vector PyObject constructor"
        epetraVector = Epetra.Vector([[0,1],[2,3]])
        self.assertEqual(epetraVector.MyLength(), 4)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                k = nj*i + j
                self.assertEqual(epetraVector[i,j], float(k))

    def testNonContiguous(self):
        "Test the 2D Vector noncontiguous Numeric array constructor"
        nonContigArray = swapaxes(self.numPyArray,0,1)
        self.assertEqual(nonContigArray.iscontiguous(), False)
        epetraVector = Epetra.Vector(nonContigArray)
        self.assertEqual(epetraVector.iscontiguous(), True)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                self.assertEqual(epetraVector[i,j], self.numPyArray[j,i])

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
        "Test the 3D Vector BlockMap constructor"
        epetraVector = Epetra.Vector(self.map)
        self.assertEqual(len(epetraVector), self.length)
        for i in range(self.length):
            self.assertEqual(epetraVector[i], 0.0)
        epetraVector[:] = self.numPyArray.flat
        for i in range(self.length):
            self.assertEqual(epetraVector[i], i*self.scale)

    def testBlockMap_Numeric(self):
        "Test the 3D Vector BlockMap, Numeric array constructor"
        epetraVector = Epetra.Vector(self.map, self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    self.assertEqual(epetraVector[i,j,k], self.numPyArray[i,j,k])

    def testNumeric(self):
        "Test the 3D Vector Numeric array constructor"
        epetraVector = Epetra.Vector(self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    self.assertEqual(epetraVector[i,j,k], self.numPyArray[i,j,k])

    def testPyObject(self):
        "Test the 3D Vector PyObject constructor"
        epetraVector = Epetra.Vector([[[0,1],[2,3]],[[4,5],[6,7]]])
        self.assertEqual(epetraVector.MyLength(), 8)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    l = nk*nj*i + nk*j + k
                    self.assertEqual(epetraVector[i,j,k], float(l))

    def testNonContiguous01(self):
        "Test the 3D Vector noncontig(0,1) Numeric array constructor"
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
        "Test the 3D Vector noncontig(0,2) Numeric array constructor"
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
        "Test the 3D Vector noncontig(1,2) Numeric array constructor"
        nonContigArray = swapaxes(self.numPyArray,1,2)
        self.assertEqual(nonContigArray.iscontiguous(), False)
        epetraVector = Epetra.Vector(nonContigArray)
        self.assertEqual(epetraVector.iscontiguous(), True)
        (ni,nj,nk) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    self.assertEqual(epetraVector[i,j,k], self.numPyArray[i,k,j])

class EpetraSerialDenseTestCase(unittest.TestCase):
    "TestCase class for Epetra SerialDense Vectors, Matrices and Solvers"

    def setUp(self):
        self.size = 4
        self.rows = 2
        self.cols = 4

    def testDefaultVectorConstructor(self):
        "Test the SerialDenseVector default constructor"
        sdv = Epetra.SerialDenseVector()
        self.assertEqual(sdv.CV(), Epetra.Copy)
        self.assertEqual(sdv.Length(), 0)

    def testVectorSizeResize(self):
        "Test the SerialDenseVector Size and Resize methods"
        sdv = Epetra.SerialDenseVector()
        sdv.Size(3*self.size)
        self.assertEqual(sdv.Length(), 3*self.size)
        sdv.Resize(self.size)
        self.assertEqual(sdv.Length(), self.size)

    def testVectorSizedConstructor(self):
        "Test the SerialDenseVector sized constructor"
        sdv = Epetra.SerialDenseVector(self.size)
        self.assertEqual(sdv.CV(), Epetra.Copy)
        self.assertEqual(sdv.Length(), self.size)

    def testVectorCopyConstructor(self):
        "Test the SerialDenseVector copy constructor"
        sdv1 = Epetra.SerialDenseVector(self.size)
        sdv2 = Epetra.SerialDenseVector(sdv1)
        self.assertEqual(sdv1.Length(), sdv2.Length())

    def testVectorIndexErrors(self):
        "Test the SerialDenseVector index errors "
        sdv = Epetra.SerialDenseVector(self.size)
        self.assertRaises(TypeError, sdv.__getitem__, 0,1)
        self.assertRaises(TypeError, sdv.__setitem__, 0,1,3.14)

    def testMatrixDefaultConstructor(self):
        "Test the SerialDenseMatrix default constructor"
        sdm = Epetra.SerialDenseMatrix()
        self.assertEqual(sdm.CV(), Epetra.Copy)
        self.assertEqual(sdm.M(), 0)
        self.assertEqual(sdm.N(), 0)

    def testMatrixShapeReshape(self):
        "Test the SerialDenseMatrix Shape and Reshape methods"
        sdm = Epetra.SerialDenseMatrix()
        sdm.Shape(self.rows,self.cols)
        self.assertEqual(sdm.M(), self.rows)
        self.assertEqual(sdm.N(), self.cols)
        sdm.Reshape(2*self.rows,2*self.cols)
        self.assertEqual(sdm.M(), 2*self.rows)
        self.assertEqual(sdm.N(), 2*self.cols)

    def testMatrixSizedConstructor(self):
        "Test the SerialDenseMatrix sized constructor"
        sdm = Epetra.SerialDenseMatrix(self.rows,self.cols)
        self.assertEqual(sdm.M(), self.rows)
        self.assertEqual(sdm.N(), self.cols)

    def testMatrixCopyConstructor(self):
        "Test the SerialDenseMatrix copy constructor"
        sdm1 = Epetra.SerialDenseMatrix(self.rows,self.cols)
        sdm2 = Epetra.SerialDenseMatrix(sdm1)
        self.assertEqual(sdm1.M(), sdm2.M())
        self.assertEqual(sdm1.N(), sdm2.N())

    def testMatrixNormOneNormInf(self):
        "Test the SerialDenseMatrix NormOne and NormInf methods"
        scalar = 2.0
        size   = self.size
        sdm    = Epetra.SerialDenseMatrix(size,size)
        for i in range(size):
            for j in range(size):
                sdm[i,j] = scalar
        self.assertEqual(sdm.NormOne(), scalar*size)
        self.assertEqual(sdm.NormInf(), scalar*size)

    def testMatrixIndexErrors(self):
        "Test the SerialDenseMatrix index errors "
        sdm = Epetra.SerialDenseMatrix(self.size,self.size)
        self.assertRaises(TypeError, sdm.__getitem__, 0,1,2)
        self.assertRaises(TypeError, sdm.__setitem__, 0,1,2,3.14)

    def testSolver(self):
        "Test the SerialDenseSolver"
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
        for i in range(size):
            for j in range(size):
                if i==j: self.assertAlmostEqual(idty[i,j],1.0,10)
                else:    self.assertAlmostEqual(idty[i,j],0.0,10)


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
    suite.addTest(unittest.makeSuite(EpetraSerialDenseTestCase))

    # Run the test suite
    unittest.TextTestRunner(verbosity=2).run(suite)
