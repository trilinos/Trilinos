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

## class EpetraCompObjectTestCase(unittest.TestCase):
##     "TestCase class for ComObjects"

##     def setUp(self):
##         self.compObject = Epetra.CompObject()

##     def test(self):

class EpetraVectorTestCase(unittest.TestCase):
    "TestCase class for Vector objects"

    def setUp(self):
        self.length     = 3 * 3
        self.scale      = 1.0 / (self.length-1)
        self.comm       = Epetra.SerialComm()
        self.map        = Epetra.Map(self.length,0,self.comm)
        self.numPyArray = arange(self.length) * self.scale
        self.numPyArray.shape = (3,self.length/3)
        
    def testBlockMap(self):
        "Test the Vector BlockMap constructor"
        epetraVector = Epetra.Vector(self.map)
        self.assertEqual(len(epetraVector), self.length)
        for i in range(self.length):
            self.assertEqual(epetraVector[i], 0.0)
        epetraVector[:] = self.numPyArray.flat
        for i in range(self.length):
            self.assertEqual(epetraVector[i], i*self.scale)

    def testBlockMap_Numeric(self):
        "Test the Vector BlockMap, Numeric array constructor"
        epetraVector = Epetra.Vector(self.map, self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                self.assertEqual(epetraVector[i,j], self.numPyArray[i,j])

    def testNumeric(self):
        "Test the Vector Numeric array constructor"
        epetraVector = Epetra.Vector(self.numPyArray)
        self.assertEqual(epetraVector.MyLength(), self.length)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                self.assertEqual(epetraVector[i,j], self.numPyArray[i,j])

    def testPyObject(self):
        "Test the Vector PyObject constructor"
        epetraVector = Epetra.Vector([[0,1],[2,3]])
        self.assertEqual(epetraVector.MyLength(), 4)
        (ni,nj) = epetraVector.shape
        for i in range(ni):
            for j in range(nj):
                k = nj*i + j
                self.assertEqual(epetraVector[i,j], float(k))

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(EpetraObjectTestCase    ))
    suite.addTest(unittest.makeSuite(EpetraSerialCommTestCase))
    suite.addTest(unittest.makeSuite(EpetraBlockMapTestCase  ))
    suite.addTest(unittest.makeSuite(EpetraMapTestCase       ))
    suite.addTest(unittest.makeSuite(EpetraVectorTestCase    ))
    unittest.TextTestRunner(verbosity=2).run(suite)
