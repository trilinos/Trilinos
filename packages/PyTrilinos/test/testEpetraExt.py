#! /usr/bin/env python

# System imports
from   Numeric  import *
import unittest

# Local imports
import setpath
from   PyTrilinos import Epetra
from   PyTrilinos import EpetraExt


class EpetraExtTestCase(unittest.TestCase):
    "TestCase class for EpetraExt objects"

    def setUp(self):
        self.size = 11
        self.comm = Epetra.SerialComm()
        self.map  = Epetra.Map(self.size,0,self.comm)
        self.crsg = Epetra.CrsGraph(Epetra.Copy, self.map, 3)
        self.crsg.InsertGlobalIndices(0,2,array([0,1]))
        for i in range(1,self.size-1):
            self.crsg.InsertGlobalIndices(i,3,array([i-1,i,i+1]))
        self.crsg.InsertGlobalIndices(self.size-1,2,array([self.size-2,self.size-1]))
        self.crsg.FillComplete()
        colors = array([3,2,1])
        self.colors = resize(colors,(self.size,))
        self.numEl  = [0,0,0,0]
        for i in self.colors:
            self.numEl[i] += 1

    def testMapColoring(self):
        "Test the EpetraExt CrsGraph-to-MapColoring transform"
        mapColoring  = EpetraExt.CrsGraph_MapColoring(False)
        colorMap     = mapColoring(self.crsg)
        numColors    = colorMap.NumColors()
        defaultColor = colorMap.DefaultColor()
        self.assertEqual(numColors   , max(self.colors))
        self.assertEqual(defaultColor, 0               )
        for c in range(numColors):
            self.assertEqual(colorMap.NumElementsWithColor(c+1), self.numEl[c+1])

    def testColorMapIndex(self):
        "Test the EpetraExt MapColoring-to-ColorMapIndex transform"
        mapColoring   = EpetraExt.CrsGraph_MapColoring(False)
        colorMap      = mapColoring(self.crsg)
        colorMapIndex = EpetraExt.CrsGraph_MapColoringIndex(colorMap)
        columns       = colorMapIndex(self.crsg)
        # No assert test here; just executing has to be enough

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the tets cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraExtTestCase))

    # Run the test suite
    unittest.TextTestRunner(verbosity=2).run(suite)
