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

class EpetraCrsGraphTestCase(unittest.TestCase):
    "TestCase class for Epetra CrsGraphs"

    def setUp(self):
        self.comm      = Epetra.PyComm()
        self.myPID     = self.comm.MyPID()
        self.numProc   = self.comm.NumProc()
        self.mySize    = 11
        self.size      = self.mySize * self.numProc
        self.indexBase = 0
        self.rowMap    = Epetra.Map(self.size, self.indexBase, self.comm)
        mge            = list(self.rowMap.MyGlobalElements())
        if self.indexBase not in mge: mge.insert(0,mge[0] -1)
        if self.size-1    not in mge: mge.append(  mge[-1]+1)
        self.colMap    = Epetra.Map(-1, mge, self.indexBase, self.comm)
        self.nipr      = ones(self.mySize)
        self.nipr[1:-1] += 1
        self.nipr[2:-2] += 1

    def tearDown(self):
        self.comm.Barrier()

    def fillGraphGlobal(self,graph,complete=True):
        n = self.size
        for lrid in range(graph.NumMyRows()):
            grid = graph.GRID(lrid)
            if   grid == 0  : indices = [0,1]
            elif grid == n-1: indices = [n-2,n-1]
            else            : indices = [grid-1,grid,grid+1]
            graph.InsertGlobalIndices(grid,indices)
        if complete: graph.FillComplete()

    def fillGraphLocal(self,graph,complete=True):
        n = self.size
        o = 0
        if self.myPID > 0: o = 1
        for lrid in range(graph.NumMyRows()):
            grid = graph.GRID(lrid)
            if   grid == 0  : indices = [0,1]
            elif grid == n-1: indices = [lrid+o-1,lrid+o]
            else            : indices = [lrid+o-1,lrid+o,lrid+o+1]
            graph.InsertMyIndices(lrid,indices)
        if complete: graph.FillComplete()

    def testConstructor01(self):
        "Test Epetra.CrsGraph constructor, no colMap w/fixed # of indices/row"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.NumMyRows(), self.mySize   )
        self.assertEqual(crsg.IndexBase(), self.indexBase)

    def testConstructor02(self):
        "Test Epetra.CrsGraph constructor, no colMap w/fixed # of indices/row, static profile"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3, True)
        self.assertEqual(crsg.NumMyRows(), self.mySize   )
        self.assertEqual(crsg.IndexBase(), self.indexBase)

    def testConstructor03(self):
        "Test Epetra.CrsGraph constructor, no colMap w/specified # of indices/row"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.nipr)
        self.assertEqual(crsg.NumMyRows(), self.mySize   )
        self.assertEqual(crsg.IndexBase(), self.indexBase)

    def testConstructor04(self):
        "Test Epetra.CrsGraph constructor, no colMap w/specified # of indices/row, static profile"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.nipr, True)
        self.assertEqual(crsg.NumMyRows(), self.mySize   )
        self.assertEqual(crsg.IndexBase(), self.indexBase)

    def testConstructor05(self):
        "Test Epetra.CrsGraph constructor, w/colMap w/fixed # of indices/row"
        try:
            crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.colMap, 3)
            self.assertEqual(crsg.NumMyRows(), self.mySize   )
            self.assertEqual(crsg.IndexBase(), self.indexBase)
        except TypeError:
            # A TypeError is raised for the wrapper generated by older versions
            # of SWIG (1.3.28 and earlier).  This is expected, so just ignore it
            pass

    def testConstructor06(self):
        "Test Epetra.CrsGraph constructor, w/colMap w/fixed # of indices/row, static profile"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.colMap, 3, True)
        self.assertEqual(crsg.NumMyRows(), self.mySize   )
        self.assertEqual(crsg.IndexBase(), self.indexBase)

    def testConstructor07(self):
        "Test Epetra.CrsGraph constructor, w/colMap w/specified # of indices/row"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.colMap, self.nipr)
        self.assertEqual(crsg.NumMyRows(), self.mySize   )
        self.assertEqual(crsg.IndexBase(), self.indexBase)

    def testConstructor08(self):
        "Test Epetra.CrsGraph constructor, w/colMap w/specified # of indices/row, static profile"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.colMap, self.nipr, True)
        self.assertEqual(crsg.NumMyRows(), self.mySize   )
        self.assertEqual(crsg.IndexBase(), self.indexBase)

    def testConstructor09(self):
        "Test Epetra.CrsGraph copy constructor"
        crsg1 = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.nipr)
        crsg2 = Epetra.CrsGraph(crsg1)
        self.assertEqual(crsg2.NumMyRows(), self.mySize   )
        self.assertEqual(crsg2.IndexBase(), self.indexBase)

    def testConstructor10(self):
        "Test Epetra.CrsGraph constructor with too-short array"
        self.assertRaises(ValueError, Epetra.CrsGraph, Epetra.Copy, self.rowMap,
                          self.nipr[1:-1])

    def testConstructor11(self):
        "Test Epetra.CrsGraph constructor with too-long array"
        nipr = list(self.nipr)
        nipr.append(1)
        self.assertRaises(ValueError, Epetra.CrsGraph, Epetra.Copy, self.rowMap,
                          nipr)

    def testInsertGlobalIndices(self):
        "Test Epetra.CrsGraph InsertGlobalIndices method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.fillGraphGlobal(crsg)  # This calls crsg.InsertGlobalIndices()
        n = self.size
        for lrid in range(crsg.NumMyRows()):
            grid = crsg.GRID(lrid)
            if grid in (0,n-1): numIndices = 2
            else:               numIndices = 3
            self.assertEqual(crsg.NumMyIndices(lrid), numIndices)

    def testInsertGlobalIndicesBad(self):
        "Test Epetra.CrsGraph InsertGlobalIndices method for bad indices"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        grid = crsg.GRID(0)
        self.assertRaises(ValueError, crsg.InsertGlobalIndices, grid, [0,"e","pi"])

    def testRemoveGlobalIndices1(self):
        "Test Epetra.CrsGraph RemoveGlobalIndices method w/specified indices"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.fillGraphGlobal(crsg,False)  # Cannot FillComplete() if we want to remove
        grid = crsg.GRID(0)
        indices = crsg.ExtractGlobalRowCopy(grid)
        crsg.RemoveGlobalIndices(grid,indices[1:])
        self.assertEqual(crsg.NumMyIndices(0), 1)

    def testRemoveGlobalIndices2(self):
        "Test Epetra.CrsGraph RemoveGlobalIndices method w/o indices"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.fillGraphGlobal(crsg,False)  # Cannot FillComplete() if we want to remove
        grid = crsg.GRID(0)
        crsg.RemoveGlobalIndices(grid)
        self.assertEqual(crsg.NumMyIndices(0), 0)

    def testInsertMyIndices(self):
        "Test Epetra.CrsGraph InsertMyIndices method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillGraphLocal(crsg)  # This calls crsg.InsertMyIndices()
        n = self.size
        for lrid in range(crsg.NumMyRows()):
            grid = crsg.GRID(lrid)
            if grid in (0,n-1): numIndices = 2
            else:               numIndices = 3
            self.assertEqual(crsg.NumMyIndices(lrid), numIndices)

    def testRemoveMyIndices1(self):
        "Test Epetra.CrsGraph RemoveMyIndices method w/specified indices"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillGraphLocal(crsg,False)  # Cannot FillComplete() if we want to remove
        indices = crsg.ExtractMyRowCopy(0)
        crsg.RemoveMyIndices(0,indices[1:])
        self.assertEqual(crsg.NumMyIndices(0), 1)

    def testRemoveMyIndices2(self):
        "Test Epetra.CrsGraph RemoveMyIndices method w/o indices"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillGraphLocal(crsg,False)  # Cannot FillComplete() if we want to remove
        lrid = self.mySize / 2
        result = crsg.RemoveMyIndices(lrid)
        self.assertEqual(crsg.NumMyIndices(lrid), 0)

    def testFillComplete1(self):
        "Test Epetra.CrsGraph FillComplete method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.Filled(), False)
        result = crsg.FillComplete()
        self.assertEqual(result,        0   )
        self.assertEqual(crsg.Filled(), True)

    def testFillComplete2(self):
        "Test Epetra.CrsGraph FillComplete method w/specified maps"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.Filled(), False)
        result = crsg.FillComplete(self.rowMap, self.colMap)
        self.assertEqual(result,        0   )
        self.assertEqual(crsg.Filled(), True)

    def testOptimizeStorage(self):
        "Test Epetra.CrsGraph OptimizeStorage method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.fillGraphGlobal(crsg,False)
        grid = crsg.GRID(0)
        numIndices = crsg.NumMyIndices(0)
        crsg.InsertGlobalIndices(grid,[grid])  # Duplicate, thus non-optimal
        self.assertEqual(crsg.NumMyIndices(0), numIndices+1)
        crsg.FillComplete()
        crsg.OptimizeStorage()
        self.assertEqual(crsg.NumMyIndices(0), numIndices)

    def testExtractGlobalRowCopy(self):
        "Test Epetra.CrsGraph ExtractGlobalRowCopy method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.fillGraphGlobal(crsg)
        n = crsg.GRID(self.mySize-1)
        indices = crsg.ExtractGlobalRowCopy(n)
        ncol = 3
        if n in (0,self.size-1): ncol -= 1
        self.assertEqual(len(indices), ncol)
        self.assertEqual(indices[0]  , n-1 )
        self.assertEqual(indices[1]  , n   )

    def testExtractGlobalRowCopyBad(self):
        "Test Epetra.CrsGraph ExtractGlobalRowCopy method, bad index"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.fillGraphGlobal(crsg)
        self.assertRaises(ValueError, crsg.ExtractGlobalRowCopy, self.size)

    def testExtractMyRowCopy(self):
        "Test Epetra.CrsGraph ExtractMyRowCopy method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.fillGraphGlobal(crsg)
        n = self.mySize-1
        indices = crsg.ExtractMyRowCopy(n)
        ncol = 3
        if crsg.GRID(n) in (0,self.size-1): ncol -= 1
        self.assertEqual(len(indices), ncol)
        self.assertEqual(indices[0]  , n-1 )
        self.assertEqual(indices[1]  , n   )

    def testExtractMyRowCopyBad(self):
        "Test Epetra.CrsGraph ExtractMyRowCopy method, bad index"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.fillGraphGlobal(crsg)
        self.assertRaises(ValueError, crsg.ExtractMyRowCopy, self.mySize)

    def testFilled(self):
        "Test Epetra.CrsGraph Filled method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.Filled(), False)
        self.fillGraphGlobal(crsg)
        self.assertEqual(crsg.Filled(), True)

    def testStorageOptimized(self):
        "Test Epetra.CrsGraph StorageOptimized method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.StorageOptimized(), False)
        self.fillGraphGlobal(crsg)
        self.assertEqual(crsg.StorageOptimized(), False)
        crsg.OptimizeStorage()
        self.assertEqual(crsg.StorageOptimized(), True)

    def testIndicesAreGlobal(self):
        "Test Epetra.CrsGraph IndicesAreGlobal method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.IndicesAreGlobal(), False)
        self.fillGraphGlobal(crsg,False)
        self.assertEqual(crsg.IndicesAreGlobal(), True)

    def testIndicesAreLocal(self):
        "Test Epetra.CrsGraph IndicesAreLocal method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.assertEqual(crsg.IndicesAreLocal(), False)
        self.fillGraphLocal(crsg)
        self.assertEqual(crsg.IndicesAreLocal(), True)

    def testLowerTriangular(self):
        "Test Epetra.CrsGraph LowerTriangular method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.LowerTriangular(), True)
        self.fillGraphGlobal(crsg)
        self.assertEqual(crsg.LowerTriangular(), False)

    def testUpperTriangular(self):
        "Test Epetra.CrsGraph UpperTriangular method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.UpperTriangular(), True)
        self.fillGraphGlobal(crsg)
        self.assertEqual(crsg.UpperTriangular(), False)

    def testNoDiagonal(self):
        "Test Epetra.CrsGraph NoDiagonal method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.NoDiagonal(), True)
        self.fillGraphGlobal(crsg)
        self.assertEqual(crsg.NoDiagonal(), False)

    def testMyGlobalRow(self):
        "Test Epetra.CrsGraph MyGlobalRow method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        myRows = self.rowMap.MyGlobalElements()
        for grid in range(self.size):
            self.assertEqual(crsg.MyGlobalRow(grid), grid in myRows)

    def testHaveColMap(self):
        "Test Epetra.CrsGraph HaveColMap method"
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsg.HaveColMap(), False)
        crsg = Epetra.CrsGraph(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.assertEqual(crsg.HaveColMap(), True)

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraCrsGraphTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n***********************\nTesting Epetra.CrsGraph\n***********************\n"
    verbosity = 2 * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))[0]
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
