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

from   Numeric  import *
import unittest

##########################################################################

class EpetraImportExportTestCase(unittest.TestCase):
    "TestCase for Epetra_Import and Epetra_Export objects"

    def setUp(self):
        self.comm       = Epetra.PyComm()
        self.numProc    = self.comm.NumProc()
        self.myPID      = self.comm.MyPID()
        self.myN        = 2
        self.globalN    = self.myN * self.numProc
        self.globalSize = self.globalN*self.globalN
        self.gids       = arange(self.globalSize)
        self.gids.shape = (self.globalN,self.globalN)
        self.start      = self.myN * self.myPID
        self.end        = self.start + self.myN
        self.map1       = Epetra.Map(-1,self.gids[:,self.start:self.end],0,self.comm)
        self.map2       = Epetra.Map(-1,self.gids[self.start:self.end,:],0,self.comm)
        self.importer   = Epetra.Import(self.map1,self.map2)
        self.exporter   = Epetra.Export(self.map2,self.map1)
        self.numSameIDs = 0
        if self.myPID   == 0: self.numSameIDs = self.myN
        if self.numProc == 1: self.numSameIDs = self.myN * self.myN
        self.comm.Barrier()

    def tearDown(self):
        self.comm.Barrier()

    def testConstructors1(self):
        "Test Epetra.Import/Export constructors"
        source1 = self.importer.SourceMap()
        source2 = self.exporter.SourceMap()
        target1 = self.importer.TargetMap()
        target2 = self.exporter.TargetMap()
        self.assertEqual(source1.SameAs(source2), True )
        self.assertEqual(target1.SameAs(target2), True )
        if self.numProc > 1:
            self.assertEqual(target1.SameAs(source1), False)
            self.assertEqual(target2.SameAs(source2), False)

    def testConstructors2(self):
        "Test Epetra.Import/Export copy constructors"
        importer = Epetra.Import(self.importer)
        exporter = Epetra.Export(self.exporter)
        source1  = importer.SourceMap()
        source2  = exporter.SourceMap()
        target1  = importer.TargetMap()
        target2  = exporter.TargetMap()
        self.assertEqual(source1.SameAs(source2), True )
        self.assertEqual(target1.SameAs(target2), True )
        if self.numProc > 1:
            self.assertEqual(target1.SameAs(source1), False)
            self.assertEqual(target2.SameAs(source2), False)

    def testNumSameIDs(self):
        "Test Epetra.Import/Export NumSameIDs method"
        self.assertEqual(self.importer.NumSameIDs(), self.numSameIDs)
        self.assertEqual(self.exporter.NumSameIDs(), self.numSameIDs)

    def testNumPermuteIDs(self):
        "Test Epetra.Import/Export NumPermuteIDs method"
        numPermuteIDs = self.myN * self.myN - self.numSameIDs
        self.assertEqual(self.importer.NumPermuteIDs(), numPermuteIDs)
        self.assertEqual(self.exporter.NumPermuteIDs(), numPermuteIDs)

    def testNumRemoteIDs(self):
        "Test Epetra.Import/Export NumRemoteIDs method"
        numRemoteIDs  = self.myN * (self.globalN - self.myN)
        self.assertEqual(self.importer.NumRemoteIDs(), numRemoteIDs)
        self.assertEqual(self.exporter.NumRemoteIDs(), numRemoteIDs)

    def testNumExportIDs(self):
        "Test Epetra.Import/Export NumExportIDs method"
        numExportIDs  = self.myN * (self.globalN - self.myN)
        self.assertEqual(self.importer.NumExportIDs(), numExportIDs)
        self.assertEqual(self.exporter.NumExportIDs(), numExportIDs)

    def testNumSend(self):
        "Test Epetra.Import/Export NumSend method"
        numSend  = self.myN * (self.globalN - self.myN)
        self.assertEqual(self.importer.NumSend(), 0      )
        self.assertEqual(self.exporter.NumSend(), numSend)

    def testNumRecv(self):
        "Test Epetra.Import/Export NumRecv method"
        numRecv  = self.myN * (self.globalN - self.myN)
        self.assertEqual(self.importer.NumRecv(), numRecv)
        self.assertEqual(self.exporter.NumRecv(), numRecv)

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraImportExportTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
          "\n****************************\nTesting Epetra.Import/Export\n" + \
          "****************************\n"
    verbosity = 2 * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))[0]
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
