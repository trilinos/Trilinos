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

class EpetraMapTestCase(unittest.TestCase):
    "TestCase class for Map objects"

    def setUp(self):
        self.comm        = Epetra.PyComm()
        self.numLocalEl  = 10
        self.numGlobalEl = self.numLocalEl * self.comm.NumProc()
        self.map         = Epetra.Map(self.numGlobalEl,0,self.comm)

    def tearDown(self):
        self.comm.Barrier()

    def testNumGlobalElements(self):
        "Test Epetra.Map NumGlobalElements method"
        self.assertEqual(self.map.NumGlobalElements(), self.numGlobalEl)

    def testNumMyElements(self):
        "Test Epetra.Map NumMyElements method"
        self.assertEqual(self.map.NumMyElements(), self.numLocalEl)

    def testElementSize(self):
        "Test Epetra.Map ElementSize method"
        self.assertEqual(self.map.ElementSize(), 1)

    def testIndexBase(self):
        "Test Epetra.Map IndexBase method"
        self.assertEqual(self.map.IndexBase(), 0)

    def testNumGlobalPoints(self):
        "Test Epetra.Map NumGlobalPoints method"
        self.assertEqual(self.map.NumGlobalPoints(), self.numGlobalEl)

    def testNumMyPoints(self):
        "Test Epetra.Map NumMyPoints method"
        self.assertEqual(self.map.NumMyPoints(), self.numLocalEl)

    def testMinMyElementSize(self):
        "Test Epetra.Map MinMyElementSize method"
        self.assertEqual(self.map.MinMyElementSize(), 1)

    def testMaxMyElementSize(self):
        "Test Epetra.Map MaxMyElementSize method"
        self.assertEqual(self.map.MaxMyElementSize(), 1)

    def testMinElementSize(self):
        "Test Epetra.Map MinElementSize method"
        self.assertEqual(self.map.MinElementSize(), 1)

    def testMaxElementSize(self):
        "Test Epetra.Map MaxElementSize method"
        self.assertEqual(self.map.MaxElementSize(), 1)

    def testConstantElementSize(self):
        "Test Epetra.Map ConstantElementSize method"
        self.assertEqual(self.map.ConstantElementSize(), True)

    def testSameAs(self):
        "Test Epetra.Map SameAs method"
        self.assertEqual(self.map.SameAs(self.map), True)

    def testDistributedGlobal(self):
        "Test Epetra.Map DistributedGlobal method"
        distributedGlobal = (self.comm.Label() == "Epetra::MpiComm" and
                             self.comm.NumProc() > 1)
        self.assertEqual(self.map.DistributedGlobal(), distributedGlobal)

    def testMinAllGID(self):
        "Test Epetra.Map MinAllGID method"
        self.assertEqual(self.map.MinAllGID(), 0)

    def testMaxAllGID(self):
        "Test Epetra.Map MaxAllGID method"
        self.assertEqual(self.map.MaxAllGID(), self.numGlobalEl-1)

    def testMinMyGID(self):
        "Test Epetra.Map MinMyGID method"
        self.assertEqual(self.map.MinMyGID(), self.comm.MyPID()*self.numLocalEl)

    def testMaxMyGID(self):
        "Test Epetra.Map MaxMyGID method"
        self.assertEqual(self.map.MaxMyGID(), (self.comm.MyPID()+1)*self.numLocalEl-1)

    def testMinLID(self):
        "Test Epetra.Map MinLID method"
        self.assertEqual(self.map.MinLID(), 0)

    def testMaxLID(self):
        "Test Epetra.Map MaxLID method"
        self.assertEqual(self.map.MaxLID(), self.numLocalEl-1)

    def testIDs(self):
        "Test Epetra.Map global and local IDs"
        myMinID = self.comm.MyPID() * self.numLocalEl
        myMaxID = myMinID + self.numLocalEl - 1
        for gid in range(self.map.NumGlobalElements()):
            if myMinID <= gid and gid <= myMaxID:
                lid = gid % self.numLocalEl
            else:
                lid = -1
            self.assertEqual(self.map.LID(gid)  , lid        )
            self.assertEqual(self.map.MyGID(gid), (lid != -1))
            self.assertEqual(self.map.MyLID(lid), (lid != -1))

    def testStr(self):
        "Test Epetra.Map __str__ method"
        lines   = 7 + self.numLocalEl
        if self.comm.MyPID() == 0: lines += 7
        s = str(self.map)
        s = s.splitlines()
        self.assertEquals(len(s), lines)

    def testPrint(self):
        "Test Epetra.Map Print method"
        myPID = self.comm.MyPID()
        filename = "testMap%d.dat" % myPID
        f = open(filename, "w")
        self.map.Print(f)
        f.close()
        f = open(filename, "r")
        s = f.readlines()
        f.close()
        lines = 7 + self.numLocalEl
        if myPID == 0: lines += 7
        self.assertEquals(len(s), lines)

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraMapTestCase))

    # Create a communicator
    comm = Epetra.PyComm()

    # Run the test suite
    if comm.MyPID() == 0: print >>sys.stderr, \
          "\n******************\nTesting Epetra.Map\n******************\n"
    verbosity = 2 * int(comm.MyPID() == 0)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    sys.exit(len(result.errors) + len(result.failures))
