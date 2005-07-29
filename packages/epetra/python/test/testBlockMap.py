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

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraBlockMapTestCase))

    # Run the test suite
    print >>sys.stderr, \
          "\n***********************\nTesting Epetra.BlockMap\n***********************\n"
    unittest.TextTestRunner(verbosity=2).run(suite)
