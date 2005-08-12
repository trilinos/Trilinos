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

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraMapTestCase))

    # Run the test suite
    print >>sys.stderr, \
          "\n******************\nTesting Epetra.Map\n******************\n"
    result = unittest.TextTestRunner(verbosity=2).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    sys.exit(len(result.errors) + len(result.failures))
