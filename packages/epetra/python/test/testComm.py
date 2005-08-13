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

import os
import unittest
from   Numeric  import *

##########################################################################

class EpetraPyCommTestCase(unittest.TestCase):
    "TestCase class for SerialComm communicator objects"

    def setUp(self):
        self.comm  = Epetra.PyComm()
        if self.comm.Label() == "Epetra::MpiComm":
            self.commType = "MPI"
            self.output   = "  Processor %d of %d total processors" % \
                            (self.comm.MyPID(), self.comm.NumProc()) 
        else:
            self.commType = "Serial"
            self.output   = "::Processor 0 of 1 total processors."
        self.comm.Barrier()

    def tearDown(self):
        self.comm.Barrier()

    def testMyPID(self):
        "Test Epetra.PyComm MyPID method"
        pid = self.comm.MyPID()
        if self.commType == "Serial":
            self.assertEqual(pid, 0)
        else:
            self.assert_(pid >= 0)

    def testNumProc(self):
        "Test Epetra.PyComm NumProc method"
        n = self.comm.NumProc()
        if self.commType == "Serial":
            self.assertEqual(n, 1)
        else:
            self.assert_(n > 0)

    def testStr(self):
        "Test Epetra.PyComm __str__ method"
        self.assertEqual(str(self.comm), self.output)

    def testPrint(self):
        "Test Epetra.PyComm Print method"
        filename = "testComm%d.dat" % self.comm.MyPID()
        f = open(filename, "w")
        self.comm.Print(f)
        f.close()
        f = open(filename, "r")
        self.assertEqual(f.read(), self.output)
        f.close()

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraPyCommTestCase ))

    # Create a communicator
    comm = Epetra.PyComm()
        
    # Run the test suite
    if comm.MyPID() == 0: print >>sys.stderr, \
          "\n*******************\nTesting Epetra.Comm\n*******************\n"
    verbosity = 2 * int(comm.MyPID() == 0)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    sys.exit(len(result.errors) + len(result.failures))
