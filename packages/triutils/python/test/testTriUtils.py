#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#         PyTrilinos.Triutils: Python Interface to Triutils
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
from   optparse import *
import sys
import unittest

parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
options,args = parser.parse_args()
if options.testharness:
    import setpath
    import Epetra, Triutils
else:
    try:
        import setpath
        import Epetra, Triutils
    except ImportError:
        from PyTrilinos import Epetra, Triutils
        print >>sys.stderr, "Using system-installed Epetra, Triutils"

####################################################################

class TriutilsTestCase(unittest.TestCase):
    "TestCase class for Triutils module"

    def setUp(self):
        self.comm    = Epetra.PyComm()
        self.numProc = self.comm.NumProc()
        self.gallery = Triutils.CrsMatrixGallery("laplace_2d",self.comm)
        self.gallery.Set("nx",10*self.numProc)
        self.gallery.Set("ny",10*self.numProc)
        self.comm.Barrier()

    def tearDown(self):
        # This will help tame the printing
        self.comm.Barrier()

    def testVersion(self):
        "Test Triutils Triutils_Version function"
        front   = "Triutils Version "
        version = Triutils.Triutils_Version()
        self.assertEquals(version[:len(front)], front)

    def testGetRHS(self):
        "Test Triutils.CrsMatrixGallery GetRHS method"
        rhs = self.gallery.GetRHS()
        # This is the result if swig version < 1.3.28
        if hasattr(Epetra,'NumPyMultiVectorPtr'):
            self.assertEqual(isinstance(rhs,Epetra.NumPyMultiVectorPtr), True)
        # This is the result if swig version >= 1.3.28
        else:
            self.assertEqual(isinstance(rhs,Epetra.MultiVector), True)

####################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(TriutilsTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n*******************\nTesting Triutils\n*******************\n"
    verbosity = options.verbosity * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
