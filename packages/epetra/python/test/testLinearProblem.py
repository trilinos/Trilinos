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
from   numpy    import *
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
    import Epetra
else:
    try:
        import setpath
        import Epetra
    except ImportError:
        from PyTrilinos import Epetra
        print >>sys.stderr, "Using system-installed Epetra"

##########################################################################

class EpetraLinearProblemTestCase(unittest.TestCase):
    "TestCase for Epetra_LinearProblems"

    def setUp(self):
        self.comm          = Epetra.PyComm()
        self.map           = Epetra.Map(9,0,self.comm)
        self.A             = Epetra.CrsMatrix(Epetra.Copy,self.map,3)
        self.x             = Epetra.Vector(self.map)
        self.b             = Epetra.Vector(self.map)
        self.linearProblem = Epetra.LinearProblem(self.A,self.x,self.b)
        self.comm.Barrier()

    def tearDown(self):
        self.comm.Barrier()

    # The testing here is minimal so far, but it does test the return type when
    # an Epetra_MultiVector * is returned.  This is significant because it tests
    # a special typemap converting Epetra_MultiVector -> Epetra_NumPyMultiVector
    # and the special swig register that converts an Epetra_NumPyMultiVector to
    # a python Epetra.MultiVector (when using swig 1.3.28 or better).

    def testGetLHS(self):
        "Test Epetra.LinearProblem GetLHS method"
        lhs = self.linearProblem.GetLHS()
        if hasattr(Epetra,'NumPyMultiVectorPtr'):
            self.assertEqual(isinstance(lhs,Epetra.NumPyMultiVectorPtr), True)
        else:
            self.failUnless(isinstance(lhs,Epetra.MultiVector))
            self.assertEqual(lhs[0], self.x[0])

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraLinearProblemTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
          "\n****************************\nTesting Epetra.LinearProblem\n" \
          "****************************\n"
    verbosity = options.verbosity * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
