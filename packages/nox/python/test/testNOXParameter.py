#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                 PyTrilinos.NOX: Python Interface to NOX
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
    import Epetra
    import NOX
else:
    try:
        import setpath
        import NOX
        import Epetra
    except ImportError:
        from PyTrilinos import Epetra
        from PyTrilinos import NOX
        print >>sys.stderr, "Using system-installed Epetra, NOX"

################################################################################

class ParameterTestCase(unittest.TestCase):
    "TestCase for NOX.Parameter.Lists"

    def setUp(self):
        self.comm = Epetra.PyComm()
        self.pl   = NOX.Parameter.List()
        self.pl.setParameter("Integer Parameter", 2003    )
        self.pl.setParameter("Float Parameter",   3.14    )
        self.pl.setParameter("String Parameter",  "Sandia")

    def tearDown(self):
        self.comm.Barrier()

    def testIntExistence(self):
        "Test NOX.Parameter.List integer existence"
        self.assertEqual(self.pl.isParameter("Integer Parameter" ), True )

    def testFloatExistence(self):
        "Test NOX.Parameter.List float existence"
        self.assertEqual(self.pl.isParameter("Float Parameter"   ), True )

    def testStringExistence(self):
        "Test NOX.Parameter.List string existence"
        self.assertEqual(self.pl.isParameter("String Parameter"  ), True )

    def testFakeExistence(self):
        "Test NOX.Parameter.List fake existence"
        self.assertEqual(self.pl.isParameter("Fake Parameter"    ), False)

    def testIntValues(self):
        "Test NOX.Parameter.List integer value"
        self.assertEqual(self.pl.getParameter("Integer Parameter"), 2003)

    def testFloatValues(self):
        "Test NOX.Parameter.List float value"
        self.assertEqual(self.pl.getParameter("Float Parameter"  ), 3.14)

    def testStringValues(self):
        "Test NOX.Parameter.List string value"
        self.assertEqual(self.pl.getParameter("String Parameter" ), "Sandia")

    def testSubList(self):
        "Test NOX.Parameter.List sublists"
        sl1 = self.pl.sublist("Sublist 1")
        self.assertEquals(self.pl.isParameter("Sublist 1"), True)
        sl2 = sl1.sublist("Sublist 2")
        self.assertEquals(sl1.isParameter("Sublist 2"), True)
        outputInfo = NOX.Utils.Warning                  + \
                     NOX.Utils.OuterIteration           + \
                     NOX.Utils.InnerIteration           + \
                     NOX.Utils.Parameters               + \
                     NOX.Utils.Details                  + \
                     NOX.Utils.OuterIterationStatusTest
        sl2.setParameter("New Parameter", outputInfo)
        self.assertEquals(sl2.getParameter("New Parameter"), outputInfo)

################################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(ParameterTestCase ))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n*********************\nTesting NOX.Parameter\n*********************\n"
    verbosity = options.verbosity * int(iAmRoot)
    result    = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
