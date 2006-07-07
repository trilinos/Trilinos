#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#              PyTrilinos.Galeri: Python Interface to Galeri
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
    import Galeri
else:
    try:
        import setpath
        import Epetra
        import Galeri
    except ImportError:
        from PyTrilinos import Epetra, Galeri
        print >>sys.stderr, "Using system-installed Epetra, Galeri"

##########################################################################

class CreateMapTestCase(unittest.TestCase):
    "TestCase for Galeri CreateMap function"

    def setUp(self):
        self.comm   = Epetra.PyComm()
        nx = 4 * self.comm.NumProc()
        ny = nz = 16
        mx = self.comm.NumProc()
        my = mz = 1
        self.param1 = {"n"  : nx }
        self.param2 = {"n"  : nx*ny,
                       "nx" : nx,
                       "ny" : ny,
                       "mx" : mx,
                       "my" : my    }
        self.param3 = {"n"  : nx*ny*nz,
                       "nx" : nx,
                       "ny" : ny,
                       "nz" : nz,
                       "mx" : mx,
                       "my" : my,
                       "mz" : mz }
        self.comm.Barrier()

    def tearDown(self):
        self.comm.Barrier()

    def testLinear(self):
        "Test Galeri CreateMap for type 'Linear'"
        map = Galeri.CreateMap("Linear", self.comm, self.param1)
        self.assertEqual(isinstance(map, Epetra.Map), True)

    def testCartesian2D(self):
        "Test Galeri CreateMap for type 'Cartesian2D'"
        map = Galeri.CreateMap("Cartesian2D", self.comm, self.param2)
        self.assertEqual(isinstance(map, Epetra.Map), True)

    def testCartesian3D(self):
        "Test Galeri CreateMap for type 'Cartesian3D'"
        map = Galeri.CreateMap("Cartesian3D", self.comm, self.param3)
        self.assertEqual(isinstance(map, Epetra.Map), True)

    def testRandom(self):
        "Test Galeri CreateMap for type 'Random'"
        map = Galeri.CreateMap("Random", self.comm, self.param1)
        self.assertEqual(isinstance(map, Epetra.Map), True)

    def testInterlaced(self):
        "Test Galeri CreateMap for type 'Interlaced'"
        map = Galeri.CreateMap("Interlaced", self.comm, self.param1)
        self.assertEqual(isinstance(map, Epetra.Map), True)

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(CreateMapTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, "\n************************" \
       "\nTesting Galeri.CreateMap\n************************\n"
    verbosity = 2 * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
