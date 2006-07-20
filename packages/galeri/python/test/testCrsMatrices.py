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
        from PyTrilinos import Teuchos, Epetra, Galeri
        print >>sys.stderr, "Using system-installed Epetra, Galeri"

##########################################################################

class CreateCrsMatrixTestCase(unittest.TestCase):
    "TestCase for Galeri CreateCrsMatrix function"

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
        self.map1   = Galeri.CreateMap("Interlaced", self.comm, self.param1)
        self.map2   = Galeri.CreateMap("Random"    , self.comm, self.param2)
        self.map3   = Galeri.CreateMap("Random"    , self.comm, self.param3)
        self.comm.Barrier()

    def tearDown(self):
        self.comm.Barrier()

    def testDiag(self):
        "Test Galeri CreateCrsMatrix for problem 'Diag'"
        matrix = Galeri.CreateCrsMatrix("Diag", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testTridiag(self):
        "Test Galeri CreateCrsMatrix for problem 'Tridiag'"
        matrix = Galeri.CreateCrsMatrix("Tridiag", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testLaplace1D(self):
        "Test Galeri CreateCrsMatrix for problem 'Laplace1D'"
        matrix = Galeri.CreateCrsMatrix("Laplace1D", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testLaplace1DNeumann(self):
        "Test Galeri CreateCrsMatrix for problem 'Laplace1DNeumann'"
        matrix = Galeri.CreateCrsMatrix("Laplace1DNeumann", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testCauchy(self):
        "Test Galeri CreateCrsMatrix for problem 'Cauchy'"
        matrix = Galeri.CreateCrsMatrix("Cauchy", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testFielder(self):
        "Test Galeri CreateCrsMatrix for problem 'Fielder'"
        matrix = Galeri.CreateCrsMatrix("Fielder", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testHanowa(self):
        "Test Galeri CreateCrsMatrix for problem 'Hanowa'"
        matrix = Galeri.CreateCrsMatrix("Hanowa", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testHilbert(self):
        "Test Galeri CreateCrsMatrix for problem 'Hilbert'"
        matrix = Galeri.CreateCrsMatrix("Hilbert", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testJordanBlock(self):
        "Test Galeri CreateCrsMatrix for problem 'JordanBlock'"
        matrix = Galeri.CreateCrsMatrix("JordanBlock", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testKMS(self):
        "Test Galeri CreateCrsMatrix for problem 'KMS'"
        matrix = Galeri.CreateCrsMatrix("KMS", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testLehmer(self):
        "Test Galeri CreateCrsMatrix for problem 'Lehmer'"
        matrix = Galeri.CreateCrsMatrix("Lehmer", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testOnes(self):
        "Test Galeri CreateCrsMatrix for problem 'Ones'"
        matrix = Galeri.CreateCrsMatrix("Ones", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testPei(self):
        "Test Galeri CreateCrsMatrix for problem 'Pei'"
        matrix = Galeri.CreateCrsMatrix("Pei", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testRis(self):
        "Test Galeri CreateCrsMatrix for problem 'Ris'"
        matrix = Galeri.CreateCrsMatrix("Ris", self.map1, self.param1)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testCross2D(self):
        "Test Galeri CreateCrsMatrix for problem 'Cross2D'"
        matrix = Galeri.CreateCrsMatrix("Cross2D", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testBigCross2D(self):
        "Test Galeri CreateCrsMatrix for problem 'BigCross2D'"
        matrix = Galeri.CreateCrsMatrix("BigCross2D", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testStar2D(self):
        "Test Galeri CreateCrsMatrix for problem 'Star2D'"
        matrix = Galeri.CreateCrsMatrix("Star2D", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testBigStar2D(self):
        "Test Galeri CreateCrsMatrix for problem 'BigStar2D'"
        matrix = Galeri.CreateCrsMatrix("BigStar2D", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testLaplace2D(self):
        "Test Galeri CreateCrsMatrix for problem 'Laplace2D'"
        matrix = Galeri.CreateCrsMatrix("Laplace2D", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testStretched2D(self):
        "Test Galeri CreateCrsMatrix for problem 'Stretched2D'"
        matrix = Galeri.CreateCrsMatrix("Stretched2D", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testUniFlow2D(self):
        "Test Galeri CreateCrsMatrix for problem 'UniFlow2D'"
        matrix = Galeri.CreateCrsMatrix("UniFlow2D", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testRecirc2D(self):
        "Test Galeri CreateCrsMatrix for problem 'Recirc2D'"
        matrix = Galeri.CreateCrsMatrix("Recirc2D", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testBiharmonic2D(self):
        "Test Galeri CreateCrsMatrix for problem 'Biharmonic2D'"
        matrix = Galeri.CreateCrsMatrix("Biharmonic2D", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testLaplace2DFourthOrder(self):
        "Test Galeri CreateCrsMatrix for problem 'Laplace2DFourthOrder'"
        matrix = Galeri.CreateCrsMatrix("Laplace2DFourthOrder", self.map2, self.param2)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

    def testLaplace3D(self):
        "Test Galeri CreateCrsMatrix for problem 'Laplace3D'"
        matrix = Galeri.CreateCrsMatrix("Laplace3D", self.map3, self.param3)
        self.assertEqual(isinstance(matrix, Epetra.CrsMatrix), True)

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(CreateCrsMatrixTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, "\n******************************" \
       "\nTesting Galeri.CreateCrsMatrix\n******************************\n"
    verbosity = options.verbosity * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
