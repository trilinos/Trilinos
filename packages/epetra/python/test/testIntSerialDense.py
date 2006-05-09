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
from   numpy    import *

##########################################################################

class EpetraIntSerialDenseVectorTestCase(unittest.TestCase):
    "TestCase class for Epetra IntSerialDense Vectors"

    def setUp(self):
        self.comm = Epetra.PyComm()
        self.size = 4
        self.rows = 2
        self.cols = 4
        self.list = [0, -1, 2, -3]
        # Query whether Copy objects are actually Copy objects.  For older
        # versions of SWIG, they are View objects, but newer versions of SWIG
        # get it right
        isdv = Epetra.IntSerialDenseVector(self.size)
        self.CV = isdv.CV()
        if self.CV == Epetra.View:
            self.CVStr = "View"
            self.YNStr = "no"
        else:
            self.CVStr = "Copy"
            self.YNStr = "yes"

    def tearDown(self):
        self.comm.Barrier()

    def testVectorConstructor0(self):
        "Test Epetra.IntSerialDenseVector default constructor"
        isdv = Epetra.IntSerialDenseVector()
        self.assertEqual(isdv.CV(), Epetra.Copy)
        self.assertEqual(isdv.Length(), 0)

    def testVectorConstructor1(self):
        "Test Epetra.IntSerialDenseVector sized constructor"
        isdv = Epetra.IntSerialDenseVector(self.size)
        self.assertEqual(isdv.CV(), self.CV)
        self.assertEqual(isdv.Length(), self.size)

    def testVectorConstructor2(self):
        "Test Epetra.IntSerialDenseVector 1D list constructor"
        isdv = Epetra.IntSerialDenseVector(self.list)
        self.assertEqual(isdv.CV(), Epetra.View)
        self.assertEqual(isdv.Length(), len(self.list))
        self.failUnless((isdv == self.list).all())

    def testVectorConstructor3(self):
        "Test Epetra.IntSerialDenseVector 2D list constructor"
        list = [[0, -1, 2], [-3, 4, -5]]
        isdv = Epetra.IntSerialDenseVector(list)
        self.assertEqual(isdv.CV(), Epetra.View)
        self.assertEqual(isdv.Length(), len(list)*len(list[0]))
        self.failUnless((isdv == list).all())

    def testVectorConstructor4(self):
        "Test Epetra.IntSerialDenseVector copy constructor"
        isdv1 = Epetra.IntSerialDenseVector(self.size)
        isdv2 = Epetra.IntSerialDenseVector(isdv1)
        self.assertEqual(isdv2.CV(), self.CV)
        self.assertEqual(isdv1.Length(), isdv2.Length())

    def testVectorConstructor5(self):
        "Test Epetra.IntSerialDenseVector bad-list constructor"
        list = [0,1.0,"e","pi"]
        self.assertRaises(TypeError,Epetra.IntSerialDenseVector,list)

    def testVectorPrint(self):
        "Test Epetra.IntSerialDenseVector Print method"
        isdv     = Epetra.IntSerialDenseVector(self.size)
        isdv[:]  = 0
        filename = "testIntSerialDense%d.dat" % self.comm.MyPID()
        f = open(filename, "w")
        isdv.Print(f)
        f.close()
        out = "Data access mode: %s\nA_Copied: %s\nLength(M): %d\n" % \
              (self.CVStr,self.YNStr,self.size) + self.size * "0 " + "\n"
        f = open(filename, "r")
        self.assertEqual(f.read(), out)
        f.close()

    def testVectorStr(self):
        "Test Epetra.IntSerialDenseVector __str__ method"
        isdv = Epetra.IntSerialDenseVector(self.size)
        isdv[:] = 0
        out = "[" + self.size * "0 "
        out = out [:-1] + "]"
        self.assertEquals(str(isdv), out)

    def testVectorSize(self):
        "Test Epetra.IntSerialDenseVector Size method"
        isdv = Epetra.IntSerialDenseVector()
        isdv.Size(3*self.size)
        self.assertEqual(isdv.Length(), 3*self.size)
        self.assertEqual(isdv.Length(), len(isdv))
        self.failUnless((isdv == 0).all())

    def testVectorResize1(self):
        "Test Epetra.IntSerialDenseVector Resize method to smaller"
        isdv = Epetra.IntSerialDenseVector(3*self.size)
        self.assertEqual(isdv.Length(), 3*self.size)
        self.assertEqual(isdv.Length(), len(isdv))
        isdv[:] = range(len(isdv))
        isdv.Resize(self.size)
        self.assertEqual(isdv.Length(), self.size)
        self.assertEqual(isdv.Length(), len(isdv))
        for i in range(len(isdv)):
            self.assertEqual(isdv[i], i)

    def testVectorResize2(self):
        "Test Epetra.IntSerialDenseVector Resize method to larger"
        isdv = Epetra.IntSerialDenseVector(self.size)
        self.assertEqual(isdv.Length(), self.size)
        self.assertEqual(isdv.Length(), len(isdv))
        isdv[:] = range(len(isdv))
        isdv.Resize(3*self.size)
        self.assertEqual(isdv.Length(), 3*self.size)
        self.assertEqual(isdv.Length(), len(isdv))
        for i in range(len(isdv)):
            if i < self.size:
                self.assertEqual(isdv[i], i)
            else:
                self.assertEqual(isdv[i], 0.0)

    def testVectorGetitem1(self):
        "Test Epetra.IntSerialDenseVector __getitem__ method"
        isdv = Epetra.IntSerialDenseVector(self.list)
        for i in range(len(self.list)):
            self.assertEqual(isdv[i],self.list[i])

    def testVectorGetitem2(self):
        "Test Epetra.IntSerialDenseVector __call__ method"
        isdv = Epetra.IntSerialDenseVector(self.list)
        for i in range(len(self.list)):
            self.assertEqual(isdv(i),self.list[i])

    def testVectorSetitem(self):
        "Test Epetra.IntSerialDenseVector __setitem__ method"
        isdv = Epetra.IntSerialDenseVector(len(self.list))
        isdv[:] = 0
        self.failUnless((isdv == 0).all())
        for i in range(len(isdv)):
            isdv[i] = self.list[i]
        self.failUnless((isdv == self.list).all())

    def testVectorIndexErrors(self):
        "Test Epetra.IntSerialDenseVector index errors"
        isdv = Epetra.IntSerialDenseVector(self.size)
        self.assertRaises(TypeError, isdv.__getitem__, 0,1)
        self.assertRaises(TypeError, isdv.__setitem__, 0,1,3.14)

    def testVectorRandom(self):
        "Test Epetra.IntSerialDenseVector Random method"
        isdv    = Epetra.IntSerialDenseVector(self.size)
        isdv[:] = 0
        self.failUnless((isdv == 0).all())
        result = isdv.Random()
        self.assertEqual(result, 0)
        maxInt = 2**31 - 1
        for i in range(self.size):
            self.failUnless(isdv[i] >= 0     )
            self.failUnless(isdv[i] <= maxInt)

    def testVectorLength(self):
        "Test Epetra.IntSerialDenseVector Length method"
        isdv = Epetra.IntSerialDenseVector()
        self.assertEqual(isdv.Length(),0)
        self.assertEqual(isdv.Length(),len(isdv))
        isdv.Size(self.size)
        self.assertEqual(isdv.Length(),self.size)
        self.assertEqual(isdv.Length(),len(isdv))
        isdv.Resize(2*self.size)
        self.assertEqual(isdv.Length(),2*self.size)
        self.assertEqual(isdv.Length(),len(isdv))
        isdv.Resize(self.size)
        self.assertEqual(isdv.Length(),self.size)
        self.assertEqual(isdv.Length(),len(isdv))

    def testVectorValues(self):
        "Test Epetra.IntSerialDenseVector Values method"
        isdv  = Epetra.IntSerialDenseVector(self.list)
        vals = isdv.Values()
        self.assertNotEqual(type(isdv), type(vals))
        self.assertEqual(isdv.Length(), len(vals))
        for i in range(len(vals)):
            self.assertEqual(isdv[i], vals[i])

    def testVectorInPlaceAdd(self):
        "Test Epetra.IntSerialDenseVector += operator"
        isdv = Epetra.IntSerialDenseVector(self.list)
        self.failUnless((isdv == self.list).all())
        isdv += isdv
        for i in range(isdv.Length()):
            self.assertEqual(isdv[i],2*self.list[i])

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraIntSerialDenseVectorTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n*****************************" \
       "\nTesting Epetra.IntSerialDense" \
       "\n*****************************\n"
    verbosity = 2 * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
