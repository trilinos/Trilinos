#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                PyTrilinos: Python Interface to Trilinos
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

class EpetraVectorTestCase(unittest.TestCase):
    "TestCase class for Vector objects"

    def setUp(self):
        self.length      = 9
        self.scale       = 1.0 / (self.length-1)
        self.comm        = Epetra.PyComm()
        self.map         = Epetra.Map(self.length*self.comm.NumProc(),0,self.comm)
        self.numPyArray1 = arange(self.length) * self.scale
        self.numPyArray2 = array([0,-1,2.17,-3.14,4,-5,6,-7,6.28])
        
    def tearDown(self):
        self.comm.Barrier()

    def testConstructor01(self):
        "Test Epetra.FEVector (BlockMap,bool) constructor"
        ev = Epetra.FEVector(self.map,True)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        self.assertEquals((ev == 0.0).all(), True)

    def testConstructor02(self):
        "Test Epetra.FEVector (BlockMap) constructor"
        ev = Epetra.FEVector(self.map)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())

    def testConstructor03(self):
        "Test Epetra.FEVector copy constructor"
        ev1 = Epetra.FEVector(self.map)
        ev2 = Epetra.FEVector(ev1)
        self.assertEquals(ev2.NumVectors(),   ev1.NumVectors()  )
        self.assertEquals(ev2.MyLength(),     ev1.MyLength()    )
        self.assertEquals(ev2.GlobalLength(), ev1.GlobalLength())
        for i in range(len(ev1)):
            self.assertEquals(ev1[i], ev2[i])

    def testReplaceMap1(self):
        "Test Epetra.FEVector ReplaceMap method with good map"
        blockMap = Epetra.BlockMap(3*self.comm.NumProc(),3,0,self.comm)
        ev = Epetra.FEVector(self.map)
        result = ev.ReplaceMap(blockMap)
        self.assertEquals(result, 0)
        newMap = ev.Map()
        self.assertEquals(newMap.ElementSize(), blockMap.ElementSize())

    def testReplaceMap2(self):
        "Test Epetra.FEVector ReplaceMap method with bad map"
        blockMap = Epetra.BlockMap(2*self.comm.NumProc(),5,0,self.comm)
        ev = Epetra.FEVector(self.map)
        result = ev.ReplaceMap(blockMap)
        self.assertEquals(result, -1)
        newMap = ev.Map()
        self.assertEquals(newMap.ElementSize(), self.map.ElementSize())

    def testReplaceGlobalValue1(self):
        "Test Epetra.FEVector ReplaceGlobalValue method"
        ev = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        gid = 4
        lid = self.map.LID(gid)
        self.assertEquals(ev[gid], 0.5)
        result = ev.ReplaceGlobalValue(gid,0,5.0)
        if lid >= 0:
            self.assertEquals(result, 0)
            self.assertEquals(ev[lid], 5.0)
        else:
            self.assertEquals(result, 1)

    def testReplaceGlobalValues1(self):
        "Test Epetra.FEVector ReplaceGlobalValues method"
        ev = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        gids = [2,3,6]
        lids = [self.map.LID(gid) for gid in gids]
        result = ev.ReplaceGlobalValues(gids,[2.2,3.3,6.6])
        for i in range(len(gids)):
            gid = gids[i]
            lid = lids[i]
            if lid >= 0:
                self.assertEquals(result, 0)
                self.assertAlmostEquals(ev[lid], 1.1*gid)
            else:
                self.assertEquals(result, 1)

    def testReplaceMyValue1(self):
        "Test Epetra.FEVector ReplaceMyValue method"
        ev = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        lid = 4
        self.assertEquals(ev[lid], 0.5)
        result = ev.ReplaceMyValue(lid,0,5.0)
        self.assertEquals(result, 0)
        self.assertEquals(ev[lid], 5.0)

    def testSumIntoGlobalValue1(self):
        "Test Epetra.FEVector SumIntoGlobalValue method"
        ev = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        gid = 4
        lid = self.map.LID(gid)
        self.assertEquals(ev[gid], 0.5)
        result = ev.SumIntoGlobalValue(gid,0,0.5)
        if lid >= 0:
            self.assertEquals(result, 0)
            self.assertEquals(ev[lid], 1.0)
        else:
            self.assertEquals(result, 1)

    def testSumIntoGlobalValues1(self):
        "Test Epetra.FEVector SumIntoGlobalValues method"
        ev = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        gids = [2,3,6]
        lids = [self.map.LID(gid) for gid in gids]
        result = ev.SumIntoGlobalValues(gids,[2.2,3.3,6.6])
        for i in range(len(gids)):
            gid = gids[i]
            lid = lids[i]
            if lid >= 0:
                self.assertEquals(result, 0)
                self.assertAlmostEquals(ev[lid], 0.125*lid+1.1*gid)
            else:
                self.assertEquals(result, 1)

    def testSumIntoMyValue1(self):
        "Test Epetra.FEVector SumIntoMyValue method"
        ev = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        lid = 4
        self.assertEquals(ev[lid], 0.5)
        result = ev.SumIntoMyValue(lid,0,0.5)
        self.assertEquals(result, 0)
        self.assertEquals(ev[lid], 1.0)

    def testPutScalar(self):
        "Test Epetra.FEVector PutScalar method"
        ev = Epetra.FEVector(self.map)
        for i in range(self.map.NumMyPoints()):
            self.assertEquals(ev[i], 0.0)
        scalar = 3.14
        ev.PutScalar(scalar)
        for i in range(self.map.NumMyPoints()):
            self.assertEquals(ev[i], scalar)

    def testRandom(self):
        "Test Epetra.FEVector Random method"
        ev = Epetra.FEVector(self.map)
        scalar = 3.14
        ev.PutScalar(scalar)
        for i in range(self.map.NumMyPoints()):
            self.assertEquals(ev[i], scalar)
        ev.Random()
        for i in range(self.map.NumMyPoints()):
            self.assertEquals(ev[i]>-1.0, True)
            self.assertEquals(ev[i]< 1.0, True)

    def testExtractCopy(self):
        "Test Epetra.FEVector ExtractCopy method"
        ev    = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        array = ev.ExtractCopy()
        self.assertEquals(type(array), ndarray)
        self.failUnless((ev.array == array).all())
        self.assertEquals(ev.array is array, False)

    def testExtractView(self):
        "Test Epetra.FEVector ExtractView method"
        ev    = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        array = ev.ExtractView()
        self.assertEquals(type(array), ndarray)
        self.failUnless((ev.array == array).all())
        self.assertEquals(ev.array is array, True)

    def testNumVectors(self):
        "Test Epetra.FEVector NumVectors method"
        ev = Epetra.FEVector(self.map)
        self.assertEquals(ev.NumVectors(), 1)

    def testMyLength(self):
        "Test Epetra.FEVector MyLength method"
        ev = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        self.assertEquals(ev.MyLength(), self.length)

    def testGlobalLength(self):
        "Test Epetra.FEVector GlobalLength method"
        ev = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray1
        self.assertEquals(ev.GlobalLength(), self.length*self.comm.NumProc())

    def testSeed(self):
        "Test Epetra.FEVector Seed method"
        ev   = Epetra.FEVector(self.map)
        seed = ev.Seed()
        max  = 2**31 - 1
        self.assertEquals(seed>0,   True)
        self.assertEquals(seed<max, True)

    def testSetSeed1(self):
        "Test Epetra.FEVector SetSeed method"
        ev     = Epetra.FEVector(self.map)
        seed   = 2005
        result = ev.SetSeed(seed)
        self.assertEquals(result,    0   )
        self.assertEquals(ev.Seed(), seed)

    def testSetSeed2(self):
        "Test Epetra.FEVector SetSeed method with negative seed"
        ev   = Epetra.FEVector(self.map)
        seed = -2005
        # The exception type depends on the version of SWIG used to generate the
        # wrapper
        self.assertRaises((TypeError,OverflowError),ev.SetSeed,seed)

    def testPrint(self):
        "Test Epetra.FEVector Print method"
        output = ""
        if self.comm.MyPID() == 0:
            output += "%10s%14s%20s  \n" % ("MyPID","GID","Value")
        for lid in range(self.length):
            gid = self.map.GID(lid)
            output += "%10d%14d%24d\n" % (self.comm.MyPID(),gid,0)
        ev = Epetra.FEVector(self.map)
        filename = "testVector%d.dat" % self.comm.MyPID()
        f = open(filename,"w")
        ev.Print(f)
        f.close()
        self.assertEqual(open(filename,"r").read(), output)

    def testDot(self):
        "Test Epetra.FEVector Dot method"
        map    = Epetra.Map(8*self.comm.NumProc(),0,self.comm)
        array1 = [-1, 2,-3, 4, 5, 1,-8,-7]
        array2 = [ 9, 0,-1,-2,-7,-8, 1, 5]
        ev1    = Epetra.FEVector(map)
        ev1[:] = array1
        ev2    = Epetra.FEVector(map)
        ev2[:] = array2
        dot    = ev1.Dot(ev2)
        result = -100*self.comm.NumProc()
        self.assertEqual(dot, result)

#     def testAbs(self):
#         "Test Epetra.FEVector Abs method"
#         ev1 = Epetra.FEVector(self.map,self.numPyArray1)
#         ev2 = Epetra.FEVector(self.map)
#         self.assertEquals(ev2[:],0.0)
#         result = ev2.Abs(ev1)
#         self.assertEquals(result,0)
#         self.assertEquals(ev2[:],abs(self.numPyArray1))

    def testReciprocal(self):
        "Test Epetra.FEVector Reciprocal method"
        a    = self.numPyArray1
        a[0] = 1.0  # Don't want to invert zero
        ev1  = Epetra.FEVector(self.map)
        ev1[:] = a
        ev2  = Epetra.FEVector(self.map)
        self.failUnless((ev2.array == 0.0).all())
        result = ev2.Reciprocal(ev1)
        self.assertEquals(result,0)
        self.failUnless((ev2.array == 1.0/a).all())

    def testScale1(self):
        "Test Epetra.FEVector Scale method in-place"
        a  = self.numPyArray2.copy()
        ev = Epetra.FEVector(self.map)
        ev[:] = self.numPyArray2
        result = ev.Scale(2.0)
        self.assertEquals(result,0)
        self.failUnless((abs(ev.array - 2.0*a) < 1e-10).all())

    def testScale2(self):
        "Test Epetra.FEVector Scale method with replace"
        a   = self.numPyArray1
        ev1 = Epetra.FEVector(self.map)
        ev1[:] = a
        ev2 = Epetra.FEVector(self.map)
        self.failUnless((ev2.array == 0.0).all())
        result = ev2.Scale(pi,ev1)
        self.assertEquals(result,0)
        self.failUnless((ev2.array == pi*a).all())

    def testUpdate1(self):
        "Test Epetra.FEVector Update method with one Vector"
        a   = self.numPyArray1.copy()
        b   = self.numPyArray2.copy()
        ev1 = Epetra.FEVector(self.map)
        ev1[:] = self.numPyArray1
        ev2 = Epetra.FEVector(self.map)
        ev2[:] = self.numPyArray2
        result = ev2.Update(2.0,ev1,3.0)
        self.assertEquals(result,0)
        self.failUnless((abs(ev2.array - (2.0*a + 3.0*b)) < 1e-10).all())

    def testUpdate2(self):
        "Test Epetra.FEVector Update method with two Vectors"
        ev0 = Epetra.FEVector(self.map)
        ev1 = Epetra.FEVector(self.map)
        ev1[:] = self.numPyArray1
        ev2 = Epetra.FEVector(self.map)
        ev2[:] = self.numPyArray2
        result = ev0.Update(2.0,ev1,3.0,ev2,pi)
        self.assertEquals(result,0)
        self.failUnless((ev0.array == 2.0*self.numPyArray1 + 3.0*self.numPyArray2).all())

    def testNorm1(self):
        "Test Epetra.FEVector Norm1 method"
        a      = self.numPyArray1
        ev     = Epetra.FEVector(self.map)
        ev[:]  = a
        result = sum(self.numPyArray1) * self.comm.NumProc()
        norm1  = ev.Norm1()
        self.assertEquals(norm1, result)

    def testNorm2(self):
        "Test Epetra.FEVector Norm2 method"
        a      = self.numPyArray2
        ev     = Epetra.FEVector(self.map)
        ev[:] = a
        result = sqrt(sum(self.numPyArray2*self.numPyArray2) * self.comm.NumProc())
        norm2  = ev.Norm2()
        self.assertEquals(norm2, result)

    def testNormInf(self):
        "Test Epetra.FEVector NormInf method"
        a       = self.numPyArray1
        ev      = Epetra.FEVector(self.map)
        ev[:]   = a
        result  = max(abs(self.numPyArray1))
        normInf = ev.NormInf()
        self.assertEquals(normInf, result)

    def testNormWeighted(self):
        "Test Epetra.FEVector NormWeighted method"
        a       = self.numPyArray2
        ev      = Epetra.FEVector(self.map)
        ev[:] = a
        wts     = sin(pi*(arange(self.length) + 0.5) / self.length)
        weights = Epetra.FEVector(self.map)
        weights[:] = wts
        result  = sqrt(sum((a/wts)**2)/self.length)
        norm    = ev.NormWeighted(weights)
        self.assertEquals(norm, result)

    def testMinValue(self):
        "Test Epetra.FEVector MinValue method"
        a        = self.numPyArray1
        ev       = Epetra.FEVector(self.map)
        ev[:]    = a
        result   = min(self.numPyArray1)
        minValue = ev.MinValue()
        self.assertEquals(minValue, result)

    def testMaxValue(self):
        "Test Epetra.FEVector MaxValue method"
        a        = self.numPyArray2
        ev       = Epetra.FEVector(self.map)
        ev[:]    = a
        result   = max(self.numPyArray2)
        maxValue = ev.MaxValue()
        self.assertEquals(maxValue, result)

    def testMeanValue(self):
        "Test Epetra.FEVector MeanValue method"
        a         = self.numPyArray1
        ev        = Epetra.FEVector(self.map)
        ev[:]     = a
        result    = sum(self.numPyArray1)/self.length
        meanValue = ev.MeanValue()
        self.assertEquals(meanValue, result)

    def testSetArray(self):
        "Test Epetra.FEVector __setattr__ method for 'array'"
        ev = Epetra.FEVector(self.map)
        self.assertRaises(AttributeError, ev.__setattr__, "array", "junk")

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraVectorTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n***********************\nTesting Epetra.FEVector\n***********************\n"
    verbosity = options.verbosity * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
