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

class EpetraMultiVectorTestCase(unittest.TestCase):
    "TestCase class for MultiVector objects"

    def setUp(self):
        self.length     = 9
        self.scale      = 1.0 / (self.length-1)
        self.comm       = Epetra.PyComm()
        self.map        = Epetra.Map(self.length*self.comm.NumProc(),0,self.comm)
        self.numPyArray = arange(self.length) * self.scale
        
    def tearDown(self):
        self.comm.Barrier()

    def testConstructor01(self):
        "Test Epetra.MultiVector (BlockMap,int,bool) constructor"
        emv = Epetra.MultiVector(self.map,2,True)
        self.assertEquals(emv.NumVectors(),2)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())

    def testConstructor02(self):
        "Test Epetra.MultiVector (BlockMap,int) constructor"
        emv = Epetra.MultiVector(self.map,3)
        self.assertEquals(emv.NumVectors(),3)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())

    def testConstructor03(self):
        "Test Epetra.MultiVector (BlockMap,bad-list) constructor"
        list = [0, 1.0, "e", "pi"]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),1)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv.shape[1]):
            self.assertEquals(emv[0,i], 0.0)

    def testConstructor04(self):
        "Test Epetra.MultiVector (BlockMap,1D-small-list) constructor"
        list = [0, 1.0, 2, 3.14]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),1)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            self.assertEquals(emv[0,i], list[i])

    def testConstructor05(self):
        "Test Epetra.MultiVector (BlockMap,1D-correct-list) constructor"
        list = [0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),1)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            self.assertEquals(emv[0,i], list[i])

    def testConstructor06(self):
        "Test Epetra.MultiVector (BlockMap,1D-big-list) constructor"
        list = [0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0, -1, -2]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),1)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv.shape[1]):
            self.assertEquals(emv[0,i], list[i])

    def testConstructor07(self):
        "Test Epetra.MultiVector (BlockMap,2D-small-list) constructor"
        list = [[0,    1.0, 2, 3.14],
                [2.17, 1.5, 1, 0.5 ]]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),2)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv.shape[0]):
            for j in range(len(list[i])):
                self.assertEquals(emv[i,j], list[i][j])

    def testConstructor08(self):
        "Test Epetra.MultiVector (BlockMap,2D-correct-list) constructor"
        list = [[0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0],
                [0.0, 0.5, 1, 1.5, 2.17, 3.14, 2, 1.0, 0]]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),2)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv.shape[0]):
            for j in range(len(list[i])):
                self.assertEquals(emv[i,j], list[i][j])

    def testConstructor09(self):
        "Test Epetra.MultiVector (BlockMap,2D-big-list) constructor"
        list = [[0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0, -1, -2],
                [0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0, -1, -2]]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),2)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv.shape[0]):
            for j in range(emv.shape[1]):
                self.assertEquals(emv[i,j], list[i][j])

    def testConstructor10(self):
        "Test Epetra.MultiVector (BlockMap,3D-small-list) constructor"
        list = [[[0   , 1.0], [2, 3.14]],
                [[2.17, 1.5], [1, 0.5 ]]]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),2)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            for j in range(len(list[i])):
                for k in range(len(list[i][j])):
                    self.assertEquals(emv[i,j*2+k], list[i][j][k])

    def testConstructor11(self):
        "Test Epetra.MultiVector (BlockMap,3D-correct-list) constructor"
        list = [[[0  , 1.0, 2], [3.14, 2.17, 1.5 ], [1, 0.5, 0.0]],
                [[0.0, 0.5, 1], [1.5 , 2.17, 3.14], [2, 1.0, 0  ]]]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),2)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            for j in range(len(list[i])):
                for k in range(len(list[i][j])):
                    self.assertEquals(emv[i,j,k], list[i][j][k])

    def testConstructor12(self):
        "Test Epetra.MultiVector (BlockMap,3D-big-list) constructor"
        list = [[[0, 1.0, 2, 3.14, 2.17], [1.5, 1, 0.5, 0.0, -1]],
                [[0, 1.0, 2, 3.14, 2.17], [1.5, 1, 0.5, 0.0, -1]],
                [[0, 1.0, 2, 3.14, 2.17], [1.5, 1, 0.5, 0.0, -1]]]
        emv = Epetra.MultiVector(self.map,list)
        self.assertEquals(emv.NumVectors(),3)
        self.assertEquals(emv.MyLength(), self.length)
        self.assertEquals(emv.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv.shape[0]):
            for j in range(len(list[i])):
                self.assertEquals(emv[i,j], list[i][0][j])

    def testConstructor13(self):
        "Test Epetra.MultiVector (1D-list) constructor"
        list = [0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0]
        emv = Epetra.MultiVector(list)
        self.assertEquals(emv.NumVectors(),1)
        self.assertEquals(emv.MyLength(), len(list))
        self.assertEquals(emv.GlobalLength(), len(list))
        for i in range(len(list)):
            self.assertEquals(emv[0,i], list[i])

    def testConstructor14(self):
        "Test Epetra.MultiVector (2D-list) constructor"
        list = [[0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0],
                [0.0, 0.5, 1, 1.5, 2.17, 3.14, 2, 1.0, 0]]
        emv = Epetra.MultiVector(list)
        self.assertEquals(emv.NumVectors(),2)
        self.assertEquals(emv.MyLength(), 2*len(list[0]))
        self.assertEquals(emv.GlobalLength(), 2*len(list[0]))
        for i in range(emv.shape[0]):
            for j in range(len(list[i])):
                self.assertEquals(emv[i,j], list[i][j])

    def testConstructor15(self):
        "Test Epetra.MultiVector (3D-list) constructor"
        list = [[[0   ,1.0,2  ], [3.14,2.17,1.5 ], [1,0.5,0.0 ]],
                [[0.0 ,0.5,1  ], [1.5 ,2.17,3.14], [2,1.0,0   ]],
                [[3.14,2  ,1.0], [0   ,0.0 ,0.5 ], [1,1.5,2.17]]]
        emv = Epetra.MultiVector(list)
        self.assertEquals(emv.NumVectors(),3)
        self.assertEquals(emv.MyLength(), 3*3*len(list[0][0]))
        self.assertEquals(emv.GlobalLength(), 3*3*len(list[0][0]))
        for i in range(emv.shape[0]):
            for j in range(len(list[i])):
                self.assertEquals(emv[i,j], list[i][j])

    def testConstructor16(self):
        "Test Epetra.MultiVector (bad-list) constructor"
        list = [0, 1.0, "e", "pi"]
        emv  = Epetra.MultiVector(list)
        self.assertEquals(emv.NumVectors(),1)
        self.assertEquals(emv.MyLength(),len(list))
        self.assertEquals(emv.GlobalLength(), len(list))
        for i in range(emv.shape[1]):
            self.assertEquals(emv[0,i], 0.0)

    def testConstructor17(self):
        "Test Epetra.MultiVector (Copy,MultiVector,range-of-1) constructor"
        list = [[0, 1, 2, 3, 4, 5, 5, 4, 3],
                [2, 1, 0, 0, 1, 2, 3, 4, 5]]
        emv1  = Epetra.MultiVector(self.map,list)
        emv2  = Epetra.MultiVector(Epetra.Copy,emv1,(1,))
        self.assertEquals(emv2.NumVectors(),1)
        self.assertEquals(emv2.MyLength(),self.length)
        self.assertEquals(emv2.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv2.shape[1]):
            self.assertEquals(emv2[0,i], emv1[1,i])

    def testConstructor18(self):
        "Test Epetra.MultiVector (Copy,MultiVector,range-of-4) constructor"
        squareArray = [[-1.2,  3.4, -5.6],
                       [ 7.8, -9.0,  1.2],
                       [-3.4,  5.6, -7.8]]
        multiArray = [squareArray] * 8
        emv1  = Epetra.MultiVector(self.map,multiArray)
        emv2  = Epetra.MultiVector(Epetra.Copy,emv1,range(4))
        self.assertEquals(emv2.NumVectors(),4)
        self.assertEquals(emv2.MyLength(),self.length)
        self.assertEquals(emv2.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv2.shape[0]):
            self.assertEquals(emv2[i], emv1[i])

    def testConstructor19(self):
        "Test Epetra.MultiVector (Copy,MultiVector,discontinuous-range) constructor"
        squareArray = [[-1.2,  3.4, -5.6],
                       [ 7.8, -9.0,  1.2],
                       [-3.4,  5.6, -7.8]]
        multiArray = [squareArray] * 8
        emv1  = Epetra.MultiVector(self.map,multiArray)
        indexes = (0,2,3,7)
        emv2  = Epetra.MultiVector(Epetra.Copy,emv1,indexes)
        self.assertEquals(emv2.NumVectors(),4)
        self.assertEquals(emv2.MyLength(),self.length)
        self.assertEquals(emv2.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv2.shape[0]):
            self.assertEquals(emv2[i], emv1[indexes[i]])

    def testConstructor20(self):
        "Test Epetra.MultiVector (Copy,MultiVector,bad-list) constructor"
        squareArray = [[-1.2,  3.4, -5.6],
                       [ 7.8, -9.0,  1.2],
                       [-3.4,  5.6, -7.8]]
        multiArray = [squareArray] * 8
        emv1  = Epetra.MultiVector(self.map,multiArray)
        indexes = [0, 1.0, "e", "pi"]
        emv2  = Epetra.MultiVector(Epetra.Copy,emv1,indexes)
        # This is to compensate for a bug I do not understand
        shape = list(emv1.shape)
        shape[0] = emv2.shape[0]
        emv2.shape = tuple(shape)
        # Done with the shape-shifting
        self.assertEquals(emv2.NumVectors(),len(indexes))
        self.assertEquals(emv2.MyLength(),self.length)
        self.assertEquals(emv2.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv2.shape[0]):
            self.assertEquals(emv2[i], emv1[i])

    def testConstructor21(self):
        "Test Epetra.MultiVector (View,MultiVector,range-of-1) constructor"
        list = [[0, 1, 2, 3, 4, 5, 5, 4, 3],
                [2, 1, 0, 0, 1, 2, 3, 4, 5]]
        emv1  = Epetra.MultiVector(self.map,list)
        emv2  = Epetra.MultiVector(Epetra.View,emv1,(1,))
        self.assertEquals(emv2.NumVectors(),1)
        self.assertEquals(emv2.MyLength(),self.length)
        self.assertEquals(emv2.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv2.shape[1]):
            self.assertEquals(emv2[0,i], emv1[1,i])

    def testConstructor22(self):
        "Test Epetra.MultiVector (View,MultiVector,range-of-4) constructor"
        squareArray = [[-1.2,  3.4, -5.6],
                       [ 7.8, -9.0,  1.2],
                       [-3.4,  5.6, -7.8]]
        multiArray = [squareArray] * 8
        emv1  = Epetra.MultiVector(self.map,multiArray)
        emv2  = Epetra.MultiVector(Epetra.View,emv1,range(4))
        self.assertEquals(emv2.NumVectors(),4)
        self.assertEquals(emv2.MyLength(),self.length)
        self.assertEquals(emv2.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv2.shape[0]):
            self.assertEquals(emv2[i], emv1[i])

    def testConstructor23(self):
        "Test Epetra.MultiVector (View,MultiVector,discontinuous-range) constructor"
        squareArray = [[-1.2,  3.4, -5.6],
                       [ 7.8, -9.0,  1.2],
                       [-3.4,  5.6, -7.8]]
        multiArray = [squareArray] * 8
        emv1  = Epetra.MultiVector(self.map,multiArray)
        indexes = (0,2,3,7)
        emv2  = Epetra.MultiVector(Epetra.View,emv1,indexes)
        self.assertEquals(emv2.NumVectors(),4)
        self.assertEquals(emv2.MyLength(),self.length)
        self.assertEquals(emv2.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv2.shape[0]):
            self.assertEquals(emv2[i], emv1[indexes[i]])

    def testConstructor24(self):
        "Test Epetra.MultiVector (View,MultiVector,bad-list) constructor"
        squareArray = [[-1.2,  3.4, -5.6],
                       [ 7.8, -9.0,  1.2],
                       [-3.4,  5.6, -7.8]]
        multiArray = [squareArray] * 8
        emv1  = Epetra.MultiVector(self.map,multiArray)
        indexes = [0, 1.0, "e", "pi"]
        emv2  = Epetra.MultiVector(Epetra.View,emv1,indexes)
        # This is to compensate for a bug I do not understand
        shape = list(emv1.shape)
        shape[0] = emv2.shape[0]
        emv2.shape = tuple(shape)
        # Done with the shape-shifting
        self.assertEquals(emv2.NumVectors(),len(indexes))
        self.assertEquals(emv2.MyLength(),self.length)
        self.assertEquals(emv2.GlobalLength(), self.length*comm.NumProc())
        for i in range(emv2.shape[0]):
            self.assertEquals(emv2[i], emv1[i])

    def testConstructor25(self):
        "Test Epetra.MultiVector copy constructor"
        emv1 = Epetra.MultiVector(self.map,4)
        emv2 = Epetra.MultiVector(emv1)
        self.assertEquals(emv2.NumVectors(),   emv1.NumVectors()  )
        self.assertEquals(emv2.MyLength(),     emv1.MyLength()    )
        self.assertEquals(emv2.GlobalLength(), emv1.GlobalLength())

    def testReplaceMap1(self):
        "Test Epetra.MultiVector ReplaceMap method with good map"
        blockMap = Epetra.BlockMap(3*self.comm.NumProc(),3,0,self.comm)
        emv = Epetra.MultiVector(self.map,3)
        result = emv.ReplaceMap(blockMap)
        self.assertEquals(result, 0)
        newMap = emv.Map()
        self.assertEquals(newMap.ElementSize(), blockMap.ElementSize())

    def testReplaceMap2(self):
        "Test Epetra.MultiVector ReplaceMap method with bad map"
        blockMap = Epetra.BlockMap(2*self.comm.NumProc(),5,0,self.comm)
        emv = Epetra.MultiVector(self.map,4)
        result = emv.ReplaceMap(blockMap)
        self.assertEquals(result, -1)
        newMap = emv.Map()
        self.assertEquals(newMap.ElementSize(), self.map.ElementSize())

    def testReplaceGlobalValue1(self):
        "Test Epetra.MultiVector ReplaceGlobalValue method"
        emv = Epetra.MultiVector(self.map,self.numPyArray)
        gid = 4
        lid = self.map.LID(gid)
        print "pid = %d, lid = %d" % (self.comm.MyPID(), lid)
        self.assertEquals(emv[0,gid], 0.5)
        print emv
        emv.ReplaceGlobalValue(gid,0,5.0)
        print emv
        if lid >= 0:
            self.assertEquals(emv[0,lid], 5.0)

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraMultiVectorTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n**************************\nTesting Epetra.MultiVector\n**************************\n"
    verbosity = 2 * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))[0]
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
