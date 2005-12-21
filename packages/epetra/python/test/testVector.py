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
        "Test Epetra.Vector (BlockMap,bool) constructor"
        ev = Epetra.Vector(self.map,True)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())

    def testConstructor02(self):
        "Test Epetra.Vector (BlockMap,int) constructor"
        ev = Epetra.Vector(self.map)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())

    def testConstructor03(self):
        "Test Epetra.Vector (BlockMap,bad-list) constructor"
        list = [0, 1.0, "e", "pi"]
        ev = Epetra.Vector(self.map,list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(ev.shape[0]):
            self.assertEquals(ev[i], 0.0)

    def testConstructor04(self):
        "Test Epetra.Vector (BlockMap,1D-small-list) constructor"
        list = [0, 1.0, 2, 3.14]
        ev = Epetra.Vector(self.map,list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            self.assertEquals(ev[i], list[i])

    def testConstructor05(self):
        "Test Epetra.Vector (BlockMap,1D-correct-list) constructor"
        list = [0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0]
        ev = Epetra.Vector(self.map,list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            self.assertEquals(ev[i], list[i])

    def testConstructor06(self):
        "Test Epetra.Vector (BlockMap,1D-big-list) constructor"
        list = [0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0, -1, -2]
        ev = Epetra.Vector(self.map,list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(ev.shape[0]):
            self.assertEquals(ev[i], list[i])

    def testConstructor07(self):
        "Test Epetra.Vector (BlockMap,2D-small-list) constructor"
        list = [[0,    1.0, 2],
                [2.17, 1.5, 1]]
        ev = Epetra.Vector(self.map,list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            for j in range(len(list[i])):
                self.assertEquals(ev[i*3+j], list[i][j])

    def testConstructor08(self):
        "Test Epetra.Vector (BlockMap,2D-big-list) constructor"
        list = [[0  , 1.0, 2  , 3.14, 2.17, 1.5],
                [1.5, 1  , 0.5, 0.0 , -1  , -2 ]]
        ev = Epetra.Vector(self.map,list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            for j in range(len(list[i])):
                ii = i * 6 + j
                if (ii < ev.MyLength()):
                    self.assertEquals(ev[ii], list[i][j])

    def testConstructor09(self):
        "Test Epetra.Vector (BlockMap,3D-small-list) constructor"
        list = [[[0   , 1.0], [2, 3.14]],
                [[2.17, 1.5], [1, 0.5 ]]]
        ev = Epetra.Vector(self.map,list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            for j in range(len(list[i])):
                for k in range(len(list[i][j])):
                    self.assertEquals(ev[4*i+2*j+k], list[i][j][k])

    def testConstructor10(self):
        "Test Epetra.Vector (BlockMap,3D-correct-list) constructor"
        list = [[[0  , 1.0, 2], [3.14, 2.17, 1.5 ], [1, 0.5, 0.0]]]
        ev = Epetra.Vector(self.map,list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            for j in range(len(list[i])):
                for k in range(len(list[i][j])):
                    self.assertEquals(ev[i,j,k], list[i][j][k])

    def testConstructor11(self):
        "Test Epetra.Vector (BlockMap,3D-big-list) constructor"
        list = [[[  0, 1.0, 2  ], [3.14, 2.17, 1.5 ], [1   , 0.5, 0.0]],
                [[ -1,   0, 1.0], [2   , 3.14, 2.17], [1.5 , 1  , 0.5]],
                [[0.0,  -1, 0  ], [1.0 , 2   , 3.14], [2.17, 1.5, 1  ]]]
        ev = Epetra.Vector(self.map,list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(len(list)):
            for j in range(len(list[i])):
                for k in range(len(list[i][j])):
                    ii = i*9 + j*3 + k
                    if (ii < ev.MyLength()):
                        self.assertEquals(ev[ii], list[i][j][k])

    def testConstructor12(self):
        "Test Epetra.Vector (1D-list) constructor"
        list = [0, 1.0, 2, 3.14, 2.17, 1.5, 1, 0.5, 0.0]
        ev = Epetra.Vector(list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), len(list))
        self.assertEquals(ev.GlobalLength(), len(list))
        for i in range(len(list)):
            self.assertEquals(ev[i], list[i])

    def testConstructor13(self):
        "Test Epetra.Vector (2D-list) constructor"
        list = [[0  , 1.0, 2, 3.14, 2.17, 1.5 , 1, 0.5, 0.0],
                [0.0, 0.5, 1, 1.5 , 2.17, 3.14, 2, 1.0, 0  ]]
        ev = Epetra.Vector(list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), 2*len(list[0]))
        self.assertEquals(ev.GlobalLength(), 2*len(list[0]))
        for i in range(ev.shape[0]):
            for j in range(len(list[i])):
                self.assertEquals(ev[i,j], list[i][j])

    def testConstructor14(self):
        "Test Epetra.Vector (3D-list) constructor"
        list = [[[0   ,1.0,2  ], [3.14,2.17,1.5 ], [1,0.5,0.0 ]],
                [[0.0 ,0.5,1  ], [1.5 ,2.17,3.14], [2,1.0,0   ]],
                [[3.14,2  ,1.0], [0   ,0.0 ,0.5 ], [1,1.5,2.17]]]
        ev = Epetra.Vector(list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(), 3*3*len(list[0][0]))
        self.assertEquals(ev.GlobalLength(), 3*3*len(list[0][0]))
        for i in range(ev.shape[0]):
            for j in range(len(list[i])):
                self.assertEquals(ev[i,j], list[i][j])

    def testConstructor15(self):
        "Test Epetra.Vector (bad-list) constructor"
        list = [0, 1.0, "e", "pi"]
        ev  = Epetra.Vector(list)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(),len(list))
        self.assertEquals(ev.GlobalLength(), len(list))
        for i in range(ev.shape[0]):
            self.assertEquals(ev[i], 0.0)

    def testConstructor16(self):
        "Test Epetra.Vector (Copy,MultiVector,int) constructor"
        list = [[0, 1, 2, 3, 4, 5, 5, 4, 3],
                [2, 1, 0, 0, 1, 2, 3, 4, 5]]
        emv  = Epetra.MultiVector(self.map,list)
        ev   = Epetra.Vector(Epetra.Copy,emv,1)
        self.assertEquals(ev.NumVectors(),1)
        self.assertEquals(ev.MyLength(),self.length)
        self.assertEquals(ev.GlobalLength(), self.length*comm.NumProc())
        for i in range(ev.shape[0]):
            self.assertEquals(ev[i], emv[1,i])

    def testConstructor17(self):
        "Test Epetra.Vector copy constructor"
        ev1 = Epetra.Vector(self.map)
        ev2 = Epetra.Vector(ev1)
        self.assertEquals(ev2.NumVectors(),   ev1.NumVectors()  )
        self.assertEquals(ev2.MyLength(),     ev1.MyLength()    )
        self.assertEquals(ev2.GlobalLength(), ev1.GlobalLength())
        for i in range(len(ev1)):
            self.assertEquals(ev1[i], ev2[i])

    def testReplaceMap1(self):
        "Test Epetra.Vector ReplaceMap method with good map"
        blockMap = Epetra.BlockMap(3*self.comm.NumProc(),3,0,self.comm)
        ev = Epetra.Vector(self.map)
        result = ev.ReplaceMap(blockMap)
        self.assertEquals(result, 0)
        newMap = ev.Map()
        self.assertEquals(newMap.ElementSize(), blockMap.ElementSize())

    def testReplaceMap2(self):
        "Test Epetra.Vector ReplaceMap method with bad map"
        blockMap = Epetra.BlockMap(2*self.comm.NumProc(),5,0,self.comm)
        ev = Epetra.Vector(self.map)
        result = ev.ReplaceMap(blockMap)
        self.assertEquals(result, -1)
        newMap = ev.Map()
        self.assertEquals(newMap.ElementSize(), self.map.ElementSize())

#     def testReplaceGlobalValue1(self):
#         "Test Epetra.Vector ReplaceGlobalValue method"
#         ev = Epetra.Vector(self.map,self.numPyArray1)
#         gid = 4
#         lid = self.map.LID(gid)
#         self.assertEquals(ev[gid], 0.5)
#         result = ev.ReplaceGlobalValue(gid,5.0)
#         if lid >= 0:
#             self.assertEquals(result, 0)
#             self.assertEquals(ev[lid], 5.0)
#         else:
#             self.assertEquals(result, 1)

#     def testReplaceGlobalValue2(self):
#         "Test Epetra.Vector ReplaceGlobalValue method for BlockMaps"
#         map = Epetra.BlockMap(3*self.comm.NumProc(),3,0,self.comm)
#         self.numPyArray1.shape = (1,3,3)  # 1 vector, 3 elements, 3 points per element
#         ev = Epetra.Vector(map,self.numPyArray1)
#         gid = 1
#         lid = self.map.LID(gid)
#         self.assertEquals(ev[0,gid,1], 0.5)
#         result = ev.ReplaceGlobalValue(gid,1,0,5.0)
#         if lid >= 0:
#             self.assertEquals(result, 0)
#             self.assertEquals(ev[0,lid,1], 5.0)
#         else:
#             self.assertEquals(result, 1)

#     def testSumIntoGlobalValue1(self):
#         "Test Epetra.Vector SumIntoGlobalValue method"
#         ev = Epetra.Vector(self.map,self.numPyArray1)
#         gid = 4
#         lid = self.map.LID(gid)
#         self.assertEquals(ev[0,gid], 0.5)
#         result = ev.SumIntoGlobalValue(gid,0,0.5)
#         if lid >= 0:
#             self.assertEquals(result, 0)
#             self.assertEquals(ev[0,lid], 1.0)
#         else:
#             self.assertEquals(result, 1)

#     def testSumIntoGlobalValue2(self):
#         "Test Epetra.Vector SumIntoGlobalValue method for BlockMaps"
#         map = Epetra.BlockMap(3*self.comm.NumProc(),3,0,self.comm)
#         self.numPyArray1.shape = (1,3,3)  # 1 vector, 3 elements, 3 points per element
#         ev = Epetra.Vector(map,self.numPyArray1)
#         gid = 1
#         lid = self.map.LID(gid)
#         self.assertEquals(ev[0,gid,1], 0.5)
#         result = ev.SumIntoGlobalValue(gid,1,0,0.5)
#         if lid >= 0:
#             self.assertEquals(result, 0)
#             self.assertEquals(ev[0,lid,1], 1.0)
#         else:
#             self.assertEquals(result, 1)

#     def testReplaceMyValue1(self):
#         "Test Epetra.Vector ReplaceMyValue method"
#         ev = Epetra.Vector(self.map,self.numPyArray1)
#         lid = 4
#         self.assertEquals(ev[0,lid], 0.5)
#         result = ev.ReplaceMyValue(lid,0,5.0)
#         self.assertEquals(result, 0)
#         self.assertEquals(ev[0,lid], 5.0)

#     def testReplaceMyValue2(self):
#         "Test Epetra.Vector ReplaceMyValue method for BlockMaps"
#         map = Epetra.BlockMap(3*self.comm.NumProc(),3,0,self.comm)
#         self.numPyArray1.shape = (1,3,3)  # 1 vector, 3 elements, 3 points per element
#         ev = Epetra.Vector(map,self.numPyArray1)
#         lid = 1
#         self.assertEquals(ev[0,lid,1], 0.5)
#         result = ev.ReplaceMyValue(lid,1,0,5.0)
#         self.assertEquals(result, 0)
#         self.assertEquals(ev[0,lid,1], 5.0)

#     def testSumIntoMyValue1(self):
#         "Test Epetra.Vector SumIntoMyValue method"
#         ev = Epetra.Vector(self.map,self.numPyArray1)
#         lid = 4
#         self.assertEquals(ev[0,lid], 0.5)
#         result = ev.SumIntoMyValue(lid,0,0.5)
#         self.assertEquals(result, 0)
#         self.assertEquals(ev[0,lid], 1.0)

#     def testSumIntoMyValue2(self):
#         "Test Epetra.Vector SumIntoMyValue method for BlockMaps"
#         map = Epetra.BlockMap(3*self.comm.NumProc(),3,0,self.comm)
#         self.numPyArray1.shape = (1,3,3)  # 1 vector, 3 elements, 3 points per element
#         ev = Epetra.Vector(map,self.numPyArray1)
#         lid = 1
#         self.assertEquals(ev[0,lid,1], 0.5)
#         result = ev.SumIntoMyValue(lid,1,0,0.5)
#         self.assertEquals(result, 0)
#         self.assertEquals(ev[0,lid,1], 1.0)

#     def testPutScalar(self):
#         "Test Epetra.Vector PutScalar method"
#         numVec = 3
#         ev    = Epetra.Vector(self.map,numVec)
#         for v in range(numVec):
#             for i in range(self.map.NumMyPoints()):
#                 self.assertEquals(ev[v,i], 0.0)
#         scalar = 3.14
#         ev.PutScalar(scalar)
#         for v in range(numVec):
#             for i in range(self.map.NumMyPoints()):
#                 self.assertEquals(ev[v,i], scalar)

#     def testRandom(self):
#         "Test Epetra.Vector Random method"
#         numVec = 2
#         ev    = Epetra.Vector(self.map,numVec)
#         scalar = 3.14
#         ev.PutScalar(scalar)
#         for v in range(numVec):
#             for i in range(self.map.NumMyPoints()):
#                 self.assertEquals(ev[v,i], scalar)
#         ev.Random()
#         for v in range(numVec):
#             for i in range(self.map.NumMyPoints()):
#                 self.assertEquals(ev[v,i]>-1.0, True)
#                 self.assertEquals(ev[v,i]< 1.0, True)

#     def testExtractCopy(self):
#         "Test Epetra.Vector ExtractCopy method"
#         a     = [self.numPyArray1] * 4
#         ev   = Epetra.Vector(self.map,a)
#         array = ev.ExtractCopy()
#         self.assertEquals(type(array), ArrayType)
#         self.assertEquals(ev[:], array[:])
#         self.assertEquals(ev.array is array, False)

#     def testExtractView(self):
#         "Test Epetra.Vector ExtractView method"
#         a     = [self.numPyArray1] * 4
#         ev   = Epetra.Vector(self.map,a)
#         array = ev.ExtractView()
#         self.assertEquals(type(array), ArrayType)
#         self.assertEquals(ev[:], array[:])
#         self.assertEquals(ev.array is array, True)

#     def testNumVectors(self):
#         "Test Epetra.Vector NumVectors method"
#         for i in range(1,8):
#             ev = Epetra.Vector(self.map,i)
#             self.assertEquals(ev.NumVectors(), i)

#     def testMyLength(self):
#         "Test Epetra.Vector MyLength method"
#         a   = [self.numPyArray1] * 3
#         ev = Epetra.Vector(self.map,a)
#         self.assertEquals(ev.MyLength(), self.length)

#     def testGlobalLength(self):
#         "Test Epetra.Vector GlobalLength method"
#         a   = [self.numPyArray1] * 3
#         ev = Epetra.Vector(self.map,a)
#         self.assertEquals(ev.GlobalLength(), self.length*self.comm.NumProc())

#     def testConstantStride(self):
#         "Test Epetra.Vector ConstantStride method"
#         squareArray = [[-1.2,  3.4, -5.6],
#                        [ 7.8, -9.0,  1.2],
#                        [-3.4,  5.6, -7.8]]
#         multiArray = [squareArray] * 8
#         ev1  = Epetra.Vector(self.map,multiArray)
#         indexes = (0,2,3,7)
#         ev2  = Epetra.Vector(Epetra.View,ev1,indexes)
#         self.assertEquals(ev1.ConstantStride(), True )
#         self.assertEquals(ev2.ConstantStride(), False)

#     def testStride(self):
#         "Test Epetra.Vector Stride method"
#         squareArray = [[-1.2,  3.4, -5.6],
#                        [ 7.8, -9.0,  1.2],
#                        [-3.4,  5.6, -7.8]]
#         multiArray = [squareArray] * 8
#         ev  = Epetra.Vector(self.map,multiArray)
#         self.assertEquals(ev.Stride(), 9)

#     def testSeed(self):
#         "Test Epetra.Vector Seed method"
#         ev  = Epetra.Vector(self.map,1)
#         seed = ev.Seed()
#         max  = 2**31 - 1
#         self.assertEquals(seed>0,   True)
#         self.assertEquals(seed<max, True)

#     def testSetSeed1(self):
#         "Test Epetra.Vector SetSeed method"
#         ev    = Epetra.Vector(self.map,1)
#         seed   = 2005
#         result = ev.SetSeed(seed)
#         self.assertEquals(result,     0   )
#         self.assertEquals(ev.Seed(), seed)

#     def testSetSeed2(self):
#         "Test Epetra.Vector SetSeed method with negative seed"
#         ev  = Epetra.Vector(self.map,1)
#         seed = -2005
#         self.assertRaises(TypeError,ev.SetSeed,seed)

#     def testPrint(self):
#         "Test Epetra.Vector Print method"
#         output = ""
#         if self.comm.MyPID() == 0:
#             output += "%10s%14s%20s%20s  \n" % ("MyPID","GID","Value","Value")
#         for lid in range(self.length):
#             gid = self.map.GID(lid)
#             output += "%10d%14d%24d%20d\n" % (self.comm.MyPID(),gid,0,0)
#         ev = Epetra.Vector(self.map,2)
#         filename = "testVector%d.dat" % self.comm.MyPID()
#         f = open(filename,"w")
#         ev.Print(f)
#         f.close()
#         self.assertEqual(open(filename,"r").read(), output)

#     def testDot(self):
#         "Test Epetra.Vector Dot method"
#         map    = Epetra.Map(4*self.comm.NumProc(),0,self.comm)
#         array1 = [[-1, 2,-3, 4],
#                   [ 5, 1,-8,-7]]
#         array2 = [[ 9, 0,-1,-2],
#                   [-7,-8, 1, 5]]
#         ev1   = Epetra.Vector(map,array1)
#         ev2   = Epetra.Vector(map,array2)
#         dot    = ev1.Dot(ev2)
#         result = array([-14,-86])*self.comm.NumProc()
#         self.assertEqual(dot[:], result)

#     def testAbs(self):
#         "Test Epetra.Vector Abs method"
#         a    = array([self.numPyArray1,self.numPyArray2])
#         ev1 = Epetra.Vector(self.map,a)
#         ev2 = Epetra.Vector(self.map,2)
#         self.assertEquals(ev2[:],0.0)
#         result = ev2.Abs(ev1)
#         self.assertEquals(result,0)
#         self.assertEquals(ev2[:],abs(a))

#     def testReciprocal(self):
#         "Test Epetra.Vector Reciprocal method"
#         a    = array([self.numPyArray1,self.numPyArray2])
#         a[0,0] = a[1,0] = 1.0  # Don't want to invert zero
#         ev1 = Epetra.Vector(self.map,a)
#         ev2 = Epetra.Vector(self.map,2)
#         self.assertEquals(ev2[:],0.0)
#         result = ev2.Reciprocal(ev1)
#         self.assertEquals(result,0)
#         self.assertEquals(ev2[:],1.0/a)

#     def testScale1(self):
#         "Test Epetra.Vector Scale method in-place"
#         a   = array([self.numPyArray1,self.numPyArray2])
#         ev = Epetra.Vector(self.map,a)
#         self.assertEquals(ev[:],0.0)
#         result = ev.Scale(2.0)
#         self.assertEquals(result,0)
#         self.assertEquals(ev[:],2.0*a)

#     def testScale2(self):
#         "Test Epetra.Vector Scale method with replace"
#         a    = array([self.numPyArray1,self.numPyArray2])
#         ev1 = Epetra.Vector(self.map,a)
#         ev2 = Epetra.Vector(self.map,2)
#         self.assertEquals(ev2[:],0.0)
#         result = ev2.Scale(pi,ev1)
#         self.assertEquals(result,0)
#         self.assertEquals(ev2[:],pi*a)

#     def testUpdate1(self):
#         "Test Epetra.Vector Update method with one Vector"
#         ev1 = Epetra.Vector(self.map,self.numPyArray1)
#         ev2 = Epetra.Vector(self.map,self.numPyArray2)
#         result = ev2.Update(2.0,ev1,3.0)
#         self.assertEquals(result,0)
#         self.assertEquals(ev2[:],2.0*self.numPyArray1 + 3.0*self.numPyArray2)

#     def testUpdate2(self):
#         "Test Epetra.Vector Update method with two Vectors"
#         ev0 = Epetra.Vector(self.map,1               )
#         ev1 = Epetra.Vector(self.map,self.numPyArray1)
#         ev2 = Epetra.Vector(self.map,self.numPyArray2)
#         result = ev0.Update(2.0,ev1,3.0,ev2,pi)
#         self.assertEquals(result,0)
#         self.assertEquals(ev0[:],2.0*self.numPyArray1 + 3.0*self.numPyArray2)

#     def testNorm1(self):
#         "Test Epetra.Vector Norm1 method"
#         a      = [self.numPyArray1,self.numPyArray2]
#         ev    = Epetra.Vector(self.map,a)
#         result = [sum(self.numPyArray1) * self.comm.NumProc(),
#                   sum(self.numPyArray2) * self.comm.NumProc()]
#         norm1  = ev.Norm1()
#         self.assertEquals(len(norm1),2     )
#         self.assertEquals(norm1[:],  result)

#     def testNorm2(self):
#         "Test Epetra.Vector Norm2 method"
#         a      = [self.numPyArray1,self.numPyArray2]
#         ev    = Epetra.Vector(self.map,a)
#         result = [sqrt(sum(self.numPyArray1*self.numPyArray1) * self.comm.NumProc()),
#                   sqrt(sum(self.numPyArray2*self.numPyArray2) * self.comm.NumProc())]
#         norm2  = ev.Norm2()
#         self.assertEquals(len(norm2),2     )
#         self.assertEquals(norm2[:],  result)

#     def testNormInf(self):
#         "Test Epetra.Vector NormInf method"
#         a       = [self.numPyArray1,self.numPyArray2]
#         ev     = Epetra.Vector(self.map,a)
#         result  = [max(abs(self.numPyArray1)),max(abs(self.numPyArray2))]
#         normInf = ev.NormInf()
#         self.assertEquals(len(normInf),2     )
#         self.assertEquals(normInf[:],  result)

#     def testNormWeighted(self):
#         "Test Epetra.Vector NormWeighted method"
#         a       = array([self.numPyArray1,self.numPyArray2])
#         ev     = Epetra.Vector(self.map,a)
#         wts     = sin(pi*(arange(self.length) + 0.5) / self.length)
#         weights = Epetra.Vector(self.map,wts)
#         result  = sqrt(sum((a/wts)**2,1)/self.length)
#         norm    = ev.NormWeighted(weights)
#         self.assertEquals(len(norm),2     )
#         self.assertEquals(norm[:],  result)

#     def testMinValue(self):
#         "Test Epetra.Vector MinValue method"
#         a        = [self.numPyArray1,self.numPyArray2]
#         ev      = Epetra.Vector(self.map,a)
#         result   = [min(self.numPyArray1),min(self.numPyArray2)]
#         minValue = ev.MinValue()
#         self.assertEquals(len(minValue),2     )
#         self.assertEquals(minValue[:],  result)

#     def testMaxValue(self):
#         "Test Epetra.Vector MaxValue method"
#         a        = [self.numPyArray1,self.numPyArray2]
#         ev      = Epetra.Vector(self.map,a)
#         result   = [max(self.numPyArray1),max(self.numPyArray2)]
#         maxValue = ev.MaxValue()
#         self.assertEquals(len(maxValue),2     )
#         self.assertEquals(maxValue[:],  result)

#     def testMeanValue(self):
#         "Test Epetra.Vector MeanValue method"
#         a         = [self.numPyArray1,self.numPyArray2]
#         ev       = Epetra.Vector(self.map,a)
#         result    = [sum(self.numPyArray1)/self.length,
#                      sum(self.numPyArray2)/self.length]
#         meanValue = ev.MeanValue()
#         self.assertEquals(len(meanValue),2     )
#         self.assertEquals(meanValue[:],  result)

#     def testMultiply1(self):
#         "Test Epetra.Vector Multiply method"
#         n    = 2 * self.comm.NumProc()
#         map  = Epetra.Map(n,0,self.comm)
#         ev0 = Epetra.Vector(map,n)
#         ev1 = Epetra.Vector(map,n)
#         ev2 = Epetra.Vector(map,n)
#         ev0.Random()
#         ev1.Random()
#         ev2.Random()
#         result = ev0.ply(1.0,ev1,ev2,2.0)
#         self.assertEquals(result,0)

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
       "\n*********************\nTesting Epetra.Vector\n*********************\n"
    verbosity = 2 * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))[0]
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
