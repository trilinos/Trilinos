# Copyright(C) 1999-2020 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of NTESS nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import unittest
#from vtk import *
#import vtk
from Operation.PhactoriMarkCellSurfaceStatus2 import *
import paraview.simple

class TestPhactoriMarkCellSurfaceStatus2(unittest.TestCase):

  def test_MarkFlatSurfaceCellsFromFlattendWavelet1(self):
    wv1 = paraview.simple.Wavelet()
    xfrm1 = Transform(Input=wv1)
    xfrm1.Transform = 'Transform'
    xfrm1.Transform.Scale = [1.0, 0.005, 1.0]
    xfrm1.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"
    operationParams = {}
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'markcellsurfacestatus2',
              PhactoriMarkCellSurfaceStatus2,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, xfrm1)
    newOperationBlock.GetPvFilter().UpdatePipeline()
    cellData = newOperationBlock.GetPvFilter().CellData
    #ssAry = cellData.GetArray("surfacestatus")
    #myDebugPrint3("numTuples: " + str(ssAry.GetNumberOfTuples()) + "\n")
    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().\
        GetOutputDataObject(0)
    ssAry = csdata.GetCellData().GetArray("surfacestatus")
    #myDebugPrint3("ssAry: " + str(ssAry) + "\n")
    #myDebugPrint3("tuple(0): " + str(ssAry.GetTuple1(0)) + "\n")
    #myDebugPrint3("tuple(7999): " + str(ssAry.GetTuple1(7999)) + "\n")
    #corner cells
    testStatus = ssAry.GetTuple1(0)
    self.assertEqual(testStatus, 13)
    testStatus = ssAry.GetTuple1(19)
    self.assertEqual(testStatus, 13)
    testStatus = ssAry.GetTuple1((20*20*20-1))
    self.assertEqual(testStatus, 13)
    #edges
    testStatus = ssAry.GetTuple1(1)
    self.assertEqual(testStatus, 12)
    testStatus = ssAry.GetTuple1(20)
    self.assertEqual(testStatus, 6)
    testStatus = ssAry.GetTuple1((20*20*20-1)-1)
    self.assertEqual(testStatus, 12)
    testStatus = ssAry.GetTuple1((20*20*20-1)-20)
    self.assertEqual(testStatus, 6)
    #surface cells, non corner, not edge
    testStatus = ssAry.GetTuple1(30)
    self.assertEqual(testStatus, 5)
    testStatus = ssAry.GetTuple1(402)
    self.assertEqual(testStatus, 11)
    testStatus = ssAry.GetTuple1(21)
    self.assertEqual(testStatus, 5)
    testStatus = ssAry.GetTuple1((20*20*20-1)-29)
    self.assertEqual(testStatus, 5)
    testStatus = ssAry.GetTuple1((20*20*20-1)-402)
    self.assertEqual(testStatus, 11)
    #interior
    testStatus = ssAry.GetTuple1(421)
    self.assertEqual(testStatus, 1)
    testStatus = ssAry.GetTuple1(20*20*10+35)
    self.assertEqual(testStatus, 1)
    testStatus = ssAry.GetTuple1((20*20*20-1)-421)
    self.assertEqual(testStatus, 1)

  def test_MarkNoFlatSurfaceCellsFromFlattendWavelet1(self):
    wv1 = paraview.simple.Wavelet()
    xfrm1 = Transform(Input=wv1)
    xfrm1.Transform = 'Transform'
    xfrm1.Transform.Scale = [1.0, 0.5, 1.0]
    xfrm1.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"
    operationParams = {}
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'markcellsurfacestatus2',
              PhactoriMarkCellSurfaceStatus2,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, xfrm1)
    newOperationBlock.GetPvFilter().UpdatePipeline()
    cellData = newOperationBlock.GetPvFilter().CellData
    #ssAry = cellData.GetArray("surfacestatus")
    #myDebugPrint3("numTuples: " + str(ssAry.GetNumberOfTuples()) + "\n")
    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().\
        GetOutputDataObject(0)
    ssAry = csdata.GetCellData().GetArray("surfacestatus")
    #myDebugPrint3("ssAry: " + str(ssAry) + "\n")
    #myDebugPrint3("tuple(0): " + str(ssAry.GetTuple1(0)) + "\n")
    #myDebugPrint3("tuple(7999): " + str(ssAry.GetTuple1(7999)) + "\n")
    #corner cells
    testStatus = ssAry.GetTuple1(0)
    self.assertEqual(testStatus, 10)
    testStatus = ssAry.GetTuple1(19)
    self.assertEqual(testStatus, 10)
    testStatus = ssAry.GetTuple1((20*20*20-1))
    self.assertEqual(testStatus, 10)
    #edges
    testStatus = ssAry.GetTuple1(1)
    self.assertEqual(testStatus, 9)
    testStatus = ssAry.GetTuple1(20)
    self.assertEqual(testStatus, 3)
    testStatus = ssAry.GetTuple1((20*20*20-1)-1)
    self.assertEqual(testStatus, 9)
    testStatus = ssAry.GetTuple1((20*20*20-1)-20)
    self.assertEqual(testStatus, 3)
    #surface cells, non corner, not edge
    testStatus = ssAry.GetTuple1(30)
    self.assertEqual(testStatus, 2)
    testStatus = ssAry.GetTuple1(402)
    self.assertEqual(testStatus, 8)
    testStatus = ssAry.GetTuple1(21)
    self.assertEqual(testStatus, 2)
    testStatus = ssAry.GetTuple1((20*20*20-1)-29)
    self.assertEqual(testStatus, 2)
    testStatus = ssAry.GetTuple1((20*20*20-1)-402)
    self.assertEqual(testStatus, 8)
    #interior
    testStatus = ssAry.GetTuple1(421)
    self.assertEqual(testStatus, 1)
    testStatus = ssAry.GetTuple1(20*20*10+35)
    self.assertEqual(testStatus, 1)
    testStatus = ssAry.GetTuple1((20*20*20-1)-421)
    self.assertEqual(testStatus, 1)
