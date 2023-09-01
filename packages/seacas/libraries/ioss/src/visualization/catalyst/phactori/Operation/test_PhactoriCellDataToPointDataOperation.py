# Copyright(C) 1999-2022 National Technology & Engineering Solutions
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

from phactori import *
import unittest
#from vtk import *
import vtk
from Operation.PhactoriCellDataToPointDataOperation import *
import paraview.simple

class TestPhactoriCellDataToPointDataOperation(unittest.TestCase):

  def localSetThresholdRange(self, threshOp, lowerValue, upperValue):
    global gParaViewCatalystVersionFlag
    if gParaViewCatalystVersionFlag < 51000:
      threshOp.ThresholdRange = [lowerValue, upperValue]
    else:
      threshOp.LowerThreshold = lowerValue
      threshOp.UpperThreshold = upperValue
      threshOp.ThresholdMethod = vtk.vtkThreshold.THRESHOLD_BETWEEN

  def makeTestGroupWithCellArrayWithTwoThresholdsFromWavelet(self):
    testWavelet = Wavelet()
    testWavelet.UpdatePipeline()

    testThresh1 = Threshold(registrationName='Threshold3', Input=testWavelet)
    testThresh1.Scalars = ['POINTS', 'RTData']
    self.localSetThresholdRange(testThresh1, 0.0, 100.0)
    testThresh2 = Threshold(registrationName='Threshold4', Input=testWavelet)
    testThresh2.Scalars = ['POINTS', 'RTData']
    self.localSetThresholdRange(testThresh2, 200.0, 300.0)

    testGroup = GroupDatasets(registrationName='testgroup1', Input=[testThresh1, testThresh2])

    testGroup.UpdatePipeline()
    testGroupWithCellData = PointDatatoCellData(registrationName='testGroupWithCellData', Input=testGroup)
    testGroupWithCellData.UpdatePipeline()
    return testGroupWithCellData

  def test_CellDataToPointData_Multiblock1(self):
    if PhactoriDbg(100):
      myDebugPrint3("test_CellDataToPointData1 entered\n")
    testGroupWithCellData = self.makeTestGroupWithCellArrayWithTwoThresholdsFromWavelet()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_celltopoint"
    operationParams = {
      "input": "testgroup1",
      "type": "element data to node data"
    }
    if PhactoriDbg(100):
      myDebugPrint3("about to parse and construct, testGroup: " + str(testGroup) +"\n")
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'element data to node data',
              PhactoriCellDataToPointDataOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, testGroupWithCellData)
    newOperationBlock.GetPvFilter().UpdatePipeline()

    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    blk0 = csdata.GetBlock(0)
    cellArrays = blk0.GetCellData()
    self.assertEqual(0, cellArrays.GetNumberOfArrays())
    ptArrays = blk0.GetPointData()
    self.assertEqual(1, ptArrays.GetNumberOfArrays())
    self.assertEqual("RTData", ptArrays.GetArrayName(0))
    ptdata = ptArrays.GetArray("RTData")
    self.assertNotEqual(ptdata, None)
    self.assertEqual(newOperationBlock.GetPvFilter().PieceInvariant, 0)

  def test_CellDataToPointData_Multiblock2PassPointData(self):
    if PhactoriDbg(100):
      myDebugPrint3("test_CellDataToPointData1 entered\n")
    testGroupWithCellData = self.makeTestGroupWithCellArrayWithTwoThresholdsFromWavelet()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_celltopoint"
    operationParams = {
      "input": "testgroup1",
      "type": "element data to node data",
      "pass cell data": 1,
      "piece invariant": 1
    }
    if PhactoriDbg(100):
      myDebugPrint3("about to parse and construct, testGroup: " + str(testGroup) +"\n")
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'element data to node data',
              PhactoriCellDataToPointDataOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, testGroupWithCellData)
    newOperationBlock.GetPvFilter().UpdatePipeline()

    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    blk0 = csdata.GetBlock(0)
    cellArrays = blk0.GetCellData()
    self.assertEqual(1, cellArrays.GetNumberOfArrays())
    ptArrays = blk0.GetPointData()
    self.assertEqual(1, ptArrays.GetNumberOfArrays())
    self.assertEqual("RTData", ptArrays.GetArrayName(0))
    ptdata = ptArrays.GetArray("RTData")
    self.assertNotEqual(ptdata, None)
    self.assertEqual(newOperationBlock.GetPvFilter().PieceInvariant, 1)

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


