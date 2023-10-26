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
from Operation.PhactoriPointDataToCellDataOperation import *
import paraview.simple

class TestPhactoriPointDataToCellDataOperation(unittest.TestCase):

  def localSetThresholdRange(self, threshOp, lowerValue, upperValue):
    global gParaViewCatalystVersionFlag
    if gParaViewCatalystVersionFlag < 51000:
      threshOp.ThresholdRange = [lowerValue, upperValue]
    else:
      threshOp.LowerThreshold = lowerValue
      threshOp.UpperThreshold = upperValue
      threshOp.ThresholdMethod = vtk.vtkThreshold.THRESHOLD_BETWEEN

  def makeTestGroupWithTwoThresholdsFromWavelet(self):
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
    return testGroup

  def test_PointDataToCellData_Multiblock1(self):
    if PhactoriDbg(100):
      myDebugPrint3("test_PointDataToCellData1 entered\n")
    testGroup = self.makeTestGroupWithTwoThresholdsFromWavelet()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_celltopoint"
    operationParams = {
      "input": "testgroup1",
      "type": "node data to element data"
    }
    if PhactoriDbg(100):
      myDebugPrint3("about to parse and construct, testGroup: " + str(testGroup) +"\n")
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'node data to element data',
              PhactoriPointDataToCellDataOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, testGroup)
    newOperationBlock.GetPvFilter().UpdatePipeline()

    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    blk0 = csdata.GetBlock(0)
    cellArrays = blk0.GetCellData()
    self.assertEqual(1, cellArrays.GetNumberOfArrays())
    ptArrays = blk0.GetPointData()
    self.assertEqual(0, ptArrays.GetNumberOfArrays())
    self.assertEqual("RTData", cellArrays.GetArrayName(0))
    celldata = cellArrays.GetArray("RTData")
    self.assertNotEqual(celldata, None)

  def test_PointDataToCellData_Multiblock2PassCellData(self):
    if PhactoriDbg(100):
      myDebugPrint3("test_PointDataToCellData1 entered\n")
    testGroupWithCellData = self.makeTestGroupWithTwoThresholdsFromWavelet()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_celltopoint"
    operationParams = {
      "input": "testgroup1",
      "type": "node data to element data",
      "pass point data": 1
    }
    if PhactoriDbg(100):
      myDebugPrint3("about to parse and construct, testGroup: " + str(testGroup) +"\n")
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'node data to element data',
              PhactoriPointDataToCellDataOperation,
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

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


