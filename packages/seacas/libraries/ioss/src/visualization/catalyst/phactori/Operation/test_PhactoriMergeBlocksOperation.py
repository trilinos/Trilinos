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
from Operation.PhactoriMergeBlocksOperation import *
import paraview.simple

class TestPhactoriMergeBlocksOperation(unittest.TestCase):

  def localSetThresholdRange(self, threshOp, lowerValue, upperValue):
    global gParaViewCatalystVersionFlag
    if gParaViewCatalystVersionFlag < 51000:
      threshOp.ThresholdRange = [lowerValue, upperValue]
    else:
      threshOp.LowerThreshold = lowerValue
      threshOp.UpperThreshold = upperValue
      threshOp.ThresholdMethod = vtk.vtkThreshold.THRESHOLD_BETWEEN

  def MakeTestGroup(self, threshrange1, threshrange2):
    testWavelet = Wavelet()
    testWavelet.UpdatePipeline()

    testThresh1 = Threshold(registrationName='Threshold3', Input=testWavelet)
    testThresh1.Scalars = ['POINTS', 'RTData']
    self.localSetThresholdRange(testThresh1, threshrange1[0], threshrange1[1])
    testThresh2 = Threshold(registrationName='Threshold4', Input=testWavelet)
    testThresh2.Scalars = ['POINTS', 'RTData']
    self.localSetThresholdRange(testThresh2, threshrange2[0], threshrange2[1])
    testGroup = GroupDatasets(registrationName='testgroup1', Input=[testThresh1, testThresh2])
    testGroup.UpdatePipeline()
    return testGroup

  def test_MergeBlocks1(self):
    if PhactoriDbg(100):
      myDebugPrint3("test_MergeBlocks1 entered\n")
    testGroup = self.MakeTestGroup([0.0, 100.0], [200.0, 300.0])
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_mergeblocks"
    operationParams = {
      "input": "testgroup1",
      "type": "mergeblocks"
    }
    if PhactoriDbg(100):
      myDebugPrint3("about to parse and construct, testGroup: " + str(testGroup) +"\n")
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'mergeblocks',
              PhactoriMergeBlocksOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, testGroup)
    if PhactoriDbg(100):
      myDebugPrint3("done to parse and construct, newOperationBlock: " + str(newOperationBlock) +"\n")
    newOperationBlock.GetPvFilter().UpdatePipeline()

    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    icsdClassname = csdata.GetClassName()
    if PhactoriDbg(100):
      myDebugPrint3("icsdClassname: " + str(icsdClassname) + "\n")
    self.assertEqual(icsdClassname, "vtkUnstructuredGrid")

    numPoints = csdata.GetNumberOfPoints()
    if PhactoriDbg(100):
      myDebugPrint3("csdata.GetNumberOfPoints(): " + str(numPoints) + "\n")
    self.assertEqual(numPoints, 1680)

  def test_MergeBlocksWithoutMergingPoints(self):
    if PhactoriDbg(100):
      myDebugPrint3("test_MergeBlocks1 entered\n")
    testGroup = self.MakeTestGroup([50.0, 150.0], [125.0, 200.0])

    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_mergeblocks"
    operationParams = {
      "input": "testgroup1",
      "type": "mergeblocks",
      "merge points": 1
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'mergeblocks',
              PhactoriMergeBlocksOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, testGroup)
    newOperationBlock.GetPvFilter().UpdatePipeline()
    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    icsdClassname = csdata.GetClassName()
    self.assertEqual(icsdClassname, "vtkUnstructuredGrid")
    numPoints = csdata.GetNumberOfPoints()
    self.assertEqual(numPoints, 7491)

    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_mergeblocks"
    operationParams = {
      "input": "testgroup1",
      "type": "mergeblocks",
      "merge points": 0
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'mergeblocks',
              PhactoriMergeBlocksOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, testGroup)
    newOperationBlock.GetPvFilter().UpdatePipeline()
    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    icsdClassname = csdata.GetClassName()
    self.assertEqual(icsdClassname, "vtkUnstructuredGrid")
    numPoints = csdata.GetNumberOfPoints()
    self.assertEqual(numPoints, 8339)

  def test_MergeBlocksMergePartitionsOnly(self):
    if PhactoriDbg(100):
      myDebugPrint3("test_MergeBlocks1 entered\n")
    testGroup = self.MakeTestGroup([50.0, 150.0], [125.0, 200.0])

    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_mergeblocks"
    operationParams = {
      "input": "testgroup1",
      "type": "mergeblocks",
      "merge partitions only": 1
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'mergeblocks',
              PhactoriMergeBlocksOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, testGroup)
    newOperationBlock.GetPvFilter().UpdatePipeline()
    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    icsdClassname = csdata.GetClassName()
    self.assertEqual(icsdClassname, "vtkMultiBlockDataSet")
    self.assertEqual(newOperationBlock.GetPvFilter().MergePartitionsOnly, 1)
    numPoints = csdata.GetNumberOfPoints()
    self.assertEqual(numPoints, 8339)

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


