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
from Operation.PhactoriExtractSurfaceOperation import *
import paraview.simple

class TestPhactoriExtractSurfaceOperation(unittest.TestCase):

  def localSetThresholdRange(self, threshOp, lowerValue, upperValue):
    global gParaViewCatalystVersionFlag
    if gParaViewCatalystVersionFlag < 51000:
      threshOp.ThresholdRange = [lowerValue, upperValue]
    else:
      threshOp.LowerThreshold = lowerValue
      threshOp.UpperThreshold = upperValue
      threshOp.ThresholdMethod = vtk.vtkThreshold.THRESHOLD_BETWEEN

  def test_ExtractSurface1(self):
    if PhactoriDbg(100):
      myDebugPrint3("test_ExtractSurface1 entered\n")
    testWavelet = Wavelet()
    testWavelet.UpdatePipeline()

    testThresh1 = Threshold(registrationName='Threshold3', Input=testWavelet)
    testThresh1.Scalars = ['POINTS', 'RTData']
    self.localSetThresholdRange(testThresh1, 200.0, 300.0)
    #self.localSetThresholdRange(testThresh1, -10000000.0, 1000000000.0)
    testGroup = GroupDatasets(registrationName='testgroup1', Input=[testThresh1])

    #testGroup = GroupDatasets(registrationName='testgroup1', Input=[testWavelet])

    testGroup.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_extractsurface"
    operationParams = {
      "input": "testgroup1",
      "type": "extractsurface"
    }
    if PhactoriDbg(100):
      myDebugPrint3("about to parse and construct, testGroup: " + str(testGroup) +"\n")
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'extractsurface',
              PhactoriExtractSurfaceOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, testGroup)
    if PhactoriDbg(100):
      myDebugPrint3("done to parse and construct, newOperationBlock: " + str(newOperationBlock) +"\n")
    newOperationBlock.GetPvFilter().UpdatePipeline()

    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    icsdClassname = csdata.GetClassName()
    if PhactoriDbg(100):
      myDebugPrint3("icsdClassname: " + str(icsdClassname) + "\n")
    self.assertEqual(icsdClassname, "vtkMultiBlockDataSet")

    numBlocks = csdata.GetNumberOfBlocks()
    if PhactoriDbg(100):
      myDebugPrint3("numBlocks: " + str(numBlocks) + "\n")

    blk0 = csdata.GetBlock(0)
    icsdClassname2 = blk0.GetClassName()
    if PhactoriDbg(100):
      myDebugPrint3("blk0 icsdClassname2: " + str(icsdClassname2) + "\n")
    self.assertEqual(icsdClassname2, "vtkPolyData")

    numPoints = blk0.GetNumberOfPoints()
    if PhactoriDbg(100):
      myDebugPrint3("blk0 numPoints: " + str(numPoints) + "\n")
    self.assertEqual(numPoints, 764)

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


