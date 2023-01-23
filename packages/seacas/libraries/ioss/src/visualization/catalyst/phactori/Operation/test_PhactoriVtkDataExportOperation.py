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
import vtk
from Operation.PhactoriVtkDataExportOperation import *
import paraview.simple
import subprocess
import os

class TestPhactoriVtkDataExportOperation(unittest.TestCase):

  def test_SimpleImageDataVtkMultiblockOutput(self):
    testWavelet = Wavelet()
    testWavelet.UpdatePipeline()
    #testGroup = GroupDatasets(input=[testWavelet, wavelet2])
    testGroup = GroupDatasets(registrationName='GroupDatasets1', Input=testWavelet)
    #testGroup = GroupDatasets()
    #print(dir(testGroup))
    testGroup.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_vtkdataexport_simpleimagedata"
    operationParams = {
      "type":"vtkdataexport",
      "input":"testgroup1",
      "basedirectory":".",
      "basename":"test_vtkdataexport_simpleimagedata"
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'slicewithplane',
              PhactoriVtkDataExportOperation,
              operationParams)
    newOperationBlock.mOperationSpecifics.myCopyOfInputFilter = testGroup
    newOperationBlock.mOperationSpecifics.mFilterToWriteDataFrom = testGroup
    #newOperationBlock.GetPvFilter().updatePipeline()
    newOperationBlock.ExportOperationData(None)

    vtmGeneratedFlag = os.path.exists("test_vtkdataexport_simpleimagedata_0000.vtm")
    self.assertTrue(vtmGeneratedFlag)
    subprocess.run(["rm", "-rf",
      "test_vtkdataexport_simpleimagedata_0000.vtm",
      "test_vtkdataexport_simpleimagedata_0000"])

  def localSetThresholdRange(self, threshOp, lowerValue, upperValue):
    global gParaViewCatalystVersionFlag
    if gParaViewCatalystVersionFlag < 51000:
      threshOp.ThresholdRange = [lowerValue, upperValue]
    else:
      threshOp.LowerThreshold = lowerValue
      threshOp.UpperThreshold = upperValue
      threshOp.ThresholdMethod = vtk.vtkThreshold.THRESHOLD_BETWEEN

  def test_ThresholdUnstructuredVtkMultiblockOutput(self):
    testWavelet = Wavelet()
    testWavelet.UpdatePipeline()
    testThresh1 = Threshold(registrationName='Threshold3', Input=testWavelet)
    testThresh1.Scalars = ['POINTS', 'RTData']
    self.localSetThresholdRange(testThresh1, 150.0, 300.0)
    #testGroup = GroupDatasets(input=[testWavelet, wavelet2])
    testGroup = GroupDatasets(registrationName='GroupDatasets4', Input=testThresh1)
    #testGroup = GroupDatasets()
    #print(dir(testGroup))
    testGroup.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_vtkdataexport_thresholdunstructured"
    operationParams = {
      "type":"vtkdataexport",
      "input":"testgroup1",
      "basedirectory":".",
      "basename":"test_vtkdataexport_thresholdunstructured"
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'slicewithplane',
              PhactoriVtkDataExportOperation,
              operationParams)
    newOperationBlock.mOperationSpecifics.myCopyOfInputFilter = testGroup
    newOperationBlock.mOperationSpecifics.mFilterToWriteDataFrom = testGroup
    #newOperationBlock.GetPvFilter().updatePipeline()
    newOperationBlock.ExportOperationData(None)

    vtmGeneratedFlag = os.path.exists("test_vtkdataexport_thresholdunstructured_0000.vtm")
    self.assertTrue(vtmGeneratedFlag)
    subprocess.run(["rm", "-rf",
      "test_vtkdataexport_thresholdunstructured_0000.vtm",
      "test_vtkdataexport_thresholdunstructured_0000"])

  def test_EmptyThresholdVtkMultiblockOutput(self):
    testWavelet = Wavelet()
    testWavelet.UpdatePipeline()
    testThresh1 = Threshold(registrationName='Threshold3', Input=testWavelet)
    testThresh1.Scalars = ['POINTS', 'RTData']
    self.localSetThresholdRange(testThresh1, 300.0, 400.0)
    #testGroup = GroupDatasets(input=[testWavelet, wavelet2])
    testGroup = GroupDatasets(registrationName='GroupDatasets3', Input=testThresh1)
    #testGroup = GroupDatasets()
    #print(dir(testGroup))
    testGroup.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "test_vtkdataexport_emptythreshold"
    operationParams = {
      "type":"vtkdataexport",
      "input":"testgroup1",
      "basedirectory":".",
      "basename":"test_vtkdataexport_emptythreshold"
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'slicewithplane',
              PhactoriVtkDataExportOperation,
              operationParams)
    newOperationBlock.mOperationSpecifics.myCopyOfInputFilter = testGroup
    newOperationBlock.mOperationSpecifics.mFilterToWriteDataFrom = testGroup
    #newOperationBlock.GetPvFilter().updatePipeline()
    newOperationBlock.ExportOperationData(None)

    vtmGeneratedFlag = os.path.exists("test_vtkdataexport_emptythreshold_0000.vtm")
    self.assertTrue(vtmGeneratedFlag)
    subprocess.run(["rm", "-rf",
      "test_vtkdataexport_emptythreshold_0000.vtm",
      "test_vtkdataexport_emptythreshold_0000"])

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


