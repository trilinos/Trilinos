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

from phactori import *
from paraview.simple import *

import vtk

#from .RecursivelyPrintMultiBlockStats import *
#from .PhactoriMpiUtilities import *

#phactori_combine_to_single_python_file_subpiece_begin_1

class PhactoriVtkDataExportOperation(PhactoriOperationSpecifics):
  """passes through it's input to it's output, but also does an
     ExportOperation which is a parallel vtk writer, .vtm or .vtp"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mOutputFileBasename = "myvtk1"
    self.mOutputFileBasedirectory = "CatalystVtkDataOutput"
    self.mFilterToWriteDataFrom = None
    self.mWriter = None
    self.mNumCounterDigits = 4
    #should be "vtkmultiblockdataset" or "vtkpolydata" or "vtkunstructureddata"
    self.mVtkDataType = "vtkmultiblockdataset"
    self.mVtkDataTypeExtention = ".vtm"
    self.mDataArtifactLister = None

  def SetPhactoriDataArtifactMetaDataControl(self, inDataArtifactLister):
    self.mDataArtifactLister = inDataArtifactLister

  def ParseOutputFilenameParametersFromJson(self, inJson):
    self.mOutputFileBasename = getParameterFromBlock(inJson,
      "basename", self.mPhactoriOperationBlockOwner.mName)

    self.mOutputFileBasedirectory = getParameterFromBlock(inJson,
      "basedirectory", self.mOutputFileBasedirectory)
    CreateDirectoryFromProcessZero(self.mOutputFileBasedirectory)

    self.mNumCounterDigits = getParameterFromBlock(inJson,
      'name digit count', self.mNumCounterDigits)

  def ParseParametersFromJson(self, inJson):
    self.ParseOutputFilenameParametersFromJson(inJson)

    param1 = getParameterFromBlock(inJson,
      "vtkdatatype", self.mVtkDataType)
    self.mVtkDataType = param1.lower()
    if self.mVtkDataType == "vtkmultiblockdataset":
      if PhactoriDbg(100):
        myDebugPrint3("self.mVtkDataType set to vtkmultiblockdataset\n")
      self.mVtkDataTypeExtention = ".vtm"
    elif self.mVtkDataType == "vtkpolydata":
      if PhactoriDbg(100):
        myDebugPrint3("self.mVtkDataType set to vtkpolydata\n")
      self.mVtkDataTypeExtention = ".vtp"
    elif self.mVtkDataType == "vtkunstructureddata":
      if PhactoriDbg(100):
        myDebugPrint3("self.mVtkDataType set to vtkunstructureddata\n")
      self.mVtkDataTypeExtention = ".vtu"
    else:
      myDebugPrint3AndException(
          "PhactoriVtkDataExportOperation::ParseParametersFromJson\n"
          "Error: vtkdatatype must be 'vtkmultiblockdataset' \
           or 'vtkpolydata'\n")

  def CreateParaViewFilter(self, inInputFilter):
    self.mFilterToWriteDataFrom = inInputFilter
    return self.mFilterToWriteDataFrom

  def CreateVtkMultiBlockDataWriter(self):
    mypid = SmartGetLocalProcessId()
    numproc = SmartGetNumberOfProcesses()

    if PhactoriDbg(100):
      myDebugPrint3("numproc " + str(numproc) + " mypid " + str(mypid) + "\n")

    self.mWriter = vtk.vtkXMLPMultiBlockDataWriter()
    self.mWriter.SetInputData(self.mFilterToWriteDataFrom.GetClientSideObject().GetOutputDataObject(0))
    self.mWriter.SetNumberOfPieces(numproc)
    self.mWriter.SetStartPiece(mypid)

    if PhactoriDbg(100):
      myDebugPrint3("self.mWriter: " + str(self.mWriter) + "\n")

  def CreateVtkPolyDataWriter(self):
    mypid = SmartGetLocalProcessId()

    pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
    globalController = pm.GetGlobalController()
    #gLocalProcessId = globalController.GetLocalProcessId()
    numproc = globalController.GetNumberOfProcesses()

    if PhactoriDbg(100):
      myDebugPrint3("numproc " + str(numproc) + " mypid " + str(mypid) + "\n")
    self.mWriter = vtk.vtkXMLPPolyDataWriter()
    self.mWriter.SetInputData(self.mFilterToWriteDataFrom.GetClientSideObject().GetOutputDataObject(0))
    self.mWriter.SetNumberOfPieces(numproc)
    self.mWriter.SetStartPiece(mypid)
    self.mWriter.SetEndPiece(mypid)
    self.mWriter.SetUseSubdirectory(True)

  def CreateVtkUnstructuredDataWriter(self):
    mypid = SmartGetLocalProcessId()

    #pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
    #globalController = pm.GetGlobalController()
    #gLocalProcessId = globalController.GetLocalProcessId()
    #numproc = globalController.GetNumberOfProcesses()

    #if PhactoriDbg(100):
    #  myDebugPrint3("numproc " + str(numproc) + " mypid " + str(mypid) + "\n")
    self.mWriter = vtk.vtkXMLUnstructuredGridWriter()
    self.mWriter.SetInputData(self.mFilterToWriteDataFrom.GetClientSideObject().GetOutputDataObject(0))

  def CreateCurrentCallbackOutputFilename(self):
    ftc = GetFrameTagCounter()
    outName = self.mOutputFileBasedirectory + "/" + \
      self.mOutputFileBasename + "_" + \
      str(ftc).zfill(self.mNumCounterDigits) + \
      self.mVtkDataTypeExtention
    return outName

  def ExportOperationData(self, datadescription):
    """overrides parent class method in order to output .vtm or .vtp
       (vtkMultiBlockData or vtkPolyData) files"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriVtkDataExportOperation." \
          "ExportOperationData entered\n", 100)

    self.mFilterToWriteDataFrom.UpdatePipeline()

    if self.mWriter == None:
      if(self.mVtkDataType == "vtkmultiblockdataset"):
        self.CreateVtkMultiBlockDataWriter()
      elif(self.mVtkDataType == "vtkpolydata"):
        self.CreateVtkPolyDataWriter()
      elif(self.mVtkDataType == "vtkunstructureddata"):
        self.CreateVtkUnstructuredDataWriter()

    outfilename = self.CreateCurrentCallbackOutputFilename()
    if PhactoriDbg(100):
      myDebugPrint3("outfilename: " + str(outfilename) + "\n")
    self.mWriter.SetFileName(outfilename)

    #if PhactoriDbg(100):
    #  myDebugPrint3("multiblock info right before self.mWriter.Write (start):\n")
    #  RecursivelyPrintMultiBlockStats(self.mFilterToWriteDataFrom)
    #  myDebugPrint3("multiblock info right before self.mWriter.Write (end):\n")
    #  BarrierLock("barrier lock before self.mWriter.Write")
    self.mWriter.Write()

    if self.mDataArtifactLister != None:
      self.mDataArtifactLister.AddDataExportToDataArtifactOutputList(outfilename)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriVtkDataExportOperation." \
          "ExportOperationData returning\n", 100)

#phactori_combine_to_single_python_file_subpiece_end_1

