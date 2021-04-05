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

#phactori_combine_to_single_python_file_subpiece_begin_1

class PhactoriExodusIIExportOperation(PhactoriOperationSpecifics):
  """passes through it's input to it's output, but also does an
     ExportOperation which is an exodusII writer, .exo """
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mOutputFileBasename = "myexo1"
    self.mOutputFileBasedirectory = "CatalystExodusIIOutput"
    self.mFilterToWriteDataFrom = None
    self.mWriter = None
    self.mNumCounterDigits = 4
    #should be "vtkmultiblockdataset" or "vtkpolydata"
    self.mExodusIIType = "vtkmultiblockdataset"
    self.mExodusIITypeExtention = ".exo"
    #self.mExodusIITypeExtention = ".e-s"

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

  def CreateParaViewFilter(self, inInputFilter):
    self.mFilterToWriteDataFrom = inInputFilter
    return self.mFilterToWriteDataFrom

  def CreateCurrentCallbackOutputFilename(self):
    ftc = GetFrameTagCounter()
    outName = self.mOutputFileBasedirectory + "/" + \
      self.mOutputFileBasename + "_" + \
      str(ftc).zfill(self.mNumCounterDigits) + \
      self.mExodusIITypeExtention
    return outName

  def ExportOperationData(self, datadescription):
    """overrides parent class method in order to output .vtm or .vtp
       (vtkMultiBlockData or vtkPolyData) files"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExodusIIExportOperation." \
          "ExportOperationData entered\n", 100)

    self.mFilterToWriteDataFrom.UpdatePipeline()

    #outfilename = "MySlice5." + str(GetFrameTagCounter()) + ".e-s"
    outfilename = self.CreateCurrentCallbackOutputFilename()
    if PhactoriDbg(100):
      myDebugPrint3("outfilename: " + str(outfilename) + "\n")
    if self.mWriter == None:
      self.mWriter = ExodusIIWriter(Input=self.mFilterToWriteDataFrom, IgnoreMetaDataWarning=True, FileName=outfilename)
      #exodusIIWriter1 = servermanager.writers.ExodusIIWriter(Input=slice1)
      #coprocessor.RegisterWriter(exodusIIWriter1, filename='Slice1.e-s', freq=1, paddingamount=0, DataMode='None', HeaderType='None', EncodeAppendedData=None, CompressorType='None', CompressionLevel='None')
    else:
      self.mWriter.FileName = outfilename


    ##self.mWriter.UpdatePipeline()
    UpdatePipelineWithCurrentTimeArgument(self.mWriter)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExodusIIExportOperation." \
          "ExportOperationData returning\n", 100)

#phactori_combine_to_single_python_file_subpiece_end_1

