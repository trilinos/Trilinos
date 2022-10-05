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

class PhactoriCSVExportOperation(PhactoriOperationSpecifics):
  """passes through it's input to it's output, but also does an
     ExportOperation which is an csv writer, .exo """
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mOutputFileBasename = "myexo1"
    self.mOutputFileBasedirectory = "CatalystCSVOutput"
    self.mFilterToWriteDataFrom = None
    self.mCsvCellWriter = None
    self.mCsvPointWriter = None
    self.mNumCounterDigits = 4
    self.mCSVTypeExtention = ".csv"
    self.Precision = 9

  def ParseOutputFilenameParametersFromJson(self, inJson):
    self.mOutputFileBasename = getParameterFromBlock(inJson,
      "basename", self.mPhactoriOperationBlockOwner.mName)

    self.mOutputFileBasedirectory = getParameterFromBlock(inJson,
      "basedirectory", self.mOutputFileBasedirectory)
    CreateDirectoryFromProcessZero(self.mOutputFileBasedirectory)

    self.mNumCounterDigits = getParameterFromBlock(inJson,
      "name digit count", self.mNumCounterDigits)

    self.Precision = getParameterFromBlock(inJson,
      "precision", self.Precision)
    if (self.Precision < 1) or (self.Precision > 100):
      myDebugPrint3AndException("PhactoriCSVExportOperation:\n"
        "precision must be 1-100, not " + str(self.Precision) + "\n")
    

  def ParseParametersFromJson(self, inJson):
    self.ParseOutputFilenameParametersFromJson(inJson)

  def CreateParaViewFilter(self, inInputFilter):
    self.mFilterToWriteDataFrom = inInputFilter
    return self.mFilterToWriteDataFrom

  def CreateCurrentCallbackOutputFilename(self):
    ftc = GetFrameTagCounter()
    outName_cell = self.mOutputFileBasedirectory + "/" + \
      self.mOutputFileBasename + "_cell_" + \
      str(ftc).zfill(self.mNumCounterDigits) + \
      self.mCSVTypeExtention
    outName_point = self.mOutputFileBasedirectory + "/" + \
      self.mOutputFileBasename + "_point_" + \
      str(ftc).zfill(self.mNumCounterDigits) + \
      self.mCSVTypeExtention
    return (outName_cell, outName_point)

  def ExportOperationData(self, datadescription):
    """overrides parent class method in order to output .vtm or .vtp
       (vtkMultiBlockData or vtkPolyData) files"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriCSVExportOperation." \
          "ExportOperationData entered\n", 100)

    self.mFilterToWriteDataFrom.UpdatePipeline()

    outfilenames = self.CreateCurrentCallbackOutputFilename()
    if PhactoriDbg(100):
      myDebugPrint3("outfilenames: " + str(outfilenames) + "\n")

    if self.mCsvCellWriter == None:
      self.mCsvCellWriter = CSVWriter(Input=self.mFilterToWriteDataFrom,
        Precision=self.Precision, FieldAssociation="Cell Data")
    if self.mCsvPointWriter == None:
      self.mCsvPointWriter = CSVWriter(Input=self.mFilterToWriteDataFrom,
        Precision=self.Precision, FieldAssociation="Point Data")

    if self.mCsvCellWriter != None:
      self.mCsvCellWriter.FileName = outfilenames[0]
      UpdatePipelineWithCurrentTimeArgument(self.mCsvCellWriter)
    if self.mCsvPointWriter != None:
      self.mCsvPointWriter.FileName = outfilenames[1]
      UpdatePipelineWithCurrentTimeArgument(self.mCsvPointWriter)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriCSVExportOperation." \
          "ExportOperationData returning\n", 100)

#phactori_combine_to_single_python_file_subpiece_end_1

