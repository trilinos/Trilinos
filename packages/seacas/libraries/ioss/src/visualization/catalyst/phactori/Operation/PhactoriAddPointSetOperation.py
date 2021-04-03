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

#phactori_combine_to_single_python_file_subpiece_begin_1
class PhactoriAddPointSetOperation(PhactoriOperationSpecifics):
  """filter/operation which reads in a .csv file which is a set of points and
     creates a new source from that and groups it with the input source"""
  def __init__(self):
    self.mXColumn = 'x column'
    self.mYColumn = 'y column'
    self.mZColumn = 'z column'
    self.m2dPointsFlag = False
    self.mFieldDelimeterCharacters = ','
    self.mFilename = None
    self.mInternalPvCSVReader = None
    self.mInternalPvTableToPoints = None

  def ParseParametersFromJson(self, inJson):
    if 'filename' in inJson:
      self.mFilename = inJson['filename']
    else:
      myDebugPrint3AndException(
          "PhactoriAddPointSetOperation::ParseParametersFromJson\n"
          "Error:  must have 'filename' key\n")
    if 'x column' in inJson:
      self.mXColumn = inJson['x column']
    if 'y column' in inJson:
      self.mYColumn = inJson['y column']
    if 'z column' in inJson:
      self.mZColumn = inJson['z column']
    if '2d points flag' in inJson:
      self.m2dPointsFlag = inJson['2d points flag']
    if 'field delimeter' in inJson:
      self.mFieldDelimeterCharacters = inJson['field delimeter']

  def CreateParaViewFilter(self, inInputFilter):
    """create the read pointset (and group) filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriAddPointSetOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    #csv reader
    self.mInternalPvCSVReader = CSVReader( FileName=[self.mFilename] )
    self.mInternalPvCSVReader.FieldDelimiterCharacters = self.mFieldDelimeterCharacters
    #don't know if this is necessary here
    UpdatePipelineWithCurrentTimeArgument(self.mInternalPvCSVReader)

    #table to points
    self.mInternalPvTableToPoints = TableToPoints(self.mInternalPvCSVReader)
    self.mInternalPvTableToPoints.XColumn = self.mXColumn
    self.mInternalPvTableToPoints.YColumn = self.mYColumn
    self.mInternalPvTableToPoints.ZColumn = self.mZColumn
    if self.m2dPointsFlag:
      self.mInternalPvTableToPoints.a2DPoints = 1
    else:
      self.mInternalPvTableToPoints.a2DPoints = 0
    #don't know if this is necessary here
    UpdatePipelineWithCurrentTimeArgument(self.mInternalPvTableToPoints)

    #now group this new point set source with the original one
    newParaViewFilter = GroupDatasets()
    newParaViewFilter.Input = [self.mInternalPvTableToPoints, inInputFilter]
    #newParaViewFilter = self.mInternalPvTableToPoints

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3(
          "filename: " + str(self.mInternalPvCSVReader.FileName) + "\n"
          "delimeter: -->" + str(self.mInternalPvCSVReader.FieldDelimiterCharacters) + "<--\n"
          "x column: " + str(self.mInternalPvTableToPoints.XColumn) + "\n"
          "y column: " + str(self.mInternalPvTableToPoints.YColumn) + "\n"
          "z column: " + str(self.mInternalPvTableToPoints.ZColumn) + "\n"
          "2d points flag: " + str(self.mInternalPvTableToPoints.a2DPoints) + "\n")
      myDebugPrint3("PhactoriAddPointSetOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
