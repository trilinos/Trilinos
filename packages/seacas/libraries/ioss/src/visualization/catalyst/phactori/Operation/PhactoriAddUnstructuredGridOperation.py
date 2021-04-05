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
class PhactoriAddUnstructuredGridOperation(PhactoriOperationSpecifics):
  """filter/operation which reads in a .vtu file which is an unstructured
     grid and creates a new source from that and groups it with the
     input source"""
  def __init__(self):
    self.mFilename = None
    self.mInternalPvUnstructuredGridReader = None

  def ParseParametersFromJson(self, inJson):
    if 'filename' in inJson:
      self.mFilename = inJson['filename']
    else:
      myDebugPrint3AndException(
          "PhactoriAddUnstructuredGridOperation::ParseParametersFromJson\n"
          "Error:  must have 'filename' key\n")

  def CreateParaViewFilter(self, inInputFilter):
    """create the read unstructured grid (and group) filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriAddUnstructuredGridOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    #csv reader
    self.mInternalPvUnstructuredGridReader = XMLUnstructuredGridReader( FileName=[self.mFilename] )

    #don't know if this is necessary here
    UpdatePipelineWithCurrentTimeArgument(self.mInternalPvUnstructuredGridReader)

    #now group this new point set source with the original one
    newParaViewFilter = GroupDatasets()
    newParaViewFilter.Input = [self.mInternalPvUnstructuredGridReader, inInputFilter]

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3(
          "filename: " + str(self.mInternalPvUnstructuredGridReader.FileName) + "\n")
      myDebugPrint3("PhactoriAddUnstructuredGridOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter
#phactori_combine_to_single_python_file_subpiece_end_1
