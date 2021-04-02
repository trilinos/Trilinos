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
class PhactoriWarpByVectorOperation(PhactoriOperationSpecifics):
  """manages warpbyvector filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mVectorName = ""
    self.mVariableType = "node"
    self.mScaleFactor = 1.0

  def ParseParametersFromJson(self, inJson):
    if 'variable type' in inJson:
      self.mVariableType = inJson['variable type']
    elif 'variablename' in inJson:
      self.mVectorName = inJson['variablename']
    elif 'variable vector' in inJson:
      self.mVectorName = inJson['variable vector']
    else:
      myDebugPrint3AndException(
          "PhactoriWarpByVectorOperation::ParseParametersFromJson\n"
          "Error:  must have 'vector name' key\n")
    if 'scale' in inJson:
      self.mScaleFactor = inJson['scale']

  def CreateParaViewFilter(self, inInputFilter):
    """create the warp by vector filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriWarpByVectorOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    newParaViewFilter = WarpByVector(inInputFilter)
    newParaViewFilter.ScaleFactor = self.mScaleFactor
    if self.mVariableType == 'cell':
      varTypeString = 'CELLS'
    else:
      varTypeString = 'POINTS'
    newParaViewFilter.Vectors = [varTypeString, self.mVectorName]

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3("vector name: " + self.mVectorName + "\n"
          "variable type: " + varTypeString + "\n"
          "scale factor: " + str(newParaViewFilter.ScaleFactor) + "\n")
      myDebugPrint3("PhactoriWarpByVectorOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
