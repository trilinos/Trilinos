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
class PhactoriVectorProject(PhactoriOperationSpecifics):
  """manages calculator filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.projectionList = None
    self.pvCalcs = None

  def ParseParametersFromJson(self, inJson):
    key1 = "project list"
    if key1 not in inJson:
      myDebugPrint3AndException(
        "PhactoriVectorProject::ParseParametersFromJson\n"
        "Error: must have key " + key1 + "\n")
    self.projectionList = inJson[key1]
    if len(self.projectionList) % 3 != 0:
      myDebugPrint3AndException(
        "PhactoriVectorProject::ParseParametersFromJson\n"
        "projection list must have three entries per projection")

  def CreateParaViewFilter(self, inInputFilter):
    """create the projection calculators for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriVectorProject.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    ndx = 0
    testend = len(self.projectionList) - 1
    self.pvCalcs = []
    nextInput = inInputFilter
    while ndx < testend:
      vecToProj = self.projectionList[ndx]
      vecToProjOnTo = self.projectionList[ndx+1]
      targetProjLenName = self.projectionList[ndx+2]
      ndx += 3
      functionstr = vecToProj + "." + vecToProjOnTo
      #functionstr = \
      #  vecToProj + "_X*" + vecToProjOnTo + "_X+" + \
      #  vecToProj + "_Y*" + vecToProjOnTo + "_Y+" + \
      #  vecToProj + "_Z*" + vecToProjOnTo + "_Z"
      if PhactoriDbg(100):
        myDebugPrint3("ndx " + str(ndx) + "\n"
          "vecToProj " + str(vecToProj) + "\n"
          "vecToProjOnTo " + str(vecToProjOnTo) + "\n"
          "targetProjLenName " + str(targetProjLenName) + "\n"
          "functionstr " + str(functionstr) + "\n")
      newPvCalc = Calculator(nextInput)
      newPvCalc.Function = functionstr
      newPvCalc.ResultArrayName = targetProjLenName
      #newParaViewFilter.AttributeType = self.mPointOrCell
      newPvCalc.AttributeType = "Point Data"
      self.pvCalcs.append(newPvCalc)
      nextInput = newPvCalc
      UpdatePipelineWithCurrentTimeArgument(newPvCalc)

    newParaViewFilter = self.pvCalcs[-1]
    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3(
        "PhactoriTransformOperation.CreateParaViewFilter returning\n", 100)

    return newParaViewFilter
#phactori_combine_to_single_python_file_subpiece_end_1
