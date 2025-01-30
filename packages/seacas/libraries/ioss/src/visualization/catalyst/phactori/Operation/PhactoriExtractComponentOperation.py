# Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
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
class PhactoriExtractComponentOperation(PhactoriOperationSpecifics):
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.Component = 0
    self.InputArrayName = "None"
    self.InputArrayCellsOrPoints = "POINTS"
    self.OutputArrayName = "Result"
    return

  def ParseParametersFromJson(self, inJson):
    keyval1 = "component"
    if keyval1 in inJson:
      self.Component = inJson[keyval1]

    keyval2 = "output array name"
    if keyval2 in inJson:
      self.OutputArrayName = inJson[keyval2]

    keyval3 = "input array name"
    if keyval3 in inJson:
      self.InputArrayName = inJson[keyval3]
    else:
      myDebugPrint3AndException(
        "PhactoriExtractComponentOperation:ParseParametersFromJson\n"
        "must have keys 'input array name' and 'input array type'\n")

    keyval4 = "input array type"
    if keyval4 in inJson:
      typeval = inJson[keyval4]
      if typeval == "points":
        self.InputArrayCellsOrPoints = "POINTS"
      elif typeval == "cells":
        self.InputArrayCellsOrPoints = "CELLS"
      elif typeval == "nodes":
        self.InputArrayCellsOrPoints = "POINTS"
      elif typeval == "elements":
        self.InputArrayCellsOrPoints = "CELLS"
      else:
        myDebugPrint3AndException(
          "PhactoriExtractComponentOperation:ParseParametersFromJson\n"
          "'input array type' must be 'points' 'cells' 'nodes' or 'elements'\n")
    else:
      myDebugPrint3AndException(
        "PhactoriExtractComponentOperation:ParseParametersFromJson\n"
        "must have keys 'input array name' and 'input array type'\n")

  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractComponentOperation:CreateParaViewFilter entered\n", 100)
    #info in block class should already be parsed and checked

    if PhactoriDbg(100):
      myDebugPrint3("about to call UpdatePipelineWithCurrentTimeArgument\n", 100)
    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    savedActiveSource = GetActiveSource()
    newParaViewFilter = ExtractComponent(Input = inInputFilter)
    newParaViewFilter.Component = self.Component
    newParaViewFilter.OutputArrayName = self.OutputArrayName
    newParaViewFilter.InputArray = [self.InputArrayCellsOrPoints, self.InputArrayName]
    if PhactoriDbg(100):
      myDebugPrint3(
        "newParaViewFilter.Component" + str(newParaViewFilter.Component) + "\n"
        "newParaViewFilter.OutputArrayName" + str(newParaViewFilter.OutputArrayName) + "\n"
        "newParaViewFilter.InputArray" + str(newParaViewFilter.InputArray) + "\n", 100)

    SetActiveSource(newParaViewFilter)
    if PhactoriDbg(100):
      myDebugPrint3("about to call UpdatePipelineWithCurrentTimeArgument\n", 100)
    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractComponentOperation.CreateParaViewFilter returning\n", 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
