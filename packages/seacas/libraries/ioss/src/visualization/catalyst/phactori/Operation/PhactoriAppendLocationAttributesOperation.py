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
class PhactoriAppendLocationAttributesOperation(PhactoriOperationSpecifics):
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.AppendPointLocations = 1
    self.AppendCellCenters = 1
    return

  def ParseParametersFromJson(self, inJson):
    key1 = "append point locations"
    if key1 in inJson:
      testval = inJson[key1]
      if testval:
        self.AppendPointLocations = 1
      else:
        self.AppendPointLocations = 0
    else:
      key1 = "append_point_locations"
      if key1 in inJson:
        testval = inJson[key1]
        if testval:
          self.AppendPointLocations = 1
        else:
          self.AppendPointLocations = 0
    
    key2 = "append cell centers"
    if key2 in inJson:
      testval = inJson[key2]
      if testval:
        self.AppendCellCenters = 1
      else:
        self.AppendCellCenters = 0
    else:
      key2 = "append_cell_centers"
      if key2 in inJson:
        testval = inJson[key2]
        if testval:
          self.AppendCellCenters = 1
        else:
          self.AppendCellCenters = 0

    return

  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriAppendLocationAttributesOperation:CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked
    
    if PhactoriDbg(100):
      myDebugPrint3('about to call UpdatePipelineWithCurrentTimeArgument\n', 100)
    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    savedActiveSource = GetActiveSource()
    if PhactoriDbg(100):
      myDebugPrint3('about to call AppendLocationAttributes\n', 100)
    newParaViewFilter = AppendLocationAttributes(Input = inInputFilter)

    newParaViewFilter.AppendPointLocations = self.AppendPointLocations
    newParaViewFilter.AppendCellCenters = self.AppendCellCenters

    SetActiveSource(newParaViewFilter)
    if PhactoriDbg(100):
      myDebugPrint3('about to call UpdatePipelineWithCurrentTimeArgument\n', 100)
    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriAppendLocationAttributesOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
