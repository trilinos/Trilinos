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
class PhactoriSliceOperation(PhactoriPlaneOpBase):
  """slice operation, settings and handling of slice filter"""

  def CreateParaViewFilter(self, inInputFilter):

    #don't need our own init code at this point, but this is how it would be
    #added
    #def __init__(self):
    #    MySuperClass.__init__(self)

    """create the slice plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriSliceOperation.CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()
    newParaViewFilter = Slice(Input = inInputFilter, SliceType = "Plane")
    newParaViewFilter.Triangulatetheslice = 0

    self.UpdateSlice(inInputFilter, newParaViewFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriSliceOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

  def DoUpdateDueToChangeInData(self, inIncomingPvFilter,
      outOutgoingPvFilter):
    """the PhactoriSliceOperation may need to update if the point on
       the slice plane was tied to a node, element, or variable min/max
       location"""
    if PhactoriDbg():
      myDebugPrint3("PhactoriSliceOperation::"
          "DoUpdateDueToChangeInData override executing\n")

    if self.MayChangeWithData() == False:
      if PhactoriDbg():
        myDebugPrint3("PhactoriSlicePlaneOperation::"
            "DoUpdateDueToChangeInData returning (absolute point or points)\n")
      return

    self.UpdateSlice(inIncomingPvFilter, outOutgoingPvFilter)

    if PhactoriDbg():
      myDebugPrint3("PhactoriSlicePlaneOperation::"
          "DoUpdateDueToChangeInData override returning\n")

  def UpdateSlice(self, inIncomingPvFilter, ioOutgoingPvFilter):
    """using the current info on the slice, get all the paraview stuff
       set up correctly"""

    if PhactoriDbg():
      myDebugPrint3("PhactoriSlicePlaneOperation::UpdateSlice entered\n")

    originToUse = [0,0,0]
    normalToUse = [0,1,0]
    self.CalculateUpdatedOriginAndNormal(
            inIncomingPvFilter, originToUse, normalToUse)

    if PhactoriDbg():
      myDebugPrint3('  updateslice using normal: ' + \
              str(normalToUse) + '\n')
    ioOutgoingPvFilter.SliceType.Normal = normalToUse

    if PhactoriDbg():
      myDebugPrint3('  updateslice using origin: ' + str(originToUse) + '\n')
    ioOutgoingPvFilter.SliceType.Origin = originToUse

    #these aren't changing yet
    ioOutgoingPvFilter.Crinkleslice = self.mCrinkleSetting

    if PhactoriDbg():
      myDebugPrint3("PhactoriSlicePlaneOperation::UpdateSlice returning\n")

#phactori_combine_to_single_python_file_subpiece_end_1

