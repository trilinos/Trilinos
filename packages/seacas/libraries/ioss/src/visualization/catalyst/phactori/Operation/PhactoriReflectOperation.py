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
class PhactoriReflectOperation(PhactoriOperationSpecifics):
  """manages reflect filter, also inserts group filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mPlane = "X Min"
    self.mCenter = 0.0

    #we keep internal reflect filter, because group filter is the one exposed
    #to the rest of the world
    self.mInternalReflectFilter = None

  def ParseParametersFromJson(self, inJson):
    if 'plane' in inJson:
      self.mPlane = inJson['plane']
    else:
      myDebugPrint3AndException(
          "PhactoriReflectOperation::ParseParametersFromJson\n"
          "Error:  must have 'plane' key\n")
    if 'center' in inJson:
      self.mCenter = inJson['center']

  def CreateParaViewFilter(self, inInputFilter):
    """create the reflect (and group) filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriReflectOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    self.mInternalReflectFilter = Reflect(inInputFilter)
    self.mInternalReflectFilter.Plane  = self.mPlane
    self.mInternalReflectFilter.Center  = self.mCenter
    #don't know if this is necessary here
    UpdatePipelineWithCurrentTimeArgument(self.mInternalReflectFilter)

    #now grop this new source with the original non-reflected one
    newParaViewFilter = GroupDatasets()
    newParaViewFilter.Input = [self.mInternalReflectFilter, inInputFilter]

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3("reflect plane: " + str(self.mInternalReflectFilter.Plane) + "\n"
          "reflect center: " + str(self.mInternalReflectFilter.Center) + "\n")
      myDebugPrint3("PhactoriReflectOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
