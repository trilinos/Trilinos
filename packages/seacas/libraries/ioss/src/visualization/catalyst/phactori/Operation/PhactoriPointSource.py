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
from .PhactoriMpiUtilities import *

#phactori_combine_to_single_python_file_subpiece_begin_1

class PhactoriPointSource(PhactoriOperationSpecifics):
  """Filter/operation which reads in a .json file which is a list of 3d
     coordinates and creates a new point source from that list. Resulting
     point source will have 1 element an N points. This source is intented
     to work correctly in parallel for Catalyst or pvbatch symmetric mode.
     The json file which is read in will only be read on one process and
     mpi broadcast is used to distribute the list, rather than having each
     process read the json file."""
  def __init__(self):
    self.UseDefaultRadius = True
    self.Radius = 0.0
    self.UseDefaultNumPoints = True
    self.NumPoints = 1
    self.UseDefaultCenter = True
    self.Center = [0.0,0.0,0.0]

  def ParseParametersFromJson(self, inJson):
    key1 = "radius"
    if key1 in inJson:
      self.UseDefaultRadius = False
      self.Radius = float(inJson[key1])
    key1 = "number of points"
    if key1 in inJson:
      self.UseDefaultNumPoints = False
      self.NumPoints = int(inJson[key1])
    key1 = "center"
    if key1 in inJson:
      self.UseDefaultCenter = False
      self.Center = inJson[key1]
      if len(self.Center) != 3:
        myDebugPrint3AndException("PhactoriPointSource:ParseParametersFromJson:\n"
          "bad center, need like [0.1,2.3,4.5]\n")

  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSource.CreateParaViewFilter entered\n", 100)

    savedActiveSource = GetActiveSource()

    newParaViewFilter = PointSource()
    if self.UseDefaultRadius == False:
      newParaViewFilter.Radius = self.Radius
    if self.UseDefaultNumPoints == False:
      newParaViewFilter.NumberOfPoints = self.NumPoints
    if self.UseDefaultCenter == False:
      newParaViewFilter.Center = self.Center

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("Radius: " + str(newParaViewFilter.Radius) + "\n" + \
        "NumberOfPoints: " + str(newParaViewFilter.NumberOfPoints) + "\n" + \
        "Center: " + str(newParaViewFilter.Center) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSource.CreateParaViewFilter returning\n", 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1

