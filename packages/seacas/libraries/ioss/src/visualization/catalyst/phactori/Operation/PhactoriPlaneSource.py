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

class PhactoriPlaneSource(PhactoriOperationSpecifics):
  """Phactori interface to the ParaView Plane Geometric Source. Allows user to
     specify a plane (rectangle) by 3 points in space and a tesselation
     resolution in each dimension of the rectangle"""

  def __init__(self):
    self.PlaneOrigin = [-0.5, -0.5, 0.0]
    self.PlanePoint1 = [ 0.5, -0.5, 0.0]
    self.PlanePoint2 = [-0.5,  0.5, 0.0]
    self.PlaneXResolution = 1
    self.PlaneYResolution = 1

  def ParseParametersFromJson(self, inJson):
    key1 = "origin"
    if key1 in inJson:
      self.PlaneOrigin = inJson[key1]
      try:
        if len(self.PlaneOrigin) != 3:
          myDebugPrint3AndException(
            "PhactoriPlaneSource:ParseParametersFromJson\n"
            "bad 'origin', must be 3 number list list [1.1,2.3,-4.5]\n")
      except:
        myDebugPrint3AndException(
          "PhactoriPlaneSource:ParseParametersFromJson\n"
          "bad 'origin', must be 3 number list list [1.1,2.3,-4.5]\n")

    key1 = "point1"
    if key1 in inJson:
      self.PlanePoint1 = inJson[key1]
      try:
        if len(self.PlanePoint1) != 3:
          myDebugPrint3AndException(
            "PhactoriPlaneSource:ParseParametersFromJson\n"
            "bad 'point1', must be 3 number list list [1.1,2.3,-4.5]\n")
      except:
        myDebugPrint3AndException(
          "PhactoriPlaneSource:ParseParametersFromJson\n"
          "bad 'point1', must be 3 number list list [1.1,2.3,-4.5]\n")

    key1 = "point2"
    if key1 in inJson:
      self.PlanePoint2 = inJson[key1]
      try:
        if len(self.PlanePoint2) != 3:
          myDebugPrint3AndException(
            "PhactoriPlaneSource:ParseParametersFromJson\n"
            "bad 'point2', must be 3 number list list [1.1,2.3,-4.5]\n")
      except:
        myDebugPrint3AndException(
          "PhactoriPlaneSource:ParseParametersFromJson\n"
          "bad 'point2', must be 3 number list list [1.1,2.3,-4.5]\n")

    key1 = "xresolution"
    if key1 in inJson:
      self.PlaneXResolution = inJson[key1]
      if isinstance(self.PlaneXResolution, int):
        if self.PlaneXResolution < 1:
          myDebugPrint3AndException(
            "PhactoriPlaneSource:ParseParametersFromJson\n"
            "'xresolution' must be integer >= 1\n")
      else:
        myDebugPrint3AndException(
          "PhactoriPlaneSource:ParseParametersFromJson\n"
          "'xresolution' must be integer >= 1\n")

    key1 = "yresolution"
    if key1 in inJson:
      self.PlaneYResolution = inJson[key1]
      if isinstance(self.PlaneYResolution, int):
        if self.PlaneYResolution < 1:
          myDebugPrint3AndException(
            "PhactoriPlaneSource:ParseParametersFromJson\n"
            "'yresolution' must be integer >= 1\n")
      else:
        myDebugPrint3AndException(
          "PhactoriPlaneSource:ParseParametersFromJson\n"
          "'yresolution' must be integer >= 1\n")

  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPlaneSource.CreateParaViewFilter entered\n", 100)

    savedActiveSource = GetActiveSource()

    newParaViewFilter = Plane()

    newParaViewFilter.Origin = self.PlaneOrigin
    newParaViewFilter.Point1 = self.PlanePoint1
    newParaViewFilter.Point2 = self.PlanePoint2
    newParaViewFilter.XResolution = self.PlaneXResolution
    newParaViewFilter.YResolution = self.PlaneYResolution

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3(
        "Origin: " + str(newParaViewFilter.Origin) + "\n" + \
        "Point1: " + str(newParaViewFilter.Point1) + "\n" + \
        "Point2: " + str(newParaViewFilter.Point2) + "\n" + \
        "XResolution: " + str(newParaViewFilter.XResolution) + "\n" + \
        "YResolution: " + str(newParaViewFilter.YResolution) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPlaneSource.CreateParaViewFilter returning\n", 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
