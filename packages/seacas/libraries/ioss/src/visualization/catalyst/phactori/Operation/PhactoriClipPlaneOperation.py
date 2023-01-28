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
class PhactoriClipPlaneOperation(PhactoriPlaneOpBase):
  """clip plane operation, adapter to the catalyst filter

PhactoriSliceWithPlaneOperation is the phactori manager for working with the
ParaView/Catalyst Clip() filter and its parameters, providing
access and pipeline/input/output managment via the json, lexx/yacc, or soon
yaml interface. The user may specify a named input to the filter, with the
unnamed default being the incoming data mesh. This class is confined to a
ParaView Clip operation with a clip type of plane. See
PhactoriBoxClipOperation and PhactoriCylinderClipOperation for clipping with
a type of box and cylinder, and other classes may be added for other clip
types.

For information on the Clip() filter from ParaView, see the ParaView
Documentation.

The user must define the plane with a point and normal or with three points.
There are defaults which will be used if the user does not supply some or
all of the definition.  PhactoriClipPlaneOperation is a child class of
PhactoriPlaneOpBase, along with PhactoriSliceWithPlaneOperation and
PhactoriSliceOperation. Check the documentation for PhactoriPlaneOpBase
for the many options for defining the plane point(s), including absolute or
relative 3D points, dynamic data-driven point locations, or collocating with
mesh nodes or elements (or offset therefrom).

The user may also choose to make a "crinkle" or "smooth" cut with the
"cut type" key, defaulting to
smooth. A smooth cut makes a two dimentional plane, while a crinkle cut
includes all elements cut by the plane and creates no new elements which
were not in the original mesh. The user may also choose to keep either side
of the plane by using the key "side to keep" and setting to "positive" or
"negative".

To add a PhactoriClipPlaneOperation to the incoming script, you add
a sub-block to the "operation blocks" section of the data with the "type"
key given a value of "clip". One complete but simple example
script:

::

  {
    "camera blocks":{"myclipcam1":{"type":"camera", "look direction":[1.0, 2.0, 3.0]}},
    "representation blocks":{"rep_tmprtr":{"color by scalar":"temperature"}},
    "imageset blocks":{
      "temperature_clip_1":{
        "operation":"myclipwithplane1",
        "camera":"myclipcam1",
        "representation":"rep_tmprtr",
        "image basedirectory":"CatalystOutput",
        "image basename":"clip1_temperature."
      }
    },
    "operation blocks":{
      "myclipwithplane1":{
        "type":"clip",
        "cut type":"crinkle",
        "side to keep":"negative",
        "relative point on plane":[0.1, -0.2, 0.3],
        "plane normal":[1.0, 2.0, 3.0]
      }
    }

The "point on plane" is a PhactoriUserPointInfo instance, which means it can
be an absolute point, a relative point, a data point (at min or max of a data
array), at a node/point, at an element/cell, or at a data min max or cell or
element but displaced by a vector. See the PhactoriUserPointInfo class for more
detailed info.

Further, the plane may also be define with three points instead of a point and
a normal. To use this format, omit the "plane normal" and "... point on plane"
keys and instead have keys with "point on plane A", "point and plane B", and
"point on plane C". For example:

::

  "absolute point on plane A":[0.0, 1.0, 2.0],
  "relative point on plane B":[0.0, 0.2, -0.3],
  "relative point on plane C":[0.0, 0.0, -0.0],

"""

  def CreateParaViewFilter(self, inInputFilter):

    #don't need our own init code at this point, but this is how it would be
    #added
    #def __init__(self):
    #    MySuperClass.__init__(self)

    """create the clip plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriClipPlaneOperation.CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()
    newParaViewFilter = Clip(Input = inInputFilter, ClipType = "Plane")
    self.UpdateClip(inInputFilter, newParaViewFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriClipPlaneOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

  def DoUpdateDueToChangeInData(self, inIncomingPvFilter,
      outOutgoingPvFilter):
    """the PhactoriClipPlaneOperation may need to update if the point on
       the clip plane was tied to a node, element, or variable min/max
       location"""
    if PhactoriDbg():
      myDebugPrint3("PhactoriClipPlaneOperation::"
          "DoUpdateDueToChangeInData override executing\n")

    if self.MayChangeWithData() == False:
      if PhactoriDbg():
        myDebugPrint3("PhactoriClipPlaneOperation::"
            "DoUpdateDueToChangeInData returning (absolute point or points)\n")
      return

    self.UpdateClip(inIncomingPvFilter, outOutgoingPvFilter)

    if PhactoriDbg():
      myDebugPrint3("PhactoriClipPlaneOperation::"
          "DoUpdateDueToChangeInData override returning\n")

  def UpdateClip(self, inIncomingPvFilter, ioOutgoingPvFilter):
    """using the current info on the clip, get all the paraview stuff
       set up correctly"""

    if PhactoriDbg():
      myDebugPrint3("PhactoriClipPlaneOperation::UpdateClip entered\n")

    originToUse = [0,0,0]
    normalToUse = [0,1,0]
    self.CalculateUpdatedOriginAndNormal(
            inIncomingPvFilter, originToUse, normalToUse)

    if PhactoriDbg():
      myDebugPrint3('  updateclip using normal: ' + \
              str(normalToUse) + '\n')
    ioOutgoingPvFilter.ClipType.Normal = normalToUse

    if PhactoriDbg():
      myDebugPrint3('  updateclip using origin: ' + str(originToUse) + '\n')
    ioOutgoingPvFilter.ClipType.Origin = originToUse

    #these aren't changing yet
    ioOutgoingPvFilter.Crinkleclip = self.mCrinkleSetting
    if gParaViewCatalystVersionFlag < 50502:
      ioOutgoingPvFilter.InsideOut = self.mInsideOut
    else:
      ioOutgoingPvFilter.Invert = self.mInvert

    if PhactoriDbg():
      myDebugPrint3("PhactoriClipPlaneOperation::UpdateClip returning\n")
#phactori_combine_to_single_python_file_subpiece_end_1

