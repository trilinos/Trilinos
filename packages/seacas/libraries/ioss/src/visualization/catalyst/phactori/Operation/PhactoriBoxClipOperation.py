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
class PhactoriBoxClipOperation(PhactoriOperationSpecifics):
  """clip operation with a clip type of "box", adapter to the catalyst filter

Phactori Adapter interface to the ParaView Clip filter with a clip type
of "box". Refere to the ParaView documentation for extensive
documentation of the clip filter.

To control the BoxClipOperation you have a section in the
"operation blocks" section of the phactori control script which looks like:

::

    "operation blocks":{
      "myboxclip1":{
        "type":"boxclip",
        "center at absolute point":[1.1, -2.2, 3.5],
        "absolute extents":[5.0, 6.2, 7.3]
      }
    }

"center at absolute point" may also be "center at relative point" or
"center at element" or "center at node" or "center at data point", as it is
a "phactori use point". "absolute extents" may be "relative extents" in
which case the extents of the clip box are determined relative to the size
of the bounding box of the data.

"""

  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mCenterPtInfo = PhactoriUserPointInfo()
    self.mExtentsRelAbsFlag = [True, True, True]
    self.mExtents = [0.5, 0.5, 0.5]
    self.mRotations = [0.0, 0.0, 0.0]
    self.mKeepInside = True
    self.mCrinkleSetting = None

  def ParseParametersFromJson(self, inJson):
    if PhactoriDbg(100):
      myDebugPrint3(
        'PhactoriBoxClipOperation.ParseParametersFromJson entered\n', 100)
    self.mCenterPtInfo.UserPointInfoParseParametersFromJson(
            inJson, "center at ", "")

    absExtentsKey = 'absolute extents'
    relExtentsKey = 'relative extents'
    if absExtentsKey in inJson:
      self.mExtents = inJson[absExtentsKey]
      self.mExtentsRelAbsFlag = [False, False, False]
    elif relExtentsKey in inJson:
      self.mExtents = inJson[relExtentsKey]
      self.mExtentsRelAbsFlag = [True, True, True]
    else:
      self.mExtents = [0.5, 0.5, 0.5]
      self.mExtentsRelAbsFlag = [True, True, True]

    rotationsKey = 'rotations'
    if rotationsKey in inJson:
      self.mRotations = inJson[rotationsKey]

    if 'keep inside box' in inJson:
      self.mKeepInside = inJson['keep inside box']
    else:
      self.mKeepInside = True

    if 'cut type' in inJson:
      cutType = inJson['cut type']
      self.mCrinkleSetting = 0
      if cutType == 'smooth':
        self.mCrinkleSetting = 0
      elif cutType == 'crinkle':
        self.mCrinkleSetting = 1
      else:
        if PhactoriDbg():
          myDebugPrint3('ParseParametersFromJson warning: cut type \
            should be smooth or crinkle, but is not, using smooth\n')
        self.mCrinkleSetting = 0
    else:
      self.mCrinkleSetting = 0  #debatable default

    #self.PrintSelf()
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriBoxClipOperation.ParseParametersFromJson returning\n', 100)

  def PrintSelf(self):
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriBoxClipOperation.PrintSelf entered\nCenter Point:', 100)
    self.mCenterPtInfo.PrintSelf()
    if PhactoriDbg():
      myDebugPrint3('mExtentsRelAbsFlag: ' + str(self.mExtentsRelAbsFlag) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mExtents: ' + str(self.mExtents) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mRotations: ' + str(self.mRotations) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mKeepInside: ' + str(self.mKeepInside) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mCrinkleSetting: ' + str(self.mCrinkleSetting) + '\n')

  def CreateParaViewFilter(self, inInputFilter):
    """create the box clip box for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriBoxClipOperation.CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    if PhactoriDbg(100):
      DebugPrintCellAndPointArrayInfo("data arrays coming into PhactoriBoxClipOperation.CreateParaViewFilter\n", inInputFilter, 100)

    newParaViewFilter = Clip(Input = inInputFilter, ClipType = "Box")

    viewBoundsArray = [None]
    centerToUse = self.mCenterPtInfo.GetCurrentGeometricPointWithDisplacement(
        inInputFilter, viewBoundsArray, True)

    if PhactoriDbg():
      myDebugPrint3("  about to find extents--absolute or calced from relative\n")
    myViewBounds = viewBoundsArray[0]
    if self.mExtentsRelAbsFlag[0]:
      if myViewBounds == None:
        myViewBounds = GetGlobalDataBoundsParallel(inInputFilter)
      halfXExt = self.mExtents[0] * (myViewBounds[1] - myViewBounds[0]) * 0.5
      if halfXExt <= 1e-20:
        halfXExt = 0.001
    else:
      halfXExt = self.mExtents[0] * 0.5

    if self.mExtentsRelAbsFlag[1]:
      if myViewBounds == None:
        myViewBounds = GetGlobalDataBoundsParallel(inInputFilter)
      halfYExt = self.mExtents[1] * (myViewBounds[3] - myViewBounds[2]) * 0.5
      if halfYExt <= 1e-20:
        halfYExt = 0.001
    else:
      halfYExt = self.mExtents[1] * 0.5

    if self.mExtentsRelAbsFlag[2]:
      if myViewBounds == None:
        myViewBounds = GetGlobalDataBoundsParallel(inInputFilter)
      halfZExt = self.mExtents[2] * (myViewBounds[5] - myViewBounds[4]) * 0.5
      if halfZExt <= 1e-20:
        halfZExt = 0.001
    else:
      halfZExt = self.mExtents[2] * 0.5

    if PhactoriDbg():
      myDebugPrint3("  halfExtents: " + str(halfXExt) + ", " + str(halfYExt) + ", " + str(halfZExt) + "\n")

    #newParaViewFilter.ClipType.Position = [centerToUse[0],
    #    centerToUse[1],
    #    centerToUse[2]]
    #newParaViewFilter.ClipType.Bounds = [-halfXExt, halfXExt,
    #    -halfYExt, halfYExt, -halfZExt, halfZExt]
    newParaViewFilter.ClipType.Position = [centerToUse[0] - halfXExt,
        centerToUse[1] - halfYExt,
        centerToUse[2] - halfZExt]
    newParaViewFilter.ClipType.Length = [halfXExt * 2.0, halfYExt * 2.0, halfZExt * 2.0]
    newParaViewFilter.ClipType.Rotation = [self.mRotations[0],
        self.mRotations[1],
        self.mRotations[2]]

    newParaViewFilter.Crinkleclip = self.mCrinkleSetting
    if gParaViewCatalystVersionFlag < 50502:
      if self.mKeepInside:
        newParaViewFilter.InsideOut = 1
      else:
        newParaViewFilter.InsideOut = 0
    else:
      if self.mKeepInside:
        newParaViewFilter.Invert = 1
      else:
        newParaViewFilter.Invert = 0

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      DebugPrintCellAndPointArrayInfo("data arrays coming out of PhactoriBoxClipOperation.CreateParaViewFilter\n", newParaViewFilter, 100)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriBoxClipOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

  def MayChangeWithData(self):
    if self.mCenterPtInfo.mMayChangeWithData == True:
      return True
    if self.mExtentsRelAbsFlag[0] == True:
      return True
    if self.mExtentsRelAbsFlag[1] == True:
      return True
    if self.mExtentsRelAbsFlag[2] == True:
      return True
    return False
#phactori_combine_to_single_python_file_subpiece_end_1
