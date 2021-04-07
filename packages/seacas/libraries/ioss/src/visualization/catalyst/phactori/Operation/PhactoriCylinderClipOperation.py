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
class PhactoriCylinderClipOperation(PhactoriOperationSpecifics):
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mCenterPtInfo = PhactoriUserPointInfo()
    self.mAxis = [0.0, 1.0, 0.0]
    self.mRadius = 0.5
    self.mRadiusRelAbsFlag = True
    self.mKeepInside = True
    self.mCrinkleSetting = None

  def ParseParametersFromJson(self, inJson):
    if PhactoriDbg(100):
      myDebugPrint3(
        'PhactoriCylinderClipOperation.ParseParametersFromJson entered\n', 100)
    self.mCenterPtInfo.UserPointInfoParseParametersFromJson(
            inJson, "center at ", "")

    absRadiusKey = 'absolute radius'
    relRadiusKey = 'relative radius'
    if absRadiusKey in inJson:
      self.mRadius = inJson[absRadiusKey]
      self.mRadiusRelAbsFlag = False
    elif relRadiusKey in inJson:
      self.mRadius = inJson[relRadiusKey]
      self.mRadiusRelAbsFlag = True
    else:
      self.mRadius = 0.5
      self.mRadiusRelAbsFlag = True

    axisKey = 'axis'
    if axisKey in inJson:
      self.mAxis = inJson[axisKey]

    if 'keep inside cylinder' in inJson:
      self.mKeepInside = inJson['keep inside cylinder']
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
      myDebugPrint3(
        'PhactoriCylinderClipOperation.ParseParametersFromJson returning\n',
        100)

  def PrintSelf(self):
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriCylinderClipOperation.PrintSelf entered\n" + \
        "Center Point:', 100)
    self.mCenterPtInfo.PrintSelf()
    if PhactoriDbg():
      myDebugPrint3('mRadiusRelAbsFlag: ' + str(self.mRadiusRelAbsFlag) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mRadius: ' + str(self.mRadius) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mAxis: ' + str(self.mAxis) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mKeepInside: ' + str(self.mKeepInside) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mCrinkleSetting: ' + str(self.mCrinkleSetting) + '\n')

  def CreateParaViewFilter(self, inInputFilter):
    """create the cylinder clip filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriCylinderClipOperation.CreateParaViewFilter ' +\
        'entered\n', 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    newParaViewFilter = Clip(Input = inInputFilter, ClipType = "Cylinder")
    ## Properties modified on clip4.ClipType
    #clip4.ClipType.Center = [10.93421983718872, -0.8890213072299957, 0.0]
    #clip4.ClipType.Axis = [1.0, 0.0, 0.0]
    #clip4.ClipType.Radius = 0.75

    newParaViewFilter.ClipType.Axis = self.mAxis

    viewBoundsArray = [None]
    centerToUse = self.mCenterPtInfo.GetCurrentGeometricPointWithDisplacement(
        inInputFilter, viewBoundsArray, True)

    newParaViewFilter.ClipType.Center = [centerToUse[0],
        centerToUse[1],
        centerToUse[2]]

    if PhactoriDbg():
      myDebugPrint3("  about to find radius--absolute or calced from relative\n")
    myViewBounds = viewBoundsArray[0]
    if self.mRadiusRelAbsFlag:
      if myViewBounds == None:
        myViewBounds = GetGlobalDataBoundsParallel(inInputFilter)
      dxx = myViewBounds[1] - myViewBounds[0]
      dyy = myViewBounds[3] - myViewBounds[2]
      dzz = myViewBounds[5] - myViewBounds[4]
      bounds_diagonal = math.sqrt(dxx*dxx + dyy*dyy + dzz*dzz)
      newParaViewFilter.ClipType.Radius = self.mRadius * bounds_diagonal
    else:
      newParaViewFilter.ClipType.Radius = self.mRadius

    newParaViewFilter.Crinkleclip = self.mCrinkleSetting
    if self.mKeepInside:
      newParaViewFilter.InsideOut = 1
    else:
      newParaViewFilter.InsideOut = 0

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3(
        'PhactoriCylinderClipOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

  def MayChangeWithData(self):
    if self.mCenterPtInfo.mMayChangeWithData == True:
      return True
    if self.mRadiusRelAbsFlag == True:
      return True
    return False
#phactori_combine_to_single_python_file_subpiece_end_1
