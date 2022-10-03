# Copyright(C) 1999-2021 National Technology & Engineering Solutions
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
class PhactoriGlyphOperation(PhactoriOperationSpecifics):
  """manages Glyph filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.scaleFactor = None
    self.glyphType = "Arrow"
    self.ArrowTipResolution = None
    self.ArrowTipRadius = None
    self.ArrowTipLength = None
    self.ArrowShaftResolution = None
    self.ArrowShaftRadius = None
    self.scaleArray = None
    self.orientationArray = None

  def ParseParametersFromJson(self, inJson):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGlyphOperation.ParseParametersFromJson "
          "entered\n", 100)

    key1 = "scale factor"
    if key1 in inJson:
      self.scaleFactor = inJson[key1]
      if self.scaleFactor <= 0.0:
        myDebugPrint3AndException(
          "PhactoriGlyphOperation.ParseParametersFromJson:\nbad " + \
          key1 + "\n")
        exit(-1)

    glyphSet = {}
    glyphSet["arrow"] = "Arrow"
    glyphSet["box"] = "Box"
    key1 = "glyph type"
    if key1 in inJson:
      glyphkey1 = inJson[key1]
      if glyphkey1 not in glyphSet:
        myDebugPrint3AndException(
          "PhactoriGlyphOperation.ParseParametersFromJson:\nbad " + \
          key1 + "\n")
        exit(-1)
      self.glyphType = glyphSet[glyphkey1]

    key1 = "scale array"
    if key1 in inJson:
      self.scaleArray = inJson[key1]

    key1 = "orientation array"
    if key1 in inJson:
      self.orientationArray = inJson[key1]

    if PhactoriDbg(100):
      myDebugPrint3("self.scaleArray " + str(self.scaleArray) + "\n")
      myDebugPrint3("self.orientationArray " + str(self.orientationArray) + "\n")

    if self.glyphType == "Arrow":
      key1 = "tip resolution"
      if key1 in inJson:
        self.ArrowTipResolution = inJson[key1]
        if self.ArrowTipResolution <= 0:
          myDebugPrint3AndException(
            "PhactoriGlyphOperation.ParseParametersFromJson:\nbad " + \
            key1 + "\n")
          exit(-1)
      key1 = "tip length"
      if key1 in inJson:
        self.ArrowTipLength = inJson[key1]
        if self.ArrowTipLength < 0.0:
          myDebugPrint3AndException(
            "PhactoriGlyphOperation.ParseParametersFromJson:\nbad " + \
            key1 + "\n")
          exit(-1)
      key1 = "tip radius"
      if key1 in inJson:
        self.ArrowTipRadius = inJson[key1]
        if self.ArrowTipRadius <= 0.0:
          myDebugPrint3AndException(
            "PhactoriGlyphOperation.ParseParametersFromJson:\nbad " + \
            key1 + "\n")
          exit(-1)
      key1 = "shaft resolution"
      if key1 in inJson:
        self.ArrowShaftResolution = inJson[key1]
        if self.ArrowShaftResolution <= 0:
          myDebugPrint3AndException(
            "PhactoriGlyphOperation.ParseParametersFromJson:\nbad " + \
            key1 + "\n")
          exit(-1)
      key1 = "shaft radius"
      if key1 in inJson:
        self.ArrowShaftRadius = inJson[key1]
        if self.ArrowShaftRadius <= 0.0:
          myDebugPrint3AndException(
            "PhactoriGlyphOperation.ParseParametersFromJson:\nbad " + \
            key1 + "\n")
          exit(-1)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGlyphOperation.ParseParametersFromJson "
          "returning\n", 100)

  def CreateParaViewFilter(self, inInputFilter):
    """create the Glyph filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGlyphOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    if PhactoriDbg(100):
      myDebugPrint3("self.glyphType: " + str(self.glyphType) + "\n")

    newParaViewFilter = Glyph(Input=inInputFilter, GlyphType=self.glyphType)
    #newParaViewFilter.OrientationArray = ['POINTS', 'No orientation array']
    #newParaViewFilter.ScaleArray = ['POINTS', 'No scale array']
    #newParaViewFilter.ScaleFactor = 0.1753620492294431
    newParaViewFilter.GlyphTransform = 'Transform2'
    if self.scaleArray != None:
      newParaViewFilter.ScaleArray = self.scaleArray
    else:
      newParaViewFilter.ScaleArray = ['POINTS', 'velocity_xyz']
    if self.orientationArray != None:
      newParaViewFilter.OrientationArray = self.orientationArray
    else:
      newParaViewFilter.OrientationArray = ['POINTS', 'velocity_xyz']
    if PhactoriDbg(100):
      myDebugPrint3("newParaViewFilter.ScaleArray " + str(newParaViewFilter.ScaleArray) + "\n")
      myDebugPrint3("newParaViewFilter.OrientationArray " + str(newParaViewFilter.OrientationArray) + "\n")
    #newParaViewFilter.ScaleFactor = 1e-05
    newParaViewFilter.GlyphMode = 'All Points'
    if self.scaleFactor != None:
      newParaViewFilter.ScaleFactor = self.scaleFactor
    if self.ArrowTipResolution != None:
      newParaViewFilter.GlyphType.TipResolution = self.ArrowTipResolution
    if self.ArrowShaftResolution != None:
      newParaViewFilter.GlyphType.ShaftResolution = self.ArrowShaftResolution
    if self.ArrowTipRadius != None:
      newParaViewFilter.GlyphType.TipRadius = self.ArrowTipRadius
    if self.ArrowShaftRadius != None:
      newParaViewFilter.GlyphType.ShaftRadius = self.ArrowShaftRadius
    if self.ArrowTipLength != None:
      newParaViewFilter.GlyphType.TipLength = self.ArrowTipLength

    # init the 'Arrow' selected for 'GlyphType'
    #newParaViewFilter.GlyphType.TipResolution = 6
    #newParaViewFilter.GlyphType.TipRadius = 0.05
    #newParaViewFilter.GlyphType.TipLength = 0.3
    #newParaViewFilter.GlyphType.ShaftResolution = 6
    #newParaViewFilter.GlyphType.ShaftRadius = 0.05

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)
    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGlyphOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
