# Copyright(C) 1999-2022 National Technology & Engineering Solutions
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
class PhactoriRepresentationBlock:
  def __init__(self):
    self.mName = ""
    self.mColorVariableInfo = PhactoriVariableInfo()
    self.mColorLegendFlag = True
    #self.mColorLegendPositionAndSize = ['bottom', 1.0]
    #self.mColorLegendPositionAndSize = ['right', 1.0]
    self.mColorLegendPositionAndSize = ['bottom right', 1.0]
    self.mTimeAnnotationSettings = PhactoriAnnotationViewSettings(
        'time', 'time')
    self.mDataCubeAxesFlag = False
    self.mDataCubeAxesInfo = PhactoriDataCubeAxesInfo()
    #self.mColorByBlockFlag = True
    self.mColorByBlockFlag = True
    self.mColorByBlockExplicitlySet = False
    self.mColorBySolidColorFlag = False
    self.mSolidColor = None
    self.mOpacitySetting = 1.0
    self.mOrientationAxesFlag = True
    self.mFixedColorRange = None
    self.mUseFixedColorRange = False
    self.mColorRangeMinMaxTracker = PlotValMinMaxTrkC()
    self.mUseHighlightSubranges = False
    self.mHighlightSubranges = []
    self.mUseColorSubrange = False
    self.mColorSubrange = [0.0, 1.0]
    self.mFilenameAddon = ""
    self.mColorSettings = PhactoriColorSettings()
    self.mColorMapSettings = PhactoriColorMapSettings()
    self.mPointSize = 2.0
    self.InterpretValuesAsCategories = False
    self.CategoryAndColorJsonFile = None
    self.Annotations = None
    self.IndexedColors = None
    self.IndexedOpacities = None

  def ParseSettingsFromJson(self, inRepresentationBlockJson):
    """given a python dict (presumably from json) description of a
       representation block, parse all the representation parameters out of
       the block"""

    if PhactoriDbg(100):
      myDebugPrint3('ParseSettingsFromJson entered\n', 100)

    inJsn = inRepresentationBlockJson

    #hack test for colors start
    #inJsn['surface color'] = [1.0, 0.0, 0.0]
    #inJsn['background color'] = [0.0, 0.0, 0.0]
    #inJsn['edge color'] = [0.0, 1.0, 0.0]
    ##inJsn['axes color'] = [0.0, 1.0, 1.0]
    #inJsn['text color'] = [1.0, 1.0, 0.0]
    #hack test for colors end

    #myDebugPrint3('  json: ' + str(inJsn) + '\n')

    if 'point size' in inJsn:
      self.mPointSize = inJsn['point size']

    #color by variable scalar/vector magnitude/vector component/tensor component
    self.mColorVariableInfo.\
      ParseVariableNameAndVectorOrTensorComponent(inJsn, 'color by ')

    if self.mColorVariableInfo.mVariableName != '':
      self.mColorByBlockFlag = False
      self.mColorBySolidColorFlag = False
    elif 'color by blockid' in inJsn:
      self.mColorByBlockFlag = True
      self.mColorByBlockExplicitlySet = True
      self.mColorBySolidColorFlag = False
    elif 'color by solid color' in inJsn:
      self.mColorBySolidColorFlag = True
      self.mColorByBlockFlag = False
      self.mSolidColor = inJsn['color by solid color']

    #color map range control
    if 'color legend range' in inJsn:
      self.mFixedColorRange = inJsn['color legend range']
      self.mUseFixedColorRange = True
    else:
      self.mUseFixedColorRange = False

    #highlight subranges with solid colors
    if 'highlight subrange 1' in inJsn:
      highlightSubrangeIndex = 1
      while True:
        oneSubrangeKey = 'highlight subrange ' + str(highlightSubrangeIndex)
        if oneSubrangeKey in inJsn:
          subrangeArgs = inJsn[oneSubrangeKey]
          highlightSubrangeIndex += 1
          if len(subrangeArgs) != 5:
            if PhactoriDbg():
              myDebugPrint3("highlight subrange needs 5 values\n");
            PrintOnProcessZero("highlight subrange needs 5 values, skipping " + \
                    oneSubrangeKey + "\n");
            continue
          srmin = float(subrangeArgs[0])
          srmax = float(subrangeArgs[1])
          if srmin > srmax:
            if PhactoriDbg():
              myDebugPrint3("subrange highlight min >= max: " + \
                str(srmin) + ", " + str(srmax) + "\nskipping " + \
                oneSubrangeKey + "\n", 100)
            PrintOnProcessZero("subrange highlight min >= max: " + \
              str(srmin) + ", " + str(srmax) + "\nskipping " + \
              oneSubrangeKey + "\n")
            continue
          srColor = [float(subrangeArgs[2]), float(subrangeArgs[3]),
                  float(subrangeArgs[4])]
          if (srColor[0] < 0.0) or (srColor[0] > 1.0) or \
             (srColor[1] < 0.0) or (srColor[1] > 1.0) or \
             (srColor[2] < 0.0) or (srColor[2] > 1.0):
            srColor = [1.0, 1.0, 0.0]
            if PhactoriDbg():
              myDebugPrint3(oneSubrangeKey + ": bad color "
                "(component not 0.0-1.0), using rgb 1.0, 1.0, 0.0\n", 100)
            PrintOnProcessZero(oneSubrangeKey + ": bad color "
              "(component not 0.0-1.0), using rgb 1.0, 1.0, 0.0\n")
          self.mHighlightSubranges.append(
                  [srmin, srmax, srColor])
          self.mUseHighlightSubranges = True
        else:
          break
      if PhactoriDbg():
        myDebugPrint3("parsed highlight subranges:\n" + \
                str(self.mHighlightSubranges) + "\n", 100)

    #if 'highlight subranges' in inJsn:
    #  self.mUseHighlightSubranges = True
    #  sbrngsJsn = inJsn['highlight subranges']
    #  for oneSubrange in sbrngsJsn:
    #    if 'range' in oneSubrange:
    #      srmin = oneSubrange['range'][0]
    #      srmax = oneSubrange['range'][1]
    #      if srmin < srmax:
    #        if 'color' in oneSubrange:
    #            srColor = oneSubrange['color']
    #        else:
    #            srColor = [1.0, 1.0, 0.0]
    #        self.mHighlightSubranges.append(
    #                [srmin, srmax, srColor])
    #      else:
    #        if PhactoriDbg():
    #          myDebugPrint3("highlight min >= max: " + \
    #            str(srmin) + ", " + str(srmax) + "\nskipping subrange\n", 100)
    #        PrintOnProcessZero("highlight min >= max: " + \
    #          str(srmin) + ", " + str(srmax) + "\nskipping subrange\n")
    #    else:
    #        if PhactoriDbg():
    #          myDebugPrint3("subrange is missing 'range' key; skipping\n")
    #        PrintOnProcessZero("subrange is missing 'range' key; skipping\n")
    #  if PhactoriDbg():
    #    myDebugPrint3("parsed highlight subranges:\n" + \
    #            str(self.mHighlightSubranges) + "\n", 100)

    #additional capability: use ratio-expressed subrange of
    #range that would otherwise be used to increase concentration of
    #dynamic range of color map in range of interest
    if 'color legend subrange' in inJsn:
      self.mUseColorSubrange = True
      self.mColorSubrange = inJsn['color legend subrange']
      goodSubrange = True
      if self.mColorSubrange[0] < 0.0:
        goodSubrange = False
      if self.mColorSubrange[1] > 1.0:
        goodSubrange = False
      if self.mColorSubrange[0] > \
              self.mColorSubrange[1]:
        goodSubrange = False
      if goodSubrange == False:
        myDebugPrint3AndException(
          "ParseOneRepresentationBlockC:\n"
          "bad color legend subrange, must be 0.0 <= bottom <= top <= 1.0\n")

    #choose color map by name
    self.mColorMapSettings.ParseColorMapSettings(inJsn)

    self.InterpretValuesAsCategories = \
      getParameterFromBlock(inJsn, "interpret values as categories", False)
    if self.InterpretValuesAsCategories:
      key1 = "category and color json file"
      if key1 not in inJsn:
        myDebugPrint3AndException(
          "interpret values as categories is true but '" + key1 + \
                  "' is missing\n")
      self.CategoryAndColorJsonFile = \
        getParameterFromBlock(inJsn, key1, "dummy1.json")
    elif 'color by scalar' in inJsn:
      if inJsn['color by scalar'] == "ObjectId":
        if PhactoriDbg():
          myDebugPrint3("variable is ObjectId, using categories\n")
        self.InterpretValuesAsCategories = True
        self.CategoryAndColorJsonFile = "ObjectIdBlockColors.json"

    showSurfacesFlag = getParameterFromBlock(inJsn, 'show surfaces', True)
    showEdgesFlag = getParameterFromBlock(inJsn, 'show edges', False)
    showPointsFlag = getParameterFromBlock(inJsn, 'show points', False)

    self.mOpacitySetting = getParameterFromBlock(inJsn, 'opacity',
            self.mOpacitySetting)

    #doVolumeRenderingFlag = getParameterFromBlock(inJsn, 'volume rendering', True)
    doVolumeRenderingFlag = getParameterFromBlock(inJsn, 'volume rendering', False)
    self.mScalarOpacityUnitDistance = getParameterFromBlock(
            inJsn, 'scalar opacity unit distance', -1.0)

    #self.mScalarOpacityUnitDistance = 0.01
    if PhactoriDbg():
      myDebugPrint3("doVolumeRenderingFlag: " + \
              str(doVolumeRenderingFlag) + "\n" + \
              "self.mScalarOpacityUnitDistance: " + \
              str(self.mScalarOpacityUnitDistance) + "\n")

    self.mPresetsImportFileName = getParameterFromBlock(inJsn,
            'color and opacity presets import file', None)
    if self.mPresetsImportFileName != None:
      retval = ImportPresets(self.mPresetsImportFileName)
      if retval != True:
        myDebugPrint3AndException(
          "paraview.simple.ImportPresets failed with the file:\n" + \
          str(self.mPresetsImportFileName) + "\n")
    self.mNameOfPresetToUse = getParameterFromBlock(inJsn,
            'color and opacity preset', None)

    showBoundingBoxFlag = getParameterFromBlock(inJsn, 'show bounding box', False)

    self.mMeshRenderControl = 'Surface'
    if showBoundingBoxFlag:
      self.mMeshRenderControl = 'Outline'
      if showSurfacesFlag | showEdgesFlag | showPointsFlag:
        if PhactoriDbg():
          myDebugPrint3("  warning:  when show bounding box is true, \n" + \
              "  show surfaces, show edges, and show points should be false\n")
    elif showPointsFlag:
      self.mMeshRenderControl = 'Points'
    else:
      if showSurfacesFlag:
        if showEdgesFlag:
          self.mMeshRenderControl = 'Surface With Edges'
        else:
          self.mMeshRenderControl = 'Surface'
      else:
        if showEdgesFlag:
          self.mMeshRenderControl = 'Wireframe'
          #self.mMeshRenderControl = 'Points'
        else:
          self.mMeshRenderControl = 'Outline'

    if doVolumeRenderingFlag:
      if PhactoriDbg(100):
        myDebugPrint3('doing volume rendering\n', 100)
      self.mMeshRenderControl = 'Volume'
      self.mDoingVolumeRendering = True
    else:
      self.mDoingVolumeRendering = False

    #color legend on/off
    self.mColorLegendFlag = getParameterFromBlock(inJsn,
      'show color legend', self.mColorLegendFlag)

    self.mColorLegendPositionAndSize = \
      getParameterFromBlock(inJsn, 'color legend position',
          self.mColorLegendPositionAndSize)

    if self.mColorLegendFlag == True:
      if self.mColorVariableInfo.mVariableName == '':
        self.mColorLegendFlag = False

    self.mColorRangeMinMaxTracker.PlotValMinMaxTrkCParseJson(
            inJsn, "color legend ")

    self.mTimeAnnotationSettings.ParseAvsFromJson(inJsn)

    self.mDataCubeAxesFlag = getParameterFromBlock(inJsn,
      'show axes', self.mDataCubeAxesFlag)

    self.mDataCubeAxesInfo.DcaiParseParametersFromJson(inJsn)

    self.mOrientationAxesFlag = getParameterFromBlock(inJsn,
      'show orientation axes', self.mOrientationAxesFlag)

    if "image name addon" in inJsn:
      self.mFilenameAddon = inJsn["image name addon"]

    self.mColorSettings.ParseColorSettingsFromJson(inJsn)

    if PhactoriDbg(100):
      myDebugPrint3('ParseOneRepresentationBlockC returning\n', 100)

  def SetFromRestartInfo(self, inJson):
    """given a map (json format), use the info in the map to set the
       representation state--this reads the info created in 
       GetRestartInfo"""

    if 'mColorRangeMinMaxTracker' not in inJson:
      if PhactoriDbg():
        myDebugPrint3("PhactoriRepresentationBlock::" + \
            "SetFromRestartInfo: no mColorRangeMinMaxTracker, return\n")
      return

    if PhactoriDbg():
      myDebugPrint3("Representation::SetFromRestartInfo\n" +
        "currently set to do something\n" +
        "before tracker:\n" + self.mColorRangeMinMaxTracker.SelfToStr())

    jsonItem = inJson['mColorRangeMinMaxTracker']
    self.mColorRangeMinMaxTracker.SetFromRestartInfo(jsonItem)

    if PhactoriDbg():
      myDebugPrint3("after tracker:\n" +
          self.mColorRangeMinMaxTracker.SelfToStr())

  def GetRestartInfo(self):
    """construct, in python map/json format, the information from this
       representation instance which contains the information
       which would be needed to restore the representation to the proper
       state after a simulation restart, particularly color range tracking
       information.  Return the restart info map/json"""
    newRestartInfoJson = {}
    newRestartInfoJson['mColorRangeMinMaxTracker'] = \
      self.mColorRangeMinMaxTracker.GetRestartInfo()
    return newRestartInfoJson

  def CalculateDefaultScalarOpacityUnitDistance(self, inPhactoriOperation):
    """given a phactori operation assumed to have an updated pipeline,
       calculate a default scalar opacity unit distance the same as paraview,
       with 0.05 * diagonal length of data bounding box"""
    bnds = GetGlobalDataBoundsParallel(inPhactoriOperation.mParaViewFilter)
    xxlen = bnds[1] - bnds[0]
    yylen = bnds[3] - bnds[2]
    zzlen = bnds[5] - bnds[4]
    self.mScalarOpacityUnitDistance = \
            0.05 * math.sqrt(xxlen*xxlen + yylen*yylen + zzlen*zzlen)

  def SetUpForInterpretValuesAsCateories(self, variableLookUpTable):
    if PhactoriDbg():
      myDebugPrint3("SetUpForInterpretValuesAsCateories entered\n")
    if(self.Annotations == None):
      jsonCategories = ReadAndMpiBroadcastJsonFile(self.CategoryAndColorJsonFile)
      jsonCategories = convertJsonUnicodeToStrings(jsonCategories)
      self.Annotations = jsonCategories[0]["Annotations"]
      self.IndexedColors = jsonCategories[0]["IndexedColors"]
      self.IndexedOpacities = []
      numCategories = len(self.Annotations) // 2
      for ii in range(0, numCategories):
        self.IndexedOpacities.append(1.0)
      if PhactoriDbg():
        myDebugPrint3("self.Annotations:\n" + str(self.Annotations) + "\n")
        myDebugPrint3("self.IndexedColors:\n" + str(self.IndexedColors) + "\n")
        myDebugPrint3("self.IndexedOpacities:\n" + str(self.IndexedOpacities) + "\n")
    variableLookUpTable.InterpretValuesAsCategories = 1
    variableLookUpTable.AnnotationsInitialized = 1
    variableLookUpTable.Annotations = self.Annotations
    variableLookUpTable.IndexedColors = self.IndexedColors
    variableLookUpTable.IndexedOpacities = self.IndexedOpacities
    if PhactoriDbg():
      myDebugPrint3("SetUpForInterpretValuesAsCateories returning\n")
#phactori_combine_to_single_python_file_subpiece_end_1

