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
from .PhactoriOperationBlock import *
from .PhactoriVectorLibrary import *
from paraview.simple import *
from .PhactoriSegment import *
from .PhactoriMpiUtilities import *
from .PhactoriParaviewMultiBlockRecursion import *
from .PhactoriSampledCellInfo import *
import vtk
import json

#phactori_combine_to_single_python_file_subpiece_begin_1

class GetCellsClosestToPointsInBlockRecursionParams:
  def __init__(self):
    self.testPointList = []
    self.distSqrdList = []
    self.closestList = []

  def InitializeWithPointList(self, inTestPointList):
    self.testPointList = inTestPointList
    numTestPoints = len(inTestPointList)
    for ii in range(0, numTestPoints):
      self.distSqrdList.append(sys.float_info.max)
      self.closestList.append(None)

class SetMaskValueForCellsNearSegmentsRecursionParams:
  def __init__(self):
    self.segmentList = None
    self.maskCellVariableName = None
    self.globalCellIdName = None
    self.CellDataArrayName = "V"
    self.markedCellSet = None
    self.maskTestDistanceSquared = None
    self.projectionAxis = -1
    self.leafVisitCount = 0

  def Initialize(self, segmentList, cellDataArrayName, maskCellVariableName, globalCellIdName, maskTestDistanceSquared, projectionAxis):
    self.CellDataArrayName = cellDataArrayName
    self.maskCellVariableName = maskCellVariableName
    self.globalCellIdName = globalCellIdName
    self.segmentList = segmentList
    self.maskTestDistanceSquared = maskTestDistanceSquared
    self.markedCellSet = []
    self.projectionAxis = projectionAxis
    self.leafVisitCount = 0
    if PhactoriDbg(100):
      myDebugPrint3("SetMaskValueForCellsNearSegmentsRecursionParams::Initialize\n" +\
        "self.segmentList:\n" + str(self.segmentList) + "\n")

    for oneSegment in self.segmentList:
      cellInfoListForThisSegment = []
      self.markedCellSet.append(cellInfoListForThisSegment)

class CellsListClosestToPointsListInBlockStructuredRecursionParams:
  def __init__(self):
    self.testPointList = []
    self.distSqrdList = []
    self.closestCellList = []
    self.leafVisitCount = 0

  def InitializeWithPointList(self, inTestPointList):
    self.testPointList = inTestPointList
    numTestPoints = len(inTestPointList)
    self.leafVisitCount = 0
    for ii in range(0, numTestPoints):
      self.distSqrdList.append(sys.float_info.max)
      newCellInfo = PhactoriSampledCellInfo()
      self.closestCellList.append(newCellInfo)

class GatherStructuredCellsFromSeedCellsRecursionParams:
  def __init__(self):
    self.leafVisitCount = 0
    self.dataArrayName = "noname"
    self.SeedCells = None
    self.CellsPerSeedCell = None

  def SetUpForRecursion(self, creatingSamplerInstance):
    self.SeedCells = creatingSamplerInstance.StructuredNearbyCellPointList
    self.dataArrayName = creatingSamplerInstance.CellDataArrayName
    self.CellsPerSeedCell = []
    for ii in range(0, len(self.SeedCells)):
      self.CellsPerSeedCell.append([])

class PhactoriSegmentCellSampler3(PhactoriOperationSpecifics):
  """Filter/operation which reads in a .json file which is a list of 3d
     coordinates pairs (representating segments, and then finds all cells
     in the input data that cross that segment and sets a masking value
     into a cell data array for those cells"""
  def __init__(self):
    self.myCopyOfInputFilter = None

    #[[[x1,y1],[x2.y2]],[[x3,y3],[x4.y4]],...]
    self.JsonSegmentListFileName = "PhactoriSegmentCellSampler3.json"

    self.SegmentList = None
    #self.MaskTestDistanceSquared = 0.00001
    self.MaskTestDistanceSquared = 1.6e-07

    self.SegmentDefinitionFormat = "two geometric points"

    self.PerSegmentMarkedCellInfoList = []

    #indicates whether we should conduct our segment-picking cell calculations
    #by projecting onto an ordinal plane.  -1 is no projection, 0, 1, or 2 are
    #the x, y, or z axes. This is of particular use if you are doing a slice
    #before selecting elements.  The slice does not need to be axis aligned,
    #just choose the axis plane closest to parallel to the slice plane.
    self.ProjectionAxis = -1

    #as of 2019Aug13, we are having the programmable filter core dump with
    #structured data: this allows us to turn off the programmable filter to
    #do the other operations
    self.DoProgrammableFilterToAddMaskVariable = False
    self.mProgFilterMaskSetter = None
    self.mProgFilterString = None
    self.ProgrammableFilterOutputCellVariableName = "mask1"

    self.StructuredNearbyCellPointList = None
    self.NearbyGeometryPointList = None

    #"structured" or "unstructured"
    self.collectionMethod = "segments identify cells"

    self.SampledCellsOutputFilename = "PhactoriSegmentCellSampler3_sampled_cells_"
    self.SampledCellsOutputDirectory = "PhactoriSegmentCellSampler3_Output"
    self.NumCounterDigits = 5
    self.SampledCellsOutputFilenameExtension = ".json"

    self.CellDataArrayName = "V"
    self.CellDataArrayTupleSize = 8

    #controls whether or not we output a single json file with a list of all
    #the selected cells (using MPI to combine to one process)
    self.WriteCellListToJsonFile = True

  def ValidateJsonStructuredSeedCellList(self, inJson):
    numPoints = len(inJson)
    if numPoints < 1:
      myDebugPrint3AndException(
          "PhactoriSegmentCellSampler3::ValidateJsonStructuredSeedCellList\n"
          "list must have at least one seed cell\n")

    for ptNdx, onePtJson in enumerate(inJson):
      key1 = "geometric seed point"
      if key1 not in onePtJson:
        myDebugPrint3AndException("point " + str(ptNdx) + " missing '" + key1 + "'\n" + str(inJson) + "\n")
      geompt = onePtJson[key1]
      if(len(geompt) != 3):
        myDebugPrint3AndException("point " + str(ptNdx) + " bad json format for xyz point:\n" + str(inJson) + "\n")
      key2 = "collection axis"
      if key2 in onePtJson:
        axisVal = onePtJson[key2]
        if (axisVal != "i") and (axisVal != "j") and (axisVal != "k") and \
           (axisVal != "ij") and (axisVal != "ik") and (axisVal != "jk"):
          myDebugPrint3AndException(
            "point " + str(ptNdx) + " bad json format (2):\n"\
            "needs 'i', 'j', 'k', 'ij', 'ik', or 'jk':\n" + str(inJson) + "\n")

    return True


  def ValidateJsonSegmentList2(self, testlist):
    numSegments = len(testlist)
    if numSegments < 1:
      myDebugPrint3AndException(
          "PhactoriSegmentCellSampler3::ValidateJsonSegmentList2\n"
          "list must have at least one segment\n")

    for segNdx in range(0,numSegments):
      oneSegment = testlist[segNdx]
      if len(oneSegment) != 2:
        errStr = "PhactoriSegmentCellSampler3:ValidateJsonSegmentList2\n" \
          "segment with index " + str(segNdx) + " does not have a point and a length\n"
        myDebugPrint3AndException(errStr)
      pt1 = oneSegment[0]
      if len(pt1) != 3:
        errStr = "PhactoriSegmentCellSampler3:ValidateJsonSegmentList2\n" \
          "segment with index " + str(segNdx) + " point does not have 3 components\n"
        myDebugPrint3AndException(errStr)

    return True

  def ValidateJsonSegmentList3(self, testlist):
    numSegments = len(testlist)
    if numSegments < 1:
      myDebugPrint3AndException(
          "PhactoriSegmentCellSampler3::ValidateJsonSegmentList3\n"
          "list must have at least one segment\n")

    for segNdx in range(0,numSegments):
      oneSegment = testlist[segNdx]
      if len(oneSegment) != 3:
        errStr = "PhactoriSegmentCellSampler3:ValidateJsonSegmentList3\n" \
          "segment with index " + str(segNdx) + " does not have a point, direction, and length\n"
        myDebugPrint3AndException(errStr)
      pt1 = oneSegment[0]
      dir1 = oneSegment[1]
      length1 = oneSegment[2]
      if len(pt1) != 3:
        errStr = "PhactoriSegmentCellSampler3:ValidateJsonSegmentList3\n" \
          "segment with index " + str(segNdx) + " point does not have 3 components\n"
        myDebugPrint3AndException(errStr)
      if len(dir1) != 3:
        errStr = "PhactoriSegmentCellSampler3:ValidateJsonSegmentList3\n" \
          "segment with index " + str(segNdx) + " direction does not have 3 components\n"
        myDebugPrint3AndException(errStr)

    return True

  def ValidateJsonSegmentList(self, testlist):
    numSegments = len(testlist)
    if numSegments < 1:
      myDebugPrint3AndException(
          "PhactoriSegmentCellSampler3::ValidateJsonSegmentList\n"
          "list must have at least one segment\n")

    for segNdx in range(0,numSegments):
      oneSegment = testlist[segNdx]
      if len(oneSegment) != 2:
        errStr = "PhactoriSegmentCellSampler3:ValidateJsonSegmentList\n" \
          "segment with index " + str(segNdx) + " does not have two points\n"
        myDebugPrint3AndException(errStr)
      pt1 = oneSegment[0]
      if len(pt1) != 3:
        errStr = "PhactoriSegmentCellSampler3:ValidateJsonSegmentList\n" \
          "segment with index " + str(segNdx) + " point 1 does not have 3 components\n"
        myDebugPrint3AndException(errStr)
      pt2 = oneSegment[1]
      if len(pt1) != 3:
        errStr = "PhactoriSegmentCellSampler3:ValidateJsonSegmentList\n" \
          "segment with index " + str(segNdx) + " point 2 does not have 3 components\n"
        myDebugPrint3AndException(errStr)

    return True

  def ParseParametersFromJson(self, inJson):
    keyval6 = "collection method"
    if keyval6 in inJson:
      self.collectionMethod = inJson[keyval6]
      if (self.collectionMethod != "segments identify cells") and \
         (self.collectionMethod != "seed cell with structured grid"):
        myDebugPrint3AndException("'collection method' must be one of:\n" \
          "'segments identify cells' (default)\n" \
          "'seed cell with structured grid'\n")

    if "filename" in inJson:
      self.JsonListFileName = inJson['filename']
    else:
      myDebugPrint3AndException(
          "PhactoriSegmentCellSampler3:ParseParametersFromJson\n"
          "Error:  must have 'filename' key\n")

    keyval2 = "cell center to segment test distance"
    if keyval2 not in inJson:
      myDebugPrint3AndException(
          "PhactoriSegmentCellSampler3:ParseParametersFromJson\n"
          "Error:  must have '" + keyval2 + "' key\n")
    testDist = inJson[keyval2]
    self.MaskTestDistanceSquared = testDist * testDist

    if self.collectionMethod == "seed cell with structured grid":
      keyval3 = "segment definition format"
      if keyval3 in inJson:
        myDebugPrint3AndException("with collection method '" + str(self.collectionMethod) + "'\n" \
          "you cannot have 'segment definition format' key\n")
    else:
      keyval3 = "segment definition format"
      if keyval3 in inJson:
        self.SegmentDefinitionFormat = inJson[keyval3]
        if (self.SegmentDefinitionFormat != "two geometric points") and \
           (self.SegmentDefinitionFormat != "geometric point and nearest cell") and \
           (self.SegmentDefinitionFormat != "structured cell geometric point and length") and \
           (self.SegmentDefinitionFormat != "geometric point and nearest cell and direction"):
            myDebugPrint3AndException("PhactoriSegmentCellSampler3::ParseParametersFromJson\n"
              "bad " + keyval3 + "\n")

    keyval4 = "projection axis"
    if keyval4 in inJson:
      val4 = inJson[keyval4]
      if val4 == "none":
        self.ProjectionAxis = -1
      elif val4 == "x":
        self.ProjectionAxis = 0
      elif val4 == "y":
        self.ProjectionAxis = 1
      elif val4 == "z":
        self.ProjectionAxis = 2
      elif int(val4) == 1:
        self.ProjectionAxis = 0
      elif int(val4) == 2:
        self.ProjectionAxis = 1
      elif int(val4) == 3:
        self.ProjectionAxis = 2

    keyval5 = "do programmable filter"
    if keyval5 in inJson:
      self.DoProgrammableFilterToAddMaskVariable = inJson[keyval5]

    keyval7 = "sampled cells output filename"
    if keyval7 in inJson:
      self.SampledCellsOutputFilename = inJson[keyval7]

    keyval7b = "sampled cells output directory"
    if keyval7b in inJson:
      self.SampledCellsOutputDirectory = inJson[keyval7b]
    if self.SampledCellsOutputDirectory != ".":
      CreateDirectoryFromProcessZero(self.SampledCellsOutputDirectory)

    keyval7c = "sampled cells output filename digit count"
    if keyval7c in inJson:
      self.NumCounterDigits = inJson[keyval7c]

    keyval7d = "sampled cells output filename extension"
    if keyval7d in inJson:
      self.SampledCellsOutputFilenameExtension = inJson[keyval7d]

    keyval8 = "cell data array name"
    if keyval8 in inJson:
      self.CellDataArrayName = inJson[keyval8]

    keyval9 = "cell data array tuple size"
    if keyval9 in inJson:
      self.CellDataArrayTupleSize = inJson[keyval9]

    keyval10 = "programmable filter output cell variable name"
    if keyval10 in inJson:
      self.ProgrammableFilterOutputCellVariableName = inJson[keyval10]

    keyval11 = "write cell list to json file"
    if keyval11 in inJson:
      self.WriteCellListToJsonFile = inJson[keyval11]

  def SetStructuredNearbyCellPointList(self, nearbyGeometryPoints, structuredNearbyCellPointList):
        self.NearbyGeometryPointList = nearbyGeometryPoints
        self.StructuredNearbyCellPointList = structuredNearbyCellPointList
        #self.SegmentList = []
        if PhactoriDbg(100):
          myDebugPrint3("(structured) format: " + str(self.SegmentDefinitionFormat) + "\n")
        for ii in range(0,len(self.StructuredNearbyCellPointList)):
          seedPoint = nearbyGeometryPoints[ii]
          segmentPointA = self.StructuredNearbyCellPointList[ii].cellTestPoint
          SeedToCell = vecFromAToB(seedPoint, segmentPointA)
          #SeedToCellNorm = vecNormalize(SeedToCell)
          SeedToCellNorm = vecNormalizeWithSmallCheck(SeedToCell, 0.0, [0.0, 1.0, 0.0])
          #segmentPointB = vecMultiplyAdd(segmentPointA, SeedToCellNorm, segmentLengths[ii])
          segmentPointB = vecMultiplyAdd(segmentPointA, SeedToCellNorm, 1.0)
          if PhactoriDbg(100):
            myDebugPrint3("segment " + str(ii) + "\n", 100)
            myDebugPrint3("seedPoint: " + str(seedPoint) + "\n")
            myDebugPrint3("closest cell point: " + str(segmentPointA) + "\n")
            myDebugPrint3("SeedToCellNorm: " + str(SeedToCellNorm) + "\n")
            myDebugPrint3("segmentPointB: " + str(segmentPointB) + "\n")
          newSegment = PhactoriSegment()
          newSegment.SetPoints(segmentPointA, segmentPointB)
          #self.SegmentList.append(newSegment)

  def CreateInternalStructuredSeedCellListFromJson(self, inJson):
      self.ValidateJsonStructuredSeedCellList(inJson)
      nearbyGeometryPoints = []
      collectionAxes = []
      for oneSeedCellJson in inJson:
        geometryPoint = oneSeedCellJson["geometric seed point"]
        nearbyGeometryPoints.append(list(geometryPoint))
        if "collection axis" in oneSeedCellJson:
          oneAxis = oneSeedCellJson["collection axis"]
          if oneAxis == "i":
            axisInt = 0
          elif oneAxis == "j":
            axisInt = 1
          elif oneAxis == "k":
            axisInt = 2
          elif oneAxis == "ij":
            axisInt = 3
          elif oneAxis == "ik":
            axisInt = 4
          elif oneAxis == "jk":
            axisInt = 5
        else:
          axisInt = 2
        collectionAxes.append(axisInt)
      structuredNearbyCellPointList = self.GetListOfCellTestPointsNearestListOfPointsStructured(nearbyGeometryPoints)
      for ii, oneCell in enumerate(structuredNearbyCellPointList):
        oneCell.SetCollectionAxis(collectionAxes[ii])
      self.SetStructuredNearbyCellPointList(nearbyGeometryPoints, structuredNearbyCellPointList)

  def CreateInternalSegmentListFromJson(self, segmentListJson):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.ExportOperationData entered\n", 100)

    if self.SegmentDefinitionFormat == "two geometric points":
      self.ValidateJsonSegmentList(segmentListJson)
      self.SegmentList = []
      self.NearbyGeometryPointList = []
      if PhactoriDbg(100):
        myDebugPrint3("format: " + str(self.SegmentDefinitionFormat) + "\n")
      for item in segmentListJson:
        newSegment = PhactoriSegment()
        newSegment.SetPoints(item[0], item[1])
        self.SegmentList.append(newSegment)
        self.NearbyGeometryPointList.append(item[0])
    elif self.SegmentDefinitionFormat == "geometric point and nearest cell":
      self.ValidateJsonSegmentList2(segmentListJson)
      nearbyGeometryPoints = []
      segmentLengths = []
      for item in segmentListJson:
        nearbyGeometryPoints.append(item[0])
        segmentLengths.append(item[1])
      self.NearbyGeometryPointList = nearbyGeometryPoints

      nearbyCellPoints = self.GetListOfCellTestPointsNearestListOfPoints(nearbyGeometryPoints)
      #use some vector math to construct the segments from the seed points
      #and the nearest cell to each seed point
      self.SegmentList = []
      if PhactoriDbg(100):
        myDebugPrint3("(unstructured) format: " + str(self.SegmentDefinitionFormat) + "\n")
      for ii in range(0,len(nearbyGeometryPoints)):
        seedPoint = nearbyGeometryPoints[ii]
        segmentPointA = nearbyCellPoints[ii]
        SeedToCell = vecFromAToB(seedPoint, segmentPointA)
        SeedToCellNorm = vecNormalize(SeedToCell)
        segmentPointB = vecMultiplyAdd(segmentPointA, SeedToCellNorm, segmentLengths[ii])
        if PhactoriDbg(100):
          myDebugPrint3("segment " + str(ii) + "\n", 100)
          myDebugPrint3("seedPoint: " + str(seedPoint) + "\n")
          myDebugPrint3("closest cell point: " + str(segmentPointA) + "\n")
          myDebugPrint3("SeedToCellNorm: " + str(SeedToCellNorm) + "\n")
          myDebugPrint3("segmentPointB: " + str(segmentPointB) + "\n")
        newSegment = PhactoriSegment()
        newSegment.SetPoints(segmentPointA, segmentPointB)
        self.SegmentList.append(newSegment)
    elif self.SegmentDefinitionFormat == "geometric point and nearest cell and direction":
      self.ValidateJsonSegmentList3(segmentListJson)
      nearbyGeometryPoints = []
      segmentLengths = []
      directionVecs = []
      for item in segmentListJson:
        nearbyGeometryPoints.append(item[0])
        directionVecs.append(item[1])
        segmentLengths.append(item[2])
      self.NearbyGeometryPointList = nearbyGeometryPoints

      nearbyCellPoints = self.GetListOfCellTestPointsNearestListOfPoints(nearbyGeometryPoints)
      self.SegmentList = []
      if PhactoriDbg(100):
        myDebugPrint3("format: " + str(self.SegmentDefinitionFormat) + "\n")
      for ii in range(0,len(nearbyGeometryPoints)):
        segmentPointA = nearbyCellPoints[ii]
        directionVecNorm = vecNormalize(directionVecs[ii])
        segmentPointB = vecMultiplyAdd(segmentPointA, directionVecNorm, segmentLengths[ii])
        if PhactoriDbg(100):
          myDebugPrint3("seedPoint: " + str(nearbyGeometryPoints[ii]) + "\n")
          myDebugPrint3("segmentPointA (nearest cell point): " + str(seedPoint) + "\n")
          myDebugPrint3("directionVecNorm: " + str(directionVecNorm) + "\n")
          myDebugPrint3("segmentPointB: " + str(segmentPointB) + "\n")
        newSegment = PhactoriSegment()
        newSegment.SetPoints(segmentPointA, segmentPointB)
        self.SegmentList.append(newSegment)
    else:
      myDebugPrint3AndException("bad self.SegmentDefinitionFormat\n")

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.ExportOperationData returning\n", 100)

  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.CreateParaViewFilter entered\n", 100)

    savedActiveSource = GetActiveSource()

    segmentListJson = ReadAndMpiBroadcastJsonFile(self.JsonListFileName)

    self.myCopyOfInputFilter = inInputFilter

    UpdatePipelineWithCurrentTimeArgument(self.myCopyOfInputFilter)

    if self.collectionMethod == "seed cell with structured grid":
      self.CreateInternalStructuredSeedCellListFromJson(segmentListJson)
      if PhactoriDbg(100):
        myDebugPrint3("self.StructuredNearbyCellPointList:\n", 100)
        for ii, oneCell in enumerate(self.StructuredNearbyCellPointList):
          outStr = oneCell.ToStr()
          myDebugPrint3("cell " + str(ii) + ":\n" + outStr + \
          "from geometry point: " + \
          str(self.NearbyGeometryPointList[ii]) + "\n")
    else:
      self.CreateInternalSegmentListFromJson(segmentListJson)
      if PhactoriDbg(100):
        myDebugPrint3("self.SegmentList:\n", 100)
        for oneSeg in self.SegmentList:
          myDebugPrint3(str(oneSeg.ptA) + ", " + str(oneSeg.ptB) + "\n")

    if self.collectionMethod == "seed cell with structured grid":
      perSegmentStructuredCellList = self.GatherStructuredCellsFromSeedCells()
      self.PerSegmentMarkedCellInfoList = perSegmentStructuredCellList
    else:
      perSegmentCellList = self.SetMaskValueForCellsNearSegments(
        self.myCopyOfInputFilter, self.SegmentList, "blargmask1", "Ids",
        self.MaskTestDistanceSquared, self.ProjectionAxis)
      self.PerSegmentMarkedCellInfoList = perSegmentCellList
    #UpdatePipelineWithCurrentTimeArgument(self.myCopyOfInputFilter)
    #end new work

    if self.DoProgrammableFilterToAddMaskVariable:
      if PhactoriDbg(100):
        myDebugPrint3("cell and point data arrays before filter:\n")
        RecursivelyPrintPointAndCellArrayInformation(inInputFilter)
      self.mProgFilterMaskSetter = ProgrammableFilter(
              Input = self.myCopyOfInputFilter)
      self.mProgFilterMaskSetter.CopyArrays = 1
      self.CreateProgrammableFilterString()
      self.mProgFilterMaskSetter.Script = self.mProgFilterString
      self.mProgFilterMaskSetter.UpdatePipeline()
      if PhactoriDbg(100):
        myDebugPrint3("cell and point data arrays after filter:\n")
        RecursivelyPrintPointAndCellArrayInformation(self.mProgFilterMaskSetter)
    else:
      self.mProgFilterMaskSetter = self.myCopyOfInputFilter

    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.CreateParaViewFilter returning\n", 100)
    return self.mProgFilterMaskSetter

  def GetPidWithLeastValueList(self, inLocalDistSqrdList):
    myPid = int(SmartGetLocalProcessId())
    globalDistSqrdList = UseReduceOnFloatList(inLocalDistSqrdList, 1)
    localPidList = []
    numItems = len(inLocalDistSqrdList)
    for ndx in range(0,numItems):
      if globalDistSqrdList[ndx] == inLocalDistSqrdList[ndx]:
        localPidList.append(myPid)
      else:
        localPidList.append(-1)

    pidWithDataList = UseReduceOnIntegerList(localPidList, 0)
    return pidWithDataList, globalDistSqrdList

  def UseMpiToGetGlobalCellPointsClosestStructured(self, ioCellList, inLocalDistSqrdList):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.UseMpiToGetGlobalCellPointsClosestStructured entered\n", 100)

    pidWithDataList, globalDistSqrdList = self.GetPidWithLeastValueList(inLocalDistSqrdList)
    if PhactoriDbg(100):
      myDebugPrint3("pidWithDataList:\n" + str(pidWithDataList) + "\nglobalDistSqrdList:\n" + str(globalDistSqrdList) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("local nearest cell list (ioCellList) before mpi:\n")
      for ii, oneCellInfo in enumerate(ioCellList):
        myDebugPrint3(str(ii) + ": " + \
          str(inLocalDistSqrdList[ii]) + "\n" + \
          str(oneCellInfo.cellTestPoint) + "\n"  + \
          str(oneCellInfo.pid) + "\n"  + \
          str(oneCellInfo.leafVisitCount) + "\n"  + \
          str(oneCellInfo.ijk) + "\n")

    #convert cell point list to array of doubles and ints, use mpi reduce to share
    #the values, then convert back to cell point list
    serializeFloatArray = []
    serializeIntArray = []
    myPid = int(SmartGetLocalProcessId())
    for ii, oneCell in enumerate(ioCellList):
      if pidWithDataList[ii] == myPid:
        oneCell.SerializeAppendToFloatAndIntArray(serializeFloatArray, serializeIntArray)
      else:
        oneCell.SerializeAppendToFloatAndIntArrayZeroVersion(serializeFloatArray, serializeIntArray)

    #use mpi reduce to spread array correctly
    globalSerializeFloatArray = UseReduceOnFloatList(serializeFloatArray, 2)
    globalSerializeIntArray = UseReduceOnIntegerList(serializeIntArray, 2)

    #now create return global cell point list from arrays
    for ii, oneCell in enumerate(ioCellList):
      oneCell.SerializeSetFromFloatAndIntArray(globalSerializeFloatArray, globalSerializeIntArray, ii)

    if PhactoriDbg(100):
      myDebugPrint3("global nearest cell list (ioCellList) after mpi:\n")
      for ii, oneCellInfo in enumerate(ioCellList):
        myDebugPrint3(str(ii) + ": " + \
          str(globalDistSqrdList[ii]) + "\n" + \
          str(oneCellInfo.cellTestPoint) + "\n"  + \
          str(oneCellInfo.pid) + "\n"  + \
          str(oneCellInfo.leafVisitCount) + "\n"  + \
          str(oneCellInfo.ijk) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.UseMpiToGetGlobalCellPointsClosestStructured returning\n", 100)

    return ioCellList, globalDistSqrdList

  def UseMpiToGetGlobalCellPointsClosest(self, inLocalCellPointList, inLocalDistSqrdList):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.UseMpiToGetGlobalCellPointsClosest entered\n", 100)

    pidWithDataList, globalDistSqrdList = self.GetPidWithLeastValueList(inLocalDistSqrdList)
    if PhactoriDbg(100):
      myDebugPrint3("pidWithDataList:\n" + str(pidWithDataList) + "\nglobalDistSqrdList:\n" + str(globalDistSqrdList) + "\n")

    #convert cell point list to array of doubles and ints, use mpi reduce to share
    #the values, then convert back to cell point list
    serializeFloatArray = []
    serializeIntArray = []

    #convert cell point list to array of doubles
    cellPointFloatArray = []
    myPid = SmartGetLocalProcessId()
    for ii, oneCellPoint in enumerate(inLocalCellPointList):
      if pidWithDataList[ii] == myPid:
        cellPointFloatArray.append(oneCellPoint[0])
        cellPointFloatArray.append(oneCellPoint[1])
        cellPointFloatArray.append(oneCellPoint[2])
      else:
        cellPointFloatArray.append(0.0)
        cellPointFloatArray.append(0.0)
        cellPointFloatArray.append(0.0)

    #use mpi reduce to spread array correctly
    globalCellPointFloatArray = UseReduceOnFloatList(cellPointFloatArray, 2)

    #now create return global cell point list from arrays
    numCells = len(inLocalCellPointList)
    returnGlobalCellPointList = []
    for ii in range(0,numCells):
      myndx = ii*3
      oneCellPoint = [globalCellPointFloatArray[myndx],
                      globalCellPointFloatArray[myndx+1],
                      globalCellPointFloatArray[myndx+2]]
      returnGlobalCellPointList.append(oneCellPoint)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.UseMpiToGetGlobalCellPointsClosest returning\n", 100)
    return returnGlobalCellPointList, globalDistSqrdList

  @staticmethod
  def DebugPrintCellPointsFaces(oneCell, cellIndex, inInputCsData):
    if not PhactoriDbg(100):
      return
    myDebugPrint3("cellIndex: " + str(cellIndex) + "\n");
    cellPointIds = oneCell.GetPointIds()
    numPoints = cellPointIds.GetNumberOfIds()
    numFaces = oneCell.GetNumberOfFaces()
    myDebugPrint3("numFaces: " + str(numFaces) + " numPoints: " + str(numPoints) + " pointIds: " + str(cellPointIds) + "\n")
    for ii in range(0, numFaces):
      oneFace = oneCell.GetFace(ii)
      facePoints = oneFace.GetPointIds()
      myDebugPrint3("face " + str(ii) + ": " + str(facePoints) + "\n")
      onePointXyz = [0.0, 0.0, 0.0]
      for jj in range(0, facePoints.GetNumberOfIds()):
        onePtId = facePoints.GetId(jj)
        inInputCsData.GetPoint(onePtId, onePointXyz)
        myDebugPrint3("    point " + str(jj) + ": id: " + str(onePtId) + " xyz: " + str(onePointXyz) + "\n")

  @staticmethod
  def SetMaskValueForCellsNearSegmentsInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("SetMaskValueForCellsNearSegmentsInBlock entered\n")
    inParameters.leafVisitCount += 1
    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    if (numCells != 0) and (numPoints == 0):
      if PhactoriDbg(100):
        myDebugPrint3("case where numCells is " + str(numCells) + " and numPoints is 0")
        myDebugPrint3("number of cell arrays: " + str(inInputCsData.GetCellData().GetNumberOfArrays()) + "\n")
        myDebugPrint3("number of point arrays: " + str(inInputCsData.GetPointData().GetNumberOfArrays()) + "\n")

    if (numCells == 0) or (numPoints == 0):
      if PhactoriDbg(100):
        cellData = inInputCsData.GetCellData()
        if cellData != None:
          myDebugPrint3("cellData is None, numCells is 0\n")
        else:
          myDebugPrint3("cellData has " + str(cellData.GetNumberOfArrays()) + " arrays, numCells is 0\n")
      if PhactoriDbg(100):
        myDebugPrint3("SetMaskValueForCellsNearSegmentsInBlock returning (no cells or no points)\n")
      return

    cellData = inInputCsData.GetCellData()
    maskCellData = cellData.GetArray(inParameters.maskCellVariableName)
    #if maskCellData == None:
    #  maskCellData = cellData.GetArray(inParameters.maskCellVariableName)

    globalCellIdData = cellData.GetArray(inParameters.globalCellIdName)
    outputCellArray = cellData.GetArray(inParameters.CellDataArrayName)

    if PhactoriDbg(100):
      myDebugPrint3("numCells \n" + str(numCells) + "\n")
      myDebugPrint3("inParameters.segmentList:\n" + str(inParameters.segmentList) + "\n")
      myDebugPrint3("cellData: " + str(cellData) + "\n")
      myDebugPrint3("maskCellData: " + str(maskCellData) + "\n")
      myDebugPrint3("globalCellIdData: " + str(globalCellIdData) + "\n")
      #myDebugPrint3("inInputCsData \n" + str(inInputCsData) + "\n")
      myDebugPrint3("inInputCsData.GetClassName() \n" + str(inInputCsData.GetClassName()) + "\n")
    projectionAxis = inParameters.projectionAxis
    leafVisitCount = inParameters.leafVisitCount
    myPid = int(SmartGetLocalProcessId())
    cellBounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for cellIndex in range(0,numCells):
      oneCell = inInputCsData.GetCell(cellIndex)
      cellTestPoint = GetCellTestPoint(oneCell)
      if maskCellData != None:
        maskCellData.SetValue(cellIndex, 0.0)
      for segIndex, oneSegment in enumerate(inParameters.segmentList):
        cellTestPointToSegmentDistanceSquared = oneSegment.FindDistanceSquaredToPointProjected(cellTestPoint, projectionAxis)
        #if PhactoriDbg(100):
        #  if cellIndex%5000 == 0:
        #    myDebugPrint3(str(cellIndex) + " " + str(cellTestPointToSegmentDistanceSquared) + "  " + str(cellTestPoint) + "\n")
        floatSegIndex = float(segIndex) * 5.0
        if cellTestPointToSegmentDistanceSquared < inParameters.maskTestDistanceSquared:
          #nearestPointOnSegment = oneSegment.FindNearestPointOnSegmentToPointProjected(cellTestPoint, projectionAxis)
          #if IsPointInsideCellBoundingBoxProjected(oneCell, nearestPointOnSegment, projectionAxis):
          oneCell.GetBounds(cellBounds)
          if oneSegment.IntersectsBoundingBoxProjected(cellBounds, projectionAxis):
            PhactoriSegmentCellSampler3.DebugPrintCellPointsFaces(oneCell, cellIndex, inInputCsData)
            if maskCellData != None:
              #warning--if element is near two segments, both will be marked,
              #last will win
              maskCellData.SetValue(cellIndex, floatSegIndex)
            if globalCellIdData != None:
              thisCellId = globalCellIdData.GetValue(cellIndex)
            else:
              thisCellId = None
            if outputCellArray  == None:
              thisCellDataTuple = None
            else:
              thisCellDataTuple = outputCellArray.GetTuple(cellIndex)
            newCellInfo = PhactoriSampledCellInfo()
            newCellInfo.SetFromList([cellTestPoint, [-1,-1,-1], thisCellDataTuple, myPid, leafVisitCount, cellIndex, segIndex, -1])

            inParameters.markedCellSet[segIndex].append(newCellInfo)
            if PhactoriDbg(100):
              myDebugPrint3("new cell close to segment:\n")
              myDebugPrint3(newCellInfo.ToStr())


    if PhactoriDbg(100):
      for idx, cellsForThisSegment in enumerate(inParameters.markedCellSet):
        myDebugPrint3("after this block cells for segment (index, count): " + str(idx) + ", " + str(len(cellsForThisSegment)) + "\n")
        #myDebugPrint3(str(inParameters.markedCellSet[idx]) + "\n")
      myDebugPrint3("\n")

    if PhactoriDbg(100):
      myDebugPrint3("SetMaskValueForCellsNearSegmentsInBlock returning\n")

  def SetMaskValueForCellsNearSegments(self,
        inputFilter, segmentList, maskCellVariableName, globalCellIdName, maskTestDistanceSquared, projectionAxis):
    if PhactoriDbg(100):
      myDebugPrint3("SetMaskValueForCellsNearSegments entered\n")
      numCellArrays = inputFilter.CellData.GetNumberOfArrays()
      myDebugPrint3("numCellArrays: " + str(numCellArrays) + "\n")
      for ii in range(0,numCellArrays):
        myDebugPrint3(str(inputFilter.CellData.GetArray(ii)) + "\n")
    recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    recursionObj.mParameters = SetMaskValueForCellsNearSegmentsRecursionParams()
    recursionObj.mParameters.Initialize(segmentList, self.CellDataArrayName, maskCellVariableName, globalCellIdName, maskTestDistanceSquared, projectionAxis)
    recursionObj.mOperationToDoPerBlock = self.SetMaskValueForCellsNearSegmentsInBlock
    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, inputFilter)

    #if PhactoriDbg(100):
    #  myDebugPrint3("SetMaskValueForCellsNearSegments returning\n")
    #  numCellArrays = inputFilter.CellData.GetNumberOfArrays()
    #  myDebugPrint3("numCellArrays: " + str(numCellArrays) + "\n")
    #  for ii in range(0,numCellArrays):
    #    myDebugPrint3(str(inputFilter.CellData.GetArray(ii)) + "\n")
    #  UpdatePipelineWithCurrentTimeArgument(inputFilter)
    #  numCellArrays = inputFilter.CellData.GetNumberOfArrays()
    #  myDebugPrint3("numCellArrays: " + str(numCellArrays) + "\n")
    #  for ii in range(0,numCellArrays):
    #    myDebugPrint3(str(inputFilter.CellData.GetArray(ii)) + "\n")

    return recursionObj.mParameters.markedCellSet

  @staticmethod
  def GetCellsClosestToPointsInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("GetCellsClosestToPointsInBlock entered\n")
    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    if (numCells == 0) or (numPoints == 0):
      #no cells here
      if PhactoriDbg(100):
        myDebugPrint3("GetCellsClosestToPointsInBlock returning (no cells or no points)\n")
      return

    if PhactoriDbg(100):
      myDebugPrint3(str(inParameters.testPointList) + "\n")
      myDebugPrint3(str(inParameters.distSqrdList) + "\n")

    for cellIndex in range(0,numCells):
      oneCell = inInputCsData.GetCell(cellIndex)
      cellTestPoint = GetCellTestPoint(oneCell)
      for ptndx, oneTestPt in enumerate(inParameters.testPointList):
        testDist = vecDistanceSquared(oneTestPt, cellTestPoint)
        if testDist < inParameters.distSqrdList[ptndx]:
          inParameters.closestList[ptndx] = cellTestPoint
          inParameters.distSqrdList[ptndx] = testDist

    if PhactoriDbg(100):
      myDebugPrint3(str(inParameters.testPointList) + "\n")
      myDebugPrint3(str(inParameters.distSqrdList) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("after this block:\n")
      for ii, oneCellPoint in enumerate(inParameters.closestList):
        myDebugPrint3(str(ii) + ": " + \
          str(inParameters.distSqrdList[ii]) + "\n" + \
          str(inParameters.testPointList[ii]) + "\n" + str(oneCellPoint))
      myDebugPrint3("\n")

    if PhactoriDbg(100):
      myDebugPrint3("GetCellsClosestToPointsInBlock returning\n")

  def GetCellsClosestToPointsOnThisProcessFromParaViewFilter(self, inTestPointList):
    if PhactoriDbg(100):
      myDebugPrint3("GetCellsClosestToPointsOnThisProcessFromParaViewFilter entered\n")
    recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    recursionObj.mParameters = GetCellsClosestToPointsInBlockRecursionParams()
    recursionObj.mParameters.InitializeWithPointList(inTestPointList)
    recursionObj.mOperationToDoPerBlock = self.GetCellsClosestToPointsInBlock
    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.myCopyOfInputFilter)
    if PhactoriDbg(100):
      myDebugPrint3("GetCellsClosestToPointsOnThisProcessFromParaViewFilter returning\n")
    return recursionObj.mParameters.closestList, recursionObj.mParameters.distSqrdList

  @staticmethod
  def GetCellsListClosestToPointsListInBlockStructured(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("GetCellsListClosestToPointsListInBlockStructured entered\n")
    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    inParameters.leafVisitCount += 1
    if (numCells == 0) or (numPoints == 0):
      #no cells here
      if PhactoriDbg(100):
        myDebugPrint3("GetCellsListClosestToPointsListInBlockStructured returning (no cells or no points)\n")
      return

    mypid = SmartGetLocalProcessId()
    #myExtent = [0,0,0,0,0,0]
    myExtent = inInputCsData.GetExtent()
    if PhactoriDbg(100):
      myDebugPrint3("myExtent: " + str(myExtent) + "\n")
    for ii in range(myExtent[0], myExtent[1]):
      for jj in range(myExtent[2], myExtent[3]):
        for kk in range(myExtent[4], myExtent[5]):
          #oneCell = inInputCsData.GetCell(cellIndex)
          #oneCell = inInputCsData.GetCell(ii, jj, kk)
          oneCell = localGetCellijk(ii, jj, kk, inInputCsData, myExtent)
          cellTestPoint = GetCellTestPoint(oneCell)
          if PhactoriDbg(100):
            if ((ii==myExtent[0]) and (jj == myExtent[2]) and (kk == myExtent[4])):
              myDebugPrint3("cellTestPoint " + str([ii,jj,kk]) + " " + str(cellTestPoint) + "\n")
              oneCellBounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
              oneCell.GetBounds(oneCellBounds)
              myDebugPrint3("oneCellBounds " + str([ii,jj,kk]) + " " + str(oneCellBounds) + "\n")
          for ptndx, oneTestPt in enumerate(inParameters.testPointList):
            testDist = vecDistanceSquared(oneTestPt, cellTestPoint)
            if testDist < inParameters.distSqrdList[ptndx]:
              inParameters.distSqrdList[ptndx] = testDist
              cellInfo = inParameters.closestCellList[ptndx]
              cellInfo.pid = mypid
              cellInfo.SetCellTestPoint(cellTestPoint)
              cellInfo.leafVisitCount = inParameters.leafVisitCount
              cellInfo.SetIjk(ii, jj, kk)

    if PhactoriDbg(100):
      myDebugPrint3(str(inParameters.testPointList) + "\n")
      myDebugPrint3(str(inParameters.distSqrdList) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("after this block:\n")
      for ii, oneCellInfo in enumerate(inParameters.closestCellList):
        myDebugPrint3(str(ii) + ": " + \
          str(inParameters.distSqrdList[ii]) + "\n" + \
          str(inParameters.testPointList[ii]) + "\n" + \
          str(oneCellInfo.cellTestPoint) + "\n"  + \
          str(oneCellInfo.pid) + "\n"  + \
          str(oneCellInfo.leafVisitCount) + "\n"  + \
          str(oneCellInfo.ijk) + "\n")
      myDebugPrint3("\n")

    if PhactoriDbg(100):
      myDebugPrint3("GetCellsListClosestToPointsListInBlockStructured returning\n")

  def GetListOfCellTestPointsNearestListOfPointsStructured(self, inTestPointList):
    if PhactoriDbg(100):
      myDebugPrint3("GetListOfCellTestPointsNearestListOfPointsStructured entered\n")
    #find the closest structured cells for this process, recursing into
    #multiblock dataset
    recursionParams = CellsListClosestToPointsListInBlockStructuredRecursionParams()
    recursionParams.InitializeWithPointList(inTestPointList)
    recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    recursionObj.mParameters = recursionParams
    recursionObj.mOperationToDoPerBlock = self.GetCellsListClosestToPointsListInBlockStructured
    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.myCopyOfInputFilter)

    #use mpi to find closest structured cells across all processes
    nearestCellList, distanceList = self.UseMpiToGetGlobalCellPointsClosestStructured(
      recursionObj.mParameters.closestCellList, recursionObj.mParameters.distSqrdList)

    if PhactoriDbg(100):
      myDebugPrint3("GetListOfCellTestPointsNearestListOfPointsStructured returning\n")
    return nearestCellList

  @staticmethod
  def GatherStructuredCellsInBlockFromSeedCells(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("GatherStructuredCellsInBlockFromSeedCells entered\n")

    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    inParameters.leafVisitCount += 1
    if (numCells == 0) or (numPoints == 0):
      #no cells here
      return

    cellData = inInputCsData.GetCellData()
    outputCellArray = None
    if cellData != None:
      outputCellArray = cellData.GetArray(inParameters.dataArrayName)

    if outputCellArray != None:
      dataArrayNumCmpnts = outputCellArray.GetNumberOfComponents()
      defaultTuple = []
      for ii in range(0, dataArrayNumCmpnts):
        defaultTuple.append(-1.0)
    else:
      dataArrayNumCmpnts = -1
      defaultTuple = []

    myPid = SmartGetLocalProcessId()
    myExtent = inInputCsData.GetExtent()
    myDims = inInputCsData.GetDimensions()
    if PhactoriDbg(100):
      myDebugPrint3("myExtent: " + str(myExtent) + \
                    "\nmyDims: " + str(myDims) + \
                    "\nnumCells: " + str(numCells) + "\n")
    seedCellList = inParameters.SeedCells
    for seedCellIndex, oneSeedCell in enumerate(seedCellList):
      if PhactoriDbg(100):
        myDebugPrint3("oneSeedCell: " + str(seedCellIndex) + " " + str(oneSeedCell.ijk) + " " + str([oneSeedCell.leafVisitCount, inParameters.leafVisitCount]) + " " + str(myExtent) + "\n")
      #is this seed cell in this block?
      if oneSeedCell.leafVisitCount != inParameters.leafVisitCount:
        if PhactoriDbg(100):
          myDebugPrint3("this block is not for this seed cell: " + \
            str(oneSeedCell.leafVisitCount) + " " + \
            str(inParameters.leafVisitCount) + "\n")
        continue
      if PhactoriDbg(100):
        myDebugPrint3("this block is for this seed cell: " + \
           str(oneSeedCell.leafVisitCount) + " " + \
           str(inParameters.leafVisitCount) + "\n")
      #does this process have any data for this seed cell?
      if oneSeedCell.AxisCrossesExtent(myExtent, oneSeedCell.collectionAxis) == False:
        if PhactoriDbg(100):
          myDebugPrint3("this seed cell not in myExtent: " + str(oneSeedCell.ijk) + "\n")
        continue
      ijkrange = oneSeedCell.GetIntersectingCollectionExtent(myExtent, oneSeedCell.collectionAxis)
      if PhactoriDbg(100):
        myDebugPrint3("this block has cells to collect: " + str(ijkrange) + "\n")
      #collect all the cells for this seed cell
      for ii in range(ijkrange[0], ijkrange[1]):
        for jj in range(ijkrange[2], ijkrange[3]):
          for kk in range(ijkrange[4], ijkrange[5]):
            newStructuredCell = PhactoriSampledCellInfo()
            newStructuredCell.Populate1(ii, jj, kk, myPid,
              inParameters.leafVisitCount, seedCellIndex, inInputCsData,
              myExtent, numCells,
              outputCellArray, dataArrayNumCmpnts, defaultTuple)
            inParameters.CellsPerSeedCell[seedCellIndex].append(newStructuredCell)

    if PhactoriDbg(100):
      myDebugPrint3("GatherStructuredCellsInBlockFromSeedCells returning\n")

  def GatherStructuredCellsFromSeedCells(self):
    """for each structured seed cell in self.StructuredNearbyCellPointList,
       gather cells along one u/v/w direction in the structured data (on this
       process), including the current data for the cells"""
    if PhactoriDbg(100):
      myDebugPrint3("GatherStructuredCellsFromSeedCells entered\n")

    #find the closest structured cells for this process, recursing into
    #multiblock dataset
    recursionParams = GatherStructuredCellsFromSeedCellsRecursionParams()
    recursionParams.SetUpForRecursion(self)
    recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    recursionObj.mParameters = recursionParams
    recursionObj.mOperationToDoPerBlock = self.GatherStructuredCellsInBlockFromSeedCells
    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.myCopyOfInputFilter)

    if PhactoriDbg(100):
      myDebugPrint3("GatherStructuredCellsFromSeedCells returning\n")
    return recursionParams.CellsPerSeedCell


  def GetListOfCellTestPointsNearestListOfPoints(self, pointList):
    """for each point in the list, find the cell test point (e.g. center of
       cell bounding box) which is nearest the test point.  Use MPI to work
       in parallel"""

    thisProcessNearestCellPointList, thisProcDistSqrdList = \
      self.GetCellsClosestToPointsOnThisProcessFromParaViewFilter(pointList)

    nearestCellList, distanceList = self.UseMpiToGetGlobalCellPointsClosest(
      thisProcessNearestCellPointList, thisProcDistSqrdList)

    return nearestCellList

  def WritePerSegmentCellListToFile(self, filename, perSegmentCellList):
    ff = open(filename, "w")
    if self.collectionMethod == "seed cell with structured grid":
      ff.write("{\n")
      ff.write('"format for cells in lists": ' + PhactoriSampledCellInfo.TerseOneLineJsonFormatComment() + ",\n")
      ff.write('"list of collected cells for each seed cell": [\n')
      lastIdx = len(perSegmentCellList) - 1
      for idx, cellsForThisSegment in enumerate(perSegmentCellList):
        ff.write("{\n")
        ff.write('"cells for seed cell index": ' + str(idx) + ',\n')
        thisStructuredSeedCell = self.StructuredNearbyCellPointList[idx]
        nearbyGeometryPoint = self.NearbyGeometryPointList[idx]
        ff.write('"structured seed cell xyz": ' + str(thisStructuredSeedCell.cellTestPoint) + ',\n')
        ff.write('"structured seed cell ijk": ' + str(thisStructuredSeedCell.ijk) + ',\n')
        ff.write('"structured seed cell complete": ' + thisStructuredSeedCell.ToStrTerseOneLineList() + ',\n')
        ff.write('"geometry point used to find seed cell": ' + str(nearbyGeometryPoint) + ',\n')
        ff.write('"number of cells collected for this seed cell": ' + str(len(cellsForThisSegment)) + ',\n')
        ff.write('"cell list for this seed cell": [\n')
        for cellCntr, oneCellInfo in enumerate(cellsForThisSegment):
          if cellCntr != 0:
            ff.write(",\n")
          cellStr = oneCellInfo.ToStrTerseOneLineList()
          ff.write(cellStr)
        ff.write("\n]\n")
        if idx != lastIdx:
          ff.write("},\n")
        else:
          ff.write("}\n")
      ff.write("]\n")
      ff.write("}\n")
    else:
      ff.write("{\n")
      ff.write('"format for cells in lists": ' + PhactoriSampledCellInfo.TerseOneLineJsonFormatComment() + ",\n")
      ff.write('"list of collected cells for each segment": [\n')
      lastIdx = len(perSegmentCellList) - 1
      for idx, cellsForThisSegment in enumerate(perSegmentCellList):
        ff.write("{\n")
        ff.write('"cells for segment index": ' + str(idx) + ',\n')
        thisSegment = self.SegmentList[idx]
        nearbyGeometryPoint = self.NearbyGeometryPointList[idx]
        ff.write('"segment point A": ' + str(thisSegment.ptA) + ',\n')
        ff.write('"segment point B": ' + str(thisSegment.ptB) + ',\n')
        ff.write('"geometry point used to find point A": ' + str(nearbyGeometryPoint) + ',\n')
        ff.write('"number of cells collected for this segment": ' + str(len(cellsForThisSegment)) + ',\n')
        ff.write('"cell list for this segment": [\n')
        for cellCntr, oneCellInfo in enumerate(cellsForThisSegment):
          if cellCntr != 0:
            ff.write(",\n")
          cellStr = oneCellInfo.ToStrTerseOneLineList()
          ff.write(cellStr)
        ff.write("\n]\n")
        if idx != lastIdx:
          ff.write("},\n")
        else:
          ff.write("}\n")
      ff.write("]\n")
      ff.write("}\n")
    ff.close()

  def WriteAllDataFromOneProcessUsingMPI(self, perSegmentCellList, basename = "PhactoriSegmentCellSampler3_", outputFilenameCounter = -1):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.WriteAllDataFromOneProcessUsingMPI entered\n", 100)
    myPid = SmartGetLocalProcessId()

    ##to start with, just have each process output it's data
    if PhactoriDbg(100):
      thisProcessHasData = False
      for idx, cellsForThisSegment in enumerate(perSegmentCellList):
        if len(cellsForThisSegment) != 0:
          thisProcessHasData = True
          break
      if self.collectionMethod == "seed cell with structured grid":
        structuredFlag = True
      else:
        structuredFlag = False
      if structuredFlag:
        basename = "PhactoriSegmentCellSampler3_structured_"
      else:
        basename = "PhactoriSegmentCellSampler3_"
      fname = basename + str(myPid) + self.SampledCellsOutputFilenameExtension
      if thisProcessHasData:
        self.WritePerSegmentCellListToFile(fname, perSegmentCellList)

    if PhactoriDbg(100):
      myDebugPrint3("about to serialize cell infos\n", 100)

    tupleSize = self.CellDataArrayTupleSize
    cellListSerializedFloatArray = []
    cellListSerializedIntArray = []
    for idx, cellsForThisSegment in enumerate(perSegmentCellList):
      for cellCntr, oneCellInfo in enumerate(cellsForThisSegment):
        oneCellInfo.SerializeAppendToFloatAndIntArray(cellListSerializedFloatArray, cellListSerializedIntArray, tupleSize)

    if PhactoriDbg(100):
      myDebugPrint3("done serializing " + str(len(cellListSerializedFloatArray)) + "  " + str(len(cellListSerializedIntArray)) + "\n", 100)

    pidToDoOutput = SmartGetNumberOfProcesses() - 1
    globalSerializedFloatArray = UseMpiToSendAllProcessesFloatArrayToOneProcess(cellListSerializedFloatArray, pidToDoOutput)
    globalSerializedIntArray = UseMpiToSendAllProcessesIntArrayToOneProcess(cellListSerializedIntArray, pidToDoOutput)

    if PhactoriDbg(100):
      myDebugPrint3("returned from getting globalSerializedFloatArray and globalSerializedIntArray\n")

    if myPid == pidToDoOutput:
      if PhactoriDbg(100):
       myDebugPrint3("I'm the output pid, about to serialize in\n")

      serialFloatSize, serialIntSize = PhactoriSampledCellInfo.GetSerializeFloatAndIntSize(tupleSize)
      globalNumCells = int(len(globalSerializedFloatArray) / serialFloatSize)

      globalPerSegmentCellList = []
      if self.collectionMethod == "seed cell with structured grid":
        structuredFlag = True
      else:
        structuredFlag = False
      for ii in perSegmentCellList:
        globalPerSegmentCellList.append([])
      for ii in range(0, globalNumCells):
        newCell = PhactoriSampledCellInfo()
        newCell.SerializeSetFromFloatAndIntArray(globalSerializedFloatArray, globalSerializedIntArray, ii, tupleSize)
        whichSegList = globalPerSegmentCellList[newCell.segmentIndex]
        whichSegList.append(newCell)

      if PhactoriDbg(100):
        myDebugPrint3("finished serialize, about to write to file\n")
        if PhactoriDbg(100):
         for ii, oneSegmentGlobalCellList in enumerate(globalPerSegmentCellList):
           myDebugPrint3(str(ii) + "  " + str(len(oneSegmentGlobalCellList)) + "\n")

      fname = self.CreateCurrentCallbackCellCollectionOutputFilePath(outputFilenameCounter)
      if PhactoriDbg(100):
        myDebugPrint3("output filename will be: " + str(fname) + "\n")
      self.WritePerSegmentCellListToFile(fname, globalPerSegmentCellList)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.WriteAllDataFromOneProcessUsingMPI returning\n", 100)

  def CreateCurrentCallbackCellCollectionOutputFilePath(self, outputFilenameCounter):
    if outputFilenameCounter < 0:
      ftc = GetFrameTagCounter()
    else:
      ftc = outputFilenameCounter

    outName = self.SampledCellsOutputDirectory + "/" + \
      self.SampledCellsOutputFilename + \
      str(ftc).zfill(self.NumCounterDigits) + \
      self.SampledCellsOutputFilenameExtension

    return outName

  def WriteAllDataFromOneProcessUsingMPIStructured(self, perSegmentCellList, outputFilenameCounter = -1):
    self.WriteAllDataFromOneProcessUsingMPI(perSegmentCellList, "PhactoriSegmentCellSampler3_structured_", outputFilenameCounter)

  def ExportOperationData(self, datadescription):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.ExportOperationData entered\n", 100)

    UpdatePipelineWithCurrentTimeArgument(self.myCopyOfInputFilter)

    if self.collectionMethod == "seed cell with structured grid":
      perSegmentStructuredCellList = self.GatherStructuredCellsFromSeedCells()
      if self.WriteCellListToJsonFile:
        self.WriteAllDataFromOneProcessUsingMPIStructured(perSegmentStructuredCellList)
    else:
      perSegmentCellList = self.SetMaskValueForCellsNearSegments(
        self.myCopyOfInputFilter, self.SegmentList, "blargmask1", "Ids",
        self.MaskTestDistanceSquared, self.ProjectionAxis)

      if PhactoriDbg(100):
        for idx, cellsForThisSegment in enumerate(perSegmentCellList):
          myDebugPrint3("after recursion cells for segment (index, count): " + str(idx) + ", " + str(len(cellsForThisSegment)) + "\n")
          #myDebugPrint3(str(cellsForThisSegment) + "\n")
        myDebugPrint3("\n")

      if self.WriteCellListToJsonFile:
        self.WriteAllDataFromOneProcessUsingMPI(perSegmentCellList)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSegmentCellSampler3.ExportOperationData returning\n", 100)

  def CreateProgrammableFilterString(self):

    scriptLines = []

    scriptLines.append("#print('test2')\n")

    scriptLines.append("#cellIndexesToSet = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]\n")
    scriptLines.append("#maskValuesToSet = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]\n")
    scriptLines.append("#leafVisitCount = [0,0,0,0,0,0,0,0,0,0,0]\n")

    scriptLines.append("localLeafVisitCount = 0\n")

    scriptLines.append("def flatten(input, output):\n")
    scriptLines.append("    # Copy the cells etc.\n")
    scriptLines.append("    output.ShallowCopy(input)\n")
    scriptLines.append("    celldata = output.GetCellData()\n")
    scriptLines.append("    #for now if there are no cell arrays here we don't want to add one\n")
    scriptLines.append("    #really we don't want to add one of there are actually no cells\n")
    scriptLines.append("    numCellArrays = celldata.GetNumberOfArrays()\n")
    scriptLines.append("    if numCellArrays == 0:\n")
    scriptLines.append("      return\n")
    scriptLines.append("    numCells = input.GetNumberOfCells()\n")
    scriptLines.append("    numPoints = input.GetNumberOfPoints()\n")
    scriptLines.append("    #print(str(numCells))\n")
    scriptLines.append("    ncda = vtk.vtkDoubleArray()\n")
    scriptLines.append("    ncda.SetNumberOfTuples(numCells)\n")
    scriptLines.append("    numTuples = ncda.GetNumberOfTuples()\n")
    scriptLines.append("    for ii in range(0, numTuples):\n")
    scriptLines.append("        ncda.SetValue(ii, 0.0)\n")
    scriptLines.append("    for ii, cellNdx in enumerate(cellIndexesToSet):\n")
    scriptLines.append("        if localLeafVisitCount == leafVisitCount[ii]:\n")
    scriptLines.append("          if cellNdx < numCells:\n")
    scriptLines.append("            ncda.SetValue(cellNdx, float(maskValuesToSet[ii]))\n")
    scriptLines.append("    ncda.SetName('" + self.ProgrammableFilterOutputCellVariableName + "')\n")
    scriptLines.append("    #print(ncda.GetName())\n")
    scriptLines.append("    celldata.AddArray(ncda)\n")

    scriptLines.append("input = self.GetInputDataObject(0, 0)\n")
    scriptLines.append("output = self.GetOutputDataObject(0)\n")

    scriptLines.append("if input.IsA('vtkMultiBlockDataSet'):\n")
    scriptLines.append("    output.CopyStructure(input)\n")
    scriptLines.append("    iter = input.NewIterator()\n")
    scriptLines.append("    iter.UnRegister(None)\n")
    scriptLines.append("    iter.InitTraversal()\n")
    scriptLines.append("    while not iter.IsDoneWithTraversal():\n")
    scriptLines.append("        localLeafVisitCount += 1\n")
    scriptLines.append("        curInput = iter.GetCurrentDataObject()\n")
    scriptLines.append("        curOutput = curInput.NewInstance()\n")
    scriptLines.append("        curOutput.UnRegister(None)\n")
    scriptLines.append("        output.SetDataSet(iter, curOutput)\n")
    scriptLines.append("        flatten(curInput, curOutput)\n")
    scriptLines.append("        iter.GoToNextItem();\n")
    scriptLines.append("else:\n")
    scriptLines.append("  flatten(input, output)\n")

    mainScriptString = "".join(scriptLines)

    myCellIndexesToSet = []
    myMaskValuesToSet = []
    myLeafVisitCount = []
    if self.PerSegmentMarkedCellInfoList != None:
      for idx, cellsForThisSegment in enumerate(self.PerSegmentMarkedCellInfoList):
        if cellsForThisSegment != None:
         maskValue = idx+1
         for cellIdx, oneCellInfo in enumerate(cellsForThisSegment):
           myCellIndexesToSet.append(oneCellInfo.index)
           myMaskValuesToSet.append(maskValue)
           myLeafVisitCount.append(oneCellInfo.leafVisitCount)

    newstr1 = "cellIndexesToSet = " + str(myCellIndexesToSet) + "\n"
    newstr2 = "maskValuesToSet = " + str(myMaskValuesToSet) + "\n"
    newstr3 = "leafVisitCount = " + str(myLeafVisitCount) + "\n"
    self.mProgFilterString = newstr1 + newstr2 + newstr3 + mainScriptString

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMarkCellSetsOperation constructed script:\n" + \
        self.mProgFilterString)

#phactori_combine_to_single_python_file_subpiece_end_1
