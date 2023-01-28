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
import sys
from .PhactoriParaviewMultiBlockRecursion import *
from .PhactoriVectorLibrary import *

#phactori_combine_to_single_python_file_subpiece_begin_1

global gCellStatusTrack
gCellStatusTrack = []
for ii in range(0,25):
  gCellStatusTrack.append(0)

class CellEdgeAngleMetricParameters():
  def __init__(self):
    self.numPointsPerBlock = None
    self.numCellsPerBlock = None
    self.UseSmallestAngle = True
    self.OffsetIndex = 0

  def Initialize(self, inUseSmallestAngle, inOffsetIndex):
    self.numPointsPerBlock = []
    self.numCellsPerBlock = []
    self.UseSmallestAngle = inUseSmallestAngle
    self.OffsetIndex = inOffsetIndex

class PhactoriFollowSurfaceSide1(PhactoriOperationSpecifics):
  """experimental filter to find the lengths of the edges of cells in the mesh"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.inPvFilter = None
    self.pfcelProgrammableFilter = None
    self.mProgFilterString = None
    self.ProgrammableFilterOutputCellVariableName = "surfacestatus"
    self.UseSmallestAngle = True
    self.OffsetIndex = 0

  def ParseParametersFromJson(self, inJson):

    keyval3 = "output nth smallest angle"
    keyval4 = "output nth largest angle"
    if keyval3 in inJson:
      self.UseSmallestAngle = True
      self.OffsetIndex = inJson[keyval3]
    elif keyval4 in inJson:
      self.UseSmallestAngle = False
      self.OffsetIndex = inJson[keyval3]
    else:
      self.UseSmallestAngle = True
      self.OffsetIndex = 0

    keyval10 = "output cell variable name"
    if keyval10 in inJson: 
      self.ProgrammableFilterOutputCellVariableName = inJson[keyval10]

  def CreateParaViewFilter(self, inInputFilter):
    """create the MergeBlocks filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFollowSurfaceSide1.CreateParaViewFilter "
          "entered\n", 100)
    self.inPvFilter = inInputFilter

    #savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    #newParaViewFilter = MergeBlocks(inInputFilter)

    #SetActiveSource(newParaViewFilter)
    #UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)
    #SetActiveSource(savedActiveSource)

    self.pfcelProgrammableFilter = ProgrammableFilter(Input = self.inPvFilter)
    self.pfcelProgrammableFilter.CopyArrays = 1
    self.CreateProgrammableFilterString()
    self.pfcelProgrammableFilter.Script = self.mProgFilterString
    self.pfcelProgrammableFilter.UpdatePipeline()

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFollowSurfaceSide1.CreateParaViewFilter "
          "returning\n", 100)

    return self.pfcelProgrammableFilter

  #@staticmethod
  #def IncrementItemInCellsTouchCountMap(onePt, cellsTouchCountMap):
  #  xx = onePt[0]
  #  yy = onePt[1]
  #  zz = onePt[2]
  #  if xx not in cellsTouchCountMap:
  #    cellsTouchCountMap[xx] = {}
  #  lvlx = cellsTouchCountMap[xx]
  #  if yy not in lvlx:
  #    lvlx[yy] = {}
  #  lvly = lvlx[yy]
  #  if zz not in lvly:
  #    lvly[zz] = 0
  #  lvlz = lvly[zz]
  #  lvly[zz] = lvlz + 1
    
  @staticmethod
  def IncrementTouchingCellCountForAllPointsInOneCell(inInputCsData, inParameters, inCellIndex,
        cellsTouchingPointCount):
    oneCell = inInputCsData.GetCell(inCellIndex)

    numPoints = oneCell.GetNumberOfPoints()
    pointIds = oneCell.GetPointIds()
    onePt = [0.0,0.0,0.0]
    for ii in range(0, numPoints):
      onePointId = pointIds.GetId(ii)
      cellsTouchingPointCount[onePointId] += 1
      inInputCsData.GetPoint(onePointId, onePt)

  @staticmethod
  def GetEdgeLength(inInputCsData, oneFace, edgeIndex):
    oneEdge = oneFace.GetEdge(edgeIndex)
    edgePointIds = oneEdge.GetPointIds()
    ptA = inInputCsData.GetPoint(edgePointIds.GetId(0))
    ptB = inInputCsData.GetPoint(edgePointIds.GetId(1))
    xx = ptB[0] - ptA[0]
    yy = ptB[1] - ptA[1]
    zz = ptB[2] - ptA[2]
    return math.sqrt(xx*xx + yy*yy + zz*zz)

  @staticmethod
  def GetFaceSize(inInputCsData, oneFace):
    numFaceEdges = oneFace.GetNumberOfEdges()
    if numFaceEdges != 4:
      myDebugPrint3AndException("numFaceEdges != 4 error\n")
      #use perimeter as a stand-in for area on non quads
      retVal = 0
      for ii in range(0, numFaceEdges):
        retVal += GetEdgeLength(inInputCsData, oneFace, ii)
    else:
      #assume roughly rectangular calculation for size
      side0Len = PhactoriFollowSurfaceSide1.GetEdgeLength(inInputCsData, oneFace, 0)
      side1Len = PhactoriFollowSurfaceSide1.GetEdgeLength(inInputCsData, oneFace, 1)
      side2Len = PhactoriFollowSurfaceSide1.GetEdgeLength(inInputCsData, oneFace, 2)
      side3Len = PhactoriFollowSurfaceSide1.GetEdgeLength(inInputCsData, oneFace, 3)
      retVal = (side0Len + side2Len) * (side1Len * side3Len)
    return retVal

  @staticmethod
  def GetSortedFaceSizeList(inInputCsData, oneCell):
    faceSizeList = []
    numFaces = oneCell.GetNumberOfFaces()
    for faceNdx in range(0, numFaces):
      oneFace = oneCell.GetFace(faceNdx)
      faceSizeList.append(PhactoriFollowSurfaceSide1.GetFaceSize(inInputCsData, oneFace))
    return sorted(faceSizeList)

  @staticmethod
  def IsFaceOneOfTwoLargestOnCell(inInputCsData, oneCell, oneFace):
    SortedFaceSizeList = PhactoriFollowSurfaceSide1.GetSortedFaceSizeList(inInputCsData, oneCell)
    testFaceSize = PhactoriFollowSurfaceSide1.GetFaceSize(inInputCsData, oneFace)
    if PhactoriDbg(100):
      myDebugPrint3("testFaceSize: " + str(testFaceSize) + " SortedFaceSizeList:\n" + str(SortedFaceSizeList) + "\n")
      #for ii, ss in enumerate(SortedFaceSizeList):
      #  myDebugPrint3(str(ii) + ": " + str(ss) +":\n")
    if len(SortedFaceSizeList) < 2:
      #shouldn't happen
      myDebugPrint3AndException("len(SortedFaceSizeList) < 2 error\n")
      return False
    if PhactoriDbg(100):
      myDebugPrint3("testFaceSize >= SortedFaceSizeList[-2]: " + str(testFaceSize >= SortedFaceSizeList[-2]) + "\n")
    return (testFaceSize >= SortedFaceSizeList[-2])

  @staticmethod
  def SurfaceCornerCellFaceFacesExterior(oneFace, cellsTouchingPointCount):
    numFacePoints = oneFace.GetNumberOfPoints()
    facePointIds = oneFace.GetPointIds()
    numOutsideVertices = 0
    for ii in range(0, numFacePoints):
      if cellsTouchingPointCount[ii] <= 0:
        return True
    return False

  @staticmethod
  def TestOutSideSurfaceCornerCellForMaxAreaFacing(
    inInputCsData, oneCell, pointIds, cellsTouchingPointCount):
    #faces with 4 2 2 0 need to be tested
    numFaces = oneCell.GetNumberOfFaces()
    for faceIndex in range(0, numFaces):
      oneFace = oneCell.GetFace(faceIndex)
      if PhactoriFollowSurfaceSide1.SurfaceCornerCellFaceFacesExterior(oneFace, cellsTouchingPointCount):
        if PhactoriFollowSurfaceSide1.IsFaceOneOfTwoLargestOnCell(inInputCsData, oneCell, oneFace):
          #this cell is facing the big area cell side
          return 7
    return 4

  @staticmethod
  def SurfaceEdgeCellFaceFacesExterior(oneFace, cellsTouchingPointCount):
    numFacePoints = oneFace.GetNumberOfPoints()
    facePointIds = oneFace.GetPointIds()
    numOutsideVertices = 0
    for ii in range(0, numFacePoints):
      if cellsTouchingPointCount[ii] <= 2:
        numOutsideVertices += 1
    return (numOutsideVertices >= 2)

  @staticmethod
  def TestOutSideSurfaceEdgeCellForMaxAreaFacing(
    inInputCsData, oneCell, pointIds, cellsTouchingPointCount):
    #faces with 4 4 2 2 need to be tested
    numFaces = oneCell.GetNumberOfFaces()
    for faceIndex in range(0, numFaces):
      oneFace = oneCell.GetFace(faceIndex)
      if PhactoriFollowSurfaceSide1.SurfaceEdgeCellFaceFacesExterior(oneFace, cellsTouchingPointCount):
        if PhactoriFollowSurfaceSide1.IsFaceOneOfTwoLargestOnCell(inInputCsData, oneCell, oneFace):
          #this cell is facing the big area cell side
          return 6
    return 3

  @staticmethod
  def TestOutSideSurfaceCellForMaxAreaFacing(
        inInputCsData, oneCell, pointIds, cellsTouchingPointCount):
    #find face which face faces exterior
    numFaces = oneCell.GetNumberOfFaces()
    if PhactoriDbg(100):
      for faceIndex in range(0, numFaces):
        oneFace = oneCell.GetFace(faceIndex)
        numFacePoints = oneFace.GetNumberOfPoints()
        facePtNdxs = oneFace.GetPointIds()
        nbrCntList = []
        for jj in range(0, numFacePoints):
          ptNdx = facePtNdxs.GetId(jj)
          nbrCntList.append(cellsTouchingPointCount[ptNdx])
        if PhactoriDbg(100):
          myDebugPrint3("faceIndex " + str(faceIndex) + " nbrCntList " + str(nbrCntList) + "\n")
    for faceIndex in range(0, numFaces):
      if PhactoriDbg(100):
        myDebugPrint3("testing faceIndex " + str(faceIndex) + "\n")
      oneFace = oneCell.GetFace(faceIndex)
      numFacePoints = oneFace.GetNumberOfPoints()
      facePtNdxs = oneFace.GetPointIds()
      for jj in range(0, numFacePoints):
        ptNdx = facePtNdxs.GetId(jj)
        if PhactoriDbg(100):
          myDebugPrint3("jj " + str(jj) + "  ctpc " + str(cellsTouchingPointCount[ptNdx]) + "\n")
        if cellsTouchingPointCount[ptNdx] != 4:
          #this is not the outward facing face
          if PhactoriDbg(100):
            myDebugPrint3("not 4, breaking loop\n")
          break
      if jj >= numFacePoints:
        #this is the outward facing face
        if PhactoriDbg(100):
          myDebugPrint3("all 4s, this is it, breaking loop\n")
        break
      if PhactoriDbg(100):
        myDebugPrint3("not all 4s, this not it, continuing loop\n")
    if PhactoriDbg(100):
      myDebugPrint3("faceIndex 4 4 4 4: " + str(faceIndex) + "\n")
    if faceIndex >= numFaces:
      #shouldn't happen, didn't find face with all 4's
      myDebugPrint3AndException("faceIndex >= numFaces error\n")
      return -2
    #now see if this outward facing face is one of the two largest
    if PhactoriFollowSurfaceSide1.IsFaceOneOfTwoLargestOnCell(inInputCsData, oneCell, oneFace):
      #this cell is on the surface and the outward facing side is large
      return 14 + faceIndex
    else:
      return 8 + faceIndex

  @staticmethod
  def FindSurfaceStatusForOneCell(inInputCsData, inParameters, inCellIndex,
        cellsTouchingPointCount):
    oneCell = inInputCsData.GetCell(inCellIndex)

    #cells with all 8 are interior
    #cells with 4 and 8 are outside surface
    #cells with 2 (and no 0) are outside edge
    #cells with 0 are corner

    #if cell is outside surface (4 and 8), if the face which is all 4 is the
    #largest or second largest area, then it is on the surface we want

    #if the cell is edge (2 4 8), there should be two sides which are 2 2 4 4
    #if either of these sides is the second or largest area then it is on the
    #surface we want.

    #if the cell is corner, there should be three sides which are 0 2 2 4
    #if any of thee sides is the second or largest area then it is on the
    #surface we want

    numPoints = oneCell.GetNumberOfPoints()
    pointIds = oneCell.GetPointIds()
    minCtpc = 9
    for ii in range(0, numPoints):
      testCtpc = cellsTouchingPointCount[pointIds.GetId(ii)]
      if testCtpc < minCtpc:
        minCtpc = testCtpc

    if minCtpc >= 8:
      #interior cell
      retStatus = 1
    elif minCtpc >= 4:
      #outside surface, not edge or corner
      retStatus = 2
      #retStatus = PhactoriFollowSurfaceSide1.TestOutSideSurfaceCellForMaxAreaFacing(
      #  inInputCsData, oneCell, pointIds, cellsTouchingPointCount)
    elif minCtpc >= 2:
      #outside surface edge, not corner
      retStatus = 3
      #retStatus = PhactoriFollowSurfaceSide1.TestOutSideSurfaceEdgeCellForMaxAreaFacing(
      #  inInputCsData, oneCell, pointIds, cellsTouchingPointCount)
    elif minCtpc >= 0:
      #outside surface corner
      retStatus = 4
      #retStatus = PhactoriFollowSurfaceSide1.TestOutSideSurfaceCornerCellForMaxAreaFacing(
      #  inInputCsData, oneCell, pointIds, cellsTouchingPointCount)
    else:
      #should not get here
      myDebugPrint3AndException("bad minCtpc error\n")
      retStatus = -1

    global gCellStatusTrack
    gCellStatusTrack[retStatus] += 1
    if PhactoriDbg(100):
      if inCellIndex % 1000 == 0:
        myDebugPrint3("cell index: " + str(inCellIndex) + " status: " + str(retStatus) + "\n")
        for ii, cnt in enumerate(gCellStatusTrack):
          myDebugPrint3(str(ii) + ": " + str(cnt) + "\n")

    return retStatus

  def FindSurfaceCells(surfaceStatus):
    retMap = {}
    firstSurfaceCellIndex = -1
    for oneCellIndex, oneCellStatus in enumerate(surfaceStatus):
      if oneCellStatus > 1:
        retMap[oneCellIndex] = 0
        if firstSurfaceCellIndex < 0:
          if (oneCellStatus == 2) or (oneCellStatus == 5):
            #grab first non-edge non-corner surface cell to use as seed
            firstSurfaceCellIndex = oneCellIndex
    if PhactoriDbg(100):
      myDebugPrint3("FindSurfaceCells retMap:\n")
      for cellIndex, tag in retMap
        myDebugPrint3(str(cellIndex) + " " + str(surfaceStatus[cellIndex]) + "\n")
      myDebugPrint3("FindSurfaceCells num items: " + len(retMap.keys()) + "\n")
      myDebugPrint3("firstSurfaceCellIndex: " + str(firstSurfaceCellIndex) + "\n")
    return retMap, firstSurfaceCellIndex

  def FindWhichSurfaceCellsTouchEachPoint(inInputCsData, surfaceCellMap):
    cellsTouchingPoints = {}
    for cellIndex in surfaceCellMap.keys():
      oneCell = inInputCsData.GetCell(cellIndex)
      numPoints = oneCell.GetNumberOfPoints()
      cellPointIds = oneCell.GetPointIds()
      for ii in range(0, numPoints):
        ptId = cellPointIds.GetId(ii)
        if ptId not in cellsTouchingPoints:
          cellsTouchingPoints[ptId] = []
        cellsTouchingPoints[ptId].append(cellIndex)
    if PhactoriDbg(100):
      myDebugPrint3("cellsTouchingPoints:\n")
      for ptId, cellList in cellsTouchingPoints:
        myDebugPrint(str(ptId) + " " + str(celllist) + "\n")
    return cellsTouchingPoints
    
  def FollowSurfaceSideFromSeedCell(inInputCsData, seedCellIndex, surfaceStatus):
    #get a list of all surface cells
    surfaceCellMap, seedCellIndex = FindSurfaceCells(surfaceStatus)

    #go through points, find cells attached to points
    cellsTouchingPoints = FindWhichSurfaceCellsTouchEachPoint(inInputCsData, surfaceStatus, surfaceCellMap)

    frontierCellList = [seedCellIndex]
    surfaceCellMap[seedCellIndex] = 1
    while len(frontierCellList) > 0:
      testCellIndex = frontierCellList.pop()
      if PhactoriDbg(100):
        myDebugPrint3("frontierCellList:\n" + str(frontierCellList) + "\ndoing cell: " + str(testCellIndex) + "\n")
      oneCell = inInputCsData.GetCell(testCellIndex)
      numPoints = oneCell.GetNumberOfPoints()
      cellPointIds = oneCell.GetPointIds()
      for ii in range(0, numPoints):
        ptId = cellPointIds.GetId(ii)
        cellsTouchingThisPoint = cellsTouchingPoints[ptId]
        if PhactoriDbg(100):
          myDebugPrint3(str(ii) + " of " + str(numPoints) + " id " + str(ptId) + " cttp " + str(cellsTouchingThisPoint) + "\n")
        for oneCellId, visitTag in cellsTouchingThisPoint:
          if oneCellId not in surfaceCellMap:
            if PhactoriDbg(100):
              myDebugPrint3(str(oneCellId) + " not a surface cell\n")
            continue
          if visitTag > 0:
            if PhactoriDbg(100):
              myDebugPrint3(str(oneCellId) + " already visited\n")
            continue
          surfaceCellMap[oneCellId] += 1
          #if this is not edge or corner, add to frontier
          if (surfaceStatus[oneCellId] == 2) or (surfaceStatus[oneCellId] == 5):
            if PhactoriDbg(100):
              myDebugPrint3(str(oneCellId) + " tagged visited, adding to frontier\n")
            frontierCellList.append(oneCellId)
          else:
            if PhactoriDbg(100):
              myDebugPrint3(str(oneCellId) + " tagged visited, corner or edge, not added to frontier\n")

  @staticmethod
  def MarkCellStatusInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("MarkCellStatusInBlock entered\n")

    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    if PhactoriDbg(100):
      myDebugPrint3("numCells: " + str(numCells) + \
        "  numPoints: " + str(numPoints) + "\n")
    inParameters.numCellsPerBlock.append(numCells)
    inParameters.numPointsPerBlock.append(numPoints)

    cellsTouchingPointCount = []
    for ii in range(0, numPoints):
      cellsTouchingPointCount.append(0)

    #see how many cells are touching each point
    for ii in range(0, numCells):
      #if PhactoriDbg(100):
      #  myDebugPrint3("cell index " + str(ii) + " of " + str(numCells) + "\n")
      PhactoriFollowSurfaceSide1.IncrementTouchingCellCountForAllPointsInOneCell(
        inInputCsData, inParameters, ii, cellsTouchingPointCount)

    #see which cells have how many points with 8 or other cells
    surfaceStatus = []
    for ii in range(0, numCells):
      #if PhactoriDbg(100):
      #  myDebugPrint3("(b) cell index " + str(ii) + " of " + str(numCells) + "\n")
      oneCellStatus = PhactoriFollowSurfaceSide1.FindSurfaceStatusForOneCell(
        inInputCsData, inParameters, ii, cellsTouchingPointCount)
      surfaceStatus.append(oneCellStatus)

    seedCellIndex = 0
    surfaceSideCellIndices = FollowSurfaceSideFromSeedCell(inInputCsData, seedCellIndex, surfaceStatus)

    if PhactoriDbg(100):
      myDebugPrint3("MarkCellStatusInBlock returning\n")

  def ExportOperationData(self, datadescription):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFollowSurfaceSide1.ExportOperationData entered\n",
        100)

    UpdatePipelineWithCurrentTimeArgument(self.inPvFilter)

    recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    recursionObj.mParameters = CellEdgeAngleMetricParameters()
    recursionObj.mParameters.Initialize(self.UseSmallestAngle,
      self.OffsetIndex)
    recursionObj.mOperationToDoPerBlock = self.MarkCellStatusInBlock
    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.inPvFilter)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFollowSurfaceSide1.ExportOperationData returning\n", 100)

  def CreateProgrammableFilterString(self):

    if self.mProgFilterString != None:
      return

    scriptLines = []

    scriptLines.append("import math\n")
    scriptLines.append("def IncrementTouchingCellCountForAllPointsInOneCell(inInputCsData, inCellIndex,\n")
    scriptLines.append("      cellsTouchingPointCount):\n")
    scriptLines.append("  oneCell = inInputCsData.GetCell(inCellIndex)\n")
    scriptLines.append("  numPoints = oneCell.GetNumberOfPoints()\n")
    scriptLines.append("  pointIds = oneCell.GetPointIds()\n")
    scriptLines.append("  onePt = [0.0,0.0,0.0]\n")
    scriptLines.append("  for ii in range(0, numPoints):\n")
    scriptLines.append("    onePointId = pointIds.GetId(ii)\n")
    scriptLines.append("    cellsTouchingPointCount[onePointId] += 1\n")
    scriptLines.append("    inInputCsData.GetPoint(onePointId, onePt)\n")
    scriptLines.append("def GetEdgeLength(inInputCsData, oneFace, edgeIndex):\n")
    scriptLines.append("  oneEdge = oneFace.GetEdge(edgeIndex)\n")
    scriptLines.append("  edgePointIds = oneEdge.GetPointIds()\n")
    scriptLines.append("  ptA = inInputCsData.GetPoint(edgePointIds.GetId(0))\n")
    scriptLines.append("  ptB = inInputCsData.GetPoint(edgePointIds.GetId(1))\n")
    scriptLines.append("  xx = ptB[0] - ptA[0]\n")
    scriptLines.append("  yy = ptB[1] - ptA[1]\n")
    scriptLines.append("  zz = ptB[2] - ptA[2]\n")
    scriptLines.append("  return math.sqrt(xx*xx + yy*yy + zz*zz)\n")
    scriptLines.append("def GetFaceSize(inInputCsData, oneFace):\n")
    scriptLines.append("  numFaceEdges = oneFace.GetNumberOfEdges()\n")
    scriptLines.append("  if numFaceEdges != 4:\n")
    scriptLines.append("    retVal = 0\n")
    scriptLines.append("    for ii in range(0, numFaceEdges):\n")
    scriptLines.append("      retVal += GetEdgeLength(inInputCsData, oneFace, ii)\n")
    scriptLines.append("  else:\n")
    scriptLines.append("    side0Len = GetEdgeLength(inInputCsData, oneFace, 0)\n")
    scriptLines.append("    side1Len = GetEdgeLength(inInputCsData, oneFace, 1)\n")
    scriptLines.append("    side2Len = GetEdgeLength(inInputCsData, oneFace, 2)\n")
    scriptLines.append("    side3Len = GetEdgeLength(inInputCsData, oneFace, 3)\n")
    scriptLines.append("    retVal = (side0Len + side2Len) * (side1Len * side3Len)\n")
    scriptLines.append("  return retVal\n")
    scriptLines.append("def GetSortedFaceSizeList(inInputCsData, oneCell):\n")
    scriptLines.append("  faceSizeList = []\n")
    scriptLines.append("  numFaces = oneCell.GetNumberOfFaces()\n")
    scriptLines.append("  for faceNdx in range(0, numFaces):\n")
    scriptLines.append("    oneFace = oneCell.GetFace(faceNdx)\n")
    scriptLines.append("    faceSizeList.append(GetFaceSize(inInputCsData, oneFace))\n")
    scriptLines.append("  return sorted(faceSizeList)\n")
    scriptLines.append("def IsFaceOneOfTwoLargestOnCell(inInputCsData, oneCell, oneFace):\n")
    scriptLines.append("  SortedFaceSizeList = GetSortedFaceSizeList(inInputCsData, oneCell)\n")
    scriptLines.append("  testFaceSize = GetFaceSize(inInputCsData, oneFace)\n")
    scriptLines.append("  if len(SortedFaceSizeList) < 2:\n")
    scriptLines.append("    return False\n")
    scriptLines.append("  return (testFaceSize >= SortedFaceSizeList[-2])\n")
    scriptLines.append("def SurfaceCornerCellFaceFacesExterior(oneFace, cellsTouchingPointCount):\n")
    scriptLines.append("  numFacePoints = oneFace.GetNumberOfPoints()\n")
    scriptLines.append("  facePointIds = oneFace.GetPointIds()\n")
    scriptLines.append("  numOutsideVertices = 0\n")
    scriptLines.append("  for ii in range(0, numFacePoints):\n")
    scriptLines.append("    if cellsTouchingPointCount[ii] <= 0:\n")
    scriptLines.append("      return True\n")
    scriptLines.append("  return False\n")
    scriptLines.append("def TestOutSideSurfaceCornerCellForMaxAreaFacing(\n")
    scriptLines.append("  inInputCsData, oneCell, pointIds, cellsTouchingPointCount):\n")
    scriptLines.append("  numFaces = oneCell.GetNumberOfFaces()\n")
    scriptLines.append("  for faceIndex in range(0, numFaces):\n")
    scriptLines.append("    oneFace = oneCell.GetFace(faceIndex)\n")
    scriptLines.append("    if SurfaceCornerCellFaceFacesExterior(oneFace, cellsTouchingPointCount):\n")
    scriptLines.append("      if IsFaceOneOfTwoLargestOnCell(inInputCsData, oneCell, oneFace):\n")
    scriptLines.append("        return 7\n")
    scriptLines.append("  return 4\n")
    scriptLines.append("def SurfaceEdgeCellFaceFacesExterior(oneFace, cellsTouchingPointCount):\n")
    scriptLines.append("  numFacePoints = oneFace.GetNumberOfPoints()\n")
    scriptLines.append("  facePointIds = oneFace.GetPointIds()\n")
    scriptLines.append("  numOutsideVertices = 0\n")
    scriptLines.append("  for ii in range(0, numFacePoints):\n")
    scriptLines.append("    if cellsTouchingPointCount[ii] <= 2:\n")
    scriptLines.append("      numOutsideVertices += 1\n")
    scriptLines.append("  return (numOutsideVertices >= 2)\n")
    scriptLines.append("def TestOutSideSurfaceEdgeCellForMaxAreaFacing(\n")
    scriptLines.append("  inInputCsData, oneCell, pointIds, cellsTouchingPointCount):\n")
    scriptLines.append("  numFaces = oneCell.GetNumberOfFaces()\n")
    scriptLines.append("  for faceIndex in range(0, numFaces):\n")
    scriptLines.append("    oneFace = oneCell.GetFace(faceIndex)\n")
    scriptLines.append("    if SurfaceEdgeCellFaceFacesExterior(oneFace, cellsTouchingPointCount):\n")
    scriptLines.append("      if IsFaceOneOfTwoLargestOnCell(inInputCsData, oneCell, oneFace):\n")
    scriptLines.append("        return 6\n")
    scriptLines.append("  return 3\n")
    scriptLines.append("def TestOutSideSurfaceCellForMaxAreaFacing(\n")
    scriptLines.append("      inInputCsData, oneCell, pointIds, cellsTouchingPointCount):\n")
    scriptLines.append("  numFaces = oneCell.GetNumberOfFaces()\n")
    scriptLines.append("  for faceIndex in range(0, numFaces):\n")
    scriptLines.append("    oneFace = oneCell.GetFace(faceIndex)\n")
    scriptLines.append("    numFacePoints = oneFace.GetNumberOfPoints()\n")
    scriptLines.append("    facePtNdxs = oneFace.GetPointIds()\n")
    scriptLines.append("    for jj in range(0, numFacePoints):\n")
    scriptLines.append("      ptNdx = facePtNdxs.GetId(jj)\n")
    scriptLines.append("      if cellsTouchingPointCount[ptNdx] != 4:\n")
    scriptLines.append("        break\n")
    scriptLines.append("    if jj >= numFacePoints:\n")
    scriptLines.append("      break\n")
    scriptLines.append("  if faceIndex >= numFaces:\n")
    scriptLines.append("    return -2\n")
    scriptLines.append("  if IsFaceOneOfTwoLargestOnCell(inInputCsData, oneCell, oneFace):\n")
    scriptLines.append("    return 5\n")
    scriptLines.append("  else:\n")
    scriptLines.append("    return 2\n")
    scriptLines.append("def FindSurfaceStatusForOneCell(inInputCsData, inCellIndex,\n")
    scriptLines.append("      cellsTouchingPointCount):\n")
    scriptLines.append("  oneCell = inInputCsData.GetCell(inCellIndex)\n")
    scriptLines.append("  numPoints = oneCell.GetNumberOfPoints()\n")
    scriptLines.append("  pointIds = oneCell.GetPointIds()\n")
    scriptLines.append("  minCtpc = 9\n")
    scriptLines.append("  for ii in range(0, numPoints):\n")
    scriptLines.append("    testCtpc = cellsTouchingPointCount[pointIds.GetId(ii)]\n")
    scriptLines.append("    if testCtpc < minCtpc:\n")
    scriptLines.append("      minCtpc = testCtpc\n")
    scriptLines.append("  if minCtpc >= 8:\n")
    scriptLines.append("    retStatus = 1\n")
    scriptLines.append("  elif minCtpc >= 4:\n")
    scriptLines.append("    retStatus = TestOutSideSurfaceCellForMaxAreaFacing(\n")
    scriptLines.append("      inInputCsData, oneCell, pointIds, cellsTouchingPointCount)\n")
    scriptLines.append("  elif minCtpc >= 2:\n")
    scriptLines.append("    retStatus = TestOutSideSurfaceEdgeCellForMaxAreaFacing(\n")
    scriptLines.append("      inInputCsData, oneCell, pointIds, cellsTouchingPointCount)\n")
    scriptLines.append("  elif minCtpc >= 0:\n")
    scriptLines.append("    retStatus = TestOutSideSurfaceCornerCellForMaxAreaFacing(\n")
    scriptLines.append("      inInputCsData, oneCell, pointIds, cellsTouchingPointCount)\n")
    scriptLines.append("  else:\n")
    scriptLines.append("    retStatus = -1\n")
    scriptLines.append("  return retStatus\n")

    scriptLines.append("localLeafVisitCount = 0\n")
    scriptLines.append("def flatten(input, output):\n")
    scriptLines.append("    # Copy the cells etc.\n")
    scriptLines.append("    output.ShallowCopy(input)\n")
    scriptLines.append("    numPoints = input.GetNumberOfPoints()\n")
    scriptLines.append("    celldata = output.GetCellData()\n")
    scriptLines.append("    numCells = input.GetNumberOfCells()\n")
    scriptLines.append("    ncda = vtk.vtkIntArray()\n")
    scriptLines.append("    ncda.SetNumberOfTuples(numCells)\n")
    scriptLines.append("    numTuples = ncda.GetNumberOfTuples()\n")
    scriptLines.append("    for ii in range(0, numTuples):\n")
    scriptLines.append("      ncda.SetValue(ii, 0)\n")
    scriptLines.append("    cellsTouchingPointCount = []\n")
    scriptLines.append("    for ii in range(0, numPoints):\n")
    scriptLines.append("      cellsTouchingPointCount.append(0)\n")
    scriptLines.append("    for ii in range(0, numCells):\n")
    scriptLines.append("      IncrementTouchingCellCountForAllPointsInOneCell(input, ii, cellsTouchingPointCount)\n")
    scriptLines.append("    for ii in range(0, numCells):\n")
    scriptLines.append("      oneCellSurfaceStatus = FindSurfaceStatusForOneCell(input, ii, cellsTouchingPointCount)\n")
    scriptLines.append("      ncda.SetValue(ii, oneCellSurfaceStatus)\n")
    scriptLines.append("    ncda.SetName('" + self.ProgrammableFilterOutputCellVariableName + "')\n")
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

    self.mProgFilterString = "".join(scriptLines)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFollowSurfaceSide1 constructed script:\n" + \
        self.mProgFilterString)

#phactori_combine_to_single_python_file_subpiece_end_1
