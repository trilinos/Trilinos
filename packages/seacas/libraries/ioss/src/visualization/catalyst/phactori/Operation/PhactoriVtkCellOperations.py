# Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
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
import math
from .PhactoriVectorLibrary import *
import vtk

#phactori_combine_to_single_python_file_subpiece_begin_1

def PhactoriCalculateFaceArea(inInputCsData, targetFace):
  numFaceEdges = targetFace.GetNumberOfEdges()
  if numFaceEdges >= 3:
    #assumes convex polygons, weird answer if not
    retVal = 0.0
    ptIds = targetFace.GetPointIds()
    ptA = inInputCsData.GetPoint(ptIds.GetId(0))
    ptB = inInputCsData.GetPoint(ptIds.GetId(1))
    vecAB = vecFromAToB(ptA, ptB)
    for ii in range(2, numFaceEdges):
      ptC = inInputCsData.GetPoint(ptIds.GetId(ii))
      vecAC = vecFromAToB(ptA, ptC)
      ABxAC = vecCrossProduct(vecAB, vecAC)
      retVal += vecMagnitude(ABxAC)
      vecAB = vecAC
    retVal *= 0.5
  else:
    #two or fewer vertices
    retVal = 0.0
  return retVal

def PhactoriCalculateCellFaceAreas(inInputCsData, oneCell):
  faceAreaList = []
  numFaces = oneCell.GetNumberOfFaces()
  for faceNdx in range(0, numFaces):
    gsfslTestFace = oneCell.GetFace(faceNdx)
    faceAreaList.append(PhactoriCalculateFaceArea(inInputCsData, gsfslTestFace))
  return faceAreaList

def PhactoriCalculateCellFaceAreasAndSort(inInputCsData, oneCell):
  faceAreasByIndex = PhactoriCalculateCellFaceAreas(inInputCsData, oneCell)
  faceAreasSorted = sorted(faceAreasByIndex)
  return faceAreasByIndex, faceAreasSorted

class PhactoriCountCellsTouchingFaceData:
  def __init__(self, inInputCsData):
    self.PerFaceData = {}
    self.vtkGrid = inInputCsData
    self.totalNumFaces = 0
    self.totalNumFacesWithOneCellTouching = 0
    self.totalNumFacesWithTwoCellsTouching = 0
    self.totalNumFacesWith3plCellsTouching = 0

  def DoCountForAllCellsInBlock(self):
    numCells = self.vtkGrid.GetNumberOfCells()
    for ii in range(0, numCells):
      self.DoCountForOneCell(ii)

  def DoCountForOneCell(self, cellIndex):
    oneCell = self.vtkGrid.GetCell(cellIndex)
    numFaces = oneCell.GetNumberOfFaces()
    for jj in range(0, numFaces):
      self.DoCountForOneFace(oneCell, cellIndex, jj)

  def GetFaceHash(self, oneFace):
    #for now we use the average of all face points to uniquely identify the
    #face; we could use all indices of the face but that is more complex
    #and we couldn't match up faces across block/processor boundaries
    #however, it's possible that average of vertices may fail in some
    #circumstances and we could use the x/y/z of all the vertices if necessary
    faceHash = [0.0, 0.0, 0.0]
    numFacePoints = oneFace.GetNumberOfPoints()
    facePointIds = oneFace.GetPointIds()
    facePointIds2 = []
    for kk in range(0, numFacePoints):
      facePointIds2.append(facePointIds.GetId(kk))
    facePointIds2 = sorted(facePointIds2)
    for pointId in facePointIds2:
      oneFacePoint = self.vtkGrid.GetPoint(pointId)
      faceHash[0] += oneFacePoint[0]
      faceHash[1] += oneFacePoint[1]
      faceHash[2] += oneFacePoint[2]
    return faceHash

  def DoCountForOneFace(self, oneCell, cellIndex, faceIndex):
    oneFace = oneCell.GetFace(faceIndex)
    faceHash = self.GetFaceHash(oneFace)
    self.IncrementCellTouchingFaceCount(faceHash, cellIndex, faceIndex)

  def IncrementCellTouchingFaceCount(self, faceHash, cellIndex, faceIndex):

    if faceHash[0] not in self.PerFaceData:
      self.PerFaceData[faceHash[0]] = {}
    facelvl1 = self.PerFaceData[faceHash[0]]
    if faceHash[1] not in facelvl1:
      facelvl1[faceHash[1]] = {}
    facelvl2 = facelvl1[faceHash[1]]
    if faceHash[2] not in facelvl2:
      facelvl2[faceHash[2]] = [1, cellIndex, faceIndex, -1, -1]
      self.totalNumFaces += 1
      self.totalNumFacesWithOneCellTouching += 1
    else:
      facelvl3 = facelvl2[faceHash[2]]
      if facelvl3[0] == 1:
        self.totalNumFacesWithOneCellTouching -= 1
        self.totalNumFacesWithTwoCellsTouching += 1
      elif facelvl3[0] == 2:
        self.totalNumFacesWithTwoCellsTouching -= 1
        self.totalNumFacesWith3plCellsTouching += 1
      facelvl3[0] += 1
      facelvl3[3] = cellIndex
      facelvl3[4] = faceIndex

  def GetNumberOfCellsTouchingFace(self, oneFace):
    faceHash = self.GetFaceHash(oneFace)
    returnCount = self.PerFaceData[faceHash[0]][faceHash[1]][faceHash[2]][0]
    return returnCount

  def GetNumberOfCellsTouchingEachFaceOfOneCell(self, oneCellIndex):
    oneCell = self.vtkGrid.GetCell(oneCellIndex)
    returnCountPerFace = []
    numFaces = oneCell.GetNumberOfFaces()
    for ii in range(0, numFaces):
      oneFace = oneCell.GetFace(ii)
      returnCountPerFace.append(self.GetNumberOfCellsTouchingFace(oneFace))
    return returnCountPerFace

  def GetListOfExteriorFacesOnCell(self, oneCell):
    returnExteriorFaceIndices = []
    numFaces = oneCell.GetNumberOfFaces()
    for ii in range(0, numFaces):
      oneFace = oneCell.GetFace(ii)
      if self.GetNumberOfCellsTouchingFace(oneFace) == 1:
        returnExteriorFaceIndices.append(ii)
    return returnExteriorFaceIndices

  def GetListOfExteriorFacesOnCellByIndex(self, cellIndex):
    return self.GetListOfExteriorFacesOnCell(self.vtkGrid.GetCell(cellIndex))

def PhactoriCountCellTouchingEachFace(inInputCsData):
  """goes through each cell, and for each cell goes through each face, and for
     each of those faces creates/locates a count for that face and increments
     it. Purpose is to find cells with faces that don't neighbor other cells
  """
  countCellsTouchingEachFace = PhactoriCountCellsTouchingFaceData(inInputCsData)
  countCellsTouchingEachFace.DoCountForAllCellsInBlock()
  if PhactoriDbg(100):
    cctef = countCellsTouchingEachFace
    myDebugPrint3("PhactoriCountCellTouchingEachFace:" + \
      "\nnum cells in block:      " + str(inInputCsData.GetNumberOfCells()) + \
      "\nnum faces:               " + str(cctef.totalNumFaces) + \
      "\nnum faces with 1 cell:   " + str(cctef.totalNumFacesWithOneCellTouching) + \
      "\nnum faces with 2 cells:  " + str(cctef.totalNumFacesWithTwoCellsTouching) + \
      "\nnum faces with 3+ cells: " + str(cctef.totalNumFacesWith3plCellsTouching) + "\n")

  return countCellsTouchingEachFace

def PhactoriCountCellTouchingEachPoint(inInputCsData):
  numCells = inInputCsData.GetNumberOfCells()
  numPoints = inInputCsData.GetNumberOfPoints()
  cellsTouchingPointCount = [0] * numPoints
  for ii in range(0, numCells):
    oneCell = inInputCsData.GetCell(ii)
    cellPointIds = oneCell.GetPointIds()
    numCellPoints = oneCell.GetNumberOfPoints()
    for jj in range(0, numCellPoints):
      oneCellPointId = cellPointIds.GetId(jj)
      cellsTouchingPointCount[oneCellPointId] += 1
  return cellsTouchingPointCount

def PhactoriSurfaceCornerCellFaceFacesExterior(oneFace, cellsTouchingPointCount):
  numFacePoints = oneFace.GetNumberOfPoints()
  facePointIds = oneFace.GetPointIds()
  numOutsideVertices = 0
  for ii in range(0, numFacePoints):
    if cellsTouchingPointCount[facePointIds.GetId(ii)] <= 1:
      return True
  return False

def PhactoriIsOutsideCornerHexCellFaceMaxAreaFace(oneCell, indexFaceListSize,
  sortedFaceSizeList, cellsTouchingPointCount, TestNthLargestFace = 2):
  #faces with 4 2 2 0 need to be tested
  numFaces = oneCell.GetNumberOfFaces()
  for faceIndex in range(0, numFaces):
    oneFace = oneCell.GetFace(faceIndex)
    if PhactoriSurfaceCornerCellFaceFacesExterior(oneFace, cellsTouchingPointCount):
      if indexFaceListSize[faceIndex] >= sortedFaceSizeList[-TestNthLargestFace]:
        #this cell is facing the big area cell side
        return True
  return False

def PhactoriSurfaceEdgeHexCellFaceFacesExterior(oneFace, cellsTouchingPointCount):
  numFacePoints = oneFace.GetNumberOfPoints()
  facePointIds = oneFace.GetPointIds()
  numOutsideVertices = 0
  for ii in range(0, numFacePoints):
    if cellsTouchingPointCount[facePointIds.GetId(ii)] <= 2:
      numOutsideVertices += 1
  return (numOutsideVertices >= 2)

def PhactoriIsOutsideEdgeHexCellFaceMaxAreaFace(oneCell,
  indexFaceListSize, sortedFaceSizeList, cellsTouchingPointCount,
  TestNthLargestFace = 2):
  #faces with 4 4 2 2 need to be tested
  numFaces = oneCell.GetNumberOfFaces()
  #test cell flatness
  for faceIndex in range(0, numFaces):
    oneFace = oneCell.GetFace(faceIndex)
    if PhactoriSurfaceEdgeHexCellFaceFacesExterior(oneFace, cellsTouchingPointCount):
      if indexFaceListSize[faceIndex] >= sortedFaceSizeList[-TestNthLargestFace]:
        #this cell is exterior and facing the big area cell side
        return True
  return False

def PhactoriIsOutsideSurfaceHexCellFaceMaxAreaFace(oneCell,
  indexFaceListSize, sortedFaceSizeList, cellsTouchingPointCount,
  TestNthLargestFace = 2):
  #find face which face faces exterior
  numFaces = oneCell.GetNumberOfFaces()
  for faceIndex in range(0, numFaces):
    oneFace = oneCell.GetFace(faceIndex)
    numFacePoints = oneFace.GetNumberOfPoints()
    facePtNdxs = oneFace.GetPointIds()
    outsideFaceIndex = -1
    for jj in range(0, numFacePoints):
      outsideFaceIndex = faceIndex
      ptNdx = facePtNdxs.GetId(jj)
      if cellsTouchingPointCount[ptNdx] != 4:
        outsideFaceIndex = -1
        #this is not the outward facing face
        break
    if outsideFaceIndex >= 0:
      #this is the outward facing face
      break
  if outsideFaceIndex < 0:
    #shouldn't happen, didn't find face with all 4's
    myDebugPrint3AndException("outsideFaceIndex < 0 error\n")
    return False
  #now see if this outward facing face is one of the Nth largest
  if indexFaceListSize[outsideFaceIndex] >= sortedFaceSizeList[-TestNthLargestFace]:
    return True
  else:
    return False

def PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(inInputCsData,
  inCellIndex, countCellsTouchingEachFace, flatnessTestRatio = 100.0,
  TestNthLargestFace = 2):
  oneCell = inInputCsData.GetCell(inCellIndex)
  exteriorFaceList = countCellsTouchingEachFace.GetListOfExteriorFacesOnCell(oneCell)
  numExteriorFaces = len(exteriorFaceList)
  if numExteriorFaces == 0:
    #interior cell
    retStatus = 1
  else:
    #surface cell (surface, edge, or corner)
    #test cell flatness
    indexFaceListSize, sortedFaceSizeList = PhactoriCalculateCellFaceAreasAndSort(inInputCsData, oneCell)
    flatnessRatio = sortedFaceSizeList[-1] / sortedFaceSizeList[0]
    if len(sortedFaceSizeList) <= 0:
      retStatus = 11
      return retStatus
    cellIsFlat = (flatnessRatio >= flatnessTestRatio)
    #see if any of the exterior faces are the biggest
    outsideFaceIsMaxAreaFace = False
    for outsideFaceIndex in exteriorFaceList:
      if indexFaceListSize[outsideFaceIndex] >= sortedFaceSizeList[-TestNthLargestFace]:
        outsideFaceIsMaxAreaFace = True
        break
    if numExteriorFaces == 1:
      #outside surface, not edge or corner
      if outsideFaceIsMaxAreaFace:
        if cellIsFlat:
          retStatus = 11
        else:
          retStatus = 8
      else:
        if cellIsFlat:
          retStatus = 5
        else:
          retStatus = 2
    elif numExteriorFaces == 2:
      #outside edge, not corner
      if outsideFaceIsMaxAreaFace:
        if cellIsFlat:
          retStatus = 12
        else:
          retStatus = 9
      else:
        if cellIsFlat:
          retStatus = 6
        else:
          retStatus = 3
    else:
      #outside surface corner
      if outsideFaceIsMaxAreaFace:
        if cellIsFlat:
          retStatus = 13
        else:
          retStatus = 10
      else:
        if cellIsFlat:
          retStatus = 7
        else:
          retStatus = 4
  return retStatus


def PhactoriFindSurfaceStatusForOneCell(inInputCsData, inCellIndex,
      cellsTouchingPointCount, flatnessTestRatio = 100.0):
  oneCell = inInputCsData.GetCell(inCellIndex)

  #cells with all 8 are interior
  #cells with 4 and 8 are outside surface
  #cells with 2 (and no 0) are outside edge
  #cells with 1 are corner

  #not currently handling inside edge (there are count 6 cells) or inside
  #corner (there are count 7 cells)

  #if cell is outside surface (4 and 8), if the face which is all 4 is the
  #largest or second largest area, then it is on the surface we want

  #if the cell is edge (2 4 8), there should be two sides which are 2 2 4 4
  #if either of these sides is the second or largest area then it is on the
  #surface we want.

  #if the cell is corner, there should be three sides which are 0 2 2 4
  #if any of these sides is the second or largest area then it is on the
  #surface we want

  numPoints = oneCell.GetNumberOfPoints()
  pointIds = oneCell.GetPointIds()
  minCtpc = 9
  for ii in range(0, numPoints):
    testCtpc = cellsTouchingPointCount[pointIds.GetId(ii)]
    if testCtpc < minCtpc:
      minCtpc = testCtpc

  #if inCellIndex == 172065:
  if inCellIndex == 5056:
    debugFlag = 100
  else:
    debugFlag = 25
  if minCtpc >= 8:
    #interior cell
    retStatus = 1
  else:
    #surface cell (surface, edge, or corner)
    #test cell flatness
    indexFaceListSize, sortedFaceSizeList = PhactoriCalculateCellFaceAreasAndSort(inInputCsData, oneCell)
    flatnessRatio = sortedFaceSizeList[-1] / sortedFaceSizeList[0]
    if len(sortedFaceSizeList) <= 0:
      retStatus = 11
      return retStatus
    cellIsFlat = (flatnessRatio >= flatnessTestRatio)
    if minCtpc >= 4:
      #outside surface, not edge or corner
      outsideFaceIsMaxAreaFace = PhactoriIsOutsideSurfaceHexCellFaceMaxAreaFace(
        oneCell, indexFaceListSize, sortedFaceSizeList,
        cellsTouchingPointCount)
      if outsideFaceIsMaxAreaFace:
        if cellIsFlat:
          retStatus = 11
        else:
          retStatus = 8
      else:
        if cellIsFlat:
          retStatus = 5
        else:
          retStatus = 2
    elif minCtpc >= 2:
      #outside surface edge, not corner
      outsideFaceIsMaxAreaFace = PhactoriIsOutsideEdgeHexCellFaceMaxAreaFace(
        oneCell, indexFaceListSize, sortedFaceSizeList,
        cellsTouchingPointCount)
      if outsideFaceIsMaxAreaFace:
        if cellIsFlat:
          retStatus = 12
        else:
          retStatus = 9
      else:
        if cellIsFlat:
          retStatus = 6
        else:
          retStatus = 3
    elif minCtpc >= 1:
      #outside surface corner
      outsideFaceIsMaxAreaFace = PhactoriIsOutsideCornerHexCellFaceMaxAreaFace(
        oneCell, indexFaceListSize, sortedFaceSizeList, cellsTouchingPointCount)
      if outsideFaceIsMaxAreaFace:
        if cellIsFlat:
          retStatus = 13
        else:
          retStatus = 10
      else:
        if cellIsFlat:
          retStatus = 7
        else:
          retStatus = 4
    else:
      #should not get here
      myDebugPrint3AndException("bad minCtpc error\n")
      retStatus = -1

  return retStatus

def PhactoriGetCellEdgeVector(inInputCsData, oneFaceOrCell, edgeIndex):
  oneFaceEdgePointIds = oneFaceOrCell.GetEdge(edgeIndex).GetPointIds()
  ptA = inInputCsData.GetPoint(oneFaceEdgePointIds.GetId(0))
  ptB = inInputCsData.GetPoint(oneFaceEdgePointIds.GetId(1))
  retVec = vecFromAToB(ptA,ptB)
  return retVec

def PhactoriFindFaceNormal(inInputCsData, oneFace):
  facePointIds = oneFace.GetPointIds()
  #maybe should work harder to find good normal
  ptA = inInputCsData.GetPoint(facePointIds.GetId(0))
  ptB = inInputCsData.GetPoint(facePointIds.GetId(1))
  ptC = inInputCsData.GetPoint(facePointIds.GetId(2))
  vecAB = vecFromAToB(ptA, ptB)
  vecBC = vecFromAToB(ptB, ptC)
  retNormalVec = vecCrossProduct(vecBC, vecAB)
  vecNormalize2(retNormalVec, retNormalVec)
  return retNormalVec

def PhactoriFindLargestCellFaceNormal(inInputCsData, oneCell):
  indexFaceListSize, sortedFaceSizeList = PhactoriCalculateCellFaceAreasAndSort(inInputCsData, oneCell)
  biggestFaceSize = sortedFaceSizeList[-1]
  for faceForNormalIndex in range(0, len(indexFaceListSize)):
    if indexFaceListSize[faceForNormalIndex] == biggestFaceSize:
      break
  faceForNormal = oneCell.GetFace(faceForNormalIndex)
  retNormalVec = PhactoriFindFaceNormal(inInputCsData, faceForNormal)
  return retNormalVec

def PhactoriGetAverageOfCellPoints(inInputCsData, oneCell):
  numPoints = oneCell.GetNumberOfPoints()
  retAvgPt = [0.0, 0.0, 0.0]
  if numPoints <= 0:
    return retAvgPt
  ptIds = oneCell.GetPointIds()
  for ii in range(0, numPoints):
    ptA = inInputCsData.GetPoint(ptIds.GetId(ii))
    vecAddInPlace(retAvgPt, ptA)
  avgFac = 1.0/float(numPoints)
  vecScaleInPlace(avgFac, retAvgPt)
  return retAvgPt

def PhactoriFindSelectedAngleBetweenEdgeAndCellNormal(inInputCsData, oneCell,
  compareNormal, paramUseSmallestAngle, paramOffsetIndex):
  edgeCompareDotProdList = []
  numEdges = oneCell.GetNumberOfEdges()
  for ii in range(0, numEdges):
    oneEdgeVec = PhactoriGetCellEdgeVector(inInputCsData, oneCell, ii)
    vecNormalize2(oneEdgeVec, oneEdgeVec)
    edgeCompareDotProd = abs(vecDotProduct(compareNormal, oneEdgeVec))
    edgeCompareDotProdList.append(edgeCompareDotProd)

  sortedEdgeCompareDotProdList = sorted(edgeCompareDotProdList)

  if paramUseSmallestAngle:
    outputIndex = max(0, numEdges - 1 - paramOffsetIndex)
  else:
    outputIndex = min(paramOffsetIndex, numEdges - 1)

  outputAngle = 90.0 - math.degrees(math.acos(sortedEdgeCompareDotProdList[outputIndex]))
  return outputAngle

def DebugPrintPhactoriFindHeightOfCell(inInputCsData, oneCell, compareNormal):
  if PhactoriDbg(100):
    dbght = []
    dbglen = []
    dbgretHeight = -1.0
    dbgretHeightIndex = -1
    for ii in range(0, numCellEdges):
      oneCellEdgeVec = PhactoriGetCellEdgeVector(inInputCsData, oneCell, ii)
      dp1 = vecDotProduct(oneCellEdgeVec, compareNormal)
      projectedVec = vecScale(dp1, compareNormal)
      oneCellEdgeHeight = vecMagnitude(projectedVec)
      dbght.append(oneCellEdgeHeight)
      dbglen.append(vecMagnitude(oneCellEdgeVec))
      if oneCellEdgeHeight > dbgretHeight:
        dbgretHeight = oneCellEdgeHeight
        dbgretHeightIndex = ii
    #how many edges before projection are longer than the height edge, and
    #is the height edge one of the four shortest
    myDebugPrint3("cell height: " + str(dbgretHeight) + "  index: " + \
      str(dbgretHeightIndex) + "  edgelen: " + str(dbglen[dbgretHeightIndex]) + "\n")
    myDebugPrint3("all edge length, height pairs:\n")
    numLonger = 0
    for ii in range(0, numCellEdges):
      if dbglen[ii] >= dbglen[dbgretHeightIndex]:
        numLonger += 1
      myDebugPrint3(str(ii) + ": " + str(dbglen[ii]) + ", " + str(dbght[ii]) + "\n")
    myDebugPrint3("num edges longer than cell height edge: " + str(numLonger) + "\n")

def PhactoriFindHeightOfCell(inInputCsData, oneCell, compareNormal):
  #do we need to discard eight longest sides, or will they project to less
  #than the four shortest sides?
  numCellEdges = oneCell.GetNumberOfEdges()

  retHeight = -1.0
  for ii in range(0, numCellEdges):
    #get edge vectory and project on to normal vector (which is unit length)
    #and find magnitude
    oneCellEdgeVec = PhactoriGetCellEdgeVector(inInputCsData, oneCell, ii)
    dp1 = vecDotProduct(oneCellEdgeVec, compareNormal)
    projectedVec = vecScale(dp1, compareNormal)
    oneCellEdgeHeight = vecMagnitude(projectedVec)
    if oneCellEdgeHeight > retHeight:
      retHeight = oneCellEdgeHeight
  #if PhactoriDbg(100):
  #  DebugPrintPhactoriFindHeightOfCell(inInputCsData, oneCell, compareNormal)
  return retHeight

def PhactoriFindCellEdgeAngleMetricsForOneCell(inInputCsData, inCellIndex,
      paramUseSmallestAngle, paramOffsetIndex):
  oneCell = inInputCsData.GetCell(inCellIndex)
  numEdges = oneCell.GetNumberOfEdges()
  if numEdges < 6:
    return 90.0, 0.0
  if oneCell.GetNumberOfFaces() < 4:
    return 90.0, 0.0

  compareNormal = PhactoriFindLargestCellFaceNormal(inInputCsData, oneCell)
  retAngle = PhactoriFindSelectedAngleBetweenEdgeAndCellNormal(inInputCsData,
    oneCell, compareNormal, paramUseSmallestAngle, paramOffsetIndex)
  retHeight = PhactoriFindHeightOfCell(inInputCsData, oneCell, compareNormal)
  return retAngle, retHeight

def GetPhactoriVtkCellOperationsProgrammableFilterLines(pfLns):
  pfLns.append("def PhactoriCalculateFaceArea(inInputCsData, targetFace):\n")
  pfLns.append("  numFaceEdges = targetFace.GetNumberOfEdges()\n")
  pfLns.append("  if numFaceEdges >= 3:\n")
  pfLns.append("    retVal = 0.0\n")
  pfLns.append("    ptIds = targetFace.GetPointIds()\n")
  pfLns.append("    ptA = inInputCsData.GetPoint(ptIds.GetId(0))\n")
  pfLns.append("    ptB = inInputCsData.GetPoint(ptIds.GetId(1))\n")
  pfLns.append("    vecAB = vecFromAToB(ptA, ptB)\n")
  pfLns.append("    for ii in range(2, numFaceEdges):\n")
  pfLns.append("      ptC = inInputCsData.GetPoint(ptIds.GetId(ii))\n")
  pfLns.append("      vecAC = vecFromAToB(ptA, ptC)\n")
  pfLns.append("      ABxAC = vecCrossProduct(vecAB, vecAC)\n")
  pfLns.append("      retVal += vecMagnitude(ABxAC)\n")
  pfLns.append("      vecAB = vecAC\n")
  pfLns.append("    retVal *= 0.5\n")
  pfLns.append("  else:\n")
  pfLns.append("    retVal = 0.0\n")
  pfLns.append("  return retVal\n")
  pfLns.append("def PhactoriCalculateCellFaceAreas(inInputCsData, oneCell):\n")
  pfLns.append("  faceAreaList = []\n")
  pfLns.append("  numFaces = oneCell.GetNumberOfFaces()\n")
  pfLns.append("  for faceNdx in range(0, numFaces):\n")
  pfLns.append("    gsfslTestFace = oneCell.GetFace(faceNdx)\n")
  pfLns.append("    faceAreaList.append(PhactoriCalculateFaceArea(inInputCsData, gsfslTestFace))\n")
  pfLns.append("  return faceAreaList\n")
  pfLns.append("def PhactoriCalculateCellFaceAreasAndSort(inInputCsData, oneCell):\n")
  pfLns.append("  faceAreasByIndex = PhactoriCalculateCellFaceAreas(inInputCsData, oneCell)\n")
  pfLns.append("  faceAreasSorted = sorted(faceAreasByIndex)\n")
  pfLns.append("  return faceAreasByIndex, faceAreasSorted\n")
  pfLns.append("class PhactoriCountCellsTouchingFaceData:\n")
  pfLns.append("  def __init__(self, inInputCsData):\n")
  pfLns.append("    self.PerFaceData = {}\n")
  pfLns.append("    self.vtkGrid = inInputCsData\n")
  pfLns.append("    self.totalNumFaces = 0\n")
  pfLns.append("    self.totalNumFacesWithOneCellTouching = 0\n")
  pfLns.append("    self.totalNumFacesWithTwoCellsTouching = 0\n")
  pfLns.append("    self.totalNumFacesWith3plCellsTouching = 0\n")
  pfLns.append("  def DoCountForAllCellsInBlock(self):\n")
  pfLns.append("    numCells = self.vtkGrid.GetNumberOfCells()\n")
  pfLns.append("    for ii in range(0, numCells):\n")
  pfLns.append("      self.DoCountForOneCell(ii)\n")
  pfLns.append("  def DoCountForOneCell(self, cellIndex):\n")
  pfLns.append("    oneCell = self.vtkGrid.GetCell(cellIndex)\n")
  pfLns.append("    numFaces = oneCell.GetNumberOfFaces()\n")
  pfLns.append("    for jj in range(0, numFaces):\n")
  pfLns.append("      self.DoCountForOneFace(oneCell, cellIndex, jj)\n")
  pfLns.append("  def GetFaceHash(self, oneFace):\n")
  pfLns.append("    faceHash = [0.0, 0.0, 0.0]\n")
  pfLns.append("    numFacePoints = oneFace.GetNumberOfPoints()\n")
  pfLns.append("    facePointIds = oneFace.GetPointIds()\n")
  pfLns.append("    facePointIds2 = []\n")
  pfLns.append("    for kk in range(0, numFacePoints):\n")
  pfLns.append("      facePointIds2.append(facePointIds.GetId(kk))\n")
  pfLns.append("    facePointIds2 = sorted(facePointIds2)\n")
  pfLns.append("    for pointId in facePointIds2:\n")
  pfLns.append("      oneFacePoint = self.vtkGrid.GetPoint(pointId)\n")
  pfLns.append("      faceHash[0] += oneFacePoint[0]\n")
  pfLns.append("      faceHash[1] += oneFacePoint[1]\n")
  pfLns.append("      faceHash[2] += oneFacePoint[2]\n")
  pfLns.append("    return faceHash\n")
  pfLns.append("  def DoCountForOneFace(self, oneCell, cellIndex, faceIndex):\n")
  pfLns.append("    oneFace = oneCell.GetFace(faceIndex)\n")
  pfLns.append("    faceHash = self.GetFaceHash(oneFace)\n")
  pfLns.append("    self.IncrementCellTouchingFaceCount(faceHash, cellIndex, faceIndex)\n")
  pfLns.append("  def IncrementCellTouchingFaceCount(self, faceHash, cellIndex, faceIndex):\n")
  pfLns.append("    if faceHash[0] not in self.PerFaceData:\n")
  pfLns.append("      self.PerFaceData[faceHash[0]] = {}\n")
  pfLns.append("    facelvl1 = self.PerFaceData[faceHash[0]]\n")
  pfLns.append("    if faceHash[1] not in facelvl1:\n")
  pfLns.append("      facelvl1[faceHash[1]] = {}\n")
  pfLns.append("    facelvl2 = facelvl1[faceHash[1]]\n")
  pfLns.append("    if faceHash[2] not in facelvl2:\n")
  pfLns.append("      facelvl2[faceHash[2]] = [1, cellIndex, faceIndex, -1, -1]\n")
  pfLns.append("      self.totalNumFaces += 1\n")
  pfLns.append("      self.totalNumFacesWithOneCellTouching += 1\n")
  pfLns.append("    else:\n")
  pfLns.append("      facelvl3 = facelvl2[faceHash[2]]\n")
  pfLns.append("      if facelvl3[0] == 1:\n")
  pfLns.append("        self.totalNumFacesWithOneCellTouching -= 1\n")
  pfLns.append("        self.totalNumFacesWithTwoCellsTouching += 1\n")
  pfLns.append("      elif facelvl3[0] == 2:\n")
  pfLns.append("        self.totalNumFacesWithTwoCellsTouching -= 1\n")
  pfLns.append("        self.totalNumFacesWith3plCellsTouching += 1\n")
  pfLns.append("      facelvl3[0] += 1\n")
  pfLns.append("      facelvl3[3] = cellIndex\n")
  pfLns.append("      facelvl3[4] = faceIndex\n")
  pfLns.append("  def GetNumberOfCellsTouchingFace(self, oneFace):\n")
  pfLns.append("    faceHash = self.GetFaceHash(oneFace)\n")
  pfLns.append("    returnCount = self.PerFaceData[faceHash[0]][faceHash[1]][faceHash[2]][0]\n")
  pfLns.append("    return returnCount\n")
  pfLns.append("  def GetNumberOfCellsTouchingEachFaceOfOneCell(self, oneCellIndex):\n")
  pfLns.append("    oneCell = self.vtkGrid.GetCell(oneCellIndex)\n")
  pfLns.append("    returnCountPerFace = []\n")
  pfLns.append("    numFaces = oneCell.GetNumberOfFaces()\n")
  pfLns.append("    for ii in range(0, numFaces):\n")
  pfLns.append("      oneFace = oneCell.GetFace(ii)\n")
  pfLns.append("      returnCountPerFace.append(self.GetNumberOfCellsTouchingFace(oneFace))\n")
  pfLns.append("    return returnCountPerFace\n")
  pfLns.append("  def GetListOfExteriorFacesOnCell(self, oneCell):\n")
  pfLns.append("    returnExteriorFaceIndices = []\n")
  pfLns.append("    numFaces = oneCell.GetNumberOfFaces()\n")
  pfLns.append("    for ii in range(0, numFaces):\n")
  pfLns.append("      oneFace = oneCell.GetFace(ii)\n")
  pfLns.append("      if self.GetNumberOfCellsTouchingFace(oneFace) == 1:\n")
  pfLns.append("        returnExteriorFaceIndices.append(ii)\n")
  pfLns.append("    return returnExteriorFaceIndices\n")
  pfLns.append("  def GetListOfExteriorFacesOnCellByIndex(self, cellIndex):\n")
  pfLns.append("    return self.GetListOfExteriorFacesOnCell(self.vtkGrid.GetCell(cellIndex))\n")
  pfLns.append("def PhactoriCountCellTouchingEachFace(inInputCsData):\n")
  pfLns.append("  countCellsTouchingEachFace = PhactoriCountCellsTouchingFaceData(inInputCsData)\n")
  pfLns.append("  countCellsTouchingEachFace.DoCountForAllCellsInBlock()\n")
  pfLns.append("  return countCellsTouchingEachFace\n")
  pfLns.append("def PhactoriCountCellTouchingEachPoint(inInputCsData):\n")
  pfLns.append("  numCells = inInputCsData.GetNumberOfCells()\n")
  pfLns.append("  numPoints = inInputCsData.GetNumberOfPoints()\n")
  pfLns.append("  cellsTouchingPointCount = [0] * numPoints\n")
  pfLns.append("  for ii in range(0, numCells):\n")
  pfLns.append("    oneCell = inInputCsData.GetCell(ii)\n")
  pfLns.append("    cellPointIds = oneCell.GetPointIds()\n")
  pfLns.append("    numCellPoints = oneCell.GetNumberOfPoints()\n")
  pfLns.append("    for jj in range(0, numCellPoints):\n")
  pfLns.append("      oneCellPointId = cellPointIds.GetId(jj)\n")
  pfLns.append("      cellsTouchingPointCount[oneCellPointId] += 1\n")
  pfLns.append("  return cellsTouchingPointCount\n")
  pfLns.append("def PhactoriSurfaceCornerCellFaceFacesExterior(oneFace, cellsTouchingPointCount):\n")
  pfLns.append("  numFacePoints = oneFace.GetNumberOfPoints()\n")
  pfLns.append("  facePointIds = oneFace.GetPointIds()\n")
  pfLns.append("  numOutsideVertices = 0\n")
  pfLns.append("  for ii in range(0, numFacePoints):\n")
  pfLns.append("    if cellsTouchingPointCount[facePointIds.GetId(ii)] <= 1:\n")
  pfLns.append("      return True\n")
  pfLns.append("  return False\n")
  pfLns.append("def PhactoriIsOutsideCornerHexCellFaceMaxAreaFace(oneCell, indexFaceListSize,\n")
  pfLns.append("  sortedFaceSizeList, cellsTouchingPointCount, TestNthLargestFace = 2):\n")
  pfLns.append("  numFaces = oneCell.GetNumberOfFaces()\n")
  pfLns.append("  for faceIndex in range(0, numFaces):\n")
  pfLns.append("    oneFace = oneCell.GetFace(faceIndex)\n")
  pfLns.append("    if PhactoriSurfaceCornerCellFaceFacesExterior(oneFace, cellsTouchingPointCount):\n")
  pfLns.append("      if indexFaceListSize[faceIndex] >= sortedFaceSizeList[-TestNthLargestFace]:\n")
  pfLns.append("        return True\n")
  pfLns.append("  return False\n")
  pfLns.append("def PhactoriSurfaceEdgeHexCellFaceFacesExterior(oneFace, cellsTouchingPointCount):\n")
  pfLns.append("  numFacePoints = oneFace.GetNumberOfPoints()\n")
  pfLns.append("  facePointIds = oneFace.GetPointIds()\n")
  pfLns.append("  numOutsideVertices = 0\n")
  pfLns.append("  for ii in range(0, numFacePoints):\n")
  pfLns.append("    if cellsTouchingPointCount[facePointIds.GetId(ii)] <= 2:\n")
  pfLns.append("      numOutsideVertices += 1\n")
  pfLns.append("  return (numOutsideVertices >= 2)\n")
  pfLns.append("def PhactoriIsOutsideEdgeHexCellFaceMaxAreaFace(oneCell,\n")
  pfLns.append("  indexFaceListSize, sortedFaceSizeList, cellsTouchingPointCount,\n")
  pfLns.append("  TestNthLargestFace = 2):\n")
  pfLns.append("  numFaces = oneCell.GetNumberOfFaces()\n")
  pfLns.append("  for faceIndex in range(0, numFaces):\n")
  pfLns.append("    oneFace = oneCell.GetFace(faceIndex)\n")
  pfLns.append("    if PhactoriSurfaceEdgeHexCellFaceFacesExterior(oneFace, cellsTouchingPointCount):\n")
  pfLns.append("      if indexFaceListSize[faceIndex] >= sortedFaceSizeList[-TestNthLargestFace]:\n")
  pfLns.append("        return True\n")
  pfLns.append("  return False\n")
  pfLns.append("def PhactoriIsOutsideSurfaceHexCellFaceMaxAreaFace(oneCell,\n")
  pfLns.append("  indexFaceListSize, sortedFaceSizeList, cellsTouchingPointCount,\n")
  pfLns.append("  TestNthLargestFace = 2):\n")
  pfLns.append("  numFaces = oneCell.GetNumberOfFaces()\n")
  pfLns.append("  for faceIndex in range(0, numFaces):\n")
  pfLns.append("    oneFace = oneCell.GetFace(faceIndex)\n")
  pfLns.append("    numFacePoints = oneFace.GetNumberOfPoints()\n")
  pfLns.append("    facePtNdxs = oneFace.GetPointIds()\n")
  pfLns.append("    outsideFaceIndex = -1\n")
  pfLns.append("    for jj in range(0, numFacePoints):\n")
  pfLns.append("      outsideFaceIndex = faceIndex\n")
  pfLns.append("      ptNdx = facePtNdxs.GetId(jj)\n")
  pfLns.append("      if cellsTouchingPointCount[ptNdx] != 4:\n")
  pfLns.append("        outsideFaceIndex = -1\n")
  pfLns.append("        break\n")
  pfLns.append("    if outsideFaceIndex >= 0:\n")
  pfLns.append("      break\n")
  pfLns.append("  if outsideFaceIndex < 0:\n")
  pfLns.append("    return False\n")
  pfLns.append("  if indexFaceListSize[outsideFaceIndex] >= sortedFaceSizeList[-TestNthLargestFace]:\n")
  pfLns.append("    return True\n")
  pfLns.append("  else:\n")
  pfLns.append("    return False\n")
  pfLns.append("def PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(inInputCsData,\n")
  pfLns.append("  inCellIndex, countCellsTouchingEachFace, flatnessTestRatio = 100.0,\n")
  pfLns.append("  TestNthLargestFace = 2):\n")
  pfLns.append("  oneCell = inInputCsData.GetCell(inCellIndex)\n")
  pfLns.append("  exteriorFaceList = countCellsTouchingEachFace.GetListOfExteriorFacesOnCell(oneCell)\n")
  pfLns.append("  numExteriorFaces = len(exteriorFaceList)\n")
  pfLns.append("  if numExteriorFaces == 0:\n")
  pfLns.append("    retStatus = 1\n")
  pfLns.append("  else:\n")
  pfLns.append("    indexFaceListSize, sortedFaceSizeList = PhactoriCalculateCellFaceAreasAndSort(inInputCsData, oneCell)\n")
  pfLns.append("    flatnessRatio = sortedFaceSizeList[-1] / sortedFaceSizeList[0]\n")
  pfLns.append("    if len(sortedFaceSizeList) <= 0:\n")
  pfLns.append("      retStatus = 11\n")
  pfLns.append("      return retStatus\n")
  pfLns.append("    cellIsFlat = (flatnessRatio >= flatnessTestRatio)\n")
  pfLns.append("    outsideFaceIsMaxAreaFace = False\n")
  pfLns.append("    for outsideFaceIndex in exteriorFaceList:\n")
  pfLns.append("      if indexFaceListSize[outsideFaceIndex] >= sortedFaceSizeList[-TestNthLargestFace]:\n")
  pfLns.append("        outsideFaceIsMaxAreaFace = True\n")
  pfLns.append("        break\n")
  pfLns.append("    if numExteriorFaces == 1:\n")
  pfLns.append("      if outsideFaceIsMaxAreaFace:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 11\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 8\n")
  pfLns.append("      else:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 5\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 2\n")
  pfLns.append("    elif numExteriorFaces == 2:\n")
  pfLns.append("      if outsideFaceIsMaxAreaFace:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 12\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 9\n")
  pfLns.append("      else:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 6\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 3\n")
  pfLns.append("    else:\n")
  pfLns.append("      if outsideFaceIsMaxAreaFace:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 13\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 10\n")
  pfLns.append("      else:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 7\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 4\n")
  pfLns.append("  return retStatus\n")
  pfLns.append("def PhactoriFindSurfaceStatusForOneCell(inInputCsData, inCellIndex,\n")
  pfLns.append("      cellsTouchingPointCount, flatnessTestRatio = 100.0):\n")
  pfLns.append("  oneCell = inInputCsData.GetCell(inCellIndex)\n")
  pfLns.append("  numPoints = oneCell.GetNumberOfPoints()\n")
  pfLns.append("  pointIds = oneCell.GetPointIds()\n")
  pfLns.append("  minCtpc = 9\n")
  pfLns.append("  for ii in range(0, numPoints):\n")
  pfLns.append("    testCtpc = cellsTouchingPointCount[pointIds.GetId(ii)]\n")
  pfLns.append("    if testCtpc < minCtpc:\n")
  pfLns.append("      minCtpc = testCtpc\n")
  pfLns.append("  if inCellIndex == 5056:\n")
  pfLns.append("    debugFlag = 100\n")
  pfLns.append("  else:\n")
  pfLns.append("    debugFlag = 25\n")
  pfLns.append("  if minCtpc >= 8:\n")
  pfLns.append("    retStatus = 1\n")
  pfLns.append("  else:\n")
  pfLns.append("    indexFaceListSize, sortedFaceSizeList = PhactoriCalculateCellFaceAreasAndSort(inInputCsData, oneCell)\n")
  pfLns.append("    flatnessRatio = sortedFaceSizeList[-1] / sortedFaceSizeList[0]\n")
  pfLns.append("    if len(sortedFaceSizeList) <= 0:\n")
  pfLns.append("      retStatus = 11\n")
  pfLns.append("      return retStatus\n")
  pfLns.append("    cellIsFlat = (flatnessRatio >= flatnessTestRatio)\n")
  pfLns.append("    if minCtpc >= 4:\n")
  pfLns.append("      outsideFaceIsMaxAreaFace = PhactoriIsOutsideSurfaceHexCellFaceMaxAreaFace(\n")
  pfLns.append("        oneCell, indexFaceListSize, sortedFaceSizeList,\n")
  pfLns.append("        cellsTouchingPointCount)\n")
  pfLns.append("      if outsideFaceIsMaxAreaFace:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 11\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 8\n")
  pfLns.append("      else:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 5\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 2\n")
  pfLns.append("    elif minCtpc >= 2:\n")
  pfLns.append("      outsideFaceIsMaxAreaFace = PhactoriIsOutsideEdgeHexCellFaceMaxAreaFace(\n")
  pfLns.append("        oneCell, indexFaceListSize, sortedFaceSizeList,\n")
  pfLns.append("        cellsTouchingPointCount)\n")
  pfLns.append("      if outsideFaceIsMaxAreaFace:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 12\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 9\n")
  pfLns.append("      else:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 6\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 3\n")
  pfLns.append("    elif minCtpc >= 1:\n")
  pfLns.append("      outsideFaceIsMaxAreaFace = PhactoriIsOutsideCornerHexCellFaceMaxAreaFace(\n")
  pfLns.append("        oneCell, indexFaceListSize, sortedFaceSizeList, cellsTouchingPointCount)\n")
  pfLns.append("      if outsideFaceIsMaxAreaFace:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 13\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 10\n")
  pfLns.append("      else:\n")
  pfLns.append("        if cellIsFlat:\n")
  pfLns.append("          retStatus = 7\n")
  pfLns.append("        else:\n")
  pfLns.append("          retStatus = 4\n")
  pfLns.append("    else:\n")
  pfLns.append("      retStatus = -1\n")
  pfLns.append("  return retStatus\n")
  pfLns.append("def PhactoriGetCellEdgeVector(inInputCsData, oneFaceOrCell, edgeIndex):\n")
  pfLns.append("  oneFaceEdgePointIds = oneFaceOrCell.GetEdge(edgeIndex).GetPointIds()\n")
  pfLns.append("  ptA = inInputCsData.GetPoint(oneFaceEdgePointIds.GetId(0))\n")
  pfLns.append("  ptB = inInputCsData.GetPoint(oneFaceEdgePointIds.GetId(1))\n")
  pfLns.append("  retVec = vecFromAToB(ptA,ptB)\n")
  pfLns.append("  return retVec\n")
  pfLns.append("def PhactoriFindFaceNormal(inInputCsData, oneFace):\n")
  pfLns.append("  facePointIds = oneFace.GetPointIds()\n")
  pfLns.append("  ptA = inInputCsData.GetPoint(facePointIds.GetId(0))\n")
  pfLns.append("  ptB = inInputCsData.GetPoint(facePointIds.GetId(1))\n")
  pfLns.append("  ptC = inInputCsData.GetPoint(facePointIds.GetId(2))\n")
  pfLns.append("  vecAB = vecFromAToB(ptA, ptB)\n")
  pfLns.append("  vecBC = vecFromAToB(ptB, ptC)\n")
  pfLns.append("  retNormalVec = vecCrossProduct(vecBC, vecAB)\n")
  pfLns.append("  vecNormalize2(retNormalVec, retNormalVec)\n")
  pfLns.append("  return retNormalVec\n")
  pfLns.append("def PhactoriFindLargestCellFaceNormal(inInputCsData, oneCell):\n")
  pfLns.append("  indexFaceListSize, sortedFaceSizeList = PhactoriCalculateCellFaceAreasAndSort(inInputCsData, oneCell)\n")
  pfLns.append("  biggestFaceSize = sortedFaceSizeList[-1]\n")
  pfLns.append("  for faceForNormalIndex in range(0, len(indexFaceListSize)):\n")
  pfLns.append("    if indexFaceListSize[faceForNormalIndex] == biggestFaceSize:\n")
  pfLns.append("      break\n")
  pfLns.append("  faceForNormal = oneCell.GetFace(faceForNormalIndex)\n")
  pfLns.append("  retNormalVec = PhactoriFindFaceNormal(inInputCsData, faceForNormal)\n")
  pfLns.append("  return retNormalVec\n")
  pfLns.append("def PhactoriFindSelectedAngleBetweenEdgeAndCellNormal(inInputCsData, oneCell,\n")
  pfLns.append("  compareNormal, paramUseSmallestAngle, paramOffsetIndex):\n")
  pfLns.append("  edgeCompareDotProdList = []\n")
  pfLns.append("  numEdges = oneCell.GetNumberOfEdges()\n")
  pfLns.append("  for ii in range(0, numEdges):\n")
  pfLns.append("    oneEdgeVec = PhactoriGetCellEdgeVector(inInputCsData, oneCell, ii)\n")
  pfLns.append("    vecNormalize2(oneEdgeVec, oneEdgeVec)\n")
  pfLns.append("    edgeCompareDotProd = abs(vecDotProduct(compareNormal, oneEdgeVec))\n")
  pfLns.append("    edgeCompareDotProdList.append(edgeCompareDotProd)\n")
  pfLns.append("  sortedEdgeCompareDotProdList = sorted(edgeCompareDotProdList)\n")
  pfLns.append("  if paramUseSmallestAngle:\n")
  pfLns.append("    outputIndex = max(0, numEdges - 1 - paramOffsetIndex)\n")
  pfLns.append("  else:\n")
  pfLns.append("    outputIndex = min(paramOffsetIndex, numEdges - 1)\n")
  pfLns.append("  outputAngle = 90.0 - math.degrees(math.acos(sortedEdgeCompareDotProdList[outputIndex]))\n")
  pfLns.append("  return outputAngle\n")
  pfLns.append("def PhactoriFindHeightOfCell(inInputCsData, oneCell, compareNormal):\n")
  pfLns.append("  numCellEdges = oneCell.GetNumberOfEdges()\n")
  pfLns.append("  retHeight = -1.0 \n")
  pfLns.append("  for ii in range(0, numCellEdges):\n")
  pfLns.append("    oneCellEdgeVec = PhactoriGetCellEdgeVector(inInputCsData, oneCell, ii)\n")
  pfLns.append("    dp1 = vecDotProduct(oneCellEdgeVec, compareNormal)\n")
  pfLns.append("    projectedVec = vecScale(dp1, compareNormal)\n")
  pfLns.append("    oneCellEdgeHeight = vecMagnitude(projectedVec)\n")
  pfLns.append("    if oneCellEdgeHeight > retHeight:\n")
  pfLns.append("      retHeight = oneCellEdgeHeight\n")
  pfLns.append("  return retHeight\n")
  pfLns.append("def PhactoriFindCellEdgeAngleMetricsForOneCell(inInputCsData, inCellIndex,\n")
  pfLns.append("      paramUseSmallestAngle, paramOffsetIndex):\n")
  pfLns.append("  oneCell = inInputCsData.GetCell(inCellIndex)\n")
  pfLns.append("  numEdges = oneCell.GetNumberOfEdges()\n")
  pfLns.append("  if numEdges < 6:\n")
  pfLns.append("    return 90.0, 0.0\n")
  pfLns.append("  if oneCell.GetNumberOfFaces() < 4:\n")
  pfLns.append("    return 90.0, 0.0\n")
  pfLns.append("  compareNormal = PhactoriFindLargestCellFaceNormal(inInputCsData, oneCell)\n")
  pfLns.append("  retAngle = PhactoriFindSelectedAngleBetweenEdgeAndCellNormal(inInputCsData,\n")
  pfLns.append("    oneCell, compareNormal, paramUseSmallestAngle, paramOffsetIndex)\n")
  pfLns.append("  retHeight = PhactoriFindHeightOfCell(inInputCsData, oneCell, compareNormal)\n")
  pfLns.append("  return retAngle, retHeight\n")

#phactori_combine_to_single_python_file_subpiece_end_1

