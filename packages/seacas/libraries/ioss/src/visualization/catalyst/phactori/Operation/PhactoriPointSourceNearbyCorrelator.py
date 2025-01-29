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

from .PhactoriVectorLibrary import *
from phactori import *
from paraview.simple import *
from .PhactoriMpiUtilities import *
import sys

#phactori_combine_to_single_python_file_subpiece_begin_1

class PhactoriPointSourceNearbyCorrelator(PhactoriOperationSpecifics):
  """Filter/operation which reads in a .json file which is a list of 3d
     coordinates and creates a new point source from that list. Resulting
     point source will have 1 element an N points. This source is intented
     to work correctly in parallel for Catalyst or pvbatch symmetric mode.
     The json file which is read in will only be read on one process and
     mpi broadcast is used to distribute the list, rather than having each
     process read the json file."""
  def __init__(self):
    self.JsonListFileName = "PhactoriPointSourceFromJsonList.json"
    self.JsonList = None
    self.NearestDatapointList = None
    self.DisplacedPointList = None
    self.DisplacementDistance = 0.0
    self.ParaviewPointSource = None
    self.SourceVtkPolyData = None
    self.myVtkPoints = vtk.vtkPoints()
    self.myVtkPolyVertex = vtk.vtkPolyVertex()
    self.myVtkCellArray = vtk.vtkCellArray()
    self.myVtkPolyData = vtk.vtkPolyData()

  def ValidateJsonPointList(self):
    numPoints = len(self.JsonList)
    if numPoints < 1:
      myDebugPrint3AndException(
          "PhactoriPointSourceFromJsonList::ValidateJsonPointList\n"
          "list must have at least one element\n")

    for ptNdx in range(0,numPoints):
      jsonPt = self.JsonList[ptNdx]
      if len(jsonPt) != 3:
        errStr = "PhactoriPointSourceFromJsonList::ValidateJsonPointList\n" \
          "point with index " + str(ptNdx) + "does not have three elements\n"
        myDebugPrint3AndException(errStr)


  def ParseParametersFromJson(self, inJson):
    if 'filename' in inJson:
      self.JsonListFileName = inJson['filename']
    else:
      myDebugPrint3AndException(
          "PhactoriPointSourceFromJsonList::ParseParametersFromJson\n"
          "Error:  must have 'filename' key\n")

    self.DisplacementDistance = getParameterFromBlock(inJson, "absolute displacement", self.DisplacementDistance)

  def CalculateNearestPointDistanceSquaredList(self, inSourcePointList, inTargetPointList):
    if PhactoriDbg(100):
      myDebugPrint3("CalculateNearestPointDistanceSquaredList entered\n")
    resultDistancedSquaredList = []
    resultClosestTargetIndexList = []
    localmaxfloat = sys.float_info.max
    if len(inTargetPointList) == 0:
      if PhactoriDbg(100):
        myDebugPrint3("no points on this process, returning noresult list\n")
      for idx in range(0, len(inSourcePointList)):
        resultDistancedSquaredList.append(localmaxfloat)
      return resultDistancedSquaredList

    for srcPnt in inSourcePointList:
      if PhactoriDbg(100):
        myDebugPrint3("doing point: " + str(srcPnt) + "\n")
      tgtPnt = inTargetPointList[0]
      closetsDistSqrd = vecDistanceSquared(srcPnt, tgtPnt)
      closestIndex = 0
      for idx, tgtPnt in enumerate(inTargetPointList):
        testDist = vecDistanceSquared(srcPnt, tgtPnt)
        if testDist < closetsDistSqrd:
          closetsDistSqrd = testDist
          closestIndex = idx
      if PhactoriDbg(100):
        myDebugPrint3("result: " + str(closetsDistSqrd) + "  " + str(closestIndex) + "\n")
      resultDistancedSquaredList.append(closetsDistSqrd)
      resultClosestTargetIndexList.append(closestIndex)
    if PhactoriDbg(100):
      myDebugPrint3("CalculateNearestPointDistanceSquaredList returning\n")
    return resultDistancedSquaredList, resultClosestTargetIndexList

  def UseMpiToFindGlobalNearestPointList(self, inTargetPointList, inDistanceSquaredList, inTargetPointIndexList):
    if PhactoriDbg(100):
      myDebugPrint3("UseMpiToFindGlobalNearestPointList entered\n")
    globalMinDistanceSquaredList = UseReduceOnFloatList(inDistanceSquaredList, 1)
    myPid = SmartGetLocalProcessId()
    numSrcPnts = len(inDistanceSquaredList)
    localPidWithClosestPoint = []
    for idx in range(0,numSrcPnts):
      if inDistanceSquaredList[idx] == globalMinDistanceSquaredList[idx]:
        if PhactoriDbg(100):
          myDebugPrint3("I match closest " + str(idx) + "  " + \
            str(inDistanceSquaredList[idx]) + "  " + str(myPid) + "\n")
        localPidWithClosestPoint.append(myPid)
      else:
        localPidWithClosestPoint.append(-1)

    globalPidWithClosestList = UseReduceOnIntegerList(localPidWithClosestPoint, 0)

    localXyzList = []
    localNodeIdList = []
    for idx in range(0,numSrcPnts):
      if globalPidWithClosestList[idx] == myPid:
        nearestPt = inTargetPointList[inTargetPointIndexList[idx]]
        if PhactoriDbg(100):
          myDebugPrint3("I am closest with highest pid " + str(idx) + "  " + \
            str(inTargetPointIndexList[idx]) + "  " + str(myPid) + "\n" + str(nearestPt) + "\n")
        localXyzList.append(nearestPt[0])
        localXyzList.append(nearestPt[1])
        localXyzList.append(nearestPt[2])
        localNodeIdList.append(int(nearestPt[3]))
      else:
        localXyzList.append(0.0)
        localXyzList.append(0.0)
        localXyzList.append(0.0)
        localNodeIdList.append(0)

    if PhactoriDbg(100):
      myDebugPrint3("localXyzList:\n" + str(localXyzList) + "\nlocalNodIdList:\n" + str(localNodeIdList) + "\n")
    globalXyz = UseReduceOnFloatList(localXyzList, 2)
    globalNodeId = UseReduceOnIntegerList(localNodeIdList, 2)
    if PhactoriDbg(100):
      myDebugPrint3("globalXyz:\n" + str(globalXyz) + "\nglobalNodeId:\n" + str(globalNodeId) + "\n")

    resultPoints = []
    for idx in range(0, numSrcPnts):
      xyzIdx = idx*3
      resultPoints.append([globalXyz[xyzIdx], globalXyz[xyzIdx+1], globalXyz[xyzIdx+2], globalNodeId[idx]])

    if PhactoriDbg(100):
      myDebugPrint3("resultPoints:\n" + str(resultPoints) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("UseMpiToFindGlobalNearestPointList returning\n")
    return resultPoints

  def CalculateOneDisplacedPoint(self, srcPt, tgtPt, dispDist):
      deltaVec = vecFromAToB(srcPt, tgtPt)
      moveVec = vecNormalizeWithSmallCheck(deltaVec, 1e-10, [0.0, 1.0, 0.0])
      displacedPoint = vecMultiplyAdd(tgtPt, moveVec, dispDist)
      return displacedPoint

  def CalculateDisplacedPointList(self):
    """take the nearest data point list, move each a displacment amount along
       the vector from the source point to that point to get a set of
       displaced points"""
    if PhactoriDbg(100):
      myDebugPrint3("CalculateDisplacedPointList: " + str(self.DisplacementDistance) + "\n")

    if self.DisplacementDistance == 0.0:
      self.DisplacedPointList = self.NearestDatapointList
      return

    self.DisplacedPointList = []
    numPts = len(self.JsonList)
    for ptNdx in range(0,numPts):
      srcPt = self.JsonList[ptNdx]
      tgtPt = self.NearestDatapointList[ptNdx]
      displacedPt = self.CalculateOneDisplacedPoint(srcPt, tgtPt, self.DisplacementDistance)
      self.DisplacedPointList.append(displacedPt)

  def FindCorrelatedInputOperationPoints(self, inInputFilter):
    """in parallel, using mpi to help, find the nodes from the input filter
       which are closest to the points in self.mJsonPointList, and record the
       3d location of those points"""

    #grab local process points
    thisProcessPointList = self.mPhactoriOperationBlockOwner.MakeListOfAllPointsAndNodeIdsOnThisProcessFromParaViewFilter(inInputFilter)

    #calc local nearest points
    nearestPointDistanceSquaredList, targetPointIndexList = self.CalculateNearestPointDistanceSquaredList(self.JsonList, thisProcessPointList)

    #use mpi to figure out global closest points
    self.NearestDatapointList = self.UseMpiToFindGlobalNearestPointList(thisProcessPointList, nearestPointDistanceSquaredList,
      targetPointIndexList)

    #get set of points (presumably slightly) displaced from self.NearestDatapointList
    self.CalculateDisplacedPointList()

    if PhactoriDbg(100):
      myDebugPrint3("FindCorrelatedInputOperationPoints result:\nNearest points:\n")
      for ndx, sourcePt in enumerate(self.JsonList):
        nearestTargetPoint = self.NearestDatapointList[ndx]
        distSqrd = vecDistanceSquared(sourcePt, nearestTargetPoint)
        distx = math.sqrt(distSqrd)
        myDebugPrint3(str(ndx) + ": " + str(sourcePt) + "  " + str(distSqrd) + "  " + str(distx) + "  " + str(self.NearestDatapointList[ndx]) + "\n")
      myDebugPrint3("Displaced points:\n")
      for ndx, sourcePt in enumerate(self.JsonList):
        matchingDisplacedPoint = self.DisplacedPointList[ndx]
        distSqrd = vecDistanceSquared(sourcePt, matchingDisplacedPoint)
        distx = math.sqrt(distSqrd)
        myDebugPrint3(str(ndx) + ": " + str(sourcePt) + "  " + str(distSqrd) + "  " + str(distx) + "  " + str(self.DisplacedPointList[ndx]) + "\n")


  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSourceFromJsonList.CreateParaViewFilter entered\n", 100)

    savedActiveSource = GetActiveSource()

    self.JsonList = ReadAndMpiBroadcastJsonFile(self.JsonListFileName)

    self.ValidateJsonPointList()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    self.FindCorrelatedInputOperationPoints(inInputFilter)

    numPoints = len(self.DisplacedPointList)
    self.myVtkPoints.SetNumberOfPoints(numPoints)
    self.myVtkPolyVertex.GetPointIds().SetNumberOfIds(numPoints)

    vtkPolyVertPtIds = self.myVtkPolyVertex.GetPointIds()
    usePoint = [0.0, 0.0, 0.0]
    for ptNdx in range(0,numPoints):
      #self.myVtkPoints.SetPoint(ptNdx, self.DisplacedPointList[ptNdx])
      jsonPt = self.DisplacedPointList[ptNdx]
      usePoint[0] = jsonPt[0]
      usePoint[1] = jsonPt[1]
      usePoint[2] = jsonPt[2]
      self.myVtkPoints.SetPoint(ptNdx, usePoint)
      vtkPolyVertPtIds.SetId(ptNdx, ptNdx)

    self.myVtkPolyData.SetPoints(self.myVtkPoints)
    self.myVtkCellArray.InsertNextCell(self.myVtkPolyVertex)
    self.myVtkPolyData.SetVerts(self.myVtkCellArray)

    self.ParaviewPointSource = PVTrivialProducer()
    self.ParaviewPointSource.GetClientSideObject().SetOutput(self.myVtkPolyData)

    SetActiveSource(self.ParaviewPointSource)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(self.ParaviewPointSource)

    if PhactoriDbg(100):
      pssfjlcso = self.ParaviewPointSource.GetClientSideObject()
      pssfjlodo = pssfjlcso.GetOutputDataObject(0)
      myDebugPrint3("pssfjlodo: " + str(pssfjlodo.GetClassName()) + " numcells " + str(pssfjlodo.GetNumberOfCells()) + " numpoints " + str(pssfjlodo.GetNumberOfPoints()) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSourceFromJsonList.CreateParaViewFilter returning\n", 100)
    return self.ParaviewPointSource

#phactori_combine_to_single_python_file_subpiece_end_1

