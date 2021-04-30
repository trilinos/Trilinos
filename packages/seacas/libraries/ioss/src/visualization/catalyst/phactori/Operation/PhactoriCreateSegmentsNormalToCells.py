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

class SegmentCollectionRecursionParameters():
  def __init__(self):
    self.SegmentLength = 1.0
    self.SegmentResolution = 10
    self.DirectionReferencePoint = [0.0, 0.0, 0.0]
    self.DirectionFlag = -1
    self.SegmentCollection = None
    self.workingBlockIndex = 0

  def Initialize(self, inLen, inRes, inRefPt, inDrctnFlg):
    self.SegmentLength = inLen
    self.SegmentResolution = inRes
    self.DirectionReferencePoint = inRefPt
    self.DirectionFlag = inDrctnFlg
    self.workingBlockIndex = 0
    self.SegmentCollection = []

class PcsntcSegment:
  def __init__(self, inBlockIndex, inCellIndex, inSegPtA, inSegDir):
    self.blockIndex = inBlockIndex
    self.cellIndex = inCellIndex
    self.ptA = inSegPtA
    self.segDir = inSegDir
    self.segmentLength = 1.0
    self.ptB = None

class PhactoriCreateSegmentsNormalToCells(PhactoriOperationSpecifics):
  """takes a grid with cells and for each cell creates a line segment
     which starts in the middle of the cell and goes 'normal' to the cell.
     Normal is defined by finding the largest face on the cell and using
     that face to define the cell normal"""

  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.SegmentLength = 1.0
    self.SegmentResolution = 10
    self.DirectionReferencePoint = [0.0, 0.0, 0.0]
    self.DirectionFlag = 1
    self.RecursionResults = None
    self.IncomingPvFilter = None
    self.OutgoingGroupFilter = None
    self.myVtkPoints = None
    self.myVtkCellArray = None
    self.myVtkPolyData = None
    self.addPerPointCoordinateSystem = True
    self.localCrdSysXName = "localcrdsysx"
    self.localCrdSysYName = "localcrdsysy"
    self.localCrdSysZName = "localcrdsysz"
    self.addPerPointSegmentId = True
    self.segmentIdArrayName = "segmentid"
    self.localAxisXReference = [1.0, 0.0, 0.0]
    self.localAxisXReference2 = [0.0, 1.0, 0.0]
    self.geometricSamplingMultiplier = 1.0
    self.segmentSampleMethod = "linear"

  def ParseCheckFor3ValueList(self, itemToCheck, itemKey):
    if (type(itemToCheck) is not list) or (len(itemToCheck) != 3):
      myDebugPrint3AndException(
        "PhactoriCreateSegmentsNormalToCells.ParseParametersFromJson:\n"
        "'" + itemKey + "' must be 3 value list\n")

  def ParseParametersFromJson(self, inJson):
    key1 = "segment length"
    if key1 in inJson:
      self.SegmentLength = inJson[key1]
      if self.SegmentLength <= 0.0:
        myDebugPrint3AndException(
          "PhactoriCreateSegmentsNormalToCells.ParseParametersFromJson:\n"
          "'" + key1 + "' must be > 0.0\n")
    key1 = "segment resolution"
    if key1 in inJson:
      self.SegmentResolution = inJson[key1]
      if self.SegmentResolution < 2:
        myDebugPrint3AndException(
          "PhactoriCreateSegmentsNormalToCells.ParseParametersFromJson:\n"
          "'" + key1 + "' must be >= 2\n")
    key1 = "segment direction reference point"
    if key1 in inJson:
      self.DirectionReferencePoint = inJson[key1]
      self.ParseCheckFor3ValueList(self.DirectionReferencePoint, key1)
    key1 = "segment direction flag"
    if key1 in inJson:
      self.DirectionFlag = inJson[key1]
    key1 = "geometric progression sampling"
    if key1 in inJson:
      self.segmentSampleMethod = "geometric"
      self.geometricSamplingMultiplier = inJson[key1]
      if self.geometricSamplingMultiplier <= 0.0:
        myDebugPrint3AndException(
          "PhactoriCreateSegmentsNormalToCells.ParseParametersFromJson:\n"
          "'" + key1 + "' must be > 0.0\n")
    key1 = "global x axis"
    if key1 in inJson:
      self.localAxisXReference = inJson[key1]
      self.ParseCheckFor3ValueList(self.localAxisXReference, key1)
      self.localAxisXReference = vecNormalize(self.localAxisXReference)
      key2 = "global x axis backup"
      if key2 in inJson:
        self.localAxisXReference2 = inJson[key2]
        self.ParseCheckFor3ValueList(self.localAxisXReference2, key2)
        self.localAxisXReference2 = vecNormalize(self.localAxisXReference2)
      else:
        #construct second axis for when segment is parallel to this one
        dpxa = abs(vecDotProduct(self.localAxisXReference, [1.0, 0.0, 0.0]))
        dpya = abs(vecDotProduct(self.localAxisXReference, [0.0, 1.0, 0.0]))
        dpza = abs(vecDotProduct(self.localAxisXReference, [0.0, 0.0, 1.0]))
        if (dpxa <= dpya) and (dpxa <= dpza):
          self.localAxisXReference2 = [1.0, 0.0, 0.0]
        elif (dpya <= dpxa) and (dpya <= dpza):
          self.localAxisXReference2 = [0.0, 1.0, 0.0]
        else:
          self.localAxisXReference2 = [0.0, 0.0, 1.0]
    key1 = "add per point coordinate system"
    if key1 in inJson:
      self.addPerPointCoordinateSystem = inJson[key1]
    key1 = "add per point segment id"
    if key1 in inJson:
      self.addPerPointSegmentId = inJson[key1]

  @staticmethod
  def CreateSegmentsNormalToCellsInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("FindCellEdgeAngleMetricsInBlock entered\n")

    inParameters.workingBlockIndex += 1

    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    if PhactoriDbg(100):
      myDebugPrint3("numCells: " + str(numCells) + \
        "  numPoints: " + str(numPoints) + "\n")

    for cellIndex in range(0, numCells):
      oneCell = inInputCsData.GetCell(cellIndex)
      cellFaceNormal = PhactoriFindLargestCellFaceNormal(inInputCsData, oneCell)
      cellPtAvg =  PhactoriGetAverageOfCellPoints(inInputCsData, oneCell)
      if PhactoriDbg(100):
        if cellIndex % 1000 == 0:
          myDebugPrint3("cellIndex: " + str(cellIndex) + " c: " + \
            str(cellPtAvg) + " n: " + str(cellFaceNormal) + "\n")
      newSegmentItem = PcsntcSegment(inParameters.workingBlockIndex, cellIndex,
        cellPtAvg, cellFaceNormal)
      inParameters.SegmentCollection.append(newSegmentItem)

    if PhactoriDbg(100):
      myDebugPrint3("FindCellEdgeAngleMetricsInBlock returning\n")

  def FindSecondPointForSegments(self):
    for oneLineInfo in self.RecursionResults.SegmentCollection:
      vecFromCellToRef = vecNormalize(vecFromAToB(oneLineInfo.ptA, self.DirectionReferencePoint))
      testDotProd = vecDotProduct(oneLineInfo.segDir, vecFromCellToRef)
      if testDotProd < 0.0:
        if self.DirectionFlag < 0:
          vecScaleInPlace(-1.0, oneLineInfo.segDir)
      else:
        if self.DirectionFlag >= 0:
          vecScaleInPlace(-1.0, oneLineInfo.segDir)
      oneSegLen = self.SegmentLength
      oneLineInfo.segmentLength = oneSegLen
      oneLineInfo.ptB = vecMultiplyAdd(oneLineInfo.ptA, oneLineInfo.segDir, oneSegLen)

  def AddPerPointSegmentIdIfRequested(self):
    if self.addPerPointSegmentId == False:
      return
    #find most segments out of any processes and use that to construct unique
    #for each segment
    localNumSegmentsArray = [len(self.RecursionResults.SegmentCollection)]
    maxNumSegmentsArray = UseReduceOnIntegerList(localNumSegmentsArray, 0)
    maxNumSegmentsOnAnyProcess = maxNumSegmentsArray[0]
    totalNumPts = len(self.RecursionResults.SegmentCollection)*self.SegmentResolution
    segIdArray = vtk.vtkIntArray()
    segIdArray.SetNumberOfComponents(3)
    segIdArray.SetNumberOfTuples(totalNumPts)
    segIdArray.SetName(self.segmentIdArrayName)
    localPid = SmartGetLocalProcessId()
    numPid = SmartGetNumberOfProcesses()
    ptNdx = 0
    segNdx = 0
    #using +1 guarantees sequence gap between each process so which segments
    #went with which proceses can be reconstructed if necessary
    localSegIdCntr = localPid * (maxNumSegmentsOnAnyProcess+1)
    idTuple = [localSegIdCntr, segNdx, 0]
    for oneLineInfo in self.RecursionResults.SegmentCollection:
      thisSegmentPointIndex = 0
      for resNdx in range(0,self.SegmentResolution):
        idTuple[2] = thisSegmentPointIndex
        segIdArray.SetTuple(ptNdx, idTuple)
        ptNdx += 1
        thisSegmentPointIndex += 1
      localSegIdCntr += 1
      segNdx += 1
      idTuple[0] = localSegIdCntr
      idTuple[1] = segNdx
    self.myVtkPolyData.GetPointData().AddArray(segIdArray)

  def AddPerPointCoordinateSystemIfRequested(self):
    if self.addPerPointCoordinateSystem == False:
      return
    ptNdx = 0
    totalNumPts = len(self.RecursionResults.SegmentCollection)*self.SegmentResolution
    crdSysXArray = vtk.vtkDoubleArray()
    crdSysXArray.SetNumberOfComponents(3)
    crdSysXArray.SetNumberOfTuples(totalNumPts)
    crdSysXArray.SetName(self.localCrdSysXName)
    crdSysYArray = vtk.vtkDoubleArray()
    crdSysYArray.SetNumberOfComponents(3)
    crdSysYArray.SetNumberOfTuples(totalNumPts)
    crdSysYArray.SetName(self.localCrdSysYName)
    crdSysZArray = vtk.vtkDoubleArray()
    crdSysZArray.SetNumberOfComponents(3)
    crdSysZArray.SetNumberOfTuples(totalNumPts)
    crdSysZArray.SetName(self.localCrdSysZName)
    for oneLineInfo in self.RecursionResults.SegmentCollection:
      localYAxis = vecNormalize(oneLineInfo.segDir)
      localZAxis = vecCrossProduct(self.localAxisXReference, localYAxis)
      zAxisLen = vecMagnitude(localZAxis)
      if zAxisLen < 1.0e-20:
        localZAxis = vecCrossProduct(self.localAxisXReference2, localYAxis)
      vecNormalize2(localZAxis, localZAxis)
      localXAxis = vecCrossProduct(localYAxis, localZAxis)
      vecNormalize2(localXAxis, localXAxis)
      for resNdx in range(0,self.SegmentResolution):
        crdSysXArray.SetTuple(ptNdx, localXAxis)
        crdSysYArray.SetTuple(ptNdx, localYAxis)
        crdSysZArray.SetTuple(ptNdx, localZAxis)
        ptNdx += 1
    self.myVtkPolyData.GetPointData().AddArray(crdSysXArray)
    self.myVtkPolyData.GetPointData().AddArray(crdSysYArray)
    self.myVtkPolyData.GetPointData().AddArray(crdSysZArray)

  def CalculateGeometricProgressionSampleValues(self, oneLineInfo):
    samplePos = [0.0]
    samplePt = 0.0
    sampleSpacing = self.geometricSamplingMultiplier
    for ii in range(0, self.SegmentResolution-1):
      samplePt += sampleSpacing
      samplePos.append(samplePt)
      sampleSpacing *= self.geometricSamplingMultiplier
    maxSamplePt = samplePos[-1]
    adj1 = oneLineInfo.segmentLength/maxSamplePt
    for ii in range(0, self.SegmentResolution-1):
      samplePos[ii] *= adj1
    samplePos[-1] = oneLineInfo.segmentLength
    return samplePos

  def CreateSegmentWithChosenSamplingMethod(self, oneLineInfo, ptNdx):
    myNewPolyLine = vtk.vtkPolyLine()
    myPlyLnIds = myNewPolyLine.GetPointIds()
    myPlyLnIds.SetNumberOfIds(self.SegmentResolution)
    if self.segmentSampleMethod == "linear" or self.segmentSampleMethod == "geometric":
      if self.segmentSampleMethod == "linear":
        self.geometricSamplingMultiplier = 1.0
      samplePos = self.CalculateGeometricProgressionSampleValues(oneLineInfo)
      for resNdx in range(0,self.SegmentResolution):
        ratio1 = samplePos[resNdx]
        newPtAlongLine = vecMultiplyAdd(oneLineInfo.ptA, oneLineInfo.segDir, ratio1)
        self.myVtkPoints.SetPoint(ptNdx, newPtAlongLine)
        myPlyLnIds.SetId(resNdx, ptNdx)
        ptNdx += 1
    else:
      myDebugPrint3AndException(
        "PhactoriCreateSegmentsNormalToCells.CreateSegmentWithChosenSamplingMethod:\n"
        "bad self.segmentSampleMethod: " + str(self.segmentSampleMethod) + "\n")
    return myNewPolyLine

  def CreateParaViewSourcesForSegments(self):
    pvLineList = []
    if PhactoriDbg(100):
      myDebugPrint3("CreateParaViewSourcesForSegments\n" + \
        "self.SegmentLength: " + str(self.SegmentLength) + "\n" + \
        "self.SegmentResolution: " + str(self.SegmentResolution) + "\n" + \
        "self.DirectionReferencePoint: " + str(self.DirectionReferencePoint) + "\n" + \
        "self.DirectionFlag: " + str(self.DirectionFlag) + "\n")

    self.myVtkPoints = vtk.vtkPoints()
    self.myVtkCellArray = vtk.vtkCellArray()
    self.myVtkPolyData = vtk.vtkPolyData()

    self.FindSecondPointForSegments()
    ptNdx = 0
    totalNumPts = len(self.RecursionResults.SegmentCollection)*self.SegmentResolution
    self.myVtkPoints.SetNumberOfPoints(totalNumPts)
    for oneLineInfo in self.RecursionResults.SegmentCollection:
      myNewPolyLine = self.CreateSegmentWithChosenSamplingMethod(oneLineInfo, ptNdx)
      ptNdx += self.SegmentResolution
      self.myVtkCellArray.InsertNextCell(myNewPolyLine)

    self.myVtkPolyData.SetPoints(self.myVtkPoints)
    #self.myVtkPolyData.SetVerts(self.myVtkCellArray)
    self.myVtkPolyData.SetLines(self.myVtkCellArray)

    self.AddPerPointCoordinateSystemIfRequested()
    self.AddPerPointSegmentIdIfRequested()

    newParaViewFilter = PVTrivialProducer()
    newParaViewFilter.GetClientSideObject().SetOutput(self.myVtkPolyData)
    return newParaViewFilter

  def CreateSegmentsForAllCells(self, inInputFilter):
    recursionControlObj = PhactoriParaviewMultiBlockRecursionControl()
    self.RecursionResults = SegmentCollectionRecursionParameters()
    self.RecursionResults.Initialize(
      self.SegmentLength,
      self.SegmentResolution,
      self.DirectionReferencePoint,
      self.DirectionFlag)
    recursionControlObj.mParameters = self.RecursionResults
    recursionControlObj.mOperationToDoPerBlock = self.CreateSegmentsNormalToCellsInBlock
    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionControlObj, inInputFilter)

  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriCreateSegmentsNormalToCells.CreateParaViewFilter "
          "entered\n", 100)

    self.IncomingPvFilter = inInputFilter
    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    #go through all cells and create a segment for each one.
    #(find the geometry for all segments)
    self.CreateSegmentsForAllCells(inInputFilter)

    self.OutgoingGroupFilter = self.CreateParaViewSourcesForSegments()
    UpdatePipelineWithCurrentTimeArgument(self.OutgoingGroupFilter)

    SetActiveSource(self.OutgoingGroupFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriCreateSegmentsNormalToCells.CreateParaViewFilter "
          "returning\n", 100)

    return self.OutgoingGroupFilter

#phactori_combine_to_single_python_file_subpiece_end_1
