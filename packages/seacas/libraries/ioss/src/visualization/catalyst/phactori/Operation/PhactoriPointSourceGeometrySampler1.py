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
from .QuadInfoGQC import *
from .PhactoriMpiUtilities import *

#phactori_combine_to_single_python_file_subpiece_begin_1

class PhactoriPointSourceGeometrySampler1(PhactoriOperationSpecifics):
  """Filter/operation which reads in a .json file which is a list of 3d
     coordinates and creates a new point source from that list. Resulting
     point source will have 1 element an N points. This source is intented
     to work correctly in parallel for Catalyst or pvbatch symmetric mode.
     The json file which is read in will only be read on one process and
     mpi broadcast is used to distribute the list, rather than having each
     process read the json file."""
  def __init__(self):
    self.JsonListFileName = "PhactoriPointSourceGeometrySampler1.json"
    self.JsonList = None
    self.ParaviewPointSource = None
    self.SourceVtkPolyData = None
    self.myVtkPoints = vtk.vtkPoints()
    self.myVtkPolyVertex = vtk.vtkPolyVertex()
    self.myVtkCellArray = vtk.vtkCellArray()
    self.myVtkPolyData = vtk.vtkPolyData()
    self.NearestQuadList = None
    self.DistanceList = None

  def ValidateJsonPointList(self):
    numPoints = len(self.JsonList)
    if numPoints < 1:
      myDebugPrint3AndException(
          "PhactoriPointSourceGeometrySampler1::ValidateJsonPointList\n"
          "list must have at least one element\n")

    for ptNdx in range(0,numPoints):
      jsonPt = self.JsonList[ptNdx]
      if len(jsonPt) != 3:
        errStr = "PhactoriPointSourceGeometrySampler1::ValidateJsonPointList\n" \
          "point with index " + str(ptNdx) + "does not have three elements\n"
        myDebugPrint3AndException(errStr)


  def ParseParametersFromJson(self, inJson):
    if 'filename' in inJson:
      self.JsonListFileName = inJson['filename']
    else:
      myDebugPrint3AndException(
          "PhactoriPointSourceGeometrySampler1::ParseParametersFromJson\n"
          "Error:  must have 'filename' key\n")

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

  def UseMpiToGetGlobalQuadsClosest(self, inLocalQuadList, inLocalDistSqrdList):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSourceGeometrySampler1.UseMpiToGetGlobalQuadsClosest entered\n", 100)

    pidWithDataList, globalDistSqrdList = self.GetPidWithLeastValueList(inLocalDistSqrdList)
    if PhactoriDbg(100):
      myDebugPrint3("pidWithDataList:\n" + str(pidWithDataList) + "\nglobalDistSqrdList:\n" + str(globalDistSqrdList) + "\n")

    #convert quad list to array of doubles and ints, use mpi reduce to share
    #the values, then convert back to quad list

    #convert quad list to array of doubles and ints
    quadInfoDoubleArray = []
    quadInfoIntArray = []
    myPid = SmartGetLocalProcessId()
    for ii, oneQuad in enumerate(inLocalQuadList):
      if pidWithDataList[ii] == myPid:
        oneQuad.AddToInfoArrays(quadInfoDoubleArray, quadInfoIntArray)
      else:
        oneQuad.AddZeroQuadToInfoArrays(quadInfoDoubleArray, quadInfoIntArray)

    #use mpi reduce to spread array correctly
    globalQuadInfoDoubleArray = UseReduceOnFloatList(quadInfoDoubleArray, 2)
    globalQuadInfoIntArray = UseReduceOnIntegerList(quadInfoIntArray, 2)

    #now create return global quad list from arrays
    numQuads = len(inLocalQuadList)
    returnGlobalQuadList = []
    for ii in range(0,numQuads):
      newQuad = QuadInfoGQC()
      newQuad.PopulateFromInfoArrays(ii, globalQuadInfoDoubleArray, globalQuadInfoIntArray)
      returnGlobalQuadList.append(newQuad)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSourceGeometrySampler1.UseMpiToGetGlobalQuadsClosest returning\n", 100)
    return returnGlobalQuadList, globalDistSqrdList

  def CreateSamplePoints(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSourceGeometrySampler1.CreateSamplePoints entered\n", 100)

    thisProcessNearestQuadList, thisProcDistSqrdList = self.mPhactoriOperationBlockOwner.\
      GetQuadsClosestToPointsOnThisProcessFromParaViewFilter(inInputFilter, self.JsonList)

    if PhactoriDbg(100):
      myDebugPrint3("quads before mpi\n", 100)
      for ii, oneQuad in enumerate(thisProcessNearestQuadList):
        myDebugPrint3(str(ii) + ":  " + str(thisProcDistSqrdList[ii]) + "\n" + oneQuad.ToString(), 100)

    self.NearestQuadList, self.DistanceList = self.UseMpiToGetGlobalQuadsClosest(
      thisProcessNearestQuadList, thisProcDistSqrdList)

    if PhactoriDbg(100):
      myDebugPrint3("quads after mpi\n", 100)
      for ii, oneQuad in enumerate(self.NearestQuadList):
        myDebugPrint3(str(ii) + ":  " + str(self.NearestQuadList[ii]) + "\n" + oneQuad.ToString(), 100)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSourceGeometrySampler1.CreateSamplePoints returning\n", 100)

  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSourceGeometrySampler1.CreateParaViewFilter entered\n", 100)

    savedActiveSource = GetActiveSource()

    self.JsonList = ReadAndMpiBroadcastJsonFile(self.JsonListFileName)

    self.ValidateJsonPointList()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    self.CreateSamplePoints(inInputFilter)

    allSamplePoints = []
    for ii, oneQuad in enumerate(self.NearestQuadList):
      #oneQuadSamplePoints = oneQuad.MakePerpendicularSamplePoints()
      oneQuadSamplePoints = oneQuad.MakeLogSamplePointsAlongRayFromNearbyPoint(self.JsonList[ii], 0.0000001, 1.414, 0.05)
      allSamplePoints.extend(oneQuadSamplePoints)

    numPoints = len(allSamplePoints)
    self.myVtkPoints.SetNumberOfPoints(numPoints)
    self.myVtkPolyVertex.GetPointIds().SetNumberOfIds(numPoints)

    vtkPolyVertPtIds = self.myVtkPolyVertex.GetPointIds()
    usePoint = [0.0, 0.0, 0.0]
    for ptNdx in range(0,numPoints):
      #self.myVtkPoints.SetPoint(ptNdx, self.JsonList[ptNdx])
      quadPtAvg = allSamplePoints[ptNdx]
      usePoint[0] = quadPtAvg[0]
      usePoint[1] = quadPtAvg[1]
      usePoint[2] = quadPtAvg[2]
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
      myDebugPrint3("PhactoriPointSourceGeometrySampler1.CreateParaViewFilter returning\n", 100)
    return self.ParaviewPointSource

#phactori_combine_to_single_python_file_subpiece_end_1
