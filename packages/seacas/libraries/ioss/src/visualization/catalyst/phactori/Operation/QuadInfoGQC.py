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

import vtk
from .PhactoriVectorLibrary import *
from phactori import myDebugPrint3
from phactori import PhactoriDbg

#phactori_combine_to_single_python_file_subpiece_begin_1
class QuadInfoGQC:
  def __init__(self):
    self.ptA = [0.0, 0.0, 0.0]
    self.ptB = [0.0, 0.0, 0.0]
    self.ptC = [0.0, 0.0, 0.0]
    self.ptD = [0.0, 0.0, 0.0]
    #using average of points instead of center
    self.ptAverage = [0.0, 0.0, 0.0]
    self.normal = None
    self.cellId = -1
    self.index = -1

  def PopulateFromVtkData(self, inInputCsData, globalCellIdArray, quadIndex):
    self.index = quadIndex
    cellPointIds = vtk.vtkIdList()
    inInputCsData.GetCellPoints(quadIndex, cellPointIds)
    numids = cellPointIds.GetNumberOfIds()
    if numids != 4:
      myDebugPrint3AndException("GetQuadClosestToPointsInBlock2Params: not a quad")
    inInputCsData.GetPoint(cellPointIds.GetId(0), self.ptA)
    inInputCsData.GetPoint(cellPointIds.GetId(1), self.ptB)
    inInputCsData.GetPoint(cellPointIds.GetId(2), self.ptC)
    inInputCsData.GetPoint(cellPointIds.GetId(3), self.ptD)
    self.ptAverage[0] = (self.ptA[0] + self.ptB[0] + self.ptC[0] + self.ptD[0]) * 0.25
    self.ptAverage[1] = (self.ptA[1] + self.ptB[1] + self.ptC[1] + self.ptD[1]) * 0.25
    self.ptAverage[2] = (self.ptA[2] + self.ptB[2] + self.ptC[2] + self.ptD[2]) * 0.25
    if globalCellIdArray != None:
      self.cellId = globalCellIdArray.GetValue(quadIndex)

  def PopulateFromInfoArrays(self, index, doubleArray, integerArray):
    """given arrays which were populated with the help of AddToInfoArrays(),
       set up this quad to match the one represented in the arrays at the
       given index"""
    ptIndx = index*15
    idIndx = index
    self.ptA[0] = doubleArray[ptIndx]
    self.ptA[1] = doubleArray[ptIndx+1]
    self.ptA[2] = doubleArray[ptIndx+2]
    self.ptB[0] = doubleArray[ptIndx+3]
    self.ptB[1] = doubleArray[ptIndx+4]
    self.ptB[2] = doubleArray[ptIndx+5]
    self.ptC[0] = doubleArray[ptIndx+6]
    self.ptC[1] = doubleArray[ptIndx+7]
    self.ptC[2] = doubleArray[ptIndx+8]
    self.ptD[0] = doubleArray[ptIndx+9]
    self.ptD[1] = doubleArray[ptIndx+10]
    self.ptD[2] = doubleArray[ptIndx+11]
    self.ptAverage[0] = doubleArray[ptIndx+12]
    self.ptAverage[1] = doubleArray[ptIndx+13]
    self.ptAverage[2] = doubleArray[ptIndx+14]
    self.cellId = integerArray.append(self.cellId)

  def AddToInfoArrays(self, outDoubleArray, outIntegerArray):
    """add the quad points, point average, and cellId to the argument arrays
       for purspose of doing mpi broadcast/reduce work"""
    outDoubleArray.append(self.ptA[0])
    outDoubleArray.append(self.ptA[1])
    outDoubleArray.append(self.ptA[2])
    outDoubleArray.append(self.ptB[0])
    outDoubleArray.append(self.ptB[1])
    outDoubleArray.append(self.ptB[2])
    outDoubleArray.append(self.ptC[0])
    outDoubleArray.append(self.ptC[1])
    outDoubleArray.append(self.ptC[2])
    outDoubleArray.append(self.ptD[0])
    outDoubleArray.append(self.ptD[1])
    outDoubleArray.append(self.ptD[2])
    outDoubleArray.append(self.ptAverage[0])
    outDoubleArray.append(self.ptAverage[1])
    outDoubleArray.append(self.ptAverage[2])
    outIntegerArray.append(self.cellId)

  def AddZeroQuadToInfoArrays(self, outDoubleArray, outIntegerArray):
    """add zeros to the argument arrays representing a blank quad
       for purspose of doing mpi broadcast/reduce work"""
    for ii in range(0, 15):
        outDoubleArray.append(0.0)
    outIntegerArray.append(int(0))

  def CalculateDistanceSquaredTo(self, onePt):
    return vecDistanceSquared(onePt, self.ptAverage)

  def ToString(self):
    retStr = str(self.index) + " " + str(self.cellId) + "\n" + \
      str(self.ptA) + "\n" + \
      str(self.ptB) + "\n" + \
      str(self.ptC) + "\n" + \
      str(self.ptD) + "\n" + \
      str(self.ptAverage) + "\n"
    return retStr

  def CalculateNormal(self):
    vec1 = vecFromAToB(self.ptA, self.ptC)
    vec2 = vecFromAToB(self.ptB, self.ptD)
    vecCrs = vecCrossProduct(vec1, vec2)
    self.normal = vecNormalize(vecCrs)

  def MakeLogSamplePointsAlongRayStartingAtAveragePoint(self, rayDirection, startStep, multiplier, maxDistance):
    smplList = []
    newPt = vecMultiplyAdd(self.ptAverage, rayDirection, 0.0)
    smplList.append(newPt)
    thisDistance = startStep
    count = 1
    while thisDistance < maxDistance:
      newPt = vecMultiplyAdd(self.ptAverage, rayDirection, thisDistance)
      smplList.append(newPt)
      thisDistance *= multiplier
      count += 1
      if count > 200:
        break
    if PhactoriDbg(100):
      myDebugPrint3("MakeLogSamplePointsAlongRayStartingAtAveragePoint:\n" + \
        str(len(smplList)) + "  " + \
        str(thisDistance) + "  " + \
        str(count) + \
        "\n" + str(smplList[0]) + \
        "\n" + str(smplList[-1]) + \
        "\n", 100)
    return smplList

  def MakeSamplePointsAlongRayStartingAtAveragePoint(self, rayDirection, numSamples, maxDistance):
    smplList = []
    #for ii in range(-50,51):
    sampleStep = maxDistance / (float(numSamples - 1))
    for ii in range(0, numSamples):
      rr = float(ii) * sampleStep
      newPt = vecMultiplyAdd(self.ptAverage, rayDirection, rr)
      smplList.append(newPt)
    return smplList

  def MakeSamplePointsAlongRayFromNearbyPoint(self, rayDefiningPoint):
    """create a vector from rayDefiningPoint to the average point of this ray,
       VV.  Create sample points along a ray starting at the average point of
       the quad (self.ptAverage) along the vector VV"""
    vecVUN = vecFromAToB(rayDefiningPoint, self.ptAverage)
    vecVV= vecNormalize(vecVUN)
    return self.MakeSamplePointsAlongRayStartingAtAveragePoint(vecVV, 50, 0.05)

  def MakeLogSamplePointsAlongRayFromNearbyPoint(self, rayDefiningPoint, startStep, multiplier, maxDistance):
    """create a vector from rayDefiningPoint to the average point of this ray,
       VV.  Create sample points along a ray starting at the average point of
       the quad (self.ptAverage) along the vector VV.  Space the sample points
       in a logarithmic manner controlled by startStep, multiplier, and
       maxDistance"""
    vecVUN = vecFromAToB(rayDefiningPoint, self.ptAverage)
    vecVV= vecNormalize(vecVUN)
    return self.MakeLogSamplePointsAlongRayStartingAtAveragePoint(vecVV, startStep, multiplier, maxDistance)

  def MakePerpendicularSamplePoints(self):
    """create a list of sample points based on the average quad point (later
       the center) and thr normal to the quad"""
    self.CalculateNormal()
    return self.MakeSamplePointsAlongRayStartingAtAveragePoint(self.normal, 50, 0.05)

#phactori_combine_to_single_python_file_subpiece_end_1

