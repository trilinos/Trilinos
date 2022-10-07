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

#from phactori import *
import math
from .PhactoriVectorLibrary import *
from phactori import myDebugPrint3AndException

#phactori_combine_to_single_python_file_subpiece_begin_1
class PhactoriSegment:
  def __init__(self):
    self.ptA = None
    self.ptB = None
    self.AtoB = None
    self.lengthSquared = None
    self.oneOverLenSqrd = None
    self.boundingBox = None

  def SetPoints(self, ptA, ptB):
    self.ptA = list(ptA)
    self.ptB = list(ptB)
    self.AtoB = vecFromAToB(self.ptA, self.ptB)
    self.lengthSquared = vecMagnitudeSquared(self.AtoB)
    if self.lengthSquared == 0.0:
      myDebugPrint3AndException("PhactoriSegment:SetPoints self.lengthSquared is 0.0\n")
    self.oneOverLenSqrd = 1.0 / self.lengthSquared

    self.boundingBox = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if self.ptA[0] <= self.ptB[0]:
      self.boundingBox[0] = ptA[0]
      self.boundingBox[1] = ptB[0]
    else:
      self.boundingBox[0] = ptB[0]
      self.boundingBox[1] = ptA[0]
    if self.ptA[1] <= self.ptB[1]:
      self.boundingBox[2] = ptA[1]
      self.boundingBox[3] = ptB[1]
    else:
      self.boundingBox[2] = ptB[1]
      self.boundingBox[3] = ptA[1]
    if self.ptA[2] <= self.ptB[2]:
      self.boundingBox[4] = ptA[2]
      self.boundingBox[5] = ptB[2]
    else:
      self.boundingBox[4] = ptB[2]
      self.boundingBox[5] = ptA[2]

    #if PhactoriDbg(100):
    #  myDebugPrint3("PhactoriSegment:SetPoints\n" +\
    #    str(self.ptA) + "\n" + \
    #    str(self.ptB) + "\n" + \
    #    str(self.AtoB) + "\n" + \
    #    str(self.lengthSquared) + "\n" + \
    #    str(self.oneOverLenSqrd) + "\n")

  def FindDistanceToPoint(self, testPoint):
    return math.sqrt(self.FindDistanceSquaredToPoint(testPoint))

  def FindDistanceToPointProjected(self, testPoint, projectionAxis):
    return math.sqrt(self.FindDistanceSquaredToPoint(testPoint), projectionAxis)

  def FindNearestPointOnSegmentToPoint(self, testPoint):
    if self.lengthSquared == 0.0:
      return list(self.ptA)
    AtoTest = vecFromAToB(self.ptA, testPoint)
    tt = vecDotProduct(AtoTest, self.AtoB) * self.oneOverLenSqrd
    if tt <= 0.0:
      return list(self.ptA)
    elif tt >= 1.0:
      return list(self.ptB)
    else:
      testProjectedToSegment = vecMultiplyAdd(self.ptA, self.AtoB, tt)
      return testProjectedToSegment

  def ProjectPointOntoPlane(self, pointToChange, projectionAxis):
    if projectionAxis == 0:
      pointToChange[0] = 0.0
    elif projectionAxis == 1:
      pointToChange[1] = 0.0
    elif projectionAxis == 2:
      pointToChange[2] = 0.0

  def FindNearestPointOnSegmentToPointProjected(self, testPoint, projectionAxis):
    retPoint = self.FindNearestPointOnSegmentToPoint(testPoint)
    self.ProjectPointOntoPlane(retPoint, projectionAxis)
    return retPoint

  def FindDistanceSquaredToPoint(self, testPoint):
    testPointProjectedOntoSegment = self.FindNearestPointOnSegmentToPoint(testPoint)
    return vecDistanceSquared(testPoint, testPointProjectedOntoSegment)

  def FindDistanceSquaredToPointProjected(self, testPoint, projectionAxis):
    testPointProjectedOntoSegment = self.FindNearestPointOnSegmentToPointProjected(testPoint, projectionAxis)
    localTestPoint = list(testPoint)
    self.ProjectPointOntoPlane(localTestPoint, projectionAxis)
    return vecDistanceSquared(localTestPoint, testPointProjectedOntoSegment)

  def MyBoundingBoxIntersectsTargetBoundingBox(self, testBoundingBox):
    if self.boundingBox[0] > testBoundingBox[1]:
      return False
    if self.boundingBox[1] < testBoundingBox[0]:
      return False
    if self.boundingBox[2] > testBoundingBox[3]:
      return False
    if self.boundingBox[3] < testBoundingBox[2]:
      return False
    if self.boundingBox[4] > testBoundingBox[5]:
      return False
    if self.boundingBox[5] < testBoundingBox[4]:
      return False
    return True

  def EitherEndpointIsInBoundingBox(self, testBoundingBox):
    ptIsIn = True
    if (self.ptA[0] < testBoundingBox[0]):
      ptIsIn = False
    if (self.ptA[0] > testBoundingBox[1]):
      ptIsIn = False
    if (self.ptA[1] < testBoundingBox[2]):
      ptIsIn = False
    if (self.ptA[1] > testBoundingBox[3]):
      ptIsIn = False
    if (self.ptA[2] < testBoundingBox[4]):
      ptIsIn = False
    if (self.ptA[2] > testBoundingBox[5]):
      ptIsIn = False
    if ptIsIn == True:
      return True
    ptIsIn = True
    if (self.ptB[0] < testBoundingBox[0]):
      ptIsIn = False
    if (self.ptB[0] > testBoundingBox[1]):
      ptIsIn = False
    if (self.ptB[1] < testBoundingBox[2]):
      ptIsIn = False
    if (self.ptB[1] > testBoundingBox[3]):
      ptIsIn = False
    if (self.ptB[2] < testBoundingBox[4]):
      ptIsIn = False
    if (self.ptB[2] > testBoundingBox[5]):
      ptIsIn = False
    if ptIsIn == True:
      return True
    return False

  def IntersectsXYRectangle(self, minU, maxU, minV, maxV, posW):
    #myDebugPrint3("IntersectsXYRectangle entered\n")
    #myDebugPrint3(str([minU, maxU, minV, maxV, posW]) + "\n")
    #myDebugPrint3(str(self.boundingBox) + "\n")
    if self.boundingBox[4] > posW:
      return False
    if self.boundingBox[5] < posW:
      return False
    tt = (posW - self.ptA[2])/(self.ptB[2] - self.ptA[2])
    #myDebugPrint3(str(tt) + "\n")
    uIntersect = self.ptA[0] + tt * (self.ptB[0] - self.ptA[0])
    #myDebugPrint3(str(uIntersect) + "\n")
    if uIntersect < minU:
      return False
    if uIntersect > maxU:
      return False
    vIntersect = self.ptA[1] + tt * (self.ptB[1] - self.ptA[1])
    #myDebugPrint3(str(vIntersect) + "\n")
    if vIntersect < minV:
      return False
    if vIntersect > maxV:
      return False
    return True

  def IntersectsXZRectangle(self, minU, maxU, minV, maxV, posW):
    if self.boundingBox[2] > posW:
      return False
    if self.boundingBox[3] < posW:
      return False
    tt = (posW - self.ptA[1])/(self.ptB[1] - self.ptA[1])
    uIntersect = self.ptA[0] + tt * (self.ptB[0] - self.ptA[0])
    if uIntersect < minU:
      return False
    if uIntersect > maxU:
      return False
    vIntersect = self.ptA[2] + tt * (self.ptB[2] - self.ptA[2])
    if vIntersect < minV:
      return False
    if vIntersect > maxV:
      return False
    return True

  def IntersectsYZRectangle(self, minU, maxU, minV, maxV, posW):
    if self.boundingBox[0] > posW:
      return False
    if self.boundingBox[1] < posW:
      return False
    tt = (posW - self.ptA[0])/(self.ptB[0] - self.ptA[0])
    uIntersect = self.ptA[1] + tt * (self.ptB[1] - self.ptA[1])
    if uIntersect < minU:
      return False
    if uIntersect > maxU:
      return False
    vIntersect = self.ptA[2] + tt * (self.ptB[2] - self.ptA[2])
    if vIntersect < minV:
      return False
    if vIntersect > maxV:
      return False
    return True

  def IntersectsBoundingBox(self, testBoundingBox):
    #myDebugPrint3("PhactoriSegment.IntersectsBoundingBox entered\n")
    #myDebugPrint3("ptA: " + str(self.ptA) + "\n")
    #myDebugPrint3("ptB: " + str(self.ptB) + "\n")
    #myDebugPrint3("box: " + str(testBoundingBox) + "\n")

    if self.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox) == False:
      #myDebugPrint3("PhactoriSegment.IntersectsBoundingBox returning false due to bounding box test\n")
      return False

    if self.EitherEndpointIsInBoundingBox(testBoundingBox) == True:
      #myDebugPrint3("PhactoriSegment.IntersectsBoundingBox returning false due to endpoint in box test\n")
      return True

    #myDebugPrint3("PhactoriSegment.IntersectsBoundingBox testing against all six sides\n")
    if self.IntersectsXYRectangle(testBoundingBox[0], testBoundingBox[1], testBoundingBox[2], testBoundingBox[3], testBoundingBox[4]):
      return True
    if self.IntersectsXYRectangle(testBoundingBox[0], testBoundingBox[1], testBoundingBox[2], testBoundingBox[3], testBoundingBox[5]):
      return True
    if self.IntersectsXZRectangle(testBoundingBox[0], testBoundingBox[1], testBoundingBox[4], testBoundingBox[5], testBoundingBox[2]):
      return True
    if self.IntersectsXZRectangle(testBoundingBox[0], testBoundingBox[1], testBoundingBox[4], testBoundingBox[5], testBoundingBox[3]):
      return True
    if self.IntersectsYZRectangle(testBoundingBox[2], testBoundingBox[3], testBoundingBox[4], testBoundingBox[5], testBoundingBox[4]):
      return True
    if self.IntersectsYZRectangle(testBoundingBox[2], testBoundingBox[3], testBoundingBox[4], testBoundingBox[5], testBoundingBox[5]):
      return True

    #myDebugPrint3("PhactoriSegment.IntersectsBoundingBox returning false, no side intersection\n")
    return False

  def IntersectsBoundingBoxProjected(self, boundingBox, projectionAxis):
    #myDebugPrint3("PhactoriSegment.IntersectsBoundingBoxProjected entered\n")
    if (projectionAxis < 0) or (projectionAxis > 2):
      return self.IntersectsBoundingBox(boundingBox)

    #inefficient but simple method of handling projection:  make projected
    #bounding box (rectangle) and projected segment and do 3d version
    localProjectedPointA = list(self.ptA)
    localProjectedPointB = list(self.ptB)
    localProjectedBoundingBox = list(boundingBox)
    localProjectedPointA[projectionAxis] = 0.0
    localProjectedPointB[projectionAxis] = 0.0
    localProjectedBoundingBox[2*projectionAxis] = 0.0
    localProjectedBoundingBox[2*projectionAxis+1] = 0.0
    localProjectedSegment = PhactoriSegment()
    localProjectedSegment.SetPoints(localProjectedPointA, localProjectedPointB)
    intersectionFlag = localProjectedSegment.IntersectsBoundingBox(localProjectedBoundingBox)
    return intersectionFlag

#phactori_combine_to_single_python_file_subpiece_end_1
