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
import unittest
from Operation.PhactoriSegment import *
from paraview.simple import *

class TestPhactoriSegment(unittest.TestCase):

  def test_SetPoints(self):
    theSegment = PhactoriSegment()
    theSegment.SetPoints([-1.0, -1.0, -1.0], [0.0, 1.0, 2.0])
    self.assertEqual(theSegment.AtoB[0], 1.0)
    self.assertEqual(theSegment.AtoB[1], 2.0)
    self.assertEqual(theSegment.AtoB[2], 3.0)
    self.assertEqual(theSegment.lengthSquared, 14.0)
    self.assertEqual(theSegment.oneOverLenSqrd, 1.0/14.0)

  def test_FindNearestPointOnSegmentToPoint(self):
    theSegment = PhactoriSegment()
    theSegment.SetPoints([-1.0, -1.0, -1.0], [1.0, 2.0, 3.0])
    testpt1 = theSegment.FindNearestPointOnSegmentToPoint([-10.0, -10.0, -10.0])
    self.assertEqual(testpt1, theSegment.ptA)
    testpt1 = theSegment.FindNearestPointOnSegmentToPoint([10.0, 10.0, 10.0])
    self.assertEqual(testpt1, theSegment.ptB)
    testpt1 = theSegment.FindNearestPointOnSegmentToPoint([0.0, 0.0, 0.0])
    self.assertEqual(testpt1, [-0.3793103448275862, -0.06896551724137923, 0.24137931034482762])

  def test_FindDistanceSquaredToPoint(self):
    theSegment = PhactoriSegment()
    theSegment.SetPoints([-1.0, -1.0, -1.0], [1.0, 2.0, 3.0])
    testDist = theSegment.FindDistanceSquaredToPoint([0.5, 0.6, 0.7])
    self.assertEqual(testDist, 0.3496551724137931)

  def test_MyBoundingBoxIntersectsTargetBoundingBox(self):
    theSegment = PhactoriSegment()
    testBoundingBox = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]
    theSegment.SetPoints([-3.0, -3.0, -3.0], [-2.0, -2.0, -2.0])
    self.assertFalse(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
    theSegment.SetPoints([-3.0, -3.0, -3.0], [2.0, 2.0, 2.0])
    self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
    theSegment.SetPoints([-3.0, -3.0, -3.0], [0.0, 0.0, 0.0])
    self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
    theSegment.SetPoints([0.0, 0.0, 0.0], [3.0, 2.0, 1.0])
    self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))

    for axesTestIndex in range(0,3):
      tp1 = [0.0, 0.0, 0.0]
      tp2 = [0.0, 0.0, 0.0]
      tp1[axesTestIndex] = -3.0
      tp2[axesTestIndex] = -1.01
      theSegment.SetPoints(tp1, tp2)
      self.assertFalse(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertFalse(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = 3.0
      tp2[axesTestIndex] = 1.01
      theSegment.SetPoints(tp1, tp2)
      self.assertFalse(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertFalse(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = -3.0
      tp2[axesTestIndex] = -1.0
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = 3.0
      tp2[axesTestIndex] = 1.0
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = -3.0
      tp2[axesTestIndex] = -0.99
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = 3.0
      tp2[axesTestIndex] = 0.99
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.MyBoundingBoxIntersectsTargetBoundingBox(testBoundingBox))

  def test_EitherEndpointIsInBoundingBox(self):
    theSegment = PhactoriSegment()
    testBoundingBox = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]
    theSegment.SetPoints([-3.0, -3.0, -3.0], [-2.0, -2.0, -2.0])
    self.assertFalse(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
    theSegment.SetPoints([-3.0, -3.0, -3.0], [0.1, 0.1, 0.1])
    self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
    theSegment.SetPoints([0.1, 0.1, 0.1], [3.0, 3.0, 3.0])
    self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
    theSegment.SetPoints([0.1, 0.1, 0.1], [-0.1, -0.1, -0.1])
    self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))

    for axesTestIndex in range(0,3):
      tp1 = [0.0, 0.0, 0.0]
      tp2 = [0.0, 0.0, 0.0]
      tp1[axesTestIndex] = -3.0
      tp2[axesTestIndex] = -1.01
      theSegment.SetPoints(tp1, tp2)
      self.assertFalse(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertFalse(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = 3.0
      tp2[axesTestIndex] = 1.01
      theSegment.SetPoints(tp1, tp2)
      self.assertFalse(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertFalse(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = -3.0
      tp2[axesTestIndex] = -1.0
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = 3.0
      tp2[axesTestIndex] = 1.0
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = -3.0
      tp2[axesTestIndex] = -0.99
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      tp1[axesTestIndex] = 3.0
      tp2[axesTestIndex] = 0.99
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.EitherEndpointIsInBoundingBox(testBoundingBox))

  def test_IntersectsBoundingBoxProjected(self):
    theSegment = PhactoriSegment()
    testBoundingBox = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]
    theSegment.SetPoints([-3.0, -3.0, -3.0], [-2.0, -2.0, -2.0])
    self.assertFalse(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
    theSegment.SetPoints([-2.0, -2.0, -2.0], [2.0, 2.0, 2.0])
    self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
    theSegment.SetPoints([-3.0, -3.0, -3.0], [2.0, 2.0, 2.0])
    self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))

    for axesTestIndex in range(0,3):
      tp1 = [-2.0, -2.0, -2.0]
      tp2 = [ 2.0,  2.0,  2.0]
      tp1[axesTestIndex] = 0.0
      tp2[axesTestIndex] = 0.0
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      tp1 = [-1.875, -1.925, -2.125]
      tp2 = [ 2.125,  1.875,  2.5]
      tp1[axesTestIndex] = 0.125
      tp2[axesTestIndex] = -0.125
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      tp1 = [-1.1, -1.1, -1.1]
      tp2 = [1.1, 1.1, 1.1]
      tp2[axesTestIndex] = -0.9
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      tp2[axesTestIndex] = -0.999
      theSegment.SetPoints(tp2, tp1)
      self.assertFalse(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      theSegment.SetPoints(tp1, tp2)
      self.assertFalse(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      tp1 = [-1.1, -1.1, -1.1]
      tp2 = [1.1, 1.1, 1.1]
      tp1[axesTestIndex] = 0.9
      theSegment.SetPoints(tp2, tp1)
      self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      theSegment.SetPoints(tp1, tp2)
      self.assertTrue(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      tp1[axesTestIndex] = 0.999
      theSegment.SetPoints(tp2, tp1)
      self.assertFalse(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))
      theSegment.SetPoints(tp1, tp2)
      self.assertFalse(theSegment.IntersectsBoundingBoxProjected(testBoundingBox, -1))


if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()
