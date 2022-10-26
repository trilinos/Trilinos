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

import unittest
#from vtk import *
import vtk
from Operation.PhactoriCreateSegmentsNormalToCells import *
import paraview.simple

class TestPhactoriCreateSegmentsNormalToCells(unittest.TestCase):

  def MakeOneHexahdtronGridFrom8Points1(self, pointList):

    # Create the points
    points = vtk.vtkPoints()
    points.InsertNextPoint(pointList[0])
    points.InsertNextPoint(pointList[1])
    points.InsertNextPoint(pointList[2])
    points.InsertNextPoint(pointList[3])
    points.InsertNextPoint(pointList[4])
    points.InsertNextPoint(pointList[5])
    points.InsertNextPoint(pointList[6])
    points.InsertNextPoint(pointList[7])

    # Create a hexahedron from the points
    hex = vtk.vtkHexahedron()
    hex.GetPointIds().SetId(0,0)
    hex.GetPointIds().SetId(1,1)
    hex.GetPointIds().SetId(2,2)
    hex.GetPointIds().SetId(3,3)
    hex.GetPointIds().SetId(4,4)
    hex.GetPointIds().SetId(5,5)
    hex.GetPointIds().SetId(6,6)
    hex.GetPointIds().SetId(7,7)

    # Add the hexahedron to a cell array
    hexs = vtk.vtkCellArray()
    hexs.InsertNextCell(hex)

    # Add the points and hexahedron to an unstructured grid
    uGrid = vtk.vtkUnstructuredGrid()
    uGrid.SetPoints(points)
    uGrid.InsertNextCell(hex.GetCellType(), hex.GetPointIds())
    return uGrid

  def MakeOneHexahdronGrid1(self):
    myPtList = []
    myPtList.append([0.0, 0.0, 0.0])
    myPtList.append([3.0, 0.0, 0.0])
    myPtList.append([3.0, 4.0, 0.0])
    myPtList.append([0.0, 4.0, 0.0])
    myPtList.append([0.0, 0.0, 1.0])
    myPtList.append([3.0, 0.0, 2.0])
    myPtList.append([3.0, 4.0, 3.0])
    myPtList.append([0.0, 4.0, 4.0])

    return self.MakeOneHexahdtronGridFrom8Points1(myPtList)

  def MakeOneFlatHexahdronGrid1(self):
    myPtList = []
    myPtList.append([0.0, 0.0, 0.125])
    myPtList.append([3.0, 0.0, 0.25])
    myPtList.append([3.125, 0.03125, 0.375])
    myPtList.append([0.375, 0.03125, 0.5])
    myPtList.append([0.5, 0.0, 1.0])
    myPtList.append([3.0675, 0.0, 2.125])
    myPtList.append([3.325, 0.03125, 3.25])
    myPtList.append([0.625, 0.03125, 4.375])

    return self.MakeOneHexahdtronGridFrom8Points1(myPtList)

  def MakeOneFlatHexahdronGrid2(self):
    myPtList = []
    myPtList.append([0.0, 0.0, 0.0])
    myPtList.append([1.0, 0.0, 0.0])
    myPtList.append([1.0, 0.03125, 0.0])
    myPtList.append([0.0, 0.03125, 0.0])
    myPtList.append([0.0, 0.0, 1.0])
    myPtList.append([1.0, 0.0, 1.0])
    myPtList.append([1.0, 0.03125, 1.0])
    myPtList.append([0.0, 0.03125, 1.0])

    return self.MakeOneHexahdtronGridFrom8Points1(myPtList)

  def test_Parse1(self):
    pcsntc1 = PhactoriCreateSegmentsNormalToCells()
    testJson1 = {"segment length":5.0}
    pcsntc1.ParseParametersFromJson(testJson1)
    self.assertEqual(pcsntc1.SegmentLength, 5.0)

  def test_Parse2(self):
    pcsntc1 = PhactoriCreateSegmentsNormalToCells()
    testJson1 = \
      {"global x axis":[1.25, 1.5, 1.75],
       "global x axis backup":[1.5, -1.5, 0.0]}
    pcsntc1.ParseParametersFromJson(testJson1)
    self.assertEqual(pcsntc1.localAxisXReference,
      [0.4767312946227962, 0.5720775535473553, 0.6674238124719146])
    self.assertEqual(pcsntc1.localAxisXReference2,
      [0.7071067811865476, -0.7071067811865476, 0.0])
    testJson1 = {"global x axis":[1.0, 0.5, 0.0]}
    pcsntc1.ParseParametersFromJson(testJson1)
    self.assertEqual(pcsntc1.localAxisXReference,
      [0.8944271909999159, 0.4472135954999579, 0.0])
    self.assertEqual(pcsntc1.localAxisXReference2, [0.0, 0.0, 1.0])
    testJson1 = {"global x axis":[1.0, 0.0, 0.5]}
    pcsntc1.ParseParametersFromJson(testJson1)
    self.assertEqual(pcsntc1.localAxisXReference2, [0.0, 1.0, 0.0])
    testJson1 = {"global x axis":[0.0, 1.0, 0.5]}
    pcsntc1.ParseParametersFromJson(testJson1)
    self.assertEqual(pcsntc1.localAxisXReference2, [1.0, 0.0, 0.0])

  def test_Parse3(self):
    pcsntc1 = PhactoriCreateSegmentsNormalToCells()
    testJson1 = {"global x axis":[1.25, 1.5]}
    self.assertRaises(Exception, pcsntc1.ParseParametersFromJson, testJson1)
    testJson1 = \
      {"global x axis":[1.25, 1.5, 1.75],
       "global x axis backup":[1.5, -1.5]}
    self.assertRaises(Exception, pcsntc1.ParseParametersFromJson, testJson1)
    testJson1 = {"global x axis":"blah"}
    self.assertRaises(Exception, pcsntc1.ParseParametersFromJson, testJson1)
    testJson1 = \
      {"global x axis":[1.25, 1.5, 1.75],
       "global x axis backup":1.5}
    self.assertRaises(Exception, pcsntc1.ParseParametersFromJson, testJson1)

  def test_Parse4(self):
    pcsntc1 = PhactoriCreateSegmentsNormalToCells()
    testJson1 = {"segment direction reference point":[1.0, -2.0, 3.5]}
    pcsntc1.ParseParametersFromJson(testJson1)
    self.assertEqual(pcsntc1.DirectionReferencePoint, [1.0, -2.0, 3.5])
    testJson1 = {"segment direction reference point":[-2.0, 3.5]}
    self.assertRaises(Exception, pcsntc1.ParseParametersFromJson, testJson1)
    testJson1 = {"segment direction reference point":13.5}
    self.assertRaises(Exception, pcsntc1.ParseParametersFromJson, testJson1)


  def test_CalculateGeometricProgressionSampleValues_linear(self):
    pcsntc1 = PhactoriCreateSegmentsNormalToCells()
    pcsntc1.SegmentLength = 9.0
    pcsntc1.SegmentResolution = 10
    pcsntc1.geometricSamplingMultiplier = 1.0
    pcsntc1.segmentSampleMethod = "linear"
    oneSeg = PcsntcSegment(5, 2, [1.0, 1.0, 1.0], [1.0, -1.0, 2.0])
    oneSeg.segmentLength = pcsntc1.SegmentLength
    testVals = pcsntc1.CalculateGeometricProgressionSampleValues(oneSeg)
    self.assertEqual(testVals, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])

  def test_CalculateGeometricProgressionSampleValues_geometric(self):
    pcsntc1 = PhactoriCreateSegmentsNormalToCells()
    pcsntc1.SegmentLength = 32767.0 * 0.125
    pcsntc1.SegmentResolution = 16
    pcsntc1.geometricSamplingMultiplier = 2.0
    pcsntc1.segmentSampleMethod = "geometric"
    oneSeg = PcsntcSegment(5, 2, [1.0, 1.0, 1.0], [1.0, -1.0, 2.0])
    oneSeg.segmentLength = pcsntc1.SegmentLength
    testVals = pcsntc1.CalculateGeometricProgressionSampleValues(oneSeg)
    goldVals = [0.0, 0.125, 0.375, 0.875, 1.875, 3.875, 7.875, 15.875, 31.875, 63.875, 127.875, 255.875, 511.875, 1023.875, 2047.875, 4095.875]
    self.assertEqual(testVals, goldVals)

  def test_CreateSegmentsForAllCells(self):
    wv1 = paraview.simple.Wavelet()
    wv1.WholeExtent = [-2, 2, -2, 2, -2, 2]
    wv1.UpdatePipeline()

    pcsntc1 = PhactoriCreateSegmentsNormalToCells()
    pcsntc1.CreateSegmentsForAllCells(wv1)
    firstSegment = pcsntc1.RecursionResults.SegmentCollection[0]
    #self.assertEqual((1, 0.0, [-1.5, -1.5, -1.5], [1.0, 0.0, -0.0]), firstSegment)
    self.assertEqual(1, firstSegment.blockIndex)
    self.assertEqual(0, firstSegment.cellIndex)
    self.assertEqual([-1.5, -1.5, -1.5], firstSegment.ptA)
    self.assertEqual([1.0, 0.0, -0.0], firstSegment.segDir)

    lastSegment = pcsntc1.RecursionResults.SegmentCollection[-1]
    #self.assertEqual((1, 63, [1.5, 1.5, 1.5], [1.0, 0.0, -0.0]), lastSegment)
    self.assertEqual(1, lastSegment.blockIndex)
    self.assertEqual(63, lastSegment.cellIndex)
    self.assertEqual([1.5, 1.5, 1.5], lastSegment.ptA)
    self.assertEqual([1.0, 0.0, -0.0], lastSegment.segDir)

    xfrm1 = Transform(Input=wv1)
    xfrm1.Transform = 'Transform'
    xfrm1.Transform.Scale = [1.0, 0.005, 1.0]
    xfrm1.UpdatePipeline()
    #wv1grid = xfrm1.GetClientSideObject().GetOutputDataObject(0)
    #wv1grid = wv1.GetClientSideObject().GetOutputDataObject(0)
    pcsntc1.CreateSegmentsForAllCells(xfrm1)
    firstSegment = pcsntc1.RecursionResults.SegmentCollection[0]
    #self.assertEqual((1, 0, [-1.5, -0.0075, -1.5], [0.0, 1.0, 0.0]), firstSegment)
    self.assertEqual(1, firstSegment.blockIndex)
    self.assertEqual(0, firstSegment.cellIndex)
    self.assertEqual([-1.5, -0.0075, -1.5], firstSegment.ptA)
    self.assertEqual([0.0, 1.0, 0.0], firstSegment.segDir)
    lastSegment = pcsntc1.RecursionResults.SegmentCollection[-1]
    #self.assertEqual((1, 63, [1.5, 0.0075, 1.5], [0.0, 1.0, 0.0]), lastSegment)
    self.assertEqual(1, lastSegment.blockIndex)
    self.assertEqual(63, lastSegment.cellIndex)
    self.assertEqual([1.5, 0.0075, 1.5], lastSegment.ptA)
    self.assertEqual([0.0, 1.0, 0.0], lastSegment.segDir)

  def test_CreateParaViewSourcesForSegments(self):
    wv1 = paraview.simple.Wavelet()
    wv1.WholeExtent = [-2, 2, -2, 2, -2, 2]
    wv1.UpdatePipeline()
    xfrm1 = Transform(Input=wv1)
    xfrm1.Transform = 'Transform'
    xfrm1.Transform.Scale = [1.0, 0.005, 1.0]
    xfrm1.UpdatePipeline()
    #wv1grid = xfrm1.GetClientSideObject().GetOutputDataObject(0)
    #wv1grid = wv1.GetClientSideObject().GetOutputDataObject(0)
    pcsntc1 = PhactoriCreateSegmentsNormalToCells()
    pcsntc1.CreateSegmentsForAllCells(xfrm1)
    pcsntc1.CreateParaViewSourcesForSegments()


if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()
