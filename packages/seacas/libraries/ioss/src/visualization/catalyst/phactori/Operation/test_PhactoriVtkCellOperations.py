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
from Operation.PhactoriVtkCellOperations import *
import paraview.simple

class TestPhactoriVtkCellOperations(unittest.TestCase):

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

  def test_PhactoriCalculateFaceArea(self):
    testugrid1 = self.MakeOneHexahdronGrid1()
    testhex1 = testugrid1.GetCell(0)
    testface0 = testhex1.GetFace(0)
    testarea0 = PhactoriCalculateFaceArea(testugrid1, testface0)
    self.assertEqual(testarea0, 10.0)
    testface0 = testhex1.GetFace(1)
    testarea0 = PhactoriCalculateFaceArea(testugrid1, testface0)
    self.assertEqual(testarea0, 10.0)
    testface0 = testhex1.GetFace(2)
    testarea0 = PhactoriCalculateFaceArea(testugrid1, testface0)
    self.assertEqual(testarea0, 4.5)
    testface0 = testhex1.GetFace(3)
    testarea0 = PhactoriCalculateFaceArea(testugrid1, testface0)
    self.assertEqual(testarea0, 10.5)
    testface0 = testhex1.GetFace(4)
    testarea0 = PhactoriCalculateFaceArea(testugrid1, testface0)
    self.assertEqual(testarea0, 12.0)
    testface0 = testhex1.GetFace(5)
    testarea0 = PhactoriCalculateFaceArea(testugrid1, testface0)
    self.assertEqual(testarea0, 14.262087348130013)

  def test_PhactoriCalculateCellFaceAreas(self):
    testugrid1 = self.MakeOneHexahdronGrid1()
    testhex1 = testugrid1.GetCell(0)
    cellFaceAreaList = PhactoriCalculateCellFaceAreas(testugrid1, testhex1)
    self.assertEqual(cellFaceAreaList, [10.0, 10.0, 4.5, 10.5, 12.0, 14.262087348130013])

  def test_PhactoriCountCellTouchingEachFace(self):
    wv1 = paraview.simple.Wavelet()
    wv1.UpdatePipeline()
    wv1grid = wv1.GetClientSideObject().GetOutputDataObject(0)
    cellsTouchingEachFace = PhactoriCountCellTouchingEachFace(wv1grid)
    self.assertEqual(cellsTouchingEachFace.totalNumFaces, 25200)
    #corners
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(0)
    self.assertEqual(counts1, [1, 2, 1, 2, 1, 2])
    exteriorList1 = cellsTouchingEachFace.GetListOfExteriorFacesOnCellByIndex(0)
    self.assertEqual(exteriorList1, [0,2,4])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(20*20*20-1)
    self.assertEqual(counts1, [2, 1, 2, 1, 2, 1])
    exteriorList1 = cellsTouchingEachFace.GetListOfExteriorFacesOnCellByIndex(20*20*20-1)
    self.assertEqual(exteriorList1, [1,3,5])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(19)
    self.assertEqual(counts1, [2, 1, 1, 2, 1, 2])
    #edges, not corner
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(1)
    self.assertEqual(counts1, [2, 2, 1, 2, 1, 2])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(20*20*20-2)
    self.assertEqual(counts1, [2, 2, 2, 1, 2, 1])
    #surface, not corner or edge
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(30)
    self.assertEqual(counts1, [2, 2, 2, 2, 1, 2])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(402)
    self.assertEqual(counts1, [2, 2, 1, 2, 2, 2])
    exteriorList1 = cellsTouchingEachFace.GetListOfExteriorFacesOnCellByIndex(402)
    self.assertEqual(exteriorList1, [2])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(420)
    self.assertEqual(counts1, [1, 2, 2, 2, 2, 2])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(439)
    self.assertEqual(counts1, [2, 1, 2, 2, 2, 2])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(21)
    self.assertEqual(counts1, [2, 2, 2, 2, 1, 2])
    exteriorList1 = cellsTouchingEachFace.GetListOfExteriorFacesOnCellByIndex(21)
    self.assertEqual(exteriorList1, [4])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell((20*20*20-1)-29)
    self.assertEqual(counts1, [2, 2, 2, 2, 2, 1])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell((20*20*20-1)-402)
    self.assertEqual(counts1, [2, 2, 2, 1, 2, 2])
    #interior
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(421)
    self.assertEqual(counts1, [2, 2, 2, 2, 2, 2])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(20*20*10+35)
    self.assertEqual(counts1, [2, 2, 2, 2, 2, 2])
    exteriorList1 = cellsTouchingEachFace.GetListOfExteriorFacesOnCellByIndex(20*20*10+35)
    self.assertEqual(exteriorList1, [])
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell((20*20*20-1)-422)
    self.assertEqual(counts1, [2, 2, 2, 2, 2, 2])

    #testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-29, cellsTouchingEachPoint)
    #self.assertEqual(testStatus, 2)
    #testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-402, cellsTouchingEachPoint)
    #self.assertEqual(testStatus, 8)

  def test_PhactoriCountCellTouchingEachPoint(self):
    wv1 = paraview.simple.Wavelet()
    wv1.UpdatePipeline()
    wv1grid = wv1.GetClientSideObject().GetOutputDataObject(0)
    cellsTouchingEachPoint = PhactoriCountCellTouchingEachPoint(wv1grid)

    self.assertEqual(len(cellsTouchingEachPoint), 21*21*21)
    #corners
    self.assertEqual(cellsTouchingEachPoint[0], 1)
    self.assertEqual(cellsTouchingEachPoint[-1], 1)
    self.assertEqual(cellsTouchingEachPoint[21-1], 1)
    self.assertEqual(cellsTouchingEachPoint[21*21-1], 1)
    #edges
    self.assertEqual(cellsTouchingEachPoint[1], 2)
    self.assertEqual(cellsTouchingEachPoint[-2], 2)
    #surface
    self.assertEqual(cellsTouchingEachPoint[40], 4)
    self.assertEqual(cellsTouchingEachPoint[-40], 4)
    #interior
    self.assertEqual(cellsTouchingEachPoint[1000], 8)
    self.assertEqual(cellsTouchingEachPoint[-1000], 8)

  def test_PhactoriFindSurfaceStatusForOneCell_CubicCells(self):
    wv1 = paraview.simple.Wavelet()
    wv1.UpdatePipeline()
    wv1grid = wv1.GetClientSideObject().GetOutputDataObject(0)
    cellsTouchingEachPoint = PhactoriCountCellTouchingEachPoint(wv1grid)
    #corner cells
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 0, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 10)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 19, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 10)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1), cellsTouchingEachPoint)
    self.assertEqual(testStatus, 10)
    #edges
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 1, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 9)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-2, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 9)
    #surface cells, non corner, not edge
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 30, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 8)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-30, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 8)
    #interior
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 421, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 1)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 20*20*10+35, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 1)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-422, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 1)

  def test_PhactoriFindSurfaceStatusForOneCellUsingFaceInfo_FlatCells(self):
    wv1 = paraview.simple.Wavelet()
    xfrm1 = Transform(Input=wv1)
    xfrm1.Transform = 'Transform'
    xfrm1.Transform.Scale = [1.0, 0.005, 1.0]
    xfrm1.UpdatePipeline()
    wv1grid = xfrm1.GetClientSideObject().GetOutputDataObject(0)
    cellsTouchingEachFace = PhactoriCountCellTouchingEachFace(wv1grid)
    self.assertEqual(cellsTouchingEachFace.totalNumFaces, 25200)
    #corner cells
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, 0, cellsTouchingEachFace)
    self.assertEqual(testStatus, 13)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, 19, cellsTouchingEachFace)
    self.assertEqual(testStatus, 13)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, (20*20*20-1), cellsTouchingEachFace)
    self.assertEqual(testStatus, 13)
    #edges
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, 1, cellsTouchingEachFace)
    self.assertEqual(testStatus, 12)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, 20, cellsTouchingEachFace)
    counts1 = cellsTouchingEachFace.GetNumberOfCellsTouchingEachFaceOfOneCell(20)
    extList1 = cellsTouchingEachFace.GetListOfExteriorFacesOnCellByIndex(20)
    self.assertEqual(testStatus, 6)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, (20*20*20-1)-1, cellsTouchingEachFace)
    self.assertEqual(testStatus, 12)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, (20*20*20-1)-20, cellsTouchingEachFace)
    self.assertEqual(testStatus, 6)
    #surface cells, non corner, not edge
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, 30, cellsTouchingEachFace)
    self.assertEqual(testStatus, 5)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, 402, cellsTouchingEachFace)
    self.assertEqual(testStatus, 11)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, 21, cellsTouchingEachFace)
    self.assertEqual(testStatus, 5)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, (20*20*20-1)-29, cellsTouchingEachFace)
    self.assertEqual(testStatus, 5)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, (20*20*20-1)-402, cellsTouchingEachFace)
    self.assertEqual(testStatus, 11)
    #interior
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, 421, cellsTouchingEachFace)
    self.assertEqual(testStatus, 1)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, 20*20*10+35, cellsTouchingEachFace)
    self.assertEqual(testStatus, 1)
    testStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(wv1grid, (20*20*20-1)-421, cellsTouchingEachFace)
    self.assertEqual(testStatus, 1)

  def test_PhactoriFindSurfaceStatusForOneCell_FlatCells(self):
    wv1 = paraview.simple.Wavelet()
    xfrm1 = Transform(Input=wv1)
    xfrm1.Transform = 'Transform'
    #xfrm1.Transform.Scale = [1.0, 0.05, 1.0]
    xfrm1.Transform.Scale = [1.0, 0.005, 1.0]
    xfrm1.UpdatePipeline()
    wv1grid = xfrm1.GetClientSideObject().GetOutputDataObject(0)
    cellsTouchingEachPoint = PhactoriCountCellTouchingEachPoint(wv1grid)
    #corner cells
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 0, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 13)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 19, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 13)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1), cellsTouchingEachPoint)
    self.assertEqual(testStatus, 13)
    #edges
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 1, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 12)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 20, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 6)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-1, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 12)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-20, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 6)
    #surface cells, non corner, not edge
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 30, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 5)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 402, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 11)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 21, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 5)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-29, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 5)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-402, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 11)
    #interior
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 421, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 1)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 20*20*10+35, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 1)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-421, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 1)

  def test_PhactoriFindSurfaceStatusForOneCell_SlightlyFlatCells(self):
    wv1 = paraview.simple.Wavelet()
    xfrm1 = Transform(Input=wv1)
    xfrm1.Transform = 'Transform'
    xfrm1.Transform.Scale = [1.0, 0.25, 1.0]
    xfrm1.UpdatePipeline()
    wv1grid = xfrm1.GetClientSideObject().GetOutputDataObject(0)
    cellsTouchingEachPoint = PhactoriCountCellTouchingEachPoint(wv1grid)
    #corner cells
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 0, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 10)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 19, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 10)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1), cellsTouchingEachPoint)
    self.assertEqual(testStatus, 10)
    #edges
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 1, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 9)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 20, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 3)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-1, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 9)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-20, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 3)
    #surface cells, non corner, not edge
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 30, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 2)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 402, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 8)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 21, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 2)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-29, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 2)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-402, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 8)
    #interior
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 421, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 1)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, 20*20*10+35, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 1)
    testStatus = PhactoriFindSurfaceStatusForOneCell(wv1grid, (20*20*20-1)-421, cellsTouchingEachPoint)
    self.assertEqual(testStatus, 1)

  def test_PhactoriFindFaceNormal(self):
    testugrid1 = self.MakeOneHexahdronGrid1()
    testhex1 = testugrid1.GetCell(0)
    testNormal = PhactoriFindFaceNormal(testugrid1, testhex1.GetFace(0))
    self.assertEqual(testNormal, [1.0, 0.0, 0.0])
    testNormal = PhactoriFindFaceNormal(testugrid1, testhex1.GetFace(1))
    self.assertEqual(testNormal, [-1.0, 0.0, 0.0])
    testNormal = PhactoriFindFaceNormal(testugrid1, testhex1.GetFace(2))
    self.assertEqual(testNormal, [0.0, 1.0, 0.0])
    testNormal = PhactoriFindFaceNormal(testugrid1, testhex1.GetFace(3))
    self.assertEqual(testNormal, [0.0, -1.0, 0.0])
    testNormal = PhactoriFindFaceNormal(testugrid1, testhex1.GetFace(4))
    self.assertEqual(testNormal, [0.0, 0.0, 1.0])
    testNormal = PhactoriFindFaceNormal(testugrid1, testhex1.GetFace(5))
    self.assertEqual(testNormal, [0.3076923076923077, 0.23076923076923078, -0.9230769230769231])

  def test_PhactoriFindLargestCellFaceNormal(self):
    testugrid1 = self.MakeOneHexahdronGrid1()
    testhex1 = testugrid1.GetCell(0)
    testNormal = PhactoriFindLargestCellFaceNormal(testugrid1, testhex1)
    self.assertEqual(testNormal, [0.3076923076923077, 0.23076923076923078, -0.9230769230769231])

  def test_PhactoriGetCellEdgeVector(self):
    testugrid1 = self.MakeOneHexahdronGrid1()
    testhex1 = testugrid1.GetCell(0)
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 0)
    self.assertEqual(testEdgeVec, [3.0, 0.0, 0.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 1)
    self.assertEqual(testEdgeVec, [0.0, 4.0, 0.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 2)
    self.assertEqual(testEdgeVec, [3.0, 0.0, 0.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 3)
    self.assertEqual(testEdgeVec, [0.0, 4.0, 0.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 4)
    self.assertEqual(testEdgeVec, [3.0, 0.0, 1.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 5)
    self.assertEqual(testEdgeVec, [0.0, 4.0, 1.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 6)
    self.assertEqual(testEdgeVec, [3.0, 0.0, -1.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 7)
    self.assertEqual(testEdgeVec, [0.0, 4.0, 3.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 8)
    self.assertEqual(testEdgeVec, [0.0, 0.0, 1.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 9)
    self.assertEqual(testEdgeVec, [0.0, 0.0, 2.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 10)
    self.assertEqual(testEdgeVec, [0.0, 0.0, 4.0])
    testEdgeVec = PhactoriGetCellEdgeVector(testugrid1, testhex1, 11)
    self.assertEqual(testEdgeVec, [0.0, 0.0, 3.0])
    
  def test_PhactoriFindCellEdgeAngleMetricsForOneCell(self):
    testugrid1 = self.MakeOneHexahdronGrid1()
    testAngle, testHeight = PhactoriFindCellEdgeAngleMetricsForOneCell(testugrid1, 0, True, 0)
    self.assertEqual(testAngle, 67.38013505195958)
    self.assertEqual(testHeight, 3.692307692307693) 
    testugrid2 = self.MakeOneFlatHexahdronGrid1()
    testAngle, testHeight = PhactoriFindCellEdgeAngleMetricsForOneCell(testugrid2, 0, True, 0)
    self.assertEqual(testAngle, 10.024987862075733)
    self.assertEqual(testHeight, 0.03125)
    testugrid3 = self.MakeOneFlatHexahdronGrid2()
    testAngle, testHeight = PhactoriFindCellEdgeAngleMetricsForOneCell(testugrid3, 0, True, 0)
    self.assertEqual(testAngle, 90.0)
    self.assertEqual(testHeight, 0.03125)

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


