# Copyright(C) 1999-2020, 2024, 2024 National Technology & Engineering Solutions
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
from Operation.PhactoriSampledCellInfo import *
from Operation.PhactoriSegmentCellSampler3 import *
from paraview.simple import *
import os

class TestPhactoriSegmentCellSampler3(unittest.TestCase):

  def test_ValidateJsonStructuredSeedCellList(self):
    pscs1 = PhactoriSegmentCellSampler3()
    goodSeedPointsJson = [
      {
        "geometric seed point":[0.05, 0.00001, 0.0],
        "collection axis":"i"
      },
      {
        "geometric seed point":[0.10, 0.00001, 0.0],
        "collection axis":"j"
      },
      {
        "geometric seed point":[0.20, 0.00001, 0.0]
      },
      {
        "geometric seed point":[0.50, 0.00001, 0.0],
        "collection axis":"k"
      },
      {
        "geometric seed point":[0.30, 0.00001, 0.0],
        "collection axis":"ij"
      },
      {
        "geometric seed point":[0.40, 0.00001, 0.0],
        "collection axis":"ik"
      },
      {
        "geometric seed point":[0.45, 0.00001, 0.0],
        "collection axis":"jk"
      }
    ]
    self.assertTrue(pscs1.ValidateJsonStructuredSeedCellList(goodSeedPointsJson))

    badSeedPointsJson1 = []
    with self.assertRaises(Exception):
      pscs1.ValidateJsonSegmentList(badSeedPointsJson1)

    badSeedPointsJson2 = goodSeedPointsJson
    badSeedPointsJson2 = [
      {
        "X geometric seed point":[0.05, 0.00001, 0.0],
        "collection axis":"i"
      }
    ]
    with self.assertRaises(Exception):
      pscs1.ValidateJsonSegmentList(badSeedPointsJson2)

    badSeedPointsJson3 = goodSeedPointsJson
    badSeedPointsJson3 = [
      {
        "geometric seed point":[0.00001, 0.0],
        "collection axis":"i"
      }
    ]
    with self.assertRaises(Exception):
      pscs1.ValidateJsonSegmentList(badSeedPointsJson3)

    badSeedPointsJson4 = goodSeedPointsJson
    badSeedPointsJson4 = [
      {
        "geometric seed point":[0.00001, 0.0, 1.7],
        "collection axis":"m"
      }
    ]
    with self.assertRaises(Exception):
      pscs1.ValidateJsonSegmentList(badSeedPointsJson4)

  def test_ValidateJsonSegmentList_format_two_point(self):
    pscs1 = PhactoriSegmentCellSampler3()

    goodJsonSegments = [[[0.0, 0.0, 0.0],[1.0,1.0,1.0]],[[2.0,3.0,4.0],[5.0,6.0,7.0]],[[-1.0,-2.0,17.0],[3.4,-16.7,143.9]]]
    self.assertTrue(pscs1.ValidateJsonSegmentList(goodJsonSegments))

    badJsonSegments1 = [[[0.0, 0.0, 0.0],[1.0,1.0,1.0]],[[2.0,3.0],[5.0,6.0,7.0]]]
    with self.assertRaises(Exception):
      pscs1.ValidateJsonSegmentList(badJsonSegments1)
    badJsonSegments2 = []
    with self.assertRaises(Exception):
      pscs1.ValidateJsonSegmentList(badJsonSegments2)
    badJsonSegments3 = [[[0.0, 0.0, 0.0],[1.0,1.0,1.0]],[2.0,3.0,4.0],[5.0,6.0,7.0]]
    with self.assertRaises(Exception):
      pscs1.ValidateJsonSegmentList(badJsonSegments3)

  def test_ValidateJsonSegmentList_format_seed_point_and_length(self):
    pscs1 = PhactoriSegmentCellSampler3()
    goodJsonSegments1 = [[[0.0, 1.0, 2.0], 37.6],[[-67.4,32.1,5.6], 14.0]]
    self.assertTrue(pscs1.ValidateJsonSegmentList2(goodJsonSegments1))
    badJsonSegments2 = [[[0.0, 1.0, 2.0]],[[-67.4,32.1,5.6], 14.0]]
    with self.assertRaises(Exception):
      pscs1.ValidateJsonSegmentList2(badJsonSegments2)
    badJsonSegments3 = [[[0.0, 1.0, 2.0], 37.6],[[-67.4,5.6], 14.0]]
    with self.assertRaises(Exception):
      pscs1.ValidateJsonSegmentList2(badJsonSegments3)

  def test_GetListOfCellTestPointsNearestListOfPoints(self):
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"
    operationParams = {
      "filename":"dummyfilename",
      "cell center to segment test distance":0.5,
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'segmentcellsampler3',
              PhactoriSegmentCellSampler3,
              operationParams)

    testPoints = [[0.1, 0.2, 0.3], [1.1, 1.2, 1.3], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]
    testSphere = Sphere()
    testSphere.UpdatePipeline()
    #testSphere.Update()
    #pvClientSideData = testSphere.GetClientSideObject().GetOutputDataObject(0)
    #if pvClientSideData == None:
        #myDebugPrint3('DoMethodPerBlock: pvClientSideData is None\n')
    #outstr3 = "num points in sphere: " + str(pvClientSideData.GetNumberOfPoints()) + "\n"
    #myDebugPrint3(outstr3)

    PhactoriSegmentCellSampler3Instance = newOperationBlock.mOperationSpecifics
    PhactoriSegmentCellSampler3Instance.myCopyOfInputFilter = testSphere
    thisProcessNearestCellPointList, thisProcDistSqrdList = PhactoriSegmentCellSampler3Instance.\
      GetCellsClosestToPointsOnThisProcessFromParaViewFilter(testPoints)
    #outstr1 = "thisProcessNearestCellPointList:\n" + str(thisProcessNearestCellPointList) + "\n"
    #outstr2 = "thisProcDistSqrdList:\n" + str(thisProcDistSqrdList) + "\n"
    #myDebugPrint3(outstr1)
    #myDebugPrint3(outstr2)

    goldNearestCellPointList = [
        [0.07670053094625473, 0.2721584066748619, 0.381114661693573],
        [0.30392880737781525, 0.13820958137512207, 0.381114661693573],
        [0.33366745710372925, 0.17234453558921814, 0.21150268241763115],
        [0.41607651114463806, 0.17234453558921814, 0.0]]
    goldProcDistSqrdList = [0.012329289253703989, 2.6054785017610183, 0.1857699955473856, 0.03674579092911934]

    self.assertEqual(thisProcessNearestCellPointList, goldNearestCellPointList)
    self.assertEqual(thisProcDistSqrdList, goldProcDistSqrdList)

    #myDebugPrint3("starting box now 2:")
    testBox = Box()
    testBox.UpdatePipeline()
    testPoints = [[0.45, 0.0, 0.0],[0.0, -0.55, 0.0],[0.0, 0.0, 0.2]]

    PhactoriSegmentCellSampler3Instance.myCopyOfInputFilter = testBox
    thisProcessNearestCellPointList, thisProcDistSqrdList = PhactoriSegmentCellSampler3Instance.\
      GetCellsClosestToPointsOnThisProcessFromParaViewFilter(testPoints)
    outstr1 = "box thisProcessNearestCellPointList:\n" + str(thisProcessNearestCellPointList) + "\n"
    outstr2 = "box thisProcDistSqrdList:\n" + str(thisProcDistSqrdList) + "\n"
    #myDebugPrint3(outstr1)
    #myDebugPrint3(outstr2)
    #myDebugPrint3("done box now:")
    goldBoxNearestCellPointList = [[0.5, 0.0, 0.0], [0.0, -0.5, 0.0], [0.0, 0.0, 0.5]]
    goldBoxProcDistSqrdList = [0.0024999999999999988, 0.0025000000000000044, 0.09]

    self.assertEqual(thisProcessNearestCellPointList, goldBoxNearestCellPointList)
    self.assertEqual(thisProcDistSqrdList, goldBoxProcDistSqrdList)

  def test_GetListOfCellTestPointsNearestListOfPointsStructured(self):
    testWavelet = Wavelet()
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"
    operationParams = {
      "filename":"dummyfilename",
      "cell center to segment test distance":0.5,
      "vtk grid type":"structured",
      "segment definition format":"geometric point and nearest cell"
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'segmentcellsampler3',
              PhactoriSegmentCellSampler3,
              operationParams)

    testGeometryPoints = [
      [-11.0,  -0.1,  -0.1],
      [ -0.1, -11.0,  -0.1],
      [ -0.1,  -0.1, -11.0],
      [-11.0, -11.0,  -0.1],
      [-11.0,  -0.1, -11.0],
      [ -0.1, -11.0, -11.0],
      [-11.0, -11.0, -11.0],
      [ 11.0,   0.1,   0.1],
      [  0.1,  11.0,   0.1],
      [  0.1,   0.1,  11.0],
      [ 11.0,  11.0,   0.1],
      [ 11.0,   0.1,  11.0],
      [  0.1,  11.0,  11.0],
      [ 11.0,  11.0,  11.0],
      [  0.1,   0.1,   0.1],
      [ -0.1,  -0.1,  -0.1]
      ]

    newOperationBlock.mOperationSpecifics.myCopyOfInputFilter = testWavelet
    testWavelet.UpdatePipeline()
    globalNearstPointlist = newOperationBlock.mOperationSpecifics.\
      GetListOfCellTestPointsNearestListOfPointsStructured(testGeometryPoints)

    goldNearestPointList = [
      [ -9.5,  -0.5,  -0.5],
      [ -0.5,  -9.5,  -0.5],
      [ -0.5,  -0.5,  -9.5],
      [ -9.5,  -9.5,  -0.5],
      [ -9.5,  -0.5,  -9.5],
      [ -0.5,  -9.5,  -9.5],
      [ -9.5,  -9.5,  -9.5],
      [  9.5,   0.5,   0.5],
      [  0.5,   9.5,   0.5],
      [  0.5,   0.5,   9.5],
      [  9.5,   9.5,   0.5],
      [  9.5,   0.5,   9.5],
      [  0.5,   9.5,   9.5],
      [  9.5,   9.5,   9.5],
      [  0.5,   0.5,   0.5],
      [ -0.5,  -0.5,  -0.5]
    ]
    goldIjkList = [
      [-10, -1, -1],
      [-1, -10, -1],
      [-1, -1, -10],
      [-10, -10, -1],
      [-10, -1, -10],
      [-1, -10, -10],
      [-10, -10, -10],
      [9, 0, 0],
      [0, 9, 0],
      [0, 0, 9],
      [9, 9, 0],
      [9, 0, 9],
      [0, 9, 9],
      [9, 9, 9],
      [0, 0, 0],
      [-1, -1, -1]
    ]
    outstr1 = "globalNearstPointlist begin:\n"
    #myDebugPrint3(outstr1)
    #for ii, oneCell in enumerate(globalNearstPointlist):
    #  myDebugPrint3(str(oneCell.ijk) + ",\n")
    for ii, oneCell in enumerate(globalNearstPointlist):
      #myDebugPrint3(oneCell.ToStr())
      self.assertEqual(oneCell.cellTestPoint, goldNearestPointList[ii])
      self.assertEqual(oneCell.ijk, goldIjkList[ii])
    outstr1 = "globalNearstPointlist end\n"
    #myDebugPrint3(outstr1)

  def test_GatherStructuredCellsFromSeedCells(self):
    testWavelet = Wavelet()
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"
    operationParams = {
      "filename":"dummyfilename",
      "cell center to segment test distance":0.5,
      "vtk grid type":"structured",
      "collection method":"seed cell with structured grid"
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'segmentcellsampler3',
              PhactoriSegmentCellSampler3,
              operationParams)

    PhactoriSegmentCellSampler3Instance = newOperationBlock.mOperationSpecifics
    PhactoriSegmentCellSampler3Instance.myCopyOfInputFilter = testWavelet

    #default collection axis is "k"
    testGeometryPointsJson = [
      {"geometric seed point":[0.5,  0.5,  0.5]},
      {"geometric seed point":[5.0,  5.0,  -5.0]},
    ]

    self.assertTrue(PhactoriSegmentCellSampler3Instance.ValidateJsonStructuredSeedCellList(testGeometryPointsJson))
    PhactoriSegmentCellSampler3Instance.CreateInternalStructuredSeedCellListFromJson(testGeometryPointsJson)
    testWavelet.UpdatePipeline()

    perSegmentStructuredCellList = PhactoriSegmentCellSampler3Instance.GatherStructuredCellsFromSeedCells()
    self.assertEqual(len(perSegmentStructuredCellList),len(PhactoriSegmentCellSampler3Instance.StructuredNearbyCellPointList))

    goldPoints = [[],[]]
    for ii in range(0,20):
      goldPoints[0].append([0.5, 0.5, -9.5 + float(ii)])
      goldPoints[1].append([4.5, 4.5, -9.5 + float(ii)])

    for ii, oneSeglist in enumerate(perSegmentStructuredCellList):
      localResultPoints = []
      localDataTuples = []
      #myDebugPrint3("begin cell test points for segment k axis " + str(ii) + "\n")
      for oneStructuredCell in oneSeglist:
        localResultPoints.append(oneStructuredCell.cellTestPoint)
        #myDebugPrint3(str(oneStructuredCell.cellTestPoint) + "\n")
        localDataTuples.append(oneStructuredCell.dataTuple)
      self.assertEqual(localResultPoints, goldPoints[ii])
      #myDebugPrint3("end cell test points for segment k axis " + str(ii) + "\n")

    #test j axis collection instead of k
    operationParams = {
      "filename":"dummyfilename",
      "cell center to segment test distance":0.5,
      "vtk grid type":"structured",
      "collection method":"seed cell with structured grid"
    }

    PhactoriSegmentCellSampler3Instance.ParseParametersFromJson(operationParams)

    testGeometryPointsJson = [
      {"geometric seed point":[0.5,  0.5,  0.5],
       "collection axis":"j"},
      {"geometric seed point":[5.0,  5.0,  -5.0],
       "collection axis":"j"}
    ]

    self.assertTrue(PhactoriSegmentCellSampler3Instance.ValidateJsonStructuredSeedCellList(testGeometryPointsJson))
    PhactoriSegmentCellSampler3Instance.CreateInternalStructuredSeedCellListFromJson(testGeometryPointsJson)
    testWavelet.UpdatePipeline()

    perSegmentStructuredCellList = PhactoriSegmentCellSampler3Instance.GatherStructuredCellsFromSeedCells()
    self.assertEqual(len(perSegmentStructuredCellList),len(PhactoriSegmentCellSampler3Instance.StructuredNearbyCellPointList))

    goldPoints = [[],[]]
    for ii in range(0,20):
      goldPoints[0].append([0.5, -9.5 + float(ii), 0.5])
      goldPoints[1].append([4.5, -9.5 + float(ii), -5.5])

    for ii, oneSeglist in enumerate(perSegmentStructuredCellList):
      localResultPoints = []
      #myDebugPrint3("begin cell test points for segment j axis " + str(ii) + "\n")
      for oneStructuredCell in oneSeglist:
        localResultPoints.append(oneStructuredCell.cellTestPoint)
        #myDebugPrint3(str(oneStructuredCell.cellTestPoint) + "\n")
      self.assertEqual(localResultPoints, goldPoints[ii])
      #myDebugPrint3("end cell test points for segment j axis " + str(ii) + "\n")

    #test i axis collection instead of j or k
    operationParams = {
      "filename":"dummyfilename",
      "cell center to segment test distance":0.5,
      "vtk grid type":"structured",
      "collection method":"seed cell with structured grid"
    }

    PhactoriSegmentCellSampler3Instance.ParseParametersFromJson(operationParams)

    testGeometryPointsJson = [
      {"geometric seed point":[0.5,  0.5,  0.5],
       "collection axis":"i"},
      {"geometric seed point":[5.0,  5.0,  -5.0],
       "collection axis":"i"}
    ]

    self.assertTrue(PhactoriSegmentCellSampler3Instance.ValidateJsonStructuredSeedCellList(testGeometryPointsJson))
    PhactoriSegmentCellSampler3Instance.CreateInternalStructuredSeedCellListFromJson(testGeometryPointsJson)
    testWavelet.UpdatePipeline()

    perSegmentStructuredCellList = PhactoriSegmentCellSampler3Instance.GatherStructuredCellsFromSeedCells()
    self.assertEqual(len(perSegmentStructuredCellList),len(PhactoriSegmentCellSampler3Instance.StructuredNearbyCellPointList))

    goldPoints = [[],[]]
    for ii in range(0,20):
      goldPoints[0].append([-9.5 + float(ii), 0.5, 0.5])
      goldPoints[1].append([-9.5 + float(ii), 4.5, -5.5])

    for ii, oneSeglist in enumerate(perSegmentStructuredCellList):
      localResultPoints = []
      #myDebugPrint3("begin cell test points for segment i axis " + str(ii) + "\n")
      for oneStructuredCell in oneSeglist:
        localResultPoints.append(oneStructuredCell.cellTestPoint)
        #myDebugPrint3(str(oneStructuredCell.cellTestPoint) + "\n")
      self.assertEqual(localResultPoints, goldPoints[ii])
      #myDebugPrint3("end cell test points for segment i axis " + str(ii) + "\n")

  def test_WriteAllDataFromOneProcessUsingMPI(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"

    testOutFileBasename = "test_WriteAllDataFromOneProcessUsingMPI_output_cells_"
    operationParams = {
      "filename":"dummyfilename",
      "cell center to segment test distance":0.5,
      "segment definition format":"two geometric points",
      "cell data array name":"RTData",
      "cell data array tuple size":1,
      "sampled cells output filename":testOutFileBasename,
      "sampled cells output directory":".",
      "sampled cells output filename extension":".txt",
      "collection method":"segments identify cells"
    }

    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'segmentcellsampler3',
              PhactoriSegmentCellSampler3,
              operationParams)

    PhactoriSegmentCellSampler3Instance = newOperationBlock.mOperationSpecifics
    PhactoriSegmentCellSampler3Instance.myCopyOfInputFilter = testWavelet

    testSegments = [
      [[0.5, 0.5, 0.5], [11.5, 11.5, 11.5]],
      [[5.0, 5.0, 5.0], [-8.0, -7.0, -6.0]]
    ]
    PhactoriSegmentCellSampler3Instance.ValidateJsonSegmentList(testSegments)
    PhactoriSegmentCellSampler3Instance.CreateInternalSegmentListFromJson(testSegments)

    perSegmentCellList = PhactoriSegmentCellSampler3Instance.SetMaskValueForCellsNearSegments(
      PhactoriSegmentCellSampler3Instance.myCopyOfInputFilter, PhactoriSegmentCellSampler3Instance.SegmentList, "mask1", "Ids",
      PhactoriSegmentCellSampler3Instance.MaskTestDistanceSquared, PhactoriSegmentCellSampler3Instance.ProjectionAxis)
    PhactoriSegmentCellSampler3Instance.PerSegmentMarkedCellInfoList = perSegmentCellList

    PhactoriSegmentCellSampler3Instance.WriteAllDataFromOneProcessUsingMPI(perSegmentCellList, "PhactoriSegmentCellSampler3_", 0)

    goldTestStringJson = """
    {
    "format for cells in lists": {"PhactoriSampledCellInfo output format 1 info":[
    " [cellTestPoint, ijk, dataTuple, pid, index, segmentIndex, collectionAxis]",
    " cellTestPoint is [X, Y, Z], ijk is [i, j, k], dataTuple is [c1, c2, ... cN]"]},
    "list of collected cells for each segment": [
    {
    "cells for segment index": 0,
    "segment point A": [0.5, 0.5, 0.5],
    "segment point B": [11.5, 11.5, 11.5],
    "geometry point used to find point A": [0.5, 0.5, 0.5],
    "number of cells collected for this segment": 10,
    "cell list for this segment": [
    [[0.5, 0.5, 0.5], [-1, -1, -1], [244.8736114501953], 0, 1, -1, 0, 2],
    [[1.5, 1.5, 1.5], [-1, -1, -1], [233.40476989746094], 0, 1, -1, 0, 2],
    [[2.5, 2.5, 2.5], [-1, -1, -1], [238.9958038330078], 0, 1, -1, 0, 2],
    [[3.5, 3.5, 3.5], [-1, -1, -1], [225.47100830078125], 0, 1, -1, 0, 2],
    [[4.5, 4.5, 4.5], [-1, -1, -1], [178.6822967529297], 0, 1, -1, 0, 2],
    [[5.5, 5.5, 5.5], [-1, -1, -1], [149.9278106689453], 0, 1, -1, 0, 2],
    [[6.5, 6.5, 6.5], [-1, -1, -1], [141.2411651611328], 0, 1, -1, 0, 2],
    [[7.5, 7.5, 7.5], [-1, -1, -1], [120.53302001953125], 0, 1, -1, 0, 2],
    [[8.5, 8.5, 8.5], [-1, -1, -1], [82.78067779541016], 0, 1, -1, 0, 2],
    [[9.5, 9.5, 9.5], [-1, -1, -1], [55.57729721069336], 0, 1, -1, 0, 2]
    ]
    },
    {
    "cells for segment index": 1,
    "segment point A": [5.0, 5.0, 5.0],
    "segment point B": [-8.0, -7.0, -6.0],
    "geometry point used to find point A": [5.0, 5.0, 5.0],
    "number of cells collected for this segment": 18,
    "cell list for this segment": [
    [[-7.5, -6.5, -5.5], [-1, -1, -1], [129.1289520263672], 0, 1, -1, 1, 2],
    [[-6.5, -5.5, -4.5], [-1, -1, -1], [170.07781982421875], 0, 1, -1, 1, 2],
    [[-5.5, -4.5, -3.5], [-1, -1, -1], [193.2375946044922], 0, 1, -1, 1, 2],
    [[-4.5, -3.5, -2.5], [-1, -1, -1], [199.42271423339844], 0, 1, -1, 1, 2],
    [[-3.5, -3.5, -2.5], [-1, -1, -1], [206.755859375], 0, 1, -1, 1, 2],
    [[-3.5, -2.5, -2.5], [-1, -1, -1], [217.16062927246094], 0, 1, -1, 1, 2],
    [[-2.5, -2.5, -1.5], [-1, -1, -1], [226.48875427246094], 0, 1, -1, 1, 2],
    [[-2.5, -1.5, -1.5], [-1, -1, -1], [249.02000427246094], 0, 1, -1, 1, 2],
    [[-1.5, -1.5, -0.5], [-1, -1, -1], [259.11346435546875], 0, 1, -1, 1, 2],
    [[-1.5, -0.5, -0.5], [-1, -1, -1], [260.3305358886719], 0, 1, -1, 1, 2],
    [[-0.5, -0.5, 0.5], [-1, -1, -1], [264.2397155761719], 0, 1, -1, 1, 2],
    [[-0.5, 0.5, 0.5], [-1, -1, -1], [246.28480529785156], 0, 1, -1, 1, 2],
    [[0.5, 0.5, 1.5], [-1, -1, -1], [238.22740173339844], 0, 1, -1, 1, 2],
    [[0.5, 1.5, 1.5], [-1, -1, -1], [234.4701690673828], 0, 1, -1, 1, 2],
    [[1.5, 1.5, 1.5], [-1, -1, -1], [233.40476989746094], 0, 1, -1, 1, 2],
    [[2.5, 2.5, 2.5], [-1, -1, -1], [238.9958038330078], 0, 1, -1, 1, 2],
    [[3.5, 3.5, 3.5], [-1, -1, -1], [225.47100830078125], 0, 1, -1, 1, 2],
    [[4.5, 4.5, 4.5], [-1, -1, -1], [178.6822967529297], 0, 1, -1, 1, 2]
    ]
    }
    ]
    }
    """

    goldJson = json.loads(goldTestStringJson)
    testFileName = testOutFileBasename + "00000.txt"
    ff = open(testFileName, "r")
    testJson = json.load(ff)
    ff.close()
    self.assertEqual(goldJson, testJson)

    #remove file that got made during test
    os.remove(testFileName)

  def test_GatherStructuredCellsFromSeedCells_with_ij_ik_jk(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"

    testOutFileBasename = "test_GatherStructuredCellsFromSeedCells_with_ij_ik_jk_output_cells_"
    operationParams = {
      "filename":"dummyfilename",
      "cell center to segment test distance":0.5,
      "cell data array name":"RTData",
      "cell data array tuple size":1,
      "sampled cells output filename":testOutFileBasename,
      "sampled cells output directory":".",
      "collection method":"seed cell with structured grid"
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'segmentcellsampler3',
              PhactoriSegmentCellSampler3,
              operationParams)

    PhactoriSegmentCellSampler3Instance = newOperationBlock.mOperationSpecifics
    PhactoriSegmentCellSampler3Instance.myCopyOfInputFilter = testWavelet

    #default collection axis is "k"
    testGeometryPointsJson = [
      {"geometric seed point":[0.5,  0.5,  0.5],
       "collection axis":"ij"},
      {"geometric seed point":[-0.125,  -0.125,  -0.125],
       "collection axis":"ik"},
      {"geometric seed point":[4.75,  3.75,  -2.125],
       "collection axis":"jk"},
      {"geometric seed point":[0.5,  0.5,  11.0],
       "collection axis":"ij"},
      {"geometric seed point":[0.5,  0.5,  -11.0],
       "collection axis":"ij"},
      {"geometric seed point":[11.0, 11.0, 11.0],
       "collection axis":"ik"},
      {"geometric seed point":[-11.0, -11.0, -11.0],
       "collection axis":"jk"}
    ]
    self.assertTrue(PhactoriSegmentCellSampler3Instance.ValidateJsonStructuredSeedCellList(testGeometryPointsJson))
    PhactoriSegmentCellSampler3Instance.CreateInternalStructuredSeedCellListFromJson(testGeometryPointsJson)
    testWavelet.UpdatePipeline()

    perSegmentStructuredCellList = PhactoriSegmentCellSampler3Instance.GatherStructuredCellsFromSeedCells()
    #PhactoriSegmentCellSampler3Instance.WriteAllDataFromOneProcessUsingMPIStructured(perSegmentStructuredCellList, 0)

    for oneSeedPointCellList in perSegmentStructuredCellList:
      oneSeedPointCellList[0].index = -1
      oneSeedPointCellList[-1].index = -1

    goldStr = "[[-9.5, -9.5, 0.5], [-10, -10, 0], [117.52127838134766], 0, 1, -1, 0, 2]"
    testStr = perSegmentStructuredCellList[0][0].ToStrTerseOneLineList()
    self.assertEqual(goldStr, testStr)
    goldStr = "[[9.5, 9.5, 0.5], [9, 9, 0], [91.66453552246094], 0, 1, -1, 0, 2]"
    testStr = perSegmentStructuredCellList[0][-1].ToStrTerseOneLineList()
    self.assertEqual(goldStr, testStr)
    testStr = perSegmentStructuredCellList[1][0].ToStrTerseOneLineList()
    goldStr = "[[-9.5, -0.5, -9.5], [-10, -1, -10], [114.62344360351562], 0, 1, -1, 1, 2]"
    testStr = perSegmentStructuredCellList[1][-1].ToStrTerseOneLineList()
    goldStr = "[[9.5, -0.5, 9.5], [9, -1, 9], [114.94000244140625], 0, 1, -1, 1, 2]"
    testStr = perSegmentStructuredCellList[2][0].ToStrTerseOneLineList()
    goldStr = "[[4.5, -9.5, -9.5], [4, -10, -10], [108.53125], 0, 1, -1, 2, 2]"
    testStr = perSegmentStructuredCellList[2][-1].ToStrTerseOneLineList()
    goldStr = "[[4.5, 9.5, 9.5], [4, 9, 9], [82.35794830322266], 0, 1, -1, 2, 2]"
    testStr = perSegmentStructuredCellList[3][0].ToStrTerseOneLineList()
    goldStr = "[[-9.5, -9.5, 9.5], [-10, -10, 9], [81.43404388427734], 0, 1, -1, 3, 2]"
    testStr = perSegmentStructuredCellList[3][-1].ToStrTerseOneLineList()
    goldStr = "[[9.5, 9.5, 9.5], [9, 9, 9], [55.57729721069336], 0, 1, -1, 3, 2]"
    testStr = perSegmentStructuredCellList[4][0].ToStrTerseOneLineList()
    goldStr = "[[-9.5, -9.5, -9.5], [-10, -10, -10], [81.43404388427734], 0, 1, -1, 4, 2]"
    testStr = perSegmentStructuredCellList[4][-1].ToStrTerseOneLineList()
    goldStr = "[[9.5, 9.5, -9.5], [9, 9, -10], [55.57729721069336], 0, 1, -1, 4, 2]"
    testStr = perSegmentStructuredCellList[5][0].ToStrTerseOneLineList()
    goldStr = "[[-9.5, 9.5, -9.5], [-10, 9, -10], [55.2607421875], 0, 1, -1, 5, 2]"
    testStr = perSegmentStructuredCellList[5][-1].ToStrTerseOneLineList()
    goldStr = "[[9.5, 9.5, 9.5], [9, 9, 9], [55.57729721069336], 0, 1, -1, 5, 2]"
    testStr = perSegmentStructuredCellList[6][0].ToStrTerseOneLineList()
    goldStr = "[[-9.5, -9.5, -9.5], [-10, -10, -10], [81.43404388427734], 0, 1, -1, 6, 2]"
    testStr = perSegmentStructuredCellList[6][-1].ToStrTerseOneLineList()
    goldStr = "[[-9.5, 9.5, 9.5], [-10, 9, 9], [55.2607421875], 0, 1, -1, 6, 2]"


  def test_WriteAllDataFromOneProcessUsingMPIStructured(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"

    testOutFileBasename = "test_WriteAllDataFromOneProcessUsingMPIStructured_output_cells_"
    operationParams = {
      "filename":"dummyfilename",
      "cell center to segment test distance":0.5,
      "cell data array name":"RTData",
      "cell data array tuple size":1,
      "sampled cells output filename":testOutFileBasename,
      "sampled cells output directory":".",
      "collection method":"seed cell with structured grid"
    }
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'segmentcellsampler3',
              PhactoriSegmentCellSampler3,
              operationParams)

    PhactoriSegmentCellSampler3Instance = newOperationBlock.mOperationSpecifics
    PhactoriSegmentCellSampler3Instance.myCopyOfInputFilter = testWavelet

    #default collection axis is "k"
    testGeometryPointsJson = [
      {"geometric seed point":[0.5,  0.5,  0.5],
       "collection axis":"i"},
      {"geometric seed point":[5.0,  5.0,  -5.0],
       "collection axis":"i"},
      {"geometric seed point":[0.5,  0.5,  0.5],
       "collection axis":"j"},
      {"geometric seed point":[5.0,  5.0,  -5.0],
       "collection axis":"j"},
      {"geometric seed point":[0.5,  0.5,  0.5],
       "collection axis":"k"},
      {"geometric seed point":[5.0,  5.0,  -5.0],
       "collection axis":"k"}
    ]

    self.assertTrue(PhactoriSegmentCellSampler3Instance.ValidateJsonStructuredSeedCellList(testGeometryPointsJson))
    PhactoriSegmentCellSampler3Instance.CreateInternalStructuredSeedCellListFromJson(testGeometryPointsJson)
    testWavelet.UpdatePipeline()

    perSegmentStructuredCellList = PhactoriSegmentCellSampler3Instance.GatherStructuredCellsFromSeedCells()
    PhactoriSegmentCellSampler3Instance.WriteAllDataFromOneProcessUsingMPIStructured(perSegmentStructuredCellList, 0)

    goldTestStringJson = """
    {
    "format for cells in lists": {"PhactoriSampledCellInfo output format 1 info":[
    " [cellTestPoint, ijk, dataTuple, pid, index, segmentIndex, collectionAxis]",
    " cellTestPoint is [X, Y, Z], ijk is [i, j, k], dataTuple is [c1, c2, ... cN]"]},
    "list of collected cells for each seed cell": [
    {
    "cells for seed cell index": 0,
    "structured seed cell xyz": [0.5, 0.5, 0.5],
    "structured seed cell ijk": [0, 0, 0],
    "structured seed cell complete": [[0.5, 0.5, 0.5], [0, 0, 0], [], 0, 1, -1, -1, 0],
    "geometry point used to find seed cell": [0.5, 0.5, 0.5],
    "number of cells collected for this seed cell": 20,
    "cell list for this seed cell": [
    [[-9.5, 0.5, 0.5], [-10, 0, 0], [153.8872833251953], 0, 1, -1, 0, 2],
    [[-8.5, 0.5, 0.5], [-9, 0, 0], [169.47491455078125], 0, 1, -1, 0, 2],
    [[-7.5, 0.5, 0.5], [-8, 0, 0], [183.55810546875], 0, 1, -1, 0, 2],
    [[-6.5, 0.5, 0.5], [-7, 0, 0], [198.17469787597656], 0, 1, -1, 0, 2],
    [[-5.5, 0.5, 0.5], [-6, 0, 0], [209.9029998779297], 0, 1, -1, 0, 2],
    [[-4.5, 0.5, 0.5], [-5, 0, 0], [222.11965942382812], 0, 1, -1, 0, 2],
    [[-3.5, 0.5, 0.5], [-4, 0, 0], [230.2534637451172], 0, 1, -1, 0, 2],
    [[-2.5, 0.5, 0.5], [-3, 0, 0], [238.78099060058594], 0, 1, -1, 0, 2],
    [[-1.5, 0.5, 0.5], [-2, 0, 0], [242.37562561035156], 0, 1, -1, 0, 2],
    [[-0.5, 0.5, 0.5], [-1, 0, 0], [246.28480529785156], 0, 1, -1, 0, 2],
    [[0.5, 0.5, 0.5], [0, 0, 0], [244.8736114501953], 0, 1, -1, 0, 2],
    [[1.5, 0.5, 0.5], [1, 0, 0], [243.7585906982422], 0, 1, -1, 0, 2],
    [[2.5, 0.5, 0.5], [2, 0, 0], [237.45396423339844], 0, 1, -1, 0, 2],
    [[3.5, 0.5, 0.5], [3, 0, 0], [231.49801635742188], 0, 1, -1, 0, 2],
    [[4.5, 0.5, 0.5], [4, 0, 0], [220.98251342773438], 0, 1, -1, 0, 2],
    [[5.5, 0.5, 0.5], [5, 0, 0], [210.91000366210938], 0, 1, -1, 0, 2],
    [[6.5, 0.5, 0.5], [6, 0, 0], [197.31800842285156], 0, 1, -1, 0, 2],
    [[7.5, 0.5, 0.5], [7, 0, 0], [184.24732971191406], 0, 1, -1, 0, 2],
    [[8.5, 0.5, 0.5], [8, 0, 0], [168.96693420410156], 0, 1, -1, 0, 2],
    [[9.5, 0.5, 0.5], [9, 0, 0], [154.20382690429688], 0, 1, -1, 0, 2]
    ]
    },
    {
    "cells for seed cell index": 1,
    "structured seed cell xyz": [4.5, 4.5, -5.5],
    "structured seed cell ijk": [4, 4, -6],
    "structured seed cell complete": [[4.5, 4.5, -5.5], [4, 4, -6], [], 0, 1, -1, -1, 0],
    "geometry point used to find seed cell": [5.0, 5.0, -5.0],
    "number of cells collected for this seed cell": 20,
    "cell list for this seed cell": [
    [[-9.5, 4.5, -5.5], [-10, 4, -6], [119.83075714111328], 0, 1, -1, 1, 2],
    [[-8.5, 4.5, -5.5], [-9, 4, -6], [132.06900024414062], 0, 1, -1, 1, 2],
    [[-7.5, 4.5, -5.5], [-8, 4, -6], [142.9117431640625], 0, 1, -1, 1, 2],
    [[-6.5, 4.5, -5.5], [-7, 4, -6], [154.47286987304688], 0, 1, -1, 1, 2],
    [[-5.5, 4.5, -5.5], [-6, 4, -6], [163.4069366455078], 0, 1, -1, 1, 2],
    [[-4.5, 4.5, -5.5], [-5, 4, -6], [173.16384887695312], 0, 1, -1, 1, 2],
    [[-3.5, 4.5, -5.5], [-4, 4, -6], [179.2395782470703], 0, 1, -1, 1, 2],
    [[-2.5, 4.5, -5.5], [-3, 4, -6], [186.16876220703125], 0, 1, -1, 1, 2],
    [[-1.5, 4.5, -5.5], [-2, 4, -6], [188.67095947265625], 0, 1, -1, 1, 2],
    [[-0.5, 4.5, -5.5], [-1, 4, -6], [192.0256805419922], 0, 1, -1, 1, 2],
    [[0.5, 4.5, -5.5], [0, 4, -6], [190.61447143554688], 0, 1, -1, 1, 2],
    [[1.5, 4.5, -5.5], [1, 4, -6], [190.0539093017578], 0, 1, -1, 1, 2],
    [[2.5, 4.5, -5.5], [2, 4, -6], [184.84173583984375], 0, 1, -1, 1, 2],
    [[3.5, 4.5, -5.5], [3, 4, -6], [180.48411560058594], 0, 1, -1, 1, 2],
    [[4.5, 4.5, -5.5], [4, 4, -6], [172.02670288085938], 0, 1, -1, 1, 2],
    [[5.5, 4.5, -5.5], [5, 4, -6], [164.41392517089844], 0, 1, -1, 1, 2],
    [[6.5, 4.5, -5.5], [6, 4, -6], [153.61619567871094], 0, 1, -1, 1, 2],
    [[7.5, 4.5, -5.5], [7, 4, -6], [143.60096740722656], 0, 1, -1, 1, 2],
    [[8.5, 4.5, -5.5], [8, 4, -6], [131.56103515625], 0, 1, -1, 1, 2],
    [[9.5, 4.5, -5.5], [9, 4, -6], [120.1473159790039], 0, 1, -1, 1, 2]
    ]
    },
    {
    "cells for seed cell index": 2,
    "structured seed cell xyz": [0.5, 0.5, 0.5],
    "structured seed cell ijk": [0, 0, 0],
    "structured seed cell complete": [[0.5, 0.5, 0.5], [0, 0, 0], [], 0, 1, -1, -1, 1],
    "geometry point used to find seed cell": [0.5, 0.5, 0.5],
    "number of cells collected for this seed cell": 20,
    "cell list for this seed cell": [
    [[0.5, -9.5, 0.5], [0, -10, 0], [175.404052734375], 0, 1, -1, 2, 2],
    [[0.5, -8.5, 0.5], [0, -9, 0], [179.89768981933594], 0, 1, -1, 2, 2],
    [[0.5, -7.5, 0.5], [0, -8, 0], [179.4281463623047], 0, 1, -1, 2, 2],
    [[0.5, -6.5, 0.5], [0, -7, 0], [201.8100128173828], 0, 1, -1, 2, 2],
    [[0.5, -5.5, 0.5], [0, -6, 0], [230.8294219970703], 0, 1, -1, 2, 2],
    [[0.5, -4.5, 0.5], [0, -5, 0], [235.75021362304688], 0, 1, -1, 2, 2],
    [[0.5, -3.5, 0.5], [0, -4, 0], [227.8350830078125], 0, 1, -1, 2, 2],
    [[0.5, -2.5, 0.5], [0, -3, 0], [238.86163330078125], 0, 1, -1, 2, 2],
    [[0.5, -1.5, 0.5], [0, -2, 0], [261.5865173339844], 0, 1, -1, 2, 2],
    [[0.5, -0.5, 0.5], [0, -1, 0], [262.8285217285156], 0, 1, -1, 2, 2],
    [[0.5, 0.5, 0.5], [0, 0, 0], [244.8736114501953], 0, 1, -1, 2, 2],
    [[0.5, 1.5, 0.5], [0, 1, 0], [241.0914306640625], 0, 1, -1, 2, 2],
    [[0.5, 2.5, 0.5], [0, 2, 0], [253.91702270507812], 0, 1, -1, 2, 2],
    [[0.5, 3.5, 0.5], [0, 3, 0], [250.46009826660156], 0, 1, -1, 2, 2],
    [[0.5, 4.5, 0.5], [0, 4, 0], [223.89569091796875], 0, 1, -1, 2, 2],
    [[0.5, 5.5, 0.5], [0, 5, 0], [206.5272979736328], 0, 1, -1, 2, 2],
    [[0.5, 6.5, 0.5], [0, 6, 0], [210.22640991210938], 0, 1, -1, 2, 2],
    [[0.5, 7.5, 0.5], [0, 7, 0], [204.92098999023438], 0, 1, -1, 2, 2],
    [[0.5, 8.5, 0.5], [0, 8, 0], [175.08787536621094], 0, 1, -1, 2, 2],
    [[0.5, 9.5, 0.5], [0, 9, 0], [149.2307586669922], 0, 1, -1, 2, 2]
    ]
    },
    {
    "cells for seed cell index": 3,
    "structured seed cell xyz": [4.5, 4.5, -5.5],
    "structured seed cell ijk": [4, 4, -6],
    "structured seed cell complete": [[4.5, 4.5, -5.5], [4, 4, -6], [], 0, 1, -1, -1, 1],
    "geometry point used to find seed cell": [5.0, 5.0, -5.0],
    "number of cells collected for this seed cell": 20,
    "cell list for this seed cell": [
    [[4.5, -9.5, -5.5], [4, -10, -6], [138.4343719482422], 0, 1, -1, 3, 2],
    [[4.5, -8.5, -5.5], [4, -9, -6], [139.57861328125], 0, 1, -1, 3, 2],
    [[4.5, -7.5, -5.5], [4, -8, -6], [135.86862182617188], 0, 1, -1, 3, 2],
    [[4.5, -6.5, -5.5], [4, -7, -6], [155.1950225830078], 0, 1, -1, 3, 2],
    [[4.5, -5.5, -5.5], [4, -6, -6], [181.42018127441406], 0, 1, -1, 3, 2],
    [[4.5, -4.5, -5.5], [4, -5, -6], [183.88121032714844], 0, 1, -1, 3, 2],
    [[4.5, -3.5, -5.5], [4, -4, -6], [173.90802001953125], 0, 1, -1, 3, 2],
    [[4.5, -2.5, -5.5], [4, -3, -6], [183.33624267578125], 0, 1, -1, 3, 2],
    [[4.5, -1.5, -5.5], [4, -2, -6], [204.96865844726562], 0, 1, -1, 3, 2],
    [[4.5, -0.5, -5.5], [4, -1, -6], [205.65621948242188], 0, 1, -1, 3, 2],
    [[4.5, 0.5, -5.5], [4, 0, -6], [187.70130920410156], 0, 1, -1, 3, 2],
    [[4.5, 1.5, -5.5], [4, 1, -6], [184.4735870361328], 0, 1, -1, 3, 2],
    [[4.5, 2.5, -5.5], [4, 2, -6], [198.39163208007812], 0, 1, -1, 3, 2],
    [[4.5, 3.5, -5.5], [4, 3, -6], [196.53305053710938], 0, 1, -1, 3, 2],
    [[4.5, 4.5, -5.5], [4, 4, -6], [172.02670288085938], 0, 1, -1, 3, 2],
    [[4.5, 5.5, -5.5], [4, 5, -6], [157.1180419921875], 0, 1, -1, 3, 2],
    [[4.5, 6.5, -5.5], [4, 6, -6], [163.6114044189453], 0, 1, -1, 3, 2],
    [[4.5, 7.5, -5.5], [4, 7, -6], [161.3614501953125], 0, 1, -1, 3, 2],
    [[4.5, 8.5, -5.5], [4, 8, -6], [134.768798828125], 0, 1, -1, 3, 2],
    [[4.5, 9.5, -5.5], [4, 9, -6], [112.26107025146484], 0, 1, -1, 3, 2]
    ]
    },
    {
    "cells for seed cell index": 4,
    "structured seed cell xyz": [0.5, 0.5, 0.5],
    "structured seed cell ijk": [0, 0, 0],
    "structured seed cell complete": [[0.5, 0.5, 0.5], [0, 0, 0], [], 0, 1, -1, -1, 2],
    "geometry point used to find seed cell": [0.5, 0.5, 0.5],
    "number of cells collected for this seed cell": 20,
    "cell list for this seed cell": [
    [[0.5, 0.5, -9.5], [0, 0, -10], [154.5513153076172], 0, 1, -1, 4, 2],
    [[0.5, 0.5, -8.5], [0, 0, -9], [166.3123321533203], 0, 1, -1, 4, 2],
    [[0.5, 0.5, -7.5], [0, 0, -8], [179.6851806640625], 0, 1, -1, 4, 2],
    [[0.5, 0.5, -6.5], [0, 0, -7], [198.03260803222656], 0, 1, -1, 4, 2],
    [[0.5, 0.5, -5.5], [0, 0, -6], [208.2532196044922], 0, 1, -1, 4, 2],
    [[0.5, 0.5, -4.5], [0, 0, -5], [216.9244384765625], 0, 1, -1, 4, 2],
    [[0.5, 0.5, -3.5], [0, 0, -4], [230.74717712402344], 0, 1, -1, 4, 2],
    [[0.5, 0.5, -2.5], [0, 0, -3], [236.71856689453125], 0, 1, -1, 4, 2],
    [[0.5, 0.5, -1.5], [0, 0, -2], [238.22740173339844], 0, 1, -1, 4, 2],
    [[0.5, 0.5, -0.5], [0, 0, -1], [244.8736114501953], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 0.5], [0, 0, 0], [244.8736114501953], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 1.5], [0, 0, 1], [238.22740173339844], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 2.5], [0, 0, 2], [236.71856689453125], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 3.5], [0, 0, 3], [230.74717712402344], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 4.5], [0, 0, 4], [216.9244384765625], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 5.5], [0, 0, 5], [208.2532196044922], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 6.5], [0, 0, 6], [198.03260803222656], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 7.5], [0, 0, 7], [179.6851806640625], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 8.5], [0, 0, 8], [166.3123321533203], 0, 1, -1, 4, 2],
    [[0.5, 0.5, 9.5], [0, 0, 9], [154.5513153076172], 0, 1, -1, 4, 2]
    ]
    },
    {
    "cells for seed cell index": 5,
    "structured seed cell xyz": [4.5, 4.5, -5.5],
    "structured seed cell ijk": [4, 4, -6],
    "structured seed cell complete": [[4.5, 4.5, -5.5], [4, 4, -6], [], 0, 1, -1, -1, 2],
    "geometry point used to find seed cell": [5.0, 5.0, -5.0],
    "number of cells collected for this seed cell": 20,
    "cell list for this seed cell": [
    [[4.5, 4.5, -9.5], [4, 4, -10], [128.51821899414062], 0, 1, -1, 5, 2],
    [[4.5, 4.5, -8.5], [4, 4, -9], [137.5346221923828], 0, 1, -1, 5, 2],
    [[4.5, 4.5, -7.5], [4, 4, -8], [148.25210571289062], 0, 1, -1, 5, 2],
    [[4.5, 4.5, -6.5], [4, 4, -7], [164.0957794189453], 0, 1, -1, 5, 2],
    [[4.5, 4.5, -5.5], [4, 4, -6], [172.02670288085938], 0, 1, -1, 5, 2],
    [[4.5, 4.5, -4.5], [4, 4, -5], [178.6822967529297], 0, 1, -1, 5, 2],
    [[4.5, 4.5, -3.5], [4, 4, -4], [190.81858825683594], 0, 1, -1, 5, 2],
    [[4.5, 4.5, -2.5], [4, 4, -3], [195.480224609375], 0, 1, -1, 5, 2],
    [[4.5, 4.5, -1.5], [4, 4, -2], [196.0938720703125], 0, 1, -1, 5, 2],
    [[4.5, 4.5, -0.5], [4, 4, -1], [202.28575134277344], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 0.5], [4, 4, 0], [202.28575134277344], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 1.5], [4, 4, 1], [196.0938720703125], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 2.5], [4, 4, 2], [195.480224609375], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 3.5], [4, 4, 3], [190.81858825683594], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 4.5], [4, 4, 4], [178.6822967529297], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 5.5], [4, 4, 5], [172.02670288085938], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 6.5], [4, 4, 6], [164.0957794189453], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 7.5], [4, 4, 7], [148.25210571289062], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 8.5], [4, 4, 8], [137.5346221923828], 0, 1, -1, 5, 2],
    [[4.5, 4.5, 9.5], [4, 4, 9], [128.51821899414062], 0, 1, -1, 5, 2]
    ]
    }
    ]
    }
    """

    goldJson = json.loads(goldTestStringJson)
    testFileName = testOutFileBasename + "00000.json"
    ff = open(testFileName, "r")
    testJson = json.load(ff)
    ff.close()
    self.assertEqual(goldJson, testJson)

    #remove file that got made during test
    os.remove(testFileName)


if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()
