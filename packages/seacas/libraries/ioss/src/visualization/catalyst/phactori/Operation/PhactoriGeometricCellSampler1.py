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
from .PhactoriOperationBlock import *
from .PhactoriVectorLibrary import *
from paraview.simple import *
from .PhactoriSegment import *
from .PhactoriMpiUtilities import *
from .PhactoriParaviewMultiBlockRecursion import *
from .PhactoriSampledCellInfo import *
from .PhactoriOutputFileNameAndDirectoryManager import *

import vtk
import json

#phactori_combine_to_single_python_file_subpiece_begin_1


class GatherGeometricallySampledCellsRecursionParams:
  def __init__(self):
    self.leafVisitCount = 0
    self.dataArrayNames = ["noname"]
    self.boundingBox = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]
    self.CellsMeetingGeometricCriteria = []

  def CellTestPointMeetsGeometricCriteria(self, cellTestPoint):
    if cellTestPoint[0] < self.boundingBox[0]:
      return False
    if cellTestPoint[0] > self.boundingBox[1]:
      return False
    if cellTestPoint[1] < self.boundingBox[2]:
      return False
    if cellTestPoint[1] > self.boundingBox[3]:
      return False
    if cellTestPoint[2] < self.boundingBox[4]:
      return False
    if cellTestPoint[2] > self.boundingBox[5]:
      return False
    return True

  def SetUpForRecursion(self, creatingSamplerInstance):
    self.dataArrayNames = creatingSamplerInstance.CellDataArrayNames
    self.boundingBox = creatingSamplerInstance.SamplingGeometryBoundingBox

class CollectDataOnSampledCellsRecursionParams:
  def __init__(self):
    self.leafVisitCount = 0
    self.dataArrayNames = ["noname"]
    self.GeometricallySampledCellsByRecursionLeaf = None

  def SetUpForRecursion(self, creatingSamplerInstance):
    self.dataArrayNames = creatingSamplerInstance.CellDataArrayNames
    self.GeometricallySampledCellsByRecursionLeaf = creatingSamplerInstance.GeometricallySampledCellsByRecursionLeaf


def localGetCellijk(ii, jj, kk, inInputCsData, myExtent):
  #returnCell = inInputCsData.GetCell(ii, jj, kk)
  computedIndex = vtk.vtkStructuredData.ComputeCellIdForExtent(myExtent, [ii,jj,kk])
  returnCell = inInputCsData.GetCell(computedIndex)
  return returnCell

class PhactoriGeometricCellSampler1(PhactoriOperationSpecifics):
  """Filter/operation which reads in a .json file which is a list of 3d
     coordinates pairs (representating segments, and then finds all cells
     in the input data that cross that segment and sets a masking value
     into a cell data array for those cells"""
  def __init__(self):
    self.myCopyOfInputFilter = None

    #indicates whether we should conduct our segment-picking cell calculations
    #by projecting onto an ordinal plane.  -1 is no projection, 0, 1, or 2 are
    #the x, y, or z axes. This is of particular use if you are doing a slice
    #before selecting elements.  The slice does not need to be axis aligned,
    #just choose the axis plane closest to parallel to the slice plane.
    self.ProjectionAxis = -1

    #as of 2019Aug13, we are having the programmable filter core dump with
    #structured data: this allows us to turn off the programmable filter to
    #do the other operations
    self.DoProgrammableFilterToAddMaskVariable = False
    self.mProgFilterMaskSetter = None
    self.mProgFilterString = None
    self.ProgrammableFilterOutputCellVariableName = "mask1"
    self.mainScriptString = None

    self.WriteOutGeometricallySampledCellData = True
    self.GeometricSampledCellOutputFilePathManager = PhactoriOutputFileNameAndDirectoryManager(
      "PhactoriGeometricCellSampler1_geometric_",
      "PhactoriGeometricCellSampler1_Output", 5, ".csv")

    self.WriteOutDataControlledSampledCellData = True
    self.DataControlledSampledCellOutputFilePathManager = PhactoriOutputFileNameAndDirectoryManager(
      "PhactoriGeometricCellSampler1_data_controlled_",
      "PhactoriGeometricCellSampler1_Output", 5, ".csv")

    self.CellDataArrayNames = ["V"]
    self.CellDataArrayTupleSize = 8

    self.SamplingGeometryBoundingBox = [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5]
    self.GeometricallySampledCellsForThisProcess = None
    self.GeometricallySampledCellsByRecursionLeaf = None

    self.DoPerCallbackDataControlledSampling = True
    self.DataControlledSampledCellsForThisProcess = None
    self.DataControlledSamplingMethod = "cells with ratio of min/max data value"
    #self.DataControlledSamplingMethod = "cells within distance of min/max highest data value cell"
    #self.DataControlledSamplingMethod = "cells within bounding box around min/max highest data value cell"
    self.DataControlledRatioBasis = "ratio is from zero to data maximum"
    #self.DataControlledRatioBasis = "ratio is from data minimum to data maximum"
    self.DataControlledRatio = 0.75
    self.DataControlledRatioGreaterOrLess = "cells greater/equal"
    #self.DataControlledRatioGreaterOrLess = "cells less/equal"
    #self.RatioBasis = "data range based"
    self.DataControlledDistance = 0.1
    self.DataControlledUseMinOrMax = "max"
    self.DataControlledBoundingBox = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]

  def CheckBoundingBox(self, boundingBox, parsekey, infoString):
    goodFlag = True
    if len(boundingBox) != 6:
      myDebugPrint3AndException(infoString + "\n"
        "bounding box is invalid for key: " + parseKey + "\n"
        "must be six float list, like:\n"
        "[-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]\n")
    if (boundingBox[0] > boundingBox[1]) or \
       (boundingBox[2] > boundingBox[3]) or \
       (boundingBox[4] > boundingBox[5]):
      myDebugPrint3AndException(infoString + "\n"
        "bounding box is invalid for key: " + parseKey + "\n"
        "we need max >= min for each axis, like:\n"
        "[-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]\n")

  def ParseParametersFromJson(self, inJson):
    #if "filename" in inJson:
    #  self.JsonListFileName = inJson['filename']
    #else:
    #  myDebugPrint3AndException(
    #      "PhactoriGeometricCellSampler1:ParseParametersFromJson\n"
    #      "Error:  must have 'filename' key\n")

    keyval1 = "sampling geometry bounding box"
    if keyval1 in inJson:
      newBB = inJson[keyval1]
      self.CheckBoundingBox(newBB, keyval1,
        "PhactoriGeometricCellSampler1:ParseParametersFromJson")
      self.SamplingGeometryBoundingBox = newBB

    keyval4 = "projection axis"
    if keyval4 in inJson:
      val4 = inJson[keyval4]
      if val4 == "none":
        self.ProjectionAxis = -1
      elif val4 == "x":
        self.ProjectionAxis = 0
      elif val4 == "y":
        self.ProjectionAxis = 1
      elif val4 == "z":
        self.ProjectionAxis = 2
      elif int(val4) == 1:
        self.ProjectionAxis = 0
      elif int(val4) == 2:
        self.ProjectionAxis = 1
      elif int(val4) == 3:
        self.ProjectionAxis = 2

    keyval5 = "do programmable filter"
    if keyval5 in inJson:
      self.DoProgrammableFilterToAddMaskVariable = inJson[keyval5]

    self.GeometricSampledCellOutputFilePathManager.ParseOutputFileNameControl(
      inJson, "geometric sampled cell output file ")

    self.DataControlledSampledCellOutputFilePathManager.ParseOutputFileNameControl(
      inJson, "data controlled sampled cell output file ")

    keyval7e = "export geometric sampled cell data as json file"
    if keyval7e in inJson:
      self.WriteOutGeometricallySampledCellData = inJson[keyval7e]

    keyval7f = "export data controlled sampled cell data as json file"
    if keyval7f in inJson:
      self.WriteOutDataControlledSampledCellData = inJson[keyval7f]

    keyval8 = "cell data array names"
    if keyval8 in inJson:
      self.CellDataArrayNames = inJson[keyval8]

    keyval9 = "cell data array tuple size"
    if keyval9 in inJson:
      self.CellDataArrayTupleSize = inJson[keyval9]

    keyval10 = "programmable filter output cell variable name"
    if keyval10 in inJson:
      self.ProgrammableFilterOutputCellVariableName = inJson[keyval10]

    keyval11 = "data controlled sampling method"
    if keyval11 in inJson:
      self.DataControlledSamplingMethod = inJson[keyval11]
      if self.DataControlledSamplingMethod == "cells within distance of min/max highest data value cell":
        keyval11b = "data controlled distance"
        if keyval11b not in inJson:
          myDebugPrint3AndException("must have 'data controlled distance'\n")
        self.DataControlledDistance = inJson[keyval11b]
      elif self.DataControlledSamplingMethod == "cells with ratio of min/max data value":
        keyval11b = "data controlled ratio of min/max"
        if keyval11b in inJson:
          self.DataControlledRatio = inJson[keyval11b]
          if (self.DataControlledRatio < 0.0) or (self.DataControlledRatio > 1.0):
            myDebugPrint3AndException("'data controlled ratio of min/max must be in range [0.0, 1.0]\n")
          if PhactoriDbg(100):
            myDebugPrint3("self.DataControlledRatio set to " + str(self.DataControlledRatio) + "\n", 100)
        keyval11c = "data controlled ratio basis"
        if keyval11c in inJson:
          self.DataControlledRatioBasis = inJson[keyval11c]
          if (self.DataControlledRatioBasis != "ratio is from zero to data maximum") and \
             (self.DataControlledRatioBasis != "ratio is from data minimum to data maximum"):
            myDebugPrint3("bad 'data controlled ratio basis'\n"
              "must be one of:\n"
              "'ratio is from zero to data maximum'\n"
              "'ratio is from data minimum to data maximum'\n")
        keyval11d = "collect cells relative to ratio"
        if keyval11d in inJson:
          self.DataControlledRatioGreaterOrLess = inJson[keyval11d]
          if (self.DataControlledRatioGreaterOrLess != "cells greater/equal") and \
             (self.DataControlledRatioGreaterOrLess != "cells less/equal"):
            myDebugPrint3("bad 'collect cells relative to ratio'\n"
              "must be one of:\n"
              "'cells greater/equal'\n"
              "'cells less/equal'\n")
      elif self.DataControlledSamplingMethod == "cells within bounding box around min/max highest data value cell":
        keyval11b = "data controlled bounding box"
        if keyval11b in inJson:
          self.DataControlledBoundingBox = inJson[keyval11b]
          self.CheckBoundingBox(self.DataControlledBoundingBox, keyval11b,
            "PhactoriGeometricCellSampler1:ParseParametersFromJson")
      else:
        myDebugPrint3AndException("bad 'data controlled sampling method'\n")
      if PhactoriDbg(100):
        myDebugPrint3("self.DataControlledSamplingMethod set to " + \
          str(self.DataControlledSamplingMethod) + "\n", 100)

    keyval11c = "data controlled sampling use min or max"
    if keyval11c in inJson:
      if inJson[keyval11c] == "max":
        self.DataControlledUseMinOrMax = "max"
      else:
        self.DataControlledUseMinOrMax = "min"

  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGeometricCellSampler1.CreateParaViewFilter entered\n", 100)

    savedActiveSource = GetActiveSource()

    self.myCopyOfInputFilter = inInputFilter

    UpdatePipelineWithCurrentTimeArgument(self.myCopyOfInputFilter)

    self.CreateInternalListOfGeometricallySampledCellsOnThisProcess()
    self.CreateInternalListOfDataControlledSampledCellsOnThisProcess()

    if self.DoProgrammableFilterToAddMaskVariable:
      self.mProgFilterMaskSetter = ProgrammableFilter(
              Input = self.myCopyOfInputFilter)
      self.mProgFilterMaskSetter.CopyArrays = 1
      self.CreateProgrammableFilterString()
      self.mProgFilterMaskSetter.Script = self.mProgFilterString
      self.mProgFilterMaskSetter.UpdatePipeline()
    else:
      self.mProgFilterMaskSetter = self.myCopyOfInputFilter

    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGeometricCellSampler1.CreateParaViewFilter returning\n", 100)
    return self.mProgFilterMaskSetter

  @staticmethod
  def GatherGeometricallySampledCellsInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("GatherGeometricallySampledCellsInBlock entered\n")

    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    inParameters.leafVisitCount += 1
    if (numCells == 0) or (numPoints == 0):
      #no cells here
      return

    cellData = inInputCsData.GetCellData()
    outputCellArrays = []
    dataArrayNumCmpnts = 0
    if cellData != None:
      for oneDataArrayName in inParameters.dataArrayNames:
        oneCellDataArray = cellData.GetArray(oneDataArrayName)
        outputCellArrays.append(oneCellDataArray)
        if oneCellDataArray != None:
          dataArrayNumCmpnts += oneCellDataArray.GetNumberOfComponents()
    if PhactoriDbg(100):
      myDebugPrint3("inParameters.dataArrayNames: " + str(inParameters.dataArrayNames) + "\n")
      for oneCellDataArray in outputCellArrays:
        if oneCellDataArray == None:
          myDebugPrint3("None\n")
        else:
          myDebugPrint3(str(oneCellDataArray.GetNumberOfTuples()) + "\n")
     
    if dataArrayNumCmpnts != 0:
      defaultTuple = []
      for ii in range(0, dataArrayNumCmpnts):
        defaultTuple.append(-1.0)
    else:
      dataArrayNumCmpnts = -1
      defaultTuple = []

    myPid = SmartGetLocalProcessId()
    for cellIndex in range(0,numCells):
      oneCell = inInputCsData.GetCell(cellIndex)
      cellTestPoint = GetCellTestPoint(oneCell)
      if inParameters.CellTestPointMeetsGeometricCriteria(cellTestPoint):
            if dataArrayNumCmpnts == 0:
              thisCellDataTuple = defaultTuple
            else:
              thisCellDataTuple = []
              for oneCellDataArray in outputCellArrays:
                if oneCellDataArray != None:
                  dataTupleSub1 = oneCellDataArray.GetTuple(cellIndex)
                  thisCellDataTuple.extend(dataTupleSub1)
            newCellInfo = PhactoriSampledCellInfo()
            newCellInfo.SetFromList([cellTestPoint, [-1,-1,-1], thisCellDataTuple, myPid, inParameters.leafVisitCount, cellIndex, -1, -1])
            inParameters.CellsMeetingGeometricCriteria.append(newCellInfo)


    if PhactoriDbg(100):
      myDebugPrint3("GatherGeometricallySampledCellsInBlock returning\n")

  def CreateInternalListOfGeometricallySampledCellsOnThisProcess(self):
    if PhactoriDbg(100):
      myDebugPrint3("CreateInternalListOfGeometricallySampledCellsOnThisProcess entered\n")

    recursionParams = GatherGeometricallySampledCellsRecursionParams()
    recursionParams.SetUpForRecursion(self)
    recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    recursionObj.mParameters = recursionParams
    recursionObj.mOperationToDoPerBlock = self.GatherGeometricallySampledCellsInBlock

    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.myCopyOfInputFilter)
    if PhactoriDbg(100):
      myDebugPrint3("CreateInternalListOfGeometricallySampledCellsOnThisProcess returning\n")

    self.GeometricallySampledCellsForThisProcess = recursionParams.CellsMeetingGeometricCriteria
    self.CreateLeafRecursionTrackingStructure()

  def PointIsInsideBoundingBox(self, testPoint, testBoundingBox):
    if testPoint[0] < testBoundingBox[0]:
      return False
    if testPoint[0] > testBoundingBox[1]:
      return False
    if testPoint[1] < testBoundingBox[2]:
      return False
    if testPoint[1] > testBoundingBox[3]:
      return False
    if testPoint[2] < testBoundingBox[4]:
      return False
    if testPoint[2] > testBoundingBox[5]:
      return False
    return True

  def CreateInternalListOfDataControlledSampledCellsOnThisProcess(self):
    if PhactoriDbg(100):
      myDebugPrint3("CreateInternalListOfDataControlledSampledCellsOnThisProcess entered\n")

    if self.DoPerCallbackDataControlledSampling == False:
      self.DataControlledSampledCellsForThisProcess = self.GeometricallySampledCellsForThisProcess
      if PhactoriDbg(100):
        myDebugPrint3("CreateInternalListOfDataControlledSampledCellsOnThisProcess returning (no processing)\n")
      return

    self.DataControlledSampledCellsForThisProcess = []

    #find min/max cell this process
    currentMax = -sys.float_info.max
    currentMin = sys.float_info.max
    if len(self.GeometricallySampledCellsForThisProcess) >= 1:
      for ii, oneCell in enumerate(self.GeometricallySampledCellsForThisProcess):
        testVal = oneCell.dataTuple[0]
        if testVal > currentMax:
          localMaxCell = oneCell
          currentMax = testVal
        if testVal < currentMin:
          localMinCell = oneCell
          currentMin = testVal
    else:
      myPid = SmartGetLocalProcessId()
      defaultDataTuple = []
      for ii in range(0, self.CellDataArrayTupleSize):
        defaultDataTuple.append(0.0)
      localMinCell = PhactoriSampledCellInfo()
      defaultDataTuple[0] = currentMax
      localMinCell.SetFromList([[0.0,0.0,0.0], [-1,-1,-1], defaultDataTuple, myPid, -1, -1, -1, -1])
      localMaxCell = PhactoriSampledCellInfo()
      defaultDataTuple[0] = currentMin
      localMaxCell.SetFromList([[0.0,0.0,0.0], [-1,-1,-1], defaultDataTuple, myPid, -1, -1, -1, -1])

    if PhactoriDbg(100):
      myDebugPrint3("localMinCell:\n" + localMinCell.ToStrTerseOneLineList() + "\n")
      myDebugPrint3("localMaxCell:\n" + localMaxCell.ToStrTerseOneLineList() + "\n")
    #find min/max cell all processes
    globalMinMaxCellPair = PhactoriLocalToGlobalCellsWithMinMaxDataUsingMPI([localMinCell, localMaxCell], self.CellDataArrayTupleSize)
    if PhactoriDbg(100):
      myDebugPrint3("globalMinCell:\n" + globalMinMaxCellPair[0].ToStrTerseOneLineList() + "\n")
      myDebugPrint3("globalMaxCell:\n" + globalMinMaxCellPair[1].ToStrTerseOneLineList() + "\n")
      myDebugPrint3("self.DataControlledSamplingMethod: " + self.DataControlledSamplingMethod + "\n")

    if self.DataControlledSamplingMethod == "cells with ratio of min/max data value":
      if self.DataControlledRatioBasis == "ratio is from zero to data maximum":
        maxCell = globalMinMaxCellPair[1]
        thresholdValue = self.DataControlledRatio * maxCell.dataTuple[0]
        if PhactoriDbg(100):
          myDebugPrint3("ratio is from zero to data maximum: " + str(maxCell.dataTuple[0]) + \
          ", " + str(self.DataControlledRatio) + "\n")
      else:
        minCell = globalMinMaxCellPair[0]
        maxCell = globalMinMaxCellPair[1]
        minValue = minCell.dataTuple[0]
        maxValue = maxCell.dataTuple[0]
        minMaxRange = maxValue - minValue
        thresholdValue = minValue + self.DataControlledRatio * minMaxRange
        if PhactoriDbg(100):
          myDebugPrint3("ratio is from data minimum to data maximum: " + str(minValue) + ", " + \
            str(maxValue) + ", " + str(minMaxRange) + ", " + \
            str(self.DataControlledRatio) + "\n")
      if PhactoriDbg(100):
        myDebugPrint3("thresholdValue: " + str(thresholdValue) + "\n")
      if self.DataControlledRatioGreaterOrLess == "cells greater/equal":
        for ii, oneCell in enumerate(self.GeometricallySampledCellsForThisProcess):
          if oneCell.dataTuple[0] >= thresholdValue:
            self.DataControlledSampledCellsForThisProcess.append(oneCell)
      else:
        for ii, oneCell in enumerate(self.GeometricallySampledCellsForThisProcess):
          if oneCell.dataTuple[0] <= thresholdValue:
            self.DataControlledSampledCellsForThisProcess.append(oneCell)
      if PhactoriDbg(100):
        myDebugPrint3("number of cells out of how many got sampled: " + \
          str(len(self.DataControlledSampledCellsForThisProcess)) + ", " + \
          str(len(self.GeometricallySampledCellsForThisProcess)) + "\n")
    elif self.DataControlledSamplingMethod == "cells within distance of min/max highest data value cell":
      if PhactoriDbg(100):
        myDebugPrint3("doing distance test: " + \
          str(self.DataControlledUseMinOrMax) + ", " + \
          str(self.DataControlledDistance) + "\n")
      if self.DataControlledUseMinOrMax == "max":
        centerCellPoint = globalMinMaxCellPair[1].cellTestPoint
      else:
        centerCellPoint = globalMinMaxCellPair[0].cellTestPoint
      if PhactoriDbg(100):
        myDebugPrint3("centerCellPoint: " + str(centerCellPoint) + "\n")
      for ii, oneCell in enumerate(self.GeometricallySampledCellsForThisProcess):
        cellDistance = vecDistance(centerCellPoint, oneCell.cellTestPoint)
        if cellDistance <= self.DataControlledDistance:
          self.DataControlledSampledCellsForThisProcess.append(oneCell)
    elif self.DataControlledSamplingMethod == "cells within bounding box around min/max highest data value cell":
      if PhactoriDbg(100):
        myDebugPrint3("doing bounding box test: " + \
          str(self.DataControlledUseMinOrMax) + ", " + \
          str(self.DataControlledBoundingBox) + "\n")
      if self.DataControlledUseMinOrMax == "max":
        centerCellPoint = globalMinMaxCellPair[1].cellTestPoint
      else:
        centerCellPoint = globalMinMaxCellPair[0].cellTestPoint
      if PhactoriDbg(100):
        myDebugPrint3("centerCellPoint: " + str(centerCellPoint) + "\n")
      testBoundingBox = [
        centerCellPoint[0] + self.DataControlledBoundingBox[0],
        centerCellPoint[0] + self.DataControlledBoundingBox[1],
        centerCellPoint[1] + self.DataControlledBoundingBox[2],
        centerCellPoint[1] + self.DataControlledBoundingBox[3],
        centerCellPoint[2] + self.DataControlledBoundingBox[4],
        centerCellPoint[2] + self.DataControlledBoundingBox[5]
        ]
      for ii, oneCell in enumerate(self.GeometricallySampledCellsForThisProcess):
        if self.PointIsInsideBoundingBox(oneCell.cellTestPoint, testBoundingBox):
          self.DataControlledSampledCellsForThisProcess.append(oneCell)
    else:
      myDebugPrint3AndException("CreateInternalListOfDataControlledSampledCellsOnThisProcess:\n" +\
        "bad self.DataControlledSamplingMethod\n")
 

    if PhactoriDbg(100):
      myDebugPrint3("CreateInternalListOfDataControlledSampledCellsOnThisProcess returning\n")

    #self.CreateLeafRecursionTrackingStructure()

  @staticmethod
  def CollectDataOnSampledCellsInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("CollectDataOnSampledCellsInBlock entered 1\n")

    inParameters.leafVisitCount += 1

    if inParameters.GeometricallySampledCellsByRecursionLeaf == None:
      #no cells collected for sampling on this process
      if PhactoriDbg(100):
        myDebugPrint3("CollectDataOnSampledCellsInBlock returning 1\n")
      return

    if inParameters.leafVisitCount not in inParameters.GeometricallySampledCellsByRecursionLeaf:
      #no cells collected for sampling on this block on this process
      if PhactoriDbg(100):
        myDebugPrint3("CollectDataOnSampledCellsInBlock returning 2\n")
      return

    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    if (numCells == 0) or (numPoints == 0):
      #no cells here
      if PhactoriDbg(100):
        myDebugPrint3("CollectDataOnSampledCellsInBlock returning 3\n")
      return

    cellData = inInputCsData.GetCellData()
    if cellData == None:
      if PhactoriDbg(100):
        myDebugPrint3("CollectDataOnSampledCellsInBlock returning 4\n")
      return

    testDataCellArray = None
    testDataCellArray = cellData.GetArray(inParameters.dataArrayNames[0])

    if testDataCellArray == None:
      if PhactoriDbg(100):
        myDebugPrint3("CollectDataOnSampledCellsInBlock returning 5\n")
      return

    outputCellArrays = []
    dataArrayNumCmpnts = 0
    if cellData != None:
      for oneDataArrayName in inParameters.dataArrayNames:
        oneCellDataArray = cellData.GetArray(oneDataArrayName)
        outputCellArrays.append(oneCellDataArray)
        if oneCellDataArray != None:
          dataArrayNumCmpnts += oneCellDataArray.GetNumberOfComponents()
    if PhactoriDbg(100):
      myDebugPrint3("inParameters.dataArrayNames: " + str(inParameters.dataArrayNames) + "\n")
      for oneCellDataArray in outputCellArrays:
        if oneCellDataArray == None:
          myDebugPrint3("None\n")
        else:
          myDebugPrint3(str(oneCellDataArray.GetNumberOfTuples()) + "\n")
     
    if dataArrayNumCmpnts != 0:
      defaultTuple = []
      for ii in range(0, dataArrayNumCmpnts):
        defaultTuple.append(-1.0)
    else:
      dataArrayNumCmpnts = -1
      defaultTuple = []

    myPid = SmartGetLocalProcessId()

    sampleCellList = inParameters.GeometricallySampledCellsByRecursionLeaf[inParameters.leafVisitCount]
    if PhactoriDbg(100):
      myDebugPrint3("num cells, leaf: " + str(len(sampleCellList)) + ", " + str(inParameters.leafVisitCount) + "\n")
    for oneSampleCell in sampleCellList:
      #thisCellDataTuple = outputCellArray.GetTuple(oneSampleCell.index)
      if dataArrayNumCmpnts == 0:
        thisCellDataTuple = defaultTuple
      else:
        thisCellDataTuple = []
        for oneCellDataArray in outputCellArrays:
          if oneCellDataArray != None:
            dataTupleSub1 = oneCellDataArray.GetTuple(oneSampleCell.index)
            thisCellDataTuple.extend(dataTupleSub1)
      tupleCount = min(len(thisCellDataTuple), len(oneSampleCell.dataTuple))
      if PhactoriDbg(100):
        myDebugPrint3(str(oneSampleCell.index) + ", " + str(thisCellDataTuple) + ", " + str(tupleCount) + "\n")
      for ii in range(0, tupleCount):
        oneSampleCell.dataTuple[ii] = thisCellDataTuple[ii]

    if PhactoriDbg(100):
      myDebugPrint3("CollectDataOnSampledCellsInBlock returning 6\n")

  def CollectDataOnSampledCellsOnThisProcess(self):
    """go through existing list of cells for this process, and get the data
       value for the current time into each cell sample instance"""
    if PhactoriDbg(100):
      myDebugPrint3("CollectDataOnSampledCellsOnThisProcess entered\n")

    recursionParams = CollectDataOnSampledCellsRecursionParams()
    recursionParams.SetUpForRecursion(self)
    recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    recursionObj.mParameters = recursionParams
    recursionObj.mOperationToDoPerBlock = self.CollectDataOnSampledCellsInBlock

    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.myCopyOfInputFilter)
    if PhactoriDbg(100):
      myDebugPrint3("CollectDataOnSampledCellsOnThisProcess returning\n")

  def CreateLeafRecursionTrackingStructure(self):
    if PhactoriDbg(100):
      myDebugPrint3("CreateLeafRecursionTrackingStructure entered\n")

    if len(self.GeometricallySampledCellsForThisProcess) <= 0:
      self.GeometricallySampledCellsByRecursionLeaf = None
      if PhactoriDbg(100):
       myDebugPrint3("CreateLeafRecursionTrackingStructure returning (no cells)\n")
      return

    self.GeometricallySampledCellsByRecursionLeaf = {}
    for oneCell in self.GeometricallySampledCellsForThisProcess:
      leafVal = oneCell.leafVisitCount
      if leafVal not in self.GeometricallySampledCellsByRecursionLeaf:
        self.GeometricallySampledCellsByRecursionLeaf[leafVal] = []
      listForLeaf = self.GeometricallySampledCellsByRecursionLeaf[leafVal]
      listForLeaf.append(oneCell)

    if PhactoriDbg(100):
      for leafKey in self.GeometricallySampledCellsByRecursionLeaf:
        cellListForLeaf = self.GeometricallySampledCellsByRecursionLeaf[leafKey]
        myDebugPrint3(str(len(cellListForLeaf)) + " cells for leaf " + str(leafKey) + "\n")
      
    if PhactoriDbg(100):
      myDebugPrint3("CreateLeafRecursionTrackingStructure returning\n")

  def WriteCellListToFile(self, filename, cellList):

    sumTuple = []
    for ii in range(0, self.CellDataArrayTupleSize):
      sumTuple.append(0.0)
    for oneCellInfo in cellList:
      dataTuple = oneCellInfo.dataTuple
      for ii in range(0, self.CellDataArrayTupleSize):
        sumTuple[ii] += dataTuple[ii]
    floatCellCount = float(len(cellList))
    averageTuple = []
    for oneItem in sumTuple:
      averageTuple.append(oneItem / floatCellCount)

    ff = open(filename, "w")
    ff.write("{\n")
    ff.write('"format for cells in lists": ' + PhactoriSampledCellInfo.TerseOneLineJsonFormatComment() + ",\n")
    ff.write('"sampling geometry bounding box": ' + str(self.SamplingGeometryBoundingBox) + ",\n")
    ff.write('"number of cells": ' + str(len(cellList)) + ",\n")
    ff.write('"cell variable names": ')
    ff.write('[')
    lastIndex = len(self.CellDataArrayNames) - 1
    for ii, oneArrayName in enumerate(self.CellDataArrayNames):
      if ii != lastIndex:
        ff.write('"' + str(oneArrayName) + '", ')
      else:
        ff.write('"' + str(oneArrayName) + '"')
    ff.write('],\n')
    ff.write('"data tuple size": ' + str(self.CellDataArrayTupleSize) + ",\n")
    ff.write('"sum variable values": ' + str(sumTuple) + ",\n")
    ff.write('"average variable values": ' + str(averageTuple) + ",\n")
    ff.write('"list of sampled cells": [\n')
    for cellCntr, oneCellInfo in enumerate(cellList):
      if cellCntr != 0:
        ff.write(",\n")
      cellStr = oneCellInfo.ToStrTerseOneLineList()
      ff.write(cellStr)
    ff.write("]\n")
    ff.write("}\n")

    ff.close()

  def GatherSampledCellDataFromAllProcessesUsingMPI(self, thisProcessCellList):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGeometricCellSampler1.GatherSampledCellDataFromAllProcessesUsingMPI entered\n", 100)

    myPid = SmartGetLocalProcessId()

    ##to start with, just have each process output it's data
    if PhactoriDbg(100):
      thisProcessHasData = False
      if thisProcessCellList != None:
        if len(thisProcessCellList) != 0:
          basename = "PhactoriGeometricCellSampler1_oneprocess_" + str(myPid) + ".json"
          self.WriteCellListToFile(pidfname, thisProcessCellList)

    if PhactoriDbg(100):
      myDebugPrint3("about to serialize cell infos\n", 100)

    tupleSize = self.CellDataArrayTupleSize
    cellListSerializedFloatArray = []
    cellListSerializedIntArray = []
    for idx, oneCellInfo in enumerate(thisProcessCellList):
      oneCellInfo.SerializeAppendToFloatAndIntArray(cellListSerializedFloatArray, cellListSerializedIntArray, tupleSize)

    if PhactoriDbg(100):
      myDebugPrint3("done serializing " + str(len(cellListSerializedFloatArray)) + "  " + str(len(cellListSerializedIntArray)) + "\n", 100)

    pidToDoOutput = SmartGetNumberOfProcesses() - 1
    globalSerializedFloatArray = UseMpiToSendAllProcessesFloatArrayToOneProcess(cellListSerializedFloatArray, pidToDoOutput)
    globalSerializedIntArray = UseMpiToSendAllProcessesIntArrayToOneProcess(cellListSerializedIntArray, pidToDoOutput)

    if PhactoriDbg(100):
      myDebugPrint3("returned from getting globalSerializedFloatArray and globalSerializedIntArray\n")

    globalCellList = []
    if myPid == pidToDoOutput:
      if PhactoriDbg(100):
       myDebugPrint3("I'm the output pid, about to serialize in\n")

      serialFloatSize, serialIntSize = PhactoriSampledCellInfo.GetSerializeFloatAndIntSize(tupleSize)
      globalNumCells = len(globalSerializedFloatArray) / serialFloatSize

      for ii in range(0, globalNumCells):
        newCell = PhactoriSampledCellInfo()
        newCell.SerializeSetFromFloatAndIntArray(globalSerializedFloatArray, globalSerializedIntArray, ii, tupleSize)
        globalCellList.append(newCell)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGeometricCellSampler1.GatherSampledCellDataFromAllProcessesUsingMPI returning\n", 100)

    return globalCellList


  def GatherSampledCellDataFromAllProcessesAndWriteFromOneProcessUsingMPI(self, outputFilenameCounter = -1):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGeometricCellSampler1.GatherSampledCellDataFromAllProcessesAndWriteFromOneProcessUsingMPI entered\n", 100)
    myPid = SmartGetLocalProcessId()

    if outputFilenameCounter < 0:
      ftc = GetFrameTagCounter()
    else:
      ftc = outputFilenameCounter


    pidToDoOutput = SmartGetNumberOfProcesses() - 1
    if self.WriteOutGeometricallySampledCellData == True:
      globalCellList = self.GatherSampledCellDataFromAllProcessesUsingMPI(self.GeometricallySampledCellsForThisProcess)
      if myPid == pidToDoOutput:
        if PhactoriDbg(100):
          myDebugPrint3("finished serialize, about to write to geometric sample cell file\n" +\
            "number of cells from all processes: " + str(len(globalCellList)) + "\n")

        fname = self.GeometricSampledCellOutputFilePathManager.GetOutputFilePath(ftc, "")
        if PhactoriDbg(100):
          myDebugPrint3("output filename will be: " + str(fname) + "\n")
        self.WriteCellListToFile(fname, globalCellList)

    if self.WriteOutDataControlledSampledCellData == True:
      globalCellList = self.GatherSampledCellDataFromAllProcessesUsingMPI(self.DataControlledSampledCellsForThisProcess)
      if myPid == pidToDoOutput:
        if PhactoriDbg(100):
          myDebugPrint3("finished serialize, about to write to data controlled sample cell file\n" +\
            "number of cells from all processes: " + str(len(globalCellList)) + "\n")

        fname = self.DataControlledSampledCellOutputFilePathManager.GetOutputFilePath(ftc, "")
        if PhactoriDbg(100):
          myDebugPrint3("output filename will be: " + str(fname) + "\n")
        self.WriteCellListToFile(fname, globalCellList)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGeometricCellSampler1.GatherSampledCellDataFromAllProcessesAndWriteFromOneProcessUsingMPI returning\n", 100)

  def ExportOperationData(self, datadescription):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGeometricCellSampler1.ExportOperationData entered\n", 100)

    UpdatePipelineWithCurrentTimeArgument(self.myCopyOfInputFilter)

    #if we are supposed to resample data after first timestep, we call Create again
    #self.CreateInternalListOfGeometricallySampledCellsOnThisProcess
    #the idea here is that we already have our cell collection and we want to pull the data
    self.CollectDataOnSampledCellsOnThisProcess()

    self.CreateInternalListOfDataControlledSampledCellsOnThisProcess()

    self.GatherSampledCellDataFromAllProcessesAndWriteFromOneProcessUsingMPI()

    if self.DoProgrammableFilterToAddMaskVariable:
      self.CreateProgrammableFilterString()
      self.mProgFilterMaskSetter.Script = self.mProgFilterString
      self.mProgFilterMaskSetter.UpdatePipeline()

    #if PhactoriDbg(100):
    #  minCellmaxCellPair = PhactoriFindCellWithMinMaxDataOnThisProcess(self.myCopyOfInputFilter, "eqps")
    #  myDebugPrint3("min cell this process:\n" + minCellmaxCellPair[0].ToStr())
    #  myDebugPrint3("max cell this process:\n" + minCellmaxCellPair[1].ToStr())
    #  globalMinCellmaxCellPair = PhactoriFindCellWtihMinMaxDataUsingMPI(self.myCopyOfInputFilter, "eqps")
    #  myDebugPrint3("min cell global:\n" + globalMinCellmaxCellPair[0].ToStr())
    #  myDebugPrint3("max cell global:\n" + globalMinCellmaxCellPair[1].ToStr())

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGeometricCellSampler1.ExportOperationData returning\n", 100)

  def CreateProgrammableFilterString(self):

    if self.mainScriptString == None:
      scriptLines = []

      scriptLines.append("#print('test2')\n")

      scriptLines.append("#cellIndexesToSet = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]\n")
      scriptLines.append("#maskValuesToSet = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]\n")
      scriptLines.append("#leafVisitCount = [0,0,0,0,0,0,0,0,0,0,0]\n")

      scriptLines.append("localLeafVisitCount = 0\n")

      scriptLines.append("def flatten(input, output):\n")
      scriptLines.append("    # Copy the cells etc.\n")
      scriptLines.append("    output.ShallowCopy(input)\n")
      scriptLines.append("    numPoints = input.GetNumberOfPoints()\n")
      scriptLines.append("    celldata = output.GetCellData()\n")
      scriptLines.append("    numCells = input.GetNumberOfCells()\n")
      scriptLines.append("    #print(str(numCells))\n")
      scriptLines.append("    ncda = vtk.vtkDoubleArray()\n")
      scriptLines.append("    ncda.SetNumberOfTuples(numCells)\n")
      scriptLines.append("    numTuples = ncda.GetNumberOfTuples()\n")
      scriptLines.append("    for ii in range(0, numTuples):\n")
      scriptLines.append("        ncda.SetValue(ii, 0.0)\n")
      scriptLines.append("    for ii, cellNdx in enumerate(cellIndexesToSet):\n")
      scriptLines.append("        if localLeafVisitCount == leafVisitCount[ii]:\n")
      scriptLines.append("          if cellNdx < numCells:\n")
      scriptLines.append("            ncda.SetValue(cellNdx, float(maskValuesToSet[ii]))\n")
      scriptLines.append("    ncda.SetName('" + self.ProgrammableFilterOutputCellVariableName + "')\n")
      scriptLines.append("    #print(ncda.GetName())\n")
      scriptLines.append("    celldata.AddArray(ncda)\n")

      scriptLines.append("input = self.GetInputDataObject(0, 0)\n")
      scriptLines.append("output = self.GetOutputDataObject(0)\n")

      scriptLines.append("if input.IsA('vtkMultiBlockDataSet'):\n")
      scriptLines.append("    output.CopyStructure(input)\n")
      scriptLines.append("    iter = input.NewIterator()\n")
      scriptLines.append("    iter.UnRegister(None)\n")
      scriptLines.append("    iter.InitTraversal()\n")
      scriptLines.append("    while not iter.IsDoneWithTraversal():\n")
      scriptLines.append("        localLeafVisitCount += 1\n")
      scriptLines.append("        curInput = iter.GetCurrentDataObject()\n")
      scriptLines.append("        curOutput = curInput.NewInstance()\n")
      scriptLines.append("        curOutput.UnRegister(None)\n")
      scriptLines.append("        output.SetDataSet(iter, curOutput)\n")
      scriptLines.append("        flatten(curInput, curOutput)\n")
      scriptLines.append("        iter.GoToNextItem();\n")
      scriptLines.append("else:\n")
      scriptLines.append("  flatten(input, output)\n")
  
      self.mainScriptString = "".join(scriptLines)

    myCellIndexesToSet = []
    myMaskValuesToSet = []
    myLeafVisitCount = []
    if self.GeometricallySampledCellsForThisProcess != None:
      for cellIdx, oneCellInfo in enumerate(self.GeometricallySampledCellsForThisProcess):
         maskValue = 1.0
         myCellIndexesToSet.append(oneCellInfo.index)
         myMaskValuesToSet.append(maskValue)
         myLeafVisitCount.append(oneCellInfo.leafVisitCount)
    if self.DataControlledSampledCellsForThisProcess != None:
      for cellIdx, oneCellInfo in enumerate(self.DataControlledSampledCellsForThisProcess):
         maskValue = 2.0
         myCellIndexesToSet.append(oneCellInfo.index)
         myMaskValuesToSet.append(maskValue)
         myLeafVisitCount.append(oneCellInfo.leafVisitCount)
   
    newstr1 = "cellIndexesToSet = " + str(myCellIndexesToSet) + "\n" 
    newstr2 = "maskValuesToSet = " + str(myMaskValuesToSet) + "\n" 
    newstr3 = "leafVisitCount = " + str(myLeafVisitCount) + "\n" 
    self.mProgFilterString = newstr1 + newstr2 + newstr3 + self.mainScriptString

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMarkCellSetsOperation constructed script:\n" + \
        self.mProgFilterString)

#phactori_combine_to_single_python_file_subpiece_end_1

