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
from .PhactoriVectorLibrary import *
from paraview.simple import *
from .PhactoriSampledCellInfo import *
import vtk
import json
from .PhactoriParaviewMultiBlockRecursion import *

#phactori_combine_to_single_python_file_subpiece_begin_1

def GetCellTestPoint(theCell):
  theCellBounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  theCell.GetBounds(theCellBounds)
  retXyz = [0.5 * (theCellBounds[0] + theCellBounds[1]),
            0.5 * (theCellBounds[2] + theCellBounds[3]),
            0.5 * (theCellBounds[4] + theCellBounds[5])]
  return retXyz

def localGetCellijk(ii, jj, kk, inInputCsData, myExtent):
  #returnCell = inInputCsData.GetCell(ii, jj, kk)
  computedIndex = vtk.vtkStructuredData.ComputeCellIdForExtent(myExtent, [ii,jj,kk])
  returnCell = inInputCsData.GetCell(computedIndex)
  return returnCell

class PhactoriSampledCellInfo:
  def __init__(self):
    self.cellTestPoint = [0.0,0.0,0.0]
    self.pid = -1
    self.leafVisitCount = -1
    self.ijk = [-1,-1,-1]
    self.dataTuple = []
    self.index = -1
    self.segmentIndex = -1
    self.collectionAxis = 2

  def ToStr(self):
    outStr = "PhactoriSampledCellInfo:\n" +\
             "cellTestPoint: " + str(self.cellTestPoint) + "\n" +\
             "pid: " + str(self.pid) + "\n" +\
             "leafVisitCount: " + str(self.leafVisitCount) + "\n" +\
             "ijk: " + str(self.ijk) + "\n" +\
             "index: " + str(self.index) + "\n" +\
             "segmentIndex: " + str(self.segmentIndex) + "\n" +\
             "collectionAxis: " + str(self.collectionAxis) + "\n" +\
             "dataTuple: " + str(self.dataTuple) + "\n"
    return outStr

  def SetFromList(self, newValues):
    tripletVal = newValues[0]
    self.cellTestPoint[0] = tripletVal[0]
    self.cellTestPoint[1] = tripletVal[1]
    self.cellTestPoint[2] = tripletVal[2]
    tripletVal = newValues[1]
    self.ijk[0] = tripletVal[0]
    self.ijk[1] = tripletVal[1]
    self.ijk[2] = tripletVal[2]
    if newValues[2] == None:
      self.dataTuple = None
    else:
      self.dataTuple = list(newValues[2])
    self.pid = newValues[3]
    self.leafVisitCount = newValues[4]
    self.index = newValues[5]
    self.segmentIndex = newValues[6]
    self.collectionAxis = newValues[7]

  def GetAsList(self):
    return [self.cellTestPoint, self.ijk, self.dataTuple, self.pid, self.leafVisitCount, self.index, self.segmentIndex, self.collectionAxis]

  def ToStrTerseOneLineList(self):
    return str(self.GetAsList())
    #outStr = "[" + \
    #         str(self.cellTestPoint) + "," + \
    #         str(self.ijk) + "," + \
    #         str(self.dataTuple) + "," + \
    #         str(self.pid) + "," + \
    #         str(self.index) + "," + \
    #         str(self.segmentIndex) + "," + \
    #         str(self.collectionAxis) + \
    #         "]"
    #return outStr

  @staticmethod
  def TerseOneLineJsonFormatComment():
    outStr = '{"PhactoriSampledCellInfo output format 1 info":[\n' + \
             '" [cellTestPoint, ijk, dataTuple, pid, index, segmentIndex, collectionAxis]",\n' + \
             '" cellTestPoint is [X, Y, Z], ijk is [i, j, k], dataTuple is [c1, c2, ... cN]"]}'
    return outStr

  def SetIndex(self, newIndex):
    self.index = newIndex

  def SetCollectionAxis(self, newCollectionAxis):
    self.collectionAxis = newCollectionAxis

  def AxisCrossesExtent(self, myExtent, whichAxis):
    # i -> 0
    # j -> 1
    # k -> 2
    # ij -> 3
    # ik -> 4
    # jk -> 5
    if whichAxis < 3:
      if whichAxis != 0:
        if self.ijk[0] < myExtent[0]:
          return False
        if self.ijk[0] >= myExtent[1]:
          return False
      if whichAxis != 1:
        if self.ijk[1] < myExtent[2]:
          return False
        if self.ijk[1] >= myExtent[3]:
          return False
      if whichAxis != 2:
        if self.ijk[2] < myExtent[4]:
          return False
        if self.ijk[2] >= myExtent[5]:
          return False
    else:
      if whichAxis == 3:
        if self.ijk[2] < myExtent[4]:
          return False
        if self.ijk[2] >= myExtent[5]:
          return False
      if whichAxis == 4:
        if self.ijk[1] < myExtent[2]:
          return False
        if self.ijk[1] >= myExtent[3]:
          return False
      if whichAxis == 5:
        if self.ijk[0] < myExtent[0]:
          return False
        if self.ijk[0] >= myExtent[1]:
          return False
    return True

  def GetIntersectingCollectionExtent(self, whichExtent, whichAxis):
    # i -> 0
    # j -> 1
    # k -> 2
    # ij -> 3
    # ik -> 4
    # jk -> 5
    retext = list(whichExtent)
    if self.AxisCrossesExtent(whichExtent, whichAxis) == False:
      retExt[1] = retExt[0]
      retExt[3] = retExt[2]
      retExt[5] = retExt[4]
      return retExt
    if (whichAxis != 0) and (whichAxis != 3) and (whichAxis != 4):
      retext[0] = self.ijk[0]
      retext[1] = self.ijk[0] + 1
    if (whichAxis != 1) and (whichAxis != 3) and (whichAxis != 5):
      retext[2] = self.ijk[1]
      retext[3] = self.ijk[1] + 1
    if (whichAxis != 2) and (whichAxis != 4) and (whichAxis != 5):
      retext[4] = self.ijk[2]
      retext[5] = self.ijk[2] + 1
    return retext

  def SetCellTestPoint(self, inPt):
    self.cellTestPoint[0] = inPt[0]
    self.cellTestPoint[1] = inPt[1]
    self.cellTestPoint[2] = inPt[2]

  def SetIjk(self, ii, jj, kk):
    self.ijk[0] = ii
    self.ijk[1] = jj
    self.ijk[2] = kk

  def SetDataTuple(self, inDataTuple):
    self.dataTuple = list(inDataTuple)

  def SerializeAppendToFloatAndIntArray(self, outSerialFloatArray, outSerialIntArray, inTupleSize = 0):
    outSerialFloatArray.append(self.cellTestPoint[0])
    outSerialFloatArray.append(self.cellTestPoint[1])
    outSerialFloatArray.append(self.cellTestPoint[2])
    if self.dataTuple != None:
      localTupleLen = len(self.dataTuple)
      for ii in range(0,inTupleSize):
        if ii >= localTupleLen:
          outSerialFloatArray.append(0.0)
        else:
          outSerialFloatArray.append(self.dataTuple[ii])
    else:
      for ii in range(0,inTupleSize):
          outSerialFloatArray.append(0.0)
    outSerialIntArray.append(self.pid)
    outSerialIntArray.append(self.leafVisitCount)
    outSerialIntArray.append(self.ijk[0])
    outSerialIntArray.append(self.ijk[1])
    outSerialIntArray.append(self.ijk[2])
    outSerialIntArray.append(self.segmentIndex)

  def SerializeAppendToFloatAndIntArrayZeroVersion(self, outSerialFloatArray, outSerialIntArray, inTupleSize = -1):
    outSerialFloatArray.append(0.0)
    outSerialFloatArray.append(0.0)
    outSerialFloatArray.append(0.0)
    if inTupleSize < 0:
      for ffval in self.dataTuple:
        outSerialFloatArray.append(0.0)
    else:
      for ffval in range(0, inTupleSize):
        outSerialFloatArray.append(0.0)
    outSerialIntArray.append(0)
    outSerialIntArray.append(0)
    outSerialIntArray.append(0)
    outSerialIntArray.append(0)
    outSerialIntArray.append(0)
    outSerialIntArray.append(0)

  @staticmethod
  def GetSerializeFloatAndIntSize(tupleSize = 0):
    floatsize = 3 + tupleSize
    intsize = 6
    return floatsize, intsize

  def SerializeSetFromFloatAndIntArray(self, inSerialFloatArray, inSerialIntArray, inIndex, inTupleSize = 0):
    #floatsize, intsize = GetSerializeFloatAndIntSize()
    floatsize = 3 + inTupleSize
    intsize = 6
    floatIndex = floatsize * inIndex
    intIndex = intsize * inIndex
    self.cellTestPoint[0] = inSerialFloatArray[floatIndex+0]
    self.cellTestPoint[1] = inSerialFloatArray[floatIndex+1]
    self.cellTestPoint[2] = inSerialFloatArray[floatIndex+2]
    self.dataTuple = []
    if inTupleSize > 0:
      for ii in range(0, inTupleSize):
        self.dataTuple.append(inSerialFloatArray[floatIndex+3+ii])
    self.pid = inSerialIntArray[intIndex + 0]
    self.leafVisitCount = inSerialIntArray[intIndex + 1]
    self.ijk[0] = inSerialIntArray[intIndex + 2]
    self.ijk[1] = inSerialIntArray[intIndex + 3]
    self.ijk[2] = inSerialIntArray[intIndex + 4]
    self.segmentIndex = inSerialIntArray[intIndex + 5]

  def Populate1(self, ii, jj, kk, myPid, leafVisitCount, seedCellIndex,
        inInputCsData, myExtent, numCells, outputCellArray, dataArrayNumCmpnts, defaultTuple):
    #used during recursion to set up the cell based on the ii, jj, kk indices
    self.SetIjk(ii,jj,kk)
    self.pid = myPid
    self.leafVisitCount = leafVisitCount
    self.segmentIndex = seedCellIndex
    #dataCell = inInputCsData.GetCell(ii, jj, kk)
    dataCell = localGetCellijk(ii, jj, kk, inInputCsData, myExtent)
    cellTestPoint = GetCellTestPoint(dataCell)
    self.SetCellTestPoint(cellTestPoint)
    cellId = vtk.vtkStructuredData.ComputeCellIdForExtent(myExtent, [ii,jj,kk])
    self.SetIndex(cellId)
    if dataArrayNumCmpnts > 0:
      if (cellId < numCells) and (cellId >= 0):
        dataTuple = outputCellArray.GetTuple(cellId)
        self.SetDataTuple(dataTuple)
      else:
        self.SetDataTuple(defaultTuple)

class PhactoriFindCellWithMinMaxDataOnThisProcessRecursionParams:
  def __init__(self):
    self.leafVisitCount = 0
    self.dataArrayName = "noname"
    self.minCell = sys.float_info.max
    self.maxCell = None
    self.dataTotal = 0.0
    self.dataCount = 0
    self.currentMinVal = sys.float_info.max
    self.currentMaxVal = -sys.float_info.max

  def SetUpForRecursion(self, cellDataArrayName):
    self.dataTotal = 0.0
    self.dataCount = 0
    self.currentMinVal = sys.float_info.max
    self.currentMaxVal = -sys.float_info.max
    self.dataArrayName = cellDataArrayName
    self.minCell = PhactoriSampledCellInfo()
    self.minCell.SetFromList([[0.0,0.0,0.0], [-1,-1,-1], [self.currentMinVal], -1, -1, -1, -1, -1])
    self.maxCell = PhactoriSampledCellInfo()
    self.maxCell.SetFromList([[0.0,0.0,0.0], [-1,-1,-1], [self.currentMaxVal], -1, -1, -1, -1, -1])

def PhactoriFindCellWithMinMaxDataOnThisProcessInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellWithMinMaxDataOnThisProcessInBlock entered\n")

    inParameters.leafVisitCount += 1

    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    if (numCells == 0) or (numPoints == 0):
      if PhactoriDbg(100):
        myDebugPrint3("PhactoriFindCellWithMinMaxDataOnThisProcessInBlock returning 1\n")
      #no cells here
      return

    cellData = inInputCsData.GetCellData()
    if cellData == None:
      if PhactoriDbg(100):
        myDebugPrint3("PhactoriFindCellWithMinMaxDataOnThisProcessInBlock returning 2\n")

    outputCellArray = None
    outputCellArray = cellData.GetArray(inParameters.dataArrayName)

    if outputCellArray == None:
      if PhactoriDbg(100):
        myDebugPrint3("PhactoriFindCellWithMinMaxDataOnThisProcessInBlock returning 3\n")

    myPid = SmartGetLocalProcessId()
    for cellIndex in range(0,numCells):
      thisCellDataTuple = outputCellArray.GetTuple(cellIndex)
      cellDataVal = thisCellDataTuple[0]
      inParameters.dataTotal += cellDataVal
      inParameters.dataCount += 1
      if cellDataVal < inParameters.currentMinVal:
        inParameters.currentMinVal = cellDataVal
        oneCell = inInputCsData.GetCell(cellIndex)
        cellTestPoint = GetCellTestPoint(oneCell)
        inParameters.minCell.index = cellIndex
        inParameters.minCell.dataTuple[0] = cellDataVal
        inParameters.minCell.cellTestPoint[0] = cellTestPoint[0]
        inParameters.minCell.cellTestPoint[1] = cellTestPoint[1]
        inParameters.minCell.cellTestPoint[2] = cellTestPoint[2]
        inParameters.minCell.pid = myPid
        inParameters.leafVisitCount = inParameters.leafVisitCount
      if cellDataVal > inParameters.currentMaxVal:
        inParameters.currentMaxVal = cellDataVal
        oneCell = inInputCsData.GetCell(cellIndex)
        cellTestPoint = GetCellTestPoint(oneCell)
        inParameters.currentMinVal = cellDataVal
        cellTestPoint = GetCellTestPoint(oneCell)
        inParameters.maxCell.index = cellIndex
        inParameters.maxCell.dataTuple[0] = cellDataVal
        inParameters.maxCell.cellTestPoint[0] = cellTestPoint[0]
        inParameters.maxCell.cellTestPoint[1] = cellTestPoint[1]
        inParameters.maxCell.cellTestPoint[2] = cellTestPoint[2]
        inParameters.maxCell.pid = myPid
        inParameters.leafVisitCount = inParameters.leafVisitCount

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellWithMinMaxDataOnThisProcessInBlock returning\n")

def PhactoriFindCellWithMinMaxDataOnThisProcess(paraviewFilter, cellDataArrayName):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellWithMinMaxDataOnThisProcess entered\n" + \
        "paraviewFilter: " + str(paraviewFilter) + "\n")

    recursionParams = PhactoriFindCellWithMinMaxDataOnThisProcessRecursionParams()
    recursionParams.SetUpForRecursion(cellDataArrayName)
    recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    recursionObj.mParameters = recursionParams
    recursionObj.mOperationToDoPerBlock = PhactoriFindCellWithMinMaxDataOnThisProcessInBlock

    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, paraviewFilter)
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellWithMinMaxDataOnThisProcess returning\n")

    return [recursionParams.minCell, recursionParams.maxCell]

def PhactoriFindCellWtihMinMaxDataUsingMPI(paraviewFilter, cellDataArrayName):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellWtihMinMaxDataUsingMPI entered\n" + \
        "paraviewFilter: " + str(paraviewFilter) + "\n")

    thisPidMinMaxCells = PhactoriFindCellWithMinMaxDataOnThisProcess(paraviewFilter, cellDataArrayName)
    globalMinMaxCells = PhactoriLocalToGlobalCellsWithMinMaxDataUsingMPI(thisPidMinMaxCells, 1)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellWtihMinMaxDataUsingMPI returning\n")
    return globalMinMaxCells

def PhactoriLocalToGlobalCellsWithMinMaxDataUsingMPI(localPidMinMaxCellPair, tupleSize):

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriLocalToGlobalCellsWithMinMaxDataUsingMPI entered\n")
    #find overall min/max
    minVal = localPidMinMaxCellPair[0].dataTuple[0]
    maxVal = localPidMinMaxCellPair[1].dataTuple[0]
    if PhactoriDbg(100):
      myDebugPrint3("local min/max, tupleSize: " + str([minVal, maxVal, tupleSize]) + "\n")
    localMinMax = [-minVal, maxVal]
    globalMinMax = UseReduceOnFloatList(localMinMax, 0)
    globalMinVal = -globalMinMax[0]
    globalMaxVal = globalMinMax[1]
    if PhactoriDbg(100):
      myDebugPrint3("global min/max, tupleSize: " + str([globalMinVal, globalMaxVal, tupleSize]) + "\n")

    localPidMinMax = [-1,-1]
    myPid = SmartGetLocalProcessId()
    if globalMinVal == minVal:
      localPidMinMax[0] = myPid
    if globalMaxVal == maxVal:
      localPidMinMax[1] = myPid

    if PhactoriDbg(100):
      myDebugPrint3("localPidMinMax: " + str(localPidMinMax) + "\n")
    globalPidMinMax = UseReduceOnIntegerList(localPidMinMax, 0)
    if PhactoriDbg(100):
      myDebugPrint3("globalPidMinMax: " + str(globalPidMinMax) + "\n")

    localSerializedFloatArray = []
    localSerializedIntArray = []
    if globalPidMinMax[0] == myPid:
      localPidMinMaxCellPair[0].SerializeAppendToFloatAndIntArray(localSerializedFloatArray, localSerializedIntArray, tupleSize)
    else:
      localPidMinMaxCellPair[0].SerializeAppendToFloatAndIntArrayZeroVersion(localSerializedFloatArray, localSerializedIntArray, tupleSize)
    if globalPidMinMax[1] == myPid:
      localPidMinMaxCellPair[1].SerializeAppendToFloatAndIntArray(localSerializedFloatArray, localSerializedIntArray, tupleSize)
    else:
      localPidMinMaxCellPair[1].SerializeAppendToFloatAndIntArrayZeroVersion(localSerializedFloatArray, localSerializedIntArray, tupleSize)

    globalSerializedFloatArray = UseReduceOnFloatList(localSerializedFloatArray, 2)
    globalSerializedIntArray = UseReduceOnIntegerList(localSerializedIntArray, 2)

    globalMinCell = PhactoriSampledCellInfo()
    globalMaxCell = PhactoriSampledCellInfo()
    globalMinCell.SerializeSetFromFloatAndIntArray(globalSerializedFloatArray, globalSerializedIntArray, 0, tupleSize)
    globalMaxCell.SerializeSetFromFloatAndIntArray(globalSerializedFloatArray, globalSerializedIntArray, 1, tupleSize)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriLocalToGlobalCellsWithMinMaxDataUsingMPI returning\n")
    return [globalMinCell, globalMaxCell]

#phactori_combine_to_single_python_file_subpiece_end_1
