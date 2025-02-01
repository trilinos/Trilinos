
from phactori import *
from paraview.simple import *

#utilities for doing various data grid/geometric operations in parallel
#in particular:
#GetListOfGridPointsNearestListOfPointsV5
#  takes a list of geometric points, and goes through all processors and all
#  blocks and finds the data grid point which is nearest to each point in the list
#GetListOfCellTestPointsNearestListOfPointsV5
#  takes a list of geometric points, and goes through all processors and all
#  blocks and all cells in each block and for each geometric point finds the
#  cell bounding box center which is nearest to each point in the list

#phactori_combine_to_single_python_file_subpiece_begin_1

def GetGridPointsClosestToPointsInBlockV5(recursionObject, inInputCsData, inParameters):
  if PhactoriDbg(100):
    myDebugPrint3("GetGridPointsClosestToPointsInBlockV5 entered\n")
  numCells = inInputCsData.GetNumberOfCells()
  numPoints = inInputCsData.GetNumberOfPoints()
  if (numCells == 0) or (numPoints == 0):
    #no cells here
    if PhactoriDbg(100):
      myDebugPrint3("GetGridPointsClosestToPointsInBlockV5 returning (no cells or no points)\n")
    return

  if PhactoriDbg(100):
    myDebugPrint3(str(inParameters.testPointList) + "\n")
    myDebugPrint3(str(inParameters.distSqrdList) + "\n")

  pointsArray = inInputCsData.GetPoints()
  gridPtXyz = [0.0, 0.0, 0.0]
  for gridPtNdx in range(0,numPoints):
    pointsArray.GetPoint(gridPtNdx, gridPtXyz)
    for ptndx, oneTestPt in enumerate(inParameters.testPointList):
      testDist = vecDistanceSquared(oneTestPt, gridPtXyz)
      if testDist < inParameters.distSqrdList[ptndx]:
        inParameters.closestList[ptndx] = list(gridPtXyz)
        inParameters.distSqrdList[ptndx] = testDist

  if PhactoriDbg(100):
    myDebugPrint3(str(inParameters.testPointList) + "\n")
    myDebugPrint3(str(inParameters.distSqrdList) + "\n")

  if PhactoriDbg(100):
    myDebugPrint3("after this block:\n")
    for ii, oneGridPoint in enumerate(inParameters.closestList):
      myDebugPrint3(str(ii) + ": " + \
        str(inParameters.distSqrdList[ii]) + "\n" + \
        str(inParameters.testPointList[ii]) + "\n" + str(oneGridPoint))
    myDebugPrint3("\n")

  if PhactoriDbg(100):
    myDebugPrint3("GetGridPointsClosestToPointsInBlockV5 returning\n")

def GetCellsClosestToPointsInBlockV5(recursionObject, inInputCsData, inParameters):
  if PhactoriDbg(100):
    myDebugPrint3("GetCellsClosestToPointsInBlock entered\n")
  numCells = inInputCsData.GetNumberOfCells()
  numPoints = inInputCsData.GetNumberOfPoints()
  if (numCells == 0) or (numPoints == 0):
    #no cells here
    if PhactoriDbg(100):
      myDebugPrint3("GetCellsClosestToPointsInBlock returning (no cells or no points)\n")
    return

  if PhactoriDbg(100):
    myDebugPrint3(str(inParameters.testPointList) + "\n")
    myDebugPrint3(str(inParameters.distSqrdList) + "\n")

  for cellIndex in range(0,numCells):
    oneCell = inInputCsData.GetCell(cellIndex)
    cellTestPoint = GetCellTestPoint(oneCell)
    for ptndx, oneTestPt in enumerate(inParameters.testPointList):
      testDist = vecDistanceSquared(oneTestPt, cellTestPoint)
      if testDist < inParameters.distSqrdList[ptndx]:
        inParameters.closestList[ptndx] = cellTestPoint
        inParameters.distSqrdList[ptndx] = testDist

  if PhactoriDbg(100):
    myDebugPrint3(str(inParameters.testPointList) + "\n")
    myDebugPrint3(str(inParameters.distSqrdList) + "\n")

  if PhactoriDbg(100):
    myDebugPrint3("after this block:\n")
    for ii, oneCellPoint in enumerate(inParameters.closestList):
      myDebugPrint3(str(ii) + ": " + \
        str(inParameters.distSqrdList[ii]) + "\n" + \
        str(inParameters.testPointList[ii]) + "\n" + str(oneCellPoint))
    myDebugPrint3("\n")

  if PhactoriDbg(100):
    myDebugPrint3("GetCellsClosestToPointsInBlock returning\n")

class GetCellsClosestToPointsInBlockRecursionParamsV5:
  def __init__(self):
    self.testPointList = []
    self.distSqrdList = []
    self.closestList = []

  def InitializeWithPointList(self, inTestPointList):
    self.testPointList = inTestPointList
    numTestPoints = len(inTestPointList)
    for ii in range(0, numTestPoints):
      self.distSqrdList.append(sys.float_info.max)
      self.closestList.append(None)

def GetCellsClosestToPointsOnThisProcessFromParaViewFilterV5(inInputFilter, inTestPointList):
  if PhactoriDbg(100):
    myDebugPrint3("GetCellsClosestToPointsOnThisProcessFromParaViewFilter entered\n")
  recursionObj = PhactoriParaviewMultiBlockRecursionControl()
  recursionObj.mParameters = GetCellsClosestToPointsInBlockRecursionParamsV5()
  recursionObj.mParameters.InitializeWithPointList(inTestPointList)
  recursionObj.mOperationToDoPerBlock = GetCellsClosestToPointsInBlockV5
  PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, inInputFilter)
  if PhactoriDbg(100):
    myDebugPrint3("GetCellsClosestToPointsOnThisProcessFromParaViewFilter returning\n")
  return recursionObj.mParameters.closestList, recursionObj.mParameters.distSqrdList

def GetGridPointsClosestToPointsOnThisProcessFromParaViewFilterV5(inInputFilter, inTestPointList):
  if PhactoriDbg(100):
    myDebugPrint3("GetGridPointsClosestToPointsOnThisProcessFromParaViewFilterV5 entered\n")
  recursionObj = PhactoriParaviewMultiBlockRecursionControl()
  recursionObj.mParameters = GetCellsClosestToPointsInBlockRecursionParamsV5()
  recursionObj.mParameters.InitializeWithPointList(inTestPointList)
  recursionObj.mOperationToDoPerBlock = GetGridPointsClosestToPointsInBlockV5
  PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, inInputFilter)
  if PhactoriDbg(100):
    myDebugPrint3("GetGridPointsClosestToPointsOnThisProcessFromParaViewFilterV5 returning\n")
  return recursionObj.mParameters.closestList, recursionObj.mParameters.distSqrdList

def GetPidWithLeastValueListV5(inLocalDistSqrdList):
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


def UseMpiToGetGlobalCellPointsClosestV5(inInputFilter, inLocalCellPointList, inLocalDistSqrdList):
  if PhactoriDbg(100):
    myDebugPrint3("PhactoriSegmentCellSampler3.UseMpiToGetGlobalCellPointsClosest entered\n", 100)

  if PhactoriDbg(100):
    myDebugPrint3("inLocalCellPointList:\n" + str(inLocalCellPointList) + "\ninLocalDistSqrdList:\n" + str(inLocalDistSqrdList) + "\n")
  pidWithDataList, globalDistSqrdList = GetPidWithLeastValueListV5(inLocalDistSqrdList)
  if PhactoriDbg(100):
    myDebugPrint3("pidWithDataList:\n" + str(pidWithDataList) + "\nglobalDistSqrdList:\n" + str(globalDistSqrdList) + "\n")

  #convert cell point list to array of doubles and ints, use mpi reduce to share
  #the values, then convert back to cell point list
  serializeFloatArray = []
  serializeIntArray = []

  #convert cell point list to array of doubles
  cellPointFloatArray = []
  myPid = SmartGetLocalProcessId()
  for ii, oneCellPoint in enumerate(inLocalCellPointList):
    if pidWithDataList[ii] == myPid:
      cellPointFloatArray.append(oneCellPoint[0])
      cellPointFloatArray.append(oneCellPoint[1])
      cellPointFloatArray.append(oneCellPoint[2])
    else:
      cellPointFloatArray.append(0.0)
      cellPointFloatArray.append(0.0)
      cellPointFloatArray.append(0.0)

  #use mpi reduce to spread array correctly
  globalCellPointFloatArray = UseReduceOnFloatList(cellPointFloatArray, 2)

  #now create return global cell point list from arrays
  numCells = len(inLocalCellPointList)
  returnGlobalCellPointList = []
  for ii in range(0,numCells):
    myndx = ii*3
    oneCellPoint = [globalCellPointFloatArray[myndx],
                    globalCellPointFloatArray[myndx+1],
                    globalCellPointFloatArray[myndx+2]]
    returnGlobalCellPointList.append(oneCellPoint)

  if PhactoriDbg(100):
    myDebugPrint3("returnGlobalCellPointList:\n" + str(returnGlobalCellPointList) + "\n")
  if PhactoriDbg(100):
    myDebugPrint3("PhactoriSegmentCellSampler3.UseMpiToGetGlobalCellPointsClosest returning\n", 100)
  return returnGlobalCellPointList, globalDistSqrdList

def GetListOfCellTestPointsNearestListOfPointsV5(inInputFilter, pointList):
  """for each point in the list, find the cell test point (e.g. center of
     cell bounding box) which is nearest the test point.  Use MPI to work
     in parallel"""

  thisProcessNearestCellPointList, thisProcDistSqrdList = \
    GetCellsClosestToPointsOnThisProcessFromParaViewFilterV5(inInputFilter, pointList)

  nearestCellList, distanceList = UseMpiToGetGlobalCellPointsClosestV5(
    inInputFilter, thisProcessNearestCellPointList, thisProcDistSqrdList)

  return nearestCellList

def GetListOfGridPointsNearestListOfPointsV5(inInputFilter, pointList):
  """for each point in the list, find the point in the data grid
     which is nearest the test point.  Use MPI to work
     in parallel"""

  thisProcessNearestGridPointList, thisProcDistSqrdList = \
    GetGridPointsClosestToPointsOnThisProcessFromParaViewFilterV5(inInputFilter, pointList)

  nearestCellList, distanceList = UseMpiToGetGlobalCellPointsClosestV5(
    inInputFilter, thisProcessNearestGridPointList, thisProcDistSqrdList)

  return nearestCellList



#phactori_combine_to_single_python_file_subpiece_end_1
