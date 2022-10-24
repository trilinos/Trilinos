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
from .PhactoriSegment import *
from .PhactoriMpiUtilities import *
import json

#phactori_combine_to_single_python_file_subpiece_begin_1

global FlagToTestOutgoingPvGeometryFilter
FlagToTestOutgoingPvGeometryFilter = False

global WriteEachDeadCellElementToFiles
WriteEachDeadCellElementToFiles = False

class BlockRecursionControlItem:
  """see DoMethodPerBlock(); mOperationPerBlock should be set to a method
     which takes 1 parameter, mParameters should be set to the parameter
     instance which will be passes to the mOperationPerBlock call"""
  def __init__(self):
    self.mOperationPerBlock = None
    self.mParameters = None

class PhactoriOperationBlock:
  """manages one stage of the data pipeline, analagous to ParaView Filter

  An instance of this class represents and manages one stage of the data
  pipeline which has been set up for management by phatori.  It creates and
  holds a reference to a ParaView/Catalyst/vtk Filter which is doing the
  real work--this class might be thought of as an adapter or interface from
  phactori to ParaView.  As currently implemented, an operation is assumed to
  have only one input and one output.  Multiple operations can have the same
  input, so a tree structure is allowed rather than just a linear pipe.
  Operations with multiple inputs and outputs are conceiveable, and may be
  added pending user requirements.
  The instance is presumed to contain a name unique amound the operation
  blocks and keeps a reference to the input operation (by name), the
  ParaView/Catalyst filter which is built, and some flags determining where
  we are in the construction process.
  """
  def __init__(self):
    self.mName = ""
    self.mInputOperationName = None
    self.mType = ""
    self.mHasBeenConstructed = None
    self.mParaViewFilter = None
    #for keeping track of data bounds for this operation and only updating
    #it when necessary
    self.mDataBoundsIsCurrent = False
    self.mDataBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    self.mOperationSpecifics = None

    #If this is a purely internal pipeline filter, this will remain None.
    #If somebody uses the filter for rendering, it will be created and
    #handed off for rendering--only needs creation once.
    self.mOutgoingGeometryFilter = None

    #stores annotation time source, if one gets created associated with this
    #operation
    self.mParaViewTimeAnnotationSource = None

    self.mRecursionParameterStore = {}


  def GetPvFilter(self):
    return self.mParaViewFilter

  def CreateOutgoingPvGeometryFilter(self):
    savedActiveSource = GetActiveSource()
    self.mOutgoingGeometryFilter = ExtractSurface(self.GetPvFilter())
    if PhactoriDbg():
      myDebugPrint3("PhactoriOperationBlock::CreateOutgoingPvGeometryFilter\n"
          "created self.mOutgoingGeometryFilter: " + \
      str(self.mOutgoingGeometryFilter) + "\n")

    SetActiveSource(savedActiveSource)

  def GetOutgoingPvGeometryFilter(self):
    global FlagToTestOutgoingPvGeometryFilter
    if FlagToTestOutgoingPvGeometryFilter == False:
      if PhactoriDbg():
        myDebugPrint3("PhactoriOperationBlock::GetOutgoingPvGeometryFilter\n"
            "not creating or using external geometry filter:\n")
      return self.GetPvFilter()
    if self.mOutgoingGeometryFilter == None:
      self.CreateOutgoingPvGeometryFilter()
    return self.mOutgoingGeometryFilter

  def GetListOfInputOperationNames(self):
    retList = self.mOperationSpecifics.GetListOfInputOperationNamesForThisOperationType()
    if self.mInputOperationName != None:
      retList.append(self.mInputOperationName)
    return retList

  def DoUpdateDueToChangeInData(self, ioPipeAndViewsState):
    outputPvFilter = self.mParaViewFilter
    if self.mInputOperationName == None:
      inputOperation = ioPipeAndViewsState.mIncomingDefaultOperation
    else:
      inputOperation = ioPipeAndViewsState.mOperationBlocks[
          self.mInputOperationName]
    if inputOperation == None:
      myDebugPrint3AndException(
        "PhactoriOperationBlock::DoUpdateDueToChangeInData:\n"
        "couldn't find input operation with name: " + \
        str(self.mInputOperationName) + "\n")

    self.mOperationSpecifics.DoUpdateDueToChangeInData(
        inputOperation.mParaViewFilter, outputPvFilter)

  def DoMethodPerBlockFromParaViewFilter(self, inRecursionControlItem, inPvFilter):
    pvClientSideData = inPvFilter.GetClientSideObject().GetOutputDataObject(0)
    if pvClientSideData == None:
      if PhactoriDbg(100):
        myDebugPrint3(
          'DoMethodPerBlock: pvClientSideData is None, returning',100)

    self.DoMethodPerBlockRecurse1(pvClientSideData, inRecursionControlItem)

  def DoMethodPerBlock(self, inRecursionControlItem):
    """DoMethodPerBlock is a generic method for doing recursion through the
       multiblock dataset and doing something (a callback) on a per leaf block
       basis.  Gets the clientside data and calls DoMethodPerBlockRecurse1
       which calls itself on internal nodes and calls
       inRecursionControlItem.mOperationToDoPerBlock on leaf block nodes"""
    pvFilter = self.GetPvFilter()
    self.DoMethodPerBlockFromParaViewFilter(inRecursionControlItem, pvFilter)

  def DoMethodPerBlockRecurse1(self, inInputCsData, inRecursionControlItem):
    """DoMethodPerBlockRecurse1 is a generic method for doing recursion through
       the multiblock dataset and doing something (a callback) on a per leaf block
       basis.  Called by DoMethodPerBlock which got clientside data and calls
       itself on internal nodes and calls
       inRecursionControlItem.mOperationToDoPerBlock on leaf block nodes"""
    #if PhactoriDbg(100):
    #  myDebugPrint3('DoMethodPerBlockRecurse1 entered\n', 100)

    icsdClassname = inInputCsData.GetClassName()
    if icsdClassname == "vtkMultiBlockDataSet" or \
       icsdClassname == "vtkExodusIIMultiBlockDataSet":
      #myDebugPrint3('recursing: ' + icsdClassname + '\n')
      numBlocks = inInputCsData.GetNumberOfBlocks()
      for ii in range(0, numBlocks):
        oneBlock = inInputCsData.GetBlock(ii)
        if(oneBlock != None):
          self.DoMethodPerBlockRecurse1(oneBlock, inRecursionControlItem)
    else:
      inRecursionControlItem.mOperationToDoPerBlock(inInputCsData,
              inRecursionControlItem.mParameters)

    #if PhactoriDbg(100):
    #  myDebugPrint3('DoMethodPerBlockRecurse1 returning\n', 100)

  def OutputElementListFromOneBlockToFile(self, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("OutputElementListFromOneBlockToFile entered\n")
    cellData = inInputCsData.GetCellData()
    numTuplesX = cellData.GetNumberOfTuples()
    numArrays = cellData.GetNumberOfArrays()

    global WriteEachDeadCellElementToFiles
    #if inParameters.mBlockCount == 0:
    if (inParameters.mFlag1 == 0) and (numTuplesX > 0):
      inParameters.mFlag1 = 1
      if WriteEachDeadCellElementToFiles:
        inParameters.mOutFileff.write("element index")
        for arrayIndex in range(0, numArrays):
            inParameters.mOutFileff.write(
              ", " + str(cellData.GetArray(arrayIndex).GetName()))
        inParameters.mOutFileff.write("\n")

    if WriteEachDeadCellElementToFiles:
      for cellIndex in range(0, numTuplesX):
          inParameters.mOutFileff.write(
                  str(cellIndex + inParameters.mElementCount))
          for arrayIndex in range(0, numArrays):
              inParameters.mOutFileff.write(
                ", " + str(cellData.GetArray(arrayIndex).GetTuple1(cellIndex)))
          inParameters.mOutFileff.write("\n")

    killed_criteria_array = cellData.GetArray("KILLED")
    for cellIndex in range(0, numTuplesX):
      kval = killed_criteria_array.GetTuple1(cellIndex)
      kvalindex = int(kval)
      if (kvalindex >= 0) and (kvalindex < 10):
        inParameters.mKilledByCriteriaCount[kvalindex] += 1


    inParameters.mElementCount += inInputCsData.GetNumberOfCells()
    inParameters.mNodeCount += inInputCsData.GetNumberOfPoints()
    inParameters.mBlockCount += 1

    if PhactoriDbg(100):
      myDebugPrint3("OutputElementListFromOneBlockToFile returning\n")

  class FindClosestNPointsToListParams:
    """recursion structure for FindClosestNPointsToList().  Also servers to
       store/track data for passing back answer"""
    def __init__(self, inParentOperation, inNumToFind,
                 inTargetGlobalNodeIdList, inTargetPointXyzList):
      self.mParentOperation = inParentOperation
      self.mNumToFind = inNumToFind

      #list of points we need to find closest local points to
      #(global node ids and geometry xyzs)
      self.mTargetIds = inTargetGlobalNodeIdList
      self.mTargetXyzs = inTargetPointXyzList

      self.mKdtree = None

      #closest N local points (local id and id from list and geometry xyzs
      #from local and from list and distance squared
      self.mThisProcIds = vtk.vtkIntArray()
      self.mTargetMatchIds = vtk.vtkIntArray()
      self.mThisProcXyzs = vtk.vtkDoubleArray()
      self.mTargetMatchXyzs = vtk.vtkDoubleArray()
      self.mDistSqrds = vtk.vtkDoubleArray()

      self.mThisProcIds.SetNumberOfValues(inNumToFind)
      self.mTargetMatchIds.SetNumberOfValues(inNumToFind)
      self.mThisProcXyzs.SetNumberOfValues(inNumToFind*3)
      self.mTargetMatchXyzs.SetNumberOfValues(inNumToFind*3)
      self.mDistSqrds.SetNumberOfValues(inNumToFind)

      #set default values indicating nothing found in those entries
      for ii in range(0,inNumToFind):
        self.mThisProcIds.SetValue(ii, -1)
        self.mTargetMatchIds.SetValue(ii, -1)
        self.mDistSqrds.SetValue(ii, sys.float_info.max)
      self.mMinDistSqrd = sys.float_info.max
      #index of the item that currently has the biggest distance
      self.mcfndx = 0

    def SetUpKdtree(self):
      """take the points in this instance (presumably the target object) and
         put them in a kdtree (using scipy) so we can find the closest point
         in log(n) time"""
      from scipy import spatial
      kdtreepts = []
      numpts = self.mTargetIds.GetNumberOfValues()
      for ii in range(0,numpts):
        pndx = ii*3
        kdtreepts.append([
          self.mTargetXyzs.GetValue(pndx),
          self.mTargetXyzs.GetValue(pndx+1),
          self.mTargetXyzs.GetValue(pndx+2)])
      self.mKdtree = None
      self.mKdtree = spatial.KDTree(kdtreepts)

    def TestPointWithKdtree(self, inSrcId, inSrcXyz):
      nrstdist, nrstndx = self.mKdtree.query(inSrcXyz)
      #print "nrstdist: ", str(nrstdist)
      #print "nrstndx: ", str(nrstndx)
      tgtId = self.mTargetIds.GetValue(nrstndx)
      tgtXyz = self.mKdtree.data[nrstndx]
      #self.TestPointSub1(inSrcId, inSrcXyz[0], inSrcXyz[1], inSrcXyz[2],
      #    tgtId, tgtXyz[0], tgtXyz[1], tgtXyz[2])
      self.TestPointSub2(inSrcId, inSrcXyz, tgtId, tgtXyz, nrstdist)

    def TestPoint(self, inId, inXyz):
      """given an xyz point in the local processor (and its global node id)
         see if it is closer to any of the target points than the current
         set of nearest points and, if so, put it in the set, dropping
         others out if necessary"""
      if PhactoriDbg():
        myDebugPrint3("TestPoint id " + str(inId))
      tgtxyzs = self.mTargetXyzs
      numtgtpts = self.mTargetIds.GetNumberOfValues()
      for pp in range(0, numtgtpts):
        tgtid = self.mTargetIds.GetValue(pp)
        ndx = pp*3
        self.TestPointSub1(inId, inXyz[0], inXyz[1], inXyz[2],
            tgtid, tgtxyzs.GetValue(ndx), tgtxyzs.GetValue(ndx+1),
            tgtxyzs.GetValue(ndx+2))

    def TestPointSub2(self, inSrcId, inSrcXyz, inTgtId, inTgtXyz, inNrstdist):
      dstsqd = inNrstdist*inNrstdist
      if dstsqd >= self.mMinDistSqrd:
        #we already have mNumToFind closer
        #if PhactoriDbg():
        #  myDebugPrint3("inSrcId inTgtId " + str(inSrcId) + " " + \
        #      str(inTgtId) + \
        #      " too far: " + str(dstsqd) + " >= " + \
        #      str(self.mMinDistSqrd) + "\n")
        return

      #replace the previous point that was furthest
      self.mDistSqrds.SetValue(self.mcfndx, dstsqd)
      self.mThisProcIds.SetValue(self.mcfndx, inSrcId)
      self.mTargetMatchIds.SetValue(self.mcfndx, inTgtId)
      gndx = self.mcfndx * 3
      self.mThisProcXyzs.SetValue(gndx, inSrcXyz[0])
      self.mThisProcXyzs.SetValue(gndx+1, inSrcXyz[1])
      self.mThisProcXyzs.SetValue(gndx+2, inSrcXyz[2])
      self.mTargetMatchXyzs.SetValue(gndx, inTgtXyz[0])
      self.mTargetMatchXyzs.SetValue(gndx+1, inTgtXyz[1])
      self.mTargetMatchXyzs.SetValue(gndx+2, inTgtXyz[2])
      #if PhactoriDbg():
      #  myDebugPrint3("closer point found put in index " + \
      #      str(self.mcfndx) + \
      #      "  dstsqrd: " + str(self.mDistSqrds.GetValue(self.mcfndx)) + "\n" + \
      #      "\nsource id xyz: " + \
      #      str(self.mThisProcIds.GetValue(self.mcfndx)) + "   " + \
      #      str(self.mThisProcXyzs.GetValue(gndx)) + ", " + \
      #      str(self.mThisProcXyzs.GetValue(gndx+1)) + ", " + \
      #      str(self.mThisProcXyzs.GetValue(gndx+2)) + ", " + \
      #      "\ntarget id xyz: " + \
      #      str(self.mTargetMatchIds.GetValue(self.mcfndx)) + "   " + \
      #      str(self.mTargetMatchXyzs.GetValue(gndx)) + ", " + \
      #      str(self.mTargetMatchXyzs.GetValue(gndx+1)) + ", " + \
      #      str(self.mTargetMatchXyzs.GetValue(gndx+2)) + "\n")

      #now find which in the current list has the biggest distance, as it is
      #next in line for replacement (we do this to avoid having to shift
      #elements every time
      self.mcfndx = 0
      self.mMinDistSqrd = self.mDistSqrds.GetValue(0)
      #if PhactoriDbg():
      #  myDebugPrint3("find next index to be replaced try 0 \n" + \
      #      str(self.mMinDistSqrd) + "\n")
      for kk in range(1, self.mNumToFind):
        #if PhactoriDbg():
        #  myDebugPrint3("try " + str(kk) + " " + \
        #      str(self.mDistSqrds.GetValue(kk)) + " >? " + \
        #      str(self.mMinDistSqrd)+ "\n")
        if self.mDistSqrds.GetValue(kk) > self.mMinDistSqrd:
          self.mcfndx = kk
          self.mMinDistSqrd = self.mDistSqrds.GetValue(kk)
          #if PhactoriDbg():
          #  myDebugPrint3("yes, now " + str(self.mcfndx) + " " + \
          #      str(self.mMinDistSqrd) + "\n")
      #if PhactoriDbg():
      #  myDebugPrint3("next to be replaced ndx: " + str(self.mcfndx) + \
      #      " sid: " + \
      #      str(self.mThisProcIds.GetValue(self.mcfndx)) + \
      #      " tid: " + \
      #      str(self.mTargetMatchIds.GetValue(self.mcfndx)) + \
      #      " dsqrd: " + \
      #      str(self.mDistSqrds.GetValue(self.mcfndx)) + "\n")

    def TestPointSub1(self, inId, inX, inY, inZ, tId, tX, tY, tZ):
      """given an xyz point in the local processor (and its global node id)
         see if it is closer to one of the target points than the current
         set of nearest points and, if so, put it in the set, dropping
         others out if necessary"""
      ddx = inX - tX
      ddy = inY - tY
      ddz = inZ - tZ
      dstsqd = ddx*ddx + ddy*ddy + ddz*ddz
      if dstsqd >= self.mMinDistSqrd:
        #we already have mNumToFind closer
        #if PhactoriDbg():
        #  myDebugPrint3("inId tId " + str(inId) + " " + str(tId) + \
        #      " too far: " + str(dstsqd) + " >= " + \
        #      str(self.mMinDistSqrd) + "\n")
        return

      #replace the previous point that was furthest
      self.mDistSqrds.SetValue(self.mcfndx, dstsqd)
      self.mThisProcIds.SetValue(self.mcfndx, inId)
      self.mTargetMatchIds.SetValue(self.mcfndx, tId)
      gndx = self.mcfndx * 3
      self.mThisProcXyzs.SetValue(gndx, inX)
      self.mThisProcXyzs.SetValue(gndx+1, inY)
      self.mThisProcXyzs.SetValue(gndx+2, inZ)
      self.mTargetMatchXyzs.SetValue(gndx, tX)
      self.mTargetMatchXyzs.SetValue(gndx+1, tY)
      self.mTargetMatchXyzs.SetValue(gndx+2, tZ)
      if PhactoriDbg():
        myDebugPrint3("closer point found put in index " + \
            str(self.mcfndx) + \
            "  dstsqrd: " + str(self.mDistSqrds.GetValue(self.mcfndx)) + "\n" + \
            "\nsource id xyz: " + \
            str(self.mThisProcIds.GetValue(self.mcfndx)) + "   " + \
            str(self.mThisProcXyzs.GetValue(gndx)) + ", " + \
            str(self.mThisProcXyzs.GetValue(gndx+1)) + ", " + \
            str(self.mThisProcXyzs.GetValue(gndx+2)) + ", " + \
            "\ntarget id xyz: " + \
            str(self.mTargetMatchIds.GetValue(self.mcfndx)) + "   " + \
            str(self.mTargetMatchXyzs.GetValue(gndx)) + ", " + \
            str(self.mTargetMatchXyzs.GetValue(gndx+1)) + ", " + \
            str(self.mTargetMatchXyzs.GetValue(gndx+2)) + "\n")

      #now find which in the current list has the biggest distance, as it is
      #next in line for replacement (we do this to avoid having to shift
      #elements every time
      self.mcfndx = 0
      self.mMinDistSqrd = self.mDistSqrds.GetValue(0)
      if PhactoriDbg():
        myDebugPrint3("find next index to be replaced try 0 \n" + \
            str(self.mMinDistSqrd) + "\n")
      for kk in range(1, self.mNumToFind):
        if PhactoriDbg():
          myDebugPrint3("try " + str(kk) + " " + \
              str(self.mDistSqrds.GetValue(kk)) + " >? " + \
              str(self.mMinDistSqrd)+ "\n")
        if self.mDistSqrds.GetValue(kk) > self.mMinDistSqrd:
          self.mcfndx = kk
          self.mMinDistSqrd = self.mDistSqrds.GetValue(kk)
          if PhactoriDbg():
            myDebugPrint3("yes, now " + str(self.mcfndx) + " " + \
                str(self.mMinDistSqrd) + "\n")
      if PhactoriDbg():
        myDebugPrint3("next to be replaced ndx: " + str(self.mcfndx) + \
            " sid: " + \
            str(self.mThisProcIds.GetValue(self.mcfndx)) + \
            " tid: " + \
            str(self.mTargetMatchIds.GetValue(self.mcfndx)) + \
            " dsqrd: " + \
            str(self.mDistSqrds.GetValue(self.mcfndx)) + "\n")

    def ToStr(self):
      retStr = "closest " + str(self.mNumToFind) + " points:\n" + \
      "index: source id, target id, dist sqrd: source xyz: target xyz\n"
      for ii in range(0, self.mNumToFind):
        pp = ii*3
        adst = str(ii) + ": " + \
               str(self.mThisProcIds.GetValue(ii)) + ", " + \
               str(self.mTargetMatchIds.GetValue(ii)) + ", " + \
               str(self.mDistSqrds.GetValue(ii)) + ": " + \
               str(self.mThisProcXyzs.GetValue(pp)) + ", " + \
               str(self.mThisProcXyzs.GetValue(pp+1)) + ", " + \
               str(self.mThisProcXyzs.GetValue(pp+2)) + ": " + \
               str(self.mTargetMatchXyzs.GetValue(pp)) + ", " + \
               str(self.mTargetMatchXyzs.GetValue(pp+1)) + ", " + \
               str(self.mTargetMatchXyzs.GetValue(pp+2)) + "\n"
        retStr += adst
      return retStr


  def FindClosestNPointsToList(self, inGlobalNodeIdList, inPointXyzList,
            inNumToFind):
    """given a list of node ids and xyz points, recursively find the nearest
       (geometrically) inNumToFind points in the local processor to the xyz
       points in inPointXyzList.  Returns a list of inNumToFind global node
       ids and inNumToFind xyz points
       inGlobalNodeIds is vtkIntArray
       inPointXyzs is vtkDoubleArray
       returns a FindClosestNPointsToListParams instance which has list of
       length inNumToFind which contain the node id of the closest local
       process point, the node id of the corresponding point from the
       inGlobalNodeIdList, the xyzs of each of those, and this distance
       between (squared to save computation)"""
    recursionItem = BlockRecursionControlItem()
    recursionItem.mParameters = \
        PhactoriOperationBlock.FindClosestNPointsToListParams(
            self, inNumToFind, inGlobalNodeIdList, inPointXyzList)
    recursionItem.mOperationToDoPerBlock = \
            self.FindClosestNPointsToListInBlock
    self.DoMethodPerBlock(recursionItem)
    return recursionItem.mParameters

  def FindClosestNPointsToListInBlock(self, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("FindClosestNPointsToListInBlock entered\n")
    pointsArray = inInputCsData.GetPoints()
    if pointsArray == None:
      #no points here
      return
    pointsData = inInputCsData.GetPointData()
    globalNodeIdArray = pointsData.GetArray('GlobalNodeId')
    numPoints = pointsArray.GetNumberOfPoints()
    ptXyz = [0.0, 0.0, 0.0]

    #get target to set up Kdtree for quickly finding nearest point
    inParameters.SetUpKdtree()

    for ndx in range(0,numPoints):
      #thePoint = pointsArray.GetPoint(ndx, ptXyz)
      if ndx % 100 == 0:
        if PhactoriDbg():
          myDebugPrint3("test " + str(ndx) + " of " + str(numPoints) + "\n")
      pointsArray.GetPoint(ndx, ptXyz)
      if(globalNodeIdArray == None):
        theGlobalNodeId = ndx + 1
      else:
        theGlobalNodeId = globalNodeIdArray.GetValue(ndx)

      inParameters.TestPointWithKdtree(theGlobalNodeId, ptXyz)
      #inParameters.TestPoint(theGlobalNodeId, ptXyz)
    if PhactoriDbg(100):
      myDebugPrint3("FindClosestNPointsToListInBlock returning\n")

  class MakeListOfAllPoints1Params:
    def __init__(self):
      self.mGlobalNodeIdList = vtk.vtkIntArray()
      self.mPointXYZList = vtk.vtkDoubleArray()

  def MakeListOfAllPoints1(self):
    """recursively going through multiblock setup, make a list of all the
       points in this operation output in this process.  We get a list
       of global node ids and a list of xyz geometry points"""
    recursionItem = BlockRecursionControlItem()
    recursionItem.mParameters = \
        PhactoriOperationBlock.MakeListOfAllPoints1Params()
    recursionItem.mOperationToDoPerBlock = \
            self.MakeListOfAllPointsInBlock1
    self.DoMethodPerBlock(recursionItem)
    return recursionItem.mParameters.mGlobalNodeIdList, \
           recursionItem.mParameters.mPointXYZList

  def MakeListOfAllPointsInBlock1(self, inInputCsData, inParameters):
    #if PhactoriDbg(100):
    #  myDebugPrint3("MakeListOfAllPointsInBlock1 entered\n")
    pointsArray = inInputCsData.GetPoints()
    if pointsArray == None:
      #no points here
      return
    pointsData = inInputCsData.GetPointData()
    globalNodeIdArray = pointsData.GetArray('GlobalNodeId')
    numPoints = pointsArray.GetNumberOfPoints()
    ptXyz = [0.0, 0.0, 0.0]
    for ndx in range(0,pointsArray.GetNumberOfPoints()):
      #thePoint = pointsArray.GetPoint(ndx, ptXyz)
      pointsArray.GetPoint(ndx, ptXyz)
      if globalNodeIdArray == None:
        inParameters.mGlobalNodeIdList.InsertNextValue(ndx+1)
      else:
        inParameters.mGlobalNodeIdList.InsertNextValue(globalNodeIdArray.GetValue(ndx))
      inParameters.mPointXYZList.InsertNextValue(ptXyz[0])
      inParameters.mPointXYZList.InsertNextValue(ptXyz[1])
      inParameters.mPointXYZList.InsertNextValue(ptXyz[2])
    #if PhactoriDbg(100):
    #  myDebugPrint3("MakeListOfAllPointsInBlock1 returning\n")

  class MakeListOfAllPointsAndNodeIdsInBlock2Params:
    def __init__(self):
      self.mPointXYZAndNodeIdList = []

  def MakeListOfAllPointsAndNodeIdsInBlock2(self, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("MakeListOfAllPointsAndNodeIdsInBlock2 entered\n")
    pointsArray = inInputCsData.GetPoints()
    if pointsArray == None:
      #no points here
      return
    pointsData = inInputCsData.GetPointData()
    globalNodeIdArray = pointsData.GetArray('GlobalNodeId')
    numPoints = pointsArray.GetNumberOfPoints()
    ptXyz = [0.0, 0.0, 0.0]
    for ndx in range(0,pointsArray.GetNumberOfPoints()):
      #thePoint = pointsArray.GetPoint(ndx, ptXyz)
      pointsArray.GetPoint(ndx, ptXyz)
      if globalNodeIdArray == None:
        thisNodeId = -1
      else:
        thisNodeId = int(globalNodeIdArray.GetValue(ndx))
      inParameters.mPointXYZAndNodeIdList.append([ptXyz[0], ptXyz[1], ptXyz[2], thisNodeId])
    if PhactoriDbg(100):
      numPtsThisBlock = len(inParameters.mPointXYZAndNodeIdList)
      myDebugPrint3("num points after this block: " + str(numPtsThisBlock) + "\n")
      if numPtsThisBlock > 0:
        myDebugPrint3("last point: " + str(inParameters.mPointXYZAndNodeIdList[-1]) + "\n")
    if PhactoriDbg(100):
      myDebugPrint3("MakeListOfAllPointsAndNodeIdsInBlock2 returning\n")

  def MakeListOfAllPointsAndNodeIdsOnThisProcessFromParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("MakeListOfAllPointsAndNodeIdsOnThisProcessFromParaViewFilter entered\n")
    recursionItem = BlockRecursionControlItem()
    recursionItem.mParameters = PhactoriOperationBlock.MakeListOfAllPointsAndNodeIdsInBlock2Params()
    recursionItem.mOperationToDoPerBlock = self.MakeListOfAllPointsAndNodeIdsInBlock2
    self.DoMethodPerBlockFromParaViewFilter(recursionItem, inInputFilter)
    if PhactoriDbg(100):
      numPtsThisProcess = len(recursionItem.mParameters.mPointXYZAndNodeIdList)
      myDebugPrint3("num points in this process: " + str(numPtsThisProcess) + "\n")
      if numPtsThisProcess > 0:
        myDebugPrint3("last point: " + str(recursionItem.mParameters.mPointXYZAndNodeIdList[-1]) + "\n")
    if PhactoriDbg(100):
      myDebugPrint3("MakeListOfAllPointsAndNodeIdsOnThisProcessFromParaViewFilter returning\n")
    return recursionItem.mParameters.mPointXYZAndNodeIdList

  class GetQuadClosestToPointsInBlock2Params:
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
      if PhactoriDbg(100):
        myDebugPrint3("InitializeWithPointList: " + str(numTestPoints) + " " +\
          str(len(self.distSqrdList)) + " " + str(len(self.closestList)) +\
          "\n" + str(self.testPointList) + "\n" +\
          "\n" + str(self.distSqrdList) + "\n" +\
          "\n" + str(self.closestList) + "\n")

  def GetQuadsClosestToPointsInBlock2(self, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("GetCloseQuadsInBlock2 entered\n")
    numQuads = inInputCsData.GetNumberOfCells()
    if numQuads == 0:
      #no quads here
      return

    if PhactoriDbg(100):
      myDebugPrint3(str(inParameters.testPointList) + "\n")
      myDebugPrint3(str(inParameters.distSqrdList) + "\n")

    cellsData = inInputCsData.GetCellData()
    globalCellIdArray = cellsData.GetArray('GlobalCellId')
    for quadIndex in range(0,numQuads):
      thisQuad = QuadInfoGQC()
      thisQuad.PopulateFromVtkData(inInputCsData, globalCellIdArray, quadIndex)
      for ptndx, oneTestPt in enumerate(inParameters.testPointList):
        testDist = thisQuad.CalculateDistanceSquaredTo(oneTestPt)
        if testDist < inParameters.distSqrdList[ptndx]:
          inParameters.closestList[ptndx] = thisQuad
          inParameters.distSqrdList[ptndx] = testDist

    if PhactoriDbg(100):
      myDebugPrint3(str(inParameters.testPointList) + "\n")
      myDebugPrint3(str(inParameters.distSqrdList) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("after this block:\n")
      for ii, oneQuad in enumerate(inParameters.closestList):
        myDebugPrint3(str(ii) + ": " + \
          str(inParameters.distSqrdList[ii]) + "\n" + \
          str(inParameters.testPointList[ii]) + "\n" + oneQuad.ToString())
      myDebugPrint3("\n")

    if PhactoriDbg(100):
      myDebugPrint3("GetCloseQuadsInBlock2 returning\n")

  def GetQuadsClosestToPointsOnThisProcessFromParaViewFilter(self, inInputFilter, inTestPointList):
    if PhactoriDbg(100):
      myDebugPrint3("GetQuadsClosestToPointsOnThisProcessFromParaViewFilter entered\n")
    recursionItem = BlockRecursionControlItem()
    recursionItem.mParameters = PhactoriOperationBlock.GetQuadClosestToPointsInBlock2Params()
    recursionItem.mParameters.InitializeWithPointList(inTestPointList)
    recursionItem.mOperationToDoPerBlock = self.GetQuadsClosestToPointsInBlock2
    self.DoMethodPerBlockFromParaViewFilter(recursionItem, inInputFilter)
    if PhactoriDbg(100):
      myDebugPrint3("GetQuadsClosestToPointsOnThisProcessFromParaViewFilter returning\n")
    return recursionItem.mParameters.closestList, recursionItem.mParameters.distSqrdList

  class OutputElementListToFileParams:
    def __init__(self):
      self.mOutFileff = None
      self.mBlockCount = 0
      self.mElementCount = 0
      self.mNodeCount = 0
      self.mFlag1 = 0
      self.mKilledByCriteriaCount = [0,0,0,0,0,0,0,0,0,0]

  def OutputElementListToFile(self, inFileNameToWrite):
    recursionItem = BlockRecursionControlItem()

    recursionItem.mParameters = \
      PhactoriOperationBlock.OutputElementListToFileParams()

    if WriteEachDeadCellElementToFiles:
      recursionItem.mParameters.mOutFileff = open(inFileNameToWrite, "w+b")
    recursionItem.mOperationToDoPerBlock = \
            self.OutputElementListFromOneBlockToFile
    self.DoMethodPerBlock(recursionItem)
    if WriteEachDeadCellElementToFiles:
      recursionItem.mParameters.mOutFileff.close()
    return recursionItem.mParameters

  def OutputSingleElementFromBlockToTimeHistoryFile(
          self, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("OutputSingleElementFromBlockToTimeHistoryFile entered\n")
    cellData = inInputCsData.GetCellData()
    numCellTuplesX = cellData.GetNumberOfTuples()
    numCellArrays = cellData.GetNumberOfArrays()
    pointData = inInputCsData.GetPointData()
    numPointTuplesX = cellData.GetNumberOfTuples()
    numPointArrays = pointData.GetNumberOfArrays()

    #if inParameters.mBlockCount == 0:
    if (inParameters.mFlag1 == 0) and (numCellTuplesX > 0):
      inParameters.mFlag1 = 1
      #as constructed, we should only open once, so we can do 'w' not 'a+'
      inParameters.mOutElementFileff = open(inParameters.mOutElementFileName, "w")
      inParameters.mOutElementFileff.write("step, simtime")
      #inParameters.mOutElementFileff.write(", element index")
      for arrayIndex in range(0, numCellArrays):
          inParameters.mOutElementFileff.write(
            ", " + str(cellData.GetArray(arrayIndex).GetName()))
      inParameters.mOutElementFileff.write("\n")

      inParameters.mOutNodeFileff = open(inParameters.mOutNodeFileName, "w")
      inParameters.mOutNodeFileff.write("step, simtime")
      for arrayIndex in range(0, numPointArrays):
          oneArray = pointData.GetArray(arrayIndex)
          numComp = oneArray.GetNumberOfComponents()
          oneArrayName = str(oneArray.GetName())
          if numComp == 1:
            inParameters.mOutNodeFileff.write(", " + oneArrayName)
          elif(numComp == 3):
            inParameters.mOutNodeFileff.write(
                    ", " + oneArrayName + "x" + \
                    ", " + oneArrayName + "y" + \
                    ", " + oneArrayName + "z")
          else:
            myDebugPrint3AndException(
              "OutputSingleElementFromBlockToTimeHistoryFile:\n" \
              "expecting 1 or 3 components (A)\n")
      inParameters.mOutNodeFileff.write("\n")

    #for now, only dump 1 element
    if numCellTuplesX > 0:
        numCellsToDo = 1
    else:
        numCellsToDo = 0
    for cellIndex in range(0, numCellsToDo):
        inParameters.mOutElementFileff.write(str(inParameters.mTimeStep) + ", " + \
                str(inParameters.mSimulationTime))
        #inParameters.mOutElementFileff.write(", " + \
        #        str(cellIndex + inParameters.mElementCount))
        for arrayIndex in range(0, numCellArrays):
            inParameters.mOutElementFileff.write(
              ", " + str(cellData.GetArray(arrayIndex).GetTuple1(cellIndex)))
        inParameters.mOutElementFileff.write("\n")

    #for now, only dump 1 node
    if numPointTuplesX > 0:
        numPointsToDo = 1
    else:
        numPointsToDo = 0
    for pointIndex in range(0, numPointsToDo):
        inParameters.mOutNodeFileff.write(str(inParameters.mTimeStep) + ", " + \
                str(inParameters.mSimulationTime))
        #inParameters.mOutNodeFileff.write(", " + \
        #        str(pointIndex + inParameters.mCellCount))
        for arrayIndex in range(0, numPointArrays):
            #inParameters.mOutNodeFileff.write(
            #  ", " + str(pointData.GetArray(arrayIndex).GetTuple1(0)))
            oneArray = pointData.GetArray(arrayIndex)
            numComp = oneArray.GetNumberOfComponents()
            if(numComp == 1):
              inParameters.mOutNodeFileff.write(
                ", " + str(pointData.GetArray(arrayIndex).GetTuple1(pointIndex)))
            elif(numComp == 3):
              theTuple = pointData.GetArray(arrayIndex).GetTuple3(0)
              inParameters.mOutNodeFileff.write( \
                ", " + str(theTuple[0]) + \
                ", " + str(theTuple[1]) + \
                ", " + str(theTuple[2]))
            else:
              myDebugPrint3AndException(
                "OutputSingleElementFromBlockToTimeHistoryFile:\n" \
                "expecting 1 or 3 components (B)\n")
        inParameters.mOutNodeFileff.write("\n")

    inParameters.mNodeCount += inInputCsData.GetNumberOfPoints()
    inParameters.mElementCount += inInputCsData.GetNumberOfCells()
    inParameters.mBlockCount += 1

    if PhactoriDbg(100):
      numCells = inInputCsData.GetNumberOfCells()
      if numCells > 0:
        numPoints = inInputCsData.GetNumberOfPoints()
        myDebugPrint3(
          "numCells: " + str(numCells) + "\n" \
          "numCellTuplesX: " + str(numCellTuplesX) + "\n" \
          "numPoints: " + str(numPoints) + "\n" \
          "numPointTuplesX: " + str(numPointTuplesX) + "\n")

    if(inParameters.mOutElementFileff != None):
      inParameters.mOutElementFileff.flush()
      inParameters.mOutNodeFileff.flush()

    #dump 1 node

    if PhactoriDbg(100):
      myDebugPrint3("OutputSingleElementFromBlockToTimeHistoryFile returning\n")

  class OutputSingleElementToTimeHistoryFileParams:
    def __init__(self):
      self.mOutElementFileff = None
      self.mOutElementFileName = None
      self.mOutNodeFileff = None
      self.mOutNodeFileName = None
      self.mTimeStep = 0
      self.mSimulationTime = 0.0
      self.mBlockCount = 0
      self.mElementCount = 0
      self.mNodeCount = 0
      self.mFlag1 = 0

    def ResetBlockNodeElementCount(self):
      self.mBlockCount = 0
      self.mElementCount = 0
      self.mNodeCount = 0

  def OutputSingleElementToTimeHistoryFile(self,
          inElementFileNameToWrite, inNodeFileNameToWrite,
          inTimeStep, inSimulationTime):
    recursionItem = BlockRecursionControlItem()

    if inElementFileNameToWrite in self.mRecursionParameterStore:
      recursionItem.mParameters = \
        self.mRecursionParameterStore[inElementFileNameToWrite]
      recursionItem.mParameters.ResetBlockNodeElementCount()
    else:
      recursionItem.mParameters = \
        PhactoriOperationBlock.OutputSingleElementToTimeHistoryFileParams()
      self.mRecursionParameterStore[inElementFileNameToWrite] = \
        recursionItem.mParameters

    #recursionItem.mParameters.mOutFileff = open(inFileNameToWrite, "w+b")
    recursionItem.mParameters.mOutElementFileName = inElementFileNameToWrite
    recursionItem.mParameters.mOutNodeFileName = inNodeFileNameToWrite
    recursionItem.mParameters.mTimeStep = inTimeStep
    recursionItem.mParameters.mSimulationTime = inSimulationTime
    recursionItem.mOperationToDoPerBlock = \
            self.OutputSingleElementFromBlockToTimeHistoryFile
    self.DoMethodPerBlock(recursionItem)
    #recursionItem.mParameters.mOutFileff.close()

  def ExportOperationData(self, datadescription):
    """this will be called once per callback (before WriteImages) to allow the
       operation to export any desired data which is not an image.  We call
       the operation specifics version of this method."""
    self.mOperationSpecifics.ExportOperationData(datadescription)

#phactori_combine_to_single_python_file_subpiece_end_1
