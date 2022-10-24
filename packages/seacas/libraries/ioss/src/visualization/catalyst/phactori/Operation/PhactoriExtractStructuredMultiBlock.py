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
from paraview.simple import *

#phactori_combine_to_single_python_file_subpiece_begin_1

global gStructuredGridFound
gStructuredGridFound = None
global gStructuredGridFoundExtent
gStructuredGridFoundExtent = None

def FigureBlockIndicesFromBlockListOneBlock(includeIndexList, includeBlockList,
  inMetaData, ioFlatIndexCounter, inCsdata, inForceSetting):
  """determine if this one block should have it's flat index tripped on for
     the extract block filter (leaf item of recursion)"""
  if PhactoriDbg(100):
    myDebugPrint3("FigureBlockIndicesFromBlockListOneBlock entered\n"
      "ioFlatIndexCounter " + str(ioFlatIndexCounter,) + " inForceSetting " + str(inForceSetting) + "\n"
      "2 inMetaData: " + str(inMetaData) + "\n")

  if inMetaData == None:
    thisBlockName = None
  else:
    thisBlockName = inMetaData.Get(vtk.vtkCompositeDataSet.NAME())

  if (thisBlockName == None) and (inForceSetting != 1):
    if PhactoriDbg(100):
      myDebugPrint3("block with no name " + \
          " not in include list, not + to mBlockIndices (flat index " + \
          str(ioFlatIndexCounter[0] - 1) + ")\n")
  elif (inForceSetting == 1) or (thisBlockName in includeBlockList):
    includeIndexList.append(int(ioFlatIndexCounter[0]) - 1)
    blockClassName = inCsdata.GetClassName()
    if blockClassName == "vtkStructuredGrid":
      global gStructuredGridFound
      gStructuredGridFound = inCsdata
      global gStructuredGridFoundExtent
      gStructuredGridFoundExtent = inMetaData.Get(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT())
      if PhactoriDbg(100):
        myDebugPrint3("this leaf is structured grid: " + str(thisBlockName) + "\n"
          "vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(): " + str(gStructuredGridFoundExtent) + "\n"
          "A inCsdata.GetExtent(): " + str(inCsdata.GetExtent()) + "\n")
    #global gDuplicateNameCounter
    #if thisBlockName in gDuplicateNameCounter:
    #  oldCount = gDuplicateNameCounter[thisBlockName]
    #  gDuplicateNameCounter[thisBlockName] = oldCount+1
    #else:
    #  gDuplicateNameCounter[thisBlockName] = 1
    if PhactoriDbg(100):
      myDebugPrint3("block " + str(thisBlockName) + \
          " in include list, + to mBlockIndices (flat index " + \
          str(ioFlatIndexCounter[0] - 1) + ")\n")
  else:
    if PhactoriDbg(100):
      myDebugPrint3("block " + str(thisBlockName) + \
          " not in include list, not + to mBlockIndices (flat index " + \
          str(ioFlatIndexCounter[0] - 1) + ")\n")

  if PhactoriDbg(100):
    myDebugPrint3("FigureBlockIndicesFromBlockListOneBlock returning\n")

def FigureBlockIndicesFromBlockListRecurse1(includeIndexList, includeBlockList,
  inCsdata, inMetaData, ioFlatIndexCounter, inForceSetting):
  """recursively go through multiblock dataset to determine flat indices
     of blocks to be included; we are assuming only leaf blocks are
     actually named"""

  icsdClassname = inCsdata.GetClassName()

  if PhactoriDbg(100):
    if inMetaData == None:
      thisBlockName = None
    else:
      thisBlockName = inMetaData.Get(vtk.vtkCompositeDataSet.NAME())
    myDebugPrint3("FigureBlockIndicesFromBlockListRecurse1 " + \
      str(thisBlockName) + " index: " + str(int(ioFlatIndexCounter[0])) + \
      " inForceSetting: " + str(inForceSetting) + " classname: " + \
      str(icsdClassname) + "\n")

  ioFlatIndexCounter[0] += 1

  if icsdClassname == "vtkMultiBlockDataSet":
    thisIsLeaf = False
    nonLeafType = 0
  elif icsdClassname == "vtkExodusIIMultiBlockDataSet":
    thisIsLeaf = False
    nonLeafType = 1
  elif icsdClassname == "vtkMultiPieceDataSet":
    thisIsLeaf = False
    nonLeafType = 2
  elif icsdClassname == "vtkPartitionedDataSet":
    thisIsLeaf = False
    nonLeafType = 3
  else:
    thisIsLeaf = True
    nonLeafType = -1

  if not thisIsLeaf:
    if inForceSetting == 1:
      includeIndexList.append(int(ioFlatIndexCounter[0]) - 1)
      includeOrExcludeAllSubBlocks = 1
      if PhactoriDbg(100):
        myDebugPrint3("this non-leaf block added inForceSetting == 1\n")
    elif inForceSetting == -1:
      includeOrExcludeAllSubBlocks = -1
      if PhactoriDbg(100):
        myDebugPrint3("this non-leaf block not added inForceSetting == -1\n")
    else:
      includeOrExcludeAllSubBlocks = 0
      #this is a non-leaf node, but we want to add it (or not) depending on
      #the filter settings
      if inMetaData == None:
        thisBlockName = None
      else:
        thisBlockName = inMetaData.Get(vtk.vtkCompositeDataSet.NAME())

      if (thisBlockName != None) and (thisBlockName in includeBlockList):
        if PhactoriDbg(100):
          myDebugPrint3("include list, append nonleaf " + str(thisBlockName) + \
            " ndx " + str(int(ioFlatIndexCounter[0]) - 1) + "\n")
        includeIndexList.append(int(ioFlatIndexCounter[0]) - 1)
        includeOrExcludeAllSubBlocks = 1
        if PhactoriDbg(100):
          myDebugPrint3("includeOrExcludeAllSubBlocks set to 1\n")
      else:
        if PhactoriDbg(100):
          myDebugPrint3("include list, exclude nonleaf " + str(thisBlockName) + \
            " ndx " + str(int(ioFlatIndexCounter[0]) - 1) + "\n")

    if nonLeafType == 2:
      numBlocks = inCsdata.GetNumberOfPieces()
    elif nonLeafType == 3:
      numBlocks = inCsdata.GetNumberOfPartitions()
    else:
      numBlocks = inCsdata.GetNumberOfBlocks()

    if PhactoriDbg(100):
      myDebugPrint3("recursing, flat index: " + \
          str(ioFlatIndexCounter[0]) + \
          "    num blocks: " + \
          str(numBlocks) + "\n")
    for ii in range(0, numBlocks):
      if nonLeafType == 2:
        oneBlock = inCsdata.GetPiece(ii)
      elif nonLeafType == 3:
        oneBlock = inCsdata.GetPartition(ii)
      else:
        oneBlock = inCsdata.GetBlock(ii)
      oneBlockMetaData = inCsdata.GetMetaData(ii)
      if PhactoriDbg(100):
        myDebugPrint3("oneBlockMetaData: " + str(oneBlockMetaData) + "\n")
        #myDebugPrint3("name: " + \
        #    oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME()) + "\n")
      if oneBlock != None:
        if oneBlockMetaData != None:
          theBlockName = oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME())
          if PhactoriDbg(100):
            myDebugPrint3("name for block " + str(ii) + ": " + \
              str(theBlockName) + "\n")
        else:
          if PhactoriDbg(100):
            myDebugPrint3("block " + str(ii) + " meta data was None\n")
        FigureBlockIndicesFromBlockListRecurse1(includeIndexList, includeBlockList,
          oneBlock, oneBlockMetaData, ioFlatIndexCounter,
          includeOrExcludeAllSubBlocks)
      else:
        #I think we need to count here to be okay with pruned stuff; maybe
        #we need to set extract block to no pruning (?)
        ioFlatIndexCounter[0] += 1
        if includeOrExcludeAllSubBlocks == 1:
          if PhactoriDbg(100):
            myDebugPrint3("leaf block None index " + \
            str(int(ioFlatIndexCounter[0]) - 1) + \
            " appended due to includeOrExcludeAllSubBlocks == 1\n\n")
          includeIndexList.append(int(ioFlatIndexCounter[0]) - 1)
        else:
          if PhactoriDbg(100):
            myDebugPrint3("force include is not on, so check include list\n")
          includeThisBlock = False
          if (includeBlockList != None) and (oneBlockMetaData != None):
            thisBlockName = oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME())
            if PhactoriDbg(100):
              myDebugPrint3("include list not None, thisBlockName " + str(thisBlockName) + "\n"
                "includeBlockList " + str(includeBlockList) + "\n")
            if (thisBlockName != None) and (thisBlockName in includeBlockList):
              includeThisBlock = True
              if PhactoriDbg(100):
                myDebugPrint3("leaf block None index " + \
                str(int(ioFlatIndexCounter[0]) - 1) + \
                " included due to name being in includeBlockList\n")
              includeIndexList.append(int(ioFlatIndexCounter[0]) - 1)
          if includeThisBlock == False:
            if PhactoriDbg(100):
              myDebugPrint3("leaf block None index " + \
              str(int(ioFlatIndexCounter[0]) - 1) + \
              " excluded due to includeOrExcludeAllSubBlocks != 1\n")
  else:
    if inForceSetting != -1:
      if inForceSetting == 1:
        if PhactoriDbg(100):
          myDebugPrint3("this leaf block will be added inForceSetting == 1\n")
      FigureBlockIndicesFromBlockListOneBlock(includeIndexList,
        includeBlockList, inMetaData, ioFlatIndexCounter, inCsdata, inForceSetting)
    else:
      if PhactoriDbg(100):
        myDebugPrint3("this leaf block not added inForceSetting == -1\n")
  #if icsdClassname == "vtkMultiPieceDataSet":
  #  numpieces = inCsdata.GetNumberOfPieces()
  #  # ioFlatIndexCounter[0] += numpieces - 1
  #  ioFlatIndexCounter[0] += numpieces


def FigureBlockIndicesFromBlockList(includeBlockList, inInputFilter):
  """from the includeBlockList create a indices list
     of blocks to put in filter"""
  if PhactoriDbg(100):
    myDebugPrint3("FigureBlockIndicesFromBlockList entered"
        "\nincludeBlockList: " + str(includeBlockList) + "\n")

  #global gDuplicateNameCounter
  #gDuplicateNameCounter = {}
  csdata = inInputFilter.GetClientSideObject().GetOutputDataObject(0)
  flatIndexCounter = [0]
  listOfBlockIndicesToInclude = []
  FigureBlockIndicesFromBlockListRecurse1(listOfBlockIndicesToInclude,
         includeBlockList, csdata, None, flatIndexCounter, 0)

  #if PhactoriDbg(100):
  #  myDebugPrint3("number of times existing block names found:\n")
  #  for blknm, count in gDuplicateNameCounter.items():
  #    myDebugPrint3(blknm + ": " + str(count) + "\n");
  if PhactoriDbg(100):
    myDebugPrint3("FigureBlockIndicesFromBlockList final indices list:\n" \
        + str(listOfBlockIndicesToInclude) + \
        "\nFigureBlockIndicesFromBlockList returning\n")

  return listOfBlockIndicesToInclude

def PrintMultiBlockStatsOneBlock(inCsdata, inMetaData, ioFlatIndexCounter):
  icsdClassname = inCsdata.GetClassName()
  myDebugPrint3("leaf block " + str(icsdClassname) + " index " + str(ioFlatIndexCounter[0]) + "\n")
  if icsdClassname == "vtkStructuredGrid":
    wholeExtent = inMetaData.Get(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT())
    myDebugPrint3("this leaf is structured grid: " + str(ioFlatIndexCounter[0]) + "\n"
      "vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(): " + str(wholeExtent) + "\n"
      "B inCsdata.GetExtent(): " + str(inCsdata.GetExtent()) + "\n")

def PrintMultiBlockStatsRecurse1(inCsdata, inMetaData, recursionLevel, ioFlatIndexCounter):
  myDebugPrint3("PrintMultiBlockStatsRecurse1 entered recursionLevel " + str(recursionLevel) + "\n")

  icsdClassname = inCsdata.GetClassName()
  myDebugPrint3("inCsdata.GetClassName() " + str(icsdClassname) + " ioFlatIndexCounter[0] " + str(ioFlatIndexCounter[0]) + "\n")

  ioFlatIndexCounter[0] += 1

  if icsdClassname == "vtkMultiBlockDataSet":
    thisIsLeaf = False
    nonLeafType = 0
  elif icsdClassname == "vtkExodusIIMultiBlockDataSet":
    thisIsLeaf = False
    nonLeafType = 1
  elif icsdClassname == "vtkMultiPieceDataSet":
    thisIsLeaf = False
    nonLeafType = 2
  elif icsdClassname == "vtkPartitionedDataSet":
    thisIsLeaf = False
    nonLeafType = 3
  else:
    thisIsLeaf = True
    nonLeafType = -1

  myDebugPrint3("thisIsLeaf " + str(thisIsLeaf) + " nonLeafType " + str(nonLeafType) + "\n")
  if not thisIsLeaf:
    if nonLeafType == 2:
      numBlocks = inCsdata.GetNumberOfPieces()
    elif nonLeafType == 3:
      numBlocks = inCsdata.GetNumberOfPartitions()
    else:
      numBlocks = inCsdata.GetNumberOfBlocks()

    myDebugPrint3("recursing, flat index: " + \
        str(ioFlatIndexCounter[0]) + \
        "    num blocks: " + \
        str(numBlocks) + "\n")
    for ii in range(0, numBlocks):
      myDebugPrint3("doing sub block " + str(ii) + " of " + str(numBlocks) + "\n")
      if nonLeafType == 2:
        oneBlock = inCsdata.GetPiece(ii)
      elif nonLeafType == 3:
        oneBlock = inCsdata.GetPartition(ii)
      else:
        oneBlock = inCsdata.GetBlock(ii)
      oneBlockMetaData = inCsdata.GetMetaData(ii)
      myDebugPrint3("oneBlockMetaData: " + str(oneBlockMetaData) + "\n")
      #myDebugPrint3("name: " + \
      #    oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME()) + "\n")
      if oneBlock != None:
        PrintMultiBlockStatsRecurse1(oneBlock, oneBlockMetaData, recursionLevel + 1, ioFlatIndexCounter)
      else:
        myDebugPrint3("this sub block is None, index " + str(ioFlatIndexCounter[0]) + "\n")
        ioFlatIndexCounter[0] += 1
  else:
    PrintMultiBlockStatsOneBlock(inCsdata, inMetaData, ioFlatIndexCounter)

  myDebugPrint3("PrintMultiBlockStatsRecurse1 returning recursionLevel " + str(recursionLevel) + "\n")

def RecursivelyPrintMultiBlockStats(inFilter):
  if PhactoriDbg(100) == False:
    return

  myDebugPrint3("RecursivelyPrintMultiBlockStats entered\n")
  csdata = inFilter.GetClientSideObject().GetOutputDataObject(0)
  flatIndexCounter = [0]
  recursionLevel = 0
  PrintMultiBlockStatsRecurse1(csdata, None, recursionLevel, flatIndexCounter)

  myDebugPrint3("RecursivelyPrintMultiBlockStats returning\n")

class PhactoriExtractStructuredMultiBlock(PhactoriOperationSpecifics):
  """manages PhactoriExtractStructuredMultiBlock filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.ControlJson = None
    self.SavedIncludeBlockIndexList = None
    self.ExtractBlockFilters = []
    self.ExtractSubsetFilters = []
    self.FinalGroupFilter = None
    self.mStructuredGridsIncluded = []
    self.mStructuredGridsIncludedExtent  = []

  def ParseParametersFromJson(self, inJson):
    #["blk-4":{"ijkrange":[imin,imax,jmin,jmax,kmin,kmax]]
    #["blk-5":{"ijkrange":[imin,imax,jmin,jmax,kmin,kmax],"ijkrangefunction":[-1,0,0,1,0,1]}]
    #["blk-6":{"ijkrange":[imin,imax,jmin,jmax,kmin,kmax],"ijkrangefunction":[-1,0,-1,0,0,1]}]
    if "ExtractSubsetSettingsByBlock" not in inJson:
      myDebugPrint3AndException("PhactoriExtractStructuredMultiBlock.ParseParametersFromJson\n"
        "missing key 'ExtractSubsetSettingsByBlock'\n")

    self.ControlJson = inJson["ExtractSubsetSettingsByBlock"]
    dummy = 0

  def MakeExtractBlockFilter(self, inInputFilter, blockName):
    if PhactoriDbg(100):
      myDebugPrint3("MakeExtractBlockFilter entered: " + str(blockName) + "\n")
    newPvExtractBlockFilter = ExtractBlock(inInputFilter)
    #newParaViewFilter.PruneOutput = 1
    #newParaViewFilter.MaintainStructure = 0
    #newParaViewFilter.MaintainStructure = 1
    global gStructuredGridFound
    gStructuredGridFound = None
    self.SavedIncludeBlockIndexList = FigureBlockIndicesFromBlockList([blockName], inInputFilter)
    self.mStructuredGridsIncluded.append(gStructuredGridFound)
    self.mStructuredGridsIncludedExtent.append(gStructuredGridFoundExtent)

    if PhactoriDbg(100):
      myDebugPrint3("self.SavedIncludeBlockIndexList:\n" + str(self.SavedIncludeBlockIndexList) + "\n")
    newPvExtractBlockFilter.BlockIndices = self.SavedIncludeBlockIndexList
    if PhactoriDbg(100):
      myDebugPrint3("MakeExtractBlockFilter returning: " + str(blockName) + "\n")
    return newPvExtractBlockFilter

  def CreateIjkRangeFromSettingsForThisBlock(self, ijkJson, blockExtent):

    rngFuncAllAbsolute = True
    if "ijkrangefunction" in ijkJson:
      rngFunc = ijkJson["ijkrangefunction"]
      for ijksetting in rngFunc:
        if ijksetting != 0:
          rngFuncAllAbsolute = False
          break
    else:
      #rngFuncAllAbsolute = True
      rngFunc = [0,0,0,0,0,0]

    if rngFuncAllAbsolute:
      retIjkRange = ijkJson["ijkrange"]
      if PhactoriDbg(100):
        myDebugPrint3("rngFunc is [0,0,0,0,0,0], all absolute, "
          "so just return explicit ijkrange:\nretIjkRange: \n" + \
          str(retIjkRange) + "\n")
      return retIjkRange

    #return ijkJson["ijkrange"]
    if blockExtent == None:
      if rngFuncAllAbsolute == False:
        return [0,-1,0,-1,0,-1]

    if "ijkrange" in ijkJson:
      ijkRangeFromJson = ijkJson["ijkrange"]
    else:
      ijkRangeFromJson = [0,0,0,0,0,0]

    if PhactoriDbg(100):
      myDebugPrint3("ijkrange: " + str(ijkRangeFromJson) + "\n" + \
        "ijkrangefunction: " + str(rngFunc) + "\n" + \
        "blockExtent: " + str(blockExtent) + "\n")

    retIjkRange = [0,0,0,0,0,0]
    minIndex = [0,0,2,2,4,4]
    maxIndex = [1,1,3,3,5,5]
    for ii in range(0,6):
      if rngFunc[ii] == 0:
        retIjkRange[ii] = ijkRangeFromJson[ii]
      elif rngFunc[ii] < 0:
        retIjkRange[ii] = blockExtent[minIndex[ii]] + ijkRangeFromJson[ii]
      else:
        retIjkRange[ii] = blockExtent[maxIndex[ii]] - ijkRangeFromJson[ii]

    #validate/fix (with warning) extraction extents
    hadToCorrect = False
    uncorrectedIjkRange = list(retIjkRange)
    for ii in range(0,6):
      extentMin = blockExtent[minIndex[ii]]
      extentMax = blockExtent[maxIndex[ii]]
      if retIjkRange[ii] < extentMin:
        hadToCorrect = True
        retIjkRange[ii] = extentMin
      if retIjkRange[ii] > extentMax:
        hadToCorrect = True
        retIjkRange[ii] = extentMax
    for ii in range(0,3):
      ndx1 = ii * 2
      ndx2 = ndx1 + 1
      if retIjkRange[ndx1] > retIjkRange[ndx2]:
        hadToCorrect = True
        retIjkRange[ndx1] = retIjkRange[ndx2]
    if hadToCorrect:
      if PhactoriDbg(100):
        myDebugPrint3("return ijkrange had to be corrected because it was outside block extent"
          "\nblockExtent: " + str(blockExtent) +
          "\nijkRangeFromJson: " + str(retIjkRange) +
          "\nrngFunc (from json): " + str(rngFunc) +
          "\nuncorrectedIjkRange: " + str(uncorrectedIjkRange) +
          "\nretIjkRange: " + str(retIjkRange) + "\n")
    else:
      if PhactoriDbg(100):
        myDebugPrint3("return ijkrange did not need correction"
          "\nblockExtent: " + str(blockExtent) +
          "\nijkRangeFromJson: " + str(retIjkRange) +
          "\nrngFunc (from json): " + str(rngFunc) +
          "\nretIjkRange: " + str(retIjkRange) + "\n")

    return retIjkRange


  def MakeExtractSubsetFilter(self, extractBlockFilter, oneBlockExtractSubsetJson, oneBlockExtents):
    if PhactoriDbg(100):
      myDebugPrint3("MakeOneExtractBlockExtractSubsetFilter entered\n")

    #ijkrange = oneBlockExtractSubsetJson["ijkrange"]
    ijkrange = self.CreateIjkRangeFromSettingsForThisBlock(oneBlockExtractSubsetJson, oneBlockExtents)
    UpdatePipelineWithCurrentTimeArgument(extractBlockFilter)
    extractSubset1 = ExtractSubset(Input=extractBlockFilter)
    extractSubset1.VOI = ijkrange
    UpdatePipelineWithCurrentTimeArgument(extractSubset1)
    if PhactoriDbg(100):
      myDebugPrint3("MakeOneExtractBlockExtractSubsetFilter returning\n")
    return extractSubset1

  def MakeOneExtractBlockExtractSubsetFilterPair(self, inInputFilter, blockName, oneBlockExtractSubsetJson):
    if PhactoriDbg(100):
      myDebugPrint3("MakeOneExtractBlockExtractSubsetFilterPair entered: " + str(blockName) + "\n")
    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    newExtractBlockPvFilter = self.MakeExtractBlockFilter(inInputFilter, blockName)
    UpdatePipelineWithCurrentTimeArgument(newExtractBlockPvFilter)
    self.ExtractBlockFilters.append(newExtractBlockPvFilter)

    if PhactoriDbg(100):
      myDebugPrint3("pesmb one extract block structure begin: " + blockName + "\n")
      RecursivelyPrintMultiBlockStats(newExtractBlockPvFilter)
      myDebugPrint3("pesmb one extract block structure end: " + blockName + "\n")

    oneBlockExtent = self.mStructuredGridsIncludedExtent[-1]
    newExtractSubsetPvFilter = self.MakeExtractSubsetFilter(newExtractBlockPvFilter, oneBlockExtractSubsetJson, oneBlockExtent)
    UpdatePipelineWithCurrentTimeArgument(newExtractSubsetPvFilter)
    self.ExtractSubsetFilters.append(newExtractSubsetPvFilter)

    if PhactoriDbg(100):
      myDebugPrint3("pesmb one extract subset structure begin: " + blockName + "\n")
      RecursivelyPrintMultiBlockStats(newExtractSubsetPvFilter)
      myDebugPrint3("pesmb one extract subset structure end: " + blockName + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("MakeOneExtractBlockExtractSubsetFilterPair returning: " + str(blockName) + "\n")

  def MakeGroupAllExtractSubsetsFilter(self):
    if PhactoriDbg(100):
      myDebugPrint3("MakeGroupAllExtractSubsetsFilter entered\n")
    self.FinalGroupFilter = GroupDatasets(Input = self.ExtractSubsetFilters)

    UpdatePipelineWithCurrentTimeArgument(self.FinalGroupFilter)

    if PhactoriDbg(100):
      myDebugPrint3("MakeGroupAllExtractSubsetsFilter returning\n")

  def CreateParaViewFilter(self, inInputFilter):
    """create the PhactoriExtractStructuredMultiBlock filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractStructuredMultiBlock.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    if PhactoriDbg(100):
      myDebugPrint3("pesmb incoming multiblock structure begin\n")
      RecursivelyPrintMultiBlockStats(inInputFilter)
      myDebugPrint3("pesmb incoming multiblock structure end\n")

    for blockName, oneBlockExtractSubsetJson in self.ControlJson.items():
      self.MakeOneExtractBlockExtractSubsetFilterPair(inInputFilter, blockName, oneBlockExtractSubsetJson)

    self.MakeGroupAllExtractSubsetsFilter()
    #newParaViewFilter = PhactoriExtractStructuredMultiBlock(inInputFilter)

    if PhactoriDbg(100):
      myDebugPrint3("pesmb outgoing group multiblock structure begin\n")
      RecursivelyPrintMultiBlockStats(self.FinalGroupFilter)
      myDebugPrint3("pesmb outgoing group multiblock structure end\n")


    SetActiveSource(self.FinalGroupFilter)
    UpdatePipelineWithCurrentTimeArgument(self.FinalGroupFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractStructuredMultiBlock.CreateParaViewFilter "
          "returning\n", 100)

    return self.FinalGroupFilter

  def DoUpdateDueToChangeInData(self, inIncomingPvFilter, outOutgoingPvFilter):
    UpdatePipelineWithCurrentTimeArgument(self.FinalGroupFilter)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractStructuredMultiBlock printing out multiblock this callback:\n")
      myDebugPrint3("pesmb outgoing group multiblock structure begin (per frame callback)\n")
      RecursivelyPrintMultiBlockStats(self.FinalGroupFilter)
      myDebugPrint3("pesmb outgoing group multiblock structure end (per frame callback)\n")

#phactori_combine_to_single_python_file_subpiece_end_1
