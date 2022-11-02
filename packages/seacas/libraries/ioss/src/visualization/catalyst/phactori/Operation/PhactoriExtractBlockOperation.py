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
global gDuplicateNameCounter
gDuplicateNameCounter = {}

class PhactoriExtractBlockOperation(PhactoriOperationSpecifics):
  """Extract Block operation, phactori interface adapter to the catalyst filter
For detailed information on the ExtractBlock() filter from ParaView, see the
ParaView Documentation. The basic operation of this filter is to take a list
of block names and create a new mesh which includes only the blocks in the
which have the names in the given list (and their children). Alternatively
the list can contain blocks to be excluded instead of included. The current
ParaView (5.8) implementation of this filter takes a list of block indices
which can become inaccurate if the mesh structure is unexpected, so this
Phactori adapter uses the block names to create the block index list at
run time rather than at script creation time.

To add a PhactoriExtractBlock to the incoming script, you add
a sub-block to the "operation blocks" section of the data with the "type"
key given a value of "extractblock". You also supply a list of block names
with a key of either "include blocks" or "exclude blocks". One complete but
simple example script:

::

  {
    "camera blocks":{"mycam1":{"type":"camera", "look direction":[1.0, 2.0, 3.0]}},
    "representation blocks":{"rep_tmprtr":{"color by scalar":"temperature"}},
    "imageset blocks":{
      "temperature_on_slice_1":{
        "operation":"myincludeblocks1",
        "camera":"mycam1",
        "representation":"rep_tmprtr",
        "image basedirectory":"CatalystOutput",
        "image basename":"myincludeblocks1_temperature."
      }
    },
    "operation blocks":{
      "myincludeblocks1":{
        "type":"extractblock",
        "include blocks":["blk-2", "blk-3", "blk-7"]
      }
    }
  }

As mentioned previously, you can specify "exclude blocks" instead of
"include blocks", but you may not specify both.

"""

  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    #mIncludeblockList OR mExcludeBlockList will be filled in from parsing
    self.mIncludeBlockList = None
    self.mExcludeBlockList = None
    self.mFlatBlockIndicesSpecifiedDirectly = False

    #this will be list of included block indicies, and is calculated from
    #mIncludeBlockList / mExcludeBlockList and passed directly to
    #ExcludeBlockFilter.BlockIndices
    self.mBlockIndices = []

  def ParseParametersFromJson(self, inJson):

    if 'include blocks' in inJson:
      self.mIncludeBlockList = inJson['include blocks']

    if 'exclude blocks' in inJson:
      self.mExcludeBlockList = inJson['exclude blocks']

    if 'include blocks by flat index list' in inJson:
      self.mFlatBlockIndicesSpecifiedDirectly = True
      self.mBlockIndices = inJson['include blocks by flat index list']

    if self.mIncludeBlockList == None and self.mExcludeBlockList == None and self.mFlatBlockIndicesSpecifiedDirectly == False:
      myDebugPrint3AndException(
          "PhactoriExtractBlockOperation::ParseParametersFromJson\n"
          "Error:  must have 'include block list' or 'exclude block list' or 'include blocks by flat index list'\n")

  def FigureBlockIndicesFromBlockListOneBlock(self, inMetaData,
          ioFlatIndexCounter):
    """determine if this one block should have it's flat index tripped on for
       the extract block filter (leaf item of recursion)"""

    if inMetaData == None:
      thisBlockName = None
    else:
      thisBlockName = inMetaData.Get(vtk.vtkCompositeDataSet.NAME())

    if self.mIncludeBlockList != None:
      if thisBlockName == None:
        if PhactoriDbg(100):
          myDebugPrint3("block with no name " + \
              " not in include list, not + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
      elif thisBlockName in self.mIncludeBlockList:
        self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
        global gDuplicateNameCounter
        if thisBlockName in gDuplicateNameCounter:
          oldCount = gDuplicateNameCounter[thisBlockName]
          gDuplicateNameCounter[thisBlockName] = oldCount+1
        else:
          gDuplicateNameCounter[thisBlockName] = 1
        if PhactoriDbg(100):
          myDebugPrint3("block " + str(thisBlockName) + \
              " in include list, + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
      else:
        if PhactoriDbg(100):
          myDebugPrint3("block " + str(thisBlockName) + \
              " not in include list, not + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
    elif self.mExcludeBlockList != None:
      if thisBlockName == None:
        self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
        if PhactoriDbg(100):
          myDebugPrint3("block with no name " + \
              " not in exclude list, + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
      elif thisBlockName not in self.mExcludeBlockList:
        self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
        if PhactoriDbg(100):
          myDebugPrint3("block " + str(thisBlockName) + \
              " not in exclude list, + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
      else:
        if PhactoriDbg(100):
          myDebugPrint3("block " + str(thisBlockName) + \
              " in exclude list, not + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
    else:
      myDebugPrint3AndException(
          "PhactoriExtractBlockOperation::"
          "FigureBlockIndicesFromBlockListOneBlock\n"
          "Error:  must have include block list or exclude block list\n")

  def FigureBlockIndicesFromBlockListRecurse1(self, inCsdata, inMetaData,
          ioFlatIndexCounter, inForceSetting):
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
        self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
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
        if self.mIncludeBlockList != None:
          if (thisBlockName != None) and (thisBlockName in self.mIncludeBlockList):
            if PhactoriDbg(100):
              myDebugPrint3("include list, append nonleaf " + str(thisBlockName) + \
                " ndx " + str(int(ioFlatIndexCounter[0]) - 1) + "\n")
            self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
            includeOrExcludeAllSubBlocks = 1
            if PhactoriDbg(100):
              myDebugPrint3("includeOrExcludeAllSubBlocks set to 1\n")
          else:
            if PhactoriDbg(100):
              myDebugPrint3("include list, exclude nonleaf " + str(thisBlockName) + \
                " ndx " + str(int(ioFlatIndexCounter[0]) - 1) + "\n")
        elif self.mExcludeBlockList != None:
          if (thisBlockName == None) or (thisBlockName not in self.mExcludeBlockList):
            if PhactoriDbg(100):
              myDebugPrint3("exclude list, append nonleaf " + str(thisBlockName) + \
                " ndx " +str(int(ioFlatIndexCounter[0]) - 1) + "\n")
            self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
          else:
            if PhactoriDbg(100):
              myDebugPrint3("exclude list, exclude nonleaf " + str(thisBlockName) + \
                " ndx " + str(int(ioFlatIndexCounter[0]) - 1) + "\n")
          if (thisBlockName != None) and (thisBlockName in self.mExcludeBlockList):
            includeOrExcludeAllSubBlocks = -1
            if PhactoriDbg(100):
              myDebugPrint3("includeOrExcludeAllSubBlocks set to -1\n")
        else:
          myDebugPrint3AndException(
            "PhactoriExtractBlockOperation::"
            "FigureBlockIndicesFromBlockListRecurse1\n"
            "Error:  must have include block list or exclude block list\n")

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
          self.FigureBlockIndicesFromBlockListRecurse1(oneBlock,
              oneBlockMetaData, ioFlatIndexCounter, includeOrExcludeAllSubBlocks)
        else:
          if PhactoriDbg(100):
            myDebugPrint3("this sub block is None, now handle it\n")
          #I think we need to count here to be okay with pruned stuff; maybe
          #we need to set extract block to no pruning (?)
          ioFlatIndexCounter[0] += 1
          if includeOrExcludeAllSubBlocks == 1:
            if PhactoriDbg(100):
              myDebugPrint3("leaf block None index " + \
              str(int(ioFlatIndexCounter[0]) - 1) + \
              " appended due to includeOrExcludeAllSubBlocks == 1\n")
            self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
          else:
            if PhactoriDbg(100):
              myDebugPrint3("force include is not on, so check include list\n")
            includeThisBlock = False
            if (self.mIncludeBlockList != None) and (oneBlockMetaData != None):
              thisBlockName = oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME())
              if PhactoriDbg(100):
                myDebugPrint3("include list not None, thisBlockName " + str(thisBlockName) + "\n"
                  "self.mIncludeBlockList " + str(self.mIncludeBlockList) + "\n")
                myDebugPrint3("thisBlockName != None: " + str(thisBlockName != None) + "\n")
                myDebugPrint3("in list: " + str(thisBlockName in self.mIncludeBlockList) + "\n")
              if (thisBlockName != None) and (thisBlockName in self.mIncludeBlockList):
                includeThisBlock = True
                if PhactoriDbg(100):
                  myDebugPrint3("leaf block None index " + \
                  str(int(ioFlatIndexCounter[0]) - 1) + \
                  " included due to name being in self.mIncludeBlockList\n")
                self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
            if includeThisBlock == False:
              if PhactoriDbg(100):
                myDebugPrint3("leaf block None index " + \
                str(int(ioFlatIndexCounter[0]) - 1) + \
                " excluded due to includeOrExcludeAllSubBlocks != 1\n")
    else:
      if inForceSetting == 1:
        self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
        if PhactoriDbg(100):
          myDebugPrint3("this leaf block added inForceSetting == 1\n")
      elif inForceSetting != -1:
        self.FigureBlockIndicesFromBlockListOneBlock(inMetaData,
          ioFlatIndexCounter)
      else:
        if PhactoriDbg(100):
          myDebugPrint3("this leaf block not added inForceSetting == -1\n")
    #if icsdClassname == "vtkMultiPieceDataSet":
    #  numpieces = inCsdata.GetNumberOfPieces()
    #  # ioFlatIndexCounter[0] += numpieces - 1
    #  ioFlatIndexCounter[0] += numpieces


  def FigureBlockIndicesFromBlockList(self, inInputFilter):
    """from the list of include/exclude blocks create a indices list
       of blocks to put in filter"""
    if PhactoriDbg(100):
      myDebugPrint3("FigureBlockIndicesFromBlockList entered"
          "\nself.mIncludeBlockList: " + str(self.mIncludeBlockList) + \
          "\nself.mExcludeBlockList: " + str(self.mExcludeBlockList) + "\n")

    global gDuplicateNameCounter
    gDuplicateNameCounter = {}
    csdata = inInputFilter.GetClientSideObject().GetOutputDataObject(0)
    flatIndexCounter = [0]
    self.FigureBlockIndicesFromBlockListRecurse1(csdata, None,
        flatIndexCounter, 0)

    if PhactoriDbg(100):
      myDebugPrint3("number of times existing block names found:\n")
      for blknm, count in gDuplicateNameCounter.items():
        myDebugPrint3(blknm + ": " + str(count) + "\n");
    if PhactoriDbg(100):
      myDebugPrint3("FigureBlockIndicesFromBlockList final indices list:\n" \
          + str(self.mBlockIndices) + \
          "\nFigureBlockIndicesFromBlockList returning\n")

  def CreateParaViewFilter(self, inInputFilter):
    """create the extract block filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractBlockOperation.CreateParaViewFilter "
          "entered\n", 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    if PhactoriDbg():
      myDebugPrint3("  extractblock inInputFilter point data arrays:\n")
      numArrays = inInputFilter.PointData.GetNumberOfArrays()
      for ii in range (0, numArrays):
        myDebugPrint3("  " + str(ii) + ":  " + inInputFilter.PointData.GetArray(ii).GetName() + "\n")

    if PhactoriDbg():
      myDebugPrint3("  extractblock inInputFilter cell data arrays:\n")
      numArrays = inInputFilter.CellData.GetNumberOfArrays()
      for ii in range (0, numArrays):
        myDebugPrint3("  " + str(ii) + ":  " + inInputFilter.CellData.GetArray(ii).GetName() + "\n")

    newParaViewFilter = ExtractBlock(inInputFilter)

    if self.mFlatBlockIndicesSpecifiedDirectly == False:
      self.FigureBlockIndicesFromBlockList(inInputFilter)

    #newParaViewFilter.PruneOutput = 1
    #newParaViewFilter.MaintainStructure = 0
    newParaViewFilter.MaintainStructure = 1

    newParaViewFilter.BlockIndices = self.mBlockIndices

    SetActiveSource(newParaViewFilter)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)
    if PhactoriDbg():
      myDebugPrint3("  extractblock newParaViewFilter point data arrays:\n")
      numArrays = newParaViewFilter.PointData.GetNumberOfArrays()
      for ii in range (0, numArrays):
        myDebugPrint3("  " + str(ii) + ":  " + newParaViewFilter.PointData.GetArray(ii).GetName() + "\n")

    if PhactoriDbg():
      myDebugPrint3("  extractblock newParaViewFilter cell data arrays:\n")
      numArrays = newParaViewFilter.CellData.GetNumberOfArrays()
      for ii in range (0, numArrays):
        myDebugPrint3("  " + str(ii) + ":  " + newParaViewFilter.CellData.GetArray(ii).GetName() + "\n")

    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractBlockOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter
#phactori_combine_to_single_python_file_subpiece_end_1
