# Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
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
    self.mExpandedIncludeBlockList = None
    self.mExcludeBlockList = None
    self.mRemoveBlanksAndDashes = True

  def ParseParametersFromJson(self, inJson):

    if 'include blocks' in inJson:
      self.mIncludeBlockList = inJson['include blocks']

    if 'exclude blocks' in inJson:
      self.mExcludeBlockList = inJson['exclude blocks']

    if self.mIncludeBlockList == None and self.mExcludeBlockList == None:
      myDebugPrint3AndException(
          "PhactoriExtractBlockOperation::ParseParametersFromJson\n"
          "Error:  must have 'include block list' or 'exclude block list'\n")

    #2022Aug09 only include blocks is working
    if self.mIncludeBlockList == None:
      myDebugPrint3AndException(
          "PhactoriExtractBlockOperation::ParseParametersFromJson\n"
          "Error:  must have 'include block list'\n"
          "currently (2022Aug09) exclude blocks is not working\n")

  def MakeExpandedIncludeBlockList(self):
    if self.mExpandedIncludeBlockList == None:
      self.mExpandedIncludeBlockList = list(self.mIncludeBlockList)
      for oneItem in self.mIncludeBlockList:
        #add leader, remove spaces and dashes, this might change
        if self.mRemoveBlanksAndDashes:
          baseItem = oneItem.replace(" ","")
          baseItem = baseItem.replace("-","")
        else:
          baseItem = oneItem
        extraItem = "/Root/ElementBlocks/" + baseItem
        self.mExpandedIncludeBlockList.append(extraItem)
        extraItem = "/Root/Bases/Base/Zones/" + baseItem
        self.mExpandedIncludeBlockList.append(extraItem)
    if PhactoriDbg(100):
      myDebugPrint3("self.mExpandedIncludeBlockList:\n" + \
        str(self.mExpandedIncludeBlockList) + "\n")

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
    self.MakeExpandedIncludeBlockList()
    newParaViewFilter.Selectors = self.mExpandedIncludeBlockList

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
