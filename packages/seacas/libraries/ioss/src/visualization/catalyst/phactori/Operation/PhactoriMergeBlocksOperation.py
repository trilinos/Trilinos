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
class PhactoriMergeBlocksOperation(PhactoriOperationSpecifics):
  """merge blocks (MergeBlocks) operation, adapter to the catalyst filter

PhactoriMergeBlocksOperation is the phactori manager for working with the
ParaView/Catalyst MergeBlocks() filter and its parameters, providing
access and pipeline/input/output managment via the json, lexx/yacc, or soon
yaml interface. The user may specify a named input to the filter, with the
unnamed default being the incoming data mesh.

The MergeBlocks filter combines all the elements from all blocks of a
multiblock dataset into a single block with one set of elemeeents and one set
of points. The input is, e.g., a vtkMultiBlockDataset, and the output is
usually an vtkUnstructuredGrid or vtkPoly. Sometimes other filters which are
failing will succeed when placed after a MergeBlocks filter.

For moree information on the MergeBlocks() filter from ParaView, see the
ParaView Documentation.

To control the MergeBlocks you have a section in the
"operation blocks" section of the phactori control script which looks like:

::

    "operation blocks":{
      "mymergeblocks1":{
        "type":"mergeblocks"
      }
    }

As with other pipeline operations, you can add an
"input":"other_op_name" line to place this after another operation in the
pipeline such as a clip or a slice or a contour.
"""

  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mergePoints = True
    self.mergePartitionsOnly = False

  def ParseParametersFromJson(self, inJson):
    key1 = "merge points"
    if key1 in inJson:
      self.mergePoints = inJson[key1]
    key1 = "merge partitions only"
    if key1 in inJson:
      self.mergePartitionsOnly = inJson[key1]
    dummy = 0

  def CreateParaViewFilter(self, inInputFilter):
    """create the MergeBlocks filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMergeBlocksOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    newParaViewFilter = MergeBlocks(inInputFilter)
    if self.mergePoints:
      newParaViewFilter.MergePoints = 1
    else:
      newParaViewFilter.MergePoints = 0
    if self.mergePartitionsOnly:
      newParaViewFilter.MergePartitionsOnly = 1
    else:
      newParaViewFilter.MergePartitionsOnly = 0

    SetActiveSource(newParaViewFilter)
    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMergeBlocksOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
