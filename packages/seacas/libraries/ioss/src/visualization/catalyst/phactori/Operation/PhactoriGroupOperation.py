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
class PhactoriGroupOperation(PhactoriOperationSpecifics):
  """manages group filter: group the results of two or more operations into
     a single multiblock source"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    #filled in from parsing
    self.mOperationNameList = None
    #filled in by finding operations id'd by self.mOperationNameList
    self.mOperationList = None

  def ParseParametersFromJson(self, inJson):

    if 'operation group list' in inJson:
      self.mOperationNameList = inJson['operation group list']
    else:
      myDebugPrint3AndException(
          "PhactoriGroupOperation::ParseParametersFromJson\n"
          "Error:  must have 'operation group list' token\n")

  def GetListOfInputOperationNamesForThisOperationType(self):
    """the group operation depends on all the operations in
       self.mOperationNameList"""
    retList = list(self.mOperationNameList)
    return retList

  def CreateParaViewFilter2(self, ioPipeAndViewsState):
    """create the group filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGroupOperation.CreateParaViewFilter "
          "entered\n", 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    #get direct list of operations (not through name list)
    #should be valid here
    self.mOperationList = []
    for ii in self.mOperationNameList:
      if ii not in ioPipeAndViewsState.mOperationBlocks:
        myDebugPrint3AndException(
          "PhactoriGroupOperation::CreateParaViewFilter\n"
          "Error:  operation '" + str(ii) + "' not available\n")
      inputOperationBlock = ioPipeAndViewsState.mOperationBlocks[ii]
      if inputOperationBlock.GetPvFilter() == None:
        myDebugPrint3AndException(
          "PhactoriGroupOperation::CreateParaViewFilter\n"
          "Error:  operation '" + str(ii) + "' paraview filter not "
          "constructed\n")
      self.mOperationList.append(inputOperationBlock.GetPvFilter())
      
    newParaViewFilter = GroupDatasets(Input = self.mOperationList)

    SetActiveSource(newParaViewFilter)
    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGroupOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter
#phactori_combine_to_single_python_file_subpiece_end_1

