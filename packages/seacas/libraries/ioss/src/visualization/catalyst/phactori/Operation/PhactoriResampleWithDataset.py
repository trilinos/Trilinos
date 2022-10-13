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
class PhactoriResampleWithDatasetOperation(PhactoriOperationSpecifics):
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.ResampleDatasetName = None
    self.InternalParaViewFilterPtr = None
    return

  def ParseParametersFromJson(self, inJson):
    key1 = "sample dataset"
    if key1 not in inJson:
      myDebugPrint3AndException("PhactoriResampleWithDatasetOperation:ParseParametersFromJson\n"
        "'" + key1 + "' must be in block to specify dataset as source for the resampling\n")
    self.ResampleDatasetName = inJson[key1]

  def GetListOfInputOperationNamesForThisOperationType(self):
    """PhactoriResampleWithDatasetOperation depends on the sampling source as
       well as the input source"""
    retList = []
    retList.append(self.ResampleDatasetName)
    return retList

  def DebugPrintInputPortAndOutputPortInfo(self, infoTag):
    if PhactoriDbg(100):
      myDebugPrint3(infoTag)
      stcso = self.InternalParaViewFilterPtr.GetClientSideObject()
      myDebugPrint3(str(stcso) + "\n")
      myDebugPrint3(str(stcso.GetClassName()) + "\n")
      stodo = stcso.GetOutputDataObject(0)
      myDebugPrint3(str(stodo) + "\n")
      myDebugPrint3(str(stodo.GetClassName()) + "\n")
      stido0 = stcso.GetInputDataObject(0,0)
      myDebugPrint3(str(stido0) + "\n")
      myDebugPrint3(str(stido0.GetClassName()) + "\n")
      stido1 = stcso.GetInputDataObject(1,0)
      myDebugPrint3(str(stido1) + "\n")
      myDebugPrint3(str(stido1.GetClassName()) + "\n")
    #  myDebugPrint3("stodo: " + str(stodo.GetClassName()) + " numcells " + str(stodo.GetNumberOfCells()) + " numpoints " + str(stodo.GetNumberOfPoints()) + "\n")
    #  myDebugPrint3("stido0: " + str(stido0.GetClassName()) + "\n")
    #  myDebugPrint3("stido1: " + str(stido1.GetClassName()) + " numcells " + str(stido1.GetNumberOfCells()) + " numpoints " + str(stido1.GetNumberOfPoints()) + "\n")

  def CreateParaViewFilter2(self, inInputFilter, inPipeAndViewsState):
    """create the clip plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriResampleWithDatasetOperation:CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    #get seed source operation
    paraviewResampleDataset = inPipeAndViewsState.GetParaViewSourceByOperationName(self.ResampleDatasetName)
    if PhactoriDbg(100):
      myDebugPrint3("self.ResampleDatasetName: " + str(self.ResampleDatasetName) + \
        "\nparaviewResampleDataset: " + str(paraviewResampleDataset) + "\n", 100)

    UpdatePipelineWithCurrentTimeArgument(paraviewResampleDataset)

    if PhactoriDbg(100):
      psscso = paraviewResampleDataset.GetClientSideObject()
      pssodo = psscso.GetOutputDataObject(0)
      myDebugPrint3("pssodo: " + str(pssodo.GetClassName()) + " numcells " + str(pssodo.GetNumberOfCells()) + " numpoints " + str(pssodo.GetNumberOfPoints()) + "\n")

    savedActiveSource = GetActiveSource()

    #self.InternalParaViewFilterPtr = ResampleWithDataset(Input = inInputFilter, Source=paraviewResampleDataset)
    self.InternalParaViewFilterPtr = ResampleWithDataset(SourceDataArrays = inInputFilter, DestinationMesh=paraviewResampleDataset)
    self.InternalParaViewFilterPtr.PassCellArrays = 1
    self.InternalParaViewFilterPtr.PassPointArrays = 1

    self.InternalParaViewFilterPtr.CellLocator = 'Static Cell Locator'

    self.DebugPrintInputPortAndOutputPortInfo("pre filter update 1\n")

    SetActiveSource(self.InternalParaViewFilterPtr)
    SetActiveSource(savedActiveSource)

    self.DebugPrintInputPortAndOutputPortInfo("pre filter update 2\n")

    #UpdatePipelineWithCurrentTimeArgument(self.InternalParaViewFilterPtr)

    if PhactoriDbg(100):
      myDebugPrint3("about to call self.InternalParaViewFilterPtr.UpdatePipeline()\n")

    self.InternalParaViewFilterPtr.UpdatePipeline()

    self.DebugPrintInputPortAndOutputPortInfo("post filter update 1\n")

    if PhactoriDbg(100):
      myDebugPrint3("self.InternalParaViewFilterPtr: " + str(self.InternalParaViewFilterPtr) + "\n", 100)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriResampleWithDatasetOperation.CreateParaViewFilter returning\n', 100)

    return self.InternalParaViewFilterPtr

#phactori_combine_to_single_python_file_subpiece_end_1
