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
class PhactoriStreamTracerSeedSourceOperation(PhactoriOperationSpecifics):
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.SeedSourceName = None
    self.InternalParaViewFilterPtr = None
    return

  def ParseParametersFromJson(self, inJson):
    key1 = "seed source"
    if key1 not in inJson:
      myDebugPrint3AndException("PhactoriStreamTracerSeedSourceOperation:ParseParametersFromJson\n"
        "'" + key1 + "' must be in block to specify source for seed points\n")
    self.SeedSourceName = inJson[key1]

  def GetListOfInputOperationNamesForThisOperationType(self):
    """PhactoriStreamTracerSeedSourceOperation depends on the seed sources as
       well as the input source"""
    retList = []
    retList.append(self.SeedSourceName)
    return retList

  def DebugPrintInputPortAndOutputPortInfo(self, infoTag):
    if PhactoriDbg(100):
      myDebugPrint3(infoTag)
      stcso = self.InternalParaViewFilterPtr.GetClientSideObject()
      stodo = stcso.GetOutputDataObject(0)
      myDebugPrint3("stodo: " + str(stodo.GetClassName()) + " numcells " + str(stodo.GetNumberOfCells()) + " numpoints " + str(stodo.GetNumberOfPoints()) + "\n")
      stido0 = stcso.GetInputDataObject(0,0)
      myDebugPrint3("stido0: " + str(stido0.GetClassName()) + "\n")
      stido1 = stcso.GetInputDataObject(1,0)
      myDebugPrint3("stido1: " + str(stido1.GetClassName()) + " numcells " + str(stido1.GetNumberOfCells()) + " numpoints " + str(stido1.GetNumberOfPoints()) + "\n")

  def CreateParaViewFilter2(self, inInputFilter, inPipeAndViewsState):
    """create the clip plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriStreamTracerSeedSourceOperation:CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    #get seed source operation
    paraviewSeedSource = inPipeAndViewsState.GetParaViewSourceByOperationName(self.SeedSourceName)
    if PhactoriDbg(100):
      myDebugPrint3("self.SeedSourceName: " + str(self.SeedSourceName) + \
        "\nparaviewSeedSource: " + str(paraviewSeedSource) + "\n", 100)

    UpdatePipelineWithCurrentTimeArgument(paraviewSeedSource)

    if PhactoriDbg(100):
      psscso = paraviewSeedSource.GetClientSideObject()
      pssodo = psscso.GetOutputDataObject(0)
      myDebugPrint3("pssodo: " + str(pssodo.GetClassName()) + " numcells " + str(pssodo.GetNumberOfCells()) + " numpoints " + str(pssodo.GetNumberOfPoints()) + "\n")

    savedActiveSource = GetActiveSource()

    self.InternalParaViewFilterPtr = StreamTracerWithCustomSource(Input = inInputFilter, SeedSource=paraviewSeedSource)

    #psscso = paraviewSeedSource.GetClientSideObject()
    #pssodo = psscso.GetOutputDataObject(0)
    #stcso = self.InternalParaViewFilterPtr.GetClientSideObject()
    #stcso.SetInputDataObject(1,pssodo)

    self.InternalParaViewFilterPtr.Vectors = ["CELLS", "velocity_xyz"]
    #self.InternalParaViewFilterPtr.Vectors = ["POINTS", "velocity_xyz"]
    self.InternalParaViewFilterPtr.MaximumStreamlineLength = 100.0
    self.InternalParaViewFilterPtr.IntegrationDirection = 'FORWARD'
    self.InternalParaViewFilterPtr.MaximumSteps = 5000
    #self.InternalParaViewFilterPtr.InitialStepLength = 0.01

    self.DebugPrintInputPortAndOutputPortInfo("pre filter update 1\n")

    SetActiveSource(self.InternalParaViewFilterPtr)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(self.InternalParaViewFilterPtr)

    self.DebugPrintInputPortAndOutputPortInfo("post filter update 1\n")

    if PhactoriDbg(100):
      myDebugPrint3("self.InternalParaViewFilterPtr: " + str(self.InternalParaViewFilterPtr) + "\n", 100)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriStreamTracerSeedSourceOperation.CreateParaViewFilter returning\n', 100)

    return self.InternalParaViewFilterPtr

#phactori_combine_to_single_python_file_subpiece_end_1
