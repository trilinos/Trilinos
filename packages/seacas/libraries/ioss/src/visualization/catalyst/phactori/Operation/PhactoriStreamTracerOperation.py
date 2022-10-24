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
class PhactoriStreamTracerOperation(PhactoriOperationSpecifics):
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.LineSourcePoint1 = [-0.5, -0.5, -0.5]
    self.LineSourcePoint2 = [0.5, 0.5, 0.5]
    self.Resolution = 500
    return

  def ParseParametersFromJson(self, inJson):
    key1 = "line source point 1"
    self.LineSourcePoint1 = getParameterFromBlock(inJson, key1, self.LineSourcePoint1)
    if len(self.LineSourcePoint1) != 3:
      myDebugPrint3AndException(
        "'" + key1 + "' must be list with three floats for x, y, z\n")

    key2 = "line source point 2"
    self.LineSourcePoint2 = getParameterFromBlock(inJson, key2, self.LineSourcePoint2)
    if len(self.LineSourcePoint2) != 3:
      myDebugPrint3AndException(
        "'" + key2 + "' must be list with three floats for x, y, z\n")

    key3 = "resolution"
    self.Resolution = getParameterFromBlock(inJson, key3, self.Resolution)
    if isinstance(self.Resolution, int) == False:
      myDebugPrint3AndException(
        "'" + key3 + "' must be an integer\n")


  def CreateParaViewFilter(self, inInputFilter):
    """create the clip plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriStreamTracerOperation:CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    savedActiveSource = GetActiveSource()

    newParaViewFilter = StreamTracer(Input = inInputFilter, SeedType="High Resolution Line Source")
    newParaViewFilter.Vectors = ["CELLS", "velocity_xyz"]
    newParaViewFilter.MaximumStreamlineLength = 100.0
    newParaViewFilter.SeedType.Point1 = self.LineSourcePoint1
    newParaViewFilter.SeedType.Point2 = self.LineSourcePoint2
    newParaViewFilter.SeedType.Resolution = self.Resolution
    newParaViewFilter.IntegrationDirection = 'FORWARD'
    newParaViewFilter.MaximumSteps = 5000
    #newParaViewFilter.InitialStepLength = 0.01

    if PhactoriDbg(100):
      stcso = newParaViewFilter.GetClientSideObject()
      stodo = stcso.GetOutputDataObject(0)
      myDebugPrint3("stodo st: " + str(stodo.GetClassName()) + " numcells " + str(stodo.GetNumberOfCells()) + " numpoints " + str(stodo.GetNumberOfPoints()) + "\n")
      stido0 = stcso.GetInputDataObject(0,0)
      myDebugPrint3("stido0 st: " + str(stido0.GetClassName()) + "\n")
      stido1 = stcso.GetInputDataObject(1,0)
      myDebugPrint3("stido1 st: " + str(stido1.GetClassName()) + " numcells " + str(stido1.GetNumberOfCells()) + " numpoints " + str(stido1.GetNumberOfPoints()) + "\n")

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      stcso = newParaViewFilter.GetClientSideObject()
      stodo = stcso.GetOutputDataObject(0)
      myDebugPrint3("B stodo st: " + str(stodo.GetClassName()) + " numcells " + str(stodo.GetNumberOfCells()) + " numpoints " + str(stodo.GetNumberOfPoints()) + "\n")
      stido0 = stcso.GetInputDataObject(0,0)
      myDebugPrint3("B stido0 st: " + str(stido0.GetClassName()) + "\n")
      stido1 = stcso.GetInputDataObject(1,0)
      myDebugPrint3("B stido1 st: " + str(stido1.GetClassName()) + " numcells " + str(stido1.GetNumberOfCells()) + " numpoints " + str(stido1.GetNumberOfPoints()) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriStreamTracerOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
