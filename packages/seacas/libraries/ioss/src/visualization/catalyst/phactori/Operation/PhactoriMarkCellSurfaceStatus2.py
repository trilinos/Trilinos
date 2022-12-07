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
import sys
from .PhactoriParaviewMultiBlockRecursion import *
from .PhactoriVectorLibrary import *
from .PhactoriVtkCellOperations import *

#phactori_combine_to_single_python_file_subpiece_begin_1

class CellEdgeAngleMetricParameters():
  def __init__(self):
    self.numPointsPerBlock = None
    self.numCellsPerBlock = None
    self.UseSmallestAngle = True
    self.OffsetIndex = 0
    self.OutputCellVariableName = None

  def Initialize(self, inUseSmallestAngle, inOffsetIndex, inCellVariableName):
    self.numPointsPerBlock = []
    self.numCellsPerBlock = []
    self.UseSmallestAngle = inUseSmallestAngle
    self.OffsetIndex = inOffsetIndex
    self.OutputCellVariableName = inCellVariableName

class PhactoriMarkCellSurfaceStatus2(PhactoriOperationSpecifics):
  """experimental filter to find the lengths of the edges of cells in the mesh"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.inPvFilter = None
    self.pfcelProgrammableFilter = None
    self.mProgFilterString = None
    self.OutputCellVariableName = "surfacestatus"
    self.UseSmallestAngle = True
    self.OffsetIndex = 0
    self.DoOperationInCreateParaViewMethod = True
    self.CreateProgrammableFilter = False

  def ParseParametersFromJson(self, inJson):

    keyval3 = "output nth smallest angle"
    keyval4 = "output nth largest angle"
    if keyval3 in inJson:
      self.UseSmallestAngle = True
      self.OffsetIndex = inJson[keyval3]
    elif keyval4 in inJson:
      self.UseSmallestAngle = False
      self.OffsetIndex = inJson[keyval3]
    else:
      self.UseSmallestAngle = True
      self.OffsetIndex = 0

    keyval10 = "output cell variable name"
    if keyval10 in inJson:
      self.OutputCellVariableName = inJson[keyval10]

    keyval11 = "do operation in createparaview method"
    if keyval11 in inJson:
      self.DoOperationInCreateParaViewMethod = inJson[keyval11]

    keyval12 = "create programmable filter"
    if keyval12 in inJson:
      self.ProgrammableFilterOutputCellVariableName = inJson[keyval10]

  def CreateParaViewFilter(self, inInputFilter):
    """create the MergeBlocks filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMarkCellSurfaceStatus2.CreateParaViewFilter "
          "entered\n", 100)
    self.inPvFilter = inInputFilter

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    if self.DoOperationInCreateParaViewMethod:
      if PhactoriDbg(100):
        myDebugPrint3("make variable during CreateParaViewFilter\n", 100)
      recursionObj = PhactoriParaviewMultiBlockRecursionControl()
      recursionObj.mParameters = CellEdgeAngleMetricParameters()
      recursionObj.mParameters.Initialize(self.UseSmallestAngle,
        self.OffsetIndex, self.OutputCellVariableName)
      recursionObj.mOperationToDoPerBlock = self.MarkCellStatusInBlock
      PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.inPvFilter)

      UpdatePipelineWithCurrentTimeArgument(inInputFilter)
      if PhactoriDbg(100):
        myDebugPrint3("done make variable during CreateParaViewFilter\n", 100)

    if self.CreateProgrammableFilter:
      self.pfcelProgrammableFilter = ProgrammableFilter(Input = self.inPvFilter)
      self.pfcelProgrammableFilter.CopyArrays = 1
      self.CreateProgrammableFilterString()
      self.pfcelProgrammableFilter.Script = self.mProgFilterString
      self.pfcelProgrammableFilter.UpdatePipeline()
      UpdatePipelineWithCurrentTimeArgument(self.pfcelProgrammableFilter)
      if PhactoriDbg(100):
        myDebugPrint3("self.CreateProgrammableFilter true, created programmable filter\n", 100)

    SetActiveSource(savedActiveSource)

    if self.CreateProgrammableFilter:
      if PhactoriDbg(100):
        myDebugPrint3("PhactoriMarkCellSurfaceStatus2.CreateParaViewFilter "
          "returning programmable filter\n", 100)
      return self.pfcelProgrammableFilter
    else:
      if PhactoriDbg(100):
        myDebugPrint3("PhactoriMarkCellSurfaceStatus2.CreateParaViewFilter "
          "returning passthrough filter\n", 100)
      return self.inPvFilter

  @staticmethod
  def MarkCellStatusInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("MarkCellStatusInBlock entered\n")

    numCells = inInputCsData.GetNumberOfCells()

    celldata = inInputCsData.GetCellData()

    newCellArrayName = inParameters.OutputCellVariableName

    inParameters.numCellsPerBlock.append(numCells)
    ncda = vtk.vtkIntArray()
    ncda.SetNumberOfTuples(numCells)

    if numCells <= 0:
      ncda.SetName(newCellArrayName)
      celldata.AddArray(ncda)

      if PhactoriDbg(100):
        myDebugPrint3("numCells 0, MarkCellStatusInBlock returning 2\n")

    #see how many cells are touching each point
    cellsTouchingPointCount = PhactoriCountCellTouchingEachPoint(inInputCsData)

    #see which cells have how many points with 8 or other cells
    for ii in range(0, numCells):
      #if PhactoriDbg(100):
      #  myDebugPrint3("(b) cell index " + str(ii) + " of " + str(numCells) + "\n")
      oneCellSurfaceStatus = PhactoriFindSurfaceStatusForOneCell(inInputCsData, ii,
        cellsTouchingPointCount)
      ncda.SetValue(ii, oneCellSurfaceStatus)

    ncda.SetName(newCellArrayName)
    celldata.AddArray(ncda)

    if PhactoriDbg(100):
      myDebugPrint3("MarkCellStatusInBlock returning\n")

  def CreateProgrammableFilterString(self):

    if self.mProgFilterString != None:
      return

    scriptLines = []

    scriptLines.append("import math\n")
    GetPhactoriVectorLibraryProgrammableFilterLines(scriptLines)
    GetPhactoriVtkCellOperationsProgrammableFilterLines(scriptLines)

    scriptLines.append("localLeafVisitCount = 0\n")
    scriptLines.append("def flatten(input, output):\n")
    scriptLines.append("    # Copy the cells etc.\n")
    scriptLines.append("    output.ShallowCopy(input)\n")
    scriptLines.append("    celldata = output.GetCellData()\n")
    scriptLines.append("    numCells = input.GetNumberOfCells()\n")
    scriptLines.append("    ncda = vtk.vtkIntArray()\n")
    scriptLines.append("    ncda.SetNumberOfTuples(numCells)\n")

    #the first three of the next lines use face connections to find
    #interior/surface/edge/corner cells, and the next three use point
    #connections. Point connections create 1/3 the size intermediate data as
    #face connectinos, but face connections should work with wedge cells too
    #scriptLines.append("    cellsTouchingPointCount = PhactoriCountCellTouchingEachPoint(input)\n")
    #scriptLines.append("    for ii in range(0, numCells):\n")
    #scriptLines.append("      oneCellSurfaceStatus = PhactoriFindSurfaceStatusForOneCell(input, ii, cellsTouchingPointCount)\n")
    scriptLines.append("    cellsTouchingEachFace = PhactoriCountCellTouchingEachFace(input)\n")
    scriptLines.append("    for ii in range(0, numCells):\n")
    scriptLines.append("      oneCellSurfaceStatus = PhactoriFindSurfaceStatusForOneCellUsingFaceInfo(input, ii, cellsTouchingEachFace)\n")

    scriptLines.append("      ncda.SetValue(ii, oneCellSurfaceStatus)\n")
    scriptLines.append("    ncda.SetName('" + self.OutputCellVariableName + "')\n")
    scriptLines.append("    celldata.AddArray(ncda)\n")

    scriptLines.append("input = self.GetInputDataObject(0, 0)\n")
    scriptLines.append("output = self.GetOutputDataObject(0)\n")

    scriptLines.append("if input.IsA('vtkMultiBlockDataSet'):\n")
    scriptLines.append("    output.CopyStructure(input)\n")
    scriptLines.append("    iter = input.NewIterator()\n")
    scriptLines.append("    iter.UnRegister(None)\n")
    scriptLines.append("    iter.InitTraversal()\n")
    scriptLines.append("    while not iter.IsDoneWithTraversal():\n")
    scriptLines.append("        localLeafVisitCount += 1\n")
    scriptLines.append("        curInput = iter.GetCurrentDataObject()\n")
    scriptLines.append("        curOutput = curInput.NewInstance()\n")
    scriptLines.append("        curOutput.UnRegister(None)\n")
    scriptLines.append("        output.SetDataSet(iter, curOutput)\n")
    scriptLines.append("        flatten(curInput, curOutput)\n")
    scriptLines.append("        iter.GoToNextItem();\n")
    scriptLines.append("else:\n")
    scriptLines.append("  flatten(input, output)\n")

    self.mProgFilterString = "".join(scriptLines)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMarkCellSurfaceStatus2 constructed script:\n" + \
        self.mProgFilterString)

#phactori_combine_to_single_python_file_subpiece_end_1
