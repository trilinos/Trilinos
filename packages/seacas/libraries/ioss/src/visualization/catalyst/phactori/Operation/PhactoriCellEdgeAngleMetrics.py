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
    self.AngleCellVariableName = "minedgenormalangle"
    self.HeightCellVariableName = "cellheight"

  def Initialize(self, inUseSmallestAngle, inOffsetIndex,
      inAngleCellVariableName, inHeightCellVariableName):
    self.numPointsPerBlock = []
    self.numCellsPerBlock = []
    self.UseSmallestAngle = inUseSmallestAngle
    self.OffsetIndex = inOffsetIndex
    self.AngleCellVariableName = inAngleCellVariableName
    self.HeightCellVariableName = inHeightCellVariableName

class PhactoriCellEdgeAngleMetrics(PhactoriOperationSpecifics):
  """experimental filter to find the lengths of the edges of cells in the mesh"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.inPvFilter = None
    self.pfcelProgrammableFilter = None
    self.mProgFilterString = None
    self.OutputAngleCellVariableName = "minedgenormalangle"
    self.OutputHeightCellVariableName = "cellheight"
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

    keyval10 = "output angle cell variable name"
    if keyval10 in inJson: 
      self.OutputAngleCellVariableName = inJson[keyval10]
    keyval10 = "output height cell variable name"
    if keyval10 in inJson: 
      self.OutputHeightCellVariableName = inJson[keyval10]

    keyval11 = "do operation in createparaview method"
    if keyval11 in inJson:
      self.DoOperationInCreateParaViewMethod = inJson[keyval11]

    keyval12 = "create programmable filter"
    if keyval12 in inJson:
      self.CreateProgrammableFilter = inJson[keyval12]

  def CreateParaViewFilter(self, inInputFilter):
    """create the MergeBlocks filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriCellEdgeAngleMetrics.CreateParaViewFilter "
          "entered\n", 100)
    self.inPvFilter = inInputFilter

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(self.inPvFilter)

    if self.DoOperationInCreateParaViewMethod:
      if PhactoriDbg(100):
        myDebugPrint3("point and cell array info before calculating metrics:")
        RecursivelyPrintPointAndCellArrayInformation(self.inPvFilter)

      recursionObj = PhactoriParaviewMultiBlockRecursionControl()
      recursionObj.mParameters = CellEdgeAngleMetricParameters()
      recursionObj.mParameters.Initialize(self.UseSmallestAngle,
        self.OffsetIndex,
        self.OutputAngleCellVariableName,
        self.OutputHeightCellVariableName)
      recursionObj.mOperationToDoPerBlock = self.FindCellEdgeAngleMetricsInBlock
      PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.inPvFilter)

      UpdatePipelineWithCurrentTimeArgument(self.inPvFilter)

      if PhactoriDbg(100):
        myDebugPrint3("point and cell array info after calculating metrics:")
        RecursivelyPrintPointAndCellArrayInformation(self.inPvFilter)

    #SetActiveSource(newParaViewFilter)
    #UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if self.CreateProgrammableFilter:
      self.pfcelProgrammableFilter = ProgrammableFilter(Input = self.inPvFilter)
      self.pfcelProgrammableFilter.CopyArrays = 1
      self.CreateProgrammableFilterString()
      self.pfcelProgrammableFilter.Script = self.mProgFilterString
      self.pfcelProgrammableFilter.UpdatePipeline()

    SetActiveSource(savedActiveSource)

    if self.CreateProgrammableFilter:
      if PhactoriDbg(100):
        myDebugPrint3("PhactoriCellEdgeAngleMetrics.CreateParaViewFilter "
          "returning programmable filter\n", 100)
      return self.pfcelProgrammableFilter
    else:
      if PhactoriDbg(100):
        myDebugPrint3("PhactoriCellEdgeAngleMetrics.CreateParaViewFilter "
          "returning passthrough filter\n", 100)
      return self.inPvFilter

  @staticmethod
  def FindCellEdgeAngleMetricsInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("FindCellEdgeAngleMetricsInBlock entered\n")

    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    if PhactoriDbg(100):
      myDebugPrint3("numCells: " + str(numCells) + \
        "  numPoints: " + str(numPoints) + "\n")
    inParameters.numCellsPerBlock.append(numCells)
    inParameters.numPointsPerBlock.append(numPoints)

    ncda = vtk.vtkDoubleArray()
    ncda.SetNumberOfTuples(numCells)
    ncda2 = vtk.vtkDoubleArray()
    ncda2.SetNumberOfTuples(numCells)

    outputAngles = []
    outputHeights = []
    for ii in range(0, numCells):
      #if PhactoriDbg(100):
      #  myDebugPrint3("cell index " + str(ii) + " of " + str(numCells) + "\n")
      edgeAngle, cellHeight = PhactoriFindCellEdgeAngleMetricsForOneCell(inInputCsData, ii,
        inParameters.UseSmallestAngle, inParameters.OffsetIndex)
      outputAngles.append(edgeAngle)
      outputHeights.append(cellHeight)
      ncda.SetValue(ii, edgeAngle)
      ncda2.SetValue(ii, cellHeight)

    celldata = inInputCsData.GetCellData()
    ncda.SetName(inParameters.AngleCellVariableName)
    celldata.AddArray(ncda)
    ncda2.SetName(inParameters.HeightCellVariableName)
    celldata.AddArray(ncda2)
    if PhactoriDbg(100):
      myDebugPrint3("numCells: " + str(numCells) + "\n"
        "inParameters.AngleCellVariableName: " + inParameters.AngleCellVariableName + "\n"
        "inParameters.HeightCellVariableName: " + inParameters.HeightCellVariableName + "\n")
      if numCells > 0:
        myDebugPrint3("edgeAngle 0: " + str(ncda.GetValue(0)) + "\n")

    if PhactoriDbg(100):
      minOutputAngle = 91.0
      minOutputAngleIndex = -1
      maxOutputAngle = -1.0
      maxOutputAngleIndex = -1
      minHeight = 1e15
      minHeightIndex = -1
      maxHeight = -1.0
      maxHeightIndex = -1
      for ii in range(0, numCells):
        edgeAngle = outputAngles[ii]
        cellHeight = outputAngles[ii]
        if edgeAngle < minOutputAngle:
          minOutputAngle = edgeAngle
          minOutputAngleIndex = ii
        if edgeAngle > maxOutputAngle:
          maxOutputAngle = edgeAngle
          maxOutputAngleIndex = ii
        if cellHeight < minHeight:
          minHeight = cellHeight
          minHeightIndex = ii
        if cellHeight > maxHeight:
          maxHeight = cellHeight
          maxHeightIndex = ii
      myDebugPrint3(
        "minOutputAngle, index: " + str(minOutputAngle) + ", " + str(minOutputAngleIndex) + "\n" + \
        "maxOutputAngle, index: " + str(maxOutputAngle) + ", " + str(maxOutputAngleIndex) + "\n" + \
        "minHeight, index: " + str(minHeight) + ", " + str(minHeightIndex) + "\n" + \
        "maxHeight, index: " + str(maxHeight) + ", " + str(maxHeightIndex) + "\n" + \
        "FindCellEdgeAngleMetricsInBlock returning\n")
    return outputAngles, outputHeights

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
    scriptLines.append("    numPoints = input.GetNumberOfPoints()\n")
    scriptLines.append("    celldata = output.GetCellData()\n")
    scriptLines.append("    numCells = input.GetNumberOfCells()\n")
    scriptLines.append("    #print(str(numCells))\n")
    scriptLines.append("    ncda = vtk.vtkDoubleArray()\n")
    scriptLines.append("    ncda.SetNumberOfTuples(numCells)\n")
    scriptLines.append("    ncda2 = vtk.vtkDoubleArray()\n")
    scriptLines.append("    ncda2.SetNumberOfTuples(numCells)\n")
    scriptLines.append("    numTuples = ncda.GetNumberOfTuples()\n")
    scriptLines.append("    for ii in range(0, numTuples):\n")
    scriptLines.append("        ncda.SetValue(ii, -1.0)\n")
    scriptLines.append("        ncda2.SetValue(ii, -1.0)\n")
    scriptLines.append("    for ii in range(0, numCells):\n")
    scriptLines.append("      cellAngle, cellHeight = PhactoriFindCellEdgeAngleMetricsForOneCell(input, ii, " + \
      str(self.UseSmallestAngle) + ", " + str(self.OffsetIndex) + ")\n")
    scriptLines.append("      ncda.SetValue(ii, cellAngle)\n")
    scriptLines.append("      ncda2.SetValue(ii, cellHeight)\n")
    scriptLines.append("    ncda.SetName('" + self.OutputAngleCellVariableName + "')\n")
    scriptLines.append("    celldata.AddArray(ncda)\n")
    scriptLines.append("    ncda2.SetName('" + self.OutputHeightCellVariableName + "')\n")
    scriptLines.append("    celldata.AddArray(ncda2)\n")

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
      myDebugPrint3("PhactoriCellEdgeAngleMetrics constructed script:\n" + \
        self.mProgFilterString)

#phactori_combine_to_single_python_file_subpiece_end_1
