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
import sys
from .PhactoriParaviewMultiBlockRecursion import *
from .PhactoriVectorLibrary import *

#phactori_combine_to_single_python_file_subpiece_begin_1

class CellEdgeLengthCalculationParameters():
  def __init__(self):
    self.numPointsPerBlock = None
    self.numCellsPerBlock = None
    self.OutputingFromShortest = True
    self.OutputShortestLongestIndex = 0

  def Initialize(self, inOutputingFromShortest, inOutputShortestLongestIndex):
    self.numPointsPerBlock = []
    self.numCellsPerBlock = []
    self.OutputingFromShortest = inOutputingFromShortest
    self.OutputShortestLongestIndex = inOutputShortestLongestIndex

class PhactoriFindCellEdgeLengths(PhactoriOperationSpecifics):
  """experimental filter to find the lengths of the edges of cells in the mesh"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.inPvFilter = None
    self.pfcelProgrammableFilter = None
    self.mProgFilterString = None
    self.ProgrammableFilterOutputCellVariableName = "mincelledgelength"
    self.OutputingFromShortest = True
    self.OutputShortestLongestIndex = 0

  def ParseParametersFromJson(self, inJson):

    keyval3 = "output nth shortest edge"
    keyval4 = "output nth longest edge"
    if keyval3 in inJson:
      self.OutputingFromShortest = True
      self.OutputShortestLongestIndex = inJson[keyval3]
    elif keyval4 in inJson:
      self.OutputingFromShortest = False
      self.OutputShortestLongestIndex = inJson[keyval3]
    else:
      self.OutputingFromShortest = True
      self.OutputShortestLongestIndex = 0

    keyval10 = "output cell variable name"
    if keyval10 in inJson:
      self.ProgrammableFilterOutputCellVariableName = inJson[keyval10]

  def CreateParaViewFilter(self, inInputFilter):
    """create the MergeBlocks filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellEdgeLengths.CreateParaViewFilter "
          "entered\n", 100)
    self.inPvFilter = inInputFilter

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    #SetActiveSource(newParaViewFilter)
    #UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)
    #SetActiveSource(savedActiveSource)

    self.pfcelProgrammableFilter = ProgrammableFilter(Input = self.inPvFilter)
    self.pfcelProgrammableFilter.CopyArrays = 1
    self.CreateProgrammableFilterString()
    self.pfcelProgrammableFilter.Script = self.mProgFilterString
    self.pfcelProgrammableFilter.UpdatePipeline()

    SetActiveSource(self.pfcelProgrammableFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellEdgeLengths.CreateParaViewFilter "
          "returning\n", 100)

    return self.pfcelProgrammableFilter

  @staticmethod
  def FindCellEdgeLengthsForOneCell(inInputCsData, inParameters, inCellIndex,
        paramOutputingFromShortest, paramOutputShortestLongestIndex):
    oneCell = inInputCsData.GetCell(inCellIndex)
    #if PhactoriDbg(100):
    #  myDebugPrint3("oneCell: " + str(oneCell) + "\n")
    numEdges = oneCell.GetNumberOfEdges()
    ptA = [0.0, 0.0, 0.0]
    ptB = [0.0, 0.0, 0.0]
    edgeLengthList = []
    for ii in range(0, numEdges):
      oneEdge = oneCell.GetEdge(ii)
      edgePoints = oneEdge.GetPointIds()
      ptAId = edgePoints.GetId(0)
      ptBId = edgePoints.GetId(1)
      inInputCsData.GetPoint(ptAId, ptA)
      inInputCsData.GetPoint(ptBId, ptB)
      edgeLength = vecDistance(ptA, ptB)
      edgeLengthList.append(edgeLength)
    sortedEdgeLengths = sorted(edgeLengthList)

    if paramOutputingFromShortest:
      outputIndex = min(paramOutputShortestLongestIndex, len(sortedEdgeLengths) - 1)
    else:
      outputIndex = max(0, len(sortedEdgeLengths) - 1 - paramOutputShortestLongestIndex)

    if PhactoriDbg(100):
      myDebugPrint3(
        "shortest edge length: " + str(sortedEdgeLengths[0]) + "\n"
        " longest edge length: " + str(sortedEdgeLengths[-1]) + "\n"
        "  output edge length: " + str(sortedEdgeLengths[outputIndex]) + "\n"
        " shortest/longest flag: " + str(paramOutputingFromShortest) + "\n"
        "shortest/longest index: " + str(paramOutputShortestLongestIndex) + "\n"
        "           outputIndex: " + str(outputIndex) + "\n")

  @staticmethod
  def FindCellEdgeLengthsForCellsInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("FindCellEdgeLengthsForCellsInBlock entered\n")

    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    if PhactoriDbg(100):
      myDebugPrint3("numCells: " + str(numCells) + \
        "  numPoints: " + str(numPoints) + "\n")
    inParameters.numCellsPerBlock.append(numCells)
    inParameters.numPointsPerBlock.append(numPoints)

    for ii in range(0, numCells):
      if PhactoriDbg(100):
        myDebugPrint3("cell index " + str(ii) + " of " + str(numCells) + "\n")
      PhactoriFindCellEdgeLengths.FindCellEdgeLengthsForOneCell(
        inInputCsData, inParameters, ii,
        inParameters.OutputingFromShortest,
        inParameters.OutputShortestLongestIndex)

    if PhactoriDbg(100):
      myDebugPrint3("FindCellEdgeLengthsForCellsInBlock returning\n")

  def ExportOperationData(self, datadescription):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellEdgeLengths.ExportOperationData entered\n",
        100)

    UpdatePipelineWithCurrentTimeArgument(self.inPvFilter)

    #recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    #recursionObj.mParameters = CellEdgeLengthCalculationParameters()
    #recursionObj.mParameters.Initialize(self.OutputingFromShortest,
    #  self.OutputShortestLongestIndex)
    #recursionObj.mOperationToDoPerBlock = self.FindCellEdgeLengthsForCellsInBlock
    #PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.inPvFilter)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriFindCellEdgeLengths.ExportOperationData returning\n",
        100)

  def CreateProgrammableFilterString(self):

    if self.mProgFilterString != None:
      return

    scriptLines = []

    scriptLines.append("#print('test2')\n")

    scriptLines.append("localLeafVisitCount = 0\n")

    scriptLines.append("import math\n")
    scriptLines.append("def FindEdgeLengthsForOneCell(inInputCsData, inCellIndex):\n")
    scriptLines.append("  oneCell = inInputCsData.GetCell(inCellIndex)\n")
    scriptLines.append("  numEdges = oneCell.GetNumberOfEdges()\n")
    scriptLines.append("  ptA = [0.0, 0.0, 0.0]\n")
    scriptLines.append("  ptB = [0.0, 0.0, 0.0]\n")
    scriptLines.append("  edgeLengthList = []\n")
    scriptLines.append("  for ii in range(0, numEdges):\n")
    scriptLines.append("    oneEdge = oneCell.GetEdge(ii)\n")
    scriptLines.append("    edgePoints = oneEdge.GetPointIds()\n")
    scriptLines.append("    ptAId = edgePoints.GetId(0)\n")
    scriptLines.append("    ptBId = edgePoints.GetId(1)\n")
    scriptLines.append("    inInputCsData.GetPoint(ptAId, ptA)\n")
    scriptLines.append("    inInputCsData.GetPoint(ptBId, ptB)\n")
    scriptLines.append("    dx = ptB[0] - ptA[0]\n")
    scriptLines.append("    dy = ptB[1] - ptA[1]\n")
    scriptLines.append("    dz = ptB[2] - ptA[2]\n")
    scriptLines.append("    edgeLength = dx*dx + dy*dy + dz*dz\n")
    scriptLines.append("    edgeLengthList.append(edgeLength)\n")
    scriptLines.append("  sortedEdgeLengths = sorted(edgeLengthList)\n")
    pickIndexStr = str(self.OutputShortestLongestIndex)
    if self.OutputingFromShortest:
      scriptLines.append("  outputIndex = min("+pickIndexStr+", len(sortedEdgeLengths) - 1)\n")
    else:
      scriptLines.append("  outputIndex = max(0, len(sortedEdgeLengths) - 1 - "+pickIndexStr+")\n")

    scriptLines.append("  returnEdgeLength = math.sqrt(sortedEdgeLengths[outputIndex])\n")
    scriptLines.append("  return returnEdgeLength\n")

    scriptLines.append("def flatten(input, output):\n")
    scriptLines.append("    # Copy the cells etc.\n")
    scriptLines.append("    output.ShallowCopy(input)\n")
    scriptLines.append("    numPoints = input.GetNumberOfPoints()\n")
    scriptLines.append("    celldata = output.GetCellData()\n")
    scriptLines.append("    numCells = input.GetNumberOfCells()\n")
    scriptLines.append("    #print(str(numCells))\n")
    scriptLines.append("    ncda = vtk.vtkDoubleArray()\n")
    scriptLines.append("    ncda.SetNumberOfTuples(numCells)\n")
    scriptLines.append("    numTuples = ncda.GetNumberOfTuples()\n")
    scriptLines.append("    for ii in range(0, numTuples):\n")
    scriptLines.append("        ncda.SetValue(ii, -1.0)\n")
    scriptLines.append("    for ii in range(0, numCells):\n")
    scriptLines.append("      minEdgeLengthForCell = FindEdgeLengthsForOneCell(input, ii)\n")
    scriptLines.append("      ncda.SetValue(ii, minEdgeLengthForCell)\n")
    scriptLines.append("    ncda.SetName('" + self.ProgrammableFilterOutputCellVariableName + "')\n")
    scriptLines.append("    #print(ncda.GetName())\n")
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
      myDebugPrint3("PhactoriFindCellEdgeLengths constructed script:\n" + \
        self.mProgFilterString)

#phactori_combine_to_single_python_file_subpiece_end_1
