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
from .PhactoriOperationBlock import *
from .PhactoriVectorLibrary import *
from paraview.simple import *
from .PhactoriSegment import *
from .PhactoriMpiUtilities import *
from .PhactoriParaviewMultiBlockRecursion import *
from .PhactoriSampledCellInfo import *
import vtk
import json

#phactori_combine_to_single_python_file_subpiece_begin_1
class PhactoriStructuredGridSampler(PhactoriOperationSpecifics):
  """given a multiblock dataset which has leaves that ar structured grids,
     take specified i/j/k sample ranges from specified blocks.  Optionally
     make mask variable for sampled cells with programmable filter"""
  def __init__(self):
    self.myCopyOfInputFilter = None
    self.SampleControlSet = None

  def ValidateSampleList1(self, sampleList):
    if len(sampleList) < 1:
      myDebugPrint3AndException(
        "PhactoriStructuredGridSampler:ValidateSampleList1\n" \
        "list of samples must have at least one entry\n")

    for oneSample in sampleList:
      sampleBlockName = oneSample[0]
      sampleBlockType = oneSample[1]
      validType = False
      if sampleBlockType == "imax":
        validType = True
      if sampleBlockType == "imin":
        validType = True
      if sampleBlockType == "jmax":
        validType = True
      if sampleBlockType == "jmin":
        validType = True
      if sampleBlockType == "kmax":
        validType = True
      if sampleBlockType == "kmin":
        validType = True
      if validType != True:
        myDebugPrint3AndException(
          "PhactoriStructuredGridSampler:ValidateSampleList1\n" \
          "invalid sample type: " + str(sampleBlockType) + "\n" \
          "must be 'imin' 'imax' 'jmin' 'jmax' 'kmin' or 'kmax'\n")

  def ParseParametersFromJson(self, inJson):
    keyval1 = "samples to be taken from all blocks"
    if keyval1 in inJson:
      sampleControlList = inJson[keyval1]
      self.ValidateSampleList1(sampleControlList)
      self.SampleControlSet = {}
      for oneSampleControl in sampleControlList:
        self.SampleControlSet[oneSampleControl[0]] = oneSampleControl[1]

  @staticmethod
  def GatherStructuredSampledCellsInBlock(recursionObject, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("GatherStructuredSampledCellsInBlock entered\n")

    numCells = inInputCsData.GetNumberOfCells()
    numPoints = inInputCsData.GetNumberOfPoints()
    inParameters.leafVisitCount += 1
    if (numCells == 0) or (numPoints == 0):
      #no cells here
      return

    if PhactoriDbg(100):
      myDebugPrint3("GatherStructuredSampledCellsInBlock returning\n")

  def GatherStructuredSampledCellsOnThisProcess(self):
    if PhactoriDbg(100):
      myDebugPrint3("GatherStructuredCellsFromSeedCells entered\n")

    recursionParams = GatherStructuredSampledCellsOnThisProcessRecursionParams()
    recursionParams.SetUpForRecursion(self)
    recursionObj = PhactoriParaviewMultiBlockRecursionControl()
    recursionObj.mParameters = recursionParams
    recursionObj.mOperationToDoPerBlock = self.GatherStructuredSampledCellsInBlock
    PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, self.myCopyOfInputFilter)

    if PhactoriDbg(100):
      myDebugPrint3("GatherStructuredCellsFromSeedCells returning\n")
    return recursionParams.CellsPerSeedCell

    cellData = inInputCsData.GetCellData()
    outputCellArray = None
    #if cellData != None:
    #  outputCellArray = cellData.GetArray(inParameters.dataArrayName)

    if outputCellArray != None:
      dataArrayNumCmpnts = outputCellArray.GetNumberOfComponents()
      defaultTuple = []
      for ii in range(0, dataArrayNumCmpnts):
        defaultTuple.append(-1.0)
    else:
      dataArrayNumCmpnts = -1
      defaultTuple = []


  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriStructuredGridSampler.CreateParaViewFilter entered\n", 100)

    savedActiveSource = GetActiveSource()

    ##segmentListJson = ReadAndMpiBroadcastJsonFile(self.JsonListFileName)

    self.myCopyOfInputFilter = inInputFilter

    UpdatePipelineWithCurrentTimeArgument(self.myCopyOfInputFilter)

    #need to figure global extents
    self.GatherStructuredSampledCellsOnThisProcess()

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriStructuredGridSampler.CreateParaViewFilter returning\n", 100)

#phactori_combine_to_single_python_file_subpiece_end_1
