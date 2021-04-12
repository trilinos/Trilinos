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

from paraview.simple import *
from phactori import *

#phactori_combine_to_single_python_file_subpiece_begin_1
class PhactoriParaviewMultiBlockRecursionControl:
  """see PhactoriDoParaviewMultiblockRecursion();
     mOperationPerBlock should be set to a method
     which takes 1 parameter, mParameters should be set to the parameter
     instance which will be passed to the mOperationPerBlock call"""
  def __init__(self):
    self.mOperationToDoPerBlock = None
    self.mParameters = None

def PhactroiParaviewDoMethodPerBlockRecurse1(inInputCsData, inRecursionControlItem):
  """PhactroiParaviewDoMethodPerBlockRecurse1 is a generic method for doing recursion through
     the multiblock dataset and doing something (a callback) on a per leaf block
     basis.  Called by PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter which got clientside data and calls
     itself on internal nodes and calls
     inRecursionControlItem.mOperationToDoPerBlock on leaf block nodes"""
  if PhactoriDbg(100):
    myDebugPrint3('PhactroiParaviewDoMethodPerBlockRecurse1 entered\n', 100)

  icsdClassname = inInputCsData.GetClassName()
  if icsdClassname == "vtkMultiBlockDataSet" or \
     icsdClassname == "vtkExodusIIMultiBlockDataSet":
    myDebugPrint3('recursing: ' + icsdClassname + '\n')
    numBlocks = inInputCsData.GetNumberOfBlocks()
    for ii in range(0, numBlocks):
      oneBlock = inInputCsData.GetBlock(ii)
      if PhactoriDbg(100):
        DebugPrintBlockName(inInputCsData, ii)
      if(oneBlock != None):
        PhactroiParaviewDoMethodPerBlockRecurse1(oneBlock, inRecursionControlItem)
  else:
    inRecursionControlItem.mOperationToDoPerBlock(inRecursionControlItem, inInputCsData,
            inRecursionControlItem.mParameters)

  if PhactoriDbg(100):
    myDebugPrint3('PhactroiParaviewDoMethodPerBlockRecurse1 returning\n', 100)

def DebugPrintBlockName(csData, blockIndex):
  if PhactoriDbg(100):
    oneBlock = csData.GetBlock(blockIndex)
    if oneBlock != None:
      oneBlockMetaData = csData.GetMetaData(blockIndex)
      if oneBlockMetaData != None:
        #myDebugPrint3("oneBlockMetaData: " + str(oneBlockMetaData) + "\n")
        theBlockName = oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME())
        myDebugPrint3("block index, name: " + str(blockIndex) + ", " + str(theBlockName) + "\n")
      else:
        myDebugPrint3("oneBlockMetaData is None (1)\n")
    else:
      myDebugPrint3("this block is None, now check meta data\n")
      oneBlockMetaData = csData.GetMetaData(blockIndex)
      if oneBlockMetaData != None:
        #myDebugPrint3("oneBlockMetaData: " + str(oneBlockMetaData) + "\n")
        theBlockName = oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME())
        myDebugPrint3("block index, name (2): " + str(blockIndex) + ", " + str(theBlockName) + "\n")
      else:
        myDebugPrint3("oneBlockMetaData is None (2)\n")

def PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(inRecursionControlItem,  inPvFilter):
  """grab client side data object, and use that to do recursion"""  
  pvClientSideData = inPvFilter.GetClientSideObject().GetOutputDataObject(0)
  if pvClientSideData == None:
    if PhactoriDbg(100):
      myDebugPrint3(
        'DoMethodPerBlock: pvClientSideData is None, returning',100)
    return

  PhactroiParaviewDoMethodPerBlockRecurse1(pvClientSideData, inRecursionControlItem)

class PrintPointAndCellArrayInformationRecursionParams:
  def __init__(self):
    self.leafVisitCount = 0

def PrintPointAndCellArrayInformationInBlock(recursionObject, inInputCsData, inParameters):
  inParameters.leafVisitCount += 1
  myDebugPrint3("inParameters.leafVisitCount " + str(inParameters.leafVisitCount) + "\n")

  celldata = inInputCsData.GetCellData()
  numcellarrays = celldata.GetNumberOfArrays()
  myDebugPrint3("number of cell data arrays: " + str(numcellarrays) + "\n")
  numcelltuples = celldata.GetNumberOfTuples()
  myDebugPrint3("number of cell data tuples: " + str(numcelltuples) + "\n")
  for ii in range(0, numcellarrays):
    myDebugPrint3(str(ii) + ": " + str(celldata.GetArray(ii).GetNumberOfComponents()) + ": " + str(celldata.GetArrayName(ii)) + "\n")

  pointdata = inInputCsData.GetPointData()
  numpointarrays = pointdata.GetNumberOfArrays()
  myDebugPrint3("number of point data arrays: " + str(numpointarrays) + "\n")
  numpointtuples = pointdata.GetNumberOfTuples()
  myDebugPrint3("number of point data tuples: " + str(numpointtuples) + "\n")
  for ii in range(0, numpointarrays):
    myDebugPrint3(str(ii) + ": " + str(pointdata.GetArray(ii).GetNumberOfComponents()) + ": " + str(pointdata.GetArrayName(ii)) + "\n")

def RecursivelyPrintPointAndCellArrayInformation(pvFilter):
  if PhactoriDbg(100):
    myDebugPrint3("RecursivelyPrintPointAndCellArrayInformation entered\n")

  recursionParams = PrintPointAndCellArrayInformationRecursionParams()
  recursionObj = PhactoriParaviewMultiBlockRecursionControl()
  recursionObj.mParameters = recursionParams
  recursionObj.mOperationToDoPerBlock = PrintPointAndCellArrayInformationInBlock

  PhactoriRecusivelyDoMethodPerBlockFromParaViewFilter(recursionObj, pvFilter)
  if PhactoriDbg(100):
    myDebugPrint3("RecursivelyPrintPointAndCellArrayInformation returning\n")

#phactori_combine_to_single_python_file_subpiece_end_1
