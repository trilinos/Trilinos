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
class PhactoriCalculatorOperation(PhactoriOperationSpecifics):
  """manages calculator filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mFunction = ""
    self.mResultArrayName = 'Result'
    self.mPointOrCell = 'Point Data'

  def ParseParametersFromJson(self, inJson):
    if 'function' in inJson:
      #function is supposed to be array of terms as strings; concatenate them
      #to get actual function string; this is this way due to sierra parsing
      #self.mFunction = "".join(inJson['function'])
      #
      #since numbers are coming across as floats, we need to convert them
      #to strings for the paraview calculator function as well
      functionStringItems = []
      for oneItem in inJson['function']:
          functionStringItems.append(str(oneItem))
      self.mFunction = "".join(functionStringItems)
    else:
      myDebugPrint3AndException(
          "PhactoriCalculatorOperation::ParseParametersFromJson\n"
          "Error:  must have 'function' key\n")
    if 'resultname' in inJson:
      self.mResultArrayName = inJson['resultname']
    if 'output variable name' in inJson:
      self.mResultArrayName = inJson['output variable name']
    if 'vector ' in inJson:
      self.mResultArrayName = inJson['resultname']
    if 'element or node data' in inJson:
      if inJson['element or node data'] == 'node':
        self.mPointOrCell = 'Point Data'
      elif inJson['element or node data'] == 'element':
        self.mPointOrCell = 'Cell Data'
      else:
        myDebugPrint3AndException(
          "PhactoriCalculatorOperation::ParseParametersFromJson\n"
          "element or node data must be 'node' or 'element'\n")
    elif 'nodeorelement' in inJson:
      if inJson['nodeorelement'] == 'node':
        self.mPointOrCell = 'Point Data'
      elif inJson['nodeorelement'] == 'element':
        self.mPointOrCell = 'Cell Data'
      else:
        myDebugPrint3AndException(
          "PhactoriCalculatorOperation::ParseParametersFromJson\n"
          "nodeorelement data must be 'node' or 'element'\n")
    elif 'pointorcell' in inJson:
      if inJson['pointorcell'] == 'point':
        self.mPointOrCell = 'Point Data'
      elif inJson['pointorcell'] == 'cell':
        self.mPointOrCell = 'Cell Data'
      else:
        myDebugPrint3AndException(
          "PhactoriCalculatorOperation::ParseParametersFromJson\n"
          "node or element data must be 'node' or 'point' or 'element' or 'cell'\n")

  def CreateParaViewFilter(self, inInputFilter):
    """create the calculator filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriCalculatorOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    if PhactoriDbg():
      myDebugPrint3("  calculator inInputFilter point data arrays:\n")
      numArrays = inInputFilter.PointData.GetNumberOfArrays()
      for ii in range (0, numArrays):
        myDebugPrint3("  " + str(ii) + ":  " + inInputFilter.PointData.GetArray(ii).GetName() + "\n")
  
    if PhactoriDbg():
      myDebugPrint3("  calculator inInputFilter cell data arrays:\n")
      numArrays = inInputFilter.CellData.GetNumberOfArrays()
      for ii in range (0, numArrays):
        myDebugPrint3("  " + str(ii) + ":  " + inInputFilter.CellData.GetArray(ii).GetName() + "\n")


    newParaViewFilter = Calculator(inInputFilter)
    mypid = SmartGetLocalProcessId()
    if self.mFunction.find("Hat") < 0:
      constfunc = "11.5"
    else:
      constfunc = "1.0*iHat+2.0*jHat+3.0*kHat"

    #self.mPointOrCell = 'Cell Data'

    if self.mPointOrCell == "Cell Data":
      if inInputFilter.CellData.GetNumberOfArrays() == 0:
        myDebugPrint3("process " + str(mypid) + " has no cell arrays, using function " + constfunc + "\n")
        newParaViewFilter.Function = constfunc
      #elif str(inInputFilter.CellData.GetArray(0).GetName()) != "V":
      #  myDebugPrint3("process " + str(mypid) + " 0th cell array not V, using function " + constfunc + "\n")
      #  newParaViewFilter.Function = constfunc
      else:
        myDebugPrint3("process " + str(mypid) + " has cell arrays, using function " + self.mFunction + "\n")
        newParaViewFilter.Function = self.mFunction
    else:
      if inInputFilter.PointData.GetNumberOfArrays() == 0:
        myDebugPrint3("process " + str(mypid) + " has no point arrays, using function " + constfunc + "\n")
        newParaViewFilter.Function = constfunc
      else:
        myDebugPrint3("process " + str(mypid) + " has cell arrays, using function " + self.mFunction + "\n")
        newParaViewFilter.Function = self.mFunction
      #newParaViewFilter.AttributeMode = 'Cell Data'
      #newParaViewFilter.AttributeMode = 'Point Data'

    if gParaViewCatalystVersionFlag < 50502:
      newParaViewFilter.AttributeMode = self.mPointOrCell
    else:
      newParaViewFilter.AttributeType = self.mPointOrCell

    newParaViewFilter.ResultArrayName = self.mResultArrayName

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3("function: " + str(newParaViewFilter.Function) + "\n"
        "result array name: " + str(newParaViewFilter.ResultArrayName) + "\n"
        "PhactoriTransformOperation.CreateParaViewFilter returning\n", 100)

    return newParaViewFilter
#phactori_combine_to_single_python_file_subpiece_end_1
