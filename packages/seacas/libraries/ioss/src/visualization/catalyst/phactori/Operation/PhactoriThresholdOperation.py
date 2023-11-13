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
from Operation.PhactoriMpiUtilities import *
try:
  from .PhactoriParaviewMultiBlockRecursion import *
except:
  from Operation.PhactoriParaviewMultiBlockRecursion import *

#phactori_combine_to_single_python_file_subpiece_begin_1
class PhactoriThresholdOperation(PhactoriOperationSpecifics):
  """Threshold operation, phactori interface adapter to the catalyst filter

  PhactoriThresholdOperation is the phactori manager for working with the
  paraview/catalyst Threshold() operation and its paramters, providing
  access and pipeline/input/output managment via the json, lexx/yacc, or soon
  yaml interface. The user may specify a named input to the filter, with the
  unnamed default being the incoming data mesh.

  For detailed information on the Threshold filter from ParaView, see the
  ParaView Documentation. The basic operation of the Threshold filter is to
  extract all cells and points for which cell or point mesh variable (e.g.
  pressure, temperature, vonmises, displacement) is inside a specified range.

  To add a PhactoriThresholdOperation to the incoming script, you add
  a sub-block to the "operation blocks" section of the data with the "type"
  key given a value of "threshold". You also specify one of a range of values
  to keep, a value below which to keep values above, or a value above which to
  keep points/nodes. One complete but simple example
  script:

::

  {
    "camera blocks":{"mycam1":{"type":"camera", "look direction":[1.0, 2.0, 3.0]}},
    "representation blocks":{"rep_tmprtr":{"color by scalar":"temperature"}},
    "imageset blocks":{
      "temperature_on_slice_1":{
        "operation":"mytempthrsh1",
        "camera":"mycam1",
        "representation":"rep_tmprtr",
        "image basedirectory":"CatalystOutput",
        "image basename":"threshold1_temperature."
      }
    },
    "operation blocks":{
      "mytempthrsh1":{
        "type":"threshold",
        "variable scalar":"temperature",
        "keep between":[10.0, 100.0]
      }
    }
  }

Rather than "keep between" with a min/max range you can also have "keep above" with a single value or "keep below" with a single value, e.g.

::

  "keep above": 100.0
  or
  "keep below": 10.0

"""

  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    #self.mVariableName = None
    #self.mElementOrNode = ""
    self.mVariableInfo = PhactoriVariableInfo()
    self.mRange = [0.0, 1.0]

    #flag allows system to be set up to bypass this operation, especially
    #for testing camera settings against geometry without thresholds,
    #iso surfaces messing the tests up
    self.mBypassFlag = False

  def DoUpdateDueToChangeInData(self, inIncomingPvFilter,
      outOutgoingPvFilter):
    """the PhactoriThresholdOperation may need to update if the variable type
       was not detectable in a previous callback, but is now detectable"""
    if PhactoriDbg():
      myDebugPrint3("PhactoriThresholdOperation::"
          "DoUpdateDueToChangeInData override executing\n")

    #BarrierLock("PhactoriThresholdOperation sync 700")

    if self.mVariableInfo.mVariableTypeNeedsDetection == True:
      UpdatePipelineWithCurrentTimeArgument(inIncomingPvFilter)
      self.HandleVariableNeedingDetection(inIncomingPvFilter,
          outOutgoingPvFilter)
      if self.mVariableInfo.mVariableTypeNeedsDetection == False:
        UpdatePipelineWithCurrentTimeArgument(outOutgoingPvFilter)

    
  def ParseParametersFromJson(self, inJson):

    self.mBypassFlag = getParameterFromBlock(inJson, 'bypass flag',
        self.mBypassFlag)
    if self.mBypassFlag == True:
      #this operation is bypassed: don't do other parsing
      return

    foundVariableFlag = self.mVariableInfo.\
        ParseVariableNameAndVectorOrTensorComponent(inJson, 'variable ')

    if foundVariableFlag == False:
      errStr = """error!  inJson has no variable info in
                  PhactoriThresholdOperation.ParseParametersFromJson\n"""
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    #self.mVariableName = inJson['variable']

    #if 'variable type' in inJson:
    #  self.mElementOrNode = inJson['variable type']
    #  if self.mElementOrNode != 'element' and self.mElementOrNode != 'node':
    #    errStr = "error!  inJson 'variable type' is neither node nor\nelement PhactoriThresholdOperation.ParseParametersFromJson\n"
    #    myDebugPrint3(errStr)
    #    raise Exception(errStr)
    #else:
    #  errStr = "error!  inJson has no 'variable type' key in\nPhactoriThresholdOperation.ParseParametersFromJson\n"
    #  myDebugPrint3(errStr)
    #  raise Exception(errStr)

    if 'keep between' in inJson:
      self.mRange = inJson['keep between']
    elif 'keep below' in inJson:
      rangeMin = -sys.float_info.max
      rangeMax = inJson['keep below']
      self.mRange = [rangeMin, rangeMax]
    elif 'keep above' in inJson:
      rangeMin = inJson['keep above']
      rangeMax = sys.float_info.max
      self.mRange = [rangeMin, rangeMax]
    else:
      if PhactoriDbg():
        myDebugPrint3("  no keep between/above/below, using keep above 0.0\n")
      rangeMin = 0.0
      rangeMax = sys.float_info.max
      self.mRange = [rangeMin, rangeMax]

  def HandleVariableNeedingDetection(self, inIncomingPvFilter,
      outOutgoingPvFilter):

    detectResult = self.mVariableInfo.DetectVariableType(
        inIncomingPvFilter, True)

    if self.mVariableInfo.mVariableType == 'element':
      scalarType = 'CELLS'
    elif self.mVariableInfo.mVariableType == 'node':
      scalarType = 'POINTS'
    else:
      errStr = """error!  in PhactoriThresholdOperation.CreateParaViewFilter
                  variable type is not element or node or detect but:
                  """ + \
                  str(self.mVariableInfo.mVariableType) + \
                  """\n"""
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    #hack to deal with how paraview/catalyst expects naming of threshold
    #variables which are vector components or magnitudes

    localVariableName = GetThresholdContourHackVariableNameString(self.mVariableInfo)

    if PhactoriDbg():
      myDebugPrint3('  name for variable in threshold: ' + localVariableName + '\n')
    outOutgoingPvFilter.Scalars = [scalarType, localVariableName]

  def CreateParaViewFilter(self, inInputFilter):
    """create the threshold filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriThresholdOperation.CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    #if the operation bypass flag is set to true, we want this
    #operation to do nothing and just pass the input to the output
    if self.mBypassFlag == True:
      if PhactoriDbg():
        myDebugPrint3("PhactoriThresholdOperation::CreateParaViewFilter\n" +
            "mBypassFlag was true, so paraview output filter is input filter")
      return inInputFilter

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    #BarrierLock("PhactoriThresholdOperation sync 100")
    if PhactoriDbg():
      numCellArrays = inInputFilter.CellData.GetNumberOfArrays()
      myDebugPrint3("before threshold cell numCellArrays: " + str(numCellArrays) + "\n")
      for ii in range(0, numCellArrays):
        myDebugPrint3(inInputFilter.CellData.GetArray(ii).Name + "\n")

    newParaViewFilter = Threshold(inInputFilter)
    global gParaViewCatalystVersionFlag
    if gParaViewCatalystVersionFlag < 51000:
      newParaViewFilter.ThresholdRange = self.mRange
    else:
      newParaViewFilter.LowerThreshold = self.mRange[0]
      newParaViewFilter.UpperThreshold = self.mRange[1]
      newParaViewFilter.ThresholdMethod = vtk.vtkThreshold.THRESHOLD_BETWEEN

    self.HandleVariableNeedingDetection(inInputFilter, newParaViewFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)
    #AddFilterToFilterMap(ioOperationBlock.mName, newParaViewFilter)

    if PhactoriDbg(100):
      RecursivelyPrintPointAndCellArrayInformation(inInputFilter)

    #inputSource = GetActiveSource()
    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)
    if PhactoriDbg():
      numCellArrays = newParaViewFilter.CellData.GetNumberOfArrays()
      myDebugPrint3("after threshold cell numCellArrays: " + str(numCellArrays) + "\n")
      for ii in range(0, numCellArrays):
        myDebugPrint3(newParaViewFilter.CellData.GetArray(ii).Name + "\n")

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriThresholdOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1

