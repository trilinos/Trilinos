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
class PhactoriContourOperation(PhactoriOperationSpecifics):
  """clip plane operation, adapter to the catalyst Contour filter

PhactoriContourOperation is the phactori manager for working with the
ParaView/Catalyst Contour() filter and its parameters, providing
access and pipeline/input/output managment via the json, lexx/yacc, or soon
yaml interface. The user may specify a named input to the filter, with the
unnamed default being the incoming data mesh. i

For information on the Contour() filter from ParaView, see the ParaView
Documentation. The basic operation is to take a value and a variable and
calculate an iso surface in the mesh corresponding to the value and create
a 3D surface at that iso value. Multiple values can be specified in the same
operation, resulting in multiple surface.

The user must specify the variable to be used and one or more isosurface
values.

To add a PhactoriContourOperation to the incoming script, you add
a sub-block to the "operation blocks" section of the data with the "type"
key given a value of "contour". One complete but simple example script:

::

  {
    "camera blocks":{"mycam1":{"type":"camera", "look direction":[1.0, 2.0, 3.0]}},
    "representation blocks":{"rep_tmprtr":{"color by scalar":"temperature"}},
    "imageset blocks":{
      "temperature_contour_1":{
        "operation":"myclipwithplane1",
        "camera":"mycam1",
        "representation":"rep_tmprtr",
        "image basedirectory":"CatalystOutput",
        "image basename":"contour1_temperature."
      }
    },
    "operation blocks":{
      "mycontour1":{
        "type":"contour",
        "variable scalar":"temperature",
        "contour value":[100.0, 150.0, 200.0, 1000.0]
      }
    }


As an alternative to "contour value" you are allowed to specify a key
"contour value sequence" with a start value, a step value, and an end value.
Contour values will automatically be created from start, spacing by step,
until a value greater than or equal to end is reached.  For example:


::

  "contour value sequence": [100.0, 200.0, 1100.0]


will get 100, 300, 500, 700, and 900. If it is 1101.0 for end, 1100.0 will
also be used.

"""

  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mVariableInfo = PhactoriVariableInfo()
    self.mContourValue = [0.0]

    #flag allows system to be set up to bypass this operation, especially
    #for testing camera settings against geometry without thresholds,
    #iso surfaces messing the tests up
    self.mBypassFlag = False

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
                  PhactoriContourOperation.ParseParametersFromJson\n"""
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    if 'contour value' in inJson:
      self.mContourValue = inJson['contour value']
    elif 'contour value sequence' in inJson:
      seqCntrl = inJson['contour value sequence']
      rStart = seqCntrl[0]
      rStep = seqCntrl[1]
      rEnd = seqCntrl[2]
      self.mContourValue = []
      rVal = rStart
      if rEnd > rStart:
        if rStep <= 0.0:
          errStr = """error!  bad contour value sequence in
                      PhactoriContourOperation.ParseParametersFromJson\n"""
          if PhactoriDbg():
            myDebugPrint3(errStr)
          raise Exception(errStr)
        while rVal < rEnd:
          self.mContourValue.append(rVal)
          rVal += rStep
      else:
        if rStep >= 0.0:
          errStr = """error!  bad contour value sequence in
                      PhactoriContourOperation.ParseParametersFromJson\n"""
          if PhactoriDbg():
            myDebugPrint3(errStr)
          raise Exception(errStr)
        while rVal > rEnd:
          self.mContourValue.append(rVal)
          rVal += rStep
    else:
      errStr = """error!  inJson has no 'contour value' or 'contour value sequence' key in
                  PhactoriContourOperation.ParseParametersFromJson\n"""
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)


  def CreateParaViewFilter(self, inInputFilter):
    """create the clip plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriContourOperation.CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    #if the operation bypass flag is set to true, we want this
    #operation to do nothing and just pass the input to the output
    if self.mBypassFlag == True:
      if PhactoriDbg():
        myDebugPrint3("PhactoriContourOperation::CreateParaViewFilter\n" +
            "mBypassFlag was true, so paraview output filter is input filter")
      return inInputFilter

    savedActiveSource = GetActiveSource()
    SetActiveSource(inInputFilter)
    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    newParaViewFilter = Contour(Input = inInputFilter,
        PointMergeMethod="Uniform Binning")
    newParaViewFilter.ComputeScalars = 1

    if PhactoriDbg():
      myDebugPrint3("  contour filter is: " + str(newParaViewFilter) + "\n")
    newParaViewFilter.PointMergeMethod = "Uniform Binning"

    detectResult = self.mVariableInfo.DetectVariableType(inInputFilter, True)
    #if detectResult == False:
      #PutOperationOnListToDetectVariableTypeInFutureSteps()

    #hack to deal with how paraview/catalyst expects naming of threshold
    #variables which are vector components or magnitudes

    localVariableName = GetThresholdContourHackVariableNameString(self.mVariableInfo)

    newParaViewFilter.Isosurfaces = self.mContourValue

    if PhactoriDbg():
      myDebugPrint3('  name for variable in threshold: ' + localVariableName + '\n')

    #for contour, if data type is 'element'/'CELLS', the contour filter will
    #(apparently) automatically add a cell-to-point filter
    newParaViewFilter.ContourBy = ['POINTS', localVariableName]

    if PhactoriDbg():
      myDebugPrint3("  still okay\n")
    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    #debugging contour error
    #if PhactoriDbg():
    #  newParaViewFilter.UpdatePipeline()
    #  pointData = newParaViewFilter.PointData
    #  numPointArrays = pointData.GetNumberOfArrays()
    #  myDebugPrint3("number of point arrays: " + str(numPointArrays) + "\n"
    #    "ComputeScalars: " + str(newParaViewFilter.ComputeScalars) + "\n")
    #  for ii in range(numPointArrays):
    #    onePointArray = pointData.GetArray(ii)
    #    myDebugPrint3(str(ii) + ": " + onePointArray.GetName() + ": " + \
    #        str(onePointArray.GetNumberOfTuples()) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriContourOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter
#phactori_combine_to_single_python_file_subpiece_end_1
