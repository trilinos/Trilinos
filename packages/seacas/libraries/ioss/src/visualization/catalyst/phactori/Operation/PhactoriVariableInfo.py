# Copyright(C) 1999-2022 National Technology & Engineering Solutions
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
class PhactoriVariableMinMaxAvgSumCntSave:
  def __init__(self):
    self.mStatsTestCounter = -1
    self.mMin = 0.0
    self.mMax = 0.0
    self.mSum = 0.0
    self.mCount = 0
    self.mIdsTestCounter = -1
    self.mLocalFoundMinId = False
    self.mLocalMinIdCount = -1
    self.mMinId = -1;
    self.mLocalFoundMaxId = False
    self.mLocalMaxIdCount = -1
    self.mMaxId = -1;

class PhactoriVariableInfo:
  def __init__(self):
    self.mVariableIntendedForUseFlag = False
    self.mVariableName = ""
    self.mVariableComponent = None
    self.mVariableType = None
    self.mVariableTypeNeedsDetection = False
    self.mVariableTypeWasDetected = False
    self.mVariableTypeWasCopied = False
    self.mVariableIsVectorMagnitude = False
    self.mVariableIsVectorComponent = False
    self.mVectorBaseName = None
    self.mStats = PhactoriVariableMinMaxAvgSumCntSave()
    #as of 2022Aug09 we want this to default to False
    self.mAddSeparatorToVectorVariableName = False

  def SelfToStr(self):
    retStr = "PhactoriVariableInfo:\n" +\
      "\n  mVariableName: " + str(self.mVariableName) +\
      "\n  mVariableComponent: " + str(self.mVariableComponent) +\
      "\n  mVariableType: " + str(self.mVariableType) +\
      "\n  mVariableTypeNeedsDetection: " + str(self.mVariableTypeNeedsDetection) +\
      "\n  mVariableTypeWasDetected: " + str(self.mVariableTypeWasDetected) +\
      "\n  mVariableTypeWasCopied: " + str(self.mVariableTypeWasCopied) +\
      "\n  mVariableIsVectorMagnitude: " + \
          str(self.mVariableIsVectorMagnitude) +\
      "\n  mVariableIsVectorComponent: " + \
          str(self.mVariableIsVectorComponent) +\
      "\n  mVectorBaseName: " + str(self.mVectorBaseName) +\
      "\n  mVariableIntendedForUseFlag: " + str(self.mVariableIntendedForUseFlag) +\
      "\n"
    return retStr

  def CheckAndParseScalarVariable1(self, inJson, inBaseString):
    testKey = inBaseString + 'scalar'
    if testKey in inJson:
      self.mVariableName = inJson[testKey]
      self.mVariableComponent = None
      return True
    else:
      return False

  def CheckAndParseVectorComponentNumber1(self, inJson, inBaseString):
    testKey = inBaseString + 'vector component number'
    if testKey in inJson:
      self.mVariableIsVectorComponent = True
      varCompNumPair = inJson[testKey]
      self.mVariableName = varCompNumPair[0]
      self.mVariableComponent = int(varCompNumPair[1])
      return True
    else:
      return False

  def CheckAndParseVectorComponentSuffix(self, inJson, inBaseString):
    testKey = inBaseString + 'vector component'
    if testKey in inJson:
      self.mVariableIsVectorComponent = True
      varName = inJson[testKey]
      resultPair = breakSpecialVarNameIntoBaseAndComponent(varName,
                      self.mAddSeparatorToVectorVariableName)
      self.mVariableName = resultPair[0]
      self.mVariableComponent = resultPair[1]
      self.mVectorBaseName = resultPair[2]
      return True
    else:
      return False

  def CheckAndParseVectorMagnitude(self, inJson, inBaseString):
    testKey = inBaseString + 'vector magnitude'
    if testKey in inJson:
      self.mVariableIsVectorMagnitude = True
      if self.mAddSeparatorToVectorVariableName:
        self.mVariableName = inJson[testKey] + GetSeparatorString()
      else:
        self.mVariableName = inJson[testKey]
      self.mVectorBaseName = inJson[testKey]
      return True
    else:
      return False

  def CheckAndParseTensorComponent(self, inJson, inBaseString):
    testKey = inBaseString + 'tensor component'
    if testKey in inJson:
      self.mVariableName = inJson[testKey]
      return True
    else:
      return False

  def CheckAndParseVectorVariable1(self, inJson, inBaseString):
    variableFoundFlag = False
    variableFoundFlag = self.CheckAndParseVectorComponentNumber1(inJson, inBaseString)
    if not variableFoundFlag:
      variableFoundFlag = self.CheckAndParseVectorComponentSuffix(inJson, inBaseString)
      if not variableFoundFlag:
        variableFoundFlag = self.CheckAndParseVectorMagnitude(inJson, inBaseString)
        if not variableFoundFlag:
          variableFoundFlag = self.CheckAndParseTensorComponent(inJson, inBaseString)
          if not variableFoundFlag:
            variableFoundFlag = False
            self.mVariableName = ''
            self.mComponent = None
    return variableFoundFlag

  def ParseVariableNameAndVectorOrTensorComponent(self, inJson, inBaseString):
    "take a base string such as 'y axis variable ' or 'variable ', use it\n    to construct keys to define a scalar variable, vector component, vector\n    magnitude, or tensor component, and see if the json has those keys.\n    If the json has a key, use it to grab the variable name and setup.\n    Also look to see if the type of the variable is specifically defined\n    (node or element) or needs to be detected\n    "
    variableFoundFlag = False
    self.mVariableIsVectorComponent = False
    self.mVariableIsVectorMagnitude = False

    if 'add separator to vector variable name' in inJson:
      self.mAddSeparatorToVectorVariableName = \
        inJson['add separator to vector variable name']

    variableFoundFlag = self.CheckAndParseScalarVariable1(inJson, inBaseString)
    if not variableFoundFlag:
      variableFoundFlag = self.CheckAndParseVectorVariable1(inJson, inBaseString)
    if variableFoundFlag:
      #it is now apparent variable is intended to be used, not ignored
      self.mVariableIntendedForUseFlag = True
      if 'variable type' in inJson:
        variableType = inJson['variable type']
        self.mVariableType = variableType
        if variableType != 'element' \
            and variableType != 'node' \
            and variableType != 'global' \
            and variableType != 'detect':
          errStr = 'error!  inJson has variable type is neither node nor '\
              'element nor detect nor global in ParseVariableNameAndVectorOrTensorComponent\n'
          if PhactoriDbg():
            myDebugPrint3(errStr)
          raise Exception(errStr)
        if variableType == 'detect':
          self.mVariableType = 'element'
          self.mVariableTypeNeedsDetection = True
          self.mVariableTypeWasDetected = False
        elif variableType == 'element':
          self.mVariableType = 'element'
          self.mVariableTypeNeedsDetection = False
          self.mVariableTypeWasDetected = True
        elif variableType == 'node':
          self.mVariableType = 'node'
          self.mVariableTypeNeedsDetection = False
          self.mVariableTypeWasDetected = True
        elif variableType == 'global':
          self.mVariableType = 'global'
          self.mVariableTypeNeedsDetection = False
          self.mVariableTypeWasDetected = True
      else:
        self.mVariableType = 'element'
        self.mVariableTypeNeedsDetection = True
        self.mVariableTypeWasDetected = False
    else:
      self.mVariableType = None

    if PhactoriDbg():
      str1 = "ParseVariableNameAndVectorOrTensorComponent:\n"
      str2 = "inBaseString:\n" + str(inBaseString) + "\n"
      str3 = "inJson:\n" + str(inJson) + "\n"
      str4 = "instance state after parsing:\n" + self.SelfToStr()
      myDebugPrint3(str1 + str2 + str3 + str4)
    return variableFoundFlag

  def CopyVariableTypeFrom(self, inSourceVariableInfo):
    self.mVariableType = inSourceVariableInfo.mVariableType
    self.mVariableTypeNeedsDetection = False
    self.mVariableTypeWasDetected = True
    self.mVariableTypeWasCopied = True

  def DetectVariableType(self, inInputCsData, inInputIsProxy = False,
          inAllowGlobalVariableDetection = False):

    if self.mVariableTypeNeedsDetection == False:
      return True

    #testing hack--force to cell
    #self.mVariableTypeNeedsDetection = False
    #self.mVariableTypeWasDetected = True
    #self.mVariableType = 'element'
    #return

    if PhactoriDbg():
      myDebugPrint3('PhactoriVariableInfo::DetectVariableType entered\n')

    if inInputIsProxy:
      pointData = inInputCsData.PointData
    else:
      pointData = inInputCsData.GetPointData()
    if pointData != None:
      testPointArray = pointData.GetArray(self.mVariableName)
      if testPointArray != None:
        if PhactoriDbg():
          myDebugPrint3('  type node detected!\n')
        self.mVariableType = 'node'
        self.mVariableTypeNeedsDetection = False
        self.mVariableTypeWasDetected = True
        return True

    if inInputIsProxy:
      cellData = inInputCsData.CellData
    else:
      cellData = inInputCsData.GetCellData()
    if cellData != None:
      testCellArray = cellData.GetArray(self.mVariableName)
      if testCellArray != None:
        if PhactoriDbg():
          myDebugPrint3('  type element detected!\n')
        self.mVariableType = 'element'
        self.mVariableTypeNeedsDetection = False
        self.mVariableTypeWasDetected = True
        return True

    if inAllowGlobalVariableDetection:
      if inInputIsProxy:
        fieldData = inInputCsData.FieldData
      else:
        fieldData = inInputCsData.GetFieldData()
      if fieldData != None:
        testFieldArray = fieldData.GetArray(self.mVariableName)
        if testFieldArray != None:
          if PhactoriDbg():
            myDebugPrint3('  type global detected!\n')
          self.mVariableType = 'global'
          self.mVariableTypeNeedsDetection = False
          self.mVariableTypeWasDetected = True
          return True

    if PhactoriDbg():
      myDebugPrint3('  type not detected!\n')
      #default to 'element' type knowing assumption may be wrong,
      #leave mVariableTypeNeedsDetection set to True so we know the variable
      #type has not yet been successfully detected
      self.mVariableType = 'element'
      #self.mVariableType = 'node'
    return False

  def GetXYPlotAxisTitle(self):
    """for this variable info, construct the name to put on the axis of an XY
    plot, e.g. with magnitude or component information"""

    if self.mVariableName[-1] == '_':
      axisTitle = str(self.mVariableName[0:-1])
    else:
      axisTitle = self.mVariableName

    if self.mVariableIsVectorComponent:
      if self.mVariableComponent == 0:
        axisTitle += ' X component'
      elif self.mVariableComponent == 1:
        axisTitle += ' Y component'
      elif self.mVariableComponent == 2:
        axisTitle += ' Z component'
      else:
        if PhactoriDbg():
          myDebugPrint3("  variable component is not 0, 1, or 2, using 0 (X)\n")
        localVariableName = self.mVariableName + '_X'
    elif self.mVariableIsVectorMagnitude:
      axisTitle += " Magnitude"

    return axisTitle

#phactori_combine_to_single_python_file_subpiece_end_1
