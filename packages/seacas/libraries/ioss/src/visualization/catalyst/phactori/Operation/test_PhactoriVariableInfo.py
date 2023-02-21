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
import unittest
from paraview.simple import *

class Test_PhactoriVariableInfo(unittest.TestCase):

  def test_ScalarVariableParsing(self):
    testVarInfo = PhactoriVariableInfo()
    testJson = {"color by scalar": "temperature"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "temperature")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, True)
    #defaults to variable type of element if needing detection
    self.assertEqual(testVarInfo.mVariableType, "element")
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, False)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, False)

  def test_ScalarVariableWithCellOrNodeIdentificationParsing(self):
    testVarInfo = PhactoriVariableInfo()
    testJson = {"color by scalar": "pressure", "variable type": "element"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "pressure")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, False)
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, True)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, False)
    self.assertEqual(testVarInfo.mVariableType, "element")
    testVarInfo = PhactoriVariableInfo()
    testJson = {"color by scalar": "temperature", "variable type": "node"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "temperature")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, False)
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, True)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, False)
    self.assertEqual(testVarInfo.mVariableType, "node")
    testJson = {"color by scalar": "speed", "variable type": "global"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "speed")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, False)
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, True)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, False)
    self.assertEqual(testVarInfo.mVariableType, "global")
    testJson = {"color by scalar": "temperature", "variable type": "detect"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "temperature")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, True)
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, False)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, False)
    self.assertEqual(testVarInfo.mVariableType, "element")

  def test_VectorMagnitudeVariableParsing(self):
    testVarInfo = PhactoriVariableInfo()
    testJson = {"color by vector magnitude": "velocity"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "velocity")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, True)
    #defaults to variable type of element if needing detection
    self.assertEqual(testVarInfo.mVariableType, "element")
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, False)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, True)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, False)

  def test_VectorComponentViaLetterSuffixVariableParsing1(self):
    testVarInfo = PhactoriVariableInfo()
    testJson = {"color by vector component": "velocityX"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "velocity")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, True)
    #defaults to variable type of element if needing detection
    self.assertEqual(testVarInfo.mVariableType, "element")
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, False)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, True)
    self.assertEqual(testVarInfo.mVariableComponent, 0)
    testJson = {"color by vector component": "velocityY"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableComponent, 1)
    testJson = {"color by vector component": "velocityZ"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableComponent, 2)
    testJson = {"color by vector component": "velocityY"}

  def test_VectorComponentViaNumberSuffixVariableParsing1(self):
    testVarInfo = PhactoriVariableInfo()
    testJson = {"color by vector component": "velocity0"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "velocity")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, True)
    #defaults to variable type of element if needing detection
    self.assertEqual(testVarInfo.mVariableType, "element")
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, False)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, True)
    self.assertEqual(testVarInfo.mVariableComponent, 0)
    testJson = {"color by vector component": "velocity2"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableComponent, 2)
    testJson = {"color by vector component": "velocity10"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableComponent, 10)
    testJson = {"color by vector component": "velocity11"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableComponent, 11)
    testJson = {"color by vector component": "velocity99"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableComponent, 99)

  def test_VectorComponentViaLetterSuffixVariableWithSeparatorParsing1(self):
    testVarInfo = PhactoriVariableInfo()
    testVarInfo.mAddSeparatorToVectorVariableName = True
    newRoot = PhactoriPipeAndViewsState()
    SetgPipeAndViewStateForTesting(newRoot)
    testJson = {"color by vector component": "velocity_Y"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "velocity_")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, True)
    #defaults to variable type of element if needing detection
    self.assertEqual(testVarInfo.mVariableType, "element")
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, False)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, True)
    self.assertEqual(testVarInfo.mVariableComponent, 1)
    testVarInfo = PhactoriVariableInfo()
    testJson = {"color by vector magnitude": "velocity"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "velocity")
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, True)

  def test_VectorComponentViaNumberParsing1(self):
    testVarInfo = PhactoriVariableInfo()
    testJson = {"color by vector component number": ["velocity", 0] }
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "velocity")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, True)
    #defaults to variable type of element if needing detection
    self.assertEqual(testVarInfo.mVariableType, "element")
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, False)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, True)
    self.assertEqual(testVarInfo.mVariableComponent, 0)

    testJson = {"variable vector component number": ["velocity", 1] }
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "variable ")
    self.assertEqual(testVarInfo.mVariableComponent, 1)

    testJson = {"variable vector component number": ["velocity", 77] }
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "variable ")
    self.assertEqual(testVarInfo.mVariableComponent, 77)

  def test_ScalarVariableWithCellOrNodeIdentificationParsing(self):
    testVarInfo = PhactoriVariableInfo()
    testJson = {"color by vector component number": ["velocity", 2], "variable type":"node"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableIntendedForUseFlag, True)
    self.assertEqual(testVarInfo.mVariableName, "velocity")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, False)
    self.assertEqual(testVarInfo.mVariableType, "node")
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, True)
    self.assertEqual(testVarInfo.mVariableIsVectorMagnitude, False)
    self.assertEqual(testVarInfo.mVariableIsVectorComponent, True)
    self.assertEqual(testVarInfo.mVariableComponent, 2)
    testJson = {"color by vector component number": ["velocity", 3], "variable type":"element"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableType, "element")
    self.assertEqual(testVarInfo.mVariableComponent, 3)
    testJson = {"color by vector component number": ["velocity", 0], "variable type":"global"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableType, "global")
    self.assertEqual(testVarInfo.mVariableComponent, 0)
    testJson = {"color by vector component number": ["velocity", 7], "variable type":"detect"}
    testVarInfo.ParseVariableNameAndVectorOrTensorComponent(testJson, "color by ")
    self.assertEqual(testVarInfo.mVariableTypeNeedsDetection, True)
    self.assertEqual(testVarInfo.mVariableTypeWasDetected, False)
    #defaults to variable type of element if needing detection
    self.assertEqual(testVarInfo.mVariableType, "element")
    self.assertEqual(testVarInfo.mVariableComponent, 7)

  #def __init__(self):
  #  self.mVariableIntendedForUseFlag = False
  #  self.mVariableName = ""
  #  self.mVariableComponent = None
  #  self.mVariableType = None
  #  self.mVariableTypeNeedsDetection = False
  #  self.mVariableTypeWasDetected = False
  #  self.mVariableTypeWasCopied = False
  #  self.mVariableIsVectorMagnitude = False
  #  self.mVariableIsVectorComponent = False
  #  self.mVectorBaseName = None
  #  self.mStats = PhactoriVariableMinMaxAvgSumCntSave()
  #  #as of 2022Aug09 we want this to default to False
  #  self.mAddSeparatorToVectorVariableName = False

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


