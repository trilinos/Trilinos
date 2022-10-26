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

global RollUpAllJsonForRunningDuringSetUpClass
RollUpAllJsonForRunningDuringSetUpClass = True

#from phactori import *
import unittest
##from paraview.simple import *
import os
import subprocess
import json
import importlib
import time

#exodus2catalyst_import = importlib.reload("exodus2catalyst4test")
#import exodus2catalyst4test

class TestPhactoriCameraBlock(unittest.TestCase):

  def MakeTestJsonSet1(self, testCam1Json, testImageBasename):
    #add camera type
    testCam1Json["camera type"] = "camera"
    #make this version of boilerplate json for test
    testJsonSetup1 = {
      "camera blocks":{
        "testcam1":testCam1Json
      },
      "representation blocks":{},
      "imageset blocks":{
        "imageset1":{
          "camera":"testcam1",
           "image basename":testImageBasename
        }
      },
      "operation blocks":{}
    }
    return testJsonSetup1

  @staticmethod
  def runit1(testJson1, testName = "testcatalystscript1"):
    #print("runit entered")
    #print("testName:")
    #print(testName)
    #print("testJson1:")
    #print(str(testJson1))
    testJsonFilename = testName + ".json"
    os.environ["PHACTORI_TEST_CAMERA_MODE"] = "on"
    os.environ["SNL_CATALYST_SIERRA_USAGE_LOG_FLAG"] = "disable"
    ff = open(testJsonFilename, "w")
    json.dump(testJson1, ff)
    ff.close()
    myenv = os.environ
    #"/home/jamauld/vbworkdisk1/paraview/paraview_install_2020Aug26/ParaView-5.8.0-MPI-Linux-Python3.7-64bit/bin/pvbatch",
    runresult1 = subprocess.run(["pvbatch",
        "exodus2catalyst4test.py", "-pj", testJsonFilename,
        "TestData/SimpleExodusForTest1.e"], env=myenv)

    #print("runit calling import_module\n")
    #exodus2catalyst_import = importlib.import_module("exodus2catalyst4test")
    ##exodus2catalyst_import = importlib.reload("exodus2catalyst4test")
    #print("runit returned from import_module, calling main2\n")
    #exodus2catalyst_import.main2(testJsonFilename, "SimpleExodusForTest1.e")
    #print("runit returned calling main2\n")
    #print("runit calling main2\n")
    #exodus2catalyst4test.main2(testJsonFilename, "SimpleExodusForTest1.e")
    #print("runit returned from calling main2\n")

    #print("runit returning")
    subprocess.run(["rm", testJsonFilename])

  @classmethod
  def setUpClass(cls):
    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass:
      testJson1 = cls.GetBaseTestJson1()

      testname, thisTestSubJson1 = cls.getjson_CameraLookDirectionWithDefaultLookAtPointZMinus1()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraLookDirectionWithDefaultLookAtPointZPlus1()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraLookDirectionWithDefaultLookAtPointXYZMinus1()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraLookAtAbsolutePoint1()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraLookAtRelativePoint1()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraLookAtElement1()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraLookAtNode1()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraLookDirection()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraLookAtRelativeDistance()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraLookAtAbsoluteDistance()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraAtAbsolutePoint()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraAtRelativePoint()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraAtElement()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraAtElementDisplaced()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraAtNodeDisplaced()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraUpVector()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraFOV()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraParallelProjection1()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)
      testname, thisTestSubJson1 = cls.getjson_CameraImageNameAddon()
      cls.MergeTestJson1(testJson1, thisTestSubJson1)

      TestPhactoriCameraBlock.runit1(testJson1, "test_PhactoriCameraBlockRolledUpTests")

  def CheckVec3AlmostEqual(self, vec1, vec2):
    for ii in range (0,3):
      self.assertAlmostEqual(vec1[ii], vec2[ii])

  def CheckCameraPostionAndFocalPoint(self, testOutputJsonFileName,
    baselineCamPos = None, baselineCamFocalPt = None,
    baselineCamViewUp = None, baselineCamViewAngle = None,
    baselineCamParallelProjection = None, baselineCamParallelScale = None):
    ff = open(testOutputJsonFileName,"r")
    testJson = json.load(ff)
    ff.close()
    camPos = testJson["CameraPosition"]
    camFocalPt = testJson["CameraFocalPoint"]
    camViewUp = testJson["CameraViewUp"]
    camViewAngle = testJson["CameraViewAngle"]
    camParallelProjection = testJson["CameraParallelProjection"]
    camParallelScale = testJson["CameraParallelScale"]
    if baselineCamPos != None:
      self.CheckVec3AlmostEqual(camPos, baselineCamPos)
    if baselineCamFocalPt != None:
      self.CheckVec3AlmostEqual(camFocalPt, baselineCamFocalPt)
    if baselineCamViewUp != None:
      self.CheckVec3AlmostEqual(baselineCamViewUp, camViewUp)
    if baselineCamViewAngle != None:
      self.assertEqual(baselineCamViewAngle, camViewAngle)
    if baselineCamParallelProjection != None:
      self.assertEqual(baselineCamParallelProjection, camParallelProjection)
    if baselineCamParallelScale != None:
      self.assertAlmostEqual(baselineCamParallelScale, camParallelScale)

  def RemoveTestOutputFiles2(self, testname):
    cmditems = []
    cmditems.append("rm")
    doRemoveImages = False
    for ii in range(0,3):
      cmditems.append(testname + ".000" + str(ii) + ".png.test.camera.txt")
      if doRemoveImages:
        cmditems.append(testname + ".000" + str(ii) + ".png")
    subprocess.run(cmditems)

  @staticmethod
  def GetBaseTestJson1():
    testJsonSetup1 = {
      "camera blocks":{},
      "representation blocks":{},
      "imageset blocks":{},
      "operation blocks":{}
    }
    return testJsonSetup1

  @staticmethod
  def MergeTestJson1(allTestsJson, testJsonToMerge):
    #print("MergeTestJson1 entered")
    #print("allTestsJson:")
    #print(str(allTestsJson))
    #print("testJsonToMerge:")
    #print(str(testJsonToMerge))
    cbjson = testJsonToMerge["camera blocks"]
    for cameraNameKey, cameraJsonValue in cbjson.items():
      allTestsJson["camera blocks"][cameraNameKey] = cameraJsonValue
    isjson = testJsonToMerge["imageset blocks"]
    for imagesetBlockNameKey, imagesetBlockJsonValue in isjson.items():
      allTestsJson["imageset blocks"][imagesetBlockNameKey] = imagesetBlockJsonValue
    #print("allTestsJson after merge:")
    #print(str(allTestsJson))
    #print("MergeTestJson1 returning")

  @staticmethod
  def getjson_ForTestHelper1(testname, thisTestJson):
    cameraName = testname + "_cam"
    imagesetName = testname + "_imageset"
    imageBasename = testname + "."
    myjson = {
      "camera blocks":{
        cameraName:thisTestJson
      },
      "imageset blocks":{
        imagesetName:{
          "image size":[800,450],
          "camera":cameraName,
          "image basename":imageBasename
        }
      }
    }
    return testname, myjson

  @staticmethod
  def getjson_CameraLookDirectionWithDefaultLookAtPointZMinus1():
    testname = "CameraLookDirectionWithDefaultLookAtPointZMinus1"
    thisTestJson = {"look direction":[0.0,0.0,-1.0],"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookDirectionWithDefaultLookAtPointZMinus1(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookDirectionWithDefaultLookAtPointZMinus1()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)
    #testCam1 =  {"look direction":[0.0,0.0,-1.0]}
    #testJson1 = self.MakeTestJsonSet1(testCam1, testImageBasename)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [5.0, -0.001999974250793457, 11.69152674091426]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname);

  @staticmethod
  def getjson_CameraLookDirectionWithDefaultLookAtPointZPlus1():
    testname = "CameraLookDirectionWithDefaultLookAtPointZPlus1"
    thisTestJson = {"look direction":[0.0,0.0,1.0],"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookDirectionWithDefaultLookAtPointZPlus1(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookDirectionWithDefaultLookAtPointZPlus1()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [5.0, -0.001999974250793457, -12.94152674091426]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraLookDirectionWithDefaultLookAtPointXYZMinus1():
    testname = "CameraLookDirectionWithDefaultLookAtPointXYZMinus1"
    thisTestJson = {"look direction":[-1.0,-1.0,-1.0],"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookDirectionWithDefaultLookAtPointXYZMinus1(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookDirectionWithDefaultLookAtPointXYZMinus1()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [15.955234275935393, 10.9532343016846, 10.330234275935393]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraLookAtAbsolutePoint1():
    testname = "CameraLookAtAbsolutePoint1"
    thisTestJson = {"look at absolute point":[1.1, 2, 3e-8],"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookAtAbsolutePoint1(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookAtAbsolutePoint1()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [20.04062128767651, 20.94062128767651, 18.94062131767651]
    baselineCamFocalPt = [1.1, 2.0, 3e-08]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraLookAtRelativePoint1():
    testname = "CameraLookAtRelativePoint1"
    thisTestJson = {"look at relative point":[1.5, 0.5, 0.5],"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookAtRelativePoint1(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookAtRelativePoint1()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [36.283676885212415, 18.279676936710832, 16.28367688521242]
    baselineCamFocalPt = [20.0, 1.996000051498413, 0.0]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraLookAtElement1():
    testname = "CameraLookAtElement1"
    thisTestJson = {"look at element":17,"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookAtElement1(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookAtElement1()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [17.57686849224281, 13.693803655979154, 12.070801931497572]
    baselineCamFocalPt = [4.881066560745239, 0.9980017244815826, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraLookAtNode1():
    testname = "CameraLookAtNode1"
    thisTestJson = {"look at node":20,"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookAtNode1(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookAtNode1()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [16.37710985121342, 9.618462409349407, 8.368462409349407]
    baselineCamFocalPt = [6.758647441864014, 0.0, -1.25]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraLookDirection():
    testname = "CameraLookDirection"
    thisTestJson = {"look direction": [1, 2, 3],"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookDirection(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookDirection()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [1.3139980597431027, -7.374003854764588, -11.683005820770692]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraLookAtRelativeDistance():
    testname = "CameraLookAtRelativeDistance"
    thisTestJson = {"look at relative distance": 2.0,"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookAtRelativeDistance(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookAtRelativeDistance()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [26.910468551870785, 21.908468577619992, 21.285468551870785]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraLookAtAbsoluteDistance():
    testname = "CameraLookAtAbsoluteDistance"
    thisTestJson = {"look at absolute distance": 15.0,"camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraLookAtAbsoluteDistance(self):
    testname, thisTestSubJson1 = self.getjson_CameraLookAtAbsoluteDistance()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [13.660254037844387, 8.658254063593594, 8.035254037844387]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraAtAbsolutePoint():
    testname = "CameraAtAbsolutePoint"
    thisTestJson = {"camera at absolute point": [-2.0, 3.0, 30.0], "camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraAtAbsolutePoint(self):
    testname, thisTestSubJson1 = self.getjson_CameraAtAbsolutePoint()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [-2.0, 3.0, 30.0]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraAtRelativePoint():
    testname = "CameraAtRelativePoint"
    thisTestJson = {"camera at relative point": [-0.5, 1.5, 20.0], "camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraAtRelativePoint(self):
    testname, thisTestSubJson1 = self.getjson_CameraAtRelativePoint()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [0.0, 5.992000102996826, 24.375]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraAtElement():
    testname = "CameraAtElement"
    thisTestJson = {"camera at element": 1, "camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraAtElement(self):
    testname, thisTestSubJson1 = self.getjson_CameraAtElement()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [0.37538909912109375, 0.3326896131038666, -0.625]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraAtElementDisplaced():
    testname = "CameraAtElementDisplaced"
    thisTestJson = {"camera at element displaced": [1,-3.0,10.0,20.0], "look at element": 1, "camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraAtElementDisplaced(self):
    testname, thisTestSubJson1 = self.getjson_CameraAtElementDisplaced()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [-2.6246109008789062, 10.332689613103867, 19.375]
    baselineCamFocalPt = [0.37538909912109375, 0.3326896131038666, -0.625]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraAtNodeDisplaced():
    testname = "CameraAtNodeDisplaced"
    thisTestJson = {"camera at node displaced": [1,-3.0,10.0,20.0], "look at node": 1, "camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraAtNodeDisplaced(self):
    testname, thisTestSubJson1 = self.getjson_CameraAtNodeDisplaced()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [-3.0, 10.0, 20.0]
    baselineCamFocalPt = [0.0, 0.0, 0.0]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraUpVector():
    testname = "CameraUpVector"
    thisTestJson = {"up vector": [0,1,2], "camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraUpVector(self):
    testname, thisTestSubJson1 = self.getjson_CameraUpVector()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [16.67302481767549, 11.671024843424698, 11.048024817675492]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    baselineCamViewUp = [0.0, 1.0, 2.0]
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt, baselineCamViewUp)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraFOV():
    testname = "CameraFOV"
    thisTestJson = {"camera fov": 45,
      "look at relative point": [0.0, 0.0, 0.0],
      "look at absolute distance": 10.0,
      "camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraFOV(self):
    testname, thisTestSubJson1 = self.getjson_CameraFOV()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [10.773502691896258, 5.771502717645465, 5.148502691896258]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    baselineCamAngle = 45.0
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos, baselineCamFocalPt, None, baselineCamAngle)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraParallelProjection1():
    testname = "CameraParallelProjection1"
    thisTestJson = {"projection type": "parallel", "look direction": [-5,-1,-1], "camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraParallelProjection1(self):
    testname, thisTestSubJson1 = self.getjson_CameraParallelProjection1()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    baselineCamPos = [21.08689812331651, 3.2153796504125087, 2.592379624663302]
    baselineCamFocalPt = [5.0, -0.001999974250793457, -0.625]
    baselineCamParallelProjection = 1
    baselineCamParallelScale = 5.69160334053072
    self.CheckCameraPostionAndFocalPoint(testOutFile, baselineCamPos,
      baselineCamFocalPt, None, None,
      baselineCamParallelProjection, baselineCamParallelScale)
    self.RemoveTestOutputFiles2(testname)

  @staticmethod
  def getjson_CameraImageNameAddon():
    testname = "CameraImageNameAddon"
    thisTestJson = {"image name addon": "_foo_", "camera type":"camera"}
    return TestPhactoriCameraBlock.getjson_ForTestHelper1(testname, thisTestJson)

  def test_CameraImageNameAddon(self):
    testname, thisTestSubJson1 = self.getjson_CameraImageNameAddon()
    testOutFileEnding = ".png.test.camera.txt"
    testOutFile = testname + ".0002" + testOutFileEnding

    testJson1 = self.GetBaseTestJson1()
    self.MergeTestJson1(testJson1, thisTestSubJson1)

    global RollUpAllJsonForRunningDuringSetUpClass
    if RollUpAllJsonForRunningDuringSetUpClass == False:
      self.runit1(testJson1, testname)

    #name of image file should be different, so test for that.
    imageNameGotAddon = os.path.exists("CameraImageNameAddon._foo_0002.png.test.camera.txt")
    self.assertTrue(imageNameGotAddon)
    subprocess.run(["rm",
      "CameraImageNameAddon._foo_0000.png.test.camera.txt",
      "CameraImageNameAddon._foo_0001.png.test.camera.txt",
      "CameraImageNameAddon._foo_0002.png.test.camera.txt",])


if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()
