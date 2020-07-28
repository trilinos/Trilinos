# Copyright(C) 1999-2020 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details







#these settings allow us to specify a particular json format (or bccolli.py)
#text file to use as the phactori visualization setup file.
#if gUseHardwiredJsonScriptFile is None, operation is normal, but if it
#is set to a file string, we ignore any other parsing and use the specified
#file
global gUseHardwiredJsonScriptFile
gUseHardwiredJsonScriptFile = None
#gUseHardwiredJsonScriptFile = "myscript1.json"

global WriteEachDeadCellElementToFiles
WriteEachDeadCellElementToFiles = False
global WriteDeadCellSummaryFile
WriteDeadCellSummaryFile = True

global gPipeAndViewsState
gPipeAndViewsState = None

try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing
import numpy as np

import vtk

global gParaViewCatalystVersionFlag
#gParaViewCatalystVersionFlag = 40300
#gParaViewCatalystVersionFlag = 40100
#gParaViewCatalystVersionFlag = 50400
gParaViewCatalystVersionFlag = 50502

global gPointsString
if gParaViewCatalystVersionFlag <= 40100:
  gPointsString = 'POINT_DATA'
  gCellsString = 'CELL_DATA'
else:
  gPointsString = 'POINTS'
  gCellsString = 'CELLS'

import math
import sys
import datetime

global gBypassUserData
#gBypassUserData = True
gBypassUserData = False

global gUserDataBypassStrings
gUserDataBypassStrings = None

global gBypassUserDataJson
gBypassUserDataJson = None

global TestUserDataForBypassScriptCompletedFlag
TestUserDataForBypassScriptCompletedFlag = False

global bccolli_controls

def TestUserDataForBypassScript(datadescription):

  global TestUserDataForBypassScriptCompletedFlag
  if TestUserDataForBypassScriptCompletedFlag == True:
    #already tested before, just return
    return

  TestUserDataForBypassScriptCompletedFlag = True

  global gBypassUserData
  if gBypassUserData:
    #already true, just return
    return

  if PhactoriDbg(100):
    myDebugPrint3("TestUserDataForBypassScript entered\n")

  global gUseHardwiredJsonScriptFile
  if gUseHardwiredJsonScriptFile == None:
    fd = datadescription.GetUserData()

    if fd == None:
      if PhactoriDbg(100):
        myDebugPrint3("no user data, no catalyst_script_extra_file\n")
      return

    sa = fd.GetAbstractArray(0)

    if sa == None:
      if PhactoriDbg(100):
        myDebugPrint3("no user data, no catalyst_script_extra_file (2)\n")
      return

    if(sa.GetNumberOfValues() > 8):
      catalyst_script_extra_file = sa.GetValue(8)
      if PhactoriDbg(100):
        myDebugPrint3("  catalyst_script_extra_file: ->" + \
                catalyst_script_extra_file + "<-\n")
    else:
      if PhactoriDbg(100):
        myDebugPrint3("  catalyst_script_extra_file: ->NONE<- (no bypass)\n")
      return

    if catalyst_script_extra_file == "":
      if PhactoriDbg(100):
        myDebugPrint3("empty catalyst_script_extra_file, returning\n")
      return;
  else:
    catalyst_script_extra_file = gUseHardwiredJsonScriptFile

  global gUserDataBypassStrings
  global gBypassUserDataJson

  if catalyst_script_extra_file == "bccolli_controls.py":
    #due to unique name for catalyst_script_extra_file, do python import
    #to obtain json description of visualization

    if PhactoriDbg(100):
      myDebugPrint3("using bypass script bccolli_controls.py to specify\n" \
              "vis instead of catalyst block\n")
    gBypassUserData = True
    global bccolli_controls
    import bccolli_controls
    gUserDataBypassStrings = [
      bccolli_controls.catalystSierraInputInJsonFormatStr,
      "_",
      "simple",
      "",
      "",
      "0",
      "simple.e",
      "",
      catalyst_script_extra_file
    ]
    gBypassUserDataJson = bccolli_controls.catalystSierraInputInJsonFormat
  else:
    #since it's not the unique name bccolli_controls, treat it as a text file
    #containing a json format data set.  Read in this json and use it as
    #the phactori json format visualization description

    if PhactoriDbg(100):
        myDebugPrint3("using json data file: " + catalyst_script_extra_file + \
              "\nto specify vis instead of catalyst block\n")
    gBypassUserData = True

    import json

    try:
      with open(catalyst_script_extra_file) as data_file:
        if PhactoriDbg(100):
          myDebugPrint3("open successful: \n")
        visDescpJson = json.load(data_file)
        visDescpJson = convertJsonUnicodeToStrings(visDescpJson)
      if PhactoriDbg(100):
        myDebugPrint3("load successful:\n" + str(visDescpJson) + "\n")
      gUserDataBypassStrings = [
        str(visDescpJson),
        "_",
        "simple",
        "",
        "",
        "0",
        "simple.e",
        "",
        catalyst_script_extra_file
      ]
      gBypassUserDataJson =  visDescpJson
    except:
      import traceback
      excpt_str = "exception while loading json from file: " + \
              catalyst_script_extra_file + ":\n" +  \
              traceback.format_exc() + "\n"
      if PhactoriDbg(100):
        myDebugPrint3("exception while loading json from file: " + \
          catalyst_script_extra_file + ":\n" + excpt_str)
      myDebugPrint3AndException(excpt_str)


def GetBypassUserDataFlag():
  global gBypassUserData
  return gBypassUserData


#test at creation time if this global was preset (above) to bypass user data
#if so, read in info by importing bccolli_controls.py
#this is used if we want to hardwire this phactori.py file to read from
#bccolli_controls.py, rather than setting from the input deck
if gBypassUserData:
  global bccolli_controls
  import bccolli_controls
  myBypassJsonString = bccolli_controls.catalystSierraInputInJsonFormatStr
  gBypassUserDataJson = bccolli_controls.catalystSierraInputInJsonFormat
else:
  myBypassJsonString = """{"camera blocks":{},"representation blocks":{},"imageset blocks":{"is1":{"color by scalar": "mixture_fraction","operation":"fooOperation"},"is2":{"color by scalar": "density_ra","operation":"contour1"}},"operation blocks":{"fooOperation":{"type":"clip","plane normal":[0.0,0.0,1.0],"relative point on plane":[0.0,0.0,0.0]},"contour1":{"type":"contour","variable scalar":"mixture_fraction","contour value":[0.7]}},"scatter plot blocks":{},"plot over time blocks":{}}"""

gUserDataBypassStrings = [
  myBypassJsonString,
  "_",
  "simple",
  "",
  "",
  "0",
  "simple.e",
  "",
  ""
]

def SetViewMapCBypassJsonString(inViewMapCString):
  global gUserDataBypassStrings
  gUserDataBypassStrings[0] = inViewMapCString
  if PhactoriDbg():
    myDebugPrint3("bypass json view map string changed to:\n" + gUserDataBypassStrings[0] + "\n")

global gEnableTemporaryExtractBlockTest
global gExtractBlockAList
global gExtractBlockBList
global gExtractBlockAOperationName
global gExtractBlockBOperationName
gEnableTemporaryExtractBlockTest = True
gExtractBlockAList = ['block_1']
gExtractBlockBList = ['block_2']
gExtractBlockAOperationName = "ReplaceWithExtractBlockA_Test"
gExtractBlockBOperationName = "ReplaceWithExtractBlockB_Test"

global gEnableTemporaryMultiOpViewTest
global gMultiOpViewTestOp2Name
global gMultiOpViewTestRep2Name
global gMultiOpViewTestImagesetName
gEnableTemporaryMultiOpViewTest = True
gMultiOpViewTestOp2Name = "ReplaceWithExtractBlockB_Test"
gMultiOpViewTestRep2Name = "repB"
gMultiOpViewTestImagesetName = "multioptestimageset"

class PhactoriPhysicalEyeAndScreenSetup:
  """contains physical dimensions of screen (in whatever screen units)
     and location of between-the-eyes position relative to the screen
     in screen units--(0.0, 0.0) is right at middle of screen and distance
     from the betwee-the-eyes position and the screen"""

  def __init__(self):
    """inScreenWidth, inScreenHeight obvious, in whatever units (feet, meters,
       inches); inIpdInScreenUnits interpupilary distance in same units as
       screen width and height, inEyeToScreenDistance distance from the
       between-the-eye position to the screen, in the same units;
       inBetweenEyeX, inBetweenEyeY where the eye is looking straight ahead
       at the screen, with (0.0, 0.0) being the center and 0.5 * width and
       height being at the corners"""

    #we operates by setting eye at (0,0,0) and setting the screen corners
    #appropriately
    self.mLeftEyeScreenBottomLeft = None
    self.mLeftEyeScreenBottomRight = None
    self.mLeftEyeScreenTopRight = None
    self.mRightEyeScreenBottomLeft = None
    self.mRightEyeScreenBottomRight = None
    self.mRightEyeScreenTopRight = None

  def ParseSettingsFromJson(self, inJsn):
    jsnGoodFlag = True
    if 'screen size' not in inJsn:
      jsnGoodFlag = False
    if 'screen ipd' not in inJsn:
      jsnGoodFlag = False
    if 'eye to screen distance' not in inJsn:
      jsnGoodFlag = False
    if 'between eyes projected screen position' not in inJsn:
      jsnGoodFlag = False
    if jsnGoodFlag == False:
      myDebugPrint3AndException("PhactoriPhysicalEyeAndScreenSetup:"
        "ParseSettingsFromJson:\n"
        "incoming json needs all of 'screen size' 'screen ipd'\n"
        "'eye to screen distance' and "
        "'between eyes projected screen position'\n")

    screenSize = inJsn['screen size']
    halfIpd = inJsn['screen ipd'] * 0.5
    halfW = screenSize[0] * 0.5
    halfH = screenSize[1] * 0.5
    eyeToScreenDistance = inJsn['eye to screen distance']
    eyePositionOnScreen = inJsn['between eyes projected screen position']
    betweenEyeX = eyePositionOnScreen[0]
    betweenEyeY = eyePositionOnScreen[1]
    self.mLeftEyeScreenBottomLeft = [
        -halfW + halfIpd - betweenEyeX,
        -halfH - betweenEyeY,
        -eyeToScreenDistance
      ]
    self.mLeftEyeScreenBottomRight = [
         halfW + halfIpd - betweenEyeX,
        -halfH - betweenEyeY,
        -eyeToScreenDistance
      ]
    self.mLeftEyeScreenTopRight = [
         halfW + halfIpd - betweenEyeX,
         halfH - betweenEyeY,
        -eyeToScreenDistance
      ]
    self.mRightEyeScreenBottomLeft = [
        -halfW - halfIpd - betweenEyeX,
        -halfH - betweenEyeY,
        -eyeToScreenDistance
      ]
    self.mRightEyeScreenBottomRight = [
         halfW - halfIpd - betweenEyeX,
        -halfH - betweenEyeY,
        -eyeToScreenDistance
      ]
    self.mRightEyeScreenTopRight = [
         halfW - halfIpd - betweenEyeX,
         halfH - betweenEyeY,
        -eyeToScreenDistance
      ]

  def GetScreenBottomLeft(self, inLeftEyeFlag):
    if inLeftEyeFlag:
      return self.mLeftEyeScreenBottomLeft
    else:
      return self.mRightEyeScreenBottomLeft

  def GetScreenBottomRight(self, inLeftEyeFlag):
    if inLeftEyeFlag:
      return self.mLeftEyeScreenBottomRight
    else:
      return self.mRightEyeScreenBottomRight

  def GetScreenTopRight(self, inLeftEyeFlag):
    if inLeftEyeFlag:
      return self.mLeftEyeScreenTopRight
    else:
      return self.mRightEyeScreenTopRight

  def GetEyePosition(self, inLeftEyeFlag):
    return [0.0, 0.0, 0.0]

global g1000x1000in1920x1080desktop30inch
#multiply factor to get, e.g. 800x800
global gPixFactor
gPixFactor = 0.8
#units are centimeters
#myXFactor = 10.0
myXFactor = 1.0
g1000x1000in1920x1080desktop30inch = {
  'screen size' : [33.866*myXFactor*gPixFactor, 33.866*myXFactor*gPixFactor],
  'screen ipd' : 6.2*myXFactor,
  'eye to screen distance': 55.88*myXFactor,
  'between eyes projected screen position': [0.0, 0.0],
  'ipd in model units': 6.2,
  'virtual self size multiplier': 0.25
  #'virtual self size multiplier': 1.0
}

global g1008x1792in1920x1200desktop30inch
screendiag1 = 30.0
screenaspectratio1 = 9.0/16.0
screendiagangle1 = math.atan(screenaspectratio1)
screenheight1 = screendiag1 * math.sin(screendiagangle1)
screenwidth1 = screendiag1 * math.cos(screendiagangle1)
screenxpixels1 = 1920.0
screenypixels1 = 1200.0
imagexpixels1 = 1792.0
imageypixels1 = 1008.0
imagewidth1 = (imagexpixels1 / screenxpixels1) * screenwidth1
imageheight1 = (imageypixels1 / screenypixels1) * screenheight1
#print "g1008x1792in1920x1200desktop30inch calculcations:"
#print "screendiag1: " + str(screendiag1)
#print "screenwidth1: " + str(screenwidth1)
#print "screenheight1: " + str(screenheight1)
#print "imagewidth1: " + str(imagewidth1)
#print "imageheight1: " + str(imageheight1)

g1008x1792in1920x1200desktop30inch = {
  'screen size' : [imagewidth1, imageheight1],
  'screen ipd' : 6.2,
  'eye to screen distance': 55.88,
  'between eyes projected screen position': [0.0, 0.0],
  'ipd in model units': 6.2,
  #'virtual self size multiplier': 1.0,
  'virtual self size multiplier': 0.001,
  'auto size 1': True,
  'auto size 1 distance ratio': 1.0,
  'auto size 1 view angle delta in degrees': 0.0
}


global gWiff9x16FromAudienceBig
#units are in feet
gTestFactorX = 1.0
gWiff9x16FromAudienceBig = {
  'screen size' : [15.33333*gTestFactorX, 8.75*gTestFactorX],
  'screen ipd' : 0.2034121*gTestFactorX,
  #screen to chair back is 22.41666667 feet, knock off .4 for head width
  'eye to screen distance': 22.0*gTestFactorX,
  #bottom of screen 3 feet above floor
  #eyes in front row 3.75 feet above floor
  #bottom of screen at -4.375 feet--below center
  'between eyes projected screen position': [0.0*gTestFactorX, (3.75 - 3.0 - 4.375)*gTestFactorX],
  #'ipd in model units': 0.2034121*gTestFactorX,
  'ipd in model units': 0.0,
  #'virtual self size multiplier': 0.25
  #'virtual self size multiplier': 0.001
  'virtual self size multiplier': 1.0
}

global gWiff9x16FromAudienceSmall
#units are in feet
gWiff9x16FromAudienceSmall = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 22.0,
  'between eyes projected screen position': [0.0, (3.75 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  'virtual self size multiplier': 1.0
}

global gWiff9x16FromAudienceAutoSizeA
#units are in feet
gWiff9x16FromAudienceAutoSizeA = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 22.0,
  'between eyes projected screen position': [0.0, (3.75 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  #'virtual self size multiplier': 0.001,
  'virtual self size multiplier': 0.25,
  #'virtual self size multiplier': 2.0,
  #'auto size 1 distance ratio': 1.0,
  #'auto size 1 distance ratio': 0.75,
  'auto size 1 distance ratio': 0.7,
  #'auto size 1 distance ratio': 0.001,
  #'auto size 1 distance ratio': 0.25,
  'auto size 1': True,
  'auto size 1 lock nth eye position': 1,
  'auto size 1 view angle delta in degrees': 11.0
  #'auto size 1 view angle delta in degrees': 9.0
  #'auto size 1 view angle delta in degrees': 0.0
}
global gWiff9x16FromAudienceAutoSizeA2
#units are in feet
gWiff9x16FromAudienceAutoSizeA2 = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 22.0,
  'between eyes projected screen position': [0.0, (3.75 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  #'virtual self size multiplier': 0.25,
  #'virtual self size multiplier': 0.5,
  'virtual self size multiplier': 0.4,
  #'auto size 1 distance ratio': 0.7,
  'auto size 1 distance ratio': 0.85,
  #'auto size 1 distance ratio': 0.35,
  'auto size 1': True,
  'auto size 1 lock nth eye position': 1,
  'auto size 1 view angle delta in degrees': 11.0
  #'auto size 1 view angle delta in degrees': 9.0
  #'auto size 1 view angle delta in degrees': 0.0
}


global gWiff9x16FromAudienceAutoSizeB
#units are in feet
gWiff9x16FromAudienceAutoSizeB = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 22.0,
  'between eyes projected screen position': [0.0, (3.75 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  #'virtual self size multiplier': 1.0,
  #'auto size 1 distance ratio': 1.0,
  #'virtual self size multiplier': 0.25,
  #'virtual self size multiplier': 0.125,
  #'virtual self size multiplier': 1.0,
  #'virtual self size multiplier': 1.5,
  #'virtual self size multiplier': 2.5,
  'virtual self size multiplier': 1.75,
  #'auto size 1 distance ratio': 1.0,
  #'auto size 1 distance ratio': 0.5,
  'auto size 1 distance ratio': 0.125,
  'auto size 1 lock nth eye position': 1,
  'auto size 1': True,
  #'auto size 1 view angle delta in degrees': 9.0
  #'auto size 1 view angle delta in degrees': 9.0
  'auto size 1 view angle delta in degrees': 9.0
  #'auto size 1 view angle delta in degrees': 11.0
}

#this one has locked eye position
global gWiff9x16FromAudienceAutoSizeC
#units are in feet
gWiff9x16FromAudienceAutoSizeC = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 22.0,
  'between eyes projected screen position': [0.0, (3.75 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  #'virtual self size multiplier': 1.0,
  #'auto size 1 distance ratio': 0.6,
  'virtual self size multiplier': 2.0,
  #'auto size 1 distance ratio': 0.25,
  #'auto size 1 distance ratio': 0.125,
  'auto size 1 distance ratio': 0.1,
  'auto size 1': True,
  'auto size 1 lock nth eye position': 5,
  'auto size 1 view angle delta in degrees': 9.0
}

global gWiff9x16FromAudienceAutoSizeD
#units are in feet
gWiff9x16FromAudienceAutoSizeD = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 22.0,
  'between eyes projected screen position': [0.0, (3.75 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  #'virtual self size multiplier': 1.0,
  #'auto size 1 distance ratio': 0.9,
  'virtual self size multiplier': 4.0,
  'auto size 1 distance ratio': 0.45,
  'auto size 1': True,
  #'auto size 1 distance ratio': 1.0,
  'auto size 1 lock nth eye position': 5,
  'auto size 1 view angle delta in degrees': 9.0
}

global gWiff9x16StandingCloseAutoSizeE
gWiff9x16StandingCloseAutoSizeE = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 5.0,
  'between eyes projected screen position': [0.0, (5.25 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  #'virtual self size multiplier': 1.0,
  'virtual self size multiplier': 0.4,
  #'auto size 1 distance ratio': 1.0,
  'auto size 1 distance ratio': 0.4,
  'auto size 1': True,
  #'auto size 1 distance ratio': 1.0,
  #'auto size 1 lock nth eye position': 5,
  #'auto size 1 view angle delta in degrees': 0.0
  #'auto size 1 view angle delta in degrees': 10.0
  'auto size 1 view angle delta in degrees': 5.0
}

global gWiff9x16StandingCloseAutoSizeF
gWiff9x16StandingCloseAutoSizeF = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 5.0,
  'between eyes projected screen position': [0.0, (5.25 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  #'virtual self size multiplier': 1.0,
  'virtual self size multiplier': 0.8,
  #'auto size 1 distance ratio': 1.0,
  'auto size 1 distance ratio': 0.25,
  'auto size 1': True,
  #'auto size 1 distance ratio': 1.0,
  #'auto size 1 lock nth eye position': 5,
  #'auto size 1 view angle delta in degrees': 0.0
  #'auto size 1 view angle delta in degrees': 10.0
  'auto size 1 view angle delta in degrees': 5.0
}

global gWiff9x16feetCloseJson
global gWiff9x16StandingCloseAutoSizeG
gWiff9x16StandingCloseAutoSizeG = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 5.0,
  'between eyes projected screen position': [0.0, (5.25 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  'virtual self size multiplier': 0.6,
  'auto size 1 distance ratio': 0.4,
  'auto size 1': True,
  #'auto size 1 distance ratio': 1.0,
  #'auto size 1 lock nth eye position': 5,
  'auto size 1 view angle delta in degrees': 0.0
}

global gWiff9x16feetCloseJson
#units are in feet
gWiff9x16feetCloseJson = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 5.0,
  #bottom of screen 3 feet above floor
  #eyes standing in front row 5.25 feet above floor
  #bottom of screen at -4.375 feet--below center
  'between eyes projected screen position': [0.0, (5.25 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  'virtual self size multiplier': 1.0
}

global gWiff9x16feetCloseBigJson
#units are in feet
gWiff9x16feetCloseBigJson = {
  'screen size' : [15.33333, 8.75],
  'screen ipd' : 0.2034121,
  'eye to screen distance': 5.0,
  'between eyes projected screen position': [0.0, (5.25 - 3.0 - 4.375)],
  'ipd in model units': 0.2034121,
  'virtual self size multiplier': 0.1
}

global gModelCoordinatesIpd
global gEnableTemporaryOffAxisProjectionTest

gEnableTemporaryOffAxisProjectionTest = True

#gOaptCameraNameList = [
#  'OaptDesktop1',
#  'OaptWiffBigInFrontOfScreen',
#  'OaptWiffBigAtScreen',
#  'OaptWiffBigBehindScreen',
#  'OaptWiffBigWayBehindScreen',
#  'OaptWiffSmallAtScreen',
#  'OaptWiffSmallInFrontOfScreen',
#  'OaptWiffSmallWayInFrontOfScreen',
#  'OaptWiffStandingClose',
#  'OaptWiffStandingCloseBig',
#  'OaptWiffStandingCloseBig2',
#  'OaptWiffAutoSize1AtScreen',
#  'OaptWiffAutoSize2AtScreen',
#  'OaptWiffAutoSize3AtScreen',
#  'OaptWiffAutoSize4AtScreen',
#  'OaptWiffAutoSize5AtScreen',
#  'OaptWiffAutoSize6AtScreen']

global gOaptPhysicalSettingsForCamera
gOaptPhysicalSettingsForCamera = {
  'OaptDesktop1': g1008x1792in1920x1200desktop30inch,
  'OaptWiffBigInFrontOfScreen': gWiff9x16FromAudienceBig,
  'OaptWiffBigAtScreen': gWiff9x16FromAudienceBig,
  'OaptWiffBigBehindScreen': gWiff9x16FromAudienceBig,
  'OaptWiffBigWayBehindScreen': gWiff9x16FromAudienceBig,
  'OaptWiffSmallAtScreen': gWiff9x16FromAudienceSmall,
  'OaptWiffSmallInFrontOfScreen': gWiff9x16FromAudienceSmall,
  'OaptWiffSmallWayInFrontOfScreen': gWiff9x16FromAudienceSmall,
  'OaptWiffStandingClose': gWiff9x16StandingCloseAutoSizeE,
  'OaptWiffStandingCloseBig': gWiff9x16StandingCloseAutoSizeF,
  'OaptWiffStandingCloseBig2': gWiff9x16StandingCloseAutoSizeG,
  'OaptWiffAutoSize1AtScreen': gWiff9x16FromAudienceAutoSizeA,
  #'OaptWiffAutoSize1AtScreen': gWiff9x16FromAudienceBig,
  'OaptWiffAutoSize2AtScreen': gWiff9x16FromAudienceAutoSizeB,
  'OaptWiffAutoSize3AtScreen': gWiff9x16FromAudienceAutoSizeA2,
  'OaptWiffAutoSize4AtScreen': gWiff9x16StandingCloseAutoSizeE,
  'OaptWiffAutoSize5AtScreen': gWiff9x16StandingCloseAutoSizeF,
  'OaptWiffAutoSize6AtScreen': gWiff9x16FromAudienceAutoSizeB,
}

#list of cameras from gOaptPhysicalSettingsForCamera
global gOaptCameraNameList
gOaptCameraNameList = []
for key, value in gOaptPhysicalSettingsForCamera.iteritems():
  gOaptCameraNameList.append(key)


notchScaleFactor = 250.0
#swap out operations from input deck for testing operations which aren't
#available through the input deck syntax yet

#rotation about z axis
global gSubstituteTestOperation1Json
gSubstituteTestOperation1Json = {

  'type': 'transform',
  'rotate': [0.0, 0.0, -45.0]

  #'type': 'calculator',
  #'function': 'damage*damage',

  #'scale': [notchScaleFactor, notchScaleFactor, notchScaleFactor]
  #'scale': [100.0, 100.0, 100.0]
  #'scale': [10.0, 10.0, 10.0]
  #'scale': [1.0, 1.0, 1.0]
  #'scale': [5.0, 5.0, 5.0]
}

#z displacement isolation
#rotation about x axis
global gSubstituteTestOperation2Json
gSubstituteTestOperation2Json = {
  'type': 'transform',
  'rotate': [-90.0, 0.0, 0.0]
}

#z displacement isolation
#global gSubstituteTestOperation2Json
#gSubstituteTestOperation2Json = {
#  'type': 'calculator',
#  'function': '0.0*iHat + 0.0*jHat + displ__Z*kHat',
#  'result array name': 'displ_',
#  #'result array name': 'mydisp_',
#}

global gSubstituteTestOperation3Json
gSubstituteTestOperation3Json = {
  'type': 'warpbyvector',
  'variable type': 'node',
  'variablename': 'displ_',
  #'vector name': 'mydisp_',
  #'scale factor': 20.0 * 10000.0,
  'scale': 10.0 * notchScaleFactor,
}

#x/y displacement isolation
global gSubstituteTestOperation4Json
gSubstituteTestOperation4Json = {
  'type': 'calculator',
  'function': 'displ__X*iHat + displ__Y*jHat + 0.0*kHat',
  'result array name': 'displ_',
  #'result array name': 'mydisp_',
}

global gSubstituteTestOperation5Json
gSubstituteTestOperation5Json = {
  'type': 'warpbyvector',
  'variable type': 'node',
  'variablename': 'displ_',
  #'vector name': 'mydisp_',
  #'scale factor': 20.0 * 10000.0,
  'scale': 10.0 * notchScaleFactor,
}

#all displacement exaggeration
global gSubstituteTestOperation6Json
gSubstituteTestOperation6Json = {
  'type': 'warpbyvector',
  'variable type': 'node',
  'vector name': 'displ_',
  #'vector name': 'mydisp_',
  #'scale factor': 20.0 * 10000.0,
  #'scale factor': 10.0 * 1000.0,
  'scale factor': 10.0 * notchScaleFactor,
}

#reflect operation
global gSubstituteTestOperation7Json
gSubstituteTestOperation7Json = {
  'type': 'reflect',
  'plane': 'Z Min',
  'center': 0.0,
}

#add pointset operation
global gSubstituteTestOperation8Json
gSubstituteTestOperation8Json = {
  'type': 'add unstructured grid',
  'filename': 'physicalpointsandoutline.vtu',
}

#add group operation
global gSubstituteTestOperation9Json
gSubstituteTestOperation9Json = {
  'type': 'group',
  'operation group list': ["ztestslice1", "ztestslice2"]
}


global gEnableSubstituteOperationTesting
gEnableSubstituteOperationTesting = True
global gSubstituteOperationTestingMap
gSubstituteOperationTestingMap = {
  "SubstituteTestOperation1": gSubstituteTestOperation1Json,
  "SubstituteTestOperation2": gSubstituteTestOperation2Json,
  "SubstituteTestOperation3": gSubstituteTestOperation3Json,
  "SubstituteTestOperation4": gSubstituteTestOperation4Json,
  "SubstituteTestOperation5": gSubstituteTestOperation5Json,
  "SubstituteTestOperation6": gSubstituteTestOperation6Json,
  "SubstituteTestOperation7": gSubstituteTestOperation7Json,
  "SubstituteTestOperation8": gSubstituteTestOperation8Json,
  "SubstituteTestOperation9": gSubstituteTestOperation9Json,
}

global gSubstituteOperationTestingList
gSubstituteOperationTestingList = []
for key, value in gSubstituteOperationTestingMap.iteritems():
  gSubstituteOperationTestingList.append(key)


#swap out imagesets from input deck for testing operations which aren't
#available through the input deck syntax yet

global gSubstituteTestImageset1Json
gSubstituteTestImageset1Json = {
  'color by vector component': 'mydisp_Z'
}
global gEnableSubstituteImagesetTesting
gEnableSubstituteImagesetTesting = True
global gSubstituteImagesetTestingMap
gSubstituteImagesetTestingMap = {
  "SubstituteTestImageset1": gSubstituteTestImageset1Json,
  #"SubstituteTestImageset2": gSubstituteTestImageset2Json,
  #"SubstituteTestImageset3": gSubstituteTestImageset3Json,
}

global gSubstituteImagesetTestingList
gSubstituteImagesetTestingList = []
for key, value in gSubstituteImagesetTestingMap.iteritems():
  gSubstituteImagesetTestingList.append(key)
#-----

global gDefaultTimeFormatString
gDefaultTimeFormatString = "Time: %.10e"
#gDefaultTimeFormatString = "Time: %.6e"
#gDefaultTimeFormatString = "Time: %f"

global gDefaultImageSizeX
gDefaultImageSizeX = 1920

global gDefaultImageSizeY
gDefaultImageSizeY = 1080

global gPixelBorderRatioXY
gPixelBorderRatioXY = [0.05,0.05]

class PhactoriAxisSettings:
  def __init__(self, inWhichAxis):
    #0-x axis, 1-yaxis, 2-zaxis
    self.WhichAxis = inWhichAxis
    self.TitleColor = [1.0, 1.0, 1.0]
    self.LabelColor = [1.0, 1.0, 1.0]
    self.TitleFontFamily = "Arial"
    self.TitleFontSize = 14
    self.LabelFontFamily = "Arial"
    self.LabelFontSize = 12

  def SetInfo1(self, inTitleColor, inLabelColor,
          inTitleFontFamily, inTitleFontSize,
          inLabelFontFamily, inLabelFontSize):
    self.TitleColor[0] = inTitleColor[0]
    self.TitleColor[1] = inTitleColor[1]
    self.TitleColor[2] = inTitleColor[2]
    self.LabelColor[0] = inLabelColor[0]
    self.LabelColor[1] = inLabelColor[1]
    self.LabelColor[2] = inLabelColor[2]
    if inTitleFontFamily != None:
      self.TitleFontFamily = inTitleFontFamily
    if inTitleFontSize != None:
      self.TitleFontSize = inTitleFontSize
    if inLabelFontFamily != None:
      self.LabelFontFamily = inLabelFontFamily
    if inLabelFontSize != None:
      self.LabelFontSize = inLabelFontSize

  def Apply(self, ioAxesGrid):
    if self.WhichAxis == 0:
      ioAxesGrid.XTitleColor[0] = self.TitleColor[0]
      ioAxesGrid.XTitleColor[1] = self.TitleColor[1]
      ioAxesGrid.XTitleColor[2] = self.TitleColor[2]
      ioAxesGrid.XLabelColor[0] = self.LabelColor[0]
      ioAxesGrid.XLabelColor[1] = self.LabelColor[1]
      ioAxesGrid.XLabelColor[2] = self.LabelColor[2]
      ioAxesGrid.XTitleFontFamily = self.TitleFontFamily
      ioAxesGrid.XTitleFontSize = self.TitleFontSize
      ioAxesGrid.XLabelFontFamily = self.LabelFontFamily
      ioAxesGrid.XLabelFontSize = self.LabelFontSize
    elif self.WhichAxis == 1:
      ioAxesGrid.YTitleColor[0] = self.TitleColor[0]
      ioAxesGrid.YTitleColor[1] = self.TitleColor[1]
      ioAxesGrid.YTitleColor[2] = self.TitleColor[2]
      ioAxesGrid.YLabelColor[0] = self.LabelColor[0]
      ioAxesGrid.YLabelColor[1] = self.LabelColor[1]
      ioAxesGrid.YLabelColor[2] = self.LabelColor[2]
      ioAxesGrid.YTitleFontFamily = self.TitleFontFamily
      ioAxesGrid.YTitleFontSize = self.TitleFontSize
      ioAxesGrid.YLabelFontFamily = self.LabelFontFamily
      ioAxesGrid.YLabelFontSize = self.LabelFontSize
    else:
      ioAxesGrid.ZTitleColor[0] = self.TitleColor[0]
      ioAxesGrid.ZTitleColor[1] = self.TitleColor[1]
      ioAxesGrid.ZTitleColor[2] = self.TitleColor[2]
      ioAxesGrid.ZLabelColor[0] = self.LabelColor[0]
      ioAxesGrid.ZLabelColor[1] = self.LabelColor[1]
      ioAxesGrid.ZLabelColor[2] = self.LabelColor[2]
      ioAxesGrid.ZTitleFontFamily = self.TitleFontFamily
      ioAxesGrid.ZTitleFontSize = self.TitleFontSize
      ioAxesGrid.ZLabelFontFamily = self.LabelFontFamily
      ioAxesGrid.ZLabelFontSize = self.LabelFontSize

class PhactoriColorSettings:
  def __init__(self, inUsePlotDefaultColors = False):
    self.mXAxisSettings = PhactoriAxisSettings(0)
    self.mYAxisSettings = PhactoriAxisSettings(1)
    self.mZAxisSettings = PhactoriAxisSettings(2)
    if inUsePlotDefaultColors:
      self.mBackground = [1.0, 1.0, 1.0]
      self.mEdgeColor = [0.0, 0.0, 0.5]

      #starting in paraview 5.0.1 this is the color
      #(mDiffuseColor) that affects the plot line
      #color
      #self.mDiffuseColor = [0.2, 0.2, 0.2]
      self.mDiffuseColor = [0.0, 0.0, 0.5]

      self.mBackfaceDiffuseColor = [0.2, 0.2, 0.2]
      self.mSelectionColor = [0.2, 0.2, 0.2]
      self.mAmbientColor = [1.0, 1.0, 1.0]
      self.mTextColor = [0.2, 0.2, 0.2]
      self.mTimeAnnotationColor = [0.2, 0.2, 0.2]
      self.mAxesGridColor = [0.2, 0.2, 0.2]
      #self.mGridAxesInfoHasBeenInitialized = True
      self.mXAxisSettings.SetInfo1(
              [0.2, 0.2, 0.2], [0.2, 0.2, 0.2], "Arial", 14, "Arial", 12)
      self.mYAxisSettings.SetInfo1(
              [0.2, 0.2, 0.2], [0.2, 0.2, 0.2], "Arial", 14, "Arial", 12)
      self.mZAxisSettings.SetInfo1(
              [0.2, 0.2, 0.2], [0.2, 0.2, 0.2], "Arial", 14, "Arial", 12)
      self.mAxesGridShowTicks = 1
      self.mAxesSettingsInititlized = True
    else:
      self.mBackground = [0.31999694819562063,
          0.3400015259021897, 0.4299992370489052]
      self.mEdgeColor = [0.0, 0.0, 0.5000076295109483]
      self.mDiffuseColor = [1.0, 1.0, 1.0]
      self.mBackfaceDiffuseColor = [1.0, 1.0, 1.0]
      self.mAmbientColor = [1.0, 1.0, 1.0]
      self.mSelectionColor = [1.0, 0.0, 1.0]
      self.mAxesGridColor = [1.0, 1.0, 1.0]
      self.mTextColor = [1.0, 1.0, 1.0]
      self.mTimeAnnotationColor = [1.0, 1.0, 1.0]
      #self.mGridAxesInfoHasBeenInitialized = False

    #weird test colors
    #self.mBackground = [1.0, 1.0, 1.0]
    #self.mEdgeColor = [1.0, 1.0, 0.0]
    #self.mDiffuseColor = [0.0, 1.0, 1.0]
    #self.mBackfaceDiffuseColor = [0.0, 1.0, 1.0]
    #self.mAmbientColor = [1.0, 1.0, 1.0]
    #self.mSelectionColor = [1.0, 0.0, 1.0]
    #self.mTextColor = [0.0, 1.0, 0.0]

  def SetParaviewRvRepColors(self, inParaviewRenderView,
          inParaviewRep):
    """given a paraview Representation instance, set the colors in the
       paraview representation instance from this color settings
       instance"""
    inParaviewRep.EdgeColor            = self.mEdgeColor
    inParaviewRep.DiffuseColor         = self.mDiffuseColor
    inParaviewRep.BackfaceDiffuseColor = self.mBackfaceDiffuseColor
    inParaviewRep.AmbientColor         = self.mAmbientColor
    inParaviewRep.SelectionColor       = self.mSelectionColor
    inParaviewRenderView.Background    = self.mBackground
    inParaviewRenderView.OrientationAxesLabelColor = \
        self.mTextColor
    if inParaviewRenderView.AxesGrid == None:
      PvSetupProperty(inParaviewRenderView.GetProperty("AxesGrid"))
      if inParaviewRenderView.AxesGrid == None:
        myDebugPrint3AndException("SetParaviewRvRepColors:\n"
        "could not get/create AxesGrid property\n")
    self.mXAxisSettings.Apply(inParaviewRenderView.AxesGrid)
    self.mYAxisSettings.Apply(inParaviewRenderView.AxesGrid)
    self.mZAxisSettings.Apply(inParaviewRenderView.AxesGrid)
    inParaviewRenderView.AxesGrid.GridColor = self.mAxesGridColor


  def ParsePlotColorSettingsFromJson(self, inJsn):
    """given a json structure, try to grab the color settings from it
       using plot keys rather than 3d image keys"""
    textColor = self.mXAxisSettings.TitleColor
    axesColor = self.mAxesGridColor
    dataColor = self.mEdgeColor
    surfaceColor = self.mDiffuseColor
    backfaceSurfaceColor = self.mBackfaceDiffuseColor
    if 'text color' in inJsn:
      textColor = inJsn['text color']
    if 'axes color' in inJsn:
      axesColor = inJsn['axes color']
    if 'point color' in inJsn:
      dataColor = inJsn['point color']
    if 'line color' in inJsn:
      dataColor = inJsn['line color']
    if 'surface color' in inJsn:
      surfaceColor = inJsn['surface color']
      backfaceSurfaceColor = surfaceColor
    if 'backface surface color' in inJsn:
      backfaceSurfaceColor = inJsn['backface surface color']
    self.mEdgeColor = dataColor
    self.mDiffuseColor = surfaceColor
    self.mBackfaceDiffuseColor = backfaceSurfaceColor
    self.mAxesGridColor = axesColor

    self.mOrientationAxesLabelColor = textColor
    self.mTextColor = textColor
    self.mXAxisSettings.SetInfo1(textColor, textColor, None, None, None, None)
    self.mYAxisSettings.SetInfo1(textColor, textColor, None, None, None, None)
    self.mZAxisSettings.SetInfo1(textColor, textColor, None, None, None, None)

    self.mBackground = getParameterFromBlock(inJsn,
                           'background color', self.mBackground)
    self.mTimeAnnotationColor = textColor

  def ParseColorSettingsFromJson(self, inJsn):
    """given a json structure, try to grab the color settings from it"""
    self.mDiffuseColor = getParameterFromBlock(inJsn,
                             'surface color', self.mDiffuseColor)
    #if 'back surface color' is not set, use surface color
    self.mBackfaceDiffuseColor = getParameterFromBlock(inJsn,
            'back surface color', self.mDiffuseColor)
    self.mAmbientColor = getParameterFromBlock(inJsn,
            'ambient color', self.mAmbientColor)
    self.mSelectionColor = getParameterFromBlock(inJsn,
            'selection color', self.mSelectionColor)
    self.mTextColor = getParameterFromBlock(inJsn,
            'text color', self.mTextColor)

    #for now, time annotation color will be text color
    self.mTimeAnnotationColor = self.mTextColor

    self.mBackground = getParameterFromBlock(inJsn,
            'background color', self.mBackground)
    self.mEdgeColor = getParameterFromBlock(inJsn,
            'edge color', self.mEdgeColor)

    self.mXAxisSettings.SetInfo1(self.mTextColor, self.mTextColor,
            None, None, None, None)
    self.mYAxisSettings.SetInfo1(self.mTextColor, self.mTextColor,
            None, None, None, None)
    self.mZAxisSettings.SetInfo1(self.mTextColor, self.mTextColor,
            None, None, None, None)

    #if 'axes color' is not set, use text color
    self.mAxesGrid = getParameterFromBlock(inJsn,
            'axes color', self.mTextColor)

def getParameterFromBlock(inBlock, inKey, inDefault):
  if inKey in inBlock:
    return inBlock[inKey]
  else:
    return inDefault


global gDefaultNumCounterDigits
gDefaultNumCounterDigits = 4
#gDefaultNumCounterDigits = 0

#global localCpViews
#global localCpWriters
#localCpViews = None
#localCpWriters = None

global localCoProcessorReference
localCoProcessorReference = None

global gMaxDebugPrintOutputLines
global gDebugPrintOutputLineCount
#gMaxDebugPrintOutputLines = 10000000
gMaxDebugPrintOutputLines = 100
gDebugPrintOutputLineCount = 0

#hacky debug printing to file mechanism
def myDebugPrint(outmsg):
  if PhactoriDbg():
    myDebugPrint3(outmsg)

def myDebugPrint2(outmsg):
  if PhactoriDbg():
    myDebugPrint3(outmsg)
  #global gMaxDebugPrintOutputLines
  #global gDebugPrintOutputLineCount
  #if gDebugPrintOutputLineCount > gMaxDebugPrintOutputLines:
  #  return
  #gDebugPrintOutputLineCount += 1
  ##ff = open('/space/jamauld/debug_out_dir/phactori_debug.txt', 'a+')
  ##fname = 'proc_' + str(SmartGetLocalProcessId()) + '_debug.txt'
  ##ff = open(fname, 'a+')
  ##ff.write(outmsg)
  ##ff.close()
  #print outmsg
  #if gDebugPrintOutputLineCount > gMaxDebugPrintOutputLines:
  #  print "gMaxDebugPrintOutputLines hit\n"
  #x = 1

global gMdp3RestrictToProcessListFlag
gMdp3RestrictToProcessListFlag = True
#gMdp3RestrictToProcessListFlag = False
global gMdp3ProcessIdList
gMdp3ProcessIdList = [0]
global gMdp3UseStandardOutFlag
gMdp3UseStandardOutFlag = True
#gMdp3UseStandardOutFlag = False
global gMdp3UseFileOutFlag
gMdp3UseFileOutFlag = False
#gMdp3UseFileOutFlag = True
global gMdp3PriorityRestriction
#gMdp3PriorityRestriction = 150
#gMdp3PriorityRestriction = 50
gMdp3PriorityRestriction = 500

#priority nominally
#100-very verbose info
#200-verbose info
#300-info
#400-debug
#500-warn
#600-error
#1000-error, and generate exception
def PhactoriDbg2(inPriority = 450, inOneProcessFlag = False, inOneProcessId = 0):
  return True

def PhactoriDbg(inPriority = 450, inOneProcessFlag = False, inOneProcessId = 0):
  """this is a test as to whether or not debug output should be written.
     The idea is that you do this test before printing debug output so that
     complicated strings do not need to be unnecessarily created and destroyed
     when not doing debug output"""
  global gMdp3PriorityRestriction
  if inPriority < gMdp3PriorityRestriction:
    return False

  global gMaxDebugPrintOutputLines
  global gDebugPrintOutputLineCount
  if gDebugPrintOutputLineCount > gMaxDebugPrintOutputLines:
    return False

  if inOneProcessFlag:
    if inOneProcessId != SmartGetLocalProcessId():
      return False
  #for parallel compatibility, we do everything on all processes _except_
  #do the output--so if we test for debugging then do something that requires
  #mpi communication, we go ahead and do it even if this process will have
  #no oputput
  #else:
  #  global gMdp3RestrictToProcessListFlag
  #  if gMdp3RestrictToProcessListFlag:
  #    global gMdp3ProcessIdList
  #    if SmartGetLocalProcessId() not in gMdp3ProcessIdList:
  #      return False

  return True

def PrintOnProcessZero(strToPrint):
  """prints out string if we are on process zero; typically for warnings"""
  if SmartGetLocalProcessId() == 0:
      print strToPrint

def myDebugPrint3(inMsg, inPriority = 450, inOneProcessFlag = False, inOneProcessId = 0):
  if PhactoriDbg(inPriority, inOneProcessFlag, inOneProcessId) == False:
    return

  global gMdp3RestrictToProcessListFlag
  if gMdp3RestrictToProcessListFlag:
    if SmartGetLocalProcessId() not in gMdp3ProcessIdList:
      #if we are only outputting some processes and this isn't one of those
      #we don't do any output--but we still deal with exception
      if inPriority >= 10000:
        raise Exception(inMsg)
      return

  global gMaxDebugPrintOutputLines
  global gDebugPrintOutputLineCount

  global gMdp3UseFileOutFlag
  if gMdp3UseFileOutFlag:
    fname = 'phdb_proc_' + str(SmartGetLocalProcessId()) + '.txt'
    if gDebugPrintOutputLineCount == 0:
      ff = open(fname, 'w')
    else:
      ff = open(fname, 'a+')
    ff.write(inMsg)
    if gDebugPrintOutputLineCount == gMaxDebugPrintOutputLines:
      ff.write("gMaxDebugPrintOutputLines hit\n")
    ff.close()

  global gMdp3UseStandardOutFlag
  if gMdp3UseStandardOutFlag:
    print inMsg
    if gDebugPrintOutputLineCount == gMaxDebugPrintOutputLines:
      print "gMaxDebugPrintOutputLines hit\n"

  gDebugPrintOutputLineCount += 1

  if inPriority >= 10000:
    raise Exception(inMsg)

def myDebugPrint3AndException(inMsg, inOneProcessFlag = False,
    inOneProcessId = 0):
  if PhactoriDbg(10000, inOneProcessFlag, inOneProcessId):
    myDebugPrint3(inMsg, 10000, inOneProcessFlag, inOneProcessId)

#def SetCpViewsAndCpWriters(inCpViews, inCpWriters):
#  global localCpViews
#  global localCpWriters
#  localCpViews = inCpViews
#  localCpWriters = inCpWriters

def GetSeparatorString():
  global gPipeAndViewsState
  return gPipeAndViewsState.mSeparatorString

def GetCurrentOutputResultsBlockCountId():
  global gPipeAndViewsState
  retStr = str(gPipeAndViewsState.mOutputResultsBlockCountId)
  if PhactoriDbg():
    myDebugPrint3('GetCurrentOutputResultsBlockCountId: ' + retStr + '\n')
  return retStr

def GetExtraFileString():
  global gPipeAndViewsState
  return gPipeAndViewsState.mRemeshRestartTag

def WriteImagesForCurrentPipeAndViewsState(datadescription):
  gPipeAndViewsState.WriteImages(datadescription)

def ExportOperationsDataForCurrentPipeAndViewsState(datadescription):
  gPipeAndViewsState.ExportOperationsData(datadescription)

def SetUpCoProcessor(inCoprocessor):
  global localCoProcessorReference
  localCoProcessorReference = inCoprocessor

class PhactoriRenderViewInfo:
  def __init__(self):
    self.RenderView1 = None
    self.DataRepresentation1 = None

#class PhactoriRenderViewInfo:
#  RenderView1 = None
#  DataRepresentation1 = None


global PhactoriRenderViewInfoList
PhactoriRenderViewInfoList = []

global currentPhactoriRenderViewInfo
currentPhactoriRenderViewInfo = None

phactoriFilterMap = {}

def AddFilterToFilterMap(inFilterName, inFilter):
  if PhactoriDbg(100):
    myDebugPrint3('AddFilterToFilterMap entered\n', 100)
  global phactoriFilterMap
  phactoriFilterMap[inFilterName] = inFilter
  count = 0
  for ii in phactoriFilterMap.iterkeys():
    if PhactoriDbg():
      myDebugPrint3(str(count) + ': ' + ii + '  : ' + str(phactoriFilterMap[ii]) + '\n')
    count += 1
  if PhactoriDbg(100):
    myDebugPrint3('AddFilterToFilterMap returning\n', 100)

def SetDataRepresentationToDefault(inDataRepresentation):

  #a3_vel__PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

  inDataRepresentation.SelectionPointLabelColor = [0.5, 0.5, 0.5]
  inDataRepresentation.SelectionPointFieldDataArrayName = 'displ'
  inDataRepresentation.SuppressLOD = 0
  inDataRepresentation.BlockVisibility = []
  inDataRepresentation.Position = [0.0, 0.0, 0.0]
  inDataRepresentation.BackfaceRepresentation = 'Follow Frontface'
  inDataRepresentation.SelectionOpacity = 1.0
  inDataRepresentation.SelectionPointLabelShadow = 0
  inDataRepresentation.OrientationMode = 'Direction'
  inDataRepresentation.ScaleMode = 'No Data Scaling Off'
  inDataRepresentation.Diffuse = 1.0
  inDataRepresentation.SelectionUseOutline = 0
  inDataRepresentation.SelectionPointLabelFormat = ''
  #shininess alteration, 1.0 is shiny
  inDataRepresentation.Specular = 0.1
  #inDataRepresentation.Specular = 1.0
  inDataRepresentation.SelectionVisibility = 1
  inDataRepresentation.InterpolateScalarsBeforeMapping = 1
  #inDataRepresentation.CustomRangeActive = [0, 0, 0]
  inDataRepresentation.Origin = [0.0, 0.0, 0.0]
  inDataRepresentation.Scale = [1.0, 1.0, 1.0]
  inDataRepresentation.SelectionCellLabelJustification = 'Left'
  inDataRepresentation.DiffuseColor = [1.0, 1.0, 1.0]
  inDataRepresentation.SelectionCellLabelOpacity = 1.0
  #inDataRepresentation.Source = []
  inDataRepresentation.Masking = 0
  inDataRepresentation.Opacity = 1.0
  inDataRepresentation.LineWidth = 1.0
  inDataRepresentation.MeshVisibility = 0
  inDataRepresentation.Visibility = 1
  inDataRepresentation.SelectionCellLabelFontSize = 18
  inDataRepresentation.SelectionPointLabelJustification = 'Left'
  #inDataRepresentation.OriginalBoundsRangeActive = [0, 0, 0]
  inDataRepresentation.SelectionPointLabelVisibility = 0
  inDataRepresentation.SelectOrientationVectors = ''
  inDataRepresentation.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
  inDataRepresentation.SelectionPointLabelFontFamily = 'Arial'
  inDataRepresentation.SelectScaleArray = ''

  if gParaViewCatalystVersionFlag <= 40100:
    inDataRepresentation.ColorAttributeType = gPointsString
  else:
    ColorBy(inDataRepresentation, (gPointsString,''))

  #inDataRepresentation.AxesOrigin = [0.0, 0.0, 0.0]
  inDataRepresentation.UserTransform = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  #shininess alteration, 40.0 is shiny
  inDataRepresentation.SpecularPower = 100.0
  #inDataRepresentation.SpecularPower = 40.0
  inDataRepresentation.Texture = []
  inDataRepresentation.SelectionCellLabelShadow = 0
  inDataRepresentation.AmbientColor = [1.0, 1.0, 1.0]
  inDataRepresentation.BlockOpacity = {}
  inDataRepresentation.MapScalars = 1
  inDataRepresentation.PointSize = 2.0
  inDataRepresentation.SelectionCellLabelFormat = ''
  inDataRepresentation.Scaling = 0
  inDataRepresentation.StaticMode = 0
  inDataRepresentation.SelectionCellLabelColor = [0.0, 1.0, 0.0]
  inDataRepresentation.EdgeColor = [0.0, 0.0, 0.5000076295109483]
  inDataRepresentation.SelectionCellLabelVisibility = 0
  inDataRepresentation.NonlinearSubdivisionLevel = 1
  inDataRepresentation.Representation = 'Surface'
  #inDataRepresentation.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
  #inDataRepresentation.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
  inDataRepresentation.Orientation = [0.0, 0.0, 0.0]
  #inDataRepresentation.ScalarOpacityUnitDistance = 0.8230042761351323
  inDataRepresentation.BackfaceOpacity = 1.0
  inDataRepresentation.SelectionPointLabelFontSize = 18
  inDataRepresentation.SelectionCellFieldDataArrayName = 'GlobalElementId'
  inDataRepresentation.SelectionColor = [1.0, 0.0, 1.0]
  inDataRepresentation.BlockColor = {}
  inDataRepresentation.Ambient = 0.0
  inDataRepresentation.ScaleFactor = 0.775
  inDataRepresentation.BackfaceAmbientColor = [1.0, 1.0, 1.0]
  #inDataRepresentation.ScalarOpacityFunction = a3_vel__PiecewiseFunction
  inDataRepresentation.SelectMaskArray = ''
  inDataRepresentation.SelectionLineWidth = 2.0
  inDataRepresentation.Interpolation = 'Gouraud'
  #inDataRepresentation.SelectMapper = 'Projected tetra'
  inDataRepresentation.SelectionCellLabelFontFamily = 'Arial'
  inDataRepresentation.SelectionCellLabelItalic = 0
  #inDataRepresentation.ExtractedBlockIndex = 2
  inDataRepresentation.SelectionPointLabelOpacity = 1.0
  #inDataRepresentation.UseAxesOrigin = 0
  inDataRepresentation.Pickable = 1
  #inDataRepresentation.CustomBoundsActive = [0, 0, 0]
  inDataRepresentation.SelectionRepresentation = 'Wireframe'
  inDataRepresentation.SelectionPointLabelBold = 0
  #inDataRepresentation.ColorArrayName = 'vel'
  inDataRepresentation.SelectionPointLabelItalic = 0
  #shininess alteration, 1 is shiny
  if gParaViewCatalystVersionFlag <= 40100:
    inDataRepresentation.AllowSpecularHighlightingWithScalarColoring = 0
    #inDataRepresentation.AllowSpecularHighlightingWithScalarColoring = 1
  inDataRepresentation.SpecularColor = [1.0, 1.0, 1.0]
  #inDataRepresentation.LookupTable = a3_vel__PVLookupTable
  inDataRepresentation.SelectionPointSize = 5.0
  inDataRepresentation.SelectionCellLabelBold = 0
  inDataRepresentation.Orient = 0

def PvSetupProperty(prop):
  domain = prop.FindDomain("vtkSMProxyListDomain")
  domain.SetDefaultValues(prop.SMProperty, False)
  prop.GetParent().UpdateVTKObjects()

def SetAxesGridVisibility(inPvRenderView, inSetting):
  if inPvRenderView.AxesGrid == None:
    PvSetupProperty(inPvRenderView.GetProperty("AxesGrid"))
  if inSetting:
    inPvRenderView.AxesGrid.Visibility = True
  else:
    inPvRenderView.AxesGrid.Visibility = False
  if PhactoriDbg(100):
    myDebugPrint3("SetAxesGridVisibility() AxesGrid.Visibility now " + \
            str(inPvRenderView.AxesGrid.Visibility) + "\n", 100)

def AddRenderView(inPhactoriImagesetInfo, inColorSettings,
    ImageBaseFileName, ImageType, ImageOverwriteFlag, PixelSizeX, PixelSizeY,
    numCounterDigits, inIsPlotFlag):

    global localCpViews
    global PhactoriRenderViewInfoList
    global currentPhactoriRenderViewInfo

    if PhactoriDbg(100):
      myDebugPrint3("AddRenderView entered\n", 100)

    newPhactoriRenderViewInfo = PhactoriRenderViewInfo()

    if ImageBaseFileName == None:
      image_base_name = "image_" + str(len(PhactoriRenderViewInfoList))
    else:
      image_base_name = ImageBaseFileName

    if ImageOverwriteFlag:
      image_base_name += "." + ImageType
      if PhactoriDbg():
        myDebugPrint3("image " + image_base_name + " will be overwritten repeatedly\n")
    else:
      image_base_name += "%t." + ImageType

    #image_base_name = "image_" + str(len(PhactoriRenderViewInfoList)) + "_%t.png"

    #RenderView1 = CreateCPView( CreateRenderView, image_base_name, 1, 0, 1, PixelSizeX, PixelSizeY, localCpViews )
    myDebugPrint3("about to create view and representation:\n"
      "current active source: " + str(GetActiveSource()) + "\n")

    global localCoProcessorReference

    global gSharedRenderView
    global gSharedLineChartView
    if inIsPlotFlag:
      if gSharedLineChartView == None:
        #gSharedLineChartView = localCoProcessorReference.CreateView("XYChartView",
        #        image_base_name, 1, 0, 1, PixelSizeX, PixelSizeY)
        gSharedLineChartView = CreateView("XYChartView")
        gSharedLineChartView.ViewSize = [PixelSizeX, PixelSizeY]
        #gSharedLineChartView.LeftAxisRangeMinimum = -1.0
        #gSharedLineChartView.LeftAxisRangeMaximum = 1.0
        #gSharedLineChartView.BottomAxisRangeMaximum = -2.0
        #gSharedLineChartView.RightAxisRangeMaximum = 2.0
        #gSharedLineChartView.TopAxisRangeMaximum = 2.0
      RenderView1 = gSharedLineChartView
    else:
      if gSharedRenderView == None:
      #if True:
        #import pdb
        #pdb.set_trace()

        gSharedRenderView = localCoProcessorReference.CreateView(CreateRenderView,
                image_base_name, 1, 0, 1, PixelSizeX, PixelSizeY)

      RenderView1 = gSharedRenderView

    inPhactoriImagesetInfo.mSharedPvRenderView2 = RenderView1

    RenderView1.add_attribute("cpNumCounterDigits", numCounterDigits)
    RenderView1.add_attribute("outputResultsBlockCountId", GetCurrentOutputResultsBlockCountId())
    RenderView1.add_attribute("associatedImagesetName",
        inPhactoriImagesetInfo.mName)

    if inIsPlotFlag == False:
      if gParaViewCatalystVersionFlag < 50502:
        RenderView1.LightSpecularColor = [1.0, 1.0, 1.0]
        RenderView1.LightIntensity = 1.0
        RenderView1.UseOffscreenRendering = 0
        RenderView1.UseOffscreenRenderingForScreenshots = 0
        RenderView1.LightDiffuseColor = [1.0, 1.0, 1.0]
        RenderView1.LightAmbientColor = [1.0, 1.0, 1.0]
        RenderView1.LightType = 'HeadLight'
        RenderView1.LightSwitch = 0
      RenderView1.UseOutlineForLODRendering = 0
      RenderView1.KeyLightAzimuth = 10.0
      RenderView1.UseTexturedBackground = 0
      RenderView1.UseLight = 1
      RenderView1.CameraPosition = [3.681775921856809, 3.2427490288581042, 6.445486324396935]
      RenderView1.FillLightKFRatio = 3.0
      RenderView1.Background2 = [0.0, 0.0, 0.165]
      RenderView1.FillLightAzimuth = -10.0
      RenderView1.LODResolution = 0.5
      RenderView1.BackgroundTexture = []
      RenderView1.InteractionMode = '3D'
      RenderView1.StencilCapable = 1
      RenderView1.CameraFocalPoint = [0.010957598686218262, 0.0, 3.6478042602539035e-05]
      RenderView1.ImageReductionFactor = 2
      RenderView1.CameraViewAngle = 30.0
      RenderView1.CameraParallelScale = 2.0952221328924265
      RenderView1.EyeAngle = 2.0
      RenderView1.HeadLightKHRatio = 3.0
      RenderView1.StereoRender = 0
      RenderView1.KeyLightIntensity = 0.75
      RenderView1.BackLightAzimuth = 110.0
      RenderView1.OrientationAxesInteractivity = 0
      if gParaViewCatalystVersionFlag <= 40100:
        RenderView1.UseInteractiveRenderingForSceenshots = 0
      RenderView1.Background = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]
      RenderView1.NonInteractiveRenderDelay = 0.0
      RenderView1.CenterOfRotation = [-1.9375, 0.0, 2.125]
      RenderView1.CameraParallelProjection = 0
      RenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
      RenderView1.HeadLightWarmth = 0.5
      RenderView1.MaximumNumberOfPeels = 4
      RenderView1.StereoType = 'Red-Blue'
      RenderView1.DepthPeeling = 1
      RenderView1.BackLightKBRatio = 3.5
      if gParaViewCatalystVersionFlag <= 40100:
        RenderView1.StereoCapableWindow = 1
      RenderView1.CameraViewUp = [-0.03140411320843562, 0.899940306729304, -0.43488070317911115]
      RenderView1.RemoteRenderThreshold = 20.0
      RenderView1.CacheKey = 1.0
      RenderView1.UseCache = 0
      RenderView1.KeyLightElevation = 50.0
      RenderView1.CenterAxesVisibility = 0
      RenderView1.MaintainLuminance = 0
      RenderView1.StillRenderImageReductionFactor = 1
      RenderView1.BackLightWarmth = 0.5
      RenderView1.FillLightElevation = -75.0
      RenderView1.MultiSamples = 0
      RenderView1.FillLightWarmth = 0.4
      RenderView1.AlphaBitPlanes = 1
      RenderView1.OrientationAxesVisibility = 1
      #RenderView1.CameraClippingRange = [8.69340289149996, 31.15868708101818]
      RenderView1.BackLightElevation = 0.0
      #RenderView1.ViewTime = 0.00403844503897222
      #RenderView1.ViewTime = 0.0
      RenderView1.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
      RenderView1.LODThreshold = 5.0
      RenderView1.CollectGeometryThreshold = 100.0
      RenderView1.UseGradientBackground = 0
      RenderView1.KeyLightWarmth = 0.6
      RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
      #RenderView1.OrientationAxesLabelColor = inColorSettings.mOrientationAxesLabelColor
      #rigid_body_impact_6_ff_e = CreateProducer( datadescription, "input" )

      SetAxesGridVisibility(RenderView1, 0)

    SetActiveView(RenderView1)

    if inPhactoriImagesetInfo.mPvDataRepresentation2 == None:
      if inIsPlotFlag:
        if inPhactoriImagesetInfo.m_PlotType == "PhactoriPlotOverTimeBlock":
          DataRepresentation1 = Show(UseIndexForXAxis = 0, XArrayName = "Time",
              AttributeType = "Row Data")
          #DataRepresentation1.AttributeType = "Row Data"
          if PhactoriDbg():
            myDebugPrint3("representation made for plot over time:\n" + \
                str(DataRepresentation1) + "\n")
          #DataRepresentation1.UseIndexForXAxis = 0
          #DataRepresentation1.XArrayName = "Time"
        elif inPhactoriImagesetInfo.m_PlotType == "PhactoriScatterPlotBlock":
          yvarStr = GetThresholdContourHackVariableNameString(
              inPhactoriImagesetInfo.m_YAxisVariableInfo)
          xvarStr = GetThresholdContourHackVariableNameString(
              inPhactoriImagesetInfo.m_XAxisVariableInfo)
          visvars = [yvarStr]
          atrtype = "Point Data"
          if inPhactoriImagesetInfo.m_YAxisVariableInfo.mVariableType == 'element':
            atrtype = "Cell Data"
          DataRepresentation1 = Show(inPhactoriImagesetInfo.m_producer,
              gSharedLineChartView, SeriesVisibility = visvars,
              UseIndexForXAxis = 0, XArrayName = xvarStr,
              AttributeType = atrtype)
          #DataRepresentation1.AttributeType = "Cell Data"
          #DataRepresentation1.XArrayName = 'carea'
          #DataRepresentation1.SeriesVisibility = ['vel__Magnitude (block_2)']
          #DataRepresentation1.SeriesVisibility = ['vel__Magnitude (Element Blocks)']
          if PhactoriDbg():
            myDebugPrint3("representation made for scatter plot:\n" + str(DataRepresentation1) + "\n")
            myDebugPrint3("m_producer for scatter plot:\n" + str(inPhactoriImagesetInfo.m_producer) + "\n")
            myDebugPrint3("rep input num points " + str(DataRepresentation1.Input.GetDataInformation().DataInformation.GetNumberOfPoints()) + "\n")
            myDebugPrint3("DataRepresentation1.ListProperties() " + str(DataRepresentation1.ListProperties()) + "\n")
            myDebugPrint3("DataRepresentation1.Input " + str(DataRepresentation1.Input) + "\n")
            myDebugPrint3("DataRepresentation1.AttributeType " + str(DataRepresentation1.AttributeType) + "\n")
            myDebugPrint3("DataRepresentation1.UseIndexForXAxis " + str(DataRepresentation1.UseIndexForXAxis) + "\n")
            myDebugPrint3("DataRepresentation1.Visibility " + str(DataRepresentation1.Visibility) + "\n")
            myDebugPrint3("DataRepresentation1.XArrayName " + str(DataRepresentation1.XArrayName) + "\n")
            myDebugPrint3("DataRepresentation1.SeriesVisibility " + str(DataRepresentation1.SeriesVisibility) + "\n")
            myDebugPrint3("DataRepresentation1.SeriesLabel " + str(DataRepresentation1.SeriesLabel) + "\n")
            myDebugPrint3("DataRepresentation1.SeriesColor " + str(DataRepresentation1.SeriesColor) + "\n")
            myDebugPrint3("DataRepresentation1.SeriesPlotCorner " + str(DataRepresentation1.SeriesPlotCorner) + "\n")
            myDebugPrint3("DataRepresentation1.SeriesLabelPrefix " + str(DataRepresentation1.SeriesLabelPrefix) + "\n")
            myDebugPrint3("DataRepresentation1.SeriesLineStyle " + str(DataRepresentation1.SeriesLineStyle) + "\n")
            myDebugPrint3("DataRepresentation1.SeriesLineThickness " + str(DataRepresentation1.SeriesLineThickness) + "\n")
            myDebugPrint3("DataRepresentation1.SeriesMarkerStyle " + str(DataRepresentation1.SeriesMarkerStyle) + "\n")
          for ii in range(1, len(DataRepresentation1.SeriesLineStyle), 2):
            DataRepresentation1.SeriesLineStyle[ii] = '0'
          for ii in range(1, len(DataRepresentation1.SeriesMarkerStyle), 2):
            DataRepresentation1.SeriesMarkerStyle[ii] = '4'
        else:
          myDebugPrint3AndException("bad m_PlotType")
      else:
        #UpdatePipelineWithCurrentTimeArgument(GetActiveSource())
        SetParaViewViewToCurrentTime(RenderView1)
        #myDebugPrint3(str(GetActiveView().ListProperties()) + "\n")
        #myDebugPrint3("ViewTime: " + str(GetActiveView().ViewTime) + "\n")
        DataRepresentation1 = Show()
        if PhactoriDbg():
          myDebugPrint3("representation made for 3d view:\n" + str(DataRepresentation1) + "\n")
      inPhactoriImagesetInfo.mPvDataRepresentation2 = DataRepresentation1
      inPhactoriImagesetInfo.mVisiblePvDataReps[0] = DataRepresentation1
    else:
      DataRepresentation1 = inPhactoriImagesetInfo.mPvDataRepresentation2

    if PhactoriDbg():
      myDebugPrint3("new DataRepresentation1: " + str(DataRepresentation1) + "\n")
    if PhactoriDbg():
      myDebugPrint3("new DataRepresentation1 input: " + str(DataRepresentation1.Input) + "\n")

    if inIsPlotFlag == False:
      SetDataRepresentationToDefault(DataRepresentation1)
      inColorSettings.SetParaviewRvRepColors(RenderView1, DataRepresentation1)

    newPhactoriRenderViewInfo.RenderView1 = RenderView1
    newPhactoriRenderViewInfo.DataRepresentation1 = DataRepresentation1

    currentPhactoriRenderViewInfo = newPhactoriRenderViewInfo
    PhactoriRenderViewInfoList.append(newPhactoriRenderViewInfo)

    if len(inPhactoriImagesetInfo.mVisibleOps) > 1:
      savedActiveSourceA = GetActiveSource()
      for ii in range(1,len(inPhactoriImagesetInfo.mVisibleOps)):
        if inPhactoriImagesetInfo.mVisiblePvDataReps[ii] == None:
          SetActiveSource(inPhactoriImagesetInfo.mVisibleOps[ii].\
              mParaViewFilter)
          newPvDataRep = Show()
          SetDataRepresentationToDefault(newPvDataRep)
          inPhactoriImagesetInfo.mVisiblePvDataReps[ii] = newPvDataRep
          if PhactoriDbg(100):
            myDebugPrint3("added extra newPvDataRep " + str(ii) + "\n"
            "newPvDataRep:\n" + str(newPvDataRep) + "\n")

      SetActiveSource(savedActiveSourceA)


    if PhactoriDbg(100):
      myDebugPrint3("AddRenderView returning\n", 100)

def SetCameraViewExplicitly(theRenderView, EyePosition,
                            LookAtPoint = None,
                            ClippingRange = None,
                            ViewUpVector = [0.0, 1.0, 0.0],
                            inCameraViewAngle = 30.0,
                            ViewBounds = None,
                            inPhactoriCamera = None,
                            inUseParallelProjection = False,
                            inParallelScale = 1.0):
  "Sets the camera view explicitly by specifying the camera position " \
  "(EyePosition), focal point (LookAtPoint), clipping range min/man " \
  "(ClippingRange), up view vector (ViewUpVector) " \
  "example:  SetCameraViewExplicitly(oneRenderView, [100.0, 90.0, 70.0], " \
  "  [0.0, 0.0, 0.0], [5.0, 200.0], [0.0, 1.0, 0.0]) " \
  "If ViewBounds is None and LookAtPoint is None, the global Data Bounds " \
  "Will be obtained by a GetGlobalDataBoundsParallel call to put the " \
  "LookAtPoint in the center of the data bounds.  If LookAtPoint is not " \
  "None, ViewBounds is not used; If LookAtPoint is None and ViewBounds " \
  "is not None, ViewBounds will be used and GetGlobalDataBoundsParallel " \
  "will not be called (to save mpi calls)"
  if PhactoriDbg(100):
    myDebugPrint3('phactori.SetCameraViewExplicitly entered\n', 100)

  if(LookAtPoint == None):
    if ViewBounds == None:
      ViewBounds = GetGlobalDataBoundsParallel()
    LookAtPoint = [ (ViewBounds[1] + ViewBounds[0]) * 0.5, \
                    (ViewBounds[3] + ViewBounds[2]) * 0.5, \
                    (ViewBounds[5] + ViewBounds[4]) * 0.5 ]

  if inPhactoriCamera != None:
    inPhactoriCamera.mOffAxisProjectionInfo.SetUpIfEnabled(theRenderView,
        EyePosition, LookAtPoint, ViewUpVector)

  if(ClippingRange == None):
    #if PhactoriDbg(100):
    #  myDebugPrint3('ClippingRange was None, must calculate\n', 100)
    xdelta = LookAtPoint[0] - EyePosition[0]
    ydelta = LookAtPoint[1] - EyePosition[1]
    zdelta = LookAtPoint[2] - EyePosition[2]
    delta  = math.sqrt(xdelta*xdelta + ydelta*ydelta + zdelta*zdelta)
    #ClippingRange = [delta * 0.05, delta * 2.0]
    ClippingRange = [delta * 0.01, delta * 10.0]

  #ClippingRange[0] *= 0.01
  #if PhactoriDbg(100):
  #  myDebugPrint3("CameraClippingRange (A): " + str(ClippingRange) + "\n")
  #  myDebugPrint3("CameraClippingRange (B): " + str(theRenderView.CameraClippingRange) + "\n")
  #ClippingRange[0] *= 1000.0
  #ClippingRange[1] *= 0.1

  theRenderView.CameraPosition = EyePosition
  theRenderView.CameraFocalPoint = LookAtPoint

  if inUseParallelProjection == True:
    theRenderView.CameraParallelProjection = 1
    theRenderView.CameraParallelScale = inParallelScale
    if PhactoriDbg(100):
      myDebugPrint3('CameraParallelProjection: 1\n')
      myDebugPrint3('CameraParallelScale: ' + str(inParallelScale) + '\n')
  else:
    theRenderView.CameraParallelProjection = 0
    #theRenderView.CameraClippingRange = ClippingRange

  theRenderView.CameraViewUp = ViewUpVector
  theRenderView.CenterOfRotation = LookAtPoint
  theRenderView.CameraViewAngle = inCameraViewAngle

  if PhactoriDbg(100):
    myDebugPrint3("CameraPosition: " + str(EyePosition) + "\n"
      "CameraFocalPoint: " + str(LookAtPoint) + "\n"
      "CameraClippingRange (C): " + str(ClippingRange) + "\n"
      "CameraViewUp: " + str(ViewUpVector) + "\n"
      "CenterOfRotation: " + str(LookAtPoint) + "\n"
      "CameraViewAngle: " + str(inCameraViewAngle) + "\n")

  if PhactoriDbg(100):
    myDebugPrint3('phactori.SetCameraViewExplicitly returning\n', 100)

#could do this different/better
def CheckForParallelVector(vec1, vec2):
  if PhactoriDbg():
    myDebugPrint3("CheckForParallelVector: " + str(vec1) + "  " + str(vec2) + "\n")
  maxcmpnt = 0
  maxval = abs(vec1[0])
  if(abs(vec1[1]) > maxval):
    maxcmpnt = 1
    maxval = abs(vec1[1])
  if(abs(vec1[2]) > maxval):
    maxcmpnt = 2
    maxval = abs(vec1[2])
  if PhactoriDbg():
    myDebugPrint3("maxcmpnt for vec1: " + str(maxcmpnt) + "\n")
  if vec2[maxcmpnt] == 0.0:
    if PhactoriDbg():
      myDebugPrint3("vec2 is zero there, not parallel\n")
    return False
  ratio = vec1[maxcmpnt] / vec2[maxcmpnt]
  for ii in [0,1,2]:
    tt1 = vec2[ii] * ratio
    if(abs(tt1 - vec1[ii]) > 0.000001):
      if PhactoriDbg():
        myDebugPrint3("vec2 component " + str(ii) + " is different, not parallel\n")
      return False
  if PhactoriDbg():
    myDebugPrint3("same ratio between vectors, parallel\n")
  return True

#given a bounds in the form [xmin, xmax, ymin, ymax, zmin, zmax], return the
# maximum dimension out of all those (max of xmax-xmin, ymax-ymin, and
# zmax-zmin)
def GetMaximumDimensionFromBounds(inBounds):
  maxDim = inBounds[1] - inBounds[0]
  testDim = inBounds[3] - inBounds[2]
  if(testDim > testDim):
    maxDim = testDim
  testDim = inBounds[5] - inBounds[4]
  if(testDim > testDim):
    maxDim = testDim
  return maxDim

def SetCameraLookAtPointAndLookDirection(inParaViewRenderView, \
        inLookDirection, \
        inImageSettings, \
        inFocalPoint = None, \
        inEyeToFocalPointDistance = None, \
        inEyePositionFactor = 1.0, \
        inViewUp=[0.0,1.0,0.0], \
        inCameraViewAngle = 30.0, \
        inViewBounds = None,
        inPhactoriCamera = None):
  """Set up camera with position, look at point, and positioning factor

  Sets the current camera view by looking at a focal point (inFocalPoint)
  along a specified vector (inLookDirection). The distance from the focal
  point is set using twice the largest of the data bounding box dimensions
  multiplied by inEyePositionFactor.  To see the whole of the data,
  inEyePositionFactor can be about 1.0, although trying different values
  will give you different filling of the screen and 'zooming' into the focal
  point, while values greater than 1.0 can be used to zoom out, possibly to
  bring the entire dataset into the image.  If inEyeToFocalPointDistance is
  something other than None, this distance will be used instead of the
  bounding box of the data.  If inFocalPoint is None, then
  the center of the data bounds will be used as a focal point. inViewBounds
  allows calling routine to pass along data bounds or other viewing bounds
  which are used to calculate the viewing distance from the focal point
  and the focal point if they are not supplied.  If inEyeToFocalPointDistance
  is not None and inFocalPoint is not None, than inViewBounds will not
  be used and the global data bounds will not be obtained.
  example:  SetCameraLookAtPointAndLookDirection([5.0, 4.0, 3.0],
  [1.0, 1.0, 1.0], 1.5, [0.0, 1.0, 0.0])
  """

  if PhactoriDbg(100):
    myDebugPrint3("SetCameraLookAtPointAndLookDirection entered\n", 100);

  #check for look view up same as look direction
  #need to do better check here
  if CheckForParallelVector(inLookDirection, inViewUp):
    inViewUp = [0.0, 0.0, 1.0]
    if CheckForParallelVector(inLookDirection, inViewUp):
      inViewUp = [0.0, 1.0, 0.0]

  #find maximum data dimension, then multiply by inEyePositionFactor
  if(inFocalPoint == None):
    if inViewBounds == None:
      inViewBounds = GetGlobalDataBoundsParallel()
    inFocalPoint = [ (inViewBounds[1] + inViewBounds[0]) * 0.5, \
                     (inViewBounds[3] + inViewBounds[2]) * 0.5, \
                     (inViewBounds[5] + inViewBounds[4]) * 0.5 ]

  if inEyeToFocalPointDistance != None:
    eyeToFocalPointDistance = inEyeToFocalPointDistance
  else:
    if inViewBounds == None:
      inViewBounds = GetGlobalDataBoundsParallel()
    #eyeToFocalPointDistance = GetMaximumDimensionFromBounds(inViewBounds)
    #eyeToFocalPointDistance *= inEyePositionFactor
    eyeToFocalPointDistance = CalcRelativeCameraDistance2(inFocalPoint,
            inLookDirection, inViewUp, inCameraViewAngle, inImageSettings,
            inViewBounds)
    eyeToFocalPointDistance *= inEyePositionFactor


  if PhactoriDbg():
    myDebugPrint3("  eyeToFocalPointDistance: " + str(eyeToFocalPointDistance) + "\n")
  if PhactoriDbg():
    myDebugPrint3("  inFocalPoint:            " + str(inFocalPoint) + "\n")
  if PhactoriDbg():
    myDebugPrint3("  inLookDirection:         " + str(inLookDirection) + "\n")
  if PhactoriDbg():
    myDebugPrint3("  inEyePositionFactor:     " + str(inEyePositionFactor) + "\n")
  if PhactoriDbg():
    myDebugPrint3("  inViewUp:                " + str(inViewUp) + "\n")

  #find eye position based on focal point, distance to eye point, and look direction
  xx = inLookDirection[0]
  yy = inLookDirection[1]
  zz = inLookDirection[2]
  lookDirectionMagnitude = math.sqrt(xx*xx + yy*yy + zz*zz)
  #lookDirectionMagnitude = 1.0
  if PhactoriDbg():
    myDebugPrint3("  lookDirectionMagnitude: " + str(lookDirectionMagnitude) + "\n")
  xx = xx / lookDirectionMagnitude
  yy = yy / lookDirectionMagnitude
  zz = zz / lookDirectionMagnitude
  eyePosition = [inFocalPoint[0] - xx * eyeToFocalPointDistance, \
                 inFocalPoint[1] - yy * eyeToFocalPointDistance, \
                 inFocalPoint[2] - zz * eyeToFocalPointDistance]
  if PhactoriDbg():
    myDebugPrint3("  eyePosition: " + str(eyePosition) + "\n")

  clippingRange = [0.01 * lookDirectionMagnitude, 100.0 * lookDirectionMagnitude]

  localUseParallelProjection = False
  localParallelScale = 1.0
  if inPhactoriCamera != None:
    if inPhactoriCamera.mUseParallelProjection:
      localUseParallelProjection = True
      if inPhactoriCamera.mParallelScaleAbsoluteOrRelative == 0:
        localParallelScale = inPhactoriCamera.mParallelScale
        if PhactoriDbg():
          myDebugPrint3("absolute parallel scale: " + \
                  str(localParallelScale) + "\n")
      else:
        if inViewBounds == None:
          inViewBounds = GetGlobalDataBoundsParallel()
        #for relative case, we are using the diagonal of the
        #bounding box to determine the parallel scale; other algorithms might
        #be to use the longest bounding box side, or determine all bounding
        #box corner positions from camera position and use the greatest
        #vertical or horizontal extent of those points
        dimx = inViewBounds[1] - inViewBounds[0]
        dimy = inViewBounds[3] - inViewBounds[2]
        dimz = inViewBounds[5] - inViewBounds[4]
        maxbbdim = math.sqrt(dimx*dimx + dimy*dimy + dimz*dimz)
        localParallelScale = maxbbdim * 0.5 * 1.05
        localParallelScale *= inPhactoriCamera.mParallelScale
        if PhactoriDbg():
          myDebugPrint3("relative parallel scale: " + \
                  str(localParallelScale) + "\n")

  SetCameraViewExplicitly(inParaViewRenderView, eyePosition, inFocalPoint,
      clippingRange, inViewUp, inCameraViewAngle, inViewBounds,
      inPhactoriCamera, localUseParallelProjection, localParallelScale)

  if PhactoriDbg(100):
    myDebugPrint3("SetCameraLookAtPointAndLookDirection returning\n", 100)


global gLocalProcessId
gLocalProcessId = -1

def SmartGetLocalProcessId():
  global gLocalProcessId
  #myDebugPrint3("SmartGetLocalProcessId entered\n", 100)
  if gLocalProcessId != -1:
    #print "gLocalProcessId already known: " + str(gLocalProcessId)
    return gLocalProcessId
  import vtkParallelCorePython
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()
  gLocalProcessId = globalController.GetLocalProcessId()
  #print "gLocalProcessId found: " + str(gLocalProcessId)
  #myDebugPrint3("exiting SmartGetLocalProcessId\n")
  return gLocalProcessId


def UseMPIToFillInSharedListNPerProcess(
        inThisProcVtkArray, inCountPerProc, outVtkArray, vtkArrayTypeIndex):
  """given a list of N values on each process fill in on each process a list
     which contains all the entries from each process
     (i.e. it is N * number_of_processes long)
     will work on either vtkDoubleArray types or vtkIntArray types, with
     vtkArrayTypeIndex = 0 for vtkIntArray types and 1 for vtkDoubleArray
     types"""

  mypid = SmartGetLocalProcessId()

  import vtkParallelCorePython
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()
  #gLocalProcessId = globalController.GetLocalProcessId()
  numproc = globalController.GetNumberOfProcesses()
  if PhactoriDbg(100):
    myDebugPrint3("mypid: " + str(mypid) + "  numproc: " + str(numproc) + \
    "  count: " + str(inCountPerProc) + "\n")

  outVtkArray.SetNumberOfValues(inCountPerProc * numproc)

  if vtkArrayTypeIndex == 0:
    fromRemoteProcList = vtk.vtkIntArray()
  else:
    fromRemoteProcList = vtk.vtkDoubleArray()

  for ii in range(0,numproc):
    prndx = ii * inCountPerProc
    if ii == mypid:
      globalController.Broadcast(inThisProcVtkArray, ii)
      for jj in range(0,inCountPerProc):
        outVtkArray.SetValue(prndx + jj, inThisProcVtkArray.GetValue(jj))
    else:
      fromRemoteProcList.SetNumberOfValues(inCountPerProc)
      globalController.Broadcast(fromRemoteProcList, ii)
      for jj in range(0,inCountPerProc):
        outVtkArray.SetValue(prndx + jj, fromRemoteProcList.GetValue(jj))

  return mypid, numproc


def UseMPIToCreateSharedPointList(inGlobalNodeIds, inPointXyzs):
  """used during finding nearest point (and distance) between one set of points
     and others.  The idea is that we have a fairly small set of points from one
     data source (e.g. one block of a multiblock or a decimated set of points)
     which we share to all processors (using mpi broadcast), and then each
     processor finds the nearest point and distance from all of it's points
     to the shared set, and then we use mpi reduce to find which is the
     closest and where it is.  This function takes the portion of the shared
     points which are local to this processor (obtained via
     PhactoriOperationBlock.MakeListOfAllPoints1()) and does mpi broadcasting
     to share with all the other processors and also to construct the whole
     list of shared points locally
     returns a vtkIntArray and vtkDoubleArray, where the int array is the list
     of all the global node ids and the double array is all the geometric
     pointx (stored x1y1z1x2y2z2...)
     inGlobalNodeIds is vtkIntArray
     inPointXyzs is vtkDoubleArray"""
  mypid = SmartGetLocalProcessId()

  import vtkParallelCorePython
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()
  #gLocalProcessId = globalController.GetLocalProcessId()
  numproc = globalController.GetNumberOfProcesses()
  if PhactoriDbg(100):
    myDebugPrint3("mypid: " + str(mypid) + "  numproc: " + str(numproc) + "\n")

  #first let everyone know how many points are in each process

  numptsInEachProc = vtk.vtkIntArray()
  numptsInEachProc.SetNumberOfValues(numproc)
  for ii in range(0,numproc):
    numptsInEachProc.SetValue(ii, 0)
  numptsInEachProc.SetValue(mypid, inGlobalNodeIds.GetNumberOfValues())
  #if PhactoriDbg(100):
  #  myDebugPrint3("numptsInEachProc before AllReduce:\n")
  #  for ii in range(0,numproc):
  #    myDebugPrint3(str(ii) + ": " + str(numptsInEachProc.GetValue(ii)) + "\n")

  glbNumptsInEachProc = vtk.vtkIntArray()
  glbNumptsInEachProc.SetNumberOfValues(numproc)

  globalController.AllReduce(numptsInEachProc, glbNumptsInEachProc, 0)
  #if PhactoriDbg(100):
  #  myDebugPrint3("numptsInEachProc after AllReduce:\n" + str(glbNumptsInEachProc) + "\n")
  #  for ii in range(0,numproc):
  #    myDebugPrint3(str(ii) + ": " + str(glbNumptsInEachProc.GetValue(ii)) + "\n")

  #now go through and create the shared points list
  #first the node ids
  glbNodeIdList = vtk.vtkIntArray()
  shareIdList = vtk.vtkIntArray()
  for ii in range(0,numproc):
    numPtsFromProcii = glbNumptsInEachProc.GetValue(ii)
    #if PhactoriDbg(100):
    #  myDebugPrint3("doing broadcast for process " + str(ii) + ":  numpts: " + str(numPtsFromProcii) + "\n")
    if ii == mypid:
      #if PhactoriDbg(100):
        #myDebugPrint3("I am this process, broadcasting to everyone\n")
      if(numPtsFromProcii > 0):
        globalController.Broadcast(inGlobalNodeIds, ii)
        for jj in range(0,numPtsFromProcii):
          glbNodeIdList.InsertNextValue(inGlobalNodeIds.GetValue(jj))
    else:
      #if PhactoriDbg(100):
      #  myDebugPrint3("I am not this process, receiving broadcast from this process: " + str(numPtsFromProcii) + "\n")
      if(numPtsFromProcii > 0):
        shareIdList.SetNumberOfValues(numPtsFromProcii)
        #if PhactoriDbg(100):
        #  myDebugPrint3("returned from broadcast call 2\n")
        globalController.Broadcast(shareIdList, ii)
        for jj in range(0,numPtsFromProcii):
          glbNodeIdList.InsertNextValue(shareIdList.GetValue(jj))
  #if PhactoriDbg(100):
  #  myDebugPrint3("glbNodeIdList after broadcasts:\n")
  #  for ii in range(0, glbNodeIdList.GetNumberOfValues()):
  #    myDebugPrint3(str(ii) + ": " + str(glbNodeIdList.GetValue(ii)) + "\n")

  #now the xyz geometry values for the nodes
  glbXyzList = vtk.vtkDoubleArray()
  shareXyzList = vtk.vtkDoubleArray()
  for ii in range(0,numproc):
    numPtsFromProcii = glbNumptsInEachProc.GetValue(ii)
    #if PhactoriDbg(100):
    #  myDebugPrint3("doing broadcast for process " + str(ii) + ":  numpts: " + str(numPtsFromProcii) + "\n")
    if ii == mypid:
      #if PhactoriDbg(100):
        #myDebugPrint3("I am this process, broadcasting to everyone\n")
      if(numPtsFromProcii > 0):
        globalController.Broadcast(inPointXyzs, ii)
        for jj in range(0,numPtsFromProcii*3):
          glbXyzList.InsertNextValue(inPointXyzs.GetValue(jj))
    else:
      #if PhactoriDbg(100):
      #  myDebugPrint3("I am not this process, receiving broadcast from this process: " + str(numPtsFromProcii) + "\n")
      if(numPtsFromProcii > 0):
        shareXyzList.SetNumberOfValues(numPtsFromProcii*3)
        #if PhactoriDbg(100):
        #  myDebugPrint3("returned from broadcast call 2\n")
        globalController.Broadcast(shareXyzList, ii)
        for jj in range(0,numPtsFromProcii*3):
          glbXyzList.InsertNextValue(shareXyzList.GetValue(jj))
  #if PhactoriDbg(100):
  #  myDebugPrint3("glbNodeIdList after broadcasts:\n")
  #  for ii in range(0, glbXyzList.GetNumberOfValues()):
  #    myDebugPrint3(str(ii) + ": " + str(glbXyzList.GetValue(ii)) + "\n")

  return glbNodeIdList, glbXyzList


def GetGlobalDataBoundsParallel(inFromWhichSource = None):
  if PhactoriDbg(100):
    myDebugPrint3("GetGlobalDataBoundsParallel entered\n", 100)
  import vtkParallelCorePython
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()

  if inFromWhichSource == None:
    localDataBounds = GetActiveSource().GetDataInformation().GetBounds()
  else:
    localDataBounds = inFromWhichSource.GetDataInformation().GetBounds()

  # negate so that MPI_MAX gets min instead of doing a MPI_MIN and MPI_MAX
  localarray = vtk.vtkDoubleArray()
  localarray.SetNumberOfTuples(6)
  localarray.SetValue(0, -localDataBounds[0])
  localarray.SetValue(1,  localDataBounds[1])
  localarray.SetValue(2, -localDataBounds[2])
  localarray.SetValue(3,  localDataBounds[3])
  localarray.SetValue(4, -localDataBounds[4])
  localarray.SetValue(5,  localDataBounds[5])
  globalarray = vtk.vtkDoubleArray()
  globalarray.SetNumberOfTuples(6)
  globalController.AllReduce(localarray, globalarray, 0)
  globalDataBounds = [ \
    -globalarray.GetTuple1(0), globalarray.GetTuple1(1), \
    -globalarray.GetTuple1(2), globalarray.GetTuple1(3), \
    -globalarray.GetTuple1(4), globalarray.GetTuple1(5) \
    ]
  if PhactoriDbg():
    myDebugPrint3("localDataBounds:  " + str(localDataBounds) + \
      "\nglobalDataBounds: " + str(globalDataBounds) + "\n")
  if PhactoriDbg(100):
    myDebugPrint3("GetGlobalDataBoundsParallel returning\n", 100)
  return globalDataBounds


global gDefaultLookatDistanceDatasizeRelative
gDefaultLookatDistanceDatasizeRelative = 1.0
global gDefaultLookatDistanceDatasizeRelativeMultiplier
gDefaultLookatDistanceDatasizeRelativeMultiplier = 4.0

def GetXyzForNodeOrElementParallelOneBlock(inInputCsData, inIdIsNode,
        inGlobalId, outXyz):
  """check for inGlobalId and set outXyz if present

  utility function called by GetXyzForNodeOrElementParallelRecurse1, this
  takes one unstructured grid as input, and sees if it has the node or element
  with id inGlobalId.  If it does, the method sets outXyz to the geometric
  location of the node (xyz) or center of the element bounding box (xyz) and
  returns true, otherwise it returns false without changing outXyz
  """

  if PhactoriDbg(100):
    myDebugPrint3('GetXyzForNodeOrElementParallelOneBlock entered\n', 100)

  globalIdArray = None
  if inIdIsNode:
    ptOrElData = inInputCsData.GetPointData()
    globalIdArray = ptOrElData.GetArray('GlobalNodeId')
  else:
    ptOrElData = inInputCsData.GetCellData()
    globalIdArray = ptOrElData.GetArray('GlobalElementId')

  if globalIdArray == None:
    if PhactoriDbg():
      myDebugPrint3("  this process/block has no Global Node or Element Id array to contain " + str(inGlobalId) + "\n")
    return False

  numTuples = globalIdArray.GetNumberOfTuples()
  thisProcessHasTheId = False
  idIndex = -1
  for ii in range(0, numTuples):
    #myDebugPrint3(" testing " + str(ii) + " against " + str(inGlobalNodeId) + "\n")
    #myDebugPrint3(" type array: " + str(type(globalNodeIdArray)) + "  type ii:" + str(type(ii)) + "\n")
    vv = globalIdArray.GetTuple1(ii)
    if vv == inGlobalId:
      thisProcessHasTheId = True
      idIndex = ii
      break

  if not thisProcessHasTheId:
    if PhactoriDbg():
      myDebugPrint3("  this process/block doesn't contain id" + \
        str(inGlobalId) + "\n")
    return False

  if PhactoriDbg():
    myDebugPrint3("  this process/block contains id " + str(inGlobalId) + "\n")

  if inIdIsNode:
    pointsArray = inInputCsData.GetPoints()
    numPoints = pointsArray.GetNumberOfPoints()
    if idIndex >= numPoints:
      if PhactoriDbg():
        myDebugPrint3("  this process/block has problem with index, setting xyz 0\n")
      outXyz[0] = 0.0
      outXyz[1] = 0.0
      outXyz[2] = 0.0
      return False

    thePoint = pointsArray.GetPoint(idIndex, outXyz)

    if PhactoriDbg():
      myDebugPrint3("  outXyz set to: " + str(outXyz) + "\n")
    return True
  else:
    myBounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    #myCells = inInputCsData.GetCells()
    #oneCell = myCells.GetCell(idIndex)
    oneCell = inInputCsData.GetCell(idIndex)
    oneCell.GetBounds(myBounds)
    #myCells.GetCellBounds(idIndex, myBounds)
    #ptOrElData.GetCellBounds(idIndex, myBounds)
    outXyz[0] = 0.5 * (myBounds[0] + myBounds[1])
    outXyz[1] = 0.5 * (myBounds[2] + myBounds[3])
    outXyz[2] = 0.5 * (myBounds[4] + myBounds[5])
    #xmin, xmax, ymin, ymax, zmin, zmax = myCells.GetCellBounds(idIndex)
    #outXyz[0] = 0.5 * (xmin + xmax)
    #outXyz[1] = 0.5 * (ymin + ymax)
    #outXyz[2] = 0.5 * (zmin + zmax)
    return True

def GetXyzForNodeOrElementParallelRecurse1(inInputCsData,
        inIdIsNode, inGlobalId, outXyz):
  """utility function used by GetXyzForNodeOrElementParallel

  looking for inGlobalId to potentially set outXyz, recurses through
  structure and calls GetXyzForNodeOrElementParallelOneBlock() on
  unstructured grids to do the real work of checking for the id and setting
  outXyz
  """

  #myDebugPrint3('GetXyzForNodeOrElementParallelRecurse1 entered\n', 100)

  icsdClassname = inInputCsData.GetClassName()
  if icsdClassname == "vtkMultiBlockDataSet" or \
     icsdClassname == "vtkExodusIIMultiBlockDataSet":
    #myDebugPrint3('recursing: ' + icsdClassname + '\n')
    numBlocks = inInputCsData.GetNumberOfBlocks()
    for ii in range(0, numBlocks):
      oneBlock = inInputCsData.GetBlock(ii)
      if(oneBlock != None):
        found = GetXyzForNodeOrElementParallelRecurse1(
                    oneBlock, inIdIsNode, inGlobalId, outXyz)
        if found:
          return True
  else:
    found = GetXyzForNodeOrElementParallelOneBlock(
                inInputCsData, inIdIsNode, inGlobalId, outXyz)
    if found:
      return True

  return False

def GetXyzForMinOrMaxVariable(inParaViewSource, inMinFlag,
        inVariableInfo, outXyz):
  """find the xyz location of the maximum or minimum value of the variable
     in the argument (at node or in element), recurses as necessary, does
     MPI, calls GetXyzForNodeOrElementParallel once id is found"""

  if PhactoriDbg(100):
    myDebugPrint3("GetXyzForMinOrMaxVariable entered\n"
        "  inMinFlag: " + str(inMinFlag) + "\n"
        "  variable name: " + str(inVariableInfo.mVariableName) + "\n", 100)

  #for now, disable output while doing parallel min/max search
  global gMdp3PriorityRestriction
  save_gMdp3PriorityRestriction = gMdp3PriorityRestriction
  gMdp3PriorityRestriction = 500

  #get min or max node or element id for variable
  csData = inParaViewSource.GetClientSideObject().GetOutputDataObject(0)
  FindNodeOrElementIdForMinMax(csData, inVariableInfo)
  if inMinFlag:
    theId = inVariableInfo.mStats.mMinId
  else:
    theId = inVariableInfo.mStats.mMaxId

  #decide if it is a node or an element
  detectResult = inVariableInfo.DetectVariableType(
        inParaViewSource, True)
  if detectResult == False:
    #myDebugPrint3AndException(
    #    "GetXyzForMinOrMaxVariable:\n"
    #    "could not detect element or node type\n"
    #    "variable name is: " + inVariableInfo.mVariableName + "\n")
    myDebugPrint3(
        "GetXyzForMinOrMaxVariable:\n"
        "could not detect element or node type (using 0 0 0)\n"
        "variable name is: " + inVariableInfo.mVariableName + "\n")
    outXyz[0] = 0.0
    outXyz[1] = 0.0
    outXyz[2] = 0.0
    localProcessCouldDetect = False
  if inVariableInfo.mVariableType == 'node':
    nodeFlag = True
    localProcessCouldDetect = True
  elif inVariableInfo.mVariableType == 'element':
    nodeFlag = False
    localProcessCouldDetect = True
  else:
    myDebugPrint3AndException(
        "GetXyzForMinOrMaxVariable:\n"
        "variable is not element or node type\n"
        "variable name is: " + inVariableInfo.mVariableName + "\n")

  #get xyz location of node or element
  GetXyzForNodeOrElementParallel(inParaViewSource, nodeFlag, theId, outXyz,
    localProcessCouldDetect)

  #restore output setting
  gMdp3PriorityRestriction = save_gMdp3PriorityRestriction

  if PhactoriDbg(100):
    myDebugPrint3(
      "outXyz: " + str(outXyz) + "\n"
      "GetXyzForMinOrMaxVariable returning\n", 100)


def GetXyzForNodeOrElementParallel(inParaViewSource,
        inIdIsNode, inGlobalId, outXyz, inLocalProcessCouldDetect = True):
  """given a node or element id, get item xyz location in parallel

  Recurses through all structures on all processes to find a node or element,
  then does allreduce to spread that information to all processes (for use in
  pointing cameras at a node or element with a given id; if this process has
  the element with the
  id inGlobalId, it will set outXyz and send the info around, if this
  process does not have inGlobalId, it will expect to receive the proper
  outXyz from the node that does.
  """

  if PhactoriDbg(100):
    myDebugPrint3("GetXyzForNodeOrElementParallel entered, IdIsNode: " + str(inIdIsNode) + "  id: " + str(inGlobalId) + "\n", 100)

  #theSource = GetCurrentSource()
  #theSource = GetActiveSource()
  if inLocalProcessCouldDetect:
    csData = inParaViewSource.GetClientSideObject().GetOutputDataObject(0)

    found = GetXyzForNodeOrElementParallelRecurse1(csData,
            inIdIsNode, inGlobalId, outXyz)
  else:
    found = False


  if found:
    if PhactoriDbg():
      myDebugPrint3("  this process has GlobalId " + str(inGlobalId) + ", doing broadcast of info\n")
    UseReduceToSpreadValues(outXyz)
    #GetXyzForNodeParallelBroadcast(thePointArrays, outXyz)
  else:
    if PhactoriDbg():
      myDebugPrint3("  this process does not have GlobalId " + str(inGlobalId) + ", doing receive of info\n")
    UseReduceToSpreadValues(outXyz)
    #GetXyzForNodeParallelReceive(outXyz)

  if PhactoriDbg(100):
    myDebugPrint3("GetXyzForNodeOrElementParallel returning\n", 100)


def vecDotProduct(inVecA, inVecB):
  """ returns dotproct inVecZ dot inVecB """
  return inVecA[0] * inVecB[0] + inVecA[1] * inVecB[1] + inVecA[2] * inVecB[2]

def vecCrossProduct(inVecA, inVecB):
  """ returns inVecZ X inVecB """
  return [inVecA[1] * inVecB[2] - inVecA[2] * inVecB[1],
          inVecA[2] * inVecB[0] - inVecA[0] * inVecB[2],
          inVecA[0] * inVecB[1] - inVecA[1] * inVecB[0]]

def vecCrossProduct2(outVec, inVecA, inVecB):
  """ calculates inVecA X inVecB and puts the result in outVec"""
  outVec[0] = inVecA[1] * inVecB[2] - inVecA[2] * inVecB[1]
  outVec[1] = inVecA[2] * inVecB[0] - inVecA[0] * inVecB[2]
  outVec[2] = inVecA[0] * inVecB[1] - inVecA[1] * inVecB[0]

def vecCopy(destinationVec, sourceVec):
  destinationVec[0] = sourceVec[0]
  destinationVec[1] = sourceVec[1]
  destinationVec[2] = sourceVec[2]

def vecMagnitude(inVec):
  xx = inVec[0]
  yy = inVec[1]
  zz = inVec[2]
  return math.sqrt(xx*xx + yy*yy + zz*zz)

def vecMagnitudeSquared(inVec):
  xx = inVec[0]
  yy = inVec[1]
  zz = inVec[2]
  return xx*xx + yy*yy + zz*zz

def vecNormalize(inVec):
  xx = inVec[0]
  yy = inVec[1]
  zz = inVec[2]
  mag = math.sqrt(xx*xx + yy*yy + zz*zz)
  return [xx/mag, yy/mag, zz/mag]

def vecNormalize2(outVec, inVec):
  xx = inVec[0]
  yy = inVec[1]
  zz = inVec[2]
  mag = math.sqrt(xx*xx + yy*yy + zz*zz)
  outVec[0] = xx/mag
  outVec[1] = yy/mag
  outVec[2] = zz/mag

def vecFromAToB(inVecA, inVecB):
  return [inVecB[0] - inVecA[0], inVecB[1] - inVecA[1], inVecB[2] - inVecA[2]]

def vecDistanceSquared(inPtA, inPtB):
  """calculate the square of the distance between two points"""
  ddx = inPtA[0] - inPtB[0]
  ddy = inPtA[1] - inPtB[1]
  ddz = inPtA[2] - inPtB[2]
  return ddx*ddx + ddy*ddy + ddz*ddz

def vecDistance(inPtA, inPtB):
  """calculate the distance between two points"""
  return math.sqrt(vecDistanceSquared(inPtA, inPtB))

def vecAdd(inVecA, inVecB):
  return [inVecA[0]+inVecB[0],inVecA[1]+inVecB[1],inVecA[2]+inVecB[2]]

def vecScale(inScale, inVec):
  return [inScale*inVec[0], inScale*inVec[1], inScale*inVec[2]]

def vecMultiplyAdd(inVecA, inVecB, inMM):
  """returns inVecA + (inMM * inVecB)"""
  return [inVecA[0] + inMM * inVecB[0],
          inVecA[1] + inMM * inVecB[1],
          inVecA[2] + inMM * inVecB[2]]

def CalcRelativeCameraDistance2_AA(inFocalPoint, inNormCameraDir,
        inTestPoint, inFov = 30.0):
  """finds point X which is along the inNormCameraDir vector starting from
     point inFocalPoint, such that inTestPoint is visible.  This is
     accomplished by taking a parametric value ss to move the camera position
     out from the focal point (along the inNormCameraDir vector) and finding
     the parametric value such that the angle between the vector from the
     camera position to the focal point and the camera position to inTestPoint
     is equal to (or greater than) the camera field of view.  The parametric
     value is returned.  We use the law of sines.  assumes inNormCameraDir
     is normalized and pointing in the direction from inFocalPoint towards
     where we want to place the camera (the opposite of the look direction)"""

  #angle between camera-to-focal point vector and camera-to-testpoint vector
  #is equal to half the camera fOV.  Call this angle AA.
  #The side opposite this angle has length equal to the distance between
  #the focal point and the test point.  Call this distance aa.  We also
  #know the angle between the focal point-to-camera vector and the focal
  #point-to-testpoint vector (because both vectors point in a fixed direction).
  #call this angle BB.  The remaining angle is CC, and CC = 180 - (AA + BB).
  #thus the length of the parametric side is determined by the law of sines,
  #cc/sin(CC) = aa/sin(AA) (cc is the length of interest),
  #so cc = aa * sin(CC) / sin(AA)
  #

  #need to set to 1/2 fov
  if PhactoriDbg():
    myDebugPrint3("CalcRelativeCameraDistance2_AA entered\n");
  if PhactoriDbg():
    myDebugPrint3("inTestPoint: " + str(inTestPoint) + "\n");
  angleAA = math.radians(inFov * 0.5)
  if PhactoriDbg():
    myDebugPrint3("angleAA: " + str(math.degrees(angleAA)) + "\n");

  focalToTestVec = [inTestPoint[0] - inFocalPoint[0],
                    inTestPoint[1] - inFocalPoint[1],
                    inTestPoint[2] - inFocalPoint[2]]
  aa = vecMagnitude(focalToTestVec)
  if PhactoriDbg():
    myDebugPrint3("aa: " + str(aa) + "\n");
  if aa == 0.0:
    return 0.0
  focalToTestVec[0] /= aa
  focalToTestVec[1] /= aa
  focalToTestVec[2] /= aa
  if PhactoriDbg():
    myDebugPrint3("inNormCameraDir: " + str(inNormCameraDir) + "\n");
  if PhactoriDbg():
    myDebugPrint3("focalToTestVec: " + str(focalToTestVec) + "\n");
  cosAngleBB = vecDotProduct(inNormCameraDir, focalToTestVec)
  if PhactoriDbg():
    myDebugPrint3("cosAngleBB: " + str(cosAngleBB) + "\n");
  if cosAngleBB >= 1.0:
    #camera-to-test-point in line with camera-to-focal-point
    #and already in camera view
    if PhactoriDbg():
      myDebugPrint3("cosAngleBB >= 1.0 cc is 0.0\n");
    cc = 0.0
    return cc
  elif cosAngleBB <= -1.0:
    #camera-to-test-point in line with camera-to-focal-point
    #with camera needing movement back to test point
    if PhactoriDbg():
      myDebugPrint3("cosAngleBB <= -1.0, cc is camera to test point dist\n");
    cc = aa
    return cc
  angleBB = math.acos(cosAngleBB)
  if PhactoriDbg():
    myDebugPrint3("angleBB: " + str(math.degrees(angleBB)) + "\n");
  angleCC = math.pi - (angleAA + angleBB)
  if PhactoriDbg():
    myDebugPrint3("angleCC: " + str(math.degrees(angleCC)) + "\n");
  cc = aa * math.sin(angleCC) / math.sin(angleAA)
  if PhactoriDbg():
    myDebugPrint3("cc: " + str(cc) + "\n");
  if PhactoriDbg():
    myDebugPrint3("CalcRelativeCameraDistance2_AA returning\n");
  return cc

def projectPointOntoPlane(inPoint, inPlanePoint, inPlaneNormal):
  vecA = [inPoint[0] - inPlanePoint[0],
          inPoint[1] - inPlanePoint[1],
          inPoint[2] - inPlanePoint[2]]
  dotp1 = vecDotProduct(vecA, inPlaneNormal)
  projectedPoint = [inPoint[0] - dotp1 * inPlaneNormal[0],
                    inPoint[1] - dotp1 * inPlaneNormal[1],
                    inPoint[2] - dotp1 * inPlaneNormal[2]]
  return projectedPoint

def CalcRelativeCameraDistance2_BB(inFocalPoint, inNormCameraDir,
        inNormUpVector, inNormSideVector, inFovV, inFovH, inTestPoint):
  """takes test point and makes two test points, one projected onto the
     plane defined by the up vector and one onto the plane defined by the
     up vector/look direction cross product; uses these two points to
     test vs. vertical FOV and horizontal FOV, and returns the biggest
     result.  Assumes inNormUpVector is normalized as is inNormCameraDir
     and also inNormSideVector"""
  if PhactoriDbg():
    myDebugPrint3("CalcRelativeCameraDistance2_BB entered\n");
  if PhactoriDbg():
    myDebugPrint3("inTestPoint: " + str(inTestPoint) + "\n");

  pointV = projectPointOntoPlane(inTestPoint, inFocalPoint, inNormSideVector)
  if PhactoriDbg():
    myDebugPrint3("pointV: " + str(pointV) + "\n");
  #find point projected onto plane defined by inFocalPoint and localUpVector,
  #which will be the test point to test horizontal FOV against
  pointH = projectPointOntoPlane(inTestPoint, inFocalPoint, inNormUpVector)
  if PhactoriDbg():
    myDebugPrint3("pointH: " + str(pointH) + "\n");

  dd1 = CalcRelativeCameraDistance2_AA(inFocalPoint, inNormCameraDir,
            pointH, inFovH)
  if PhactoriDbg():
    myDebugPrint3("dd1: " + str(dd1) + "\n");
  dd2 = CalcRelativeCameraDistance2_AA(inFocalPoint, inNormCameraDir,
            pointV, inFovV)
  if PhactoriDbg():
    myDebugPrint3("dd2: " + str(dd2) + "\n");
  if dd2 > dd1:
    dd1 = dd2
  if PhactoriDbg():
    myDebugPrint3("return dd: " + str(dd1) + "\n");
  return dd1


def CalcRelativeCameraDistance2(inFocalPoint, inLookDirection, inUpVector,
        inFov, inImageSettings, inBounds):
  #for each focal point, see how far back from the focal point the
  #camera position needs to be to see the focal point.  Use the largest.
  inXyPixelSize = inImageSettings.mImageSize
  inPixelBorderRatioXY = inImageSettings.mPixelBorderRatioXY
  myDebugPrint("CalcRelativeCameraDistance2 entered\n");
  myDebugPrint("inFocalPoint: " + str(inFocalPoint) + "\n");
  myDebugPrint("inLookDirection: " + str(inLookDirection) + "\n");
  myDebugPrint("inUpVector: " + str(inUpVector) + "\n");
  myDebugPrint("inFov: " + str(inFov) + "\n");
  myDebugPrint("inXyPixelSize: " + str(inXyPixelSize) + "\n");
  myDebugPrint("inPixelBorderRatioXY: " + str(inPixelBorderRatioXY) + "\n");
  normCameraDir = vecNormalize(inLookDirection)
  normCameraDir[0] = -normCameraDir[0]
  normCameraDir[1] = -normCameraDir[1]
  normCameraDir[2] = -normCameraDir[2]
  myDebugPrint("normCameraDir: " + str(normCameraDir) + "\n");

  #find orthoganol up and side vectors,
  normUpVector = vecNormalize(inUpVector)
  sideVector = vecCrossProduct(normUpVector, normCameraDir)
  normUpVector = vecCrossProduct(normCameraDir, sideVector)
  if PhactoriDbg():
    myDebugPrint3("sideVector: " + str(sideVector) + "\n");
  if PhactoriDbg():
    myDebugPrint3("normUpVector: " + str(normUpVector) + "\n");

  #find vertical fov, based on pixel border percentage
  vertPixX = 0.5 * float(inXyPixelSize[0])
  vertPixY = 0.5 * float(inXyPixelSize[1])
  vertPixBrdrX = inPixelBorderRatioXY[0] * float(inXyPixelSize[0])
  vertPixBrdrY = inPixelBorderRatioXY[1] * float(inXyPixelSize[1])

  vertPixBrdrX = math.ceil(vertPixBrdrX)
  vertPixBrdrY = math.ceil(vertPixBrdrY)

  if PhactoriDbg():
    myDebugPrint3("pixel border: " + str(vertPixBrdrX) + ", " + str(vertPixBrdrY) + "\n")
  tanInFov = math.tan(math.radians(inFov*0.5))
  if PhactoriDbg():
    myDebugPrint3("tanInFov: " + str(tanInFov) + "\n")

  #vertPixY / ss = tan(inFov*0.5)
  #1/ss = tan(inFov*0.5)/vertPixY
  #(vertPixY - vertPixBrdrY) / ss = tan(adjFov*0.5)
  #(vertPixY - vertPixBrdrY) * tan(inFov*0.5)/vertPixY = tan(adjFov*0.5)

  tanAdjustedFovV = (vertPixY - vertPixBrdrY) * tanInFov / vertPixY
  if PhactoriDbg():
    myDebugPrint3("tanAdjustedFovV: " + str(tanAdjustedFovV) + "\n")
  fovV = 2.0 * math.degrees(math.atan(tanAdjustedFovV))
  if PhactoriDbg():
    myDebugPrint3("fovV: " + str(fovV) + "\n");

  pixelSizeWithBorder = [inXyPixelSize[0] - 2 * int(vertPixBrdrX),
                         inXyPixelSize[1] - 2 * int(vertPixBrdrY)]

  #we are special casing when the pixel ratio is the same for both
  #X and Y for backwards compatibility of test images;  We should
  #take out this special case and update the test images.  The
  #difference is very minor, but detectable
  if inImageSettings.mPixelBorderRatioXY[0] == \
      inImageSettings.mPixelBorderRatioXY[1]:
    #this is the special case--use the original pixel size to calculate
    #aspect ratio, which will be very slightly different than the other
    #case due to the math.ceil clamping to an integer value
    aspectRatio = float(inXyPixelSize[0]) / float(inXyPixelSize[1])
  else:
    aspectRatio = float(pixelSizeWithBorder[0]) / float(pixelSizeWithBorder[1])


  if PhactoriDbg():
    myDebugPrint3("bordered aspectRatio: " + str(aspectRatio) + "\n");
  tanFovH = aspectRatio * tanAdjustedFovV
  fovH = 2.0 * math.degrees(math.atan(tanFovH))
  if PhactoriDbg():
    myDebugPrint3("fovH: " + str(fovH) + "\n");

  #find point projected onto plane defined by inFocalPoint and sideVector,
  #which will be the test point to test vertical FOV against
  dist1 = 0.0
  for ii in range(0,2):
    for jj in range(0,2):
      for kk in range(0,2):
        testpoint = [inBounds[ii+0], inBounds[jj+2], inBounds[kk+4]]
        testdist = CalcRelativeCameraDistance2_BB(inFocalPoint,
                     normCameraDir, normUpVector, sideVector, fovV, fovH,
                     testpoint)
        if testdist > dist1:
          dist1 = testdist
  myDebugPrint("CalcRelativeCameraDistance2 returning dist: " + str(dist1) + "\n");
  return dist1

def SetParaViewRepresentationCameraParams(inXParaViewRenderView, inCameraC,
    inLookDirection, inImageSettingsX, inParaViewSource):
  lookDirection = inLookDirection

  myViewBounds = None
  viewBoundsIo = [myViewBounds]

  #special case--if inCameraC was set up so the camera is at a specified
  #point, if the look at point is specified we need to base the look
  #direction on the camera position and the look at point, and if
  #the look direction is specified, we need to base the look at point
  #on the camera position and the look direction
  if inCameraC.mUseCameraAtPointFlag:
    if PhactoriDbg():
      myDebugPrint3("SetParaViewRepresentationCameraParams: special case of camera at point setting\n")
    #camera at point was specified--do the right stuff
    cameraPoint = inCameraC.mCameraAtPointInfo.\
        GetCurrentGeometricPointWithDisplacement(
            inParaViewSource, viewBoundsIo, True)
    if PhactoriDbg():
      myDebugPrint3("cameraPoint: " + str(cameraPoint) + "\n")
    inCameraC.mLookAtDistanceType = 'absolute'
    if inCameraC.mLookDirectionSpecifiedFlag:
      if PhactoriDbg():
        myDebugPrint3("look direction was specified\n")
      #look direction was specified, so calculate focal point based on camera
      #position and look direction and set look at distance to 1.0
      lookDirection = inCameraC.mLookDirection
      aa = lookDirection[0]
      bb = lookDirection[1]
      cc = lookDirection[2]
      lookDirectionMagnitude = math.sqrt(aa*aa + bb*bb + cc*cc)
      inCameraC.mLookAtDistance = 1.0
      focalPoint = [cameraPoint[0] + (aa / lookDirectionMagnitude),
                    cameraPoint[1] + (bb / lookDirectionMagnitude),
                    cameraPoint[2] + (cc / lookDirectionMagnitude)]
      if PhactoriDbg():
        myDebugPrint3("focalPoint: " + str(focalPoint) + "\n")
      if PhactoriDbg():
        myDebugPrint3("lookDirection: " + str(lookDirection) + "\n")
      if PhactoriDbg():
        myDebugPrint3("inCameraC.mLookAtDistance is 1.0\n")
    else:
      #look direction not specified, so we use look at point.  Calculate
      #look direction and look at distance based on camera position and
      #look at point
      if PhactoriDbg():
        myDebugPrint3("use look at point as look direction was not specified\n")
      focalPoint = inCameraC.mLookAtPointInfo. \
          GetCurrentGeometricPointWithDisplacement(
          inParaViewSource, viewBoundsIo, True)
      if PhactoriDbg():
        myDebugPrint3("focalPoint: " + str(focalPoint) + "\n")
      lookDirection = [focalPoint[0] - cameraPoint[0],
                      focalPoint[1] - cameraPoint[1],
                      focalPoint[2] - cameraPoint[2]]
      if PhactoriDbg():
        myDebugPrint3("lookDirection: " + str(lookDirection) + "\n")
      aa = lookDirection[0]
      bb = lookDirection[1]
      cc = lookDirection[2]
      inCameraC.mLookAtDistance = math.sqrt(aa*aa + bb*bb + cc*cc)
      if PhactoriDbg():
        myDebugPrint3("inCameraC.mLookAtDistance is " + str(inCameraC.mLookAtDistance) + "\n")
      #check for zero magnitude
      if inCameraC.mLookAtDistance != 0.0:
        inCameraC.mLookDirection = lookDirection
      else:
        if PhactoriDbg():
          myDebugPrint3("warning: camera position and look at point are the same\n")
        #try to use previous
        lookDirection = inCameraC.mLookDirection
        aa = lookDirection[0]
        bb = lookDirection[1]
        cc = lookDirection[2]
        inCameraC.mLookAtDistance = math.sqrt(aa*aa + bb*bb + cc*cc)
        if PhactoriDbg():
          myDebugPrint3("  coincident lookDirection: " + str(lookDirection) + "\n")
        if PhactoriDbg():
          myDebugPrint3("  coincident inCameraC.mLookAtDistance is " + str(inCameraC.mLookAtDistance) + "\n")
        if inCameraC.mLookAtDistance == 0.0:
          if PhactoriDbg():
            myDebugPrint3("  warning: previous camera look direction had zero magnitude\n")
          #previous also zero, use -1,-1,-1
          lookDirection = [-1.0, -1.0, -1.0]
          inCameraC.mLookAtDistance = math.sqrt(3.0)
          if PhactoriDbg():
            myDebugPrint3("    coincident lookDirection 2: " + str(lookDirection) + "\n")
          if PhactoriDbg():
            myDebugPrint3("    coincident inCameraC.mLookAtDistance 2 is " + str(inCameraC.mLookAtDistance) + "\n")
  else:
    UpdatePipelineWithCurrentTimeArgument(inParaViewSource)
    focalPoint = inCameraC.mLookAtPointInfo. \
        GetCurrentGeometricPointWithDisplacement(
        inParaViewSource, viewBoundsIo, True)

  if PhactoriDbg():
    myDebugPrint3('  lookDirection ' + str(lookDirection) + '\n')

  myViewBounds = viewBoundsIo[0]

  if PhactoriDbg():
    myDebugPrint3('  focalPoint ' + str(focalPoint) + '\n')

  #need to handle relative and error

  global gDefaultLookatDistanceDatasizeRelativeMultiplier
  global gDefaultLookatDistanceDatasizeRelative
  if inCameraC.mLookAtDistanceType == 'datasize relative':
    cameraPositionFactor = inCameraC.mLookAtDistance
    #cameraPositionFactor *= gDefaultLookatDistanceDatasizeRelativeMultiplier
    cameraToFocalPointDistance = None
    #we'll need view bounds for distance
    if myViewBounds == None:
      myViewBounds = GetGlobalDataBoundsParallel(inParaViewSource)
  #elif inCameraC.mLookAtDistanceType == 'absolute':
  else:
    cameraPositionFactor = gDefaultLookatDistanceDatasizeRelative
    cameraPositionFactor *= gDefaultLookatDistanceDatasizeRelativeMultiplier
    cameraToFocalPointDistance = inCameraC.mLookAtDistance


  if PhactoriDbg():
    myDebugPrint3('  cameraPositionFactor: ' + str(cameraPositionFactor) + '\n')
  if PhactoriDbg():
    myDebugPrint3('  cameraToFocalPointDistance: ' + str(cameraToFocalPointDistance) + '\n')


  SetCameraLookAtPointAndLookDirection(inParaViewRenderView = inXParaViewRenderView,
    inLookDirection = lookDirection,
    inImageSettings = inImageSettingsX,
    inFocalPoint = focalPoint,
    inEyePositionFactor = cameraPositionFactor,
    inEyeToFocalPointDistance = cameraToFocalPointDistance,
    inCameraViewAngle = inCameraC.mViewAngle,
    inViewUp = inCameraC.mViewUpVector,
    inViewBounds = myViewBounds,
    inPhactoriCamera = inCameraC)


def SetForCorrectColorByVariable(inImagesetInfo, inPhactoriOperation,
        inPvDataRepresentation, inPhactoriRepresentation,
        inInitializeColorLegendFlag):
  """sets the color by variable to a scalar, vector component/magnitude, or
     tensor component.  Assumes we've already established we are doing this
     instead of solid color or color by block"""

  if PhactoriDbg():
    myDebugPrint3("SetForCorrectColorByVariable entered\n")

  savedActiveSourceX = GetActiveSource()
  SetActiveSource(inPhactoriOperation.GetPvFilter())

  colorVarInfo = inPhactoriRepresentation.mColorVariableInfo

  if PhactoriDbg():
    myDebugPrint3("colorVarInfo.mVariableName: " + \
        str(colorVarInfo.mVariableName) + "\n")

  if colorVarInfo.mVariableName == '':
    myVariableArrayType = gCellsString
  else:
    global gPipeAndViewsState
    #the correct behavior is probably to use the imageset's filter as the
    #filter to use to determine info about the variable, such as POINT or
    #cell; however prior work used the pipeline start default incoming
    #filter.  We are using prior to get tests to pass for now, and
    #we are considering switching to 'correct.'
    useImagesetFilterToDetectVariableType = False
    if useImagesetFilterToDetectVariableType:
      paraViewSource = inPhactoriOperation.GetPvFilter()
    else:
      paraViewSource = \
        gPipeAndViewsState.mIncomingDefaultOperation.GetPvFilter()
    detectResult = colorVarInfo.DetectVariableType(
        paraViewSource, True, False)

    if detectResult == True:
      #variable type detected determine value for paraview
      if colorVarInfo.mVariableType == 'node':
        myVariableArrayType = gPointsString
      elif colorVarInfo.mVariableType == 'element':
        myVariableArrayType = gCellsString
      else:
        if PhactoriDbg():
          myDebugPrint3(errStr)
        errStr = 'CreateParaviewItemsForImagesetC error:\n'\
          'image set ' + inImagesetInfo.mName + \
          'operation ' + inPhactoriOperation.mName + \
          ' representation ' + inPhactoriRepresentation.mName + \
          '\ncolor variable name: ' + \
          colorVarInfo.mVariableName +\
          '\nvariable should be node or element at this point and is not'
        raise Exception(errStr)
    else:
      #variable type not detected, deal with it
      errStr = 'CreateParaviewItemsForImagesetC error:\n'\
          'image set ' + inImagesetInfo.mName + \
          'operation ' + inPhactoriOperation.mName + \
          ' representation ' + inPhactoriRepresentation.mName + \
          '\ncolor variable name: ' + \
          colorVarInfo.mVariableName +\
          '\nvariable is not point data or cell data, assuming cell'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      #current hack--if detection of variable type fails due to missing
      #variable; assume it is CELLS and let mpi stuff share info
      #and create odd color map otherwise
      myVariableArrayType = gCellsString

  if PhactoriDbg():
    myDebugPrint3("myVariableArrayType: " + myVariableArrayType + "\n")
      #("showColorLegend: " + str(showColorLegend) + "\n")

  colorVarName = colorVarInfo.mVariableName
  if colorVarName != '':
    if colorVarInfo.mVariableIsVectorComponent:
      if PhactoriDbg():
        myDebugPrint3("mVariableIsVectorComponent is true so we do\n"
          "ColorByVariableComponentOrMagnitudeXX (component)\n")
      ColorByVariableComponentOrMagnitudeXX(
        inPvDataRepresentation, inPhactoriRepresentation, colorVarName,
        'Component', colorVarInfo.mVariableComponent, myVariableArrayType,
        inPhactoriRepresentation.mColorMapSettings)
    elif colorVarInfo.mVariableIsVectorMagnitude:
      if PhactoriDbg():
        myDebugPrint3("mVariableIsVectorMagnitude is true so we do\n"
          "ColorByVariableComponentOrMagnitudeXX (magnitude)\n")
      ColorByVariableComponentOrMagnitudeXX(
        inPvDataRepresentation, inPhactoriRepresentation, colorVarName,
        'Magnitude', colorVarInfo.mVariableComponent, myVariableArrayType,
        inPhactoriRepresentation.mColorMapSettings)
    else:
      if PhactoriDbg():
        myDebugPrint3("not component or magnitude so we do "
            "ColorByVariableScalarXX\n")
      ColorByVariableScalarXX(inPvDataRepresentation, inPhactoriRepresentation,
          colorVarName, myVariableArrayType,
          inPhactoriRepresentation.mColorMapSettings)

    #if inInitializeColorLegendFlag:
    if True:
      if inPhactoriRepresentation.mUseFixedColorRange:
        fixedColorRange = inPhactoriRepresentation.mFixedColorRange
        SetMinimumMaximumColorValues(inPvDataRepresentation,
            fixedColorRange[0], fixedColorRange[1],
            inPhactoriOperation, inPhactoriRepresentation,
            colorVarName, myVariableArrayType)
      #else:
      #  SetMinimumMaximumColorValues(inPvDataRepresentation,
      #      -0.1, 0.1,
      #      inPhactoriRepresentation, colorVarName, myVariableArrayType)
  else:
    if gParaViewCatalystVersionFlag <= 40100:
      inPvDataRepresentation.ColorArrayName = ""
    else:
      #if PhactoriDbg():
      #  myDebugPrint3("using ColorBy() to set no color var 1\n")
      ColorBy(inPvDataRepresentation, (gPointsString,''))

  if PhactoriDbg():
    myDebugPrint3("SetForCorrectColorByVariable returning\n")

  #inPvDataRepresentation.RescaleTransferFunctionToDataRange(True)
  SetActiveSource(savedActiveSourceX)

def SetForCorrectColorBy(inImagesetInfo, inPhactoriOperation,
        inPvDataRepresentation, inPhactoriRepresentation,
        inInitializeColorLegendFlag):

  if PhactoriDbg():
    if inPhactoriOperation != None:
        opStr = inPhactoriOperation.mName
    else:
        opStr = "None"
    myDebugPrint3("SetForCorrectColorBy entered\n"
      "inImagesetInfo: " + inImagesetInfo.mName + "\n"
      "inPhactoriOperation: " + opStr + "\n"
      "inPhactoriRepresentation: " + str(inPhactoriRepresentation) + "\n"
      "inInitializeColorLegendFlag: " + str(inInitializeColorLegendFlag) + "\n")
  #color by variable
  if inPhactoriRepresentation != None:
    if inPvDataRepresentation.Opacity != inPhactoriRepresentation.mOpacitySetting:
        inPvDataRepresentation.Opacity = inPhactoriRepresentation.mOpacitySetting
    SetForCorrectColorByVariable(inImagesetInfo, inPhactoriOperation,
        inPvDataRepresentation, inPhactoriRepresentation,
        inInitializeColorLegendFlag)
  else:
    if PhactoriDbg():
      myDebugPrint3("inPhactoriRepresentation is None, so set to no color array\n")
    if gParaViewCatalystVersionFlag <= 40100:
      inPvDataRepresentation.ColorArrayName = ""
    else:
      #if PhactoriDbg():
      #  myDebugPrint3("using ColorBy() to set no color var 2\n")
      ColorBy(inPvDataRepresentation, (gPointsString,''))

  #color by block id or solid color
  if inPhactoriRepresentation != None:
    if inPhactoriRepresentation.mColorByBlockFlag == False:
      #clear block color
      inPvDataRepresentation.BlockColor = {}
    if inPhactoriRepresentation.mColorByBlockFlag:
      if PhactoriDbg():
        myDebugPrint3("inPhactoriRepresentation.mColorByBlockFlag is true\n")

      if gParaViewCatalystVersionFlag <= 40100:
        inPvDataRepresentation.ColorArrayName = ""
      else:
        #if PhactoriDbg():
        #  myDebugPrint3("using ColorBy() to set no color var 3\n")
        ColorBy(inPvDataRepresentation, (gPointsString,''))

      ColorByBlock(GetActiveSource(),
          inPvDataRepresentation,
          inPhactoriRepresentation.mColorByBlockExplicitlySet)
    elif inPhactoriRepresentation.mColorBySolidColorFlag:
      if PhactoriDbg():
        myDebugPrint3("inPhactoriRepresentation.mColorBySolidColorFlag is true\n")
      if gParaViewCatalystVersionFlag <= 40100:
        inPvDataRepresentation.ColorArrayName = ""
      else:
        #if PhactoriDbg():
        #  myDebugPrint3("using ColorBy() to set no color var 4\n")
        ColorBy(inPvDataRepresentation, (gPointsString,''))

      inPvDataRepresentation.DiffuseColor = \
          inPhactoriRepresentation.mSolidColor
    else:
      if PhactoriDbg():
        myDebugPrint3("no color by block or solid color flag\n")

  if PhactoriDbg():
    myDebugPrint3("SetForCorrectColorBy returning\n")

def SetUpOneParaViewRepresentationAndViewC(inCameraC, inLookDirection,
  inImagesetInfo,
  inColorSettingsX,
  inMeshRenderControl,
  inShowDataCubeAxes, inShowDataCubeAxesInfo, inShowOrientationAxes,
  inFixedColorRange,
  inIsPlotFlag,
  inViewBounds = None,
  inRepresentationFilenameAddon = "",
  inLookDirectionFilenameAddon = "",
  inPhactoriRepresentation = None):

  inImageSettings = inImagesetInfo.mImageSettings
  lclImageSize = inImageSettings.mImageSize
  lclImageFormat = inImageSettings.mImageFormat
  inImageBasename = inImageSettings.mImageBasename
  inImageBasedirectory = inImageSettings.mImageBasedirectory
  inNumCounterDigits = inImageSettings.mNumCounterDigits

  #myDebugPrint3('inCameraC: ' + str(inCameraC) + '\n')
  #myDebugPrint3('inLookDirection: ' + str(inLookDirection) + '\n')
  #myDebugPrint3('inImageBasename: ' + str(inImageBasename) + '\n')
  #myDebugPrint3('inImageBasedirectory: ' + str(inImageBasedirectory) + '\n')
  #myDebugPrint3('inNumCounterDigits: ' + str(inNumCounterDigits) + '\n')
  #myDebugPrint3('lclImageFormat: ' + str(lclImageFormat) + '\n')
  #myDebugPrint3('lclImageSize: ' + str(lclImageSize) + '\n')
  #myDebugPrint3('inMeshRenderControl: ' + str(inMeshRenderControl) + '\n')
  #myDebugPrint3('inShowDataCubeAxes: ' + str(inShowDataCubeAxes) + '\n')
  #myDebugPrint3('inShowOrientationAxes: ' + str(inShowOrientationAxes) + '\n')
  #myDebugPrint3('inViewBounds: ' + str(inViewBounds) + '\n')
  #myDebugPrint3('inRepresentationFilenameAddon: ' + str(inRepresentationFilenameAddon) + '\n')
  #myDebugPrint3('inLookDirectionFilenameAddon: ' + str(inLookDirectionFilenameAddon) + '\n')


  #need to handle absolute and error

  nameSuffix = inRepresentationFilenameAddon + inCameraC.mFilenameAddon + inLookDirectionFilenameAddon

  if inImageBasedirectory == None:
    fileBaseName = inImageBasename + nameSuffix
  elif inImageBasedirectory == "":
    fileBaseName = inImageBasename + nameSuffix
  else:
    import os
    fileBaseName = inImageBasedirectory + os.sep + inImageBasename + nameSuffix

  AddRenderView(inPhactoriImagesetInfo = inImagesetInfo,
    inColorSettings = inColorSettingsX,
    ImageBaseFileName = fileBaseName, ImageType = lclImageFormat,
    ImageOverwriteFlag = False,
    PixelSizeX = lclImageSize[0], PixelSizeY = lclImageSize[1],
    numCounterDigits = inNumCounterDigits,
    inIsPlotFlag = inIsPlotFlag)

  if inIsPlotFlag == False:
    SetParaViewRepresentationCameraParams(inImagesetInfo.mSharedPvRenderView2,
        inCameraC, inLookDirection, inImageSettings, GetActiveSource())
    SetForCorrectColorBy(inImagesetInfo,
        inImagesetInfo.GetInputPhactoriOperation(),
        inImagesetInfo.mPvDataRepresentation2, inPhactoriRepresentation, True)

    inImagesetInfo.mPvDataRepresentation2.Representation = inMeshRenderControl
    if(inShowOrientationAxes):
      inImagesetInfo.mSharedPvRenderView2.OrientationAxesVisibility = 1
    else:
      inImagesetInfo.mSharedPvRenderView2.OrientationAxesVisibility = 0

    UpdateRepresentationColorBy(inImagesetInfo)
    #if inShowColorLegend:
    #  if PhactoriDbg():
    #    myDebugPrint3('  CreateOneCameraViewFromViewMapCInfo color legend on\n')
    #  ShowDataColorLegendXX(inImagesetInfo, 'on',
    #      inColorLegendPositionAndSize, inColorSettingsX)
    #else:
    #  if PhactoriDbg():
    #    myDebugPrint3('  CreateOneCameraViewFromViewMapCInfo color legend off\n')
    #  ShowDataColorLegendXX(inImagesetInfo, 'off',
    #      inColorLegendPositionAndSize, inColorSettingsX)

    if inShowDataCubeAxes:
      ShowCubeAxesXX(inImagesetInfo.mSharedPvRenderView2, 'on', inShowDataCubeAxesInfo)
    else:
      ShowCubeAxesXX(inImagesetInfo.mSharedPvRenderView2, 'off')

def GetLookDirectionListFromCamera(inCamera):
  if inCamera.mType == 'multicamera8':
    directionsList = [ [[-1.0,  0.0,  0.0], 'x1.'],
                       [[ 1.0,  0.0,  0.0], 'x2.'],
                       [[ 0.0, -1.0,  0.0], 'y1.'],
                       [[ 0.0,  1.0,  0.0], 'y2.'],
                       [[ 0.0,  0.0, -1.0], 'z1.'],
                       [[ 0.0,  0.0,  1.0], 'z2.'],
                       [[-1.0, -1.0, -1.0], 'xyz1.'],
                       [[ 1.0,  1.0,  1.0], 'xyz2.'] ]
  elif inCamera.mType == 'camera':
    directionsList = [[inCamera.mLookDirection, '']]
  else:
    errStr = 'error!  GetLookDirectionListFromCamera has illegal camera type\n'
    if PhactoriDbg():
      myDebugPrint3(errStr)
    raise Exception(errStr)

  return directionsList


def ParseOneFilterTypeFromViewMapOperation(ioOperationBlock, inTypeString, inOperationClass, inOperationParamsJson):
  ioOperationBlock.mType = inTypeString

  ioOperationBlock.mOperationSpecifics = inOperationClass()
  ioOperationBlock.mOperationSpecifics.mPhactoriOperationBlockOwner = \
    ioOperationBlock
  ioOperationBlock.mOperationSpecifics.ParseParametersFromJson(inOperationParamsJson)

  #if 'input operation' in inOperationParamsJson:
  if 'input' in inOperationParamsJson:
    ioOperationBlock.mInputOperationName = inOperationParamsJson['input']
    if PhactoriDbg(100):
      myDebugPrint3("ParseOneFilterTypeFromViewMapOperation: input is: \n" + ioOperationBlock.mInputOperationName, 100)
  else:
    noticeStr = 'notice!  inOperationParamsJson has no input key, using default pipeline input\n'
    ioOperationBlock.mInputOperationName = None
    if PhactoriDbg():
      myDebugPrint3(noticeStr)


def ConstructClipPlaneOperationFromParsedOperationBlockC(ioPipeAndViewsState, ioOperationBlock):
  return


def ConstructPipelineOperationFromParsedOperationBlockC(ioPipeAndViewsState, ioOperationBlock):

  if PhactoriDbg(100):
    myDebugPrint3("ConstructPipelineOperationFromParsedOperationBlockC entered\n", 100)
  if PhactoriDbg():
    myDebugPrint3("  constructing operation named " + ioOperationBlock.mName + "\n")

  particularOperation = ioOperationBlock.mOperationSpecifics

  if ioOperationBlock.mType == 'group':
    if PhactoriDbg(100):
      myDebugPrint3("type was 'group', different pv filter construction\n",
      100)
    newParaViewFilter = particularOperation.CreateParaViewFilter2(
        ioPipeAndViewsState)
  else:
    inputSource = None
    if ioOperationBlock.mInputOperationName == None:
      inputSource = GetActiveSource()
      if PhactoriDbg():
        myDebugPrint3("mInputOperationName was none so using" + str(inputSource) + '\n')
    else:
      inputOperationBlock = ioPipeAndViewsState.mOperationBlocks[
          ioOperationBlock.mInputOperationName]
      inputSource = inputOperationBlock.GetPvFilter()
      if PhactoriDbg():
        myDebugPrint3("mInputOperationName was " + ioOperationBlock.mInputOperationName + " so using" + str(inputSource) + '\n')

    newParaViewFilter = particularOperation.CreateParaViewFilter(inputSource)

  ioOperationBlock.mParaViewFilter = newParaViewFilter

  if PhactoriDbg(100):
    myDebugPrint3("ConstructPipelineOperationFromParsedOperationBlockC returning\n", 100)


def MakeFiltersFromViewMapOperationsC(ioPipeAndViewsState, inOperationBlocksJson):
  if PhactoriDbg(100):
    myDebugPrint3('MakeFiltersFromViewMapOperationsC entered\n', 100)

  #loop through set of filter blocks (in json dict) and parse out the filter objects
  for operationName, operationParams in inOperationBlocksJson.iteritems():

    #hack to insert filter images in/out
    if operationName == "ImageFilteringStartOperationOverride":
      if PhactoriDbg():
        myDebugPrint3("overriding operation to do image set on/off filtering\n")
      operationParams['type'] = 'image_filtering_start'
      ParseOneImageStartStopFilterFromViewMap(operationParams)
      #then go on to next operation
      continue

    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = operationName
    if PhactoriDbg():
      myDebugPrint3('  creating operation named ' + newOperationBlock.mName + '\n')

    if 'type' not in operationParams:
      errStr = 'error!  operation block with name ' + str(operationName) + ' has no type key in MakeFiltersFromViewMapOperationsC\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    if PhactoriDbg():
        myDebugPrint3('  operation type: ' + str(operationParams['type']) + '\n')

    if operationParams['type'] == 'threshold':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'threshold',
              PhactoriThresholdOperation,
              operationParams)
    elif operationParams['type'] == 'extractblock':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'extractblock',
              PhactoriExtractBlockOperation,
              operationParams)
    elif operationParams['type'] == 'generatesurfacenormals':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'generatesurfacenormals',
              PhactoriGenerateSurfaceNormalsOperation,
              operationParams)
    elif operationParams['type'] == 'extractsurface':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'extractsurface',
              PhactoriExtractSurfaceOperation,
              operationParams)
    elif operationParams['type'] == 'ghostcellsgenerator':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'ghostcellsgenerator',
              PhactoriGhostCellsGeneratorOperation,
              operationParams)
    elif operationParams['type'] == 'mergeblocks':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'mergeblocks',
              PhactoriMergeBlocksOperation,
              operationParams)
    elif operationParams['type'] == 'subdivide':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'subdivide',
              PhactoriSubdivideOperation,
              operationParams)
    elif operationParams['type'] == 'triangulate':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'triangulate',
              PhactoriTriangulateOperation,
              operationParams)
    elif operationParams['type'] == 'group':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'group',
              PhactoriGroupOperation,
              operationParams)
    elif operationParams['type'] == 'transform':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'transform',
              PhactoriTransformOperation,
              operationParams)
    elif operationParams['type'] == 'calculator':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'calculator',
              PhactoriCalculatorOperation,
              operationParams)
    elif operationParams['type'] == 'add point set':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'add point set',
              PhactoriAddPointSetOperation,
              operationParams)
    elif operationParams['type'] == 'add unstructured grid':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'add unstructured grid',
              PhactoriAddUnstructuredGridOperation,
              operationParams)
    elif operationParams['type'] == 'reflect':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'reflect',
              PhactoriReflectOperation,
              operationParams)
    elif operationParams['type'] == 'warpbyvector':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'warpbyvector',
              PhactoriWarpByVectorOperation,
              operationParams)
    elif operationParams['type'] == 'clip':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'clip',
              PhactoriClipPlaneOperation,
              operationParams)
    elif operationParams['type'] == 'slice':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'slice',
              PhactoriSliceOperation,
              operationParams)
    elif operationParams['type'] == 'boxclip':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'boxclip',
              PhactoriBoxClipOperation,
              operationParams)
    elif operationParams['type'] == 'cylinderclip':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'cylinderclip',
              PhactoriCylinderClipOperation,
              operationParams)
    elif operationParams['type'] == 'element data to node data':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'element data to node data',
              PhactoriCellDataToPointDataOperation,
              operationParams)
    elif operationParams['type' ] == 'nearestpoints':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'nearestpoints',
              PhactoriNearestPointsOperation,
              operationParams)
    elif operationParams['type'] == 'castnormalrays':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'castnormalrays',
              PhactoriIntersectNodeNormalsWithSurface,
              operationParams)
    elif operationParams['type'] == 'exportvtp':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'exportvtp',
              PhactoriVtpDataExporterOperation,
              operationParams)
    elif operationParams['type'] == 'exportvtm':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'exportvtm',
              PhactoriVtmDataExporterOperation,
              operationParams)
    elif operationParams['type'] == 'contour':
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'contour',
              PhactoriContourOperation,
              operationParams)
    else:
      errStr = 'error!  in MakeFiltersFromViewMapOperationsC inOperationBlocksJson operation type is unrecognized ' + str(operationParams['type']) + '\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    ioPipeAndViewsState.mOperationBlocks[newOperationBlock.mName] = \
        newOperationBlock

  #check to make sure no operation block calls for an input which doesn't exist
  for operationName, operationBlock in ioPipeAndViewsState.mOperationBlocks.iteritems():
    if operationBlock.mInputOperationName != None:
      if operationBlock.mInputOperationName not in ioPipeAndViewsState.mOperationBlocks:
        errStr = 'error! in MakeFiltersFromViewMapOperationsC operation block with name ' + str(operationName) + ' calls for nonexistent input operation with name ' + str(operationBlock.mInputOperationName) + '\n'
        if PhactoriDbg():
          myDebugPrint3(errStr)
        raise Exception(errStr)

  #now construct each filter in pipeline, being careful to construct them in a valid order, which causes multiple passes through the set
  #first, mark all as uncreated
  for operationName, operationBlock in ioPipeAndViewsState.mOperationBlocks.iteritems():
    operationBlock.mHasBeenConstructed = False

  #we need to keep looping through the operation blocks until all have been constructed
  keepConstructing = True
  while keepConstructing:
    keepConstructing = False  #at the beginning of the loop, we haven't seen any unconstructed operations
    for operationName, operationBlock in ioPipeAndViewsState.mOperationBlocks.iteritems():
      if operationBlock.mHasBeenConstructed == False:
        #this one still needs constructing, and we'll have to loop again to
        #make sure all are constructed
        keepConstructing = True

        #determine whether or not we can construct this operation yet, or if
        #we need to wait for something else to be constructed
        canBeConstructedNow = False
        if operationBlock.mType != "group":
          if operationBlock.mInputOperationName == None:
            canBeConstructedNow = True
          else:
            inputBlock = ioPipeAndViewsState.mOperationBlocks[operationBlock.mInputOperationName]
            if inputBlock.mHasBeenConstructed == True:
              canBeConstructedNow = True
        else:
          #special case--for group operation to be constructed, we need all
          #input operations to be instructed
          canBeConstructedNow = True
          for opName in operationBlock.mOperationSpecifics.mOperationNameList:
            inputBlock = ioPipeAndViewsState.mOperationBlocks[opName]
            if inputBlock.mHasBeenConstructed == False:
              canBeConstructedNow = False
              break

        if canBeConstructedNow == True:
          #this operation's input has been constructed (or is None), so we can construct it
          ConstructPipelineOperationFromParsedOperationBlockC(ioPipeAndViewsState, operationBlock)
          operationBlock.mHasBeenConstructed = True



#helper method:  given a block A which potentially contains a key B with
#value BV, see if block set CS has key BV which indicates a refered-to
#block.  Return None if we can't get referred block, warn if A had key
#but BV didn't exist in CS (A not having B merely returns None with no
#complaint
def FindReferredBlockC(referringBlock, referringKey, referredBlockset):
  #myDebugPrint3('  FindReferredBlockC entered\n',100)
  #myDebugPrint3('  referringBlock: ' + str(referringBlock) + '\n')
  #myDebugPrint3('  referringKey: ' + str(referringKey) + '\n')
  #myDebugPrint3('  referredBlockset: ' + str(referringKey) + '\n')
  if referringKey not in referringBlock:
    return None
  referredName = referringBlock[referringKey]
  if referredName in referredBlockset:
    return referredBlockset[referredName]

  if PhactoriDbg():
    myDebugPrint3('  warning:  ' + referringKey + ' was expected to be in referredBlockset but was not, returning None\n')
  return None

def UpdatePipelineWithCurrentTimeArgument(inParaViewFilter):
  global gPipeAndViewsState
  thisTime = gPipeAndViewsState.CurrentDatadescription.GetTime()
  inParaViewFilter.UpdatePipeline(thisTime)

def SetParaViewViewToCurrentTime(inParaViewView):
  global gPipeAndViewsState
  thisTime = gPipeAndViewsState.CurrentDatadescription.GetTime()
  inParaViewView.ViewTime = thisTime

def CreateParaViewRepresentationAndViewFromInfoC(inImageset, inLookDirection, inLookDirectionFilenameAddonL):
  global gPipeAndViewsState

  if PhactoriDbg(100):
    myDebugPrint3("CreateParaViewRepresentationAndViewFromInfoC entered\n",100)
  if PhactoriDbg():
    myDebugPrint3("  imageset: " + inImageset.mName + "\n")
  if PhactoriDbg():
    myDebugPrint3("  lookdir : " + str(inLookDirection) + "  file name addon: " + str(inLookDirectionFilenameAddonL) + "\n")

  newLookDirection = list(inLookDirection)
  inImageset.mLookDirectionList.append(newLookDirection)
  inImageset.mLookDirectionFilenameAddon.append(inLookDirectionFilenameAddonL)

  theRepresentation = inImageset.mRepresentation
  theOperation = inImageset.mOperation
  inCamera = inImageset.mCamera

  #we are setting this up, for now, so that you MUST have a representation at this
  #point--default should have been added and referenced earlier
  if theRepresentation == None:
    errStr = 'error!  theRepresentation is None in CreateParaviewItemsForImagesetC\n'
    if PhactoriDbg():
      myDebugPrint3(errStr)
    raise Exception(errStr)

  meshRenderControl = theRepresentation.mMeshRenderControl

  showColorLegend = theRepresentation.mColorLegendFlag
  colorLegendPositionAndSize = theRepresentation.mColorLegendPositionAndSize
  showDataCubeAxes = theRepresentation.mDataCubeAxesFlag
  showOrientationAxes = theRepresentation.mOrientationAxesFlag
  if theRepresentation.mUseFixedColorRange == True:
    fixedColorRange = theRepresentation.mFixedColorRange
  else:
    fixedColorRange = None

  savedActiveSource = GetActiveSource()
  if PhactoriDbg():
    myDebugPrint3("  operation is " + theOperation.mName + \
        " with ParaView Filter " + str(theOperation.GetPvFilter()) + "\n")
  pvPvGeomFilterFromOp = theOperation.GetOutgoingPvGeometryFilter()
  if PhactoriDbg():
    myDebugPrint3("CreateParaViewRepresentationAndViewFromInfoC:\n"
      "setting active source to exterior geometry filter from operation\n"
      "active source A: " + str(GetActiveSource()) + "\n")
  SetActiveSource(pvPvGeomFilterFromOp)
  if PhactoriDbg():
    myDebugPrint3("active source B1: " + str(pvPvGeomFilterFromOp) + "\n")
    myDebugPrint3("active source B: " + str(GetActiveSource()) + "\n")
  UpdatePipelineWithCurrentTimeArgument(pvPvGeomFilterFromOp)

  if PhactoriDbg():
    myDebugPrint3("  operation point data arrays:\n")
    numArrays = pvPvGeomFilterFromOp.PointData.GetNumberOfArrays()
    for ii in range (0, numArrays):
      myDebugPrint3("  " + str(ii) + ":  " + pvPvGeomFilterFromOp.PointData.GetArray(ii).GetName() + "\n")

  if PhactoriDbg():
    myDebugPrint3("  operation cell data arrays:\n")
    numArrays = pvPvGeomFilterFromOp.CellData.GetNumberOfArrays()
    for ii in range (0, numArrays):
      myDebugPrint3("  " + str(ii) + ":  " + pvPvGeomFilterFromOp.CellData.GetArray(ii).GetName() + "\n")

  if theRepresentation.mColorVariableInfo.mVariableName == '':
    myVariableArrayType = gCellsString
  else:
    #paraViewSource = \
    #    gPipeAndViewsState.mIncomingDefaultOperation.GetPvFilter()
    paraViewSource = pvPvGeomFilterFromOp
    detectResult = theRepresentation.mColorVariableInfo.DetectVariableType(
        paraViewSource, True, False)

    if detectResult == True:
      #variable type detected determine value for paraview
      if theRepresentation.mColorVariableInfo.mVariableType == 'node':
        myVariableArrayType = gPointsString
      elif theRepresentation.mColorVariableInfo.mVariableType == 'element':
        myVariableArrayType = gCellsString
      else:
        if PhactoriDbg():
          myDebugPrint3(errStr)
        errStr = 'CreateParaviewItemsForImagesetC error:\n'\
          'image set ' + inImageset.mName + \
          ' representation ' + theRepresentation.mName + \
          '\ncolor variable name: ' + \
          theRepresentation.mColorVariableInfo.mVariableName +\
          '\nvariable should be node or element at this point and is not'
        raise Exception(errStr)
    else:
      #variable type not detected, deal with it
      errStr = 'CreateParaviewItemsForImagesetC error:\n'\
          'image set ' + inImageset.mName + \
          ' representation ' + theRepresentation.mName + \
          '\ncolor variable name: ' + \
          theRepresentation.mColorVariableInfo.mVariableName +\
          '\nvariable is not point data or cell data, assuming cell'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      #current hack--if detection of variable type fails due to missing
      #variable; assume it is CELLS and let mpi stuff share info
      #and create odd color map otherwise
      myVariableArrayType = gCellsString

  global gCameraFilenameSuffixCounter
  gCameraFilenameSuffixCounter = 0
  SetUpOneParaViewRepresentationAndViewC(inImageset.mCamera,
    inLookDirection,
    inImagesetInfo = inImageset,
    inColorSettingsX = inImageset.mRepresentation.mColorSettings,
    inMeshRenderControl = meshRenderControl,
    inShowDataCubeAxes = showDataCubeAxes,
    inShowDataCubeAxesInfo = theRepresentation.mDataCubeAxesInfo,
    inShowOrientationAxes = showOrientationAxes,
    inFixedColorRange = fixedColorRange,
    inIsPlotFlag = False,
    inRepresentationFilenameAddon = theRepresentation.mFilenameAddon,
    inLookDirectionFilenameAddon = inLookDirectionFilenameAddonL,
    inPhactoriRepresentation = theRepresentation)

  global currentPhactoriRenderViewInfo
  #newParaViewRenderInfoC.mParaViewInfo = currentPhactoriRenderViewInfo

  if PhactoriDbg(100):
    myDebugPrint3("CreateParaViewRepresentationAndViewFromInfoC returning\n", 100)

  SetActiveSource(savedActiveSource)

  if 1:
    if PhactoriDbg(100):
      myDebugPrint3("trying time annotation stuff\n", 100)
    timeAnnStngs = inImageset.mRepresentation.mTimeAnnotationSettings
    #since we are now using one RenderView, we can share a single time
    #annotation source representation.  However, it is conceivable that
    #we may eventually want to show different times for some reason, and
    #in that case we will need to have each imageset have its own time
    #annotation paraview stuff
    if timeAnnStngs.mVisible:
      if gPipeAndViewsState.mTimeAnnotationPv == None:
        gPipeAndViewsState.mTimeAnnotationPv = \
            PhactoriAnnotationPv(timeAnnStngs)
        gPipeAndViewsState.mTimeAnnotationPv.CreateParaViewStuff(
            inImageset.mRepresentation.mColorSettings.mTimeAnnotationColor,
            inImageset)

      #if newParaViewRenderInfoC.mTimeAnnotationPv == None:
      #  newParaViewRenderInfoC.mTimeAnnotationPv = \
      #      PhactoriAnnotationPv(timeAnnStngs)
      #  newParaViewRenderInfoC.mTimeAnnotationPv.CreateParaViewStuff(
      #      inImageset.mRepresentation.mColorSettings.mTimeAnnotationColor,
      #      inImageset)


def CreateParaviewItemsForImagesetC(inImageset):
  """given a json block structure as discussed in the catalyst sierra
     insitu wiki, create the corresponding one imageset (view) based on
     a particular imageset block and the availa"""

  theCamera = inImageset.mCamera

  #we are setting this up, for now, so that you MUST have a camera at this
  #point--default should have been added and referenced earlier
  if theCamera == None:
    errStr = 'error!  theCamera is None in "\
        "CreateParaviewItemsForImagesetC\n'
    if PhactoriDbg():
      myDebugPrint3(errStr)
    raise Exception(errStr)

  cameraLookDirectionAndFilenameAddonList = GetLookDirectionListFromCamera(theCamera)

  #counter = 0
  for oneItem in cameraLookDirectionAndFilenameAddonList:
    #print oneItem
    oneDirection = oneItem[0]
    oneDirectionFilenameAddon = oneItem[1]
    CreateParaViewRepresentationAndViewFromInfoC(inImageset, oneDirection,
            oneDirectionFilenameAddon)
    #counter = counter + 1

def ParseOneCriteriaBlockC(ioCriteriaBlock, ioCriteriaBlockJson,
        inPipeAndViewsState):
  ioCriteriaBlock.ParseParametersFromJson(ioCriteriaBlockJson)

def ParseOneMarkerBlockC(ioMarkerBlock, ioBlockJson,
        inPipeAndViewsState):
  ioMarkerBlock.ParseMarkerSettingsFromJson(ioBlockJson)

def ParseOneTextAnnotationBlockC(ioTextAnnotationBlock, ioBlockJson,
        inPipeAndViewsState):
  ioTextAnnotationBlock.ParseTextAnnotationSettingsFromJson(ioBlockJson)

def ParseOneCameraBlockC(ioCameraBlock, ioCameraBlockJson, inPipeAndViewsState):
  #parse camera type
  if 'camera type' not in ioCameraBlockJson:
    errStr = 'error!  camera block with name ' + str(ioCameraBlock.mName) + \
      ' has no camera type key in ParseCameraBlocksC\n'
    if PhactoriDbg():
      myDebugPrint3(errStr)
    raise Exception(errStr)

  cameraType = ioCameraBlockJson['camera type']

  if cameraType != 'camera' and cameraType != 'multicamera8':
    errStr = 'camera type entry is not recognized in camera ' + ioCameraBlock.mName + ' ParseCameraBlocksC\n'
    if PhactoriDbg():
      myDebugPrint3(errStr)
    raise Exception(errStr)

  ioCameraBlock.mType = cameraType

  #check for camera at node|element displaced settings
  ioCameraBlock.mCameraAtPointInfo.UserPointInfoParseParametersFromJson( \
      ioCameraBlockJson, "camera at ", " displaced", True)
  if ioCameraBlock.mCameraAtPointInfo.mParsingDetectedAtLeastOneSetting == False:
    #check for camera at settings (not displaced) if there weren't displaced ones
    ioCameraBlock.mCameraAtPointInfo.UserPointInfoParseParametersFromJson( \
        ioCameraBlockJson, "camera at ", "", False)
  if ioCameraBlock.mCameraAtPointInfo.mParsingDetectedAtLeastOneSetting:
    ioCameraBlock.mUseCameraAtPointFlag = True
  else:
    ioCameraBlock.mUseCameraAtPointFlag = False

  if PhactoriDbg():
    myDebugPrint3("camera parsing look at point:\n")
  ioCameraBlock.mLookAtPointInfo.UserPointInfoParseParametersFromJson(ioCameraBlockJson,
      "look at ", "")

  #parse look at distance
  if 'look at relative distance' in ioCameraBlockJson:
    ioCameraBlock.mLookAtDistanceType = 'datasize relative'
    ioCameraBlock.mLookAtDistance = ioCameraBlockJson['look at relative distance']
    #if 'camera at point' was specified, look at distance should not
    #be specified
  elif 'look at absolute distance' in ioCameraBlockJson:
    ioCameraBlock.mLookAtDistanceType = 'absolute'
    ioCameraBlock.mLookAtDistance = ioCameraBlockJson['look at absolute distance']
    #if 'camera at point' was specified, look at distance should not
    #be specified
  else:
    ioCameraBlock.mLookAtDistanceType = 'datasize relative'
    ioCameraBlock.mLookAtDistance = 1.0

  #parse look direction (single camera, not multicamera8)
  if ioCameraBlock.mType == 'camera':
    ioCameraBlock.mLookDirection = [-1.0, -1.0, -1.0]
    if 'look direction' in ioCameraBlockJson:
      ioCameraBlock.mLookDirection = ioCameraBlockJson['look direction']
      ioCameraBlock.mLookDirectionSpecifiedFlag = True
    else:
      ioCameraBlock.mLookDirectionSpecifiedFlag = False

  #get text to add to each image created with this camera, if any
  if "image name addon" in ioCameraBlockJson:
    ioCameraBlock.mFilenameAddon = ioCameraBlockJson["image name addon"]

  #camera field of view (for default use current value which was set at
  #creation time--see class definition for default
  ioCameraBlock.mViewAngle = getParameterFromBlock(ioCameraBlockJson,
    'camera fov', ioCameraBlock.mViewAngle)

  #up vector for the view (for default use current value which was set at
  #creation time--see class definition for default)
  ioCameraBlock.mViewUpVector = getParameterFromBlock(ioCameraBlockJson,
    'up vector', ioCameraBlock.mViewUpVector)

  #parse off-axis projection info
  ioCameraBlock.mOffAxisProjectionInfo.ParseParametersFromJson(
    ioCameraBlockJson)

  #parse parallel projection info
  perspectiveOrParallel = getParameterFromBlock(ioCameraBlockJson,
    'projection type', 'perspective')
  if perspectiveOrParallel == 'parallel':
    ioCameraBlock.mUseParallelProjection = True
    if PhactoriDbg(100):
      myDebugPrint3("projection type is now parallel\n", 100)
    if 'absolute parallel scale' in ioCameraBlockJson:
      ioCameraBlock.mParallelScaleAbsoluteOrRelative = 0
      ioCameraBlock.mParallelScale = \
        ioCameraBlockJson['absolute parallel scale']
      if PhactoriDbg(100):
        myDebugPrint3("absolute parallel scale: " + \
                str(ioCameraBlock.mParallelScale) + "\n", 100)
    elif 'relative parallel scale' in ioCameraBlockJson:
      ioCameraBlock.mParallelScaleAbsoluteOrRelative = 1
      ioCameraBlock.mParallelScale = \
        ioCameraBlockJson['relative parallel scale']
      if PhactoriDbg(100):
        myDebugPrint3("relative parallel scale: " + \
                str(ioCameraBlock.mParallelScale) + "\n", 100)
    else:
      ioCameraBlock.mParallelScale = 1.0
      ioCameraBlock.mParallelScaleAbsoluteOrRelative = 1
      if PhactoriDbg(100):
        myDebugPrint3("default relative parallel scale: " + \
                str(ioCameraBlock.mParallelScale) + "\n", 100)


def localGet1or0(inJsn, inKey, inDefault):
  value = getParameterFromBlock(inJsn, inKey, inDefault)
  if value:
    return 1
  else:
    return 0


def ParseOneRepresentationBlockC(ioRepresentationBlock, inRepresentationBlockJson, inPipeAndViewsState):
  """given a python dict (presumably from json) description of a representation block, parse all the
     representation parameters out of the block"""

  if PhactoriDbg(100):
    myDebugPrint3('ParseOneRepresentationBlockC entered\n', 100)

  inJsn = inRepresentationBlockJson

  #hack test for colors start
  #inJsn['surface color'] = [1.0, 0.0, 0.0]
  #inJsn['background color'] = [0.0, 0.0, 0.0]
  #inJsn['edge color'] = [0.0, 1.0, 0.0]
  ##inJsn['axes color'] = [0.0, 1.0, 1.0]
  #inJsn['text color'] = [1.0, 1.0, 0.0]
  #hack test for colors end

  #myDebugPrint3('  json: ' + str(inJsn) + '\n')

  if 'point size' in inJsn:
    ioRepresentationBlock.mPointSize = inJsn['point size']

  #color by variable scalar/vector magnitude/vector component/tensor component
  ioRepresentationBlock.mColorVariableInfo.\
    ParseVariableNameAndVectorOrTensorComponent(inJsn, 'color by ')

  if ioRepresentationBlock.mColorVariableInfo.mVariableName != '':
    ioRepresentationBlock.mColorByBlockFlag = False
    ioRepresentationBlock.mColorBySolidColorFlag = False
  elif 'color by blockid' in inJsn:
    ioRepresentationBlock.mColorByBlockFlag = True
    ioRepresentationBlock.mColorByBlockExplicitlySet = True
    ioRepresentationBlock.mColorBySolidColorFlag = False
  elif 'color by solid color' in inJsn:
    ioRepresentationBlock.mColorBySolidColorFlag = True
    ioRepresentationBlock.mColorByBlockFlag = False
    ioRepresentationBlock.mSolidColor = inJsn['color by solid color']

  #color map range control
  if 'color legend range' in inJsn:
    ioRepresentationBlock.mFixedColorRange = inJsn['color legend range']
    ioRepresentationBlock.mUseFixedColorRange = True
  else:
    ioRepresentationBlock.mUseFixedColorRange = False

  #highlight subranges with solid colors
  if 'highlight subrange 1' in inJsn:
    highlightSubrangeIndex = 1
    while True:
      oneSubrangeKey = 'highlight subrange ' + str(highlightSubrangeIndex)
      if oneSubrangeKey in inJsn:
        subrangeArgs = inJsn[oneSubrangeKey]
        highlightSubrangeIndex += 1
        if len(subrangeArgs) != 5:
          if PhactoriDbg():
            myDebugPrint3("highlight subrange needs 5 values\n");
          PrintOnProcessZero("highlight subrange needs 5 values, skipping " + \
                  oneSubrangeKey + "\n");
          continue
        srmin = float(subrangeArgs[0])
        srmax = float(subrangeArgs[1])
        if srmin > srmax:
          if PhactoriDbg():
            myDebugPrint3("subrange highlight min >= max: " + \
              str(srmin) + ", " + str(srmax) + "\nskipping " + \
              oneSubrangeKey + "\n", 100)
          PrintOnProcessZero("subrange highlight min >= max: " + \
            str(srmin) + ", " + str(srmax) + "\nskipping " + \
            oneSubrangeKey + "\n")
          continue
        srColor = [float(subrangeArgs[2]), float(subrangeArgs[3]),
                float(subrangeArgs[4])]
        if (srColor[0] < 0.0) or (srColor[0] > 1.0) or \
           (srColor[1] < 0.0) or (srColor[1] > 1.0) or \
           (srColor[2] < 0.0) or (srColor[2] > 1.0):
          srColor = [1.0, 1.0, 0.0]
          if PhactoriDbg():
            myDebugPrint3(oneSubrangeKey + ": bad color "
              "(component not 0.0-1.0), using rgb 1.0, 1.0, 0.0\n", 100)
          PrintOnProcessZero(oneSubrangeKey + ": bad color "
            "(component not 0.0-1.0), using rgb 1.0, 1.0, 0.0\n")
        ioRepresentationBlock.mHighlightSubranges.append(
                [srmin, srmax, srColor])
        ioRepresentationBlock.mUseHighlightSubranges = True
      else:
        break
    if PhactoriDbg():
      myDebugPrint3("parsed highlight subranges:\n" + \
              str(ioRepresentationBlock.mHighlightSubranges) + "\n", 100)

  #if 'highlight subranges' in inJsn:
  #  ioRepresentationBlock.mUseHighlightSubranges = True
  #  sbrngsJsn = inJsn['highlight subranges']
  #  for oneSubrange in sbrngsJsn:
  #    if 'range' in oneSubrange:
  #      srmin = oneSubrange['range'][0]
  #      srmax = oneSubrange['range'][1]
  #      if srmin < srmax:
  #        if 'color' in oneSubrange:
  #            srColor = oneSubrange['color']
  #        else:
  #            srColor = [1.0, 1.0, 0.0]
  #        ioRepresentationBlock.mHighlightSubranges.append(
  #                [srmin, srmax, srColor])
  #      else:
  #        if PhactoriDbg():
  #          myDebugPrint3("highlight min >= max: " + \
  #            str(srmin) + ", " + str(srmax) + "\nskipping subrange\n", 100)
  #        PrintOnProcessZero("highlight min >= max: " + \
  #          str(srmin) + ", " + str(srmax) + "\nskipping subrange\n")
  #    else:
  #        if PhactoriDbg():
  #          myDebugPrint3("subrange is missing 'range' key; skipping\n")
  #        PrintOnProcessZero("subrange is missing 'range' key; skipping\n")
  #  if PhactoriDbg():
  #    myDebugPrint3("parsed highlight subranges:\n" + \
  #            str(ioRepresentationBlock.mHighlightSubranges) + "\n", 100)

  #additional capability: use ratio-expressed subrange of
  #range that would otherwise be used to increase concentration of
  #dynamic range of color map in range of interest
  if 'color legend subrange' in inJsn:
    ioRepresentationBlock.mUseColorSubrange = True
    ioRepresentationBlock.mColorSubrange = inJsn['color legend subrange']
    goodSubrange = True
    if ioRepresentationBlock.mColorSubrange[0] < 0.0:
      goodSubrange = False
    if ioRepresentationBlock.mColorSubrange[1] > 1.0:
      goodSubrange = False
    if ioRepresentationBlock.mColorSubrange[0] > \
            ioRepresentationBlock.mColorSubrange[1]:
      goodSubrange = False
    if goodSubrange == False:
      myDebugPrint3AndException(
        "ParseOneRepresentationBlockC:\n"
        "bad color legend subrange, must be 0.0 <= bottom <= top <= 1.0\n")

  #choose color map by name
  ioRepresentationBlock.mColorMapSettings.ParseColorMapSettings(inJsn)

  showSurfacesFlag = getParameterFromBlock(inJsn, 'show surfaces', True)
  showEdgesFlag = getParameterFromBlock(inJsn, 'show edges', False)
  showPointsFlag = getParameterFromBlock(inJsn, 'show points', False)

  ioRepresentationBlock.mOpacitySetting = getParameterFromBlock(inJsn, 'opacity',
          ioRepresentationBlock.mOpacitySetting)

  #doVolumeRenderingFlag = getParameterFromBlock(inJsn, 'volume rendering', True)
  doVolumeRenderingFlag = getParameterFromBlock(inJsn, 'volume rendering', False)
  ioRepresentationBlock.mScalarOpacityUnitDistance = getParameterFromBlock(
          inJsn, 'scalar opacity unit distance', -1.0)

  #ioRepresentationBlock.mScalarOpacityUnitDistance = 0.01
  if PhactoriDbg():
      myDebugPrint3("doVolumeRenderingFlag: " + \
              str(doVolumeRenderingFlag) + "\n" + \
              "ioRepresentationBlock.mScalarOpacityUnitDistance: " + \
              str(ioRepresentationBlock.mScalarOpacityUnitDistance) + "\n")

  ioRepresentationBlock.mPresetsImportFileName = getParameterFromBlock(inJsn,
          'color and opacity presets import file', None)
  if ioRepresentationBlock.mPresetsImportFileName != None:
    retval = ImportPresets(ioRepresentationBlock.mPresetsImportFileName)
    if retval != True:
      myDebugPrint3AndException(
        "paraview.simple.ImportPresets failed with the file:\n" + \
        str(ioRepresentationBlock.mPresetsImportFileName) + "\n")
  ioRepresentationBlock.mNameOfPresetToUse = getParameterFromBlock(inJsn,
          'color and opacity preset', None)

  showBoundingBoxFlag = getParameterFromBlock(inJsn, 'show bounding box', False)

  ioRepresentationBlock.mMeshRenderControl = 'Surface'
  if showBoundingBoxFlag:
    ioRepresentationBlock.mMeshRenderControl = 'Outline'
    if showSurfacesFlag | showEdgesFlag | showPointsFlag:
      if PhactoriDbg():
        myDebugPrint3("  warning:  when show bounding box is true, \n" + \
            "  show surfaces, show edges, and show points should be false\n")
  elif showPointsFlag:
    ioRepresentationBlock.mMeshRenderControl = 'Points'
  else:
    if showSurfacesFlag:
      if showEdgesFlag:
        ioRepresentationBlock.mMeshRenderControl = 'Surface With Edges'
      else:
        ioRepresentationBlock.mMeshRenderControl = 'Surface'
    else:
      if showEdgesFlag:
        ioRepresentationBlock.mMeshRenderControl = 'Wireframe'
        #ioRepresentationBlock.mMeshRenderControl = 'Points'
      else:
        ioRepresentationBlock.mMeshRenderControl = 'Outline'

  if doVolumeRenderingFlag:
    if PhactoriDbg(100):
      myDebugPrint3('doing volume rendering\n', 100)
    ioRepresentationBlock.mMeshRenderControl = 'Volume'
    ioRepresentationBlock.mDoingVolumeRendering = True
  else:
    ioRepresentationBlock.mDoingVolumeRendering = False

  #color legend on/off
  ioRepresentationBlock.mColorLegendFlag = getParameterFromBlock(inJsn,
    'show color legend', ioRepresentationBlock.mColorLegendFlag)

  ioRepresentationBlock.mColorLegendPositionAndSize = \
    getParameterFromBlock(inJsn, 'color legend position',
        ioRepresentationBlock.mColorLegendPositionAndSize)

  if ioRepresentationBlock.mColorLegendFlag == True:
    if ioRepresentationBlock.mColorVariableInfo.mVariableName == '':
      ioRepresentationBlock.mColorLegendFlag = False

  ioRepresentationBlock.mColorRangeMinMaxTracker.PlotValMinMaxTrkCParseJson(
          inJsn, "color legend ")

  ioRepresentationBlock.mTimeAnnotationSettings.ParseAvsFromJson(inJsn)

  ioRepresentationBlock.mDataCubeAxesFlag = getParameterFromBlock(inJsn,
    'show axes', ioRepresentationBlock.mDataCubeAxesFlag)

  ioRepresentationBlock.mDataCubeAxesInfo.DcaiParseParametersFromJson(inJsn)

  ioRepresentationBlock.mOrientationAxesFlag = getParameterFromBlock(inJsn,
    'show orientation axes', ioRepresentationBlock.mOrientationAxesFlag)

  if "image name addon" in inJsn:
    ioRepresentationBlock.mFilenameAddon = inJsn["image name addon"]

  ioRepresentationBlock.mColorSettings.ParseColorSettingsFromJson(inJsn)

  if PhactoriDbg(100):
    myDebugPrint3('ParseOneRepresentationBlockC returning\n', 100)

global gDefaultImageBasename
gDefaultImageBasename = "csierra.view."

def SetDefaultImageBasename(inNewImageBasename):
  global gDefaultImageBasename
  gDefaultImageBasename = inNewImageBasename

def breakSpecialVarNameIntoBaseAndComponent(inSpecialVarName,
      inAddSeparatorFlag):

    varNameLen = len(inSpecialVarName)
    if varNameLen < 3:
      errStr = '  in breakSpecialVarNameIntoBaseAndComponent representation bad color by vector component name (too short)\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)
    lastVarChar = inSpecialVarName[varNameLen - 1]
    if PhactoriDbg():
      myDebugPrint3('  lastVarChar: ' + str(lastVarChar) + '\n')
    if lastVarChar == 'X' or lastVarChar == 'x':
      component = 0
    elif lastVarChar == 'Y' or lastVarChar == 'y':
      component = 1
    elif lastVarChar == 'Z' or lastVarChar == 'z':
      component = 2
    else:
      errStr = '  in breakSpecialVarNameIntoBaseAndComponent representation bad color by vector component name (does not end in X Y or Z)\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)
    #baseVarName = inSpecialVarName[0:(varNameLen-2)] + GetSeparatorString()
    lenToGrab = varNameLen - 1 - len(GetSeparatorString())
    varNameWithoutSeparator = inSpecialVarName[0:lenToGrab]
    if inAddSeparatorFlag:
      baseVarName = varNameWithoutSeparator + GetSeparatorString()
    else:
      baseVarName = varNameWithoutSeparator

    return [baseVarName, component, varNameWithoutSeparator]


class PhactoriImageSettings:
  def __init__(self):
    self.mImageBasename = None
    global gPipeAndViewsState
    self.mImageBasedirectory = gPipeAndViewsState.mDefaultBasedirectory
    global gDefaultImageSizeX
    global gDefaultImageSizeY
    self.mImageSize = [gDefaultImageSizeX, gDefaultImageSizeY]
    self.mImageFormat = 'png'
    global gDefaultNumCounterDigits
    self.mNumCounterDigits = gDefaultNumCounterDigits
    self.mUsingDefaultGeneratedImageBasename = False
    self.mPixelBorderRatioXY = gPixelBorderRatioXY
  def GetAspectRatio(self):
    return float(self.mImageSize[0]) / float(self.mImageSize[1])
  def GetAspectRatioInsidePixelBorder(self):
    vertPixBrdrX = self.mPixelBorderRatioXY[0] * float(self.mImageSize[0])
    vertPixBrdrY = self.mPixelBorderRatioXY[1] * float(self.mImageSize[1])
    vertPixBrdrX = math.ceil(vertPixBrdrX)
    vertPixBrdrY = math.ceil(vertPixBrdrY)
    insideVertPixX = float(self.mImageSize[0]) - 2.0 * vertPixBrdrX
    insideVertPixY = float(self.mImageSize[1]) - 2.0 * vertPixBrdrY
    insideAspectRatio = insideVertPixX / insideVertPixY
    return insideAspectRatio

  def ParseImageSettingsInfo(self, inImageSettingsJson,
        inImageBasenameKey, inImageBaseDirectoryKey):
    self.mImageBasename = getParameterFromBlock(inImageSettingsJson,
      inImageBasenameKey, self.mImageBasename)
    if self.mImageBasename == None:
      global gDefaultImageBasename
      global gPipeAndViewsState
      blockCountIdStr = GetCurrentOutputResultsBlockCountId()
      imageSetCount = gPipeAndViewsState.mImageSetCounter
      gPipeAndViewsState.mImageSetCounter += 1
      self.mImageBasename = gDefaultImageBasename + \
          ".b-" + blockCountIdStr + "-is-" + str(imageSetCount) + "."
      self.mUsingDefaultGeneratedImageBasename = True

    self.mImageBasedirectory = getParameterFromBlock(inImageSettingsJson,
      inImageBaseDirectoryKey, self.mImageBasedirectory)

    if self.mImageBasedirectory != '':
      #test for directory existence, try to create if necessary
      if SmartGetLocalProcessId() == 0:
        import os
        if PhactoriDbg():
          myDebugPrint3('find or create self.mImageBasedirectory\n')
        if os.path.exists(self.mImageBasedirectory):
          if os.path.isdir(self.mImageBasedirectory) != True:
            errStr = \
              '  path to base directory exists and is not a directory\n' + \
              ' (' + self.mImageBasedirectory + ')\n'
            if PhactoriDbg():
              myDebugPrint3(errStr)
            raise Exception(errStr)
          else:
            if PhactoriDbg():
              myDebugPrint3('self.mImageBasedirectory found and is directory\n')
        else:
          if PhactoriDbg():
            myDebugPrint3('self.mImageBasedirectory to be created\n')
          os.makedirs(self.mImageBasedirectory)
          if os.path.exists(self.mImageBasedirectory) != True:
            errStr = \
              '  failed to create directory for images\n' + \
              ' (' + self.mImageBasedirectory + ')\n'
            if PhactoriDbg():
              myDebugPrint3(errStr)
            raise Exception(errStr)

    global gDefaultImageSizeX
    global gDefaultImageSizeY
    self.mImageSize = getParameterFromBlock(inImageSettingsJson,
      'image size', self.mImageSize)

    self.mImageFormat = getParameterFromBlock(inImageSettingsJson,
      'image format', self.mImageFormat)

    self.mNumCounterDigits = getParameterFromBlock(\
      inImageSettingsJson,
      'image digit count', self.mNumCounterDigits)
    if PhactoriDbg():
      myDebugPrint3('  mNumCounterDigits is ' + str(self.mNumCounterDigits) + '\n')


def ParseOneImagesetBlockC(ioImagesetBlock, ioImagesetBlockJson,
        ioPipeAndViewsState):

    #hack--if we have an imageset with a particular name, we enable
    #interaction
    if ioImagesetBlock.mName == 'PhactoriGoInteractive':
      ioPipeAndViewsState.mInteractionEnabled = True

    #parse on/off criteria
    if 'onoff criteria' in ioImagesetBlockJson:
      inCriteriaList = ioImagesetBlockJson['onoff criteria']
      for oneCriteriaName in inCriteriaList:
        if oneCriteriaName not in \
                ioPipeAndViewsState.mImagesetOnOffCriteriaBlocks:
          myDebugPrint3AndException("ParseOneImagesetBlockC:"
            "imageset with name: " + ioImagesetBlock.mName + "\n"
            "calls for nonexistent onoff criteria with name: " + \
            oneCriteriaName + "\n")
        ioImagesetBlock.mImagesetOnOffFilter.AddStartCriteria( \
          ioPipeAndViewsState.mImagesetOnOffCriteriaBlocks[oneCriteriaName])

    if 'camera' not in ioImagesetBlockJson:
      #we have to construct and use a default camera, including parsing
      #commands in the imageset for the camera
      if PhactoriDbg():
        myDebugPrint3("  ParseOneImagesetBlockC: for imageset " + \
            ioImagesetBlock.mName + \
            " there is no camera, so we must add and reference default\n")
      if 'camera type' not in ioImagesetBlockJson:
        ioImagesetBlockJson['camera type'] = 'multicamera8'
      defaultCameraName = ioImagesetBlock.mName + '_default_camera'
      defaultCameraBlockAndWrapper = {defaultCameraName: ioImagesetBlockJson}
      ParseBlocksC(ioPipeAndViewsState.mCameraBlocks,
          defaultCameraBlockAndWrapper,
          PhactoriCameraBlock,
          ParseOneCameraBlockC,
          ioPipeAndViewsState)
      ioImagesetBlockJson['camera'] = defaultCameraName
      #myDebugPrint3(  "done adding camera, here it is:\n")
      #ioPipeAndViewsState.mCameraBlocks[defaultCameraName].PrintSelf()


    #there is a camera (if there wasn't, we added it), so get the reference
    cameraName = ioImagesetBlockJson['camera']
    if cameraName not in ioPipeAndViewsState.mCameraBlocks:
      errStr = '  in ParseOneImagesetBlockC imageset (' + str(ioImagesetBlock.mName) + ') calls for camera (' + str(cameraName) + ') which does not exist\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    ioImagesetBlock.mCamera = ioPipeAndViewsState.mCameraBlocks[cameraName]
    if PhactoriDbg():
      myDebugPrint3("  image set " + ioImagesetBlock.mName + " is using the following camera:\n")
    ioImagesetBlock.mCamera.PrintSelf()

    if 'markers' in ioImagesetBlockJson:
      ioImagesetBlock.mVisibleMarkerNames = ioImagesetBlockJson['markers']
      for oneMarkerName in ioImagesetBlock.mVisibleMarkerNames:
        if oneMarkerName in ioPipeAndViewsState.mMarkerBlocks:
            ioImagesetBlock.mVisibleMarkers.append(
                    ioPipeAndViewsState.mMarkerBlocks[oneMarkerName])
        else:
          myDebugPrint3AndException(
            "ParseOneImagesetBlockC:\n" \
            "image set with name: " + ioImagesetBlock.mName + "\n"
            "calls for nonexistent marker named: " + oneMarkerName + "\n")

    if 'text annotations' in ioImagesetBlockJson:
      ioImagesetBlock.mTextAnnotationNames = \
              ioImagesetBlockJson['text annotations']
      for oneMarkerName in ioImagesetBlock.mTextAnnotationNames:
        if oneMarkerName in ioPipeAndViewsState.mTextAnnotationBlocks:
            ioImagesetBlock.mTextAnnotations.append(
                    ioPipeAndViewsState.mTextAnnotationBlocks[oneMarkerName])
        else:
          myDebugPrint3AndException(
            "ParseOneImagesetBlockC:\n" \
            "image set with name: " + ioImagesetBlock.mName + "\n"
            "calls for nonexistent marker named: " + oneMarkerName + "\n")

    #if 'representation' not in ioImagesetBlockJson:
    #  #we have to construct and use a default representation, including parsing
    #  #commands in the imageset for the representation
    #  if PhactoriDbg():
    #    myDebugPrint3("  ParseOneImagesetBlockC: for imageset " + \
    #        ioImagesetBlock.mName + \
    #        " there is no representation, " +
    #        "so we must add and reference default\n")
    #  defaultRepName = ioImagesetBlock.mName + '_default_representation'
    #  defaultRepBlockAndWrapper = {defaultRepName: ioImagesetBlockJson}
    #  ParseBlocksC(ioPipeAndViewsState.mRepresentationBlocks,
    #      defaultRepBlockAndWrapper,
    #      PhactoriRepresentationBlock,
    #      ParseOneRepresentationBlockC,
    #      ioPipeAndViewsState)
    #  ioImagesetBlockJson['representation'] = defaultRepName

    #if 'representation' in ioImagesetBlockJson:
    #  representationName = ioImagesetBlockJson['representation']
    #  if representationName not in ioPipeAndViewsState.mRepresentationBlocks:
    #    errStr = '  in ParseOneImagesetBlockC imageset (' + str(ioImagesetBlock.mName) + ') calls for representation (' + str(representationName) + ') which does not exist\n'
    #    if PhactoriDbg():
    #      myDebugPrint3(errStr)
    #    raise Exception(errStr)
    #  ioImagesetBlock.mRepresentation = ioPipeAndViewsState.mRepresentationBlocks[representationName]
    #else:
    #  ioImagesetBlock.mRepresentation = None #need to use default representation
    #ioImagesetBlock.mVisibleReps.append(ioImagesetBlock.mRepresentation)

    ##get the operation referred to by this image set (or the default incoming
    ##data)
    #ioImagesetBlock.mOperation =\
    #        ioPipeAndViewsState.GetOperationReferredByJson(
    #                'operation', ioImagesetBlockJson)
    #ioImagesetBlock.mVisibleOps.append(ioImagesetBlock.mOperation)

    ioImagesetBlock.ParseOperationAndRepresentationPair(ioPipeAndViewsState,
      ioImagesetBlockJson, 'operation', False, 'representation', False, False)
    ioImagesetBlock.mOperation = ioImagesetBlock.mVisibleOps[0]
    ioImagesetBlock.mRepresentation = ioImagesetBlock.mVisibleReps[0]

    #initial implementation of having more than one operation visible:
    #have up to 3 visible with particular parsing names
    ioImagesetBlock.ParseOperationAndRepresentationPair(ioPipeAndViewsState,
      ioImagesetBlockJson, 'operation2', True, 'representation2', True, True)
    ioImagesetBlock.ParseOperationAndRepresentationPair(ioPipeAndViewsState,
      ioImagesetBlockJson, 'operation3', True, 'representation3', True, True)

    ioImagesetBlock.mImageSettings.ParseImageSettingsInfo(
        ioImagesetBlockJson, 'image basename', 'image basedirectory')
    ioImagesetBlock.mImageFileNameCountSettings.\
        ParseImageFileNameCountSettings(ioImagesetBlockJson)


def SetUpPlotAxisNameDetails(inVariableInfo, ioDataCubeAxisInfo):
  if inVariableInfo.mVariableIsVectorComponent:
    if PhactoriDbg():
      myDebugPrint3("  vector is component, setting y axis name\n")
    varNameForPlot = inVariableInfo.mVectorBaseName
    if inVariableInfo.mVariableComponent == 0:
      varNameForPlot = varNameForPlot + '[x]'
    elif inVariableInfo.mVariableComponent == 1:
      varNameForPlot = varNameForPlot + '[y]'
    elif inVariableInfo.mVariableComponent == 2:
      varNameForPlot = varNameForPlot + '[z]'
    else:
      varNameForPlot = varNameForPlot + str(inVariableInfo.mVariableComponent)
    ioDataCubeAxisInfo.mUseLabelFlag = True
    ioDataCubeAxisInfo.mAxisLabel = varNameForPlot
    ioDataCubeAxisInfo.mBaseAxisLabel = varNameForPlot
    if PhactoriDbg():
      myDebugPrint3("  y axis name is " + ioDataCubeAxisInfo.mAxisLabel + "\n")
  elif inVariableInfo.mVariableIsVectorMagnitude:
    if PhactoriDbg():
      myDebugPrint3("  vector is magnitude, setting y axis name\n")
    varNameForPlot = inVariableInfo.mVectorBaseName
    varNameForPlot = varNameForPlot + " Magnitude"
    ioDataCubeAxisInfo.mUseLabelFlag = True
    ioDataCubeAxisInfo.mAxisLabel = varNameForPlot
    ioDataCubeAxisInfo.mBaseAxisLabel = varNameForPlot
    if PhactoriDbg():
      myDebugPrint3("  y axis name is " + ioDataCubeAxisInfo.mAxisLabel + "\n")
  else:
    ioDataCubeAxisInfo.mUseLabelFlag = True
    ioDataCubeAxisInfo.mAxisLabel = inVariableInfo.mVariableName
    ioDataCubeAxisInfo.mBaseAxisLabel = inVariableInfo.mVariableName

def ParseOnePlotOverTimeBlockC(ioPlotOverTimeBlock, inPlotOverTimeBlockJson, inPipeAndViewsState):
  if PhactoriDbg(100):
    myDebugPrint3('ParseOnePlotOverTimeBlockC entered\n', 100)

  #hack test for colors start
  #inPlotOverTimeBlockJson['background color'] = [0.0, 0.0, 0.0]
  ##inJsn['axes color'] = [0.0, 1.0, 1.0]
  #inPlotOverTimeBlockJson['text color'] = [1.0, 1.0, 0.0]
  #inPlotOverTimeBlockJson['line color'] = [1.0, 1.0, 1.0]
  #hack test for colors end

  ioPlotOverTimeBlock.mImageSettings.ParseImageSettingsInfo(
      inPlotOverTimeBlockJson, 'plot basename', 'plot basedirectory')
  ioPlotOverTimeBlock.mImageFileNameCountSettings.\
      ParseImageFileNameCountSettings(inPlotOverTimeBlockJson)

  if PhactoriDbg():
    myDebugPrint3(str(inPlotOverTimeBlockJson) + '\n')

  #parse y axis variable
  ioPlotOverTimeBlock.m_YAxisVariableInfo.\
    ParseVariableNameAndVectorOrTensorComponent(inPlotOverTimeBlockJson,
        'variable ')

  #handle name of axis, particularly in vector component/magnitude case
  SetUpPlotAxisNameDetails(ioPlotOverTimeBlock.m_YAxisVariableInfo,
      ioPlotOverTimeBlock.m_DataCubeAxesInfo.mYAxisInfo)

  #parse whether to draw min, max mean
  if 'plot minimum' in inPlotOverTimeBlockJson:
    ioPlotOverTimeBlock.mPlotMinimumFlag = inPlotOverTimeBlockJson['plot minimum']
  else:
    ioPlotOverTimeBlock.mPlotMinimumFlag = True
  if 'plot maximum' in inPlotOverTimeBlockJson:
    ioPlotOverTimeBlock.mPlotMaximumFlag = inPlotOverTimeBlockJson['plot maximum']
  else:
    ioPlotOverTimeBlock.mPlotMaximumFlag = True
  if 'plot mean' in inPlotOverTimeBlockJson:
    ioPlotOverTimeBlock.mPlotMeanFlag = inPlotOverTimeBlockJson['plot mean']
  else:
    ioPlotOverTimeBlock.mPlotMeanFlag = True

  if 'plot id' in inPlotOverTimeBlockJson:
    idList = inPlotOverTimeBlockJson['plot id']
    for ii in idList:
      newIdPlotLine = PhactoriPlotOverTimeIdLine(ii)
      ioPlotOverTimeBlock.m_IdPlotLineList.append(newIdPlotLine)

  ioPlotOverTimeBlock.m_xyzMinMaxTrkC.PlotXYZMinMaxTrkCParseJson(inPlotOverTimeBlockJson, 'axis ')

  #hack to test missing data situations
  #if ioPlotOverTimeBlock.mName == "fooPlot":
  #  myDebugPrint3("hack to test missing data, found fooPlot, using stressthresh\n")
  #  ioPlotOverTimeBlock.mInputOperation = inPipeAndViewsState.mOperationBlocks["stressthresh"]
  #else:
  #  myDebugPrint3("hack to test missing data, did not find fooPlot\n")
  #  ioPlotOverTimeBlock.mInputOperation = inPipeAndViewsState.GetOperationReferredByJson('operation', inJsn)

  ioPlotOverTimeBlock.mInputOperation = \
      inPipeAndViewsState.GetOperationReferredByJson('operation',
          inPlotOverTimeBlockJson)

  ioPlotOverTimeBlock.mColorSettings.\
       ParsePlotColorSettingsFromJson(inPlotOverTimeBlockJson)

  if PhactoriDbg(100):
    myDebugPrint3('ParseOnePlotOverTimeBlockC returning:\n', 100)
  return

def ParseBlocksC2(ioBlockSet, inBlocksJson, inBlockClass, ioPipeAndViewsState):
  if PhactoriDbg(100):
    myDebugPrint3('ParseBlocksC2 entered\n', 100)
  count = 0
  for blockName, blockParamsJson in inBlocksJson.iteritems():
    count = count + 1
    newBlockInstance = inBlockClass()
    newBlockInstance.mName = blockName
    if PhactoriDbg():
      myDebugPrint3('  creating block named ' + newBlockInstance.mName + '\n')
    if PhactoriDbg():
      myDebugPrint3(str(blockParamsJson) + '\n')
    newBlockInstance.ParseParametersFromJson(blockParamsJson, ioPipeAndViewsState)
    ioBlockSet[newBlockInstance.mName] = newBlockInstance
  if PhactoriDbg(100):
    myDebugPrint3('ParseBlocksC2 returning\n', 100)
  return count

def ParseBlocksC(ioBlockSet, inBlocksJson, inBlockClass, inParseOneBlockMethod,
    ioPipeAndViewsState):
  if PhactoriDbg(100):
    myDebugPrint3('ParseBlocksC entered\n', 100)
  count = 0
  for blockName, blockParamsJson in inBlocksJson.iteritems():
    count = count + 1
    newBlockInstance = inBlockClass()
    newBlockInstance.mName = blockName
    if PhactoriDbg():
      myDebugPrint3('  creating block named ' + newBlockInstance.mName + '\n')
    if PhactoriDbg():
      myDebugPrint3(str(blockParamsJson) + '\n')
    inParseOneBlockMethod(newBlockInstance, blockParamsJson, ioPipeAndViewsState)
    ioBlockSet[newBlockInstance.mName] = newBlockInstance
  if PhactoriDbg(100):
    myDebugPrint3('ParseBlocksC returning\n', 100)
  return count

class PhactoriUserPointInfo:
  def __init__(self):
    #'absolute' or 'datasize relative' or 'element' or 'node'
    self.mPointType = "datasize relative"
    self.mElementOrNodeId = None
    self.mXyz = [0.0, 0.0, 0.0]
    self.mMayChangeWithData = False
    #flag tracking whether any information was filled in due to parsing
    #or if default is simply being used because of no parsed items
    self.mParseHadKeys = None
    self.mHasDisplacedInfoFlag = False
    self.mDisplacement = [0.0, 0.0, 0.0]
    self.mParsingDetectedAtLeastOneSetting = False
    self.mReferenceVariableInfo = None

    #for tracking which source input to use if point is at min or max
    #variable
    self.mInputNameToUseForMinMaxVariable = None
    self.mInputPvSourceToUseForMinMaxVariable = None

  def UserPointHasValidDisplacementInfo():
    return self.mHasDisplacedInfoFlag

  def ParseIdAndDisplacement(self, inDisplacementFlag, inJson, inKey):
    if inDisplacementFlag:
      params = inJson[inKey]
      self.mElementOrNodeId = params[0]
      self.mHasDisplacedInfoFlag = True
      self.mDisplacement[0] = params[1]
      self.mDisplacement[1] = params[2]
      self.mDisplacement[2] = params[3]
      if PhactoriDbg():
        myDebugPrint3(" displacement is " + str(self.mDisplacement) + "\n")
    else:
      self.mElementOrNodeId = inJson[inKey]

  def UserPointInfoParseParametersFromJson(self, inJson, inKeyPrefix, inKeySuffix,
          inDisplacementFlag = False):
    self.mParseHadKeys = False
    keyAbsolute = inKeyPrefix + "absolute point" + inKeySuffix
    keyRelative = inKeyPrefix + "relative point" + inKeySuffix
    keyNode = inKeyPrefix + "node" + inKeySuffix
    keyElement = inKeyPrefix + "element" + inKeySuffix
    keyMaxVar = inKeyPrefix + "max variable point" + inKeySuffix
    keyMinVar = inKeyPrefix + "min variable point" + inKeySuffix
    keyDataPt = inKeyPrefix + "data point" + inKeySuffix

    if keyAbsolute in inJson:
      self.mParsingDetectedAtLeastOneSetting = True
      self.mPointType = "absolute"
      self.mXyz = inJson[keyAbsolute]
      self.mMayChangeWithData = False
    elif keyRelative in inJson:
      self.mParsingDetectedAtLeastOneSetting = True
      self.mPointType = "datasize relative"
      self.mXyz = inJson[keyRelative]
      self.mMayChangeWithData = True
    elif keyNode in inJson:
      self.mParsingDetectedAtLeastOneSetting = True
      if PhactoriDbg():
        myDebugPrint3("PhactoriUserPointInfo.ParseParametersFromJson node found\n")
      self.mPointType = "node"
      self.ParseIdAndDisplacement(inDisplacementFlag, inJson, keyNode)
      self.mMayChangeWithData = True
      if PhactoriDbg():
        myDebugPrint3(" node id is " + str(self.mElementOrNodeId) + "\n")
    elif keyElement in inJson:
      self.mParsingDetectedAtLeastOneSetting = True
      if PhactoriDbg():
        myDebugPrint3("PhactoriUserPointInfo.ParseParametersFromJson element found\n")
      self.mPointType = "element"
      self.ParseIdAndDisplacement(inDisplacementFlag, inJson, keyElement)
      self.mMayChangeWithData = True
      if PhactoriDbg():
        myDebugPrint3(" element id is " + str(self.mElementOrNodeId) + "\n")
    elif keyMinVar in inJson:
      self.mParsingDetectedAtLeastOneSetting = True
      self.mPointType = "min variable"
      #self.ParseIdAndDisplacement(inDisplacementFlag, inJson, keyMinNode):
      self.mReferenceVariableInfo = PhactoriVariableInfo()
      self.mReferenceVariableInfo.\
          ParseVariableNameAndVectorOrTensorComponent(inJson, 'variable ')
      self.mMayChangeWithData = True
      if PhactoriDbg():
        myDebugPrint3("PhactoriUserPointInfo.ParseParametersFromJson "
            "min variable found\n"
            "var info: " + self.mReferenceVariableInfo.SelfToStr() + "\n")
      varSourceKey = inKeyPrefix + "input for min variable point" + inKeySuffix
      if varSourceKey in inJson:
        self.mInputNameToUseForMinMaxVariable = inJson[varSourceKey]
        myDebugPrint3("input name for min variable point: " + \
          self.mInputNameToUseForMinMaxVariable + "\n")
      varDispKey = \
        inKeyPrefix + "displacement for min variable point" + inKeySuffix
      if varDispKey in inJson:
        self.mDisplacement = inJson[varDispKey]
        myDebugPrint3("displacement for min variable point: " + \
          str(self.mDisplacement) + "\n")
    elif keyDataPt in inJson:
      self.mParsingDetectedAtLeastOneSetting = True
      if PhactoriDbg():
          myDebugPrint3("PhactoriUserPointInfo.ParseParametersFromJson: "
                  "parsing data point")
      params = inJson[keyDataPt]
      if inDisplacementFlag:
        #better have 7 items
        if len(params) != 7:
          myDebugPrint3AndException(
            "PhactoriUserPointInfo.ParseParametersFromJson:\n"
            "data point displaced expected 7 parameters")
      else:
        #better have 4 items
        if len(params) != 4:
          myDebugPrint3AndException(
            "PhactoriUserPointInfo.ParseParametersFromJson:\n"
            "data point expected 4 parameters")
      param0 = params[0].lower()
      if param0 == "min":
        self.mPointType = "min variable"
      if param0 == "max":
        self.mPointType = "max variable"
      if param0 == "center":
        myDebugPrint3AndException(
          "PhactoriUserPointInfo.ParseParametersFromJson:\n"
          "data point 'center' not yet implemented")
      self.mReferenceVariableInfo = PhactoriVariableInfo()
      tempJson = {}
      tempJson[params[1].lower()] = params[2]
      if PhactoriDbg():
        print "tempJson: " + str(tempJson)
      self.mReferenceVariableInfo.\
          ParseVariableNameAndVectorOrTensorComponent(tempJson, '')
      self.mMayChangeWithData = True
      if PhactoriDbg():
        myDebugPrint3("variable found\n"
            "var info: " + self.mReferenceVariableInfo.SelfToStr() + "\n")
      if params[3] != "default":
        self.mInputNameToUseForMinMaxVariable = params[3]
        if PhactoriDbg():
          myDebugPrint3("name of input operation: " +
              self.mInputNameToUseForMinMaxVariable + "\n")
      if inDisplacementFlag:
        self.mHasDisplacedInfoFlag = True
        self.mDisplacement[0] = params[4]
        self.mDisplacement[1] = params[5]
        self.mDisplacement[2] = params[6]

      else:
        if PhactoriDbg():
          myDebugPrint3("using default input operation\n")
      if PhactoriDbg():
          myDebugPrint3("PhactoriUserPointInfo.ParseParametersFromJson: "
                  "done parsing data point")
    elif keyMaxVar in inJson:
      self.mParsingDetectedAtLeastOneSetting = True
      self.mPointType = "max variable"
      #self.ParseIdAndDisplacement(inDisplacementFlag, inJson, keyMinNode):
      self.mReferenceVariableInfo = PhactoriVariableInfo()
      self.mReferenceVariableInfo.\
          ParseVariableNameAndVectorOrTensorComponent(inJson, 'variable ')
      self.mMayChangeWithData = True
      if PhactoriDbg():
        myDebugPrint3("PhactoriUserPointInfo.ParseParametersFromJson "
            "max variable found\n"
            "var info: " + self.mReferenceVariableInfo.SelfToStr() + "\n")
      varSourceKey = inKeyPrefix + "input for max variable point" + inKeySuffix
      if varSourceKey in inJson:
        self.mInputNameToUseForMinMaxVariable = inJson[varSourceKey]
        myDebugPrint3("input name for max variable point: " + \
          self.mInputNameToUseForMinMaxVariable + "\n")
      varDispKey = \
        inKeyPrefix + "displacement for max variable point" + inKeySuffix
      if varDispKey in inJson:
        self.mDisplacement = inJson[varDispKey]
        myDebugPrint3("displacement for max variable point: " + \
          str(self.mDisplacement) + "\n")
    else:
      self.mParsingDetectedAtLeastOneSetting = False
      if PhactoriDbg():
        myDebugPrint3(" no PhactoriUserPointInfo parse keys, using default\n")

    #myDebugPrint3("PhactoriUserPoint Parsed from json:\n")
    #self.PrintSelf()

  def GetCurrentGeometricPoint(self, inParaViewSource, ioViewBounds,
      inUpdateInputPipelineFlag):
    """returns the current absolute xyz of this point, calculating relative
       position or element or node position if necessary; if it obtains
       its view bounds during operation it returns them, otherwise it
       returns None
    """
    if inUpdateInputPipelineFlag:
      UpdatePipelineWithCurrentTimeArgument(inParaViewSource)
    if PhactoriDbg(100):
      myDebugPrint3("GetCurrentGeometricPoint entered\n", 100)

    #import pdb
    #pdb.set_trace()

    #self.PrintSelf()
    if self.mPointType == 'absolute':
      #focal point is exact point in coordinate system of data
      returnXyz = list(self.mXyz)
    elif self.mPointType == 'datasize relative':
      #focal point is offset from center of data bounding box, with offset
      #relative in size to the size of the overall bounding box
      relativeXyz = self.mXyz
      myViewBounds = ioViewBounds[0]
      if myViewBounds == None:
        myViewBounds = GetGlobalDataBoundsParallel(inParaViewSource)
        ioViewBounds[0] = myViewBounds
      #old style:  use maximum bounding box dimension to define
      #relative space
      #maxDim = GetMaximumDimensionFromBounds(myViewBounds)
      #returnXyz = [ (myViewBounds[1] + myViewBounds[0]) * 0.5 + \
      #                     maxDim * relativeXyz[0], \
      #                 (myViewBounds[3] + myViewBounds[2]) * 0.5 + \
      #                     maxDim * relativeXyz[1], \
      #                 (myViewBounds[5] + myViewBounds[4]) * 0.5 + \
      #                     maxDim * relativeXyz[2] ]

      #new style:  do relative in each dimension
      returnXyz = [0.0,0.0,0.0]
      returnXyz[0] = myViewBounds[0] + (myViewBounds[1] - myViewBounds[0]) * \
          (relativeXyz[0] + 0.5)
      returnXyz[1] = myViewBounds[2] + (myViewBounds[3] - myViewBounds[2]) * \
          (relativeXyz[1] + 0.5)
      returnXyz[2] = myViewBounds[4] + (myViewBounds[5] - myViewBounds[4]) * \
          (relativeXyz[2] + 0.5)
    elif self.mPointType == 'node':
      nodeId = self.mElementOrNodeId
      returnXyz = [0.0, 0.0, 0.0]
      GetXyzForNodeOrElementParallel(inParaViewSource, True, nodeId, returnXyz)
    elif self.mPointType == 'element':
      elementId = self.mElementOrNodeId
      returnXyz = [0.0, 0.0, 0.0]
      GetXyzForNodeOrElementParallel(inParaViewSource, False, elementId, returnXyz)
    elif self.mPointType == 'min variable':
      returnXyz = [0.0, 0.0, 0.0]
      if self.mInputNameToUseForMinMaxVariable == None:
        pvFilter = inParaViewSource
      else:
        if self.mInputPvSourceToUseForMinMaxVariable == None:
          global gPipeAndViewsState
          operationBlock = gPipeAndViewsState.mOperationBlocks[
              self.mInputNameToUseForMinMaxVariable]
          pvFilter = operationBlock.GetPvFilter()
          self.mInputPvSourceToUseForMinMaxVariable = pvFilter
        else:
          pvFilter = self.mInputPvSourceToUseForMinMaxVariable
      if(pvFilter == None):
        #no pvFilter for min/max yet, just use 0,0,0
        returnXyz = [0.0, 0.0, 0.0]
      else:
        UpdatePipelineWithCurrentTimeArgument(pvFilter)
        GetXyzForMinOrMaxVariable(pvFilter, True,
            self.mReferenceVariableInfo, returnXyz)
    elif self.mPointType == 'max variable':
      returnXyz = [0.0, 0.0, 0.0]
      if self.mInputNameToUseForMinMaxVariable == None:
        pvFilter = inParaViewSource
        if PhactoriDbg():
          myDebugPrint3("max var point has no source, using incoming\n")
      else:
        if PhactoriDbg():
          myDebugPrint3("max var point has source, using instead\n")
        if self.mInputPvSourceToUseForMinMaxVariable == None:
          if PhactoriDbg():
            myDebugPrint3("not previously found, must locate\n")
          #global gPipeAndViewsState
          operationBlock = gPipeAndViewsState.mOperationBlocks[
              self.mInputNameToUseForMinMaxVariable]
          pvFilter = operationBlock.GetPvFilter()
          self.mInputPvSourceToUseForMinMaxVariable = pvFilter
        else:
          if PhactoriDbg():
            myDebugPrint3("previously found.\n")
          pvFilter = self.mInputPvSourceToUseForMinMaxVariable
        if PhactoriDbg():
          myDebugPrint3("incoming: " + str(inParaViewSource) + "\n"
            "using: " + str(pvFilter) + "\n")
      if(pvFilter == None):
        #no pvFilter for min/max yet, just use 0,0,0
        returnXyz = [0.0, 0.0, 0.0]
      else:
        UpdatePipelineWithCurrentTimeArgument(pvFilter)
        GetXyzForMinOrMaxVariable(pvFilter, False,
            self.mReferenceVariableInfo, returnXyz)
    else:
      errStr = 'GetCurrentGeometricPoint error! bad mPointType\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    if PhactoriDbg(100):
      myDebugPrint3("GetCurrentGeometricPoint returning: " + \
          str(returnXyz) + "\n", 100)

    return returnXyz

  def GetCurrentGeometricPointWithDisplacement(self, inParaViewSource,
      ioViewBounds, inUpdateInputPipelineFlag):
    undisplacedXyz = self.GetCurrentGeometricPoint(inParaViewSource,
        ioViewBounds, inUpdateInputPipelineFlag)
    returnXyz = [undisplacedXyz[0] + self.mDisplacement[0],
                 undisplacedXyz[1] + self.mDisplacement[1],
                 undisplacedXyz[2] + self.mDisplacement[2]]
    return returnXyz

  def SelfToStr(self):
      return "  PhactoriUserPointInfo:\n" \
        "    mPointType: " + self.mPointType + \
        "    mElementOrNodeId: " + str(self.mElementOrNodeId) + "\n" + \
        "    mMayChangeWithData: " + str(self.mMayChangeWithData) + "\n" \
        "    mParseHadKeys: " + str(self.mParseHadKeys) + "\n" \
        "    mHasDisplacedInfoFlag: " + str(self.mHasDisplacedInfoFlag) + "\n" \
        "    mDisplacement: " + str(self.mDisplacement) + "\n" \
        "    mReferenceVariableInfo: " + str(self.mReferenceVariableInfo) + "\n"

  def PrintSelf(self):
    if PhactoriDbg():
      myDebugPrint3(self.SelfToStr())

class PhactoriOffAxisProjectionInfo:
  def __init__(self):
    self.mUseOffAxisProjection = False
    self.mPhysicalEyeAndScreenSettings = PhactoriPhysicalEyeAndScreenSetup()
    self.mLeftEyeFlag = True
    self.mIpdInModelUnits = 6.2
    self.mVirtualSelfSizeMultiplier = 1.0
    self.mAutoSize1Enabled = False
    self.mAutoSize1LockEyePosition = False
    self.mAutoSize1LockCountdown = 1
    self.mAutoSize1DistanceRatio = 1.0
    self.mAutoSize1LockedEyePosition = None
    self.mAutoSize1LockedFocalPoint = None
    self.mAutoSize1ViewAngleDeltaInRadians = 0.0

  def ParseParametersFromJson(self, inJsn):
    if 'use off axis projection' in inJsn:
      self.mUseOffAxisProjection = inJsn['use off axis projection']
      if PhactoriDbg():
        myDebugPrint3("PhactoriOffAxisProjectionInfo:ParseParametersFromJson\n"
            "mUseOffAxisProjection set to " + str(self.mUseOffAxisProjection))
    if self.mUseOffAxisProjection == True:
      self.mPhysicalEyeAndScreenSettings.ParseSettingsFromJson(inJsn)
      if 'which eye' in inJsn:
        if inJsn['which eye'] == 'left':
          self.mLeftEyeFlag = True
        else:
          self.mLeftEyeFlag = False
      else:
        if PhactoriDbg():
          myDebugPrint3("PhactoriOffAxisProjectionInfo:"
            "ParseParametersFromJson\n"
            "warning: 'which eye' not in json, defaulting to left\n")
      if 'ipd in model units' not in inJsn:
        myDebugPrint3AndException(
            "PhactoriOffAxisProjectionInfo::ParseParametersFromJson\n"
            "Error:  incoming json must have 'ipd in model units'\n")
      self.mIpdInModelUnits = inJsn['ipd in model units']
      if 'virtual self size multiplier' in inJsn:
        self.mVirtualSelfSizeMultiplier = inJsn['virtual self size multiplier']
      if 'auto size 1' in inJsn:
        self.mAutoSize1Enabled = inJsn['auto size 1']
        if 'auto size 1 lock nth eye position' in inJsn:
          self.mAutoSize1LockEyePosition = True
          self.mAutoSize1LockCountdown = inJsn['auto size 1 lock nth eye position']
        if 'auto size 1 distance ratio' in inJsn:
          self.mAutoSize1DistanceRatio = inJsn['auto size 1 distance ratio']
        if 'auto size 1 view angle delta in degrees' in inJsn:
          angleInDegrees = inJsn['auto size 1 view angle delta in degrees']
          self.mAutoSize1ViewAngleDeltaInRadians = math.radians(angleInDegrees)


  def SetUpIfEnabled(self, inRenderView, ioCameraPosition, ioLookAtPoint,
      inViewUpVector):
    theCamera = GetActiveCamera()
    if self.mUseOffAxisProjection == False:
      theCamera.SetUseOffAxisProjection(0)
      return

    projectionEyePosition = self.mPhysicalEyeAndScreenSettings.GetEyePosition(
      self.mLeftEyeFlag)
    scrnBttmLeft = self.mPhysicalEyeAndScreenSettings.GetScreenBottomLeft(
      self.mLeftEyeFlag)
    scrnBttmRight = self.mPhysicalEyeAndScreenSettings.GetScreenBottomRight(
      self.mLeftEyeFlag)
    scrnTopRight = self.mPhysicalEyeAndScreenSettings.GetScreenTopRight(
      self.mLeftEyeFlag)

    if PhactoriDbg():
      myDebugPrint3("SetUpIfEnabled setting camera parameters:"
          "\n  mProjectionEyePosition: " + str(projectionEyePosition) + \
          "\n  mScreenBottomLeft:      " + str(scrnBttmLeft) + \
          "\n  mScreenBottomRight:     " + str(scrnBttmRight) + \
          "\n  mScreenTopRight:        " + str(scrnTopRight) + \
          "\n")
    SetActiveView(inRenderView)
    theCamera.SetUseOffAxisProjection(1)
    #we are assuming camera is explicitly set up where you want the left or right eye
    #and so we don't use eye separation
    theCamera.SetEyeSeparation(0.0)
    theCamera.SetScreenBottomLeft(scrnBttmLeft)
    theCamera.SetScreenBottomRight(scrnBttmRight)
    theCamera.SetScreenTopRight(scrnTopRight)
    theCamera.SetEyePosition(projectionEyePosition)

    self.ChangeEyeAndFocusInModelSpaceForStereo(ioCameraPosition,
      ioLookAtPoint, inViewUpVector)

    if PhactoriDbg():
      myDebugPrint3("SetUpIfEnabled done setting camera parameters\n")

  def ChangeEyeAndFocusInModelSpaceForStereo(self,
      ioEyePosition, ioFocalPoint, inViewUpVector):

    #place eye and look at point so object is at screen in AutoSize1 case
    if self.mAutoSize1Enabled:
        if self.mAutoSize1LockEyePosition and self.mAutoSize1LockedEyePosition != None:
          ioEyePosition[0] = self.mAutoSize1LockedEyePosition[0]
          ioEyePosition[1] = self.mAutoSize1LockedEyePosition[1]
          ioEyePosition[2] = self.mAutoSize1LockedEyePosition[2]
          ioFocalPoint[0] = self.mAutoSize1LockedFocalPoint[0]
          ioFocalPoint[1] = self.mAutoSize1LockedFocalPoint[1]
          ioFocalPoint[2] = self.mAutoSize1LockedFocalPoint[2]
          if PhactoriDbg():
            myDebugPrint3("AutoSize1: reusing first eye position and focal point:\n" + \
                "  eye position: " + str(ioEyePosition) + "\n"
                "  focal point : " + str(ioFocalPoint) + "\n")
        else:
          fpToEyeVec = vecFromAToB(ioFocalPoint, ioEyePosition)
          fpToEyeVecNorm = vecNormalize(fpToEyeVec)
          dist = -self.mPhysicalEyeAndScreenSettings.mLeftEyeScreenBottomLeft[2]
          dist *= self.mAutoSize1DistanceRatio
          dist *= self.mVirtualSelfSizeMultiplier

          if PhactoriDbg():
            myDebugPrint3("AutoSize1:\n"
              "  before eye position: " +  str(ioEyePosition) + "\n"
              "  before focal point : " +  str(ioFocalPoint) + "\n")
          ioEyePosition[0] = ioFocalPoint[0] + fpToEyeVecNorm[0] * dist
          ioEyePosition[1] = ioFocalPoint[1] + fpToEyeVecNorm[1] * dist
          ioEyePosition[2] = ioFocalPoint[2] + fpToEyeVecNorm[2] * dist
          #angle look direction for AutoSize1 by moving focal point up or down
          #simple tangent work, but need to work it out to see it (we want
          #same distance to screen plane but looking more angled up or down)
          dy1 = -dist * math.tan(self.mAutoSize1ViewAngleDeltaInRadians)
          ioFocalPoint[0] += inViewUpVector[0] * dy1
          ioFocalPoint[1] += inViewUpVector[1] * dy1
          ioFocalPoint[2] += inViewUpVector[2] * dy1

          if PhactoriDbg():
            myDebugPrint3(
              "  calculated eye position: " + str(ioEyePosition) + "\n"
              "  calculated focal point : " + str(ioFocalPoint) + "\n")
          if self.mAutoSize1LockEyePosition:
            self.mAutoSize1LockCountdown -= 1
            if PhactoriDbg():
              myDebugPrint3("self.mAutoSize1LockCountdown now " + \
                  str(self.mAutoSize1LockCountdown) + "\n")
            if self.mAutoSize1LockCountdown == 0:
              self.mAutoSize1LockedEyePosition = list(ioEyePosition)
              self.mAutoSize1LockedFocalPoint = list(ioFocalPoint)
              if PhactoriDbg():
                myDebugPrint3("(eye position saved for locking for next use)\n")

    #perpendicular screen distance or put look at point on screen plane?
    #i think put look at point on screen plane.

    lookDirection = vecFromAToB(ioEyePosition, ioFocalPoint)
    axisBetweenEyes = vecCrossProduct(lookDirection, inViewUpVector)
    axisBetweenEyesNorm = vecNormalize(axisBetweenEyes)

    if self.mLeftEyeFlag:
      modelCsHalfIpd = -self.mIpdInModelUnits * \
          self.mVirtualSelfSizeMultiplier * 0.5
    else:
      modelCsHalfIpd = self.mIpdInModelUnits * \
          self.mVirtualSelfSizeMultiplier * 0.5
    dx = axisBetweenEyesNorm[0] * modelCsHalfIpd
    dy = axisBetweenEyesNorm[1] * modelCsHalfIpd
    dz = axisBetweenEyesNorm[2] * modelCsHalfIpd
    if PhactoriDbg():
      myDebugPrint3("ChangeEyeAndFocusInModelSpaceForStereo entered\n"
        "axisBetweenEyesNorm: " + str(axisBetweenEyesNorm) + "\n"
        "modelCsHalfIpd: " + str(modelCsHalfIpd) + "\n"
        "[dx,dy,dz]: " + str([dx,dy,dz]) + "\n"
        "prior ioEyePosition: " + str(ioEyePosition) + "\n"
        "prior ioFocalPoint: " + str(ioFocalPoint) + "\n")
    ioEyePosition[0] += dx
    ioEyePosition[1] += dy
    ioEyePosition[2] += dz
    ioFocalPoint[0] += dx
    ioFocalPoint[1] += dy
    ioFocalPoint[2] += dz
    if PhactoriDbg():
      myDebugPrint3("post ioEyePosition: " + str(ioEyePosition) + "\n"
        "post ioFocalPoint: " + str(ioFocalPoint) + "\n"
        "ChangeEyeAndFocusInModelSpaceForStereo returning\n")


class PhactoriCameraBlock:
  def __init__(self):
    self.mName = ""
    self.mType = ""
    self.mFilenameAddon = ""
    self.mLookAtPointInfo = PhactoriUserPointInfo()
    self.mUseCameraAtPointFlag = False
    self.mCameraAtPointInfo = PhactoriUserPointInfo()
    self.mUseParallelProjection = False
    #0 is absolute, 1, is relative
    self.mParallelScaleAbsoluteOrRelative = 1
    self.mParallelScale = 1.0

    #'absolute' or 'datasize relative'
    self.mLookAtDistanceType = ""
    self.mLookAtDistance = 1.0
    self.mLookDirection = [-1.0, -1.0, -1.0]
    self.mLookDirectionSpecifiedFlag = False

    self.mViewAngle = 30.0
    self.mViewUpVector = [0.0, 1.0, 0.0]

    self.mOffAxisProjectionInfo = PhactoriOffAxisProjectionInfo()

  def PrintSelf(self):
    if PhactoriDbg():
      myDebugPrint3("PhactoriCameraBlock: name " + self.mName + \
        "  type " + self.mType + "\n" +
        "  look at distance type " + self.mLookAtDistanceType + \
        "  setting " + str(self.mLookAtDistance) + "\n" +
        "  look direction " + str(self.mLookDirection) + "\n")
    self.mLookAtPointInfo.PrintSelf()

  def MayChangeWithData(self):
    if self.mLookAtDistanceType != 'absolute':
      return True
    if self.mLookAtPointInfo.mMayChangeWithData == True:
      return True
    return False

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
    self.mAddSeparatorToVectorVariableName = True

  def SelfToStr(self):
    retStr = "PhactoriVariableInfo:\n" +\
      "\n  mVariableName: " + str(self.mVariableName) +\
      "\n  mVariableComponent: " + str(self.mVariableComponent) +\
      "\n  mVariableType: " + str(self.mVariableType) +\
      "\n  mVariableIsVectorMagnitude: " + \
          str(self.mVariableIsVectorMagnitude) +\
      "\n  mVariableIsVectorComponent: " + \
          str(self.mVariableIsVectorComponent) +\
      "\n  mVectorBaseName: " + str(self.mVectorBaseName) +\
      "\n"
    return retStr

  def ParseVariableNameAndVectorOrTensorComponent(self, inJson, inBaseString):
    """take a base string such as 'y axis variable ' or 'variable ', use it
    to construct keys to define a scalar variable, vector component, vector
    magnitude, or tensor component, and see if the json has those keys.
    If the json has a key, use it to grab the variable name and setup.
    Also look to see if the type of the variable is specifically defined
    (node or element) or needs to be detected
    """
    variableFoundFlag = False
    self.mVariableIsVectorComponent = False
    self.mVariableIsVectorMagnitude = False

    if 'add separator to vector variable name' in inJson:
      self.mAddSeparatorToVectorVariableName = \
        inJson['add separator to vector variable name']

    testKey = inBaseString + 'scalar'
    if testKey in inJson:
      variableFoundFlag = True
      if PhactoriDbg():
        myDebugPrint3('  scalar found\n')
      self.mVariableName = inJson[testKey]
      self.mVariableComponent = None
    else:
      testKey = inBaseString + 'vector component'
      if testKey in inJson:
        variableFoundFlag = True
        self.mVariableIsVectorComponent = True
        if PhactoriDbg():
          myDebugPrint3('  vector component found\n')
        varName = inJson[testKey]
        resultPair = breakSpecialVarNameIntoBaseAndComponent(varName,
                        self.mAddSeparatorToVectorVariableName)
        self.mVariableName = resultPair[0]
        self.mVariableComponent = resultPair[1]
        self.mVectorBaseName = resultPair[2]
        if PhactoriDbg():
          myDebugPrint3('  non comp name: ' + str(self.mVariableName) + '\n')
      else:
        testKey = inBaseString + 'vector magnitude'
        if testKey in inJson:
          variableFoundFlag = True
          self.mVariableIsVectorMagnitude = True
          if PhactoriDbg():
            myDebugPrint3('  vector magnitude found\n')
          if self.mAddSeparatorToVectorVariableName:
            self.mVariableName = inJson[testKey] + GetSeparatorString()
          else:
            self.mVariableName = inJson[testKey]
          self.mVectorBaseName = inJson[testKey]
        else:
          testKey = inBaseString + 'tensor component'
          if testKey in inJson:
            variableFoundFlag = True
            if PhactoriDbg():
              myDebugPrint3('  tensor component found\n')
            self.mVariableName = inJson[testKey]
          else:
            variableFoundFlag = False
            if PhactoriDbg():
              myDebugPrint3('  no variable found\n')
            self.mVariableName = ''
            self.mComponent = None

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



class PhactoriDataCubeOneAxisInfo:
  def __init__(self):
    self.mUseLabelFlag = False
    self.mAxisLabel = ""
    self.mBaseAxisLabel = None

  def DcoaiParseParametersFromJson(self, parentAxesInfo, inJsn, inXyzKey):

    theKey = 'show ' + inXyzKey + ' axis tic marks'
    parentAxesInfo.mShowTicks = localGet1or0(inJsn, theKey, True)
    #parentAxesInfo.mShowEdges = 1
    #theKey = 'show ' + inXyzKey + ' axis minor tic marks'

    theKey = 'show ' + inXyzKey + ' axis label'
    value = getParameterFromBlock(inJsn, theKey, True)
    if value == False:
      self.mUseLabelFlag = True
      self.mAxisLabel = ""
    else:
      self.mUseLabelFlag = False

    theKey = inXyzKey + ' axis label name'
    if theKey in inJsn:
      self.mUseLabelFlag = True
      self.mAxisLabel = getParameterFromBlock(inJsn, theKey, "")

    if PhactoriDbg():
      myDebugPrint3('DcoaiParseParametersFromJson  ' + inXyzKey + '  ' \
          + str(parentAxesInfo.mShowTicks) + '   ' \
          + str(self.mUseLabelFlag) + '   ' \
          + str(self.mAxisLabel) + '\n')

class PhactoriDataCubeAxesInfo:
  def __init__(self):
    self.mShowEdges = 1
    self.mShowTicks = 1
    self.mXAxisInfo = PhactoriDataCubeOneAxisInfo()
    self.mYAxisInfo = PhactoriDataCubeOneAxisInfo()
    self.mZAxisInfo = PhactoriDataCubeOneAxisInfo()

  def DcaiParseParametersFromJson(self, inJsn):
    self.mXAxisInfo.DcoaiParseParametersFromJson(self, inJsn, 'x')
    self.mYAxisInfo.DcoaiParseParametersFromJson(self, inJsn, 'y')
    self.mZAxisInfo.DcoaiParseParametersFromJson(self, inJsn, 'z')


class PhactoriColorMapSettings:
  def __init__(self):
    self.mColorMapNameId = "Default"
    self.mInvertColorMap = False
    #self.mColorMapNameId = "Blue to Yellow"
    #self.mColorMapNameId = "Blue to Red Rainbow"

  def ParseColorMapSettings(self, inJsn):
    self.mColorMapNameId = getParameterFromBlock(inJsn,
        'preset color scale', self.mColorMapNameId)
    self.mInvertColorMap = getParameterFromBlock(inJsn,
        'invert color scale', self.mInvertColorMap)

class PhactoriRepresentationBlock:
  def __init__(self):
    self.mName = ""
    self.mColorVariableInfo = PhactoriVariableInfo()
    self.mColorLegendFlag = True
    #self.mColorLegendPositionAndSize = ['bottom', 1.0]
    #self.mColorLegendPositionAndSize = ['right', 1.0]
    self.mColorLegendPositionAndSize = ['bottom right', 1.0]
    self.mTimeAnnotationSettings = PhactoriAnnotationViewSettings(
        'time', 'time')
    self.mDataCubeAxesFlag = False
    self.mDataCubeAxesInfo = PhactoriDataCubeAxesInfo()
    #self.mColorByBlockFlag = True
    self.mColorByBlockFlag = True
    self.mColorByBlockExplicitlySet = False
    self.mColorBySolidColorFlag = False
    self.mSolidColor = None
    self.mOpacitySetting = 1.0
    self.mOrientationAxesFlag = True
    self.mFixedColorRange = None
    self.mUseFixedColorRange = False
    self.mColorRangeMinMaxTracker = PlotValMinMaxTrkC()
    self.mUseHighlightSubranges = False
    self.mHighlightSubranges = []
    self.mUseColorSubrange = False
    self.mColorSubrange = [0.0, 1.0]
    self.mFilenameAddon = ""
    self.mColorSettings = PhactoriColorSettings()
    self.mColorMapSettings = PhactoriColorMapSettings()
    self.mPointSize = 2.0

  def SetFromRestartInfo(self, inJson):
    """given a map (json format), use the info in the map to set the
       representation state--this reads the info created in
       GetRestartInfo"""

    if 'mColorRangeMinMaxTracker' not in inJson:
      if PhactoriDbg():
        myDebugPrint3("PhactoriRepresentationBlock::" + \
            "SetFromRestartInfo: no mColorRangeMinMaxTracker, return\n")
      return

    if PhactoriDbg():
      myDebugPrint3("Representation::SetFromRestartInfo\n" +
        "currently set to do something\n" +
        "before tracker:\n" + self.mColorRangeMinMaxTracker.SelfToStr())

    jsonItem = inJson['mColorRangeMinMaxTracker']
    self.mColorRangeMinMaxTracker.SetFromRestartInfo(jsonItem)

    if PhactoriDbg():
      myDebugPrint3("after tracker:\n" +
          self.mColorRangeMinMaxTracker.SelfToStr())

  def GetRestartInfo(self):
    """construct, in python map/json format, the information from this
       representation instance which contains the information
       which would be needed to restore the representation to the proper
       state after a simulation restart, particularly color range tracking
       information.  Return the restart info map/json"""
    newRestartInfoJson = {}
    newRestartInfoJson['mColorRangeMinMaxTracker'] = \
      self.mColorRangeMinMaxTracker.GetRestartInfo()
    return newRestartInfoJson

  def CalculateDefaultScalarOpacityUnitDistance(self, inPhactoriOperation):
    """given a phactori operation assumed to have an updated pipeline,
       calculate a default scalar opacity unit distance the same as paraview,
       with 0.05 * diagonal length of data bounding box"""
    bnds = GetGlobalDataBoundsParallel(inPhactoriOperation.mParaViewFilter)
    xxlen = bnds[1] - bnds[0]
    yylen = bnds[3] - bnds[2]
    zzlen = bnds[5] - bnds[4]
    self.mScalarOpacityUnitDistance = \
            0.05 * math.sqrt(xxlen*xxlen + yylen*yylen + zzlen*zzlen)

def GetThresholdContourHackVariableNameString(inVariableInfo):
  """gives altered variable name for magnitude or component
     in paraview, coloring by a vector component or magnitude is handled
     differently than thresholding or contouring by a vector component or
     magnitude.  For the Contour and Threshold case, we add _X, _Y, or _Z to
     the variable name for component and _Magnitude for magnitude.  In
     the coloring case, we specify component arguments and another argument
     as to whether to use magnitude or component.  This routine takes the
     given variable info (name, whether or not vector component or magnitude)
     and creates a properly munged name to do the right thing"""

  if inVariableInfo.mVariableIsVectorComponent:
    if inVariableInfo.mVariableComponent == 0:
      localVariableName = inVariableInfo.mVariableName + '_X'
    elif inVariableInfo.mVariableComponent == 1:
      localVariableName = inVariableInfo.mVariableName + '_Y'
    elif inVariableInfo.mVariableComponent == 2:
      localVariableName = inVariableInfo.mVariableName + '_Z'
    else:
      if PhactoriDbg():
        myDebugPrint3("  variable component is not 0, 1, or 2, using 0 (X)\n")
      localVariableName = inVariableInfo.mVariableName + '_X'
  elif inVariableInfo.mVariableIsVectorMagnitude:
    localVariableName = inVariableInfo.mVariableName + '_Magnitude'
  else:
    localVariableName = inVariableInfo.mVariableName

  return localVariableName

class PhactoriOperationSpecifics:
  #base class for operation specifics; probably need to reconfigure so
  #that we have a base operation abstract class rather than an operation
  #base class that has an operation specifics member which is a child
  #class of this class
  def __init__(self):
    if PhactoriDbg():
      myDebugPrint3("PhactoriOperationSpecifics::__init__ executed\n")

    #which PhactoriOperationBlock instance owns/created this specific
    #operation.  Set during ParseOneFilterTypeFromViewMapOperation()
    self.mPhactoriOperationBlockOwner = None

  def DoUpdateDueToChangeInData(self, inIncomingPvFilter,
      outOutgoingPvFilter):
    #default behavior is to do nothing on new data
    if PhactoriDbg():
      myDebugPrint3("PhactoriOperationSpecifics::"
          "DoUpdateDueToChangeInData executed\n")
      return False

  def ExportOperationData(self, datadescription):
    """this will be called once per callback (before WriteImages) to allow the
       operation to export any desired data which is not an image. The child
       class should override this method if it wants so do something"""
    return


class PhactoriThresholdOperation(PhactoriOperationSpecifics):
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
    #    errStr = """error!  inJson 'variable type' is neither node nor element
    #                PhactoriThresholdOperation.ParseParametersFromJson\n"""
    #    myDebugPrint3(errStr)
    #    raise Exception(errStr)
    #else:
    #  errStr = """error!  inJson has no 'variable type' key in
    #              PhactoriThresholdOperation.ParseParametersFromJson\n"""
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
                  """ + str(self.mVariableInfo.mVariableType) + """\n"""
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

    #CellDataList = []
    #for ii in range(inputSource.CellData.GetNumberOfArrays()):
    #    CellDataList.append(inputSource.CellData.GetArray(ii).Name)
    #myDebugPrint3('before threshold cell data items: ' + str(CellDataList) + '\n');

    newParaViewFilter = Threshold(inInputFilter)
    newParaViewFilter.ThresholdRange = self.mRange

    self.HandleVariableNeedingDetection(inInputFilter, newParaViewFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)
    #AddFilterToFilterMap(ioOperationBlock.mName, newParaViewFilter)

    #inputSource = GetActiveSource()
    #CellDataList = []
    #for ii in range(inputSource.CellData.GetNumberOfArrays()):
    #    CellDataList.append(inputSource.CellData.GetArray(ii).Name)
    #myDebugPrint3('after threshold cell data items: ' + str(CellDataList) + '\n');

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriThresholdOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

class PhactoriAddUnstructuredGridOperation(PhactoriOperationSpecifics):
  """filter/operation which reads in a .vtu file which is an unstructured
     grid and creates a new source from that and groups it with the
     input source"""
  def __init__(self):
    self.mFilename = None
    self.mInternalPvUnstructuredGridReader = None

  def ParseParametersFromJson(self, inJson):
    if 'filename' in inJson:
      self.mFilename = inJson['filename']
    else:
      myDebugPrint3AndException(
          "PhactoriAddUnstructuredGridOperation::ParseParametersFromJson\n"
          "Error:  must have 'filename' key\n")

  def CreateParaViewFilter(self, inInputFilter):
    """create the read unstructured grid (and group) filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriAddUnstructuredGridOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    #csv reader
    self.mInternalPvUnstructuredGridReader = XMLUnstructuredGridReader( FileName=[self.mFilename] )

    #don't know if this is necessary here
    UpdatePipelineWithCurrentTimeArgument(self.mInternalPvUnstructuredGridReader)

    #now group this new point set source with the original one
    newParaViewFilter = GroupDatasets()
    newParaViewFilter.Input = [self.mInternalPvUnstructuredGridReader, inInputFilter]

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3(
          "filename: " + str(self.mInternalPvUnstructuredGridReader.FileName) + "\n")
      myDebugPrint3("PhactoriAddUnstructuredGridOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

class PhactoriAddPointSetOperation(PhactoriOperationSpecifics):
  """filter/operation which reads in a .csv file which is a set of points and
     creates a new source from that and groups it with the input source"""
  def __init__(self):
    self.mXColumn = 'x column'
    self.mYColumn = 'y column'
    self.mZColumn = 'z column'
    self.m2dPointsFlag = False
    self.mFieldDelimeterCharacters = ','
    self.mFilename = None
    self.mInternalPvCSVReader = None
    self.mInternalPvTableToPoints = None

  def ParseParametersFromJson(self, inJson):
    if 'filename' in inJson:
      self.mFilename = inJson['filename']
    else:
      myDebugPrint3AndException(
          "PhactoriAddPointSetOperation::ParseParametersFromJson\n"
          "Error:  must have 'filename' key\n")
    if 'x column' in inJson:
      self.mXColumn = inJson['x column']
    if 'y column' in inJson:
      self.mYColumn = inJson['y column']
    if 'z column' in inJson:
      self.mZColumn = inJson['z column']
    if '2d points flag' in inJson:
      self.m2dPointsFlag = inJson['2d points flag']
    if 'field delimeter' in inJson:
      self.mFieldDelimeterCharacters = inJson['field delimeter']

  def CreateParaViewFilter(self, inInputFilter):
    """create the read pointset (and group) filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriAddPointSetOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    #csv reader
    self.mInternalPvCSVReader = CSVReader( FileName=[self.mFilename] )
    self.mInternalPvCSVReader.FieldDelimiterCharacters = self.mFieldDelimeterCharacters
    #don't know if this is necessary here
    UpdatePipelineWithCurrentTimeArgument(self.mInternalPvCSVReader)

    #table to points
    self.mInternalPvTableToPoints = TableToPoints(self.mInternalPvCSVReader)
    self.mInternalPvTableToPoints.XColumn = self.mXColumn
    self.mInternalPvTableToPoints.YColumn = self.mYColumn
    self.mInternalPvTableToPoints.ZColumn = self.mZColumn
    if self.m2dPointsFlag:
      self.mInternalPvTableToPoints.a2DPoints = 1
    else:
      self.mInternalPvTableToPoints.a2DPoints = 0
    #don't know if this is necessary here
    UpdatePipelineWithCurrentTimeArgument(self.mInternalPvTableToPoints)

    #now group this new point set source with the original one
    newParaViewFilter = GroupDatasets()
    newParaViewFilter.Input = [self.mInternalPvTableToPoints, inInputFilter]
    #newParaViewFilter = self.mInternalPvTableToPoints

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3(
          "filename: " + str(self.mInternalPvCSVReader.FileName) + "\n"
          "delimiter: -->" + str(self.mInternalPvCSVReader.FieldDelimiterCharacters) + "<--\n"
          "x column: " + str(self.mInternalPvTableToPoints.XColumn) + "\n"
          "y column: " + str(self.mInternalPvTableToPoints.YColumn) + "\n"
          "z column: " + str(self.mInternalPvTableToPoints.ZColumn) + "\n"
          "2d points flag: " + str(self.mInternalPvTableToPoints.a2DPoints) + "\n")
      myDebugPrint3("PhactoriAddPointSetOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

class PhactoriReflectOperation(PhactoriOperationSpecifics):
  """manages reflect filter, also inserts group filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mPlane = "X Min"
    self.mCenter = 0.0

    #we keep internal reflect filter, because group filter is the one exposed
    #to the rest of the world
    self.mInternalReflectFilter = None

  def ParseParametersFromJson(self, inJson):
    if 'plane' in inJson:
      self.mPlane = inJson['plane']
    else:
      myDebugPrint3AndException(
          "PhactoriReflectOperation::ParseParametersFromJson\n"
          "Error:  must have 'plane' key\n")
    if 'center' in inJson:
      self.mCenter = inJson['center']

  def CreateParaViewFilter(self, inInputFilter):
    """create the reflect (and group) filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriReflectOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    self.mInternalReflectFilter = Reflect(inInputFilter)
    self.mInternalReflectFilter.Plane  = self.mPlane
    self.mInternalReflectFilter.Center  = self.mCenter
    #don't know if this is necessary here
    UpdatePipelineWithCurrentTimeArgument(self.mInternalReflectFilter)

    #now group this new source with the original non-reflected one
    newParaViewFilter = GroupDatasets()
    newParaViewFilter.Input = [self.mInternalReflectFilter, inInputFilter]

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3("reflect plane: " + str(self.mInternalReflectFilter.Plane) + "\n"
          "reflect center: " + str(self.mInternalReflectFilter.Center) + "\n")
      myDebugPrint3("PhactoriReflectOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

class PhactoriWarpByVectorOperation(PhactoriOperationSpecifics):
  """manages warpbyvector filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mVectorName = ""
    self.mVariableType = "node"
    self.mScaleFactor = 1.0

  def ParseParametersFromJson(self, inJson):
    if 'variable type' in inJson:
      self.mVariableType = inJson['variable type']
    elif 'variablename' in inJson:
      self.mVectorName = inJson['variablename']
    elif 'variable vector' in inJson:
      self.mVectorName = inJson['variable vector']
    else:
      myDebugPrint3AndException(
          "PhactoriWarpByVectorOperation::ParseParametersFromJson\n"
          "Error:  must have 'vector name' key\n")
    if 'scale' in inJson:
      self.mScaleFactor = inJson['scale']

  def CreateParaViewFilter(self, inInputFilter):
    """create the warp by vector filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriWarpByVectorOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    newParaViewFilter = WarpByVector(inInputFilter)
    newParaViewFilter.ScaleFactor = self.mScaleFactor
    if self.mVariableType == 'cell':
      varTypeString = 'CELLS'
    else:
      varTypeString = 'POINTS'
    newParaViewFilter.Vectors = [varTypeString, self.mVectorName]

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3("vector name: " + self.mVectorName + "\n"
          "variable type: " + varTypeString + "\n"
          "scale factor: " + str(newParaViewFilter.ScaleFactor) + "\n")
      myDebugPrint3("PhactoriWarpByVectorOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

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

    newParaViewFilter = Calculator(inInputFilter)
    newParaViewFilter.Function = self.mFunction
    #newParaViewFilter.AttributeMode = 'Cell Data'
    #newParaViewFilter.AttributeMode = 'Point Data'
    newParaViewFilter.AttributeMode = self.mPointOrCell
    newParaViewFilter.ResultArrayName = self.mResultArrayName

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    if PhactoriDbg(100):
      myDebugPrint3("function: " + str(newParaViewFilter.Function) + "\n"
        "result array name: " + str(newParaViewFilter.ResultArrayName) + "\n"
        "PhactoriTransformOperation.CreateParaViewFilter returning\n", 100)

    return newParaViewFilter

class PhactoriTransformOperation(PhactoriOperationSpecifics):
  """manages transform filter--scale, rotate, translate"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mTranslate = [0.0, 0.0, 0.0]
    self.mRotate = [0.0, 0.0, 0.0]
    self.mScale = [1.0, 1.0, 1.0]

  def ParseParametersFromJson(self, inJson):
    if 'scale' in inJson:
      self.mScale = inJson['scale']
    if 'translate' in inJson:
      self.mTranslate = inJson['translate']
    if 'rotate' in inJson:
      self.mRotate = inJson['rotate']

  def CreateParaViewFilter(self, inInputFilter):
    """create the transform filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriTransformOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    #newParaViewFilter = Transform(inInputFilter, Transform = "Transform")
    SetActiveSource(inInputFilter)
    inInputFilter.UpdatePipeline()
    newParaViewFilter = Transform(Transform = "Transform")
    newParaViewFilter.Transform = "Transform"
    newParaViewFilter.Transform.Scale = self.mScale
    newParaViewFilter.Transform.Translate = self.mTranslate
    newParaViewFilter.Transform.Rotate = self.mRotate

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3(
        "translate: " + str(newParaViewFilter.Transform.Translate) + "\n"
        "rotate: " + str(newParaViewFilter.Transform.Rotate) + "\n"
        "scale: " + str(newParaViewFilter.Transform.Scale) + "\n"
        "input: " + str(newParaViewFilter.Input) + "\n")
      myDebugPrint3("PhactoriTransformOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

class SharedTriangleSet():
  """idea is that this encapsulates a set of triangles which are the triangles
     from the overall surface of interest which reside on a single processor.
     You can store your local processor triangles in this class and broadcast
     from them and you can receive a broadcast of this setup from other
     processes"""

  def __init__(self):
    self.numPoints = -1
    self.numTriangles = -1
    #triangle points in xyzxyzxyzxyz order
    self.PointXyzs = vtk.vtkDoubleArray()
    #global node id for each point to help avoid self intersection
    self.NodeIds = vtk.vtkIntArray()
    #triplets of indexes into self.PointXyzs to represent triangles, in
    #abcabcabcabc order
    self.Triangles = vtk.vtkIntArray()
    self.fromPid = -1

    #to share num points and num tris
    self.shareSizes = vtk.vtkIntArray()
    self.shareSizes.SetNumberOfValues(2)

    self.BspTree = None
    self.UseBSPTreeForRayIntersection = True

  def CreateFromLocalProcess(self, inPhactoriOperation):
    """figure out the set of triangles from the target surface which are on
       this process (if any): we assume we have a triangle mesh or this code
       won't work"""

    #obtain pointer to the local geometry data
    csdata = inPhactoriOperation.mParaViewFilter.GetClientSideObject().\
        GetOutputDataObject(0)

    self.numPoints = csdata.GetNumberOfPoints()
    if PhactoriDbg():
      myDebugPrint3(str(dir(csdata)) + "\n")
      myDebugPrint3("num points: " + str(self.numPoints) + "\n")
      myDebugPrint3(str(dir(vtk)) + "\n")

    pntData = csdata.GetPointData()
    cellData = csdata.GetCellData()
    numCells = csdata.GetNumberOfCells()
    gNodeIdArray = pntData.GetArray('GlobalNodeId')
    #gElmtIdArray = cellData.GetArray('GlobalElementId')
    #pntGeometryArray = csdata.GetPoints()

    self.PointXyzs.SetNumberOfValues(self.numPoints*3)
    self.NodeIds = vtk.vtkIntArray()
    self.NodeIds.SetNumberOfValues(self.numPoints)

    #this is stupid, there is probably a much faster way to do this
    ptxyz = [0.0,0.0,0.0]
    for ii in xrange(0, self.numPoints):
      ndx = ii*3
      csdata.GetPoint(ii,ptxyz)
      self.PointXyzs.SetValue(ndx, ptxyz[0])
      self.PointXyzs.SetValue(ndx+1, ptxyz[1])
      self.PointXyzs.SetValue(ndx+2, ptxyz[2])
      if(gNodeIdArray == None):
        self.NodeIds.SetValue(ii, ii)
      else:
        self.NodeIds.SetValue(ii, gNodeIdArray.GetValue(ii))

    self.Triangles.SetNumberOfValues(0)
    cellPointIds = vtk.vtkIdList()
    for ii in xrange(0, numCells):
      csdata.GetCellPoints(ii, cellPointIds)
      #numpoints should be 3
      numids = cellPointIds.GetNumberOfIds()
      #we are only doing triangles
      if numids != 3:
        if numids < 3:
          #degenerate ? try just skipping
          if PhactoriDbg():
            myDebugPrint3AndException(str(ii) + " degenerate 2 point\n")
          continue
        if True: #for now we consider this fatal error
          myDebugPrint3AndException(
            "PhactoriIntersectNodeNormalsWithSurface::CreateFromLocalProcess\n"
            "encountered non-triangle\n")
        continue
      self.Triangles.InsertNextValue(cellPointIds.GetId(0))
      self.Triangles.InsertNextValue(cellPointIds.GetId(1))
      self.Triangles.InsertNextValue(cellPointIds.GetId(2))
    self.numTriangles = self.Triangles.GetNumberOfValues() // 3

  def SendBroadcast(self, inLocalPid, globalController):
    """broadcast the triangle set from the local process to all the other
       processes"""
    if (self.numPoints < 0) or (self.numTriangles < 0):
      myDebugPrint3AndException(
        "PhactoriIntersectNodeNormalsWithSurface::SendBroadcast\n"
        "can't broadcast without being filled in first\n")

    if PhactoriDbg():
      myDebugPrint3("SharedTriangleSet::SendBroadcast entered, inLocalPid " + \
          str(inLocalPid) + \
          "\nself.numPoints " + str(self.numPoints) + \
          "  self.numTriangles " + str(self.numTriangles) + "\n")

    #first broadcast number of points and number of triangles
    self.shareSizes.SetValue(0, self.numPoints)
    self.shareSizes.SetValue(1, self.numTriangles)

    globalController.Broadcast(self.shareSizes, inLocalPid)

    self.fromPid = inLocalPid

    if self.numTriangles == 0:
      if PhactoriDbg():
        myDebugPrint3(
          "SharedTriangleSet::SendBroadcast returning, nothing to send\n")
      return

    globalController.Broadcast(self.PointXyzs, inLocalPid)
    globalController.Broadcast(self.NodeIds, inLocalPid)
    globalController.Broadcast(self.Triangles, inLocalPid)

    if PhactoriDbg():
      myDebugPrint3("SharedTriangleSet::SendBroadcast returning, inLocaPid " \
        + str(inLocalPid) + "\n")

  def ReceiveBroadcast(self, inFromPid, globalController):
    """receive a triangle set broadcast from the another process to the local
       process"""
    if PhactoriDbg():
      myDebugPrint3(
        "SharedTriangleSet::ReceiveBroadcast entered, inFromPid " + \
        str(inFromPid) + "\n")

    globalController.Broadcast(self.shareSizes, inFromPid)
    self.fromPid = inFromPid

    self.numPoints = self.shareSizes.GetValue(0)
    self.numTriangles = self.shareSizes.GetValue(1)

    if self.numTriangles == 0:
      if PhactoriDbg():
        myDebugPrint3(
          "SharedTriangleSet::ReceiveBroadcast returning, nothing to get\n")
      return

    self.PointXyzs.SetNumberOfValues(3*self.numPoints)
    self.NodeIds.SetNumberOfValues(self.numPoints)
    self.Triangles.SetNumberOfValues(3*self.numTriangles)

    globalController.Broadcast(self.PointXyzs, inFromPid)
    globalController.Broadcast(self.NodeIds, inFromPid)
    globalController.Broadcast(self.Triangles, inFromPid)

    if PhactoriDbg():
      myDebugPrint3(
          "SharedTriangleSet::ReceiveBroadcast returning, inFromPid " + \
          str(inFromPid) + \
          "\nself.numPoints " + str(self.numPoints) + \
          "  self.numTriangles " + str(self.numTriangles) + "\n")

  def IntersectWithOneRayViaBSPTree(self, inIndex, inRayStart, inRayDirection,
        ioRaycastOperation, inPointOrCellData,
        inRayEnd1, inRayEnd2, inBothDirFlag):
    if PhactoriDbg():
      myDebugPrint3(
        "bspt ray: " + str(inIndex) + " " + str(inRayStart) + "  " + \
        str(inRayDirection) + " " + str(inPointOrCellData) + "\n")

    #BSPTree interect with line
    #V.IntersectWithLine((float, float, float), (float, float, float),
    #    float, vtkPoints, vtkIdList) -> int
    #C++: virtual int IntersectWithLine(const double p1[3],
    #    const double p2[3], const double tol, vtkPoints *points,
    #    vtkIdList *cellIds)
    hitPoints = vtk.vtkPoints()
    hitIds = vtk.vtkIdList()
    segmentStart = vecMultiplyAdd(inRayStart, inRayDirection, inRayEnd1)
    segmentEnd = vecMultiplyAdd(inRayStart, inRayDirection, inRayEnd2)
    self.BspTree.IntersectWithLine(segmentStart, segmentEnd, 0.0001,
      hitPoints, hitIds)
    numhits = hitPoints.GetNumberOfPoints()
    if PhactoriDbg():
      myDebugPrint3("results from self.BspTree.IntersectWithLine:\n"
        "num hits: " + str(numhits) + "\n"
        "index, id, xyz\n")
      for hhpp in range(0, numhits):
        myDebugPrint3(str(hhpp) + ", " + str(hitIds.GetId(hhpp)) + ", " + \
          str(hitPoints.GetPoint(hhpp)) + "\n")

    numhitsB = 0
    if inBothDirFlag:
      hitPointsB = vtk.vtkPoints()
      hitIdsB = vtk.vtkIdList()
      segmentStart = vecMultiplyAdd(inRayStart, inRayDirection, -inRayEnd1)
      segmentEnd = vecMultiplyAdd(inRayStart, inRayDirection, -inRayEnd2)
      self.BspTree.IntersectWithLine(segmentStart, segmentEnd, 0.0001,
        hitPointsB, hitIdsB)
      numhitsB = hitPointsB.GetNumberOfPoints()
      if PhactoriDbg():
        myDebugPrint3("other direction self.BspTree.IntersectWithLine:\n"
          "num hits B: " + str(numhitsB) + "\n"
          "index, id, xyz B\n")
        for hhpp in range(0, numhitsB):
          myDebugPrint3(str(hhpp) + ", " + str(hitIdsB.GetId(hhpp)) + ", " + \
            str(hitPointsB.GetPoint(hhpp)) + "\n")

    #find nearest hit point
    if (numhits > 0) or (numhitsB > 0):
      hitxyz = [0.0,0.0,0.0]
      hit2xyz = [0.0,0.0,0.0]
      distsqrd = sys.float_info.max
      for ii in range(0,numhits):
        hitPoints.GetPoint(ii, hit2xyz)
        distsqrd2 = vecDistanceSquared(inRayStart, hit2xyz)
        if(distsqrd2 < distsqrd):
          distsqrd = distsqrd2
          vecCopy(hitxyz, hit2xyz)
      for ii in range(0,numhitsB):
        hitPointsB.GetPoint(ii, hit2xyz)
        distsqrd2 = vecDistanceSquared(inRayStart, hit2xyz)
        if(distsqrd2 < distsqrd):
          distsqrd = distsqrd2
          vecCopy(hitxyz, hit2xyz)

      if PhactoriDbg():
        myDebugPrint3(
          "nearest hit: " + str(distsqrd) + "  " + str(hitxyz) + "\n")

      if ioRaycastOperation.mCreateIntersectionSegmentGroup:
        if inPointOrCellData <= 0:
          ioRaycastOperation.SegmentsFromPointRays.TestAndSet(
            inIndex, distsqrd, inRayStart, hitxyz)
        else:
          ioRaycastOperation.SegmentsFromTriangleRays.TestAndSet(
            inIndex, distsqrd, inRayStart, hitxyz)
      else:
        hitdist = math.sqrt(distsqrd)
        if inPointOrCellData <= 0:
          oldvv = ioRaycastOperation.mPointNormalRayIntersectDistanceArray.\
            GetValue(inIndex)
          if (oldvv < 0.0) or (hitdist < oldvv):
            ioRaycastOperation.mPointNormalRayIntersectDistanceArray.\
              SetValue(inIndex, hitdist)
        else:
          oldvv = ioRaycastOperation.mCellNormalRayIntersectDistanceArray.\
            GetValue(inIndex)
          if (oldvv < 0.0) or (hitdist < oldvv):
            ioRaycastOperation.mCellNormalRayIntersectDistanceArray.\
              SetValue(inIndex, hitdist)

  def IntersectWithOneRay(self, inIndex, inRayStart, inRayDirection,
        ioRaycastOperation, inPointOrCellData):
    """this is not used, but we are leaving it in as a reference for how to
       do a ray/triangle intersection directly"""
    if PhactoriDbg():
      myDebugPrint3(
        "ray: " + str(inIndex) + " " + str(inRayStart) + "  " + \
        str(inRayDirection) + " " + str(inPointOrCellData) + "\n")

    ptA = [0.0,0.0,0.0]
    ptB = [0.0,0.0,0.0]
    ptC = [0.0,0.0,0.0]
    hitPoint = [0.0, 0.0, 0.0, 0.0]
    saveHit = [0.0, 0.0, 0.0, -1.0]
    #if PhactoriDbg():
    #  myDebugPrint3(
    #    "self.Triangles.GetNumberOfValues() " + str(self.Triangles.GetNumberOfValues()) + "\n")
    #  myDebugPrint3(
    #    "self.PointXyzs.GetNumberOfValues() " + str(self.PointXyzs.GetNumberOfValues()) + "\n")
    hitCount = 0
    for ii in xrange(0, self.numTriangles):
      ndx = ii * 3
      ptndxA = self.Triangles.GetValue(ndx)
      ptndxB = self.Triangles.GetValue(ndx+1)
      ptndxC = self.Triangles.GetValue(ndx+2)
      #if PhactoriDbg():
      #  myDebugPrint3("ptndxABC: " + str(ptndxA) + " " + str(ptndxB) + " " + str(ptndxC) + "\n");
      ptndxA *= 3
      ptA[0] = self.PointXyzs.GetValue(ptndxA)
      ptA[1] = self.PointXyzs.GetValue(ptndxA+1)
      ptA[2] = self.PointXyzs.GetValue(ptndxA+2)
      ptndxB *= 3
      ptB[0] = self.PointXyzs.GetValue(ptndxB)
      ptB[1] = self.PointXyzs.GetValue(ptndxB+1)
      ptB[2] = self.PointXyzs.GetValue(ptndxB+2)
      ptndxC *= 3
      ptC[0] = self.PointXyzs.GetValue(ptndxC)
      ptC[1] = self.PointXyzs.GetValue(ptndxC+1)
      ptC[2] = self.PointXyzs.GetValue(ptndxC+2)
      #if PhactoriDbg():
      #  myDebugPrint3("try " + str(ii) + " " + str(ptA) + "\n")
      hitFlag = self.IntersectRayTriangle(
                  inRayStart, inRayDirection, ptA, ptB, ptC, hitPoint)
      if hitFlag:
        hitCount += 1
        if PhactoriDbg():
          myDebugPrint3("Hit " + str(hitCount) + "! " + str(hitPoint) + "\n")
          vvx = vecFromAToB(inRayStart, hitPoint)
          vvxdist = math.sqrt(vecDotProduct(vvx,vvx))
          myDebugPrint3("calculated distance: " + str(vvxdist) + "\n")
        if ioRaycastOperation.mCreateIntersectionSegmentGroup:
          #onelinesource = Line()
          #onelinesource.Point1 = inRayStart
          #onelinesource.Point2 = [hitPoint[0], hitPoint[1], hitPoint[2]]
          #ioRaycastOperation.mGroupLineSource.Input.append(onelinesource)
          if (saveHit[3] < 0.0) or (hitPoint[3] < saveHit[3]):
            saveHit[0] = hitPoint[0]
            saveHit[1] = hitPoint[1]
            saveHit[2] = hitPoint[2]
            saveHit[3] = hitPoint[3]
        else:
          if inPointOrCellData <= 0:
            oldvv = ioRaycastOperation.mPointNormalRayIntersectDistanceArray.\
              GetValue(inIndex)
            newvv = hitPoint[3]
            if (oldvv < 0.0) or (newvv < oldvv):
              ioRaycastOperation.mPointNormalRayIntersectDistanceArray.\
                SetValue(inIndex, hitPoint[3])
          else:
            oldvv = ioRaycastOperation.mCellNormalRayIntersectDistanceArray.\
              GetValue(inIndex)
            newvv = hitPoint[3]
            if (oldvv < 0.0) or (newvv < oldvv):
              ioRaycastOperation.mCellNormalRayIntersectDistanceArray.\
                SetValue(inIndex, hitPoint[3])

    if ioRaycastOperation.mCreateIntersectionSegmentGroup:
      if saveHit[3] >= 0.0:
        onelinesource = Line()
        onelinesource.Point1 = inRayStart
        onelinesource.Point2 = [saveHit[0], saveHit[1], saveHit[2]]
        ioRaycastOperation.mGroupLineSource.Input.append(onelinesource)


  def IntersectRayTriangleNumpy(self, inRayStart, inRayDirection, ptA, ptB, ptC, outPt):
    edge1 = ptB - ptA
    edge2 = ptC - ptA
    EPSILON = 0.0000001
    hh = np.cross(inRayDirection, edge2)
    aa = np.dot(edge1, hh)
    if (aa > -EPSILON) and (aa < EPSILON):
      return False
    ff = 1.0/aa
    ss = inRayStart - ptA
    uu = np.dot(ss, hh)
    uu *= ff
    if (uu < 0.0) or (uu > 1.0):
      return False;
    qq = np.cross(ss, edge1)
    vv = np.dot(inRayDirection, qq)
    vv *= ff
    if (vv < 0.0) or (uu + vv > 1.0):
      return False;
    #At this stage we can compute t to find out where the intersection point
    #is on the line.
    tt = np.dot(edge2, qq)
    tt *= ff
    if (tt > EPSILON):
      #ray intersection
      outPt[0] = inRayStart[0] + inRayDirection[0] * tt
      outPt[1] = inRayStart[1] + inRayDirection[1] * tt
      outPt[2] = inRayStart[2] + inRayDirection[2] * tt
      outPt[3] = tt
      return True
    else:
      #This means that there is a line intersection but not a ray intersection
      return False

  def IntersectRayTriangle(self, inRayStart, inRayDirection, ptA, ptB, ptC, outPt):
    #if PhactoriDbg():
    #  myDebugPrint3("IntersectRayTriangle entered " + str(ptA) + "\n")
    #c function from wikipedia translated to python
    edge1 = vecFromAToB(ptA, ptB)
    edge2 = vecFromAToB(ptA, ptC)
    EPSILON = 0.0000001
    hh = vecCrossProduct(inRayDirection, edge2)
    aa = vecDotProduct(edge1, hh)
    if (aa > -EPSILON) and (aa < EPSILON):
      return False
    ff = 1.0/aa
    ss = vecFromAToB(ptA, inRayStart)
    uu = vecDotProduct(ss, hh)
    uu *= ff
    if (uu < 0.0) or (uu > 1.0):
      return False;
    qq = vecCrossProduct(ss, edge1)
    vv = vecDotProduct(inRayDirection, qq)
    vv *= ff
    if (vv < 0.0) or (uu + vv > 1.0):
      return False;
    #At this stage we can compute t to find out where the intersection point
    #is on the line.
    tt = vecDotProduct(edge2, qq)
    tt *= ff
    if (tt > EPSILON):
      #ray intersection
      outPt[0] = inRayStart[0] + inRayDirection[0] * tt
      outPt[1] = inRayStart[1] + inRayDirection[1] * tt
      outPt[2] = inRayStart[2] + inRayDirection[2] * tt
      outPt[3] = tt
      return True
    else:
      #This means that there is a line intersection but not a ray intersection
      return False

  def CreateBspTree(self):
    if PhactoriDbg():
      myDebugPrint3("CreateBspTree entered, creating\n")

    if self.BspTree != None:
      #clear/delete bsp tree to avoid memory leak; may be unnecessary
      if PhactoriDbg():
        myDebugPrint3("clearing old bsp tree\n")
      self.BspTree.FreeSearchStructure()
      #?Delete(self.BspTree)
      self.BspTree = None

    self.BspTreePolyData = vtk.vtkPolyData()

    tmpPoints = vtk.vtkPoints()
    ptxyz = self.PointXyzs
    numvls = ptxyz.GetNumberOfValues()
    for ii in xrange(0, numvls, 3):
      tmpPoints.InsertNextPoint(ptxyz.GetValue(ii), ptxyz.GetValue(ii+1),
                             ptxyz.GetValue(ii+2))

    tmpTriangles = vtk.vtkCellArray()
    tris = self.Triangles
    numtripts = tris.GetNumberOfValues()
    for ii in xrange(0, numtripts, 3):
      newtri = vtk.vtkTriangle()
      newtri.GetPointIds().SetId(0, tris.GetValue(ii))
      newtri.GetPointIds().SetId(1, tris.GetValue(ii+1))
      newtri.GetPointIds().SetId(2, tris.GetValue(ii+2))
      tmpTriangles.InsertNextCell(newtri)

    self.BspTreePolyData.SetPoints(tmpPoints)
    self.BspTreePolyData.SetPolys(tmpTriangles)
    self.BspTree = vtk.vtkFiltersFlowPaths.vtkModifiedBSPTree()
    self.BspTree.SetDataSet(self.BspTreePolyData);
    self.BspTree.BuildLocator();
    if PhactoriDbg():
      myDebugPrint3("dir for self.BspTree:\n")
      myDebugPrint3(str(dir(self.BspTree)) + "\n")
      myDebugPrint3("str(self.BspTree):\n")
      myDebugPrint3(str(self.BspTree) + "\n")
      myDebugPrint3("CreateBspTree returning\n")

  def IntersectWithNormalRaysFromSourceTriangles(self,
        inSourcePhactoriOp, ioRaycastOperation, inSourceSharedTriangles):
    if PhactoriDbg():
      myDebugPrint3(
        "IntersectWithNormalRaysFromSourceTriangles entered\n" \
        "target triangle set from pid " + str(self.fromPid) + "\n"\
        "source triangle set from pid " + \
          str(inSourceSharedTriangles.fromPid) + "\n"\
        "target num triangles, num points: " + str(self.numTriangles) + ", " +\
        str(self.numPoints) + "\n" + \
        "target num triangles, num points: " + str(self.numTriangles) + ", " +\
        str(self.numPoints) + "\n")

    csdata = inSourcePhactoriOp.mParaViewFilter.GetClientSideObject().\
        GetOutputDataObject(0)

    numSourceCells = csdata.GetNumberOfCells()
    if PhactoriDbg():
      myDebugPrint3("inSourcePhactoriOp.mParaViewFilter: " + \
        str(inSourcePhactoriOp.mParaViewFilter) + "\n")
      myDebugPrint3("csdata: " + str(csdata) + "\n")
      myDebugPrint3("numSourceCells: " + str(numSourceCells) + "\n")

    if numSourceCells <= 0:
      #nothing here, but maybe need to create empty item
      return

    cellData = csdata.GetPointData()

    nrmls = cellData.GetArray("Normals")
    if nrmls == None:
      myDebugPrint3AndException("IntersectWithNormalRaysFromSourceTriangles: " \
        "no Normals cell data\n")

    #pointXyza = [0.0,0.0,0.0]
    #pointXyzb = [0.0,0.0,0.0]
    #pointXyzc = [0.0,0.0,0.0]
    cntrd = [0.0, 0.0, 0.0]
    onethird = 1.0/3.0
    cellPointIds = vtk.vtkIdList()
    for ii in range(0,numSourceCells):
      theNormalTuple = nrmls.GetTuple3(ii)
      ndx = ii*3
      ##csdata.GetPoint(ii, pointXyz)
      #ptndxa = inSourceSharedTriangles.Triangles.GetValue(ndx) * 3
      #ptndxb = inSourceSharedTriangles.Triangles.GetValue(ndx+1) * 3
      #ptndxc = inSourceSharedTriangles.Triangles.GetValue(ndx+2) * 3
      #cntrd[0] = inSourceSharedTriangles.PointXyzs.GetValue(ptndxa)
      #cntrd[0] += inSourceSharedTriangles.PointXyzs.GetValue(ptndxb)
      #cntrd[0] += inSourceSharedTriangles.PointXyzs.GetValue(ptndxc)
      #cntrd[0] *= onethird
      #cntrd[1] = inSourceSharedTriangles.PointXyzs.GetValue(ptndxa+1)
      #cntrd[1] += inSourceSharedTriangles.PointXyzs.GetValue(ptndxb+1)
      #cntrd[1] += inSourceSharedTriangles.PointXyzs.GetValue(ptndxc+1)
      #cntrd[1] *= onethird
      #cntrd[2] = inSourceSharedTriangles.PointXyzs.GetValue(ptndxa+2)
      #cntrd[2] += inSourceSharedTriangles.PointXyzs.GetValue(ptndxb+2)
      #cntrd[2] += inSourceSharedTriangles.PointXyzs.GetValue(ptndxc+2)
      #cntrd[2] *= onethird
      csdata.GetCellPoints(ii, cellPointIds)
      #numpoints should be 3
      numids = cellPointIds.GetNumberOfIds()
      if numids != 3:
        if numids < 3:
          #degenerate ? try just skipping
          if PhactoriDbg():
            myDebugPrint3AndException(str(ii) + " degenerate 2 point\n")
          continue
        if True: #for now we consider this fatal error
          myDebugPrint3AndException(
            "PhactoriIntersectNodeNormalsWithSurface::CreateFromLocalProcess\n"
            "encountered non-triangle\n")
        continue
      ptndxa = cellPointIds.GetId(0)
      ptndxb = cellPointIds.GetId(1)
      ptndxc = cellPointIds.GetId(2)
      #trying to get it to work just use first point on triangle
      csdata.GetPoint(ptndxa, cntrd)

      if self.UseBSPTreeForRayIntersection:
        self.IntersectWithOneRayViaBSPTree(ii, cntrd, theNormalTuple,
          ioRaycastOperation, 1, -0.1, 10.0, False)
          #ioRaycastOperation, 1, 0.0001, 4.0, False)
      else:
        self.IntersectWithOneRay(
          ii, cntrd, theNormalTuple, ioRaycastOperation, 1)

    if PhactoriDbg():
      myDebugPrint3("IntersectWithNormalRaysFromSourceTriangles returning\n")

  def IntersectWithNormalRaysFromSource(self, inSourcePhactoriOp,
        ioRaycastOperation):
    if PhactoriDbg():
      myDebugPrint3("IntersectWithNormalRaysFromSource entered\n" \
        "triangle set from pid " + str(self.fromPid) + "\n"\
        "num triangles, num points: " + str(self.numTriangles) + ", " +
        str(self.numPoints) + "\n")

    csdata = inSourcePhactoriOp.mParaViewFilter.GetClientSideObject().\
        GetOutputDataObject(0)

    if self.UseBSPTreeForRayIntersection:
      self.CreateBspTree()

    numSourcePoints = csdata.GetNumberOfPoints()
    if PhactoriDbg():
      myDebugPrint3("inSourcePhactoriOp.mParaViewFilter: " + str(inSourcePhactoriOp.mParaViewFilter) + "\n")
      myDebugPrint3("csdata: " + str(csdata) + "\n")
      myDebugPrint3("numSourcePoints: " + str(numSourcePoints) + "\n")

    if numSourcePoints <= 0:
      #nothing here, but maybe need to create empty item
      return

    pntData = csdata.GetPointData()

    nrmls = pntData.GetArray("Normals")
    if nrmls == None:
      myDebugPrint3AndException(
          "IntersectWithNormalRaysFromSource: no Normals point data\n")

    pointXyz = [0.0,0.0,0.0]
    #pointXyz = np.array([0.0,0.0,0.0])

    for ii in range(0,numSourcePoints):
      theNormalTuple = nrmls.GetTuple3(ii)
      #theNormalTuple = np.array(nrmls.GetTuple3(ii))
      csdata.GetPoint(ii, pointXyz)
      if self.UseBSPTreeForRayIntersection:
        self.IntersectWithOneRayViaBSPTree(ii, pointXyz, theNormalTuple,
          ioRaycastOperation, 0, -0.1, 10.0, False)
          #ioRaycastOperation, 0, 0.0001, 4.0, False)
      else:
        self.IntersectWithOneRay(
          ii, pointXyz, theNormalTuple, ioRaycastOperation, 0)

    if PhactoriDbg():
      myDebugPrint3("IntersectWithNormalRaysFromSource returning\n")


class PhactoriSegmentGroup1Item:
  """one segment in PhactoriSegmentGroup1"""
  def __init__(self):
    self.ParaViewLineSource = None
    self.Point1 = [0.0, 0.0, 0.0]
    self.Point2 = [0.0, 0.0, 0.0]
    self.LengthSquared = sys.float_info.max
    self.mVtkLine = vtk.vtkLine()

  def UpdateVtkLine(self, inIdx, ioVtkPoints):
    idx = inIdx*2
    self.mVtkLine.GetPointIds().SetId(0,idx)
    self.mVtkLine.GetPointIds().SetId(1,idx+1)
    ioVtkPoints.SetPoint(idx, self.Point1)
    ioVtkPoints.SetPoint(idx+1, self.Point2)

  def TestAndSet(self, inLengthSquared, inPt1, inPt2):
    if inLengthSquared < self.LengthSquared:
      self.LengthSquared = inLengthSquared
      vecCopy(self.Point1, inPt1)
      vecCopy(self.Point2, inPt2)

class PhactoriSegmentGroup1:
  """manage a set of segments for purposes of showing raycast results"""
  def __init__(self):
    self.Segments = []
    self.mVtkPoints = vtk.vtkPoints()
    self.mVtkPolyData = vtk.vtkPolyData()
    self.mVtkCellArray = vtk.vtkCellArray()

  def SetupWithNumberOfItems(self, inCount):
    #we need to get arrays the right length and delete extra line sources
    while len(self.Segments) < inCount:
      self.Segments.append(PhactoriSegmentGroup1Item())
    while len(self.Segments) > inCount:
      pp = self.Segments.pop()
      #do we need this?
      #Delete(self.mVtkLine)
      #or this?
      #del self.mVtkLine
      exit(-1)

    self.mVtkPoints.SetNumberOfPoints(inCount * 2)

    for ii in self.Segments:
      ii.LengthSquared = sys.float_info.max

  def UpdateParaViewSegments(self, inSegmentMinLength):
    self.mVtkCellArray.SetNumberOfCells(0)
    for idx, ii in enumerate(self.Segments):
      if PhactoriDbg():
        myDebugPrint3("UpdateParaViewSegments " + str(idx) + "\n")
      ii.UpdateVtkLine(idx, self.mVtkPoints)
      tstminsqrd = inSegmentMinLength * inSegmentMinLength
      if ii.LengthSquared != sys.float_info.max:
        if tstminsqrd <= ii.LengthSquared:
          self.mVtkCellArray.InsertNextCell(ii.mVtkLine)

    self.mVtkPolyData.SetPoints(self.mVtkPoints)
    self.mVtkPolyData.SetLines(self.mVtkCellArray)

  def TestAndSet(self, inIndex, inLengthSquared, inPt1, inPt2):
    self.Segments[inIndex].TestAndSet(inLengthSquared, inPt1, inPt2)

class PhactoriVtpDataExporterOperation(PhactoriOperationSpecifics):
  """passes through it's input to it's output, but also does an
     ExportOperation which is a parallel vtk writer (.vtp)"""
  def __init__(self):
    self.mOutputFileBasename = "myvtp1"
    self.mFilterToWriteDataFrom = None
    self.mWriter = None

  def ParseParametersFromJson(self, inJson):
    dummy = 1

  def CreateParaViewFilter(self, inInputFilter):
    self.mFilterToWriteDataFrom = inInputFilter
    return self.mFilterToWriteDataFrom

  def ExportOperationData(self, datadescription):
    """this will be called once per callback (before WriteImages) to allow the
       operation to export any desired data which is not an image. The child
       class should override this method if it wants so do something.
       For PhactoriIntersectNodeNormalsWithSurface we will output information
       about the nearest and furthest intersections. In this case we dump
       a .vtp or .vtm file."""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriVtpDataExporterOperation." \
          "ExportOperationData entered\n", 100)
    self.mOutputFileBasename = self.mPhactoriOperationBlockOwner.mName
    global gPipeAndViewsState
    outfilename = self.mOutputFileBasename + "_" + \
        str(gPipeAndViewsState.mFrameTagCounter).zfill(4) + ".vtp"
    self.mFilterToWriteDataFrom.UpdatePipeline()

    #icsdClassname = inInputCsData.GetClassName()
    #if icsdClassname == "vtkMultiBlockDataSet" or \
    #   icsdClassname == "vtkExodusIIMultiBlockDataSet":


    if PhactoriDbg(100):
      myDebugPrint3("outfilename: " + outfilename + "\n" \
      "self.mFilterToWriteDataFrom: " + str(self.mFilterToWriteDataFrom) + "\n")

    if self.mWriter == None:
      mypid = SmartGetLocalProcessId()

      import vtkParallelCorePython
      pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
      globalController = pm.GetGlobalController()
      #gLocalProcessId = globalController.GetLocalProcessId()
      numproc = globalController.GetNumberOfProcesses()

      if PhactoriDbg(100):
        myDebugPrint3("numproc " + str(numproc) + " mypid " + str(mypid) + "\n")
      self.mWriter = vtk.vtkXMLPPolyDataWriter()
      self.mWriter.SetInputData(self.mFilterToWriteDataFrom.GetClientSideObject().GetOutputDataObject(0))
      self.mWriter.SetNumberOfPieces(numproc)
      self.mWriter.SetStartPiece(mypid)
      self.mWriter.SetEndPiece(mypid)
      self.mWriter.SetUseSubdirectory(True)
      #if mypid == 0:
      #  self.mWriter.SetWriteSummaryFile(1)
      #else:
      #  self.mWriter.SetWriteSummaryFile(0)

    self.mWriter.SetFileName(outfilename)
    self.mWriter.Write()

      #self.mWriter = CreateWriter(outfilename, self.mFilterToWriteDataFrom)
    #writer = XMLMultiBlockDataWriter(FileName = outfilename,
    #    Input = self.mFilterToWriteDataFrom)
    #self.mWriter.UpdatePipeline()
    #del writer
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriVtpDataExporterOperation." \
          "ExportOperationData returning\n", 100)

class PhactoriVtmDataExporterOperation(PhactoriOperationSpecifics):
  """passes through it's input to it's output, but also does an
     ExportOperation which is a parallel vtk writer (.vtm)"""
  def __init__(self):
    self.mOutputFileBasename = "myvtm1"
    self.mFilterToWriteDataFrom = None
    self.mWriter = None

  def ParseParametersFromJson(self, inJson):
    dummy = 1

  def CreateParaViewFilter(self, inInputFilter):
    self.mFilterToWriteDataFrom = inInputFilter
    return self.mFilterToWriteDataFrom

  def ExportOperationData(self, datadescription):
    """this will be called once per callback (before WriteImages) to allow the
       operation to export any desired data which is not an image. The child
       class should override this method if it wants so do something.
       For PhactoriIntersectNodeNormalsWithSurface we will output information
       about the nearest and furthest intersections. In this case we dump
       a .vtm file."""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriVtmDataExporterOperation." \
          "ExportOperationData entered\n", 100)
    self.mOutputFileBasename = self.mPhactoriOperationBlockOwner.mName
    global gPipeAndViewsState
    outfilename = self.mOutputFileBasename + "_" + \
        str(gPipeAndViewsState.mFrameTagCounter).zfill(4) + ".vtm"
    self.mFilterToWriteDataFrom.UpdatePipeline()

    #icsdClassname = inInputCsData.GetClassName()
    #if icsdClassname == "vtkMultiBlockDataSet" or \
    #   icsdClassname == "vtkExodusIIMultiBlockDataSet":


    if PhactoriDbg(100):
      myDebugPrint3("outfilename: " + outfilename + "\n" \
      "self.mFilterToWriteDataFrom: " + str(self.mFilterToWriteDataFrom) + "\n")

    if self.mWriter == None:
      mypid = SmartGetLocalProcessId()

      import vtkParallelCorePython
      pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
      globalController = pm.GetGlobalController()
      #gLocalProcessId = globalController.GetLocalProcessId()
      numproc = globalController.GetNumberOfProcesses()

      if PhactoriDbg(100):
        myDebugPrint3("numproc " + str(numproc) + " mypid " + str(mypid) + "\n")
      self.mWriter = vtk.vtkXMLPMultiBlockDataWriter()
      self.mWriter.SetInputData(self.mFilterToWriteDataFrom.GetClientSideObject().GetOutputDataObject(0))
      self.mWriter.SetNumberOfPieces(numproc)
      self.mWriter.SetStartPiece(mypid)
      #self.mWriter.SetEndPiece(mypid)
      #self.mWriter.SetUseSubdirectory(True)

      #if mypid == 0:
      #  self.mWriter.SetWriteSummaryFile(1)
      #else:
      #  self.mWriter.SetWriteSummaryFile(0)

    self.mWriter.SetFileName(outfilename)
    self.mWriter.Write()

      #self.mWriter = CreateWriter(outfilename, self.mFilterToWriteDataFrom)
    #writer = XMLMultiBlockDataWriter(FileName = outfilename,
    #    Input = self.mFilterToWriteDataFrom)
    #self.mWriter.UpdatePipeline()
    #del writer
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriVtmDataExporterOperation." \
          "ExportOperationData returning\n", 100)

class PhactoriIntersectNodeNormalsWithSurface(PhactoriOperationSpecifics):
  """takes each node ND from surface SS and its corresponding normal vector
     NN, and casts a ray RR towards the second surface TT which must be a set
     of triangles and finds if and where RR intersects TT at point XX.  For
     now it creates a new nodal variable NormalRayDistance which is the
     length of RR between ND and XX.  We may add capability in the future
     to store the XX as well, although this can be presumably be recovered
     from the normal and the distance"""

  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    #filled in from parsing
    self.mSourceOperationName = None
    self.mTargetOperationName = None

    #filled in by finding operation id'd by self.mSourceOperationName
    self.mSourceOperation = None

    #filled in by finding operation id'd by self.mTargetOperationName
    self.mTargetOperation = None

    self.mCalculationFrameTag = -1

    self.mLocalTriangleSet = SharedTriangleSet()
    self.mReceiveTriangleSet = SharedTriangleSet()

    self.mPointRayLineSource = None
    self.mTriangleRayLineSource = None
    self.mGroupLineSource = None

    self.mCreateIntersectionSegmentGroup = False
    self.mIntersectionSegmentGroupMinLength = 0.0
    self.mPointNormalRayIntersectDistanceArray = None
    self.mCellNormalRayIntersectDistanceArray = None

    self.SegmentsFromPointRays = PhactoriSegmentGroup1()
    self.SegmentsFromTriangleRays = PhactoriSegmentGroup1()

  def ParseParametersFromJson(self, inJson):
    if "target operation" in inJson:
      self.mTargetOperationName = inJson["target operation"]
    elif "target_operation" in inJson:
      self.mTargetOperationName = inJson["target_operation"]
    else:
      myDebugPrint3AndException(
          "PhactoriIntersectNodeNormalsWithSurface::ParseParametersFromJson\n"
          "Error:  must have 'target operation' or 'target_operation' token\n")
    if "input" in inJson:
      self.mSourceOperationName = inJson["input"]
    else:
      #will use default input for source
      self.mSourceOperationName = None

    if "create intersection segment group" in inJson:
      self.mCreateIntersectionSegmentGroup = \
        inJson["create intersection segment group"]
    elif "create_intersection_segment_group" in inJson:
      self.mCreateIntersectionSegmentGroup = \
        inJson["create_intersection_segment_group"]

    isgmlkey1 = "intersection segment group minimum length"
    isgmlkey2 = "intersection_segment_group_minimum_length"
    if isgmlkey1 in inJson:
      self.mIntersectionSegmentGroupMinLength = inJson[isgmlkey1]
    elif isgmlkey2 in inJson:
      self.mIntersectionSegmentGroupMinLength = inJson[isgmlkey2]

  def CreateParaViewFilter(self, inInputFilter):
    """create the filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriIntersectNodeNormalsWithSurface." \
          "CreateParaViewFilter entered\n", 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    if self.mCreateIntersectionSegmentGroup == True:

      self.mPointRayLineSource = PVTrivialProducer()
      self.mPointRayLineSource.GetClientSideObject().SetOutput(
          self.SegmentsFromPointRays.mVtkPolyData)
      self.mTriangleRayLineSource = PVTrivialProducer()
      self.mTriangleRayLineSource.GetClientSideObject().SetOutput(
          self.SegmentsFromTriangleRays.mVtkPolyData)

      #newParaViewFilter = GroupDatasets(Input = [])
      #self.mGroupLineSource = newParaViewFilter

      self.mGroupLineSource = GroupDatasets(
          Input = [self.mPointRayLineSource, self.mTriangleRayLineSource])
      self.mTubeFilter = Tube(Input = self.mGroupLineSource)
      #self.mTubeFilter.NumberofSides = 8
      self.mTubeFilter.Radius = 0.01
      #self.mTubeFilter.VaryRadius = 'By Scalar'
      newParaViewFilter = self.mTubeFilter

    else:
      newParaViewFilter = PVTrivialProducer()
      passthru = inInputFilter.GetClientSideObject().GetOutputDataObject(0)

      numpts = passthru.GetNumberOfPoints()

      newvar = vtk.vtkDoubleArray()
      newvar.SetNumberOfComponents(1)
      newvar.SetName("PointNormalRayIntersectDistance")
      for ii in xrange(0, numpts):
        newvar.InsertNextValue(float(-1.0))
      passthru.GetPointData().AddArray(newvar)
      self.mPointNormalRayIntersectDistanceArray = \
        passthru.GetPointData().GetArray("PointNormalRayIntersectDistance")

      numcells = passthru.GetNumberOfCells()

      newvar2 = vtk.vtkDoubleArray()
      newvar2.SetNumberOfComponents(1)
      newvar2.SetName("CellNormalRayIntersectDistance")
      for ii in xrange(0, numcells):
        newvar2.InsertNextValue(float(-1.0))
      passthru.GetCellData().AddArray(newvar2)
      self.mCellNormalRayIntersectDistanceArray = \
        passthru.GetCellData().GetArray("CellNormalRayIntersectDistance")

      if PhactoriDbg(100):
        numpts = passthru.GetNumberOfPoints()
        myDebugPrint3("numpts: " + str(numpts) + "\n")
        numptarrays = passthru.GetPointData().GetNumberOfArrays()
        myDebugPrint3("numptarrays: " + str(numptarrays) + "\n")
        numcells = passthru.GetNumberOfCells()
        myDebugPrint3("numcells: " + str(numcells) + "\n")
        numcellarrays = passthru.GetCellData().GetNumberOfArrays()
        myDebugPrint3("numcellarrays: " + str(numcellarrays) + "\n")

      newParaViewFilter.GetClientSideObject().SetOutput(passthru)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3(str(self.mGroupLineSource))
      myDebugPrint3("PhactoriIntersectNodeNormalsWithSurface." \
          "CreateParaViewFilter returning\n", 100)

    return newParaViewFilter

  def RunCalculationToCastRays(self, ioPipeAndViewsState):
    """our technique is as follows:  we loop through the pids of all the
       processes and each broadcasts the set triangles TTPID that pid holds
       from the target source to all the other processes.  Then each process
       casts rays from the nodes on the source surface that it holds (along
       the node normals) towards TTPID and finds the intersection if any.
       If a closer intersection is found than the existing one, the closer is
       used.  At the end of this loop the distance is known for each node on
       its own process (and if there was no intersection, which can be NAN or
       some arbitrary value).
       Initial implementation doesn't worry about self-intersection, but we
       can add additional features to expect certain number of 0 distance
       hits and ignore those or to remove self node and surrounding tris from
       the list of those to intersect"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriIntersectNodeNormalsWithSurface." \
          "RunCalculationToCastRays entered\n", 100)
    import time
    tmstrt = time.clock()

    if self.mCalculationFrameTag >= ioPipeAndViewsState.mFrameTagCounter:
      if PhactoriDbg(100):
        myDebugPrint3(\
          "self.mCalculationFrameTag indicates we have already run calc\n")
      return
    self.mCalculationFrameTag = ioPipeAndViewsState.mFrameTagCounter

    #get the source and target operations if necessary
    if(self.mTargetOperation == None):
      self.mTargetOperation = ioPipeAndViewsState.mOperationBlocks[
          self.mTargetOperationName]
      if(self.mTargetOperation == None):
        myDebugPrint3AndException("PhactoriIntersectNodeNormalsWithSurface:" \
            "RunCalculationToCastRays\ncouldn't find target " \
            "operation named " + str(self.mTargetOperationName) + "\n")
    if self.mSourceOperation == None:
      if self.mSourceOperationName == None:
        self.mSourceOperation = ioPipeAndViewsState.mIncomingDefaultOperation
      else:
        self.mSourceOperation = ioPipeAndViewsState.mOperationBlocks[
            self.mSourceOperationName]
      if self.mSourceOperation == None:
        myDebugPrint3AndException("PhactoriIntersectNodeNormalsWithSurface:" \
            "RunCalculationToCastRays\ncouldn't find source " \
            "operation named " + str(self.mSourceOperationName) + "\n")

    #first figure out the set of triangles from the target surface are on the
    #local process
    self.mLocalTriangleSet.CreateFromLocalProcess(self.mTargetOperation)

    #next figure out the set of nodes from the source surface are on the local
    #process

    #set up to hold the results

    #set up for line segments, if we are doing them
    if self.mCreateIntersectionSegmentGroup:
      csdata = self.mSourceOperation.mParaViewFilter.GetClientSideObject().\
        GetOutputDataObject(0)
      numSourcePoints = csdata.GetNumberOfPoints()
      self.SegmentsFromPointRays.SetupWithNumberOfItems(numSourcePoints)
      #self.SegmentsFromTriangleRays.SetupWithNumberOfItems(
      #  self.mLocalTriangleSet.numTriangles)

    #now loop through all processes; one process at a time share the triangles
    #from each process to the others and intersect with the local surface
    #portion

    mypid = SmartGetLocalProcessId()

    import vtkParallelCorePython
    pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
    globalController = pm.GetGlobalController()
    #gLocalProcessId = globalController.GetLocalProcessId()
    numproc = globalController.GetNumberOfProcesses()
    if PhactoriDbg(100):
      myDebugPrint3(
        "mypid: " + str(mypid) + "  numproc: " + str(numproc) + "\n")

    import vtkParallelCorePython
    pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
    globalController = pm.GetGlobalController()
    for ii in range(0,numproc):
      if ii == mypid:
        self.mLocalTriangleSet.SendBroadcast(ii, globalController)
        #self.mLocalTriangleSet.IntersectWithNormalRaysFromSourceTriangles(
        #  self.mSourceOperation, self, self.mLocalTriangleSet)
        self.mLocalTriangleSet.IntersectWithNormalRaysFromSource(
          self.mSourceOperation, self)
      else:
        self.mReceiveTriangleSet.ReceiveBroadcast(ii, globalController)
        #self.mReceiveTriangleSet.IntersectWithNormalRaysFromSourceTriangles(
        #  self.mSourceOperation, self, self.mLocalTriangleSet)
        self.mReceiveTriangleSet.IntersectWithNormalRaysFromSource(
          self.mSourceOperation, self)

    tmend = time.clock()

    if PhactoriDbg(100):
      myDebugPrint3("time for RunCalculationToCastRays: " + \
        str(tmend - tmstrt) + "\n")
    else:
      print "time for RunCalculationToCastRays raycasting: " + str(tmend - tmstrt)

    #create/handle line segments, if we are doing them
    if self.mCreateIntersectionSegmentGroup:
      self.SegmentsFromPointRays.UpdateParaViewSegments(
        self.mIntersectionSegmentGroupMinLength)
      #self.SegmentsFromTriangleRays.UpdateParaViewSegments()

    tmend2 = time.clock()
    if PhactoriDbg(100):
      myDebugPrint3(
        "time for RunCalculationToCastRays line segment making: " + \
        str(tmend2 - tmend) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriIntersectNodeNormalsWithSurface." \
          "RunCalculationToCastRays returning\n", 100)

  def ExportOperationData(self, datadescription):
    """this will be called once per callback (before WriteImages) to allow the
       operation to export any desired data which is not an image. The child
       class should override this method if it wants so do something.
       For PhactoriIntersectNodeNormalsWithSurface we will output information
       about the nearest and furthest intersections"""

    if PhactoriDbg(100):
      myDebugPrint3(
        "PhactoriIntersectNodeNormalsWithSurface::ExportOperationData "
        "entered\n", 100)

    UpdatePipelineWithCurrentTimeArgument(
      self.mPhactoriOperationBlockOwner.mParaViewFilter)
    self.RunCalculationToCastRays(gPipeAndViewsState)
    UpdatePipelineWithCurrentTimeArgument(
      self.mPhactoriOperationBlockOwner.mParaViewFilter)

    if self.mPointNormalRayIntersectDistanceArray == None:
      if PhactoriDbg(100):
        myDebugPrint3("no data values, returning\n")
      return

    numpts = self.mPointNormalRayIntersectDistanceArray.GetNumberOfValues()
    maxval = -1.0
    minval = -1.0
    minval2 = -1.0
    if numpts > 0:
      vv = self.mPointNormalRayIntersectDistanceArray.GetValue(0)
      maxval = vv
      minval = vv
      if vv >= -0.5:
        minval2 = vv
      for ii in xrange(1, numpts):
        vv = self.mPointNormalRayIntersectDistanceArray.GetValue(ii)
        if vv > maxval:
          maxval = vv
        if vv < minval:
          minval = vv
        if vv >= -0.5:
          if (minval2 < -0.5) or (vv < minval2):
            minval2 = vv

    fname = "RayIntersectionMinMax_" + \
      self.mPhactoriOperationBlockOwner.mName + "_process_" + \
      str(SmartGetLocalProcessId()) + ".txt"
    try:
      ff = open(fname, "a+b")
      ff.write("minimum, maximum, nonnegative minimum\n")
      ff.write(str(minval) + ", " + str(maxval) + ", " + str(minval2) + "\n")
    except:
      myDebugPrint3AndException(
          "PhactoriIntersectNodeNormalsWithSurface::ExportOperationData\n"
          "Error writing file: " + str(fname) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3(
        "PhactoriIntersectNodeNormalsWithSurface::ExportOperationData "
        "returning\n", 100)
    return


class PhactoriNearestPointsOperation(PhactoriOperationSpecifics):
  """finds the set of N nearest nodes between a source operation/input and
     a target operation input.  Also outputs a table of the set of nodes
     and distances at render time.  This operation as currently constructed
     pulls all the target nodes into each processor and thus may cause memory
     issues.  It may be advisible to use, e.g. quadrature decimation to
     reduce the number of points in the target and/or source first.  We may
     eventually make this a multi-stage operation to first operate on
     a reduced set of points to find a subset of the original sets of points
     to then do an exact comparison on.  We may also eventually change to
     operated based on a space partitioning tree to avoid checking every
     point against every point"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    #filled in from parsing
    self.mSourceOperationName = None
    self.mTargetOperationName = None

    #filled in by finding operation id'd by self.mSourceOperationName
    self.mSourceOperation = None

    #filled in by finding operation id'd by self.mTargetOperationName
    self.mTargetOperation = None

    #controls how many points we look for in each process, and then how many
    #points out of those we select as nearest
    #if self.mNumPointsToFind > self.mNumPointsPerProcess it is possible to
    #have incorrect points beyond self.mNumPointsPerProcess when all the
    #nearest points happen to reside on one process
    self.mNumPointsPerProcess = 5
    self.mNumPointsToFind = 5

    #class to maintain parallel shared information
    self.mClosestNPointsFromEachProcess = None

    self.mLineSource = None

    self.mCalculationFrameTag = -1

    self.mWriteCsvResults = True
    self.mCsvFileBasedirectory = gPipeAndViewsState.mDefaultBasedirectory
    self.mCsvFileName = None

    #if not None, this will be a string which identifies a text annotation
    #which the user has to create normally.  If said annotation exists, this
    #operation will overwrite the text for the annotation with information
    #on the nearest point
    self.mOutputToTextAnnotationName = None
    self.mTextAnnotationOutputTarget = None

    self.mTxtAntnFormatString = None
    self.mCsvFormatString = None

  def ParseParametersFromJson(self, inJson):

    self.mTxtAntnFormatString = getParameterFromBlock(inJson,
        'text annotation float format string', self.mTxtAntnFormatString)
    self.mTxtAntnFormatString = getParameterFromBlock(inJson,
        'text_annotation_float_format_string', self.mTxtAntnFormatString)
    self.mCsvFormatString = getParameterFromBlock(inJson,
        'csv float format string', self.mCsvFormatString)
    self.mCsvFormatString = getParameterFromBlock(inJson,
        'csv_float_format_string', self.mCsvFormatString)
    self.mOutputToTextAnnotationName = getParameterFromBlock(inJson,
        'output to text annnotation', self.mOutputToTextAnnotationName)
    self.mOutputToTextAnnotationName = getParameterFromBlock(inJson,
        'output_to_text_annnotation', self.mOutputToTextAnnotationName)
    self.mWriteCsvResults = getParameterFromBlock(inJson,
        'write_results_csv_file', self.mWriteCsvResults)
    self.mWriteCsvResults = getParameterFromBlock(inJson,
        'write results csv file', self.mWriteCsvResults)
    self.mCsvFileName = getParameterFromBlock(inJson,
        'csv_results_filename', self.mCsvFileName)
    self.mCsvFileName = getParameterFromBlock(inJson,
        'csv results filename', self.mCsvFileName)
    self.mCsvFileBasedirectory = getParameterFromBlock(inJson,
        'csv_results_basedirectory', self.mCsvFileBasedirectory)
    self.mCsvFileBasedirectory = getParameterFromBlock(inJson,
        'csv results basedirectory', self.mCsvFileBasedirectory)
    if self.mCsvFileName == None:
      #make default name including operation name and datetime, without
      #microseconds
      import datetime
      nn = datetime.datetime.now()
      self.mCsvFileName = self.mPhactoriOperationBlockOwner.mName + \
          ".nearestpoints." +  (nn.isoformat()[0:-7]) + ".csv"
    if 'target operation' in inJson:
      self.mTargetOperationName = inJson['target operation']
    elif 'target_operation' in inJson:
      self.mTargetOperationName = inJson['target_operation']
    else:
      myDebugPrint3AndException(
          "PhactoriNearestPointsOperation::ParseParametersFromJson\n"
          "Error:  must have 'target operation' or 'target_operation' token\n")
    if 'input' in inJson:
      self.mSourceOperationName = inJson['input']
    else:
      #will use default input for source
      self.mSourceOperationName = None

    if 'number of points per process' in inJson:
      self.mNumPointsPerProcess = inJson['number of points per process']
    elif 'number of points per process' in inJson:
      self.mNumPointsPerProcess = inJson['number_of_points_per_process']
    if 'number of points to find' in inJson:
      self.mNumPointsToFind = inJson['number of points to find']
      #unless user purposely set self.mNumPointsPerProcess different, we want
      #self.mNumPointsPerProcess == self.mNumPointsToFind
      if 'number of points per process' not in inJson:
        self.mNumPointsPerProcess = self.mNumPointsToFind
    elif 'number_of_points_to_find' in inJson:
      self.mNumPointsToFind = inJson['number_of_points_to_find']
      #unless user purposely set self.mNumPointsPerProcess different, we want
      #self.mNumPointsPerProcess == self.mNumPointsToFind
      if 'number_of_points per_process' not in inJson:
        self.mNumPointsPerProcess = self.mNumPointsToFind

  def CreateParaViewFilter(self, inInputFilter):
    """create the group filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriNearestPointsOperation.CreateParaViewFilter "
          "entered\n", 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    #newParaViewFilter = Sphere()
    newParaViewFilter = Line()
    newParaViewFilter.Point1 = [0.0, 0.0, 0.0]
    newParaViewFilter.Point2 = [0.01, 0.01, 0.01]

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriNearestPointsOperation.CreateParaViewFilter "
          "returning\n", 100)

    self.mLineSource = newParaViewFilter
    return newParaViewFilter

  def RunCalculationToFindNearestPoints(self, ioPipeAndViewsState):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriNearestPointsOperation." \
          "RunCalculationToFindNearestPoints entered\n", 100)

    if self.mCalculationFrameTag >= ioPipeAndViewsState.mFrameTagCounter:
      if PhactoriDbg(100):
        myDebugPrint3(\
          "self.mCalculationFrameTag indicates we have already run calc\n")
      return
    self.mCalculationFrameTag = ioPipeAndViewsState.mFrameTagCounter

    #get the source and target operations if necessary
    if(self.mTargetOperation == None):
      self.mTargetOperation = ioPipeAndViewsState.mOperationBlocks[
          self.mTargetOperationName]
      if(self.mTargetOperation == None):
        myDebugPrint3AndException("PhactoriNearestPointsOperation:" \
            "RunCalculationToFindNearestPoints\ncouldn't find target " \
            "operation named " + str(self.mTargetOperationName) + "\n")
    if self.mSourceOperation == None:
      if self.mSourceOperationName == None:
        self.mSourceOperation = ioPipeAndViewsState.mIncomingDefaultOperation
      else:
        self.mSourceOperation = ioPipeAndViewsState.mOperationBlocks[
            self.mSourceOperationName]
      if self.mSourceOperation == None:
        myDebugPrint3AndException("PhactoriNearestPointsOperation:" \
            "RunCalculationToFindNearestPoints\ncouldn't find source " \
            "operation named " + str(self.mSourceOperationName) + "\n")


    #get the list of points from the target that reside on this process
    nodeIds, pointXyzs = self.mTargetOperation.MakeListOfAllPoints1()

    if PhactoriDbg(100):
      myDebugPrint3("In this process " + str(SmartGetLocalProcessId()) + \
          " the target has " + str(nodeIds.GetNumberOfValues()) + " points\n")
      #myDebugPrint3("nodeIds in " + self.mName + \
      #    " (count " + str(len(nodeIds)) + "):\n" + str(nodeIds) + "\n")
      #for jj in range(0, nodeIds.GetNumberOfValues()):
      #  myDebugPrint3(str(jj) + ": " + str(nodeIds.GetValue(jj)) + "\n")
      #myDebugPrint3("points in " + self.mName + ":\n" + str(pointXyzs) + "\n")

    #get the list of points from the target on all processes (copied to this process)
    allProcTgtNodeIds, allProcTgtXyzs = \
        UseMPIToCreateSharedPointList(nodeIds, pointXyzs)

    if PhactoriDbg(100):
      myDebugPrint3("from all processes the target has " + \
          str(allProcTgtNodeIds.GetNumberOfValues()) + " points\n")

    #get nearest N points from the source points on this process to all the
    #points on the target (which we have from all processes)
    nearestset = self.mSourceOperation.FindClosestNPointsToList(
        allProcTgtNodeIds, allProcTgtXyzs, self.mNumPointsPerProcess)

    if PhactoriDbg(100):
        myDebugPrint3("nearest points from local process source points to " \
            "all target points:\n")
        myDebugPrint3(nearestset.ToStr())

    if self.mClosestNPointsFromEachProcess == None:
      self.mClosestNPointsFromEachProcess = ClosestNPointsFromEachProcess()

    #now share those N points from each process to all processes so for this and
    #every other process we have a list of N * numprocesses points which
    #represent the closestst N points from each process to the target
    self.mClosestNPointsFromEachProcess.\
        ParallelSetupFromFindClosestNPointsToList(nearestset)

    #now find the closest M points out of that list
    self.mClosestNPointsFromEachProcess.FindNClosestPointsNdxInList(
        self.mNumPointsToFind)

    if PhactoriDbg(100):
      myDebugPrint3( self.mClosestNPointsFromEachProcess.\
          MakeNClosestPointsNdxInListTableString())

    ndx = self.mClosestNPointsFromEachProcess.mClosestNRefs[0][1]
    srcndid = self.mClosestNPointsFromEachProcess.mThisProcIds.GetValue(ndx)
    tgtndid = self.mClosestNPointsFromEachProcess.mTargetMatchIds.GetValue(ndx)
    distsqrd = self.mClosestNPointsFromEachProcess.mDistSqrds.GetValue(ndx)
    dist = math.sqrt(distsqrd)
    pndx = ndx*3
    srcx = self.mClosestNPointsFromEachProcess.mThisProcXyzs.GetValue(pndx)
    srcy = self.mClosestNPointsFromEachProcess.mThisProcXyzs.GetValue(pndx+1)
    srcz = self.mClosestNPointsFromEachProcess.mThisProcXyzs.GetValue(pndx+2)
    tgtx = self.mClosestNPointsFromEachProcess.mTargetMatchXyzs.GetValue(pndx)
    tgty = self.mClosestNPointsFromEachProcess.mTargetMatchXyzs.GetValue(pndx+1)
    tgtz = self.mClosestNPointsFromEachProcess.mTargetMatchXyzs.GetValue(pndx+2)

    self.mLineSource.Point1 = [srcx, srcy, srcz]
    self.mLineSource.Point2 = [tgtx, tgty, tgtz]
    self.mLineSource.UpdatePipeline()

    self.WriteClosestNPointsToCsvFile(ioPipeAndViewsState)
    self.OutputToTextAnnotation(ioPipeAndViewsState)

  def OutputToTextAnnotation(self, ioPipeAndViewsState):
    if PhactoriDbg():
      myDebugPrint3("PhactoriNearestPointsOperation::OutputToTextAnnotation" \
        " called\n")
    if self.mOutputToTextAnnotationName == None:
      return

    if self.mTextAnnotationOutputTarget == None:
      self.mTextAnnotationOutputTarget = \
        ioPipeAndViewsState.mTextAnnotationBlocks[ \
          self.mOutputToTextAnnotationName]
      if self.mTextAnnotationOutputTarget == None:
        if PhactoriDbg():
          myDebugPrint3(
            "PhactoriNearestPointsOperation::OutputToTextAnnotation\n" \
            "could not find PhactoriTextAnnotationBlock named:\n" \
            + str(self.mOutputToTextAnnotationName) + "\n")
        return

    if PhactoriDbg():
      myDebugPrint3("found annotation " + \
        str(self.mOutputToTextAnnotationName) + "\n")
    clstpts = self.mClosestNPointsFromEachProcess
    ndx = clstpts.mClosestNRefs[0][1]
    srcndid = clstpts.mThisProcIds.GetValue(ndx)
    tgtndid = clstpts.mTargetMatchIds.GetValue(ndx)
    distsqrd = clstpts.mDistSqrds.GetValue(ndx)
    dist = math.sqrt(distsqrd)
    pndx = ndx*3
    srcx = clstpts.mThisProcXyzs.GetValue(pndx)
    srcy = clstpts.mThisProcXyzs.GetValue(pndx+1)
    srcz = clstpts.mThisProcXyzs.GetValue(pndx+2)
    tgtx = clstpts.mTargetMatchXyzs.GetValue(pndx)
    tgty = clstpts.mTargetMatchXyzs.GetValue(pndx+1)
    tgtz = clstpts.mTargetMatchXyzs.GetValue(pndx+2)

    if self.mTxtAntnFormatString == None:
      self.mTextAnnotationOutputTarget.mTextString = \
        "nearest: " + str(dist) + "\n" \
        "node id 1: " + str(srcndid) + "\n" \
        "x1: " + str(srcx) + "\n" \
        "y1: " + str(srcy) + "\n" \
        "z1: " + str(srcz) + "\n" \
        "node id 2: " + str(tgtndid) + "\n" \
        "x2: " + str(tgtx) + "\n" \
        "y2: " + str(tgty) + "\n" \
        "z2: " + str(tgtz)
    else:
      diststr = self.mTxtAntnFormatString%dist
      id1str = str(srcndid)
      x1str = self.mTxtAntnFormatString%srcx
      y1str = self.mTxtAntnFormatString%srcy
      z1str = self.mTxtAntnFormatString%srcz
      id2str = str(tgtndid)
      x2str = self.mTxtAntnFormatString%tgtx
      y2str = self.mTxtAntnFormatString%tgty
      z2str = self.mTxtAntnFormatString%tgtz
      self.mTextAnnotationOutputTarget.mTextString = \
        "nearest: " + diststr + "\n" \
        "node id 1: " + id1str + "\n" \
        "x1: " + x1str + "\n" \
        "y1: " + y1str + "\n" \
        "z1: " + z1str + "\n" \
        "node id 2: " + id2str + "\n" \
        "x2: " + x2str + "\n" \
        "y2: " + y2str + "\n" \
        "z2: " + z2str

    self.mTextAnnotationOutputTarget.mParaViewSource.Text = \
      self.mTextAnnotationOutputTarget.mTextString
    self.mTextAnnotationOutputTarget.mParaViewSource.UpdatePipeline()

  def WriteClosestNPointsToCsvFile(self, ioPipeAndViewsState):
    if SmartGetLocalProcessId() != 0:
      return
    if self.mWriteCsvResults == False:
      return

    if PhactoriDbg():
      myDebugPrint3("WriteClosestNPointsToCsvFile entered and writing\n")

    if self.mCsvFileBasedirectory == None:
      ffname = self.mCsvFileName
    elif self.mCsvFileBasedirectory == "":
      ffname = self.mCsvFileName
    else:
      import os
      ffname = self.mCsvFileBasedirectory + os.sep + self.mCsvFileName
    try:
      cvsff = open(ffname, "a+b")
      if self.mCalculationFrameTag == 1:
        cvsff.write("call index, simulation time, position, " \
          "distance, distance squared, " \
          "source node id, source x, source y, source z, " \
          "target node id, target x, target y, target z\n")
      clstpts = self.mClosestNPointsFromEachProcess
      try:
        simtime = ioPipeAndViewsState.CurrentDatadescription.GetTime()
      except:
        simtime = 0.0
      for ii in range(self.mNumPointsToFind):
        if PhactoriDbg():
          myDebugPrint3("doing point\n" + str(ii) + "\n")
        ndx = clstpts.mClosestNRefs[ii][1]
        srcndid = clstpts.mThisProcIds.GetValue(ndx)
        tgtndid = clstpts.mTargetMatchIds.GetValue(ndx)
        distsqrd = clstpts.mDistSqrds.GetValue(ndx)
        dist = math.sqrt(distsqrd)
        pndx = ndx*3
        srcx = clstpts.mThisProcXyzs.GetValue(pndx)
        srcy = clstpts.mThisProcXyzs.GetValue(pndx+1)
        srcz = clstpts.mThisProcXyzs.GetValue(pndx+2)
        tgtx = clstpts.mTargetMatchXyzs.GetValue(pndx)
        tgty = clstpts.mTargetMatchXyzs.GetValue(pndx+1)
        tgtz = clstpts.mTargetMatchXyzs.GetValue(pndx+2)
        if self.mCsvFormatString == None:
          cvsff.write(
            str(self.mCalculationFrameTag) + ", " + \
            str(simtime) + ", " + \
            str(ii) + ", " + \
            str(dist) + ", " + \
            str(distsqrd) + ", " + \
            str(srcndid) + ", " + \
            str(srcx) + ", " + \
            str(srcy) + ", " + \
            str(srcz) + ", " + \
            str(tgtndid) + ", " + \
            str(tgtx) + ", " + \
            str(tgty) + ", " + \
            str(tgtz) + "\n"
          )
        else:
          tagstr = str(self.mCalculationFrameTag)
          simtimestr = self.mCsvFormatString%simtime
          posstr = str(ii)
          diststr = self.mCsvFormatString%dist
          distsqrdstr = self.mCsvFormatString%distsqrd
          id1str = str(srcndid)
          x1str = self.mCsvFormatString%srcx
          y1str = self.mCsvFormatString%srcy
          z1str = self.mCsvFormatString%srcz
          id2str = str(tgtndid)
          x2str = self.mCsvFormatString%tgtx
          y2str = self.mCsvFormatString%tgty
          z2str = self.mCsvFormatString%tgtz
          cvsff.write(
            tagstr + ", " + \
            simtimestr + ", " + \
            posstr + ", " + \
            diststr + ", " + \
            distsqrdstr + ", " + \
            id1str + ", " + \
            x1str + ", " + \
            y1str + ", " + \
            z1str + ", " + \
            id2str + ", " + \
            x2str + ", " + \
            y2str + ", " + \
            z2str + "\n"
          )
      cvsff.close()
    except:
      myDebugPrint3AndException(
        "PhactoriNearestPointsOperation::" \
        "RunCalculationToFindNearestPoints\n" \
        "error opening or writing to .csv file:\n" \
        + str(ffname) + "\n")

    if PhactoriDbg():
      myDebugPrint3("WriteClosestNPointsToCsvFile returning\n")



class PhactoriGroupOperation(PhactoriOperationSpecifics):
  """manages group filter: group the results of two or more operations into
     a single multiblock source"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    #filled in from parsing
    self.mOperationNameList = None
    #filled in by finding operations id'd by self.mOperationNameList
    self.mOperationList = None

  def ParseParametersFromJson(self, inJson):

    if 'operation group list' in inJson:
      self.mOperationNameList = inJson['operation group list']
    else:
      myDebugPrint3AndException(
          "PhactoriGroupOperation::ParseParametersFromJson\n"
          "Error:  must have 'operation group list' token\n")

  def CreateParaViewFilter2(self, ioPipeAndViewsState):
    """create the group filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGroupOperation.CreateParaViewFilter "
          "entered\n", 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    #get direct list of operations (not through name list)
    #should be valid here
    self.mOperationList = []
    for ii in self.mOperationNameList:
      if ii not in ioPipeAndViewsState.mOperationBlocks:
        myDebugPrint3AndException(
          "PhactoriGroupOperation::CreateParaViewFilter\n"
          "Error:  operation '" + str(ii) + "' not available\n")
      inputOperationBlock = ioPipeAndViewsState.mOperationBlocks[ii]
      if inputOperationBlock.GetPvFilter() == None:
        myDebugPrint3AndException(
          "PhactoriGroupOperation::CreateParaViewFilter\n"
          "Error:  operation '" + str(ii) + "' paraview filter not "
          "constructed\n")
      self.mOperationList.append(inputOperationBlock.GetPvFilter())

    newParaViewFilter = GroupDatasets(Input = self.mOperationList)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGroupOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

class PhactoriGenerateSurfaceNormalsOperation(PhactoriOperationSpecifics):
  """manages ExtractSurface filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.ComputeCellNormals = True

  def ParseParametersFromJson(self, inJson):
    if "compute cell normals" in inJson:
      self.ComputeCellNormals = inJson["compute cell normals"]
    elif "compute_cell_normals" in inJson:
      self.ComputeCellNormals = inJson["compute_cell_normals"]

  def CreateParaViewFilter(self, inInputFilter):
    """create the GenerateSurfaceNormals filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGenerateSurfaceNormalsOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    newParaViewFilter = GenerateSurfaceNormals(inInputFilter)
    if self.ComputeCellNormals:
      newParaViewFilter.ComputeCellNormals = 1
    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      pvsDi = newParaViewFilter.GetDataInformation()
      numCells = pvsDi.GetNumberOfCells()
      numPoints = pvsDi.GetNumberOfPoints()
      myDebugPrint3("numCells, numPoints  " + str(numCells) + " " + str(numPoints) + "\n")
      myDebugPrint3("newParaViewFilter: " + str(newParaViewFilter) + "\n")


    if PhactoriDbg(100):
      myDebugPrint3("PhactoriGenerateSurfaceNormalsOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter


class PhactoriExtractSurfaceOperation(PhactoriOperationSpecifics):
  """manages ExtractSurface filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)

  def ParseParametersFromJson(self, inJson):
    dummy = 0

  def CreateParaViewFilter(self, inInputFilter):
    """create the ExtractSurface filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractSurfaceOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    newParaViewFilter = ExtractSurface(inInputFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractSurfaceOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter


class PhactoriGhostCellsGeneratorOperation(PhactoriOperationSpecifics):
  """manages GhostCellsGenerator filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)

  def ParseParametersFromJson(self, inJson):
    dummy = 0

  def CreateParaViewFilter(self, inInputFilter):
    """create the GhostCellsGenerator filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3(
          "PhactoriGhostCellsGeneratorOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    newParaViewFilter = GhostCellsGenerator(inInputFilter)
    #newParaViewFilter = D3(inInputFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3(
          "PhactoriGhostCellsGeneratorOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter


class PhactoriMergeBlocksOperation(PhactoriOperationSpecifics):
  """manages MergeBlocks filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)

  def ParseParametersFromJson(self, inJson):
    dummy = 0

  def CreateParaViewFilter(self, inInputFilter):
    """create the MergeBlocks filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMergeBlocksOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    newParaViewFilter = MergeBlocks(inInputFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMergeBlocksOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter


class PhactoriSubdivideOperation(PhactoriOperationSpecifics):
  """manages Subdivide filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.NumberOfSubdivisions = 1

  def ParseParametersFromJson(self, inJson):
    if "number of subdivisions" in inJson:
      self.NumberOfSubdivisions = inJson["number of subdivisions"]
    elif "number_of_subdivisions" in inJson:
      self.NumberOfSubdivisions = inJson["number_of_subdivisions"]

  def CreateParaViewFilter(self, inInputFilter):
    """create the Triangulate filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSubdivideOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    newParaViewFilter = Subdivide(inInputFilter)
    if self.NumberOfSubdivisions != 1:
      newParaViewFilter.NumberofSubdivisions = self.NumberOfSubdivisions

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriSubdivideOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter


class PhactoriTriangulateOperation(PhactoriOperationSpecifics):
  """manages Triangulate filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)

  def ParseParametersFromJson(self, inJson):
    dummy = 0

  def CreateParaViewFilter(self, inInputFilter):
    """create the Triangulate filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriTriangulateOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    newParaViewFilter = Triangulate(inInputFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriTriangulateOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter


class PhactoriExtractBlockOperation(PhactoriOperationSpecifics):
  """manages extract block filter, including creating flat block
     indices list from list of block names"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    #mIncludeblockList OR mExcludeBlockList will be filled in from parsing
    self.mIncludeBlockList = None
    self.mExcludeBlockList = None

    #this will be list of included block indices, and is calculated from
    #mIncludeBlockList / mExcludeBlockList and passed directly to
    #ExcludeBlockFilter.BlockIndices
    self.mBlockIndices = []

  def ParseParametersFromJson(self, inJson):

    if 'include blocks' in inJson:
      self.mIncludeBlockList = inJson['include blocks']

    if 'exclude blocks' in inJson:
      self.mExcludeBlockList = inJson['exclude blocks']

    if self.mIncludeBlockList == None and self.mExcludeBlockList == None:
      myDebugPrint3AndException(
          "PhactoriExtractBlockOperation::ParseParametersFromJson\n"
          "Error:  must have include block list or exclude block list\n")

  def FigureBlockIndicesFromBlockListOneBlock(self, inMetaData,
          ioFlatIndexCounter):
    """determine if this one block should have it's flat index tripped on for
       the extract block filter (leaf item of recursion)"""

    if inMetaData == None:
      thisBlockName = None
    else:
      thisBlockName = inMetaData.Get(vtk.vtkCompositeDataSet.NAME())

    if self.mIncludeBlockList != None:
      if thisBlockName == None:
        if PhactoriDbg(100):
          myDebugPrint3("block with no name " + \
              " not in include list, not + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
      elif thisBlockName in self.mIncludeBlockList:
        self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
        if PhactoriDbg(100):
          myDebugPrint3("block " + str(thisBlockName) + \
              " in include list, + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
      else:
        if PhactoriDbg(100):
          myDebugPrint3("block " + str(thisBlockName) + \
              " not in include list, not + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
    elif self.mExcludeBlockList != None:
      if thisBlockName == None:
        self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
        if PhactoriDbg(100):
          myDebugPrint3("block with no name " + \
              " not in exclude list, + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
      elif thisBlockName not in self.mExcludeBlockList:
        self.mBlockIndices.append(int(ioFlatIndexCounter[0]) - 1)
        if PhactoriDbg(100):
          myDebugPrint3("block " + str(thisBlockName) + \
              " not in exclude list, + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
      else:
        if PhactoriDbg(100):
          myDebugPrint3("block " + str(thisBlockName) + \
              " in exclude list, not + to mBlockIndices (flat index " + \
              str(ioFlatIndexCounter[0] - 1) + ")\n")
    else:
      myDebugPrint3AndException(
          "PhactoriExtractBlockOperation::"
          "FigureBlockIndicesFromBlockListOneBlock\n"
          "Error:  must have include block list or exclude block list\n")

  def FigureBlockIndicesFromBlockListRecurse1(self, inCsdata, inMetaData,
          ioFlatIndexCounter):
    """recursively go through multiblock dataset to determine flat indices
       of blocks to be included; we are assuming only leaf blocks are
       actually named"""

    ioFlatIndexCounter[0] += 1L

    icsdClassname = inCsdata.GetClassName()
    if icsdClassname == "vtkMultiBlockDataSet" or \
       icsdClassname == "vtkExodusIIMultiBlockDataSet":
      numBlocks = inCsdata.GetNumberOfBlocks()
      if PhactoriDbg(100):
        myDebugPrint3("recursing, flat index: " + \
            str(ioFlatIndexCounter[0]) + \
            "    num blocks: " + \
            str(numBlocks) + "\n")
      for ii in range(0, numBlocks):
        oneBlock = inCsdata.GetBlock(ii)
        oneBlockMetaData = inCsdata.GetMetaData(ii)
        if PhactoriDbg(100):
          myDebugPrint3("oneBlockMetaData: " + str(oneBlockMetaData) + "\n")
          #myDebugPrint3("name: " + \
          #    oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME()) + "\n")
        if oneBlock != None:
          if oneBlockMetaData != None:
            theBlockName = oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME())
            if PhactoriDbg(100):
              myDebugPrint3("name for block " + str(ii) + ":  " + \
                str(theBlockName) + "\n")
          else:
            if PhactoriDbg(100):
              myDebugPrint3("block " + str(ii) + " meta data was None\n")
          self.FigureBlockIndicesFromBlockListRecurse1(oneBlock,
              oneBlockMetaData, ioFlatIndexCounter)
        else:
          #I think we need to count here to be okay with pruned stuff; maybe
          #we need to set extract block to no pruning (?)
          ioFlatIndexCounter[0] += 1L
    else:
      self.FigureBlockIndicesFromBlockListOneBlock(inMetaData,
          ioFlatIndexCounter)
    if icsdClassname == "vtkMultiPieceDataSet":
      numpieces = inCsdata.GetNumberOfPieces()
      # ioFlatIndexCounter[0] += numpieces - 1
      ioFlatIndexCounter[0] += numpieces



  def FigureBlockIndicesFromBlockList(self, inInputFilter):
    """from the list of include/exclude blocks create a indices list
       of blocks to put in filter"""
    if PhactoriDbg(100):
      myDebugPrint3("FigureBlockIndicesFromBlockList entered\n"
          "self.mIncludeBlockList: " + str(self.mIncludeBlockList) + "\n")

    csdata = inInputFilter.GetClientSideObject().GetOutputDataObject(0)
    flatIndexCounter = [0L]
    self.FigureBlockIndicesFromBlockListRecurse1(csdata, None,
        flatIndexCounter)
    if PhactoriDbg(100):
      myDebugPrint3("FigureBlockIndicesFromBlockList final indices list:\n" \
          + str(self.mBlockIndices) + \
          "\nFigureBlockIndicesFromBlockList returning\n")

  def CreateParaViewFilter(self, inInputFilter):
    """create the extract block filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractBlockOperation.CreateParaViewFilter "
          "entered\n", 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    newParaViewFilter = ExtractBlock(inInputFilter)
    self.FigureBlockIndicesFromBlockList(inInputFilter)

    #newParaViewFilter.PruneOutput = 1
    #newParaViewFilter.MaintainStructure = 0
    newParaViewFilter.MaintainStructure = 1

    newParaViewFilter.BlockIndices = self.mBlockIndices

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractBlockOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter


class PhactoriCellDataToPointDataOperation(PhactoriOperationSpecifics):
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    return

  def ParseParametersFromJson(self, inJson):
    return

  def CreateParaViewFilter(self, inInputFilter):
    """create the clip plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriCellDataToPointDataOperation:CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()
    newParaViewFilter = CellDatatoPointData(Input = inInputFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriCellDataToPointDataOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

class PhactoriContourOperation(PhactoriOperationSpecifics):
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

class PhactoriPlaneOpBase(PhactoriOperationSpecifics):
  """base operation type for slice operation and clip operation.  They share
     parsing and settings, and just differ by which filter is created"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mPtOnPlaneInfo = PhactoriUserPointInfo()
    self.mPlaneNormal = [0.0, 1.0, 0.0]
    self.mUseThreePointsOnPlane = False
    self.mPtOnPlaneInfoA = None
    self.mPtOnPlaneInfoB = None
    self.mPtOnPlaneInfoC = None
    self.mCrinkleSetting = None
    if gParaViewCatalystVersionFlag < 50502:
      self.mInsideOut = 1
    else:
      self.mInvert = 1

  def ParseParametersFromJson(self, inJson):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPlaneOpBase.ParseParametersFromJson entered\n", 100)

    if "plane specification" in inJson:
        if inJson["plane specification"] == "point and normal":
          self.mUseThreePointsOnPlane = False
        elif inJson["plane specification"] == "three points":
          self.mUseThreePointsOnPlane = True
        else:
          myDebugPrint3AndException(
            "PhactoriPlaneOpBase.ParseParametersFromJson:\n"
            "'plane specification' must be 'point and normal' or "
            "'three points'\n")

    if self.mUseThreePointsOnPlane:
      if PhactoriDbg():
        myDebugPrint3("using three points on plane instead of plane/normal:\n")
      self.mPtOnPlaneInfoA = PhactoriUserPointInfo()
      self.mPtOnPlaneInfoA.mXyz = [0.0, 0.0, 0.0]
      self.mPtOnPlaneInfoB = PhactoriUserPointInfo()
      self.mPtOnPlaneInfoB.mXyz = [0.0, 1.0, 0.0]
      self.mPtOnPlaneInfoC = PhactoriUserPointInfo()
      self.mPtOnPlaneInfoC.mXyz = [0.0, 0.0, 1.0]
      self.mPtOnPlaneInfoA.UserPointInfoParseParametersFromJson(
              inJson, "", " on plane A")
      self.mPtOnPlaneInfoB.UserPointInfoParseParametersFromJson(
              inJson, "", " on plane B")
      self.mPtOnPlaneInfoC.UserPointInfoParseParametersFromJson(
              inJson, "", " on plane C")
    else:
      if PhactoriDbg():
        myDebugPrint3(
          "using point/normal: PhactoriPlaneOpBase parsing point onplane:\n")
      self.mPtOnPlaneInfo.UserPointInfoParseParametersFromJson(
              inJson, "", " on plane")

      if 'plane normal' in inJson:
        self.mPlaneNormal = inJson['plane normal']
      else:
        self.mPlaneNormal = [0.0, 1.0, 0.0]

    if gParaViewCatalystVersionFlag < 50502:
      if 'side to keep' in inJson:
        whichSide = inJson['side to keep']
        if whichSide == 'positive':
          self.mInsideOut = 0
        elif whichSide == 'negative':
          self.mInsideOut = 1
        else:
          if PhactoriDbg():
            myDebugPrint3('ParseParametersFromJson warning: side to keep \
              should be positive or negative, but is not, using positive\n')
      else:
        self.mInsideOut = 0
    else:
      if 'side to keep' in inJson:
        whichSide = inJson['side to keep']
        if whichSide == 'positive':
          self.mInvert = 0
        elif whichSide == 'negative':
          self.mInvert = 1
        else:
          if PhactoriDbg():
            myDebugPrint3('ParseParametersFromJson warning: side to keep \
              should be positive or negative, but is not, using positive\n')
      else:
        self.mInvert = 0

    if 'cut type' in inJson:
      cutType = inJson['cut type']
      self.mCrinkleSetting = 0
      if cutType == 'smooth':
        self.mCrinkleSetting = 0
      elif cutType == 'crinkle':
        self.mCrinkleSetting = 1
      else:
        if PhactoriDbg():
          myDebugPrint3('ParseParametersFromJson warning: cut type \
            should be smooth or crinkle, but is not, using smooth\n')
        self.mCrinkleSetting = 0
    else:
      self.mCrinkleSetting = 0  #debatable default

  def CreateParaViewFilter(self, inInputFilter):
    errStr = """error!  this is a base class which should have
                CreateParaViewFilter method overridden and not
                used"""
    if PhactoriDbg():
      myDebugPrint3(errStr)
    raise Exception(errStr)

  def CalculateUpdatedOriginAndNormal(
          self, inIncomingPvFilter, outCalculatedOrigin, outCalculatedNormal):
    """find a point on the plane and a normal to the plane, allowing for
       need to update location of points due to the point or points on the
       plane being relative points or data points"""

    if PhactoriDbg():
      myDebugPrint3(
        "PhactoriPlaneOpBase::CalculateUpdatedOriginAndNormal entered\n")

    viewBoundsArray = [None]
    if self.mUseThreePointsOnPlane:
      #defining plane by three points on plane
      originToUse = \
          self.mPtOnPlaneInfoA.GetCurrentGeometricPointWithDisplacement(
          inIncomingPvFilter, viewBoundsArray, True)
      outCalculatedOrigin[0] = originToUse[0]
      outCalculatedOrigin[1] = originToUse[1]
      outCalculatedOrigin[2] = originToUse[2]
      if PhactoriDbg():
        myDebugPrint3('  calculating normal from 3 points\n')
      tempPtB = self.mPtOnPlaneInfoB.GetCurrentGeometricPointWithDisplacement(
          inIncomingPvFilter, viewBoundsArray, False)
      tempPtC = self.mPtOnPlaneInfoC.GetCurrentGeometricPointWithDisplacement(
          inIncomingPvFilter, viewBoundsArray, False)
      vecAB = vecFromAToB(outCalculatedOrigin, tempPtB)
      vecAC = vecFromAToB(outCalculatedOrigin, tempPtC)
      vecCrossProduct2(outCalculatedNormal, vecAB, vecAC)
      #if we have 2 collinear points, our cross product will go to zero, and
      #we have to pick an arbitrary normal; we should be more selective
      #about this test
      if vecMagnitudeSquared(outCalculatedNormal,) < 0.000000000000000000001:
        if PhactoriDbg():
          myDebugPrint3('  two or more points are colinear, using 0,1,0\n')
        outCalculatedNormal[0] = 0.0
        outCalculatedNormal[1] = 1.0
        outCalculatedNormal[2] = 0.0
      else:
        vecNormalize2(outCalculatedNormal, outCalculatedNormal)
    else:
      originToUse = self.mPtOnPlaneInfo. \
        GetCurrentGeometricPointWithDisplacement(
        inIncomingPvFilter, viewBoundsArray, True)
      outCalculatedOrigin[0] = originToUse[0]
      outCalculatedOrigin[1] = originToUse[1]
      outCalculatedOrigin[2] = originToUse[2]
      outCalculatedNormal[0] = self.mPlaneNormal[0]
      outCalculatedNormal[1] = self.mPlaneNormal[1]
      outCalculatedNormal[2] = self.mPlaneNormal[2]

    if PhactoriDbg():
      myDebugPrint3(
        "PhactoriPlaneOpBase::CalculateUpdatedOriginAndNormal returning\n")


  def MayChangeWithData(self):
    if self.mUseThreePointsOnPlane:
      return (self.mPtOnPlaneInfoA.mMayChangeWithData or \
              self.mPtOnPlaneInfoB.mMayChangeWithData or \
              self.mPtOnPlaneInfoC.mMayChangeWithData)
    else:
      return self.mPtOnPlaneInfo.mMayChangeWithData


class PhactoriSliceOperation(PhactoriPlaneOpBase):
  """slice operation, settings and handling of slice filter"""

  def CreateParaViewFilter(self, inInputFilter):

    #don't need our own init code at this point, but this is how it would be
    #added
    #def __init__(self):
    #    MySuperClass.__init__(self)

    """create the slice plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriSliceOperation.CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()
    newParaViewFilter = Slice(Input = inInputFilter, SliceType = "Plane")
    newParaViewFilter.Triangulatetheslice = 0

    self.UpdateSlice(inInputFilter, newParaViewFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriSliceOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

  def DoUpdateDueToChangeInData(self, inIncomingPvFilter,
      outOutgoingPvFilter):
    """the PhactoriSliceOperation may need to update if the point on
       the slice plane was tied to a node, element, or variable min/max
       location"""
    if PhactoriDbg():
      myDebugPrint3("PhactoriSliceOperation::"
          "DoUpdateDueToChangeInData override executing\n")

    if self.MayChangeWithData() == False:
      if PhactoriDbg():
        myDebugPrint3("PhactoriSlicePlaneOperation::"
            "DoUpdateDueToChangeInData returning (absolute point or points)\n")
      return

    self.UpdateSlice(inIncomingPvFilter, outOutgoingPvFilter)

    if PhactoriDbg():
      myDebugPrint3("PhactoriSlicePlaneOperation::"
          "DoUpdateDueToChangeInData override returning\n")

  def UpdateSlice(self, inIncomingPvFilter, ioOutgoingPvFilter):
    """using the current info on the slice, get all the paraview stuff
       set up correctly"""

    if PhactoriDbg():
      myDebugPrint3("PhactoriSlicePlaneOperation::UpdateSlice entered\n")

    originToUse = [0,0,0]
    normalToUse = [0,1,0]
    self.CalculateUpdatedOriginAndNormal(
            inIncomingPvFilter, originToUse, normalToUse)

    if PhactoriDbg():
      myDebugPrint3('  updateslice using normal: ' + \
              str(normalToUse) + '\n')
    ioOutgoingPvFilter.SliceType.Normal = normalToUse

    if PhactoriDbg():
      myDebugPrint3('  updateslice using origin: ' + str(originToUse) + '\n')
    ioOutgoingPvFilter.SliceType.Origin = originToUse

    #these aren't changing yet
    ioOutgoingPvFilter.Crinkleslice = self.mCrinkleSetting

    if PhactoriDbg():
      myDebugPrint3("PhactoriSlicePlaneOperation::UpdateSlice returning\n")



class PhactoriClipPlaneOperation(PhactoriPlaneOpBase):
  """clip plane operation, settings and handling of slice filter"""

  def CreateParaViewFilter(self, inInputFilter):

    #don't need our own init code at this point, but this is how it would be
    #added
    #def __init__(self):
    #    MySuperClass.__init__(self)

    """create the clip plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriClipPlaneOperation.CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()
    newParaViewFilter = Clip(Input = inInputFilter, ClipType = "Plane")
    self.UpdateClip(inInputFilter, newParaViewFilter)

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriClipPlaneOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

  def DoUpdateDueToChangeInData(self, inIncomingPvFilter,
      outOutgoingPvFilter):
    """the PhactoriClipPlaneOperation may need to update if the point on
       the clip plane was tied to a node, element, or variable min/max
       location"""
    if PhactoriDbg():
      myDebugPrint3("PhactoriClipPlaneOperation::"
          "DoUpdateDueToChangeInData override executing\n")

    if self.MayChangeWithData() == False:
      if PhactoriDbg():
        myDebugPrint3("PhactoriClipPlaneOperation::"
            "DoUpdateDueToChangeInData returning (absolute point or points)\n")
      return

    self.UpdateClip(inIncomingPvFilter, outOutgoingPvFilter)

    if PhactoriDbg():
      myDebugPrint3("PhactoriClipPlaneOperation::"
          "DoUpdateDueToChangeInData override returning\n")

  def UpdateClip(self, inIncomingPvFilter, ioOutgoingPvFilter):
    """using the current info on the clip, get all the paraview stuff
       set up correctly"""

    if PhactoriDbg():
      myDebugPrint3("PhactoriClipPlaneOperation::UpdateClip entered\n")

    originToUse = [0,0,0]
    normalToUse = [0,1,0]
    self.CalculateUpdatedOriginAndNormal(
            inIncomingPvFilter, originToUse, normalToUse)

    if PhactoriDbg():
      myDebugPrint3('  updateclip using normal: ' + \
              str(normalToUse) + '\n')
    ioOutgoingPvFilter.ClipType.Normal = normalToUse

    if PhactoriDbg():
      myDebugPrint3('  updateclip using origin: ' + str(originToUse) + '\n')
    ioOutgoingPvFilter.ClipType.Origin = originToUse

    #these aren't changing yet
    ioOutgoingPvFilter.Crinkleclip = self.mCrinkleSetting
    if gParaViewCatalystVersionFlag < 50502:
      ioOutgoingPvFilter.InsideOut = self.mInsideOut
    else:
      ioOutgoingPvFilter.Invert = self.mInvert

    if PhactoriDbg():
      myDebugPrint3("PhactoriClipPlaneOperation::UpdateClip returning\n")


class PhactoriCylinderClipOperation(PhactoriOperationSpecifics):
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mCenterPtInfo = PhactoriUserPointInfo()
    self.mAxis = [0.0, 1.0, 0.0]
    self.mRadius = 0.5
    self.mRadiusRelAbsFlag = True
    self.mKeepInside = True
    self.mCrinkleSetting = None

  def ParseParametersFromJson(self, inJson):
    if PhactoriDbg(100):
      myDebugPrint3(
        'PhactoriCylinderClipOperation.ParseParametersFromJson entered\n', 100)
    self.mCenterPtInfo.UserPointInfoParseParametersFromJson(
            inJson, "center at ", "")

    absRadiusKey = 'absolute radius'
    relRadiusKey = 'relative radius'
    if absRadiusKey in inJson:
      self.mRadius = inJson[absRadiusKey]
      self.mRadiusRelAbsFlag = False
    elif relRadiusKey in inJson:
      self.mRadius = inJson[relRadiusKey]
      self.mRadiusRelAbsFlag = True
    else:
      self.mRadius = 0.5
      self.mRadiusRelAbsFlag = True

    axisKey = 'axis'
    if axisKey in inJson:
      self.mAxis = inJson[axisKey]

    if 'keep inside cylinder' in inJson:
      self.mKeepInside = inJson['keep inside cylinder']
    else:
      self.mKeepInside = True

    if 'cut type' in inJson:
      cutType = inJson['cut type']
      self.mCrinkleSetting = 0
      if cutType == 'smooth':
        self.mCrinkleSetting = 0
      elif cutType == 'crinkle':
        self.mCrinkleSetting = 1
      else:
        if PhactoriDbg():
          myDebugPrint3('ParseParametersFromJson warning: cut type \
            should be smooth or crinkle, but is not, using smooth\n')
        self.mCrinkleSetting = 0
    else:
      self.mCrinkleSetting = 0  #debatable default

    #self.PrintSelf()
    if PhactoriDbg(100):
      myDebugPrint3(
        'PhactoriCylinderClipOperation.ParseParametersFromJson returning\n',
        100)

  def PrintSelf(self):
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriCylinderClipOperation.PrintSelf entered\n" + \
        "Center Point:', 100)
    self.mCenterPtInfo.PrintSelf()
    if PhactoriDbg():
      myDebugPrint3('mRadiusRelAbsFlag: ' + str(self.mRadiusRelAbsFlag) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mRadius: ' + str(self.mRadius) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mAxis: ' + str(self.mAxis) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mKeepInside: ' + str(self.mKeepInside) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mCrinkleSetting: ' + str(self.mCrinkleSetting) + '\n')

  def CreateParaViewFilter(self, inInputFilter):
    """create the cylinder clip filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriCylinderClipOperation.CreateParaViewFilter ' +\
        'entered\n', 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    newParaViewFilter = Clip(Input = inInputFilter, ClipType = "Cylinder")
    ## Properties modified on clip4.ClipType
    #clip4.ClipType.Center = [10.93421983718872, -0.8890213072299957, 0.0]
    #clip4.ClipType.Axis = [1.0, 0.0, 0.0]
    #clip4.ClipType.Radius = 0.75

    newParaViewFilter.ClipType.Axis = self.mAxis

    viewBoundsArray = [None]
    centerToUse = self.mCenterPtInfo.GetCurrentGeometricPointWithDisplacement(
        inInputFilter, viewBoundsArray, True)

    newParaViewFilter.ClipType.Center = [centerToUse[0],
        centerToUse[1],
        centerToUse[2]]

    if PhactoriDbg():
      myDebugPrint3("  about to find radius--absolute or calced from relative\n")
    myViewBounds = viewBoundsArray[0]
    if self.mRadiusRelAbsFlag:
      if myViewBounds == None:
        myViewBounds = GetGlobalDataBoundsParallel(inInputFilter)
      dxx = myViewBounds[1] - myViewBounds[0]
      dyy = myViewBounds[3] - myViewBounds[2]
      dzz = myViewBounds[5] - myViewBounds[4]
      bounds_diagonal = math.sqrt(dxx*dxx + dyy*dyy + dzz*dzz)
      newParaViewFilter.ClipType.Radius = self.mRadius * bounds_diagonal
    else:
      newParaViewFilter.ClipType.Radius = self.mRadius

    newParaViewFilter.Crinkleclip = self.mCrinkleSetting
    if self.mKeepInside:
      newParaViewFilter.InsideOut = 1
    else:
      newParaViewFilter.InsideOut = 0

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3(
        'PhactoriCylinderClipOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

  def MayChangeWithData(self):
    if self.mCenterPtInfo.mMayChangeWithData == True:
      return True
    if self.mRadiusRelAbsFlag == True:
      return True
    return False


class PhactoriBoxClipOperation(PhactoriOperationSpecifics):
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.mCenterPtInfo = PhactoriUserPointInfo()
    self.mExtentsRelAbsFlag = [True, True, True]
    self.mExtents = [0.5, 0.5, 0.5]
    self.mRotations = [0.0, 0.0, 0.0]
    self.mKeepInside = True
    self.mCrinkleSetting = None

  def ParseParametersFromJson(self, inJson):
    if PhactoriDbg(100):
      myDebugPrint3(
        'PhactoriBoxClipOperation.ParseParametersFromJson entered\n', 100)
    self.mCenterPtInfo.UserPointInfoParseParametersFromJson(
            inJson, "center at ", "")

    absExtentsKey = 'absolute extents'
    relExtentsKey = 'relative extents'
    if absExtentsKey in inJson:
      self.mExtents = inJson[absExtentsKey]
      self.mExtentsRelAbsFlag = [False, False, False]
    elif relExtentsKey in inJson:
      self.mExtents = inJson[relExtentsKey]
      self.mExtentsRelAbsFlag = [True, True, True]
    else:
      self.mExtents = [0.5, 0.5, 0.5]
      self.mExtentsRelAbsFlag = [True, True, True]

    rotationsKey = 'rotations'
    if rotationsKey in inJson:
      self.mRotations = inJson[rotationsKey]

    if 'keep inside box' in inJson:
      self.mKeepInside = inJson['keep inside box']
    else:
      self.mKeepInside = True

    if 'cut type' in inJson:
      cutType = inJson['cut type']
      self.mCrinkleSetting = 0
      if cutType == 'smooth':
        self.mCrinkleSetting = 0
      elif cutType == 'crinkle':
        self.mCrinkleSetting = 1
      else:
        if PhactoriDbg():
          myDebugPrint3('ParseParametersFromJson warning: cut type \
            should be smooth or crinkle, but is not, using smooth\n')
        self.mCrinkleSetting = 0
    else:
      self.mCrinkleSetting = 0  #debatable default

    #self.PrintSelf()
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriBoxClipOperation.ParseParametersFromJson returning\n', 100)

  def PrintSelf(self):
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriBoxClipOperation.PrintSelf entered\nCenter Point:', 100)
    self.mCenterPtInfo.PrintSelf()
    if PhactoriDbg():
      myDebugPrint3('mExtentsRelAbsFlag: ' + str(self.mExtentsRelAbsFlag) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mExtents: ' + str(self.mExtents) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mRotations: ' + str(self.mRotations) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mKeepInside: ' + str(self.mKeepInside) + '\n')
    if PhactoriDbg():
      myDebugPrint3('mCrinkleSetting: ' + str(self.mCrinkleSetting) + '\n')

  def CreateParaViewFilter(self, inInputFilter):
    """create the box clip box for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriBoxClipOperation.CreateParaViewFilter entered\n', 100)
    #info in block class should already be parsed and checked

    savedActiveSource = GetActiveSource()

    newParaViewFilter = Clip(Input = inInputFilter, ClipType = "Box")

    viewBoundsArray = [None]
    centerToUse = self.mCenterPtInfo.GetCurrentGeometricPointWithDisplacement(
        inInputFilter, viewBoundsArray, True)

    newParaViewFilter.ClipType.Position = [centerToUse[0],
        centerToUse[1],
        centerToUse[2]]

    if PhactoriDbg():
      myDebugPrint3("  about to find extents--absolute or calced from relative\n")
    myViewBounds = viewBoundsArray[0]
    if self.mExtentsRelAbsFlag[0]:
      if myViewBounds == None:
        myViewBounds = GetGlobalDataBoundsParallel(inInputFilter)
      halfXExt = self.mExtents[0] * (myViewBounds[1] - myViewBounds[0]) * 0.5
    else:
      halfXExt = self.mExtents[0] * 0.5

    if self.mExtentsRelAbsFlag[1]:
      if myViewBounds == None:
        myViewBounds = GetGlobalDataBoundsParallel(inInputFilter)
      halfYExt = self.mExtents[1] * (myViewBounds[3] - myViewBounds[2]) * 0.5
    else:
      halfYExt = self.mExtents[1] * 0.5

    if self.mExtentsRelAbsFlag[2]:
      if myViewBounds == None:
        myViewBounds = GetGlobalDataBoundsParallel(inInputFilter)
      halfZExt = self.mExtents[2] * (myViewBounds[5] - myViewBounds[4]) * 0.5
    else:
      halfZExt = self.mExtents[2] * 0.5

    if PhactoriDbg():
      myDebugPrint3("  halfExtents: " + str(halfXExt) + ", " + str(halfYExt) + ", " + str(halfZExt) + "\n")

    newParaViewFilter.ClipType.Bounds = [-halfXExt, halfXExt,
        -halfYExt, halfYExt, -halfZExt, halfZExt]
    newParaViewFilter.ClipType.Rotation = [self.mRotations[0],
        self.mRotations[1],
        self.mRotations[2]]

    newParaViewFilter.Crinkleclip = self.mCrinkleSetting
    if gParaViewCatalystVersionFlag < 50502:
      if self.mKeepInside:
        newParaViewFilter.InsideOut = 1
      else:
        newParaViewFilter.InsideOut = 0
    else:
      if self.mKeepInside:
        newParaViewFilter.Invert = 1
      else:
        newParaViewFilter.Invert = 0

    SetActiveSource(newParaViewFilter)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3('PhactoriBoxClipOperation.CreateParaViewFilter returning\n', 100)

    return newParaViewFilter

  def MayChangeWithData(self):
    if self.mCenterPtInfo.mMayChangeWithData == True:
      return True
    if self.mExtentsRelAbsFlag[0] == True:
      return True
    if self.mExtentsRelAbsFlag[1] == True:
      return True
    if self.mExtentsRelAbsFlag[2] == True:
      return True
    return False

global FlagToTestOutgoingPvGeometryFilter
FlagToTestOutgoingPvGeometryFilter = False

class ClosestNPointsFromEachProcess:
  """class to manage parallel operations where we obtain a set of N points on
     each process and then share those points from each process to all
     process.  used for doing things like finding closest N points from one
     operation to another on each process and then sharing and sorting and
     finding closest overall"""

  def __init__(self):
    self.mPointsPerProcess = None
    self.mMyPid = -1
    self.mNumProcesses = -1
    self.mClosestNRefs = None

    #closest N local points (local id and id from list and geometry xyzs
    #from local and from list and distance squared
    self.mThisProcIds = vtk.vtkIntArray()
    self.mTargetMatchIds = vtk.vtkIntArray()
    self.mThisProcXyzs = vtk.vtkDoubleArray()
    self.mTargetMatchXyzs = vtk.vtkDoubleArray()
    self.mDistSqrds = vtk.vtkDoubleArray()

  def ParallelSetupFromFindClosestNPointsToList(self, inMyProcData):
    self.mPointsPerProcess = inMyProcData.mNumToFind
    mypid, numproc = UseMPIToFillInSharedListNPerProcess(
      inMyProcData.mThisProcIds, self.mPointsPerProcess, self.mThisProcIds, 0)
    UseMPIToFillInSharedListNPerProcess(
      inMyProcData.mTargetMatchIds, self.mPointsPerProcess, self.mTargetMatchIds, 0)
    UseMPIToFillInSharedListNPerProcess(
      inMyProcData.mDistSqrds, self.mPointsPerProcess, self.mDistSqrds, 1)
    UseMPIToFillInSharedListNPerProcess(
      inMyProcData.mThisProcXyzs, self.mPointsPerProcess * 3, self.mThisProcXyzs, 1)
    UseMPIToFillInSharedListNPerProcess(
      inMyProcData.mTargetMatchXyzs, self.mPointsPerProcess * 3, self.mTargetMatchXyzs, 1)

    self.mMyPid = mypid
    self.mNumProcesses = numproc

  def FindNClosestPointsNdxInList(self, inHowMany):
    import heapq

    numPoints = self.mPointsPerProcess * self.mNumProcesses
    if PhactoriDbg():
      myDebugPrint3(str(inHowMany) + " " + str(numPoints) + "\n")
    if numPoints < inHowMany:
      exit(-1)
    topNHeap = []
    for ii in range(0, inHowMany):
      heapq.heappush(topNHeap, (-self.mDistSqrds.GetValue(ii), ii))
    for ii in range(inHowMany, numPoints):
      testval = -self.mDistSqrds.GetValue(ii)
      if testval > topNHeap[0][0]:
         heapq.heappushpop(topNHeap, (testval, ii))

    self.mClosestNRefs = []
    for ii in range(0, inHowMany):
      item = heapq.heappop(topNHeap)
      self.mClosestNRefs.append((-item[0], item[1]))

    self.mClosestNRefs.reverse()

    if PhactoriDbg():
      myDebugPrint3("FindNClosestPointsNdxInList:\n" + str(self.mClosestNRefs) + "\n")

  def MakeNClosestPointsNdxInListTableString(self, inMakeHeaderLine = True):
    if inMakeHeaderLine:
      retStr = "index, distance, distance squared, " \
               "source node id, source x, source y, source z, " \
               "target node id, target x, target y, target z\n"
    else:
      retStr = ""
    for ii in range(0, len(self.mClosestNRefs)):
      ndx = self.mClosestNRefs[ii][1]
      srcndid = self.mThisProcIds.GetValue(ndx)
      tgtndid = self.mTargetMatchIds.GetValue(ndx)
      distsqrd = self.mDistSqrds.GetValue(ndx)
      dist = math.sqrt(distsqrd)
      pndx = ndx*3
      srcx = self.mThisProcXyzs.GetValue(pndx)
      srcy = self.mThisProcXyzs.GetValue(pndx+1)
      srcz = self.mThisProcXyzs.GetValue(pndx+2)
      tgtx = self.mTargetMatchXyzs.GetValue(pndx)
      tgty = self.mTargetMatchXyzs.GetValue(pndx+1)
      tgtz = self.mTargetMatchXyzs.GetValue(pndx+2)
      retStr += str(ii) + ", " + \
          str(dist) + ", " + \
          str(distsqrd) + ", " + \
          str(srcndid) + ", " + \
          str(srcx) + ", " + \
          str(srcy) + ", " + \
          str(srcz) + ", " + \
          str(tgtndid) + ", " + \
          str(tgtx) + ", " + \
          str(tgty) + ", " + \
          str(tgtz) + "\n"
    return retStr

  def ToStr(self):
    retStr = "parallel shared closest " + str(self.mPointsPerProcess) + " points:\n" + \
    "index:process:pndx source id, target id, dist sqrd: source xyz: target xyz\n"
    for ii in range(0, self.mPointsPerProcess * self.mNumProcesses):
      pp = ii*3
      adst = str(ii) + ":" + \
             str(ii/self.mPointsPerProcess) + ":" + \
             str(ii%self.mPointsPerProcess) + ": " + \
             str(self.mThisProcIds.GetValue(ii)) + ", " + \
             str(self.mTargetMatchIds.GetValue(ii)) + ", " + \
             str(self.mDistSqrds.GetValue(ii)) + ": " + \
             str(self.mThisProcXyzs.GetValue(pp)) + ", " + \
             str(self.mThisProcXyzs.GetValue(pp+1)) + ", " + \
             str(self.mThisProcXyzs.GetValue(pp+2)) + ": " + \
             str(self.mTargetMatchXyzs.GetValue(pp)) + ", " + \
             str(self.mTargetMatchXyzs.GetValue(pp+1)) + ", " + \
             str(self.mTargetMatchXyzs.GetValue(pp+2)) + "\n"
      retStr += adst
    return retStr


class BlockRecursionControlItem:
  """see DoMethodPerBlock(); mOperationPerBlock should be set to a method
     which takes 1 parameter, mParameters should be set to the parameter
     instance which will be passes to the mOperationPerBlock call"""
  def __init__(self):
    self.mOperationPerBlock = None
    self.mParameters = None

class PhactoriOperationBlock:
  """manages one stage of the data pipeline, analogous to ParaView Filter

  An instance of this class represents and manages one stage of the data
  pipeline which has been set up for management by phatori.  It creates and
  holds a reference to a ParaView/Catalyst/vtk Filter which is doing the
  real work--this class might be thought of as an adapter or interface from
  phactori to ParaView.  As currently implemented, an operation is assumed to
  have only one input and one output.  Multiple operations can have the same
  input, so a tree structure is allowed rather than just a linear pipe.
  Operations with multiple inputs and outputs are conceiveable, and may be
  added pending user requirements.
  The instance is presumed to contain a name unique amount the operation
  blocks and keeps a reference to the input operation (by name), the
  ParaView/Catalyst filter which is built, and some flags determining where
  we are in the construction process.
  """
  def __init__(self):
    self.mName = ""
    self.mInputOperationName = None
    self.mType = ""
    self.mHasBeenConstructed = False
    self.mParaViewFilter = None
    #for keeping track of data bounds for this operation and only updating
    #it when necessary
    self.mDataBoundsIsCurrent = False
    self.mDataBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    self.mOperationSpecifics = None

    #If this is a purely internal pipeline filter, this will remain None.
    #If somebody uses the filter for rendering, it will be created and
    #handed off for rendering--only needs creation once.
    self.mOutgoingGeometryFilter = None

    #stores annotation time source, if one gets created associated with this
    #operation
    self.mParaViewTimeAnnotationSource = None

    self.mRecursionParameterStore = {}


  def GetPvFilter(self):
    return self.mParaViewFilter

  def CreateOutgoingPvGeometryFilter(self):
    savedActiveSource = GetActiveSource()
    self.mOutgoingGeometryFilter = ExtractSurface(self.GetPvFilter())
    if PhactoriDbg():
      myDebugPrint3("PhactoriOperationBlock::CreateOutgoingPvGeometryFilter\n"
          "created self.mOutgoingGeometryFilter: " + \
      str(self.mOutgoingGeometryFilter) + "\n")

    SetActiveSource(savedActiveSource)

  def GetOutgoingPvGeometryFilter(self):
    global FlagToTestOutgoingPvGeometryFilter
    if FlagToTestOutgoingPvGeometryFilter == False:
      if PhactoriDbg():
        myDebugPrint3("PhactoriOperationBlock::GetOutgoingPvGeometryFilter\n"
            "not creating or using external geometry filter: ")
      return self.GetPvFilter()
    if self.mOutgoingGeometryFilter == None:
      self.CreateOutgoingPvGeometryFilter()
    return self.mOutgoingGeometryFilter

  def DoUpdateDueToChangeInData(self, ioPipeAndViewsState):
    outputPvFilter = self.mParaViewFilter
    if self.mInputOperationName == None:
      inputOperation = ioPipeAndViewsState.mIncomingDefaultOperation
    else:
      inputOperation = ioPipeAndViewsState.mOperationBlocks[
          self.mInputOperationName]
    if inputOperation == None:
      myDebugPrint3AndException(
        "PhactoriOperationBlock::DoUpdateDueToChangeInData:\n"
        "couldn't find input operation with name: " + \
        str(self.mInputOperationName) + "\n")

    self.mOperationSpecifics.DoUpdateDueToChangeInData(
        inputOperation.mParaViewFilter, outputPvFilter)

  def DoMethodPerBlock(self, inRecursionControlItem):
    """DoMethodPerBlock is a generic method for doing recursion through the
       multiblock dataset and doing something (a callback) on a per leaf block
       basis.  Gets the clientside data and calls DoMethodPerBlockRecurse1
       which calls itself on internal nodes and calls
       inRecursionControlItem.mOperationToDoPerBlock on leaf block nodes"""
    pvClientSideData = self.GetPvFilter().GetClientSideObject().\
            GetOutputDataObject(0)
    if pvClientSideData == None:
      if PhactoriDbg(100):
        myDebugPrint3(
          'DoMethodPerBlock: pvClientSideData is None, returning',100)

    self.DoMethodPerBlockRecurse1(pvClientSideData, inRecursionControlItem)

  def DoMethodPerBlockRecurse1(self, inInputCsData, inRecursionControlItem):
    """DoMethodPerBlockRecurse1 is a generic method for doing recursion through
       the multiblock dataset and doing something (a callback) on a per leaf block
       basis.  Called by DoMethodPerBlock which got clientside data and calls
       itself on internal nodes and calls
       inRecursionControlItem.mOperationToDoPerBlock on leaf block nodes"""
    if PhactoriDbg(100):
      myDebugPrint3('DoMethodPerBlockRecurse1 entered\n', 100)

    icsdClassname = inInputCsData.GetClassName()
    if icsdClassname == "vtkMultiBlockDataSet" or \
       icsdClassname == "vtkExodusIIMultiBlockDataSet":
      #myDebugPrint3('recursing: ' + icsdClassname + '\n')
      numBlocks = inInputCsData.GetNumberOfBlocks()
      for ii in range(0, numBlocks):
        oneBlock = inInputCsData.GetBlock(ii)
        if(oneBlock != None):
          self.DoMethodPerBlockRecurse1(oneBlock, inRecursionControlItem)
    else:
      inRecursionControlItem.mOperationToDoPerBlock(inInputCsData,
              inRecursionControlItem.mParameters)

    if PhactoriDbg(100):
      myDebugPrint3('DoMethodPerBlockRecurse1 returning\n', 100)

  def OutputElementListFromOneBlockToFile(self, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("OutputElementListFromOneBlockToFile entered\n")
    cellData = inInputCsData.GetCellData()
    numTuplesX = cellData.GetNumberOfTuples()
    numArrays = cellData.GetNumberOfArrays()

    global WriteEachDeadCellElementToFiles
    global WriteDeadCellSummaryFile

    #if inParameters.mBlockCount == 0:
    if (inParameters.mFlag1 == 0) and (numTuplesX > 0):
      inParameters.mFlag1 = 1
      if WriteEachDeadCellElementToFiles:
        inParameters.mOutFileff.write("element index")
        for arrayIndex in range(0, numArrays):
            inParameters.mOutFileff.write(
              ", " + str(cellData.GetArray(arrayIndex).GetName()))
        inParameters.mOutFileff.write("\n")

    if WriteEachDeadCellElementToFiles:
      for cellIndex in range(0, numTuplesX):
          inParameters.mOutFileff.write(
                  str(cellIndex + inParameters.mElementCount))
          for arrayIndex in range(0, numArrays):
              inParameters.mOutFileff.write(
                ", " + str(cellData.GetArray(arrayIndex).GetTuple1(cellIndex)))
          inParameters.mOutFileff.write("\n")

    killed_criteria_array = cellData.GetArray("KILLED")
    for cellIndex in range(0, numTuplesX):
      kval = killed_criteria_array.GetTuple1(cellIndex)
      kvalindex = int(kval)
      if (kvalindex >= 0) and (kvalindex < 10):
        inParameters.mKilledByCriteriaCount[kvalindex] += 1


    inParameters.mElementCount += inInputCsData.GetNumberOfCells()
    inParameters.mNodeCount += inInputCsData.GetNumberOfPoints()
    inParameters.mBlockCount += 1

    if PhactoriDbg(100):
      myDebugPrint3("OutputElementListFromOneBlockToFile returning\n")

  class FindClosestNPointsToListParams:
    """recursion structure for FindClosestNPointsToList().  Also servers to
       store/track data for passing back answer"""
    def __init__(self, inParentOperation, inNumToFind,
                 inTargetGlobalNodeIdList, inTargetPointXyzList):
      self.mParentOperation = inParentOperation
      self.mNumToFind = inNumToFind

      #list of points we need to find closest local points to
      #(global node ids and geometry xyzs)
      self.mTargetIds = inTargetGlobalNodeIdList
      self.mTargetXyzs = inTargetPointXyzList

      self.mKdtree = None

      #closest N local points (local id and id from list and geometry xyzs
      #from local and from list and distance squared
      self.mThisProcIds = vtk.vtkIntArray()
      self.mTargetMatchIds = vtk.vtkIntArray()
      self.mThisProcXyzs = vtk.vtkDoubleArray()
      self.mTargetMatchXyzs = vtk.vtkDoubleArray()
      self.mDistSqrds = vtk.vtkDoubleArray()

      self.mThisProcIds.SetNumberOfValues(inNumToFind)
      self.mTargetMatchIds.SetNumberOfValues(inNumToFind)
      self.mThisProcXyzs.SetNumberOfValues(inNumToFind*3)
      self.mTargetMatchXyzs.SetNumberOfValues(inNumToFind*3)
      self.mDistSqrds.SetNumberOfValues(inNumToFind)

      #set default values indicating nothing found in those entries
      for ii in range(0,inNumToFind):
        self.mThisProcIds.SetValue(ii, -1)
        self.mTargetMatchIds.SetValue(ii, -1)
        self.mDistSqrds.SetValue(ii, sys.float_info.max)
      self.mMinDistSqrd = sys.float_info.max
      #index of the item that currently has the biggest distance
      self.mcfndx = 0

    def SetUpKdtree(self):
      """take the points in this instance (presumably the target object) and
         put them in a kdtree (using scipy) so we can find the closest point
         in log(n) time"""
      from scipy import spatial
      kdtreepts = []
      numpts = self.mTargetIds.GetNumberOfValues()
      for ii in range(0,numpts):
        pndx = ii*3
        kdtreepts.append([
          self.mTargetXyzs.GetValue(pndx),
          self.mTargetXyzs.GetValue(pndx+1),
          self.mTargetXyzs.GetValue(pndx+2)])
      self.mKdtree = None
      self.mKdtree = spatial.KDTree(kdtreepts)

    def TestPointWithKdtree(self, inSrcId, inSrcXyz):
      nrstdist, nrstndx = self.mKdtree.query(inSrcXyz)
      #print "nrstdist: ", str(nrstdist)
      #print "nrstndx: ", str(nrstndx)
      tgtId = self.mTargetIds.GetValue(nrstndx)
      tgtXyz = self.mKdtree.data[nrstndx]
      #self.TestPointSub1(inSrcId, inSrcXyz[0], inSrcXyz[1], inSrcXyz[2],
      #    tgtId, tgtXyz[0], tgtXyz[1], tgtXyz[2])
      self.TestPointSub2(inSrcId, inSrcXyz, tgtId, tgtXyz, nrstdist)

    def TestPoint(self, inId, inXyz):
      """given an xyz point in the local processor (and its global node id)
         see if it is closer to any of the target points than the current
         set of nearest points and, if so, put it in the set, dropping
         others out if necessary"""
      if PhactoriDbg():
        myDebugPrint3("TestPoint id " + str(inId))
      tgtxyzs = self.mTargetXyzs
      numtgtpts = self.mTargetIds.GetNumberOfValues()
      for pp in range(0, numtgtpts):
        tgtid = self.mTargetIds.GetValue(pp)
        ndx = pp*3
        self.TestPointSub1(inId, inXyz[0], inXyz[1], inXyz[2],
            tgtid, tgtxyzs.GetValue(ndx), tgtxyzs.GetValue(ndx+1),
            tgtxyzs.GetValue(ndx+2))

    def TestPointSub2(self, inSrcId, inSrcXyz, inTgtId, inTgtXyz, inNrstdist):
      dstsqd = inNrstdist*inNrstdist
      if dstsqd >= self.mMinDistSqrd:
        #we already have mNumToFind closer
        #if PhactoriDbg():
        #  myDebugPrint3("inSrcId inTgtId " + str(inSrcId) + " " + \
        #      str(inTgtId) + \
        #      " too far: " + str(dstsqd) + " >= " + \
        #      str(self.mMinDistSqrd) + "\n")
        return

      #replace the previous point that was furthest
      self.mDistSqrds.SetValue(self.mcfndx, dstsqd)
      self.mThisProcIds.SetValue(self.mcfndx, inSrcId)
      self.mTargetMatchIds.SetValue(self.mcfndx, inTgtId)
      gndx = self.mcfndx * 3
      self.mThisProcXyzs.SetValue(gndx, inSrcXyz[0])
      self.mThisProcXyzs.SetValue(gndx+1, inSrcXyz[1])
      self.mThisProcXyzs.SetValue(gndx+2, inSrcXyz[2])
      self.mTargetMatchXyzs.SetValue(gndx, inTgtXyz[0])
      self.mTargetMatchXyzs.SetValue(gndx+1, inTgtXyz[1])
      self.mTargetMatchXyzs.SetValue(gndx+2, inTgtXyz[2])
      #if PhactoriDbg():
      #  myDebugPrint3("closer point found put in index " + \
      #      str(self.mcfndx) + \
      #      "  dstsqrd: " + str(self.mDistSqrds.GetValue(self.mcfndx)) + "\n" + \
      #      "\nsource id xyz: " + \
      #      str(self.mThisProcIds.GetValue(self.mcfndx)) + "   " + \
      #      str(self.mThisProcXyzs.GetValue(gndx)) + ", " + \
      #      str(self.mThisProcXyzs.GetValue(gndx+1)) + ", " + \
      #      str(self.mThisProcXyzs.GetValue(gndx+2)) + ", " + \
      #      "\ntarget id xyz: " + \
      #      str(self.mTargetMatchIds.GetValue(self.mcfndx)) + "   " + \
      #      str(self.mTargetMatchXyzs.GetValue(gndx)) + ", " + \
      #      str(self.mTargetMatchXyzs.GetValue(gndx+1)) + ", " + \
      #      str(self.mTargetMatchXyzs.GetValue(gndx+2)) + "\n")

      #now find which in the current list has the biggest distance, as it is
      #next in line for replacement (we do this to avoid having to shift
      #elements every time
      self.mcfndx = 0
      self.mMinDistSqrd = self.mDistSqrds.GetValue(0)
      #if PhactoriDbg():
      #  myDebugPrint3("find next index to be replaced try 0 \n" + \
      #      str(self.mMinDistSqrd) + "\n")
      for kk in range(1, self.mNumToFind):
        #if PhactoriDbg():
        #  myDebugPrint3("try " + str(kk) + " " + \
        #      str(self.mDistSqrds.GetValue(kk)) + " >? " + \
        #      str(self.mMinDistSqrd)+ "\n")
        if self.mDistSqrds.GetValue(kk) > self.mMinDistSqrd:
          self.mcfndx = kk
          self.mMinDistSqrd = self.mDistSqrds.GetValue(kk)
          #if PhactoriDbg():
          #  myDebugPrint3("yes, now " + str(self.mcfndx) + " " + \
          #      str(self.mMinDistSqrd) + "\n")
      #if PhactoriDbg():
      #  myDebugPrint3("next to be replaced ndx: " + str(self.mcfndx) + \
      #      " sid: " + \
      #      str(self.mThisProcIds.GetValue(self.mcfndx)) + \
      #      " tid: " + \
      #      str(self.mTargetMatchIds.GetValue(self.mcfndx)) + \
      #      " dsqrd: " + \
      #      str(self.mDistSqrds.GetValue(self.mcfndx)) + "\n")

    def TestPointSub1(self, inId, inX, inY, inZ, tId, tX, tY, tZ):
      """given an xyz point in the local processor (and its global node id)
         see if it is closer to one of the target points than the current
         set of nearest points and, if so, put it in the set, dropping
         others out if necessary"""
      ddx = inX - tX
      ddy = inY - tY
      ddz = inZ - tZ
      dstsqd = ddx*ddx + ddy*ddy + ddz*ddz
      if dstsqd >= self.mMinDistSqrd:
        #we already have mNumToFind closer
        #if PhactoriDbg():
        #  myDebugPrint3("inId tId " + str(inId) + " " + str(tId) + \
        #      " too far: " + str(dstsqd) + " >= " + \
        #      str(self.mMinDistSqrd) + "\n")
        return

      #replace the previous point that was furthest
      self.mDistSqrds.SetValue(self.mcfndx, dstsqd)
      self.mThisProcIds.SetValue(self.mcfndx, inId)
      self.mTargetMatchIds.SetValue(self.mcfndx, tId)
      gndx = self.mcfndx * 3
      self.mThisProcXyzs.SetValue(gndx, inX)
      self.mThisProcXyzs.SetValue(gndx+1, inY)
      self.mThisProcXyzs.SetValue(gndx+2, inZ)
      self.mTargetMatchXyzs.SetValue(gndx, tX)
      self.mTargetMatchXyzs.SetValue(gndx+1, tY)
      self.mTargetMatchXyzs.SetValue(gndx+2, tZ)
      if PhactoriDbg():
        myDebugPrint3("closer point found put in index " + \
            str(self.mcfndx) + \
            "  dstsqrd: " + str(self.mDistSqrds.GetValue(self.mcfndx)) + "\n" + \
            "\nsource id xyz: " + \
            str(self.mThisProcIds.GetValue(self.mcfndx)) + "   " + \
            str(self.mThisProcXyzs.GetValue(gndx)) + ", " + \
            str(self.mThisProcXyzs.GetValue(gndx+1)) + ", " + \
            str(self.mThisProcXyzs.GetValue(gndx+2)) + ", " + \
            "\ntarget id xyz: " + \
            str(self.mTargetMatchIds.GetValue(self.mcfndx)) + "   " + \
            str(self.mTargetMatchXyzs.GetValue(gndx)) + ", " + \
            str(self.mTargetMatchXyzs.GetValue(gndx+1)) + ", " + \
            str(self.mTargetMatchXyzs.GetValue(gndx+2)) + "\n")

      #now find which in the current list has the biggest distance, as it is
      #next in line for replacement (we do this to avoid having to shift
      #elements every time
      self.mcfndx = 0
      self.mMinDistSqrd = self.mDistSqrds.GetValue(0)
      if PhactoriDbg():
        myDebugPrint3("find next index to be replaced try 0 \n" + \
            str(self.mMinDistSqrd) + "\n")
      for kk in range(1, self.mNumToFind):
        if PhactoriDbg():
          myDebugPrint3("try " + str(kk) + " " + \
              str(self.mDistSqrds.GetValue(kk)) + " >? " + \
              str(self.mMinDistSqrd)+ "\n")
        if self.mDistSqrds.GetValue(kk) > self.mMinDistSqrd:
          self.mcfndx = kk
          self.mMinDistSqrd = self.mDistSqrds.GetValue(kk)
          if PhactoriDbg():
            myDebugPrint3("yes, now " + str(self.mcfndx) + " " + \
                str(self.mMinDistSqrd) + "\n")
      if PhactoriDbg():
        myDebugPrint3("next to be replaced ndx: " + str(self.mcfndx) + \
            " sid: " + \
            str(self.mThisProcIds.GetValue(self.mcfndx)) + \
            " tid: " + \
            str(self.mTargetMatchIds.GetValue(self.mcfndx)) + \
            " dsqrd: " + \
            str(self.mDistSqrds.GetValue(self.mcfndx)) + "\n")

    def ToStr(self):
      retStr = "closest " + str(self.mNumToFind) + " points:\n" + \
      "index: source id, target id, dist sqrd: source xyz: target xyz\n"
      for ii in range(0, self.mNumToFind):
        pp = ii*3
        adst = str(ii) + ": " + \
               str(self.mThisProcIds.GetValue(ii)) + ", " + \
               str(self.mTargetMatchIds.GetValue(ii)) + ", " + \
               str(self.mDistSqrds.GetValue(ii)) + ": " + \
               str(self.mThisProcXyzs.GetValue(pp)) + ", " + \
               str(self.mThisProcXyzs.GetValue(pp+1)) + ", " + \
               str(self.mThisProcXyzs.GetValue(pp+2)) + ": " + \
               str(self.mTargetMatchXyzs.GetValue(pp)) + ", " + \
               str(self.mTargetMatchXyzs.GetValue(pp+1)) + ", " + \
               str(self.mTargetMatchXyzs.GetValue(pp+2)) + "\n"
        retStr += adst
      return retStr


  def FindClosestNPointsToList(self, inGlobalNodeIdList, inPointXyzList,
            inNumToFind):
    """given a list of node ids and xyz points, recursively find the nearest
       (geometrically) inNumToFind points in the local processor to the xyz
       points in inPointXyzList.  Returns a list of inNumToFind global node
       ids and inNumToFind xyz points
       inGlobalNodeIds is vtkIntArray
       inPointXyzs is vtkDoubleArray
       returns a FindClosestNPointsToListParams instance which has list of
       length inNumToFind which contain the node id of the closest local
       process point, the node id of the corresponding point from the
       inGlobalNodeIdList, the xyzs of each of those, and this distance
       between (squared to save computation)"""
    recursionItem = BlockRecursionControlItem()
    recursionItem.mParameters = \
        PhactoriOperationBlock.FindClosestNPointsToListParams(
            self, inNumToFind, inGlobalNodeIdList, inPointXyzList)
    recursionItem.mOperationToDoPerBlock = \
            self.FindClosestNPointsToListInBlock
    self.DoMethodPerBlock(recursionItem)
    return recursionItem.mParameters

  def FindClosestNPointsToListInBlock(self, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("FindClosestNPointsToListInBlock entered\n")
    pointsArray = inInputCsData.GetPoints()
    if pointsArray == None:
      #no points here
      return
    pointsData = inInputCsData.GetPointData()
    globalNodeIdArray = pointsData.GetArray('GlobalNodeId')
    numPoints = pointsArray.GetNumberOfPoints()
    ptXyz = [0.0, 0.0, 0.0]

    #get target to set up Kdtree for quickly finding nearest point
    inParameters.SetUpKdtree()

    for ndx in range(0,numPoints):
      #thePoint = pointsArray.GetPoint(ndx, ptXyz)
      if ndx % 100 == 0:
        if PhactoriDbg():
          myDebugPrint3("test " + str(ndx) + " of " + str(numPoints) + "\n")
      pointsArray.GetPoint(ndx, ptXyz)
      if(globalNodeIdArray == None):
        theGlobalNodeId = ndx + 1
      else:
        theGlobalNodeId = globalNodeIdArray.GetValue(ndx)

      inParameters.TestPointWithKdtree(theGlobalNodeId, ptXyz)
      #inParameters.TestPoint(theGlobalNodeId, ptXyz)
    if PhactoriDbg(100):
      myDebugPrint3("FindClosestNPointsToListInBlock returning\n")

  class MakeListOfAllPoints1Params:
    def __init__(self):
      self.mGlobalNodeIdList = vtk.vtkIntArray()
      self.mPointXYZList = vtk.vtkDoubleArray()

  def MakeListOfAllPoints1(self):
    """recursively going through multiblock setup, make a list of all the
       points in this operation output in this process.  We get a list
       of global node ids and a list of xyz geometry points"""
    recursionItem = BlockRecursionControlItem()
    recursionItem.mParameters = \
        PhactoriOperationBlock.MakeListOfAllPoints1Params()
    recursionItem.mOperationToDoPerBlock = \
            self.MakeListOfAllPointsInBlock1
    self.DoMethodPerBlock(recursionItem)
    return recursionItem.mParameters.mGlobalNodeIdList, \
           recursionItem.mParameters.mPointXYZList

  def MakeListOfAllPointsInBlock1(self, inInputCsData, inParameters):
    #if PhactoriDbg(100):
    #  myDebugPrint3("MakeListOfAllPointsInBlock1 entered\n")
    pointsArray = inInputCsData.GetPoints()
    if pointsArray == None:
      #no points here
      return
    pointsData = inInputCsData.GetPointData()
    globalNodeIdArray = pointsData.GetArray('GlobalNodeId')
    numPoints = pointsArray.GetNumberOfPoints()
    ptXyz = [0.0, 0.0, 0.0]
    for ndx in range(0,pointsArray.GetNumberOfPoints()):
      #thePoint = pointsArray.GetPoint(ndx, ptXyz)
      pointsArray.GetPoint(ndx, ptXyz)
      if globalNodeIdArray == None:
        inParameters.mGlobalNodeIdList.InsertNextValue(ndx+1)
      else:
        inParameters.mGlobalNodeIdList.InsertNextValue(globalNodeIdArray.GetValue(ndx))
      inParameters.mPointXYZList.InsertNextValue(ptXyz[0])
      inParameters.mPointXYZList.InsertNextValue(ptXyz[1])
      inParameters.mPointXYZList.InsertNextValue(ptXyz[2])
    #if PhactoriDbg(100):
    #  myDebugPrint3("MakeListOfAllPointsInBlock1 returning\n")


  class OutputElementListToFileParams:
    def __init__(self):
      self.mOutFileff = None
      self.mBlockCount = 0
      self.mElementCount = 0
      self.mNodeCount = 0
      self.mFlag1 = 0
      self.mKilledByCriteriaCount = [0,0,0,0,0,0,0,0,0,0]

  def OutputElementListToFile(self, inFileNameToWrite):
    recursionItem = BlockRecursionControlItem()

    recursionItem.mParameters = \
      PhactoriOperationBlock.OutputElementListToFileParams()

    if WriteEachDeadCellElementToFiles:
      recursionItem.mParameters.mOutFileff = open(inFileNameToWrite, "w+b")
    recursionItem.mOperationToDoPerBlock = \
            self.OutputElementListFromOneBlockToFile
    self.DoMethodPerBlock(recursionItem)
    if WriteEachDeadCellElementToFiles:
      recursionItem.mParameters.mOutFileff.close()
    return recursionItem.mParameters

  def OutputSingleElementFromBlockToTimeHistoryFile(
          self, inInputCsData, inParameters):
    if PhactoriDbg(100):
      myDebugPrint3("OutputSingleElementFromBlockToTimeHistoryFile entered\n")
    cellData = inInputCsData.GetCellData()
    numCellTuplesX = cellData.GetNumberOfTuples()
    numCellArrays = cellData.GetNumberOfArrays()
    pointData = inInputCsData.GetPointData()
    numPointTuplesX = cellData.GetNumberOfTuples()
    numPointArrays = pointData.GetNumberOfArrays()

    #if inParameters.mBlockCount == 0:
    if (inParameters.mFlag1 == 0) and (numCellTuplesX > 0):
      inParameters.mFlag1 = 1
      #as constructed, we should only open once, so we can do 'w' not 'a+'
      inParameters.mOutElementFileff = open(inParameters.mOutElementFileName, "w")
      inParameters.mOutElementFileff.write("step, simtime")
      #inParameters.mOutElementFileff.write(", element index")
      for arrayIndex in range(0, numCellArrays):
          inParameters.mOutElementFileff.write(
            ", " + str(cellData.GetArray(arrayIndex).GetName()))
      inParameters.mOutElementFileff.write("\n")

      inParameters.mOutNodeFileff = open(inParameters.mOutNodeFileName, "w")
      inParameters.mOutNodeFileff.write("step, simtime")
      for arrayIndex in range(0, numPointArrays):
          oneArray = pointData.GetArray(arrayIndex)
          numComp = oneArray.GetNumberOfComponents()
          oneArrayName = str(oneArray.GetName())
          if numComp == 1:
            inParameters.mOutNodeFileff.write(", " + oneArrayName)
          elif(numComp == 3):
            inParameters.mOutNodeFileff.write(
                    ", " + oneArrayName + "x" + \
                    ", " + oneArrayName + "y" + \
                    ", " + oneArrayName + "z")
          else:
            myDebugPrint3AndException(
              "OutputSingleElementFromBlockToTimeHistoryFile:\n" \
              "expecting 1 or 3 components (A)\n")
      inParameters.mOutNodeFileff.write("\n")

    #for now, only dump 1 element
    if numCellTuplesX > 0:
        numCellsToDo = 1
    else:
        numCellsToDo = 0
    for cellIndex in range(0, numCellsToDo):
        inParameters.mOutElementFileff.write(str(inParameters.mTimeStep) + ", " + \
                str(inParameters.mSimulationTime))
        #inParameters.mOutElementFileff.write(", " + \
        #        str(cellIndex + inParameters.mElementCount))
        for arrayIndex in range(0, numCellArrays):
            inParameters.mOutElementFileff.write(
              ", " + str(cellData.GetArray(arrayIndex).GetTuple1(cellIndex)))
        inParameters.mOutElementFileff.write("\n")

    #for now, only dump 1 node
    if numPointTuplesX > 0:
        numPointsToDo = 1
    else:
        numPointsToDo = 0
    for pointIndex in range(0, numPointsToDo):
        inParameters.mOutNodeFileff.write(str(inParameters.mTimeStep) + ", " + \
                str(inParameters.mSimulationTime))
        #inParameters.mOutNodeFileff.write(", " + \
        #        str(pointIndex + inParameters.mCellCount))
        for arrayIndex in range(0, numPointArrays):
            #inParameters.mOutNodeFileff.write(
            #  ", " + str(pointData.GetArray(arrayIndex).GetTuple1(0)))
            oneArray = pointData.GetArray(arrayIndex)
            numComp = oneArray.GetNumberOfComponents()
            if(numComp == 1):
              inParameters.mOutNodeFileff.write(
                ", " + str(pointData.GetArray(arrayIndex).GetTuple1(pointIndex)))
            elif(numComp == 3):
              theTuple = pointData.GetArray(arrayIndex).GetTuple3(0)
              inParameters.mOutNodeFileff.write( \
                ", " + str(theTuple[0]) + \
                ", " + str(theTuple[1]) + \
                ", " + str(theTuple[2]))
            else:
              myDebugPrint3AndException(
                "OutputSingleElementFromBlockToTimeHistoryFile:\n" \
                "expecting 1 or 3 components (B)\n")
        inParameters.mOutNodeFileff.write("\n")

    inParameters.mNodeCount += inInputCsData.GetNumberOfPoints()
    inParameters.mElementCount += inInputCsData.GetNumberOfCells()
    inParameters.mBlockCount += 1

    if PhactoriDbg(100):
      numCells = inInputCsData.GetNumberOfCells()
      if numCells > 0:
        numPoints = inInputCsData.GetNumberOfPoints()
        myDebugPrint3(
          "numCells: " + str(numCells) + "\n" \
          "numCellTuplesX: " + str(numCellTuplesX) + "\n" \
          "numPoints: " + str(numPoints) + "\n" \
          "numPointTuplesX: " + str(numPointTuplesX) + "\n")

    if(inParameters.mOutElementFileff != None):
      inParameters.mOutElementFileff.flush()
      inParameters.mOutNodeFileff.flush()

    #dump 1 node

    if PhactoriDbg(100):
      myDebugPrint3("OutputSingleElementFromBlockToTimeHistoryFile returning\n")

  class OutputSingleElementToTimeHistoryFileParams:
    def __init__(self):
      self.mOutElementFileff = None
      self.mOutElementFileName = None
      self.mOutNodeFileff = None
      self.mOutNodeFileName = None
      self.mTimeStep = 0
      self.mSimulationTime = 0.0
      self.mBlockCount = 0
      self.mElementCount = 0
      self.mNodeCount = 0
      self.mFlag1 = 0

    def ResetBlockNodeElementCount(self):
      self.mBlockCount = 0
      self.mElementCount = 0
      self.mNodeCount = 0

  def OutputSingleElementToTimeHistoryFile(self,
          inElementFileNameToWrite, inNodeFileNameToWrite,
          inTimeStep, inSimulationTime):
    recursionItem = BlockRecursionControlItem()

    if inElementFileNameToWrite in self.mRecursionParameterStore:
      recursionItem.mParameters = \
        self.mRecursionParameterStore[inElementFileNameToWrite]
      recursionItem.mParameters.ResetBlockNodeElementCount()
    else:
      recursionItem.mParameters = \
        PhactoriOperationBlock.OutputSingleElementToTimeHistoryFileParams()
      self.mRecursionParameterStore[inElementFileNameToWrite] = \
        recursionItem.mParameters

    #recursionItem.mParameters.mOutFileff = open(inFileNameToWrite, "w+b")
    recursionItem.mParameters.mOutElementFileName = inElementFileNameToWrite
    recursionItem.mParameters.mOutNodeFileName = inNodeFileNameToWrite
    recursionItem.mParameters.mTimeStep = inTimeStep
    recursionItem.mParameters.mSimulationTime = inSimulationTime
    recursionItem.mOperationToDoPerBlock = \
            self.OutputSingleElementFromBlockToTimeHistoryFile
    self.DoMethodPerBlock(recursionItem)
    #recursionItem.mParameters.mOutFileff.close()

  def ExportOperationData(self, datadescription):
    """this will be called once per callback (before WriteImages) to allow the
       operation to export any desired data which is not an image.  We call
       the operation specifics version of this method."""
    self.mOperationSpecifics.ExportOperationData(datadescription)

class PhactoriTextAnnotationBlock(PhactoriOperationSpecifics):
  """represents manages one text annotation item (paraview Text source).
     The idea is that we print a 2d text item on the window front to give
     information to the user.  Notionally the text can change depending on
     data (e.g. min/max/average values from an operation), but this is not
     yet implemented.  Imagesets can selectively turn on any desired subset
     of existing text annotations.  The position of the text annotation
     on the window, the color can be set, the font size and family can be
     chosen, and the bold/italic/shadow flags can be set"""
  def __init__(self):
    self.mName = ""
    self.mTextString = ""
    #not implemented yet self.minputOperationName = ""
    #not implemented yet self.minputOperationBlock = None
    self.mParaViewSource = None
    self.mParaViewRepresentation = None
    self.mWindowLocation = None
    self.mPosition = None
    self.mColor = None
    self.mOpacity = None
    self.mBoldFlag = None
    self.mShadowFlag = None
    self.mItalicFlag = None
    self.mFontSize = None
    self.mFontFamily = None

  def UpdateAndMakeVisible(self):
    if PhactoriDbg(100):
      myDebugPrint3(
          "PhactoriTextAnnotationBlock::UpdateAndMakeVisible entered:\n" + \
          self.mName + "\n", 100)

    if self.mParaViewSource == None:
      #this text annotation needs to have its paraview components created
      if PhactoriDbg(100):
        myDebugPrint3("(paraview items will now be created)\n")
      self.CreateParaViewItems()

    #here is where we'd get input operation and update text if it
    #was data dependent; see PhactoriMarkerBlock:UpdateAndMakeVisible()

    #make it visible
    self.mParaViewRepresentation.Visibility = 1

    if PhactoriDbg(100):
      myDebugPrint3(
          "PhactoriTextAnnotationBlock::UpdateAndMakeVisible returning:\n" + \
          self.mName + "\n", 100)

  def MakeInvisible(self):
    if PhactoriDbg(100):
      myDebugPrint3("making text annotation invisible: " + \
              self.mName + "\n", 100)
    if self.mParaViewRepresentation != None:
      self.mParaViewRepresentation.Visibility = 0
    else:
      myDebugPrint3("(paraview items not yet created)\n")

  def CreateParaViewItems(self):
    """creates the paraview text source and representation for this
       PhactoriTextAnnotationBlock instance"""
    if PhactoriDbg(100):
      myDebugPrint3(
          "PhactoriTextAnnotationBlock::CreateParaViewItems entered:\n" + \
          self.mName + "\n", 100)


    savedActiveSource = GetActiveSource()

    #UpdatePipelineWithCurrentTimeArgument(inInputFilter)
    self.mParaViewSource = Text()
    self.mParaViewSource.Text = self.mTextString

    self.mParaViewRepresentation = Show()

    if self.mWindowLocation != None:
      self.mParaViewRepresentation.WindowLocation = self.mWindowLocation
    if self.mPosition != None:
      self.mParaViewRepresentation.Position = self.mPosition
    if self.mFontSize != None:
      self.mParaViewRepresentation.FontSize = self.mFontSize
    if self.mFontFamily != None:
      self.mParaViewRepresentation.FontFamily = self.mFontFamily
    if self.mColor != None:
      self.mParaViewRepresentation.Color = self.mColor
    if self.mOpacity != None:
      self.mParaViewRepresentation.Opacity = self.mOpacity
    if self.mBoldFlag != None:
      self.mParaViewRepresentation.Bold = self.mBoldFlag
    if self.mItalicFlag != None:
      self.mParaViewRepresentation.Italic = self.mItalicFlag
    if self.mShadowFlag != None:
      self.mParaViewRepresentation.Shadow = self.mShadowFlag


    #SetActiveSource(self.mParaViewSource)
    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3(
          "PhactoriTextAnnotationBlock::CreateParaViewItems returning:\n" + \
          self.mName + "\n", 100)

  def  SelfToStr(self):
    return "PhactoriTextAnnotationBlock:\n" \
      "mName: " + self.mName + "\n" \
      "mTextString: " + self.mTextString + "\n" \
      "self.mWindowLocation: " + str(self.mWindowLocation) + "\n" \
      "self.mPosition: " + str(self.mPosition) + "\n" \
      "self.mColor: " + str(self.mColor) + "\n" \
      "self.mOpacity: " + str(self.mOpacity) + "\n" \
      "self.mBoldFlag: " + str(self.mBoldFlag) + "\n" \
      "self.mShadowFlag: " + str(self.mShadowFlag) + "\n" \
      "self.mItalicFlag: " + str(self.mItalicFlag) + "\n" \
      "self.mFontSize: " + str(self.mFontSize) + "\n" \
      "self.mFontFamily: " + str(self.mFontFamily) + "\n"
      #"mInputOperationName: " + str(self.mInputOperationName) + "\n" \

  def ParseTextAnnotationSettingsFromJson(self, inJson):
    if PhactoriDbg(100):
      myDebugPrint3(
        "PhactoriTextAnnotationBlock::ParseTextAnnotationSettingsFromJson " \
                "entered\n", 100)

    if "line1" in inJson:
      self.mTextString = inJson["line1"]
    else:
      myDebugPrint3AndException(
          "PhactoriTextAnnotationBlock::ParseParametersFromJson\n"
          "Error:  must have 'line1' command\n")

    for lineIndex in range(2, 10):
      linetoken = "line" + str(lineIndex)
      if linetoken not in inJson:
        break
      self.mTextString += "\n" + inJson[linetoken]

    if PhactoriDbg(100):
        myDebugPrint3("final text string:\n" + self.mTextString + "\n", 100)

    if "windowlocation" in inJson:
      self.mWindowLocation = inJson["windowlocation"]

    if "position" in inJson:
      self.mWindowLocation = "AnyLocation"
      self.mPosition = inJson["position"]

    if "fontsize" in inJson:
      self.mFontSize = int(inJson["fontsize"])

    if "fontfamily" in inJson:
      self.mFontFamily = inJson["fontfamily"]

    if "color" in inJson:
      self.mColor = inJson["color"]

    if "opacity" in inJson:
      self.mOpacity = float(inJson["opacity"])

    if "boldflag" in inJson:
      self.mBoldFlag = inJson["boldflag"]

    if "italicflag" in inJson:
      self.mItalicFlag = inJson["italicflag"]

    if "shadowflag" in inJson:
      self.mShadowFlag = inJson["shadowflag"]

    if PhactoriDbg(100):
      myDebugPrint3("after parsing:\n" + self.SelfToStr())
      myDebugPrint3(
        "PhactoriTextAnnotationBlock::ParseTextAnnotationSettingsFromJson " \
                "returning\n", 100)


class PhactoriMarkerBlock:
  """represents and manages one 'marker' item.  The idea of a phactori marker
     is that we can place a sphere (or maybe cube or other geometric object)
     at a specific point in the scene to provide some kind of information to
     the user.  Markers can be placed at absolute or relative geometric
     points, at an element or node, or (perhaps most valuably) at a min or
     max of a particular variable (with the min or max determined at the
     pipeline point of any selected operation, e.g. filter out block X
     and place a marker at the max of Temperature on block X).  Imagesets
     can selectively turn on any desired subset of existing markers. In
     Paraview the markers are simply sphere (or maybe cube) sources.
     Additionally, the color of the marker can be set as well as the size
     (relative or absolute)"""

  def __init__(self):
    self.mName = ""
    self.mInputOperationName = None
    self.mInputOperationBlock = None
    self.mShape = "sphere"  #we will add others
    self.mParaViewSource = None
    self.mParaViewRepresentation = None
    self.mLocationPoint = PhactoriUserPointInfo()
    self.mSizeType = "absolute"  #"or datasize relative", not handled yet
    self.mSize = [0.5, 0.5, 0.5]
    self.mOrientation = [0.0, 0.0, 0.0]
    self.mColor = [0.0, 1.0, 0.0]
    self.mResolution = 16

  def UpdateAndMakeVisible(self):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMarkerBlock::UpdateAndMakeVisible entered: " + \
              self.mName + "\n", 100)

    if self.mParaViewSource == None:
      #this marker needs to have its paraview components created
      if PhactoriDbg(100):
        myDebugPrint3("(paraview items will now be created)\n")
      self.CreateParaViewItems()

    #get the current location of the point (which may change for a point
    #at a variable min/max or a relative point)
    global gPipeAndViewsState
    myViewBounds = None
    viewBoundsIo = [myViewBounds]

    if self.mInputOperationBlock == None:
      if self.mInputOperationName == None:
        self.mInputOperationBlock = \
          gPipeAndViewsState.mIncomingDefaultOperation
        if PhactoriDbg(100):
          myDebugPrint3("first calc, using default incoming operation\n")
      else:
        self.mInputOperationBlock = gPipeAndViewsState.mOperationBlocks[
          self.mInputOperationName]
        if self.mInputOperationBlock == None:
          myDebugPrint3AndException(
            "PhactoriMarkerBlock::UpdateAndMakeVisible:\n"
            "couldn't find input operation with name: " + \
            str(self.mInputOperationName) + "\n")
        if PhactoriDbg(100):
            myDebugPrint3("first calc, using operation named: " + \
                    self.mInputOperationName + "\n")

    checkPvFilter = self.mInputOperationBlock.GetPvFilter()

    updatedMarkerLocation = self.mLocationPoint.\
      GetCurrentGeometricPointWithDisplacement(checkPvFilter,
            viewBoundsIo, True)

    #update the paraview component to reflect the current marker location
    #self.mParaViewSource.Center = updatedMarkerLocation
    self.mParaViewRepresentation.Position = updatedMarkerLocation
    if PhactoriDbg(100):
      #myDebugPrint3("updatedMarkerLocation: " + \
      #        str(self.mParaViewSource.Center)+ "\n", 100)
      myDebugPrint3("updatedMarkerLocation: " + \
              str(self.mParaViewRepresentation.Position)+ "\n", 100)

    #make it visible
    self.mParaViewRepresentation.Visibility = 1

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriMarkerBlock::UpdateAndMakeVisible returning: " + \
              self.mName + "\n", 100)

  def MakeInvisible(self):
    if PhactoriDbg(100):
      myDebugPrint3("making marker invisible: " + self.mName + "\n", 100)
    if self.mParaViewRepresentation != None:
      self.mParaViewRepresentation.Visibility = 0
    else:
      if PhactoriDbg(100):
        myDebugPrint3("(paraview items not yet created\n")

  def CreateParaViewItems(self):
    """creates the paraview source (e.g box, sphere, cone, arrow) and representation
       for the marker"""
    if PhactoriDbg(100):
        myDebugPrint3("PhactoriMarkerBlock::CreateParaViewItems entered\n", 100)
    savedActiveSource = GetActiveSource()
    if self.mShape == "sphere":
      self.mParaViewSource = Sphere()
      #self.mParaViewRepresentation.SelectionPointFieldDataArrayName = 'Normals'
      self.mParaViewSource.Radius = 0.5
      self.mParaViewSource.ThetaResolution = self.mResolution
      self.mParaViewSource.PhiResolution = self.mResolution
    elif self.mShape == "box":
      self.mParaViewSource = Box()
      self.mParaViewRepresentation = Show()
      #self.mParaViewSource.ZLength = self.mSize[0]
      #self.mParaViewSource.XLength = self.mSize[1]
      #self.mParaViewSource.YLength = self.mSize[2]
    elif self.mShape == "arrow":
      self.mParaViewSource = Arrow()
      self.mParaViewSource.Invert = 1
      self.mParaViewSource.ShaftResolution = self.mResolution
      self.mParaViewSource.TipResolution = self.mResolution
    elif self.mShape == "cone":
      self.mParaViewSource = Cone()
      self.mParaViewSource.Direction = [-1.0, 0.0, 0.0]
      self.mParaViewSource.Center = [0.5, 0.0, 0.0]
      self.mParaViewSource.Resolution = self.mResolution
      self.mParaViewSource.Capping = 1

    self.mParaViewRepresentation = Show()

    self.mParaViewRepresentation.Position = [0.0, 0.0, 0.0]
    self.mParaViewRepresentation.Scale = self.mSize
    self.mParaViewRepresentation.Orientation = self.mOrientation

    if gParaViewCatalystVersionFlag <= 40100:
      self.mParaViewRepresentation.\
              AllowSpecularHighlightingWithScalarColoring = 1
    self.mParaViewRepresentation.ScaleFactor = 0.1
    self.mParaViewRepresentation.EdgeColor = [0.0, 0.0, 0.5000076295109483]
    self.mParaViewRepresentation.DiffuseColor = self.mColor

    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
        myDebugPrint3("PhactoriMarkerBlock::CreateParaViewItems returning\n", 100)


  def SelfToStr(self):
    return "marker:\n" \
      "mName: " + self.mName + "\n" \
      "mInputOperationName: " + str(self.mInputOperationName) + "\n" \
      "mShape: " + self.mShape + "\n" \
      "mSizeType: " + self.mSizeType + "\n" \
      "mSize: " + str(self.mSize) + "\n" \
      "mOrientation: " + str(self.mOrientation) + "\n" \
      "mColor: " + str(self.mColor) + "\n" \
      "mResolution: " + str(self.mResolution) + "\n" \
      "mLocationPoint:\n" + self.mLocationPoint.SelfToStr()

  def ParseMarkerSettingsFromJson(self, inJsn):
    if PhactoriDbg(100):
      myDebugPrint3(
        "PhactoriMarkerBlock::ParseMarkerSettingsFromJson entered\n", 100)

    if 'input' in inJsn:
      self.mInputOperationName = inJsn['input']
    else:
      noticeStr = \
          'notice!  inJsn has no input key, using default pipeline input\n'
      self.mInputOperationName = None
      if PhactoriDbg():
        myDebugPrint3(noticeStr)

    if "color" in inJsn:
      self.mColor = inJsn["color"]
      if len(self.mColor) != 3:
          myDebugPrint3AndException(
            "'color' item was not three element list")
    if "relative size" in inJsn:
      myDebugPrint3AndException(
        "PhactoriMarkerBlock::ParseMarkerSettingsFromJson:\n"
        "'relative size' not yet implemented")
    if "absolute size" in inJsn:
      self.mSize = inJsn["absolute size"]
    if "orientation" in inJsn:
      self.mOrientation= inJsn["orientation"]
    if "resolution" in inJsn:
      self.mResolution = inJsn["resolution"]
    if "shape" in inJsn:
      self.mShape = inJsn["shape"]
      if self.mShape == "sphere":
        if PhactoriDbg(100):
          myDebugPrint3("shape is sphere\n", 100)
      elif self.mShape == "box":
        if PhactoriDbg(100):
          myDebugPrint3("shape is box\n", 100)
      elif self.mShape == "cone":
        if PhactoriDbg(100):
          myDebugPrint3("shape is cone\n", 100)
      elif self.mShape == "arrow":
        if PhactoriDbg(100):
          myDebugPrint3("shape is arrow\n", 100)
      else:
        myDebugPrint3AndException(
          "PhactoriMarkerBlock::ParseMarkerSettingsFromJson:\n"
          "'shape' must be 'sphere' or 'box' or 'cone' or 'arrow'\n"
          "got: " + str(self.mShape) + "\n")

    #get point location for marker settings, e.g. marker at max variable point
    self.mLocationPoint.UserPointInfoParseParametersFromJson(inJsn,
      "marker at ", "")

    if self.mLocationPoint.mParsingDetectedAtLeastOneSetting == False:
      myDebugPrint3AndException(
        "PhactoriMarkerBlock::ParseMarkerSettingsFromJson:\n"
        "must have 'marker at ... point' parameters")

    if PhactoriDbg(100):
      myDebugPrint3(self.SelfToStr())
      myDebugPrint3(
        "PhactoriMarkerBlock::ParseMarkerSettingsFromJson returning\n", 100)



#class ParaViewRenderInfoC:
#  def __init__(self):
#    self.mParaViewInfo = None
#    self.mLookDirection = None

class PhactoriAnnotationPv:
  '''class manages the ParaView stuff for annotation.  The settings are
     managed in PhactoriAnnotationViewSettings, and instances of this
     class are used in PhactoriImageset to manage the actual paraview
     objects, which might be different from image set to image set,
     while still using the same PhactoriAnnotationViewSettings'''
  def __init__(self, inSettings):
    #paraview stuff not created yet
    self.mParaViewTextAnnotationSource = None
    self.mParaViewRepresentation = None
    self.mInputOperation = None
    self.mSettings = inSettings

  def CreateParaViewStuff(self, inColor, inImageset = None):
    '''assumes inImageset is valid; creates ParaView source and
       representation based on internal parameters'''
    if PhactoriDbg():
      myDebugPrint3("PhactoriAnnotationPv::CreateParaViewStuff entered\n")
    if self.mSettings.mType == 'time':
      if inImageset == None:
        if PhactoriDbg():
          myDebugPrint3("PhactoriAnnotationPv::CreateParaViewStuff\n" + \
              "this is a time annotation, and it needs a valid imageset\n" + \
              "to work off, and None was supplied\n")
        return
      self.mInputOperation = inImageset.mOperation

      savedActiveSource = GetActiveSource()
      inOprtn = self.mInputOperation
      if(inOprtn.mParaViewTimeAnnotationSource == None):
        SetActiveSource(inOprtn.GetPvFilter())
        inOprtn.mParaViewTimeAnnotationSource = AnnotateTime()
        #inOprtn.mParaViewTimeAnnotationSource.Format = 'Time: %.20e'
        inOprtn.mParaViewTimeAnnotationSource.Format = self.mSettings.mTimeFormatString
      else:
        SetActiveSource(inOprtn.mParaViewTimeAnnotationSource)
      #self.mParaViewAnnotationSource = Text()
      self.mParaViewRepresentation = Show()
      self.mParaViewRepresentation.Color = inColor
      SetActiveSource(savedActiveSource)
    else:
      self.mParaViewTextAnnotationSource = Text()
      self.mParaViewTextAnnotationSource.Text = self.mText
      self.mParaViewRepresentation = Show()
      self.mParaViewRepresentation.Color = inColor
    self.mParaViewRepresentation.FontSize = self.mSettings.mFontSize
    self.mParaViewRepresentation.WindowLocation = self.mSettings.mWindowLocation
    if PhactoriDbg():
      myDebugPrint3("window location: " +  str(self.mParaViewRepresentation.WindowLocation) + "\n")
    if PhactoriDbg():
      myDebugPrint3("PhactoriAnnotationPv::CreateParaViewStuff returning\n")

#global gTestTimeNotation
#gTestTimeNotation = 1

class PhactoriAnnotationViewSettings:
  '''class manages settings for how a text or time annotation is to be
     show in the window--size, position, alignment, visible/invisible,
     etc.  An instance is a member of a PhactoriRepresentation.  ParaView
     stuff which uses these settings is in a different class, which is
     associated with ImageSet.   Gives access to control annotation
     position, size, text, and could be extended to control all the items
     controllable in paraview (x,y position, alignment, bold/italic/shadow'''

  def __init__(self, inName, inAnnotationType):
    if PhactoriDbg():
      myDebugPrint3("ViewSettingsPhactoriAnnotation::__init__ entered\n")
    self.mName = inName
    self.mVisible = False
    if inAnnotationType == 'time':
      self.mType = 'time'
      self.mText = None
      self.mTimeFormatString = gDefaultTimeFormatString
    else:
      self.mType = 'text'
      self.mText = 'Text'
    self.mSizeScale = 1.0
    self.mFontSize = 18
    self.mFont = 'Arial'
    self.mBold = False
    self.mItalic = False
    self.mShadow = False
    self.mOpacity = 1.0
    self.mAlign = 'Center'
    self.mWindowLocation = 'UpperLeftCorner'
    if PhactoriDbg():
      myDebugPrint3("ViewSettingsPhactoriAnnotation::__init__ returning\n")

  def ParseAvsFromJson(self, inJsn):
    self.mVisible = getParameterFromBlock(inJsn, 'show time annotation',
        self.mVisible)

    #global gTestTimeNotation
    #if gTestTimeNotation == 0:
    #  self.mVisible = True
    #  gTestTimeNotation = 1
    #else:
    #  self.mVisible = False
    #  gTestTimeNotation = 0

    if 'time annotation position' in inJsn:
      windowPosSize = inJsn['time annotation position']
      winPos = windowPosSize[0]
      winSize = windowPosSize[1]
      self.SetSizeScale(winSize)
      if winPos == 'top left':
        self.SetWindowLocation('UpperLeftCorner')
      elif winPos == 'top right':
        self.SetWindowLocation('UpperRightCorner')
      elif winPos == 'top':
        self.SetWindowLocation('UpperCenter')
      elif winPos == 'bottom left':
        self.SetWindowLocation('LowerLeftCorner')
      elif winPos == 'bottom right':
        self.SetWindowLocation('LowerRightCorner')
      elif winPos == 'bottom':
        self.SetWindowLocation('LowerCenter')
      else:
        if PhactoriDbg():
          myDebugPrint3("bad time annotation position, using upper left\n")
        self.SetWindowLocation('UpperLeftCorner')

    if 'time format string' in inJsn:
      self.mTimeFormatString = inJsn['time format string']

  def SetText(self, newText):
    if self.mType != 'text':
      if PhactoriDbg():
        myDebugPrint3("ViewSettingsPhactoriAnnotationSource::SetText not text annotation\n")
      return
    self.mText = newText
    self.mParaViewTextAnnotationSource.Text = newText

  def SetSizeScale(self, inNewSizeScale):
    if inNewSizeScale < 0.001:
      if PhactoriDbg():
        myDebugPrint3(
            "ViewSettingsPhactoriAnnotationSource::SetSizeScale not in [0.001,1000]\n")
    if inNewSizeScale > 1000.0:
      if PhactoriDbg():
        myDebugPrint3(
            "ViewSettingsPhactoriAnnotationSource::SetSizeScale not in [0.001,1000]\n")
    self.mSizeScale = inNewSizeScale
    fontSizeFloat = float(self.mFontSize)
    fontSizeFloat *= self.mSizeScale
    intSizeFloat = math.floor(fontSizeFloat)
    self.mFontSize = int(intSizeFloat)
    #self.mParaViewRepresentation.FontSize = self.mFontSize

  def SetWindowLocation(self, inNewPosition):
    validPositions = ['UpperLeftCorner', 'UpperRightCorner',
        'LowerLeftCorner', 'LowerRightCorner',
        'UpperCenter', 'LowerCenter']
    if inNewPosition in validPositions:
      #self.mParaViewRepresentation.WindowLocation = inNewPosition
      self.mWindowLocation = inNewPosition
    else:
      if PhactoriDbg():
        myDebugPrint3("ViewSettingsPhactoriAnnotationSource::SetWindowLocation invalid " + \
            "position\n")

class PhactoriPolyDataForPlotLine:
  """contains time over plot data values and vtk geometry, as well as the
  routines to construct and maintain the vtk geometry
  """
  def __init__(self, inName):
    self.m_TableColumn = vtk.vtkDoubleArray()
    self.m_TableColumn.SetName(inName)

  def GetRestartInfo(self):
    """construct, in python map/json format, the information from this
       PhactoriPolyDataForPlotLine instance which contains the information
       which would be needed to restore the instance to the proper
       state after a simulation restart, particularly prior data values.
       Return the restart info map/json"""
    newJsonItem = {}
    numValues = self.m_TableColumn.GetNumberOfValues()
    if numValues > 0:
      newJsonItem["m_HasAtLeastOnePoint"] = True
    else:
      newJsonItem["m_HasAtLeastOnePoint"] = False
    if self.m_HasAtLeastOnePoint:
      #should we not make copy and just return reference?
      #we don't need deepcopy because all items are floats/non-objects
      valueList = []
      for ii in range(0, self.m_TableColumn.GetNumberOfValues()):
        valueList.append(self.m_TableColumn.GetValue(ii))
      newJsonItem["m_ColumnValues"] = valueList
    return newJsonItem

  def SetFromRestartInfo(self, inJson):
    """given a map (json format), use the info in the map to set this
       PhactoriPolyDataForPlotLine instance state--this reads the info
       created by out in GetRestartInfo.  It should only occur on a restart,
       and it fills in the data that had been generated during the
       earlier run of the system"""
    if "m_HasAtLeastOnePoint"not in inJson:
      if PhactoriDbg():
        myDebugPrint3("PhactoriPolyDataForPlotLine::SetFromRestartInfo " + \
            "m_HasAtLeastOnePoint not in inJson\n")
      return

    hasAtLeastOnePoint = inJson["m_HasAtLeastOnePoint"]
    if hasAtLeastOnePoint == False:
      return

    inJsonIsGood = True
    if "m_ColumnValues" not in inJson:
      myDebugPrint("PhactoriPolyDataForPlotLine::SetFromRestartInfo " + \
          "m_ColumnValues not in inJson\n")
      inJsonIsGood = False
    if inJsonIsGood == False:
      return

    valueList = inJson["m_ColumnValues"]
    self.m_ColumnValues = valueList[:]
    for vv in valueList:
      self.m_TableColumn.InsertNextValue(vv)

  def AppendPlotValue(self, inNewValue):
    self.m_TableColumn.InsertNextValue(inNewValue)


class PhactoriPlotOverTimeIdLine:
  """useful info for tracking a single element over time in a plot line"""
  def __init__(self, inIdToPlot):
    self.m_Id = inIdToPlot
    linename = "Id " + str(inIdToPlot)
    self.m_PlotColumn = PhactoriPolyDataForPlotLine(linename)
    self.m_FoundOnThisProcessorLastTime = False
    self.m_FoundBlockRef = None
    self.m_FoundIndex = None

class ImageFileNameCountSettings:
  """used to control par of how the image filenames are generated. The default
     is to increment a counter each time the results output block is triggered
     to call the catalyst stuff and make the count value part of the filename.
     We can optionally have a filename section that is based on the date/time,
     we can optionally include microseconds, and we can optionally convert the
     date/time to an integer which will be monotonically increasing"""
  def __init__(self):
    self.mUseDateTime = False
    self.mUseMicrosecondsWithDateTime = False
    self.mConvertDateTimeToInteger = True
    self.mUseSimulationTime = False
    #self.mSimulationTimeFormat = "simtime_%e"
    self.mSimulationTimeFormat = "simtime_%.6e"
    self.mUseCallbackCount = True
    self.mUseRemeshRestartTag = True

  def ParseImageFileNameCountSettings(self, inJson):
    self.mUseDateTime = getParameterFromBlock(inJson,
      "filename use datetime", self.mUseDateTime)
    self.mUseMicrosecondsWithDateTime = getParameterFromBlock(inJson,
      "filename datetime microseconds", self.mUseMicrosecondsWithDateTime)
    self.mConvertDateTimeToInteger = getParameterFromBlock(inJson,
      "filename datetime convert to integer", self.mConvertDateTimeToInteger)
    self.mUseSimulationTime = getParameterFromBlock(inJson,
      "filename use simulation time", self.mUseSimulationTime)
    self.mSimulationTimeFormat = getParameterFromBlock(inJson,
      "filename simulation time format", self.mSimulationTimeFormat)
    self.mUseCallbackCount = getParameterFromBlock(inJson,
      "filename use call count", self.mUseCallbackCount)
    self.mUseRemeshRestartTag = getParameterFromBlock(inJson,
      "filename use restart or remesh tag", self.mUseRemeshRestartTag)

  def GetImageFilename(self, datadescription, inImageSettings,
    inOneLookDirectionFilenameAddon, inRepresentationFilenameAddon,
    inCameraFilenameAddon):
    """given the imageset and the look direction filename addon, figure out the
       filename for the image.  Due to prior implementation, we return two
       values, a temporary name based only on the callback count, and the
       entire name we intend to have so that all processes can agree on the
       file name when it is generated and then process zero can move the image
       filename to it's final correct value.  If we return a None as the
       second value, this indicates the callback count temporary name is also
       the final name so no operating system move needs to happen"""

    inImageBasename = inImageSettings.mImageBasename
    inImageBasedirectory = inImageSettings.mImageBasedirectory
    inNumCounterDigits = inImageSettings.mNumCounterDigits
    inImageFormat = inImageSettings.mImageFormat

    timestep = datadescription.GetTimeStep()

    lastImagesetName = ""

    unusedFilenameFound = False
    extraCount = 0
    import os

    while unusedFilenameFound == False:
      imageBasename = inImageBasename
      imageBasedirectory = inImageBasedirectory
      nameSuffix = inRepresentationFilenameAddon + \
        inCameraFilenameAddon + inOneLookDirectionFilenameAddon

      if imageBasedirectory == None:
        fileBaseName = imageBasename + nameSuffix
      elif imageBasedirectory == "":
        fileBaseName = imageBasename + nameSuffix
      else:
        import os
        fileBaseName = imageBasedirectory + os.sep + imageBasename + \
            nameSuffix

      fname = fileBaseName

      fname += "%t." + inImageFormat

      extraString = GetExtraFileString()
      if extraString != '':
        fname = fname.replace("%t", extraString + "-%t")
      timestepString = str(timestep + extraCount)
      while len(timestepString) < inNumCounterDigits:
        timestepString = "0" + timestepString
      myDebugPrint2("image fname: " + fname + "\ndigit count: " + str(inNumCounterDigits) + "\ntimestepString: " + timestepString + "\n")

      rplstr1 = None
      if self.mUseDateTime:
        ttnow = gPipeAndViewsState.mCallbackDateTime
        if self.mUseMicrosecondsWithDateTime == False:
          ttnow = ttnow.replace(microsecond=0)
        rplstr1 = str(ttnow.isoformat('-'))
        if self.mConvertDateTimeToInteger:
          import re
          rplstr1 = re.sub('[-:.]', '', rplstr1)

      if self.mUseSimulationTime:
        global gSharedRenderView
        simtm = gSharedRenderView.ViewTime
        simtmstr = self.mSimulationTimeFormat % simtm
        if rplstr1 == None:
          rplstr1 = simtmstr
        else:
          rplstr1 = rplstr1 + "_" + simtmstr

      if self.mUseCallbackCount:
        #if rplstr1 is None, we just fall through and use the default
        #behavior; otherwise we need to add to the string
        if rplstr1 != None:
          rplstr1 = timestepString + "_" + rplstr1

      if rplstr1 != None:
        fnameRR = fname.replace("%t", rplstr1)
      else:
        fnameRR = None

      fname = fname.replace("%t", timestepString)

      myDebugPrint2("fname after replace: ->" + fname + "<-\n")
      #if os.path.exists(fname):
      if False:  #don't do os.path.exists; it's maybe slow
        #myDebugPrint2("  filename in use, try to find another\n")
        #unusedFilenameFound = False
        unusedFilenameFound = True
        extraCount += 1
      else:
        #myDebugPrint2("  filename is not in use, use it\n")
        unusedFilenameFound = True

    if self.mUseRemeshRestartTag == False:
      if fnameRR != None:
        import re
        #print "removing, e.g. -s0004-"
        #print "before:", fnameRR
        fnameRR = re.sub('-s[0123456789]*-', '', fnameRR)
        #print "after:", fnameRR

    return fname, fnameRR


class PhactoriPlot1Base:
  def __init__(self):
    self.mName = ""
    self.m_DataCubeAxesInfo = PhactoriDataCubeAxesInfo()

    self.m_PlotType = "PhactoriPlot1Base"
    self.mImageSettings = PhactoriImageSettings()
    self.mImageSettings.mPixelBorderRatioXY = [0.175 / self.mImageSettings.GetAspectRatio(), 0.175]
    self.mInputOperation = None
    #self.m_plotColorInfo = PlotColorInfo()
    self.mColorSettings = PhactoriColorSettings(inUsePlotDefaultColors = True)

    #used to help generate paraview View and Representation
    self.mCamera = PhactoriCameraBlock()

    #paraview filter created to hold geometric data
    self.m_producer = None
    self.m_PhactoriRenderViewInfo = None

    #reference for paraview RenderView instance which is shared among all
    #renderings
    self.mSharedPvRenderView2 = None

    #reference for paraview DataRepresentation instance for this imageset.
    #if the imageset has multiple view directions, they can share the same
    #representation
    self.mPvDataRepresentation2 = None

    #for compatibility with PhactoriImagesetBlock
    self.mLookDirectionList = []
    self.mLookDirectionFilenameAddon = []
    self.mVisibleOps = [None]
    self.mVisibleReps = [None]
    self.mVisiblePvDataReps = [None]
    self.mColorLegendRepRefs = [None]

    #used to control if we want to have image names sequentially count by
    #callback index, or if we want to use a date/time based count
    #(which can use microseconds or not, and can be converted to an integer
    #or not)
    self.mImageFileNameCountSettings = ImageFileNameCountSettings()

  def GetInputPhactoriOperation(self):
    #mainly for compatibility between plot blocks and imageset blocks, give
    #same access to phactori operation coming being visualized
    return self.mInputOperation

  def WriteImages(self, datadescription):
    """write out the .png/.jpg/whatever images associated with this plot
       block for the current timestep/state.  Must loop through camera angles
       and do a write for each one if necessary"""
    if PhactoriDbg():
      myDebugPrint3("PhactoriPlotBase1::WriteImages entered\n")

    if self.mImageSettings.mUsingDefaultGeneratedImageBasename:
      if self.m_PlotType == "PhactoriScatterPlotBlock":
        oneLookDirectionFilenameAddon = "sctr."
      elif self.m_PlotType == "PhactoriPlotOverTimeBlock":
        oneLookDirectionFilenameAddon = "plot."
      else:
        oneLookDirectionFilenameAddon = ""
    else:
      oneLookDirectionFilenameAddon = ""

    fname, fnameRR = self.mImageFileNameCountSettings.GetImageFilename(
        datadescription, self.mImageSettings,
        oneLookDirectionFilenameAddon,
        "", #no representation string for plot
        self.mCamera.mFilenameAddon)

    #used to do this:
    #view.ViewTime = datadescription.GetTime()
    #maybe need to do this?
    #UpdatePipelineWithCurrentTimeArgument()

    #SetUpViewAndRepresentationBeforeWriteImage(oneViewInfo)
    self.SetUpViewAndRepresentationBeforeWriteImage()

    #only need to update color range for first look direction, rest
    #are same
    #if self.mRepresentation.mUseFixedColorRange == False:
    #  UseDataRangeForColorValues(self.mPvDataRepresentation2,
    #      self.mRepresentation, self.mOperation)
    if PhactoriDbg(150):
      myDebugPrint3("writing plot image: " + fname + "\n" +
        str(self.mSharedPvRenderView2) + "\n")

    self.mSharedPvRenderView2.LeftAxisUseCustomRange = 0
    self.mSharedPvRenderView2.BottomAxisUseCustomRange = 0

    if self.m_PlotType == "PhactoriPlotOverTimeBlock":
      self.mPvDataRepresentation2.AttributeType = "Point Data"
      self.mPvDataRepresentation2.AttributeType = "Row Data"
      self.mSharedPvRenderView2.LeftAxisTitle = self.m_YAxisVariableInfo.GetXYPlotAxisTitle()
      self.mSharedPvRenderView2.BottomAxisTitle = "Time"
      self.mSharedPvRenderView2.ChartTitle = ""
      self.mSharedPvRenderView2.ShowLegend = 1
    else:
      self.m_producer.UpdatePipeline()
      DataRepresentation1 = self.mPvDataRepresentation2
      atrtype = "Point Data"
      if self.m_YAxisVariableInfo.mVariableType == 'element':
        atrtype = "Cell Data"
      DataRepresentation1.AttributeType = "Row Data"
      DataRepresentation1.AttributeType = atrtype
      self.mSharedPvRenderView2.LeftAxisTitle = self.m_YAxisVariableInfo.GetXYPlotAxisTitle()
      self.mSharedPvRenderView2.BottomAxisTitle = self.m_XAxisVariableInfo.GetXYPlotAxisTitle()
      #self.mSharedPvRenderView2.ChartTitle = "Scatter Plot at Time ${TIME}"
      self.mSharedPvRenderView2.ChartTitle = ""
      self.mSharedPvRenderView2.ShowLegend = 0
      if PhactoriDbg(150):
        myDebugPrint3("B representation made for scatter plot:\n" + str(DataRepresentation1) + "\n")
        myDebugPrint3("m_producer for scatter plot:\n" + str(self.m_producer) + "\n")
        myDebugPrint3("rep input num points " + str(DataRepresentation1.Input.GetDataInformation().DataInformation.GetNumberOfPoints()) + "\n")
        myDebugPrint3("DataRepresentation1.ListProperties() " + str(DataRepresentation1.ListProperties()) + "\n")
        myDebugPrint3("DataRepresentation1.Input " + str(DataRepresentation1.Input) + "\n")
        myDebugPrint3("DataRepresentation1.AttributeType " + str(DataRepresentation1.AttributeType) + "\n")
        myDebugPrint3("DataRepresentation1.UseIndexForXAxis " + str(DataRepresentation1.UseIndexForXAxis) + "\n")
        myDebugPrint3("DataRepresentation1.Visibility " + str(DataRepresentation1.Visibility) + "\n")
        myDebugPrint3("DataRepresentation1.XArrayName " + str(DataRepresentation1.XArrayName) + "\n")
        myDebugPrint3("DataRepresentation1.SeriesVisibility " + str(DataRepresentation1.SeriesVisibility) + "\n")
        myDebugPrint3("DataRepresentation1.SeriesLabel " + str(DataRepresentation1.SeriesLabel) + "\n")
        myDebugPrint3("DataRepresentation1.SeriesColor " + str(DataRepresentation1.SeriesColor) + "\n")
        myDebugPrint3("DataRepresentation1.SeriesPlotCorner " + str(DataRepresentation1.SeriesPlotCorner) + "\n")
        myDebugPrint3("DataRepresentation1.SeriesLabelPrefix " + str(DataRepresentation1.SeriesLabelPrefix) + "\n")
        myDebugPrint3("DataRepresentation1.SeriesLineStyle " + str(DataRepresentation1.SeriesLineStyle) + "\n")
        myDebugPrint3("DataRepresentation1.SeriesLineThickness " + str(DataRepresentation1.SeriesLineThickness) + "\n")
        myDebugPrint3("DataRepresentation1.SeriesMarkerStyle " + str(DataRepresentation1.SeriesMarkerStyle) + "\n")

    yaxisstng = self.m_xyzMinMaxTrkC.mXyzTrk[1]
    if yaxisstng.mUseLowestBot and yaxisstng.mUseHighestBot and yaxisstng.mUseLowestTop and yaxisstng.mUseHighestTop:
      self.mSharedPvRenderView2.LeftAxisUseCustomRange = 1
      self.mSharedPvRenderView2.LeftAxisRangeMinimum = yaxisstng.mLowestBot
      self.mSharedPvRenderView2.LeftAxisRangeMaximum = yaxisstng.mHighestTop
    xaxisstng = self.m_xyzMinMaxTrkC.mXyzTrk[0]
    if xaxisstng.mUseLowestBot and xaxisstng.mUseHighestBot and xaxisstng.mUseLowestTop and xaxisstng.mUseHighestTop:
      self.mSharedPvRenderView2.BottomAxisUseCustomRange = 1
      self.mSharedPvRenderView2.BottomAxisRangeMinimum = xaxisstng.mLowestBot
      self.mSharedPvRenderView2.BottomAxisRangeMaximum = xaxisstng.mHighestTop

    WriteImage(fname, self.mSharedPvRenderView2,
        Magnification=1)
    #handle datetime naming extra work to avoid race condition
    if fnameRR != None:
      if SmartGetLocalProcessId() == 0:
        import os
        os.rename(fname, fnameRR)

    #ClearPvViewAndPvRepAfterWriteImage(oneViewInfo)
    self.ClearPvViewAndPvRepAfterWriteImage()

  def SetUpViewAndRepresentationBeforeWriteImage(self):
    """Since it turns out to be much more memory efficient to reuse a single
       render view rather than having one per image, this routine is called
       from PhactoriDriver immediately before WriteImage in order to set the
       camera, background color, image size, etc. up appropriately for the
       WriteImage call"""
    if PhactoriDbg(150):
      myDebugPrint3("PhactoriPlot1Base:" \
        "SetUpViewAndRepresentationBeforeWriteImage entered " + \
         str(self.mName) + "\n", 150)
    self.mPvDataRepresentation2.Visibility = 1

    pvRenderView = self.mSharedPvRenderView2
    pvDataRep = self.mPvDataRepresentation2

    #image size
    #pvRenderView.ViewSize = self.mImageSettings.mImageSize
    pvRenderView.ViewSize = [int(self.mImageSettings.mImageSize[0]),
                             int(self.mImageSettings.mImageSize[1])]

    UpdatePlotViewLook(self)

  def ClearPvViewAndPvRepAfterWriteImage(self):
    """Since it turns out to be much more memory efficient to reuse a single
       render view rather than having one per image, this routine is called
       immediately after WriteImage in order to make stuff invisible again
       before the next item gets a chance to do WriteImage"""
    #if PhactoriDbg(150):
    #  myDebugPrint3("ClearPvViewAndPvRepAfterWriteImage entered\n", 150)

    #3d/pointset dataset invisible (plot or 3d viewing)
    self.mPvDataRepresentation2.Visibility = 0


class PhactoriPlotOverTimeBlock(PhactoriPlot1Base):
  """container for one time plot of variable
  """
  def __init__(self):
    PhactoriPlot1Base.__init__(self)
    self.m_PlotType = "PhactoriPlotOverTimeBlock"
    self.m_xyzMinMaxTrkC = PlotXYZMinMaxTrkC()

    self.m_XAxisVariableInfo = PhactoriVariableInfo()
    self.m_XAxisVariableInfo.mVariableName = "Time Step"

    self.m_YAxisVariableInfo = PhactoriVariableInfo()
    self.m_YAxisVariableInfo.mVariableName = None

    self.m_NumberOfEntries = 0
    self.m_vtkTable = None
    self.m_TimeColumn = PhactoriPolyDataForPlotLine("Time")
    self.m_MaxPlotLine = PhactoriPolyDataForPlotLine("Max")
    self.m_MinPlotLine = PhactoriPolyDataForPlotLine("Min")
    self.m_MeanPlotLine = PhactoriPolyDataForPlotLine("Mean")
    self.m_IdPlotLineList = []

    self.mPlotMaximumFlag = None
    self.mPlotMinimumFlag = None
    self.mPlotMeanFlag = None

  def SetFromRestartInfo(self, inJson):
    """given a map (json format), use the info in the map to set the
       plot state--this reads the info created by out in
       GetRestartInfo.  It should only occur on a restart, and
       it fills in the data that had been generated during the
       earlier run of the system"""
    if 'm_xyzMinMaxTrkC' not in inJson:
      if PhactoriDbg():
        myDebugPrint3("PhactoriPlotOverTimeBlock::SetFromRestartInfo " + \
            "m_xyzMinMaxTrkC not in inJson\n")
      return

    self.m_xyzMinMaxTrkC.SetFromRestartInfo(inJson["m_xyzMinMaxTrkC"])
    self.m_TimeColumn.SetFromRestartInfo(inJson["m_TimeColumn"])
    self.m_MaxPlotLine.SetFromRestartInfo(inJson["m_MaxPlotLine"])
    self.m_MinPlotLine.SetFromRestartInfo(inJson["m_MinPlotLine"])
    self.m_MeanPlotLine.SetFromRestartInfo(inJson["m_MeanPlotLine"])


  def GetRestartInfo(self):
    """construct, in python map/json format, the information from this
       plot over time instance which contains the information
       which would be needed to restore the plot over time to the proper
       state after a simulation restart, particularly prior plot points.
       Return the restart info map/json"""
    newJsonItem = {}
    newJsonItem["m_xyzMinMaxTrkC"] = self.m_xyzMinMaxTrkC.GetRestartInfo()
    newJsonItem["m_TimeColumn"] = self.m_TimeColumn.GetRestartInfo()
    newJsonItem["m_MaxPlotLine"] = self.m_MaxPlotLine.GetRestartInfo()
    newJsonItem["m_MinPlotLine"] = self.m_MinPlotLine.GetRestartInfo()
    newJsonItem["m_MeanPlotLine"] = self.m_MeanPlotLine.GetRestartInfo()
    return newJsonItem

    #self.m_TimeColumn = PhactoriPolyDataForPlotLine()
    #self.m_MaxPlotLine = PhactoriPolyDataForPlotLine()
    #self.m_MinPlotLine = PhactoriPolyDataForPlotLine()
    #self.m_MeanPlotLine = PhactoriPolyDataForPlotLine()

class PhactoriScatterPlotBlock(PhactoriPlot1Base):
  """container for one parallel scatterplot of one variable
  """
  def __init__(self):
    PhactoriPlot1Base.__init__(self)
    self.m_PlotType = "PhactoriScatterPlotBlock"

    self.m_YAxisVariableInfo = PhactoriVariableInfo()
    self.m_YAxisVariableInfo.mVariableName = None
    self.m_XAxisVariableInfo = PhactoriVariableInfo()
    self.m_XAxisVariableInfo.mVariableName = None

    self.m_xyzMinMaxTrkC = PlotXYZMinMaxTrkC()
    #self.m_PersistentPoints = None
    self.m_LimitNumberOfPointsFlag = False
    self.m_MaximumNumberOfPoints = 1000
    self.m_GrabMinsAndMaxesWhenLimitedFlag = True
    #self.m_Name = ""
    self.m_PolyData = None
    self.m_Points = None
    self.m_Vertex = None

    #used to keep track paraview mergeblocks filter which is automatically
    #applied to the data before doing the ScatterPlot()
    self.m_MergeBlocks = None

    #self.m_DataCubeAxesInfo = PhactoriDataCubeAxesInfo()

  def SetVarABasedOnVarB(self, ioVarInfoA, ioVarInfoB):
    if PhactoriDbg():
      myDebugPrint3("SetVarABasedOnVarB entered\n")
    #see if x axis variable type needs detecting, and if so, detect it
    detectResult = ioVarInfoB.DetectVariableType(
        self.mInputOperation.GetPvFilter(), True)
    #if detectResult == False:
      #PutPlotOnListToDetectVariableTypeInFutureSteps() ?
    if detectResult == False:
        if PhactoriDbg(10000):
          myDebugPrint3("PhactoriScatterPlotBlock." + \
            "SetVarABasedOnVarB:\n" + \
            "  error!  detect variable type failed.\n",
            10000)

    #set variable name and type for A based on B
    if ioVarInfoB.mVariableType == 'element':
      if PhactoriDbg():
        myDebugPrint3("variable B is element, setting A\n")
      ioVarInfoA.mVariableName = 'GlobalElementId'
      ioVarInfoA.mVariableType = 'element'
    else:
      if PhactoriDbg():
        myDebugPrint3("variable B is node, setting A\n")
      ioVarInfoA.mVariableName = 'GlobalNodeId'
      ioVarInfoA.mVariableType = 'node'
    if PhactoriDbg():
      myDebugPrint3("SetVarABasedOnVarB returning\n")

  def ChooseDefaultVariableIfNecessary(self):
    if PhactoriDbg():
      myDebugPrint3("ChooseDefaultVariableIfNecessary entered\n")
    if PhactoriDbg():
      myDebugPrint3("  x axis var:\n" + self.m_XAxisVariableInfo.SelfToStr())
    if PhactoriDbg():
      myDebugPrint3("  y axis var:\n" + self.m_YAxisVariableInfo.SelfToStr())

    if self.m_XAxisVariableInfo.mVariableType == None:
      if PhactoriDbg():
        myDebugPrint3("  x axis variable needs default from y axis\n")
      if self.m_YAxisVariableInfo.mVariableType == None:
        if PhactoriDbg(
            10000):
          myDebugPrint3("PhactoriScatterPlotBlock." + \
            "ChooseDefaultVariableIfNecessary:\n" + \
            "  neither x axis variable or y axis variable is specified\n",
            10000)
      self.SetVarABasedOnVarB(self.m_XAxisVariableInfo,
          self.m_YAxisVariableInfo)
      SetUpPlotAxisNameDetails(self.m_XAxisVariableInfo,
          self.m_DataCubeAxesInfo.mXAxisInfo)
    elif self.m_YAxisVariableInfo.mVariableType == None:
      if PhactoriDbg():
        myDebugPrint3("  y axis variable needs default from x axis\n")
      self.SetVarABasedOnVarB(self.m_YAxisVariableInfo,
          self.m_XAxisVariableInfo)
      SetUpPlotAxisNameDetails(self.m_YAxisVariableInfo,
          self.m_DataCubeAxesInfo.mYAxisInfo)
    if PhactoriDbg():
      myDebugPrint3("ChooseDefaultVariableIfNecessary returning\n")

  def ParseParametersFromJson(self, inJsn, inPipeAndViewsState):
    if PhactoriDbg():
      myDebugPrint3("PhactoriScatterPlotBlock.ParseParametersFromJson entered\n")
    self.mImageSettings.ParseImageSettingsInfo(inJsn,
        'plot basename', 'plot basedirectory')
    self.mImageFileNameCountSettings.ParseImageFileNameCountSettings(inJsn)

    #parse x axis variable
    self.m_XAxisVariableInfo.ParseVariableNameAndVectorOrTensorComponent(inJsn,
        'x axis variable ')

    #parse y axis variable
    self.m_YAxisVariableInfo.ParseVariableNameAndVectorOrTensorComponent(inJsn,
        'y axis variable ')

    if self.m_YAxisVariableInfo.mVariableType == None and \
      self.m_XAxisVariableInfo.mVariableType:
          if PhactoriDbg(
              10000):
            myDebugPrint3("PhactoriScatterPlotBlock.ParseParametersFromJson " + \
              ":\n  neither x axis variable or y axis variable is specified\n",
              10000)

    #handle name of each axis, particularly in vector component/magnitude case
    SetUpPlotAxisNameDetails(self.m_XAxisVariableInfo,
        self.m_DataCubeAxesInfo.mXAxisInfo)
    SetUpPlotAxisNameDetails(self.m_YAxisVariableInfo,
        self.m_DataCubeAxesInfo.mYAxisInfo)

    #hack to test missing data situations
    #if self.mName == "fooScatter":
    #  myDebugPrint3("hack to test missing data, found fooScatter, using stressthresh\n")
    #  self.mInputOperation = inPipeAndViewsState.mOperationBlocks["stressthresh"]
    #else:
    #  myDebugPrint3("hack to test missing data, did not find fooScatter\n")
    #  self.mInputOperation = inPipeAndViewsState.GetOperationReferredByJson(
    #      'operation', inJsn)

    self.mInputOperation = inPipeAndViewsState.GetOperationReferredByJson(
        'operation', inJsn)

    self.m_xyzMinMaxTrkC.PlotXYZMinMaxTrkCParseJson(inJsn, 'axis ')

    self.mColorSettings.ParsePlotColorSettingsFromJson(inJsn)

    if PhactoriDbg():
      myDebugPrint3(\
          "PhactoriScatterPlotBlock.ParseParametersFromJson returning\n")

global HandleShowAxesWithEmptyDataParaViewIssueFlag
HandleShowAxesWithEmptyDataParaViewIssueFlag = True
#HandleShowAxesWithEmptyDataParaViewIssueFlag = False

#for now, we are trying to have one paraview RenderView which is shared and
#reused for all images everywhere; we may need to change that for plots and
#multiple simulation codes in same input deck
global gSharedRenderView
gSharedRenderView = None

#shared render view to be used by all plots
global gSharedLineChartView
gSharedLineChartView = None


class PhactoriImagesetBlock:
  """contains information corresponding to an imageset block
  """
  def __init__(self):
    self.mName = ""

    #weirdness on first image for an imageset; partial fix is to render twice
    #the first time
    self.mWriteFirstImageTwiceFlag = 0

    #(for now) mOperation and mRepresentation are the 'primary' data sources
    #for this imageset, which includes directing camera positioning if
    #necessary.  Later we may update this to only have a list, not a
    #distinction between primary and other.
    self.mOperation = None
    self.mRepresentation = None

    #(for now) allowing imageset to have additional operations visible.
    #later this may become the only recording of which operations are visible
    #and mOperation and mRepresentation will be eliminated.  Also for
    #now, mVisibleOps[0] should be the same as mOperation and
    #mVisibleReps[0] should be the same as mRepresentation
    self.mVisibleOps = []
    self.mVisibleReps = []
    self.mVisiblePvDataReps = []
    self.mColorLegendRepRefs = []

    self.mVisibleMarkerNames = []
    self.mVisibleMarkers = []

    self.mTextAnnotationNames = []
    self.mTextAnnotations = []

    self.mCamera = None
    self.mImageSettings = PhactoriImageSettings()
    self.mLookDirectionList = []
    self.mLookDirectionFilenameAddon = []

    self.mHandleShowAxesWithEmptyDataParaViewIssueStatus = 0

    #reference for paraview RenderView instance which is shared among all
    #renderings
    self.mSharedPvRenderView2 = None

    #reference for paraview DataRepresentation instance for this imageset.
    #if the imageset has multiple view directions, they can share the same
    #representation
    self.mPvDataRepresentation2 = None

    #self.mOnOffCriteriaName = None
    self.mImagesetOnOffFilter = PhactoriImagesetOnOffFilter()

    self.DeadCellIoFf = None

    #used to control if we want to have image names sequentially count by
    #callback index, or if we want to use a date/time based count
    #(which can use microseconds or not, and can be converted to an integer
    #or not)
    self.mImageFileNameCountSettings = ImageFileNameCountSettings()

  def GetInputPhactoriOperation(self):
    #mainly for compatibility between plot blocks and imageset blocks, give
    #same access to phactori operation coming being visualized
    return self.mOperation

  def HandleShowAxesWithEmptyDataParaViewIssue(self):
    """we are having a paraview related issue wherein if an empty dataset
       (due to threshold or clip) has a show axes setting of true, VTK
       complains.  For now, we are turning off the cube axes if the data is
       empty and back on if data is nonempty (assuming user has cube axes
       set to on"""

    if self.mRepresentation.mDataCubeAxesFlag == False:
      #nothing to be done--it's off so we just return
      return

    #allow this functionality to be eliminated to avoid mpi communication
    global HandleShowAxesWithEmptyDataParaViewIssueFlag
    if HandleShowAxesWithEmptyDataParaViewIssueFlag != True:
      return

    if self.mHandleShowAxesWithEmptyDataParaViewIssueStatus == 2:
      #if PhactoriDbg():
      #  myDebugPrint3("HandleShowAxesWithEmptyDataParaViewIssue at status 2\n")
      return

    paraviewSource = self.mOperation.GetPvFilter()
    UpdatePipelineWithCurrentTimeArgument(paraviewSource)

    pvsDi = paraviewSource.GetDataInformation()

    #we find num cells and num points on this processor and share with others;
    #note that we use UseReduceToSpreadValues which doesn't sum values--we're
    #just trying to ascertain if there are ANY cells or points

    numCells = pvsDi.GetNumberOfCells()
    numPoints = pvsDi.GetNumberOfPoints()
    cellAndPointCount = [numCells, numPoints]
    UseReduceToSpreadValues(cellAndPointCount)

    #see if data is empty or nonempty
    if cellAndPointCount[0] + cellAndPointCount[1] != 0:
      datasetIsEmpty = False
    else:
      datasetIsEmpty = True

    #if PhactoriDbg():
    #  myDebugPrint3("HandleShowAxesWithEmptyDataParaViewIssue executing:\n")
    #  myDebugPrint3("  num cells: " + str(pvsDi.GetNumberOfCells()) + \
    #      "\n  num points: " + str(pvsDi.GetNumberOfPoints()) + "\n")
    #  myDebugPrint3("datasetIsEmpty is " + str(datasetIsEmpty) + "\n")

    #thePvDataRep = self.mParaViewRenderInfoCs[0].mParaViewInfo.DataRepresentation1

    if datasetIsEmpty:
      self.mHandleShowAxesWithEmptyDataParaViewIssueStatus = 1
      ShowCubeAxesXX(self.mSharedPvRenderView2, 'off')
    else:
      self.mHandleShowAxesWithEmptyDataParaViewIssueStatus = 2
      ShowCubeAxesXX(self.mSharedPvRenderView2, 'on')

  def WriteImages(self, datadescription):
    """write out the .png/.jpg/whatever images associated with this imageset
       block for the current timestep/state.  Must loop through camera angles
       and do a write for each one if necessary"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriImagesetBlock::WriteImages entered: " + \
          self.mName + "\n")

    global gPipeAndViewsState
    if self.mImagesetOnOffFilter.TestDrawImagesThisCallback(
            gPipeAndViewsState) == False:
      if PhactoriDbg(100):
        myDebugPrint3("WriteImages returning with no write due to filter\n")
      return

    self.WriteImagesPassedOnOffFilter(datadescription)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriImagesetBlock::WriteImages returning: " + \
          self.mName + "\n")

  def WriteImagesPassedOnOffFilter(self, datadescription):
    """write out the .png/.jpg/whatever images associated with this imageset
       block for the current timestep/state.  Must loop through camera angles
       and do a write for each one if necessary"""

    global gPipeAndViewsState

    if PhactoriDbg(100):
      myDebugPrint3(
        "PhactoriImagesetBlock::WriteImagesPassedOnOffFilter entered\n")

    for ii in range(len(self.mLookDirectionList)):
      oneLookDirection = self.mLookDirectionList[ii]
      oneLookDirectionFilenameAddon = self.mLookDirectionFilenameAddon[ii]

      fname, fnameRR = self.mImageFileNameCountSettings.GetImageFilename(
          datadescription, self.mImageSettings,
          oneLookDirectionFilenameAddon,
          self.mRepresentation.mFilenameAddon,
          self.mCamera.mFilenameAddon)

      #used to do this:
      #view.ViewTime = datadescription.GetTime()
      #maybe need to do this?
      #UpdatePipelineWithCurrentTimeArgument()

      #SetUpViewAndRepresentationBeforeWriteImage(oneViewInfo)
      self.SetUpViewAndRepresentationBeforeWriteImage(oneLookDirection, ii)

      #only need to update color range for first look direction, rest
      #are same
      if ii == 0:
        #UpdateColorRangeImmediatelyBeforeWrite(phactoriImagesetName)
        if self.mRepresentation.mUseFixedColorRange == False:
          UseDataRangeForColorValues(self.mPvDataRepresentation2,
              self.mRepresentation, self.mOperation)
        for ii in range(1, len(self.mVisibleReps)):
          oneVisOp = self.mVisibleOps[ii]
          oneVisRep = self.mVisibleReps[ii]
          oneVisPvDataRep = self.mVisiblePvDataReps[ii]
          UseDataRangeForColorValues(oneVisPvDataRep,
              oneVisRep, oneVisOp)

      if self.mName.startswith("is_element_select") == False:

        for onevisop in self.mVisibleOps:
          if onevisop.mType == "nearestpoints":
            onevisop.mOperationSpecifics.\
                RunCalculationToFindNearestPoints(gPipeAndViewsState)
          if onevisop.mType == "castnormalrays":
            onevisop.mOperationSpecifics.\
                RunCalculationToCastRays(gPipeAndViewsState)
            UpdatePipelineWithCurrentTimeArgument(onevisop.mParaViewFilter)
          if onevisop.mName == "surfaceofinterest1":
            svrng = onevisop.mParaViewFilter.ThresholdRange
            #onevisop.mParaViewFilter.ThresholdRange = [svrng[0]*0.5, svrng[1]*0.5]
            onevisop.mParaViewFilter.ThresholdRange = [1.0, 10.0]
            UpdatePipelineWithCurrentTimeArgument(onevisop.mParaViewFilter)
            onevisop.mParaViewFilter.ThresholdRange = svrng
            UpdatePipelineWithCurrentTimeArgument(onevisop.mParaViewFilter)


        if PhactoriDbg():
          myDebugPrint3("calling WriteImage() " + fname + "\n")
        WriteImage(fname, self.mSharedPvRenderView2,
            Magnification=1)
        if PhactoriDbg():
          myDebugPrint3("returned from WriteImage()\n")
        #handle datetime naming extra work to avoid race condition
        if fnameRR != None:
          if SmartGetLocalProcessId() == 0:
            import os
            os.rename(fname, fnameRR)
        #hack, double write
        #WriteImage(fname, self.mSharedPvRenderView2,
        #    Magnification=1)

      #ClearPvViewAndPvRepAfterWriteImage(oneViewInfo)
      self.ClearPvViewAndPvRepAfterWriteImage()

    if self.mWriteFirstImageTwiceFlag == 0:
      self.mWriteFirstImageTwiceFlag = 1
      if PhactoriDbg(100):
        myDebugPrint3(
          "mWriteFirstImageTwiceFlag triggers (3) re-render of first image\n")
        self.WriteImagesPassedOnOffFilter(datadescription)

    if self.mName.startswith("is_dead_cells"):
      timestep = datadescription.GetTimeStep()

      imageBasename = self.mImageSettings.mImageBasename
      imageBasedirectory = self.mImageSettings.mImageBasedirectory
      import os
      #dcfname = + str(timestep) + "_.csv"
      lpid = SmartGetLocalProcessId()
      #dcfname = imageBasedirectory + os.sep + imageBasename + \
      #        str(timestep) + "." + str(lpid) + "._.csv"
      dcfname = imageBasedirectory + os.sep + imageBasename + \
      str(timestep) + "." + str(lpid) + ".csv"
      if PhactoriDbg(100):
          myDebugPrint3("is_dead_cells: writing:\n" + dcfname + "\n")
      rcrsnParams = self.mOperation.OutputElementListToFile(dcfname)
      if PhactoriDbg(100):
          myDebugPrint3("is_dead_cells: done writing:\n" + dcfname + "\n")

      #write out dead cell count
      global WriteEachDeadCellElementToFiles
      global WriteDeadCellSummaryFile

      if WriteDeadCellSummaryFile:
        myVals = [rcrsnParams.mElementCount,
                rcrsnParams.mKilledByCriteriaCount[2],
                rcrsnParams.mKilledByCriteriaCount[3],
                rcrsnParams.mKilledByCriteriaCount[5]]
        UseReduceToSumArrayOfInts(myVals)

        if lpid == 0:
          if self.DeadCellIoFf == None:
            dcsummaryname = imageBasedirectory + os.sep + imageBasename + \
                    "dead_cell_info.csv"
            self.DeadCellIoFf = open(dcsummaryname, "w")
            self.DeadCellIoFf.write("step, simtime, number of dead cells, " \
                    "killed 2, killed 3, killed 5\n")
          timestep = datadescription.GetTimeStep()
          simTime = gPipeAndViewsState.CurrentDatadescription.GetTime()
          self.DeadCellIoFf.write(
            str(timestep) + ", " + \
            str(simTime) + "," + \
            str(myVals[0]) + ", " + \
            str(myVals[1]) + ", " + \
            str(myVals[2]) + ", " + \
            str(myVals[3]) + "\n")
          self.DeadCellIoFf.flush()

    if self.mName.startswith("is_element_select"):
      timestep = datadescription.GetTimeStep()
      #global gPipeAndViewsState
      simTime = gPipeAndViewsState.CurrentDatadescription.GetTime()
      imageBasename = self.mImageSettings.mImageBasename
      imageBasedirectory = self.mImageSettings.mImageBasedirectory
      import os
      #dcfname = imageBasedirectory + os.sep + imageBasename + ".csv"
      dcfname1 = imageBasedirectory + os.sep + imageBasename + "element.csv"
      dcfname2 = imageBasedirectory + os.sep + imageBasename + "node.csv"
      self.mOperation.OutputSingleElementToTimeHistoryFile(
              dcfname1, dcfname2, timestep, simTime)

    if PhactoriDbg(100):
      myDebugPrint3(
        "PhactoriImagesetBlock::WriteImagesPassedOnOffFilter returning\n")

  def SetUpViewAndRepresentationBeforeWriteImage(self, inLookDirection,
    inLookDirectionIndex):
    """Since it turns out to be much more memory efficient to reuse a single
       render view rather than having one per image, this routine is called
       immediately before WriteImage in order to set the
       camera, background color, image size, etc. up appropriately for the
       WriteImage call"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriImagesetBlock:" \
        "SetUpViewAndRepresentationBeforeWriteImage entered " + \
         str(self.mName) + "\n", 100)
    self.mPvDataRepresentation2.Visibility = 1

    pvRenderView = self.mSharedPvRenderView2
    pvDataRep = self.mPvDataRepresentation2

    #image size
    #pvRenderView.ViewSize = self.mImageSettings.mImageSize
    pvRenderView.ViewSize = [int(self.mImageSettings.mImageSize[0]),
                             int(self.mImageSettings.mImageSize[1])]

    #background color, text color, etc.

    theColorSettings = self.mRepresentation.mColorSettings

    theColorSettings.SetParaviewRvRepColors(pvRenderView, pvDataRep)

    phactoriRep = self.mRepresentation
    #orientation axis on or off for this 3d mesh
    if phactoriRep.mOrientationAxesFlag:
      pvRenderView.OrientationAxesVisibility = 1
    else:
      pvRenderView.OrientationAxesVisibility = 0

    #time annotation invisible (for 3d plot)
    timeAnnStngs = self.mRepresentation.mTimeAnnotationSettings
    if timeAnnStngs.mVisible:
      global gPipeAndViewsState
      if gPipeAndViewsState.mTimeAnnotationPv != None:
        gPipeAndViewsState.mTimeAnnotationPv.\
            mParaViewRepresentation.Visibility = 1

    #cube axes
    if phactoriRep.mDataCubeAxesFlag:
      ShowCubeAxesXX(pvRenderView, 'on', phactoriRep.mDataCubeAxesInfo)

    #edges/surface/wireframe
    #self.mPvDataRepresentation2.Representation = \
    #    phactoriRep.mMeshRenderControl
    self.mPvDataRepresentation2.SetRepresentationType(
            phactoriRep.mMeshRenderControl)

    if phactoriRep.mDoingVolumeRendering:
      if phactoriRep.mScalarOpacityUnitDistance < 0.0:
          phactoriRep.CalculateDefaultScalarOpacityUnitDistance(self.mOperation)
      if phactoriRep.mScalarOpacityUnitDistance != \
              self.mPvDataRepresentation2.ScalarOpacityUnitDistance:
        self.mPvDataRepresentation2.ScalarOpacityUnitDistance = \
                phactoriRep.mScalarOpacityUnitDistance

    #size of points
    self.mPvDataRepresentation2.PointSize = phactoriRep.mPointSize

    #color legend
    for oneColorLegendRepRef in self.mColorLegendRepRefs:
      if oneColorLegendRepRef != None:
        oneColorLegendRepRef.Visibility = 1

    #color by variable or solid color or color by block
    if inLookDirectionIndex == 0:
      #on multicamera 8 setups, we only do this on first (index 0)
      SetForCorrectColorBy(self, self.mOperation,
          self.mPvDataRepresentation2, phactoriRep, False)
    else:
      #for some reason, solid color doesn't stick between camera angles;
      #redo if it is solid color
      if phactoriRep.mColorBySolidColorFlag == True:
        SetForCorrectColorBy(self, self.mOperation,
            self.mPvDataRepresentation2, phactoriRep, False)

    #markers for this imageset made visible
    for oneMarker in self.mVisibleMarkers:
      oneMarker.UpdateAndMakeVisible()

    for oneTextAnnotation in self.mTextAnnotations:
      oneTextAnnotation.UpdateAndMakeVisible()

    #do extra visible operations/representations
    ii = 1
    while ii < len(self.mVisibleReps):
      oneVisOp = self.mVisibleOps[ii]
      oneVisRep = self.mVisibleReps[ii]
      oneVisPvDataRep = self.mVisiblePvDataReps[ii]
      if PhactoriDbg(100):
        myDebugPrint3("add extra visibility " + str(ii) + "  rep:\n" + \
          str(oneVisPvDataRep))
      ii += 1

      oneVisPvDataRep.Visibility = 1
      #theColorSettings = self.mRepresentation.mColorSettings
      #theColorSettings.SetParaviewRvRepColors(pvRenderView, oneVisPvDataRep)

      #edges/surface/wireframe
      oneVisPvDataRep.SetRepresentationType(oneVisRep.mMeshRenderControl)

      if oneVisRep.mDoingVolumeRendering:
        if oneVisRep.mScalarOpacityUnitDistance < 0.0:
            oneVisRep.CalculateDefaultScalarOpacityUnitDistance(oneVisOp)
        if oneVisRep.mScalarOpacityUnitDistance != \
                oneVisPvDataRep.ScalarOpacityUnitDistance:
          oneVisPvDataRep.ScalarOpacityUnitDistance = \
                  oneVisRep.mScalarOpacityUnitDistance


      #color legend
      #if self.mColorLegendRepRef != None:
      #  self.mColorLegendRepRef.Visibility = 1

      #color by variable or solid color or color by block
      if inLookDirectionIndex == 0:
        #on multicamera 8 setups, we only do this on first (index 0)
        SetForCorrectColorBy(self, oneVisOp,
                oneVisPvDataRep, oneVisRep, False)
      else:
        #for some reason, solid color doesn't stick between camera angles;
        #redo if it is solid color
        if phactoriRep.mColorBySolidColorFlag == True:
          #SetForCorrectColorBy(self, self.mPvDataRepresentation2,
          #    phactoriRep, False)
          SetForCorrectColorBy(self, oneVisOp,
                  oneVisPvDataRep, oneVisRep, False)

    if PhactoriDbg(100):
        myDebugPrint3("start setting camera parameters for imageset: " + \
                self.mName + "\n", 100)
    SetParaViewRepresentationCameraParams(
      pvRenderView,
      self.mCamera,
      inLookDirection,
      self.mImageSettings,
      self.mOperation.GetPvFilter())
    if PhactoriDbg(100):
        myDebugPrint3("done setting camera parameters for imageset: " + \
                self.mName + "\n", 100)

    if PhactoriDbg(100):
      myDebugPrint3("SetUpViewAndRepresentationBeforeWriteImage returning\n",
        100)

  def ClearPvViewAndPvRepAfterWriteImage(self):
    """Since it turns out to be much more memory efficient to reuse a single
       render view rather than having one per image, this routine is called
       immediately after WriteImage in order to make stuff invisible again
       before the next item gets a chance to do WriteImage"""
    #if PhactoriDbg(150):
    #  myDebugPrint3("ClearPvViewAndPvRepAfterWriteImage entered\n", 150)

    #cube axes invisible
    ShowCubeAxesXX(self.mSharedPvRenderView2, 'off')

    #3d/pointset dataset invisible (plot or 3d viewing)
    self.mPvDataRepresentation2.Visibility = 0

    #color legend invisible, if 3d viewing
    for oneColorLegendRepRef in self.mColorLegendRepRefs:
      if oneColorLegendRepRef != None:
        if PhactoriDbg(100):
          myDebugPrint3("C inColorLegendRepRef was " + \
            str(oneColorLegendRepRef.Visibility) + \
            " now 0: " + str(oneColorLegendRepRef) + "\n")
        oneColorLegendRepRef.Visibility = 0

    #time annotation invisible (for 3d plot)
    timeAnnStngs = self.mRepresentation.mTimeAnnotationSettings
    if timeAnnStngs.mVisible:
      global gPipeAndViewsState
      if gPipeAndViewsState.mTimeAnnotationPv != None:
        gPipeAndViewsState.mTimeAnnotationPv.\
            mParaViewRepresentation.Visibility = 0

    #markers for this imageset made invisible
    for oneMarker in self.mVisibleMarkers:
        oneMarker.MakeInvisible()

    for oneTextAnnotation in self.mTextAnnotations:
      oneTextAnnotation.MakeInvisible()

    #do extra visible operations/representations
    ii = 1
    while ii < len(self.mVisibleReps):
      oneVisPvDataRep = self.mVisiblePvDataReps[ii]
      if(oneVisPvDataRep != None):
        oneVisPvDataRep.Visibility = 0
      #this is already done above
      #oneColorLegendRepRef = self.mColorLegendRepRefs[ii]
      #if(oneColorLegendRepRef != None):
      #  oneColorLegendRepRef.Visbility = 0
      ii += 1


  def ParseOperationAndRepresentationPair(self, ioPipeAndViewsState, ioJson,
    inOperationKey, inSkipIfOperationKeyNotPresent,
    inRepresentationKey, inRepresentationKeyRequired,
    inAllowPairFromAnotherSource):
    """parse out the operation and associated representation from the
       given json; also optionally throw exception if representation
       is not given.  Used to get multiple visible operation/representation
       pairs"""
    if PhactoriDbg():
      myDebugPrint3("PhactoriImagesetBlock::"
          "ParseOperationAndRepresentationPair entered\n")
      myDebugPrint3(
          "inOperationKey: " + str(inOperationKey) + "\n"
          "inRepresentationKey: " + str(inRepresentationKey) + "\n"
          "incoming json:\n" + str(ioJson) + "\n")


    if inOperationKey not in ioJson:
      if inSkipIfOperationKeyNotPresent:
        return

    if inRepresentationKeyRequired == False and \
      inRepresentationKey not in ioJson:
      #we have to construct and use a default representation, including parsing
      #commands in the imageset for the representation
      if PhactoriDbg():
        myDebugPrint3("  ParseOperationAndRepresentationPair: for imageset " + \
            self.mName + \
            " there is no representation, " +
            "so we must add and reference default\n")
      defaultRepName = self.mName + '_default_representation'
      defaultRepBlockAndWrapper = {defaultRepName: ioJson}
      ParseBlocksC(ioPipeAndViewsState.mRepresentationBlocks,
          defaultRepBlockAndWrapper,
          PhactoriRepresentationBlock,
          ParseOneRepresentationBlockC,
          ioPipeAndViewsState)
      ioJson[inRepresentationKey] = defaultRepName

    pipeForPair = None
    if inAllowPairFromAnotherSource:
      anotherVisOperation, pipeForPair = ioPipeAndViewsState. \
          GetOperationReferredByJsonCrossPipe(inOperationKey, ioJson)
    else:
      anotherVisOperation = ioPipeAndViewsState.GetOperationReferredByJson(
                                inOperationKey, ioJson)

    self.mVisibleOps.append(anotherVisOperation)

    anotherVisRepresentation = None
    if inRepresentationKey in ioJson:
      pipeForRepresentation = ioPipeAndViewsState
      if pipeForPair != None:
        pipeForRepresentation = pipeForPair
      representationName = ioJson[inRepresentationKey]
      if representationName not in pipeForRepresentation.mRepresentationBlocks:
        errStr = "ParseOperationAndRepresentationPair::exception/error\n" + \
            "  imageset (" + str(self.mName) + \
            ") calls for representation (" + \
            str(representationName) + ") which does not exist\n"
        myDebugPrint3AndException(errStr)
      else:
        anotherVisRepresentation = \
          pipeForRepresentation.mRepresentationBlocks[representationName]

    if anotherVisRepresentation == None:
      #at this point we should have a representation--if we don't it
      #was required to be set but not set
      myDebugPrint3AndException(\
        "PhactoriImagesetBlock::ParseOperationAndRepresentationPair\n"
        "Exception/Error: operation missing required corresponding\n"
        "representation.\noperation key: " + inOperationKey + "\n"
        "missing representation key: " + inRepresentationKey + "\n")

    self.mVisibleReps.append(anotherVisRepresentation)
    self.mVisiblePvDataReps.append(None)
    self.mColorLegendRepRefs.append(None)


class PhactoriPipeAndViewsState:
  """Top State Container--Pipeline, Cameras, Representations, Imagesets, Plots

  This is essentially the top-level class for wrapping up everything about the
  phactori-controlled viewing situation.  It keeps a reference to the original
  json (assuming that was how the state was described originally), and stores
  and organizes all the other information.  Presently that means that there
  six sets, each of which stores the class instances corresponding to one
  block type from the json input format.  The block types stored are the
  camera blocks, the representation blocks, the operation blocks, the
  imageset blocks, the scatter plot blocks, and the plot over time blocks.
  A set for each type of block is kept, with the block name (assigned by the
  json author) as a key for that block.  The Operation blocks are basically
  analogous to the ParaView/Catalyst Filters, the imagesets are roughly
  analogous to the ParaView Views (but with image endpoints rather than
  interactive rendering endpoints), the Representation blocks plus the
  camera blocks are analogous to the ParaView Representations.  The
  scatter plot and plot over time blocks are simply specialized descriptions
  of requested plots, which are converted into ParaView filters and views
  to create plots which can be calculated and/or rendered in parallel at
  runtime, as opposed to bringing data to the client.  An imageset block will
  reference a camera block, a representation block, and an operation block and
  contains some additional parameters (e.g. image size and file basename)
  which will describe a view to be rendered--repeatedly at different times
  when using insitu.  The camera blocks describe 3-D viewpoints, sometimes in
  absolute 3D terms, sometimes dependent on the data.  The representation
  blocks control how the view looks, e.g. if element surface and edges are
  rendered or just surfaces and if we show axes and color legends.  The
  operation blocks can describe a potentially complex data pipeline which can
  allow rendering at various points along the pipeline.  The pipeline is
  currently assumed to have operations that only have one input and one
  output.  Multiple operations can have the same input, so a 'tree' is
  possible, not just a single linear pipe.  Operations with multiple inputs
  and outputs are conceivable but not currently implemented.  One obvious
  possible upgrade would be to allow 'geometryoutput' blocks which would
  output geometry rather than images, presumably geometry which is at the
  back end of a pipeline (e.g. to extract isosurface and the decimate).
  """
  def __init__(self):
    self.mJsonDescription = {}
    self.mCameraBlocks = {}
    self.mRepresentationBlocks = {}
    #default operation for incoming input
    self.mIncomingDefaultOperation = PhactoriOperationBlock()
    self.mIncomingDefaultOperation.mName = "default incoming input"
    self.mOperationBlocks = {}
    self.mImagesetBlocks = {}
    self.mScatterPlotBlocks = {}
    self.mPlotOverTimeBlocks = {}
    self.mImagesetOnOffCriteriaBlocks = {}
    self.mMarkerBlocks = {}
    self.mTextAnnotationBlocks = {}
    self.mCallbackDateTime = datetime.datetime.now()

    #used to essentially keep track of the number of times catalyst is
    #called back, for keying when we do or do not need to recalculate
    #things such as variable mins/maxes
    self.mFrameTagCounter = 0

    self.mImagesetOnOffFilter = PhactoriImagesetOnOffFilter()

    #the following items help id this instance to one results output block
    #in the input deck, and to a given restart or remesh instance
    self.mOutputResultsBlockId = None
    self.mCurrentDatabaseDummyFileName = None
    self.mProducer = None
    self.mPipeAndViewsState = None
    self.mRemeshRestartTag = None
    self.mSeparatorString = '_'
    self.mRemeshCount = 0
    self.mOutputResultsBlockCountId = -1
    self.mImageSetCounter = 0

    self.CurrentDatadescription = None
    self.mDefaultBasedirectory = ''

    self.mBlockIdForRestart = ''

    #self.mInteractionEnabled = True
    self.mInteractionEnabled = False
    self.mInteractionTriggerCounter = 0
    self.mInteractionTriggerTries = 7
    self.mInteractionTriggerSleep = 1
    self.mInteractionRepeatPauseSim = False

    #class instance to filter rendering of images on/off based on simulation
    #data
    self.mImagesetOnOffFilter = PhactoriImagesetOnOffFilter()

    #this item manges the paraview representation and color stuff for the
    #time annotation source.
    #since we are now using one RenderView, we can share a single time
    #annotation source representation.  However, it is conceivable that
    #we may eventually want to show different times for some reason, and
    #in that case we will need to have each imageset have its own time
    #annotation paraview stuff
    self.mTimeAnnotationPv = None

  def ExportOperationsData(self, datadescription):
    """go through all the operation blocks and call the ExportOperationData()
       method on each to allow each operation to export non-image data if
       they need to do so"""
    for operationName, operationInstance in self.mOperationBlocks.iteritems():
      operationInstance.ExportOperationData(datadescription)

  def WriteImages(self, datadescription):
    """go through all imageset blocks and plot blocks and have each of them
       write out their images"""

    #since we are sharing a paraview RenderView instance, we need to go
    #through all the plots and imagesets and turn them all to invisible
    #so they can be made visible one at a time during the write image loop

    #NOTE:  when we change this it needs to work across multiple pipes
    #i.e. all pipes need to make selves invisible
    self.SetUpForWriteImageLoop()

    for imagesetName, imagesetInstance in self.mImagesetBlocks.iteritems():
      imagesetInstance.WriteImages(datadescription)
    for splotName, splotInstance in self.mScatterPlotBlocks.iteritems():
      splotInstance.WriteImages(datadescription)
    for tplotName, tplotInstance in self.mPlotOverTimeBlocks.iteritems():
      tplotInstance.WriteImages(datadescription)


  def SetUpForWriteImageLoop(self):
    """Since it turns out to be much more memory efficient to reuse a single
    render view rather than having one per image, this routine is called
    from WriteImages immediately before we start looping through the
    imagesets and plots to render images in order to basically set all
    paraview Representations to Visibility=0 so that we can turn them
    visible appropriately as we do a WriteImage for each one."""

    for imagesetName, imagesetInstance in self.mImagesetBlocks.iteritems():
      imagesetInstance.ClearPvViewAndPvRepAfterWriteImage()
    for splotName, splotInstance in self.mScatterPlotBlocks.iteritems():
      splotInstance.ClearPvViewAndPvRepAfterWriteImage()
    for tplotName, tplotInstance in self.mPlotOverTimeBlocks.iteritems():
      tplotInstance.ClearPvViewAndPvRepAfterWriteImage()

  def GetOperationReferredByJson(self, inOperationNameKey, inJsonDict):
    """utility function to get a reference to an operation instance given a
       dict (inJsonDict) which may (or may not) contain the key
       (inOperationNameKey) which tags the name of the operation to be
       found.  If the key isn't there, we return mIncomingDefaultOperation,
       if the key grabs a name which isn't in mOperationBlocks, that is an
       exception"""
    if inOperationNameKey in inJsonDict:
      nameOfOperationToUse = inJsonDict[inOperationNameKey]
      if nameOfOperationToUse not in self.mOperationBlocks:
        errStr = '  in GetOperationReferredByJson operation with name ' +\
            str(nameOfOperationToUse) + ' does not exist\n'
        if PhactoriDbg():
          myDebugPrint3(errStr)
        raise Exception(errStr)
      returnOperation = self.mOperationBlocks[nameOfOperationToUse]
    else:
      returnOperation = self.mIncomingDefaultOperation

    return returnOperation

  def GetOperationReferredByJsonCrossPipe(self, inOperationNameKey, inJsonDict):
    """same as GetOperationReferredByJson(); however, if the operation is not
       in our self.mOperationBlocks, then we will also search all the other
       pipe and views state instances for the operation.  If we find it, we
       will return it along with the pipe we found it in"""
    if inOperationNameKey in inJsonDict:
      nameOfOperationToUse = inJsonDict[inOperationNameKey]
      if nameOfOperationToUse in self.mOperationBlocks:
        #this operation is in our own operation blocks; return it and None
        #to follow the normal course of (non cross-pipe) operation
        returnOperation = self.mOperationBlocks[nameOfOperationToUse]
        returnPipe = None
      else:
        returnOperation = None
        returnPipe = None
        for onePipeKey, onePipe in gPhactoriPipeRootMap.iteritems():
          if nameOfOperationToUse in onePipe.mOperationBlocks:
            #we found the operation in another pipe; return that and the other
            #pipe
            returnOperation = onePipe.mOperationBlocks[nameOfOperationToUse]
            returnPipe = onePipe
            break
        if returnOperation == None:
          #operation not found anywhere:  this is fatal
          errStr = "  in GetOperationReferredByJsonCrossPipe " +\
              "operation with name " +\
              str(nameOfOperationToUse) + " does not exist in any pipe\n"
          if PhactoriDbg():
            myDebugPrint3(errStr)
          raise Exception(errStr)
    else:
      #this code path requires an operation name, no implied default allowed
      errStr = "  in GetOperationReferredByJsonCrossPipe " +\
          "requires specific operation name; no implied default allowed\n"
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    return returnOperation, returnPipe


def PerRendersetInitialization(datadescription):
  """do setup necessary on a once-per insitu callback basis

  There are some things that should be done once each time everything is
  rendered, particularly incrementing a callback counter so that we can keep
  track and do some things, such as getting data bounds, as few times as
  possible (i.e. don't get data bounds for an operation stage multiple times
  when once is sufficient)
  """
  global gPipeAndViewsState
  gPipeAndViewsState.mFrameTagCounter = gPipeAndViewsState.mFrameTagCounter + 1
  if PhactoriDbg(100):
    myDebugPrint3("-X-X-X-X- mFrameTagCounter now " + \
            str(gPipeAndViewsState.mFrameTagCounter) + " -X-X-X-X-\n", 100)

  if PhactoriDbg():
    paraviewSource = gPipeAndViewsState.mIncomingDefaultOperation.GetPvFilter()
    pvsDi = paraviewSource.GetDataInformation()
    #pointDataInfo = paraviewSource.GetPointDataInformation()
    #cellDataInfo = paraviewSource.GetCellDataInformation()
    myDebugPrint3("  PerRendersetInitialization:\n  num cells: " + str(pvsDi.GetNumberOfCells()) + "\n  num points: " + str(pvsDi.GetNumberOfPoints()) + "\n")

  #we assume to begin with that all global data bounds in operations are
  #invalid as data may have changed
  for oneOperationName, oneOperation in \
          gPipeAndViewsState.mOperationBlocks.iteritems():
    oneOperation.mDataBoundsIsCurrent = False

  #update the time tag in all views
  global gSharedRenderView
  gSharedRenderView.ViewTime = datadescription.GetTime()

def UpdateRepresentationColorBy(ioPhactoriImagesetBlock):
  """given a para view representation and source and a phactori representation,
     set the color by variable in the paraview representation based on the
     phactori representation.  Used to interactively change what the
     'color by' variable is between insitu callbacks"""
  if PhactoriDbg():
    myDebugPrint3("UpdateRepresentationColorBy entered\n")

  for ii in range(0, len(ioPhactoriImagesetBlock.mVisibleReps)):
    oneVisOp = ioPhactoriImagesetBlock.mVisibleOps[ii]
    oneVisRep = ioPhactoriImagesetBlock.mVisibleReps[ii]
    oneVisPvDataRep = ioPhactoriImagesetBlock.mVisiblePvDataReps[ii]
    oneColorLegendRepRef = ioPhactoriImagesetBlock.mColorLegendRepRefs[ii]
    if oneVisRep != None:
      retVal = UpdateRepresentationColorBySub1(
          ioPhactoriImagesetBlock.mSharedPvRenderView2, oneVisPvDataRep,
          oneVisRep, oneVisOp, oneColorLegendRepRef)
      if oneColorLegendRepRef == None and retVal != None:
        if PhactoriDbg(100):
          myDebugPrint3("in PhactoriImagsetBlock named: " + \
            ioPhactoriImagesetBlock.mName + "\n" + \
            "mColorLegendRepRefs " + str(ii) +
            "is set to " + str(retVal) + "\n")
        ioPhactoriImagesetBlock.mColorLegendRepRefs[ii] = retVal

  if PhactoriDbg():
    myDebugPrint3("UpdateRepresentationColorBy returning\n")

def UpdateRepresentationColorBySub1(inPvView, inPvRep,
        inPhactoriRep, inPhactoriOp, inColorLegendRepRef):
  if PhactoriDbg():
    myDebugPrint3("UpdateRepresentationColorBySub1 entered\n")

  colorVarInfo = inPhactoriRep.mColorVariableInfo

  if PhactoriDbg():
    myDebugPrint3("setting color by to:\n" + colorVarInfo.SelfToStr() + "\n")

  colorVarName = colorVarInfo.mVariableName
  detectResult = colorVarInfo.DetectVariableType(
          inPhactoriOp.GetPvFilter(), True)
  if detectResult == False:
    if PhactoriDbg(100):
      myDebugPrint3('no detection: returning\n', 100)
    return None

  if colorVarInfo.mVariableType == 'element':
    thePvVarType = gCellsString
  elif colorVarInfo.mVariableType == 'node':
    thePvVarType = gPointsString
  else:
    if PhactoriDbg(100):
      myDebugPrint3('bad element type: returning\n', 100)
    return None

  if colorVarName != '':
    if inPhactoriRep.mColorVariableInfo.mVariableIsVectorComponent:
      ColorByVariableComponentOrMagnitudeXX(inPvRep, inPhactoriRep,
          colorVarName, 'Component', colorVarInfo.mVariableComponent,
          thePvVarType, inPhactoriRep.mColorMapSettings)
    elif inPhactoriRep.mColorVariableInfo.mVariableIsVectorMagnitude:
      ColorByVariableComponentOrMagnitudeXX(inPvRep, inPhactoriRep,
          colorVarName, 'Magnitude', colorVarInfo.mVariableComponent,
          thePvVarType, inPhactoriRep.mColorMapSettings)
    else:
      ColorByVariableScalarXX(inPvRep, inPhactoriRep, colorVarName,
          thePvVarType, inPhactoriRep.mColorMapSettings)

  #update the color legend widget
  if inPhactoriRep.mColorLegendFlag:
    onoffFlagString = 'on'
  else:
    onoffFlagString = 'off'
  retVal = ShowDataColorLegendXX(inPvView, onoffFlagString,
      inPhactoriRep.mColorLegendPositionAndSize, inPhactoriRep.mColorSettings,
      inColorLegendRepRef, inPvRep)

  if PhactoriDbg():
    myDebugPrint3("UpdateRepresentationColorBySub1 returning\n")

  return retVal

def DuringRestartUseJsonToSetUp(jsonIn, ioPipeAndViewsState):
  """used by process zero (broadcast send process) as well as other processes
     (broadcast receive processes) to actually take the info in json format
     and set the system up for proper behavior after restart, particularly
     data ranges and plots over time"""
  #go through representations and have each add it's state info to jsonOut
  if "RepresentationsRestartInfo" not in jsonIn:
    if PhactoriDbg():
      myDebugPrint3("StartOfVisualizationCallback returning, no rep info\n")
    return

  representationsJsonIn = jsonIn["RepresentationsRestartInfo"]
  for oneRepName, oneRep in \
      ioPipeAndViewsState.mRepresentationBlocks.iteritems():
    if oneRepName in representationsJsonIn:
      if PhactoriDbg():
        myDebugPrint3("restart setup callback json for rep: " + oneRepName + "\n")
      oneRep.SetFromRestartInfo(representationsJsonIn[oneRepName])
    else:
      if PhactoriDbg():
        myDebugPrint3("warning: Representation named " + oneRepName +
          " not in restart info")

  plotsOverTimeJsonIn = jsonIn["PlotsOverTimeRestartInfo"]
  for onePlotOtName, onePlotOt in \
      ioPipeAndViewsState.mPlotOverTimeBlocks.iteritems():
    if onePlotOtName in plotsOverTimeJsonIn:
      if PhactoriDbg():
        myDebugPrint3("plot over time setup callback json for rep: " + \
            onePlotOtName + "\n")
      onePlotOt.SetFromRestartInfo(plotsOverTimeJsonIn[onePlotOtName])
    else:
      if PhactoriDbg():
        myDebugPrint3("warning: plot over time named " + onePlotOtName +
          " not in restart info")


def HandleRestartUpdateProcessZero(ioPipeAndViewsState):
  if PhactoriDbg():
    myDebugPrint3("HandleRestartUpdateProcessZero entered\n")
  #import pdb
  #pdb.set_trace()

  #broadcast for full multiprocess functionality
  #color range needs no broadcast, as that will be managed when color range
  #calculated for the next frame (process 1 will share the maxes and mins
  #it has)
  import vtkParallelCorePython
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()

  #try to read file written (presumably last restart place)
  #and convert to json with expected string type
  validDataIsAvailableToBroadcast = True
  try:
    #read map/json from file
    import json
    import os

    restartInfoFilename = ioPipeAndViewsState.mDefaultBasedirectory + \
        os.sep + ioPipeAndViewsState.mBlockIdForRestart + ".phrs"

    inFile = open(restartInfoFilename, "rb")
    jsonIn = json.load(inFile)
    inFile.close()
    jsonIn = convertJsonUnicodeToStrings(jsonIn)

  except:
    if PhactoriDbg():
      myDebugPrint3("no info file to read or json failed\n")
    #need to broadcast to non-zero processes that there is no info--we
    #do this by indicating a length of zero
    validDataIsAvailableToBroadcast = False

  if validDataIsAvailableToBroadcast:
    try:
      DuringRestartUseJsonToSetUp(jsonIn, ioPipeAndViewsState)
    except:
      if PhactoriDbg():
        myDebugPrint3("exception trying to use info for setting\n")
      validDataIsAvailableToBroadcast = False

  #if no valid data is available (either no file read or bad data), send
  #indicator to other processes by sending 0 data size
  if validDataIsAvailableToBroadcast == False:
    if PhactoriDbg():
      myDebugPrint3("no valid restart data to send, send signal of no data\n")
    vtkInfoBufferSizeArray = vtk.vtkIntArray()
    vtkInfoBufferSizeArray.SetNumberOfTuples(1)
    vtkInfoBufferSizeArray.SetValue(0, 0)
    globalController.Broadcast(vtkInfoBufferSizeArray, 0)
    if PhactoriDbg():
      myDebugPrint3("HandleRestartUpdateProcessZero returning 2\n")
    return

  try:

    outStr = json.dumps(jsonIn)

    numberOfValues = len(outStr)

    if PhactoriDbg():
      myDebugPrint3("valid restart data to send.  Sending size: " + \
        str(numberOfValues) + "\n")
    #send size of data buffer first
    vtkInfoBufferSizeArray = vtk.vtkIntArray()
    vtkInfoBufferSizeArray.SetNumberOfTuples(1)
    vtkInfoBufferSizeArray.SetValue(0, numberOfValues)
    globalController.Broadcast(vtkInfoBufferSizeArray, 0)

    if PhactoriDbg():
      myDebugPrint3("now sending data\n")
    #now fill out data and send it
    vtkOutArray = vtk.vtkCharArray()
    vtkOutArray.SetNumberOfTuples(numberOfValues)
    for jj in range(0, numberOfValues):
      vtkOutArray.SetValue(jj, outStr[jj])

    globalController.Broadcast(vtkOutArray, 0)

    if PhactoriDbg():
      myDebugPrint3("done sending data\n")
  except:
    #broadcast error; not really recoverable;
    if PhactoriDbg():
      myDebugPrint3("exception trying to broadcast\n")
    if PhactoriDbg():
      myDebugPrint3("HandleRestartUpdateProcessZero returning 4\n")
    return

  if PhactoriDbg():
    myDebugPrint3("HandleRestartUpdateProcessZero returning\n")


def HandleRestartUpdateProcessNotZero(ioPipeAndViewsState):

  if PhactoriDbg():
    myDebugPrint3("HandleRestartUpdateProcessNotZero entered\n")

  try:
    #do broadcast (receive from process 0) to receive update data
    import vtkParallelCorePython
    pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
    globalController = pm.GetGlobalController()

    if PhactoriDbg():
      myDebugPrint3("trying to receive data buffer length")

    #first, get size of data.  If it is zero, we are done and there won't be
    #another broadcast
    vtkInfoBufferSizeArray = vtk.vtkIntArray()
    vtkInfoBufferSizeArray.SetNumberOfTuples(1)
    #vtkInfoBufferSizeArray.SetValue(0, 0)
    globalController.Broadcast(vtkInfoBufferSizeArray, 0)
    bufferLen = vtkInfoBufferSizeArray.GetTuple1(0)
    if bufferLen == 0:
      if PhactoriDbg():
        myDebugPrint3("broadcast indicated 0 data, so we are done\n")
      if PhactoriDbg():
        myDebugPrint3("HandleRestartUpdateProcessNotZero returning 2\n")
      return;

    if PhactoriDbg():
      myDebugPrint3("buffer length: " + str(bufferLen) + \
        ", now trying to receive data\n")
    #first go round, fixed buffer length, next expandable if necessary
    vtkOutArray = vtk.vtkCharArray()
    vtkOutArray.SetNumberOfTuples(bufferLen)
    #for jj in range(0, bufferLen):
    #  vtkOutArray.SetValue(jj, 0)

    #receive broadcast of info from process 0
    globalController.Broadcast(vtkOutArray, 0)
    if PhactoriDbg():
      myDebugPrint3("data received successfully\n")
  except:
    #broadcast error; not really recoverable;
    if PhactoriDbg():
      myDebugPrint3("exception trying to receive broadcast\n")
    if PhactoriDbg():
      myDebugPrint3("HandleRestartUpdateProcessNotZero returning 2\n")
    return

  try:
    if PhactoriDbg():
      myDebugPrint3("converting data to json format string\n")
    #convert to json-format string
    jsonCharList = []
    for ii in range(0, bufferLen):
      oneChar = vtkOutArray.GetTuple1(ii)
      jsonCharList.append(oneChar)

    jsonInString = "".join(jsonCharList)
    import json

    if PhactoriDbg():
      myDebugPrint3("json received from broadcast (converted to string):\n" +
        jsonInString)

    #convert json format string to python dict
    jsonIn = json.loads(jsonInString)

    if PhactoriDbg():
      myDebugPrint3("json received from broadcast (converted to dict):\n" +
        str(jsonIn))

    #use info to set up for post-restart (color ranges, time plots, etc.)
    DuringRestartUseJsonToSetUp(jsonIn, ioPipeAndViewsState)

  except:
    #broadcast error; not really recoverable;
    if PhactoriDbg():
      myDebugPrint3("exception trying to use broadcast info as json set info\n")
    if PhactoriDbg():
      myDebugPrint3("HandleRestartUpdateProcessNotZero returning 3\n")
    return

  if PhactoriDbg():
    myDebugPrint3("HandleRestartUpdateProcessNotZero returning\n")


def HandleRestartUpdateForOutputResultsBlock(ioPipeAndViewsState):
  """this method is a companion to SaveRestartInfoCallback.  It is called
     after the pipeline is first initialized for the given output block, and
     it will read the data stored by SaveRestartInfoCallback in a file.
     The presumption is that the pipeline will be initialized for the
     given output block only once when the system starts or only once
     when the system restarts.
     The intent is to deal with restart by allowing persistent
     information across restarts.  See SaveRestartInfoCallback."""
  if PhactoriDbg():
    myDebugPrint3("HandleRestartUpdateForOutputResultsBlock entered\n")

  #we need to behave differently for process 0 and other processes; process 0
  #needs to read the file and, if necessary, broadcast info to other
  #processes, while other processes need to receive info broadcast from 0, if
  #necessary
  if SmartGetLocalProcessId() == 0:
    HandleRestartUpdateProcessZero(ioPipeAndViewsState)
  else:
    HandleRestartUpdateProcessNotZero(ioPipeAndViewsState)


  if PhactoriDbg():
    myDebugPrint3("HandleRestartUpdateForOutputResultsBlock returning\n")


def SaveRestartInfoCallback():
  """this method is called back at the tail end of the visualization
     from the simulation.  The purpose of this function is to collect and
     save out information which is necessary for a restart to have the
     restart begin in the same state as if the simulation had run to this
     point.  Items include things such as
     recording color range mins and maxes and data values in plots over time
     so that they can be correctly restarted in the case of a simulation
     restart"""

  if PhactoriDbg():
    myDebugPrint3("SaveRestartInfoCallback entered\n")
  #only process zero gets to write anything out:  info should have already
  #been shared appropriately to process zero and we don't want every process
  #hitting the disk
  if SmartGetLocalProcessId() != 0:
    if PhactoriDbg():
      myDebugPrint3("SaveRestartInfoCallback returning, not process zero\n")
    return

  #we construct a json-format/python map style data item which we will save
  #out to a file.  We start with a blank map and allow each image set and
  #plot over time to add stuff to it.  Then we write it out to a file.

  jsonOut = {}
  global gPipeAndViewsState

  #go through representations and have each add it's state info to jsonOut
  representationsJsonOut = {}
  for oneRepName, oneRep in \
      gPipeAndViewsState.mRepresentationBlocks.iteritems():
    if PhactoriDbg():
      myDebugPrint3("making end vis callback json for rep: " + oneRepName + "\n")
    representationsJsonOut[oneRepName] = oneRep.GetRestartInfo()

  jsonOut['RepresentationsRestartInfo'] = \
      representationsJsonOut

  #go through plots over time and have each add it's state info to jsonOut
  plotsOverTimeJsonOut = {}
  for onePlotOtName, onePlotOt in \
      gPipeAndViewsState.mPlotOverTimeBlocks.iteritems():
    if PhactoriDbg():
      myDebugPrint3("making end vis callback json for plot over time: " +
        onePlotOtName + "\n")
    plotsOverTimeJsonOut[onePlotOtName] = onePlotOt.GetRestartInfo()

  jsonOut['PlotsOverTimeRestartInfo'] = plotsOverTimeJsonOut

  #write info to file
  import json

  import os
  restartInfoFilename = gPipeAndViewsState.mDefaultBasedirectory + \
      os.sep + gPipeAndViewsState.mBlockIdForRestart + ".phrs"

  outStr = json.dumps(jsonOut, sort_keys=True,
        indent = 2, separators=(',', ': '))
  outFile = open(restartInfoFilename, 'w+b')
  outFile.write(outStr)
  outFile.close()

  if PhactoriDbg():
    myDebugPrint3("SaveRestartInfoCallback returning\n")

def UpdateOneImagesetViewsWhichMayChangeWithData(ioImageset):
  if PhactoriDbg(100):
    myDebugPrint3("UpdateOneImagesetViewsWhichMayChangeWithData entered\n", 100)
  theCamera = ioImageset.mCamera

  #this is a bit hacky--we should probably update the para view sources when
  #we update operations, which we aren't doing separately yet
  #theParaViewSource = ioImageset.mOperation.GetPvFilter()
  theParaViewSource = ioImageset.mOperation.GetOutgoingPvGeometryFilter()
  UpdatePipelineWithCurrentTimeArgument(theParaViewSource)

  UpdateRepresentationColorBy(ioImageset)

  if theCamera.MayChangeWithData() == False:
    if PhactoriDbg():
      myDebugPrint3("  camera: " + theCamera.mName + "   imageset: " + ioImageset.mName + " is set to not change with data \n")
    return

  #it no longer makes sense to loop through and set stuff up here, as
  #we will do it right before we call WriteImage
  #for xx in ioImageset.mParaViewRenderInfoCs:
  #  if PhactoriDbg():
  #    myDebugPrint3("  renderview.CameraPosition before: " + str(xx.mParaViewInfo.RenderView1.CameraPosition) + "\n")
  #  lookDirectionToUse = xx.mLookDirection
  #  if theCamera.mType == 'camera':
  #    lookDirectionToUse = theCamera.mLookDirection
  #  if PhactoriDbg():
  #    myDebugPrint3("  camera: " + theCamera.mName + "   imageset: " + ioImageset.mName + "  look dir: " + str(lookDirectionToUse) + "\n")
  #  SetParaViewRepresentationCameraParams(xx.mParaViewInfo.RenderView1, theCamera, lookDirectionToUse, ioImageset.mImageSettings, theParaViewSource)
  #  if PhactoriDbg():
  #    myDebugPrint3("  renderview.CameraPosition after : " + str(xx.mParaViewInfo.RenderView1.CameraPosition) + "\n")

  if PhactoriDbg(100):
    myDebugPrint3("UpdateOneImagesetViewsWhichMayChangeWithData returning\n", 100)
  #SetActiveSource(savedActiveSource)

def UpdateAllOperationsWhichMayChangeWithData():
  if PhactoriDbg(100):
    myDebugPrint3("UpdateAllOperationsWhichMayChangeWithData entered\n", 100)
  global gPipeAndViewsState
  for oneOperationName, oneOperation in gPipeAndViewsState.mOperationBlocks.iteritems():
    oneOperation.DoUpdateDueToChangeInData(gPipeAndViewsState)
  if PhactoriDbg(100):
    myDebugPrint3("UpdateAllOperationsWhichMayChangeWithData returning\n", 100)

def UpdateAllImagesetViewsWhichMayChangeWithData():
  global gPipeAndViewsState
  for oneImagesetName, oneImageset in gPipeAndViewsState.mImagesetBlocks.iteritems():
    oneImageset.HandleShowAxesWithEmptyDataParaViewIssue()
    UpdateOneImagesetViewsWhichMayChangeWithData(oneImageset)

def HandleOperationShortcuts2(inBlockName, inJson, ioOperationBlocks, inCount):
  if PhactoriDbg(100):
    myDebugPrint3('HandleOperationShortcuts entered2\n', 100)

  lastShortcutOperationName = None

  opSuffix = str(inCount)

  if "threshold" in inJson:
    stuff = inJson["threshold"]
    if len(stuff) != 4:
      errStr = 'HandleOperationShortcuts2 threshold shortcut bad stuff, not 4 items\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)
    #stuff should have 4 items
    #[scalar|vector magnitude|vector component|tensor
        #component, component<variable name>, namekeep between|keep above|keep
        #below, [<val1>, <val2>]]
    newOperationJson = {}
    newOperationJson["type"] = "threshold"
    varkey = "variable " + stuff[0]
    newOperationJson[varkey] = stuff[1]
    thresholdValues = stuff[3]
    if stuff[2] == "keep between":
      newOperationJson["keep between"] = thresholdValues
    elif stuff[2] == "keep below":
      newOperationJson["keep below"] = thresholdValues[0]
    elif stuff[2] == "keep above":
      newOperationJson["keep above"] = thresholdValues[0]
    else:
      errStr = 'HandleOperationShortcuts2 bad stuff[2], not keep above|below|between\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    if PhactoriDbg():
      myDebugPrint3('got threshold operation shortcut ' + str(newOperationJson) + '\n')
    if lastShortcutOperationName != None:
      newOperationJson["input"] = lastShortcutOperationName
    lastShortcutOperationName = "thresholdshortcutoperation" + opSuffix
    ioOperationBlocks[lastShortcutOperationName] = newOperationJson


  if "boxclip" in inJson:
    stuff = inJson["boxclip"]
    if len(stuff) != 3:
      errStr = 'HandleOperationShortcuts2 boxclip shortcut bad stuff, not 3 items\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)
    center = stuff[0]
    extents = stuff[1]
    insideOrOutside = stuff[2]

    newOperationJson = {}
    newOperationJson["type"] = "boxclip"
    newOperationJson["center at absolute point"] = center
    newOperationJson["absolute extents"] = extents
    if insideOrOutside == "keep inside box":
      newOperationJson["keep inside box"] = True
    elif insideOrOutside == "keep outside box":
      newOperationJson["keep inside box"] = False
    else:
      newOperationJson["keep inside box"] = True

    if PhactoriDbg():
      myDebugPrint3('got boxclip operation shortcut ' + str(newOperationJson) + '\n')
    if lastShortcutOperationName != None:
      newOperationJson["input"] = lastShortcutOperationName
    lastShortcutOperationName = "boxclipshortcutoperation" + opSuffix
    ioOperationBlocks[lastShortcutOperationName] = newOperationJson

  if "clip" in inJson:
    stuff = inJson["clip"]
    if len(stuff) != 2:
      errStr = 'HandleOperationShortcuts2 clip shortcut bad stuff, not 2 items\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)
    pointOnPlane = stuff[0]
    planeNormal = stuff[1]
    if PhactoriDbg():
      myDebugPrint3('got clip operation shortcut ' + str(pointOnPlane) + str(planeNormal) + '\n')
    newOperationJson =  {"type":"clip", "absolute point on plane":pointOnPlane, "plane normal":planeNormal}
    if lastShortcutOperationName != None:
      newOperationJson["input"] = lastShortcutOperationName
    lastShortcutOperationName = "clipshortcutoperation" + opSuffix
    ioOperationBlocks[lastShortcutOperationName] = newOperationJson

  if "contour" in inJson:
    stuff = inJson["contour"]
    if len(stuff) != 4:
      errStr = 'HandleOperationShortcuts2 contour shortcut bad stuff, not 4 items\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)
    variableInfo = stuff[0]
    variableName = stuff[1]
    valueListOrSequence = stuff[2]
    contourValues = stuff[3]
    newOperationJson = {}
    newOperationJson["type"] = "contour"
    varkey = "variable " + variableInfo
    newOperationJson[varkey] = variableName
    if (valueListOrSequence.lower()=="value list"):
      newOperationJson["contour value"] = contourValues
    elif (valueListOrSequence.lower()=="value sequence"):
      newOperationJson["contour value sequence"]=contourValues
    else:
      errStr = """error! contour shortcut did not contain 'value list' or 'value sequence' in
                      HandleOperationShortcuts2\n"""
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)
    if lastShortcutOperationName != None:
      newOperationJson["input"] = lastShortcutOperationName
    lastShortcutOperationName = "contourshortcutoperation" + opSuffix
    ioOperationBlocks[lastShortcutOperationName] = newOperationJson

  if "slice" in inJson:
    stuff = inJson["slice"]
    if len(stuff) != 2:
      errStr = 'HandleOperationShortcuts2 slice shortcut bad stuff, not 2 items\n'
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)
    pointOnPlane = stuff[0]
    planeNormal = stuff[1]
    if PhactoriDbg():
      myDebugPrint3('got slice operation shortcut ' + str(pointOnPlane) + str(planeNormal) + '\n')
    newOperationJson =  {"type":"slice", "absolute point on plane":pointOnPlane, "plane normal":planeNormal}
    if lastShortcutOperationName != None:
      newOperationJson["input"] = lastShortcutOperationName
    lastShortcutOperationName = "sliceshortcutoperation" + opSuffix
    ioOperationBlocks[lastShortcutOperationName] = newOperationJson

  if lastShortcutOperationName != None:
    inJson["operation"] = lastShortcutOperationName
  if PhactoriDbg(100):
    myDebugPrint3('HandleOperationShortcuts returning\n', 100)

#go through imageset json, and locate operation shortcuts.  If found, create
#json definition of operations and add them to the structure for later
#creation, and add reference to the new operations in the imageblock json
def HandleOperationShortcuts(ioImagesetBlocks, ioOperationBlocks):
  if PhactoriDbg(100):
    myDebugPrint3('HandleOperationShortcuts entered\n', 100)
  count = 0
  for blockName, blockParamsJson in ioImagesetBlocks.iteritems():
    count = count + 1
    HandleOperationShortcuts2(blockName, blockParamsJson, ioOperationBlocks, count)
  if PhactoriDbg(100):
    myDebugPrint3('HandleOperationShortcuts returning\n', 100)
  return count

def ConstructAndAddPlotsOverTimeJson(outPlotOverTimeJsonBlocks, inJsonToCopy,
        inYAxisVariableName, inNumComponents):
  """helper function for TestForAndConstructAllPlotsOverTime, creates a
     PlotOverTime json block, not for external use"""

  global gPipeAndViewsState
  blockCountIdStr = GetCurrentOutputResultsBlockCountId()
  imageSetCount = gPipeAndViewsState.mImageSetCounter

  if 'plot basename' in inJsonToCopy:
    extraNameString = inJsonToCopy['plot basename']
  else:
    extraNameString = ''
  if inNumComponents != 3:
    newPlotOverTimeJson = inJsonToCopy.copy()
    newPlotOverTimeJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + ".b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".plot."
    newPlotOverTimeJson['plot basename'] = plotName
    newPlotOverTimeJson['variable scalar'] = inYAxisVariableName
    outPlotOverTimeJsonBlocks[plotName] = newPlotOverTimeJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newPlotOverTimeJson) + '\n')
  else:
    newPlotOverTimeJson = inJsonToCopy.copy()
    newPlotOverTimeJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + "x.b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".plot."
    newPlotOverTimeJson['plot basename'] = plotName
    newPlotOverTimeJson['variable vector component'] = inYAxisVariableName + 'x'
    outPlotOverTimeJsonBlocks[plotName] = newPlotOverTimeJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newPlotOverTimeJson) + '\n')
    newPlotOverTimeJson = inJsonToCopy.copy()
    newPlotOverTimeJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + "y.b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".plot."
    newPlotOverTimeJson['plot basename'] = plotName
    newPlotOverTimeJson['variable vector component'] = inYAxisVariableName + 'y'
    outPlotOverTimeJsonBlocks[plotName] = newPlotOverTimeJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newPlotOverTimeJson) + '\n')
    newPlotOverTimeJson = inJsonToCopy.copy()
    newPlotOverTimeJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + "z.b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".plot."
    newPlotOverTimeJson['plot basename'] = plotName
    newPlotOverTimeJson['variable vector component'] = inYAxisVariableName + 'z'
    outPlotOverTimeJsonBlocks[plotName] = newPlotOverTimeJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newPlotOverTimeJson) + '\n')
    newPlotOverTimeJson = inJsonToCopy.copy()
    newPlotOverTimeJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + "mag.b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".plot."
    newPlotOverTimeJson['plot basename'] = plotName
    lenToGrab = len(inYAxisVariableName) - len(GetSeparatorString())
    magVarName = inYAxisVariableName[0:lenToGrab]
    newPlotOverTimeJson['variable vector magnitude'] = magVarName
    outPlotOverTimeJsonBlocks[plotName] = newPlotOverTimeJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newPlotOverTimeJson) + '\n')


def ConstructAndAddScatterPlotsJson(outScatterPlotJsonBlocks, inJsonToCopy,
        inYAxisVariableName, inNumComponents, inXAxisVariableName):
  """helper function for TestForAndConstructAllScatterPlots, creates a
     scatterplot json block, not for external use"""

  global gPipeAndViewsState
  blockCountIdStr = GetCurrentOutputResultsBlockCountId()
  imageSetCount = gPipeAndViewsState.mImageSetCounter

  if 'plot basename' in inJsonToCopy:
    extraNameString = inJsonToCopy['plot basename']
  else:
    extraNameString = ''
  if inNumComponents != 3:
    newScatterPlotJson = inJsonToCopy.copy()
    newScatterPlotJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + ".b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".sctr."
    newScatterPlotJson['plot basename'] = plotName
    newScatterPlotJson['x axis variable scalar'] = inXAxisVariableName
    newScatterPlotJson['y axis variable scalar'] = inYAxisVariableName
    outScatterPlotJsonBlocks[plotName] = newScatterPlotJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newScatterPlotJson) + '\n')
  else:
    newScatterPlotJson = inJsonToCopy.copy()
    newScatterPlotJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + "x.b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".sctr."
    newScatterPlotJson['plot basename'] = plotName
    newScatterPlotJson['x axis variable scalar'] = inXAxisVariableName
    newScatterPlotJson['y axis variable vector component'] = inYAxisVariableName + 'x'
    outScatterPlotJsonBlocks[plotName] = newScatterPlotJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newScatterPlotJson) + '\n')
    newScatterPlotJson = inJsonToCopy.copy()
    newScatterPlotJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + "y.b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".sctr."
    newScatterPlotJson['plot basename'] = plotName
    newScatterPlotJson['x axis variable scalar'] = inXAxisVariableName
    newScatterPlotJson['y axis variable vector component'] = inYAxisVariableName + 'y'
    outScatterPlotJsonBlocks[plotName] = newScatterPlotJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newScatterPlotJson) + '\n')
    newScatterPlotJson = inJsonToCopy.copy()
    newScatterPlotJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + "z.b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".sctr."
    newScatterPlotJson['plot basename'] = plotName
    newScatterPlotJson['x axis variable scalar'] = inXAxisVariableName
    newScatterPlotJson['y axis variable vector component'] = inYAxisVariableName + 'z'
    outScatterPlotJsonBlocks[plotName] = newScatterPlotJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newScatterPlotJson) + '\n')
    newScatterPlotJson = inJsonToCopy.copy()
    newScatterPlotJson['all variables'] = False
    plotName = extraNameString + inYAxisVariableName + "mag.b-" + blockCountIdStr + ".is-" + str(imageSetCount) + ".sctr."
    newScatterPlotJson['plot basename'] = plotName
    newScatterPlotJson['x axis variable scalar'] = inXAxisVariableName
    lenToGrab = len(inYAxisVariableName) - len(GetSeparatorString())
    magVarName = inYAxisVariableName[0:lenToGrab]
    newScatterPlotJson['y axis variable vector magnitude'] = magVarName
    outScatterPlotJsonBlocks[plotName] = newScatterPlotJson
    if PhactoriDbg():
      myDebugPrint3('added block named ' + plotName + ': ' + \
          str(newScatterPlotJson) + '\n')

def TestForAndConstructAllPlotsOverTime(ioPlotOverTimeBlocksJson):
  """tests to see if a scatterplot block contains "all variables": true and,
     if so, construct a scatterplot block for each variable"""
  if PhactoriDbg(100):
    myDebugPrint3('TestForAndConstructAllPlotsOverTime entered\n', 100)
  potBlocks = ioPlotOverTimeBlocksJson
  blocksToPop = []
  constructAllVariablesFlag = False
  for blockName, blockParamsJson in potBlocks.iteritems():
    if PhactoriDbg():
      myDebugPrint3("testing block named " + blockName + " for all variables\n")
    if PhactoriDbg():
      myDebugPrint3("  json is: " + str(blockParamsJson) + "\n")
    testVal = getParameterFromBlock(blockParamsJson, "all variables", False)
    if testVal == True:
      constructAllVariablesFlag = True
      blockToPop = blockName
      blockToCopy = blockParamsJson
      #there should be only one like this
      break

  if constructAllVariablesFlag == False:
    #don't construct all variables
    if PhactoriDbg(100):
      myDebugPrint3('TestForAndConstructAllPlotsOverTime returning, no all variables\n', 100)
    return

  global gPipeAndViewsState
  paraviewSource = gPipeAndViewsState.mIncomingDefaultOperation.GetPvFilter()
  pointDataInfo = paraviewSource.GetPointDataInformation()
  cellDataInfo = paraviewSource.GetCellDataInformation()
  fieldDataInfo = paraviewSource.GetFieldDataInformation()

  numPointArrays = pointDataInfo.GetNumberOfArrays()
  for ii in range(numPointArrays):
    pointArray = pointDataInfo.GetArray(ii)
    pointArrayName = pointArray.GetName()
    if pointArrayName != "GlobalNodeId":
      ConstructAndAddPlotsOverTimeJson(potBlocks, blockToCopy,
          pointArrayName, pointArray.GetNumberOfComponents())
  numCellArrays = cellDataInfo.GetNumberOfArrays()
  for ii in range(numCellArrays):
    cellArray = cellDataInfo.GetArray(ii)
    cellArrayName = cellArray.GetName()
    if cellArrayName != "GlobalElementId":
      ConstructAndAddPlotsOverTimeJson(potBlocks, blockToCopy,
          cellArray.GetName(), cellArray.GetNumberOfComponents())
  numFieldArrays = fieldDataInfo.GetNumberOfArrays()
  for ii in range(numFieldArrays):
    fieldArray = fieldDataInfo.GetArray(ii)
    fieldArrayName = fieldArray.GetName()
    ConstructAndAddPlotsOverTimeJson(potBlocks, blockToCopy,
        fieldArray.GetName(), fieldArray.GetNumberOfComponents())

  gPipeAndViewsState.mImageSetCounter += 1

  potBlocks.pop(blockToPop, None)

  if PhactoriDbg(100):
    myDebugPrint3('TestForAndConstructAllPlotsOverTime returning\n', 100)


def TestForAndConstructAllScatterPlots(ioScatterPlotBlocksJson):
  """tests to see if a scatterplot block contains "all variables": true and,
     if so, construct a scatterplot block for each variable"""
  if PhactoriDbg(100):
    myDebugPrint3('TestForAndConstructAllScatterPlots entered\n', 100)
  spBlocks = ioScatterPlotBlocksJson
  blocksToPop = []
  constructAllVariablesFlag = False
  for blockName, blockParamsJson in spBlocks.iteritems():
    if PhactoriDbg():
      myDebugPrint3("testing block named " + blockName + " for all variables\n")
    if PhactoriDbg():
      myDebugPrint3("  json is: " + str(blockParamsJson) + "\n")
    testVal = getParameterFromBlock(blockParamsJson, "all variables", False)
    if testVal == True:
      constructAllVariablesFlag = True
      blockToPop = blockName
      blockToCopy = blockParamsJson
      #there should be only one like this
      break

  if constructAllVariablesFlag == False:
    #don't construct all variables
    if PhactoriDbg(100):
      myDebugPrint3('TestForAndConstructAllScatterPlots returning, no all variables\n', 100)
    return

  global gPipeAndViewsState
  paraviewSource = gPipeAndViewsState.mIncomingDefaultOperation.GetPvFilter()
  pointDataInfo = paraviewSource.GetPointDataInformation()
  cellDataInfo = paraviewSource.GetCellDataInformation()

  numPointArrays = pointDataInfo.GetNumberOfArrays()
  for ii in range(numPointArrays):
    pointArray = pointDataInfo.GetArray(ii)
    pointArrayName = pointArray.GetName()
    if pointArrayName != "GlobalNodeId":
      ConstructAndAddScatterPlotsJson(spBlocks, blockToCopy,
          pointArrayName, pointArray.GetNumberOfComponents(),
          "GlobalNodeId")
  numCellArrays = cellDataInfo.GetNumberOfArrays()
  for ii in range(numCellArrays):
    cellArray = cellDataInfo.GetArray(ii)
    cellArrayName = cellArray.GetName()
    if cellArrayName != "GlobalElementId":
      ConstructAndAddScatterPlotsJson(spBlocks, blockToCopy,
          cellArray.GetName(), pointArray.GetNumberOfComponents(),
          "GlobalElementId")

  gPipeAndViewsState.mImageSetCounter += 1

  spBlocks.pop(blockToPop, None)

  if PhactoriDbg(100):
    myDebugPrint3('TestForAndConstructAllScatterPlots returning\n', 100)

def UpdatePipeAndViewsStateFromUserA(inPipeAndViewsState, inFileName):
  """given a PhactoriPipeAndViewsState instance which has been previously
     parsed, paraview items created, etc. (and presumably used to create
     rendered images at least once), read in the json state from the
     specified file, and update the inPipeAndViewsState accordingly.
     Only do the update if the PhactoriUpdateTrigger.txt file has a new
     value since the last check.  Only read the file on processor zero,
     and use mpi to share the update with all other processors."""
  if PhactoriDbg():
    myDebugPrint3("UpdatePipeAndViewsStateFromUserA entered\n")
  if SmartGetLocalProcessId() != 0:
    if PhactoriDbg():
      myDebugPrint3("UpdatePipeAndViewsStateFromUserA returning, not process 0\n")
    return
  try:
    import json
    inFile = open(inFileName, 'rb')
    userUpdateJson = json.load(inFile)
    inFile.close()
    userUpdateJson = convertJsonUnicodeToStrings(userUpdateJson)

    if 'representation blocks' in userUpdateJson:
      representationBlocks = userUpdateJson['representation blocks']
      #could add logic to detect no change
      for repName, repJson in representationBlocks.iteritems():
        onePhactoriRep = inPipeAndViewsState.mRepresentationBlocks[repName]
        if PhactoriDbg():
          myDebugPrint3("representation name: " + str(repName) + \
              "\nrepresentation json: " + str(repJson) + "\n")
        if onePhactoriRep != None:
          ParseOneRepresentationBlockC(onePhactoriRep, repJson,
              inPipeAndViewsState)
          if PhactoriDbg():
            myDebugPrint3("new color by variable: " +
                onePhactoriRep.mColorVariableInfo.SelfToStr() + "\n")

    if 'camera blocks' in userUpdateJson:
      cameraBlocks = userUpdateJson['camera blocks']
      for cameraName, cameraJson in cameraBlocks.iteritems():
        onePhactoriCamera = inPipeAndViewsState.mCameraBlocks[cameraName]
        if PhactoriDbg():
          myDebugPrint3("camera name: " + str(cameraName) + \
              "\ncamera json: " + str(cameraJson) + "\n")
        if onePhactoriCamera != None:
          if PhactoriDbg():
            myDebugPrint3("updating camera named: " + str(cameraName) + "\n")
          ParseOneCameraBlockC(onePhactoriCamera, cameraJson,
              inPipeAndViewsState)
          if PhactoriDbg():
            myDebugPrint3("new look direction: " + \
                str(onePhactoriCamera.mLookDirection) + "\n")

  except:
    if PhactoriDbg():
      myDebugPrint3("error trying to get and use json from file\n")
  if PhactoriDbg():
    myDebugPrint3("UpdatePipeAndViewsStateFromUserA returning\n")

def SaveJsonSettingsToFile(inPipeAndViewsState, inFileName):
  """given a PhactoriPipeAndViewsState instance, dump the json that was
     used to construct it (usually called after default items have been
     added)"""
  if PhactoriDbg():
    myDebugPrint3("SaveJsonSettingsToFile entered\n")
  if SmartGetLocalProcessId() != 0:
    if PhactoriDbg():
      myDebugPrint3("SaveJsonSettingsToFile returning, not process 0\n")
    return
  try:
    import json
    outStr = json.dumps(inPipeAndViewsState.mJsonDescription, sort_keys=True,
        indent = 2, separators=(',', ': '))
    outFile = open(inFileName, 'w+b')
    outFile.write(outStr)
    outFile.close()
    outFile = open('PhactoriInteractionTrigger.txt', 'w+b')
    outStr = json.dumps({'InteractionTriggerCounter': 0,
        'InteractionState':'OncePerCallbackFromSimulation'})
    outFile.write(outStr)
    outFile.close()
  except:
    if PhactoriDbg():
      myDebugPrint3("error trying to create or write file or make json string\n")

  if PhactoriDbg():
    myDebugPrint3("SaveJsonSettingsToFile returning\n")

def PreprocessExperimentalTextAnnotationBlock(
        inViewMapC, inExpBlkname, inExpBlk):

  substringToSubstituteForSpace = None
  if 'substringtosubstituteforspace' in inExpBlk:
    substringToSubstituteForSpace = \
            str(inExpBlk['substringtosubstituteforspace'])
    if PhactoriDbg():
      myDebugPrint3("substringtosubstituteforspace: ->" + \
              substringToSubstituteForSpace + "<-\n")

  if 'line1' not in inExpBlk:
    myDebugPrint3AndException("PreprocessOneExperimentalBlock:\n" \
            "textannotation type block must have 'line1' entry")

  myNewTxtAntnBlk = {}

  for lineIndex in range(1, 10):
    oneLineString = ""
    lineKey = "line" + str(lineIndex)
    if lineKey not in inExpBlk:
      break
    oneLineString = str(inExpBlk[lineKey])
    if substringToSubstituteForSpace != None:
      oneLineString = oneLineString.replace(
              substringToSubstituteForSpace, " ")
      if PhactoriDbg():
        myDebugPrint3("line after sub: ->" + oneLineString + "\n")
    myNewTxtAntnBlk[lineKey] = oneLineString

  if 'windowlocation' in inExpBlk:
    myNewTxtAntnBlk['windowlocation'] = inExpBlk['windowlocation']

  if ('positionx' in inExpBlk) or ('positiony' in inExpBlk):
    if 'positiony' not in inExpBlk:
      myDebugPrint3AndException("PreprocessOneExperimentalBlock:\n" \
              "if you have 'positionx' you need 'positiony'\n")
    if 'positionx' not in inExpBlk:
      myDebugPrint3AndException("PreprocessOneExperimentalBlock:\n" \
              "if you have 'positiony' you need 'positionx'\n")
    posx = float(inExpBlk['positionx'])
    posy = float(inExpBlk['positiony'])
    myNewTxtAntnBlk['windowlocation'] = 'AnyLocation'
    myNewTxtAntnBlk['position'] = [posx, posy]

  if 'color' in inExpBlk:
    colorRgb = inExpBlk['color']
    myNewTxtAntnBlk['color'] = \
            [float(colorRgb[0]), float(colorRgb[1]), float(colorRgb[2])]

  if 'opacity' in inExpBlk:
    myNewTxtAntnBlk['opacity'] = float(inExpBlk['opacity'])

  if 'boldflag' in inExpBlk:
    myNewTxtAntnBlk['boldflag'] = int(inExpBlk['boldflag'])
  if 'italicflag' in inExpBlk:
    myNewTxtAntnBlk['italicflag'] = int(inExpBlk['italicflag'])
  if 'shadowflag' in inExpBlk:
    myNewTxtAntnBlk['shadowflag'] = int(inExpBlk['shadowflag'])

  if 'fontsize' in inExpBlk:
    myNewTxtAntnBlk['fontsize'] = int(inExpBlk['fontsize'])
  if 'fontfamily' in inExpBlk:
    myNewTxtAntnBlk['fontfamily'] = inExpBlk['fontfamily']

  if 'text annotation blocks' not in inViewMapC:
    inViewMapC['text annotation blocks'] = {}
  inViewMapC['text annotation blocks'][inExpBlkname] = myNewTxtAntnBlk

  if 'show_in_imageset' in inExpBlk:
    targetImageset = inExpBlk['show_in_imageset']
    imagesets = inViewMapC['imageset blocks']
    if targetImageset in imagesets:
        theImageset = imagesets[targetImageset]
        if 'text annotations' not in theImageset:
          theImageset['text annotations'] = []
        textAnnotationList = theImageset['text annotations']
        textAnnotationList.append(inExpBlkname)
    else:
      myDebugPrint3AndException("PreprocessOneExperimentalBlock:\n" \
          "textannotation target imageset is not in imageset blocks")
  else:
    myDebugPrint3AndException("PreprocessOneExperimentalBlock:\n" \
        "textannotation type block must have 'show_in_imageset' entry")


def PreprocessExperimentalAddOpacityToRepBlock(
        inViewMapC, inExpBlkname, inExpBlk):

  if 'representation blocks' not in inViewMapC:
    myDebugPrint3AndException("PreprocessExperimentalAddOpacityToRepBlock:\n" \
        "no 'representation blocks'\n")

  repNameToAddTo = inExpBlk["target_representation"]
  if repNameToAddTo == None:
    myDebugPrint3AndException("PreprocessExperimentalAddOpacityToRepBlock:\n" \
        + str(inExpBlkname) + " missing 'target_representation'\n")

  representationBlocks = inViewMapC["representation blocks"]
  targetRep =  representationBlocks[repNameToAddTo]
  if targetRep == None:
    myDebugPrint3AndException("PreprocessExperimentalAddOpacityToRepBlock:\n" \
      "no representation with name: " + str(repNameToAddTo) + "\n");

  opacitySetting = inExpBlk["opacity"]
  if opacitySetting == None:
    myDebugPrint3AndException("PreprocessExperimentalAddOpacityToRepBlock:\n" \
        + str(inExpBlkname) + " missing 'opacity'\n")

  targetRep["opacity"] = float(opacitySetting)
  if PhactoriDbg():
      myDebugPrint3("PreprocessExperimentalAddOpacityToRepBlock:\n" \
        "to representation: " + repNameToAddTo + "\n" \
        "added opacity: " + str(opacitySetting) + "\n")


def PreprocessOneExperimentalBlock(inViewMapC, inExpBlkname, inExpBlk):
  """gives us entry point to take one experimental blocks and manipulate it
     to produce other kinds of blocks to enter into inViewMapC"""
  expBlkType = inExpBlk["type"]

  if expBlkType == "textannotation":
    PreprocessExperimentalTextAnnotationBlock(
            inViewMapC, inExpBlkname, inExpBlk)

  if expBlkType == "addopacitysettingtorepresentation":
    PreprocessExperimentalAddOpacityToRepBlock(inViewMapC, inExpBlkname, inExpBlk)


def PreprocessExperimentalBlocks(inViewMapC):
  """gives us entry point to take experimental blocks and manipulate them
     to produce other kinds of blocks to enter into inViewMapC"""

  if 'experimental blocks' not in inViewMapC:
    #no experimental blocks
    return

  expBlkMap = inViewMapC['experimental blocks']
  for expBlkName, expBlk in expBlkMap.iteritems():
      PreprocessOneExperimentalBlock(inViewMapC, expBlkName, expBlk)

  if len(expBlkMap) > 0:
    if PhactoriDbg():
      myDebugPrint3(
        "altered inViewMapC after experimental block handling:\n" + \
                str(inViewMapC) + "\n")

#construct a set of views from the json-structure view/pipeline description
#specified in the sierra catalyst wiki/documentation
def CreateViewSetFromPhactoriViewMapC(inViewMapC):
  """Create data pipeline, cameras, representations, and views from data struct

  Given a dict which is a json format structure with syntax as discussed in
  the sierra catalyst paraview insitu wiki, create the corresponding data
  pipeline, representations, views, etc.  The basic idea is that from the
  incoming datastructure, we will parse out 6 (currently) sets of blocks:
  operations, cameras, representations, imagesets, scatter plots, and plots
  over time.  From these, we will construct ParaView/Catalyst data structures
  to do the data management and rendering.  See the class
  PhactoriPipeAndViewsState and the lower level related classes for more
  explanation.
  """

  if PhactoriDbg(100):
    myDebugPrint3('CreateViewSetFromPhactoriViewMapC entered\n', 100)
  if PhactoriDbg():
    myDebugPrint3(str(inViewMapC) + "\n")
  global gPipeAndViewsState
  gPipeAndViewsState.mJsonDescription = inViewMapC

  #get pointers to the various block types (currently 6, camera,
  #representation, operation, imageset, plot over time, and scatter plot)

  #check for 'experimental blocks' and handle them appropriately if they
  #are there (e.g. use them to create other blocks)
  PreprocessExperimentalBlocks(inViewMapC)

  cameraBlocks = {}
  representationBlocks = {}
  imagesetBlocks = {}
  operationBlocks = {}
  scatterplotBlocks = {}
  timeplotBlocks = {}
  criteriaBlocks = {}
  markerBlocks = {}
  textannotationBlocks = {}
  if 'camera blocks' in inViewMapC:
    cameraBlocks = inViewMapC['camera blocks']
  if 'representation blocks' in inViewMapC:
    representationBlocks = inViewMapC['representation blocks']
  if 'imageset blocks' in inViewMapC:
    imagesetBlocks = inViewMapC['imageset blocks']
  if 'operation blocks' in inViewMapC:
    operationBlocks = inViewMapC['operation blocks']
  if 'scatter plot blocks' in inViewMapC:
    scatterplotBlocks = inViewMapC['scatter plot blocks']
  if 'plot over time blocks' in inViewMapC:
    timeplotBlocks = inViewMapC['plot over time blocks']
  if 'onoff criteria blocks' in inViewMapC:
    criteriaBlocks = inViewMapC['onoff criteria blocks']
  if 'visual marker blocks' in inViewMapC:
    markerBlocks = inViewMapC['visual marker blocks']
  if 'text annotation blocks' in inViewMapC:
    textannotationBlocks = inViewMapC['text annotation blocks']
  if PhactoriDbg():
    myDebugPrint3("  cameraBlocks:\n")
    myDebugPrint3("  " + str(cameraBlocks) + "\n")
    myDebugPrint3("  representationBlocks:\n")
    myDebugPrint3("  " + str(representationBlocks) + "\n")
    myDebugPrint3("  imagesetBlocks:\n")
    myDebugPrint3("  " + str(imagesetBlocks) + "\n")
    myDebugPrint3("  operationBlocks:\n")
    myDebugPrint3("  " + str(operationBlocks) + "\n")
    myDebugPrint3("  scatterplotBlocks:\n")
    myDebugPrint3("  " + str(scatterplotBlocks) + "\n")
    myDebugPrint3("  timeplotBlocks:\n")
    myDebugPrint3("  " + str(timeplotBlocks) + "\n")
    myDebugPrint3("  onoff criteria blocks:\n")
    myDebugPrint3("  " + str(criteriaBlocks) + "\n")
    myDebugPrint3("  visual marker blocks:\n")
    myDebugPrint3("  " + str(markerBlocks) + "\n")
    myDebugPrint3("  text annotation blocks:\n")
    myDebugPrint3("  " + str(textannotationBlocks) + "\n")

  #go through imageset json, and locate operation shortcuts.  If found, create
  #json definition of operations and add them to the structure for later
  #creation, and add reference to the new operations in the imageblock json
  HandleOperationShortcuts(imagesetBlocks, operationBlocks)

  #do pipeline stuff
  gPipeAndViewsState.mIncomingDefaultOperation.mParaViewFilter = GetActiveSource()

  global gEnableTemporaryExtractBlockTest
  if gEnableTemporaryExtractBlockTest:
    #hack to test extract block filter
    if gExtractBlockAOperationName in operationBlocks:
      global gExtractBlockAList
      theBlock = operationBlocks[gExtractBlockAOperationName]
      theBlock['type'] = 'extractblock'
      theBlock['include blocks'] = gExtractBlockAList
    if gExtractBlockBOperationName in operationBlocks:
      global gExtractBlockBList
      theBlock = operationBlocks[gExtractBlockBOperationName]
      theBlock['type'] = 'extractblock'
      theBlock['include blocks'] = gExtractBlockBList

  #to test imagesets with multiple operations visible
  if gEnableTemporaryMultiOpViewTest:
    if gMultiOpViewTestImagesetName in imagesetBlocks:
      theBlock = imagesetBlocks["multioptestimageset"]
      theBlock["operation2"] = gMultiOpViewTestOp2Name
      theBlock["representation2"] = gMultiOpViewTestRep2Name

  global gEnableTemporaryOffAxisProjectionTest
  if gEnableTemporaryOffAxisProjectionTest:
    global gOaptCameraNameList

    #go through cameras and set up stereo ones
    for oneOaptCameraName in gOaptCameraNameList:
      if oneOaptCameraName in cameraBlocks:
        global gOaptPhysicalSettingsForCamera
        physicalSettingsForCamera = \
            gOaptPhysicalSettingsForCamera[oneOaptCameraName]
        theBlock = cameraBlocks[oneOaptCameraName]
        theBlock['use off axis projection'] = True
        for key, value in physicalSettingsForCamera.iteritems():
          theBlock[key] = value
        theBlock['which eye'] = 'left'
        #copy this camera to another to use for right eye
        rightCameraName = oneOaptCameraName + "Right"
        rightCameraBlock = {}
        for key, value in theBlock.iteritems():
          rightCameraBlock[key] = value
        rightCameraBlock['which eye'] = 'right'
        cameraBlocks[rightCameraName] = rightCameraBlock
        if PhactoriDbg():
          myDebugPrint3("new camera block for right eye:\n"
            + str(rightCameraBlock) + "\n")

    #go through imageset blocks and add extras for stereo
    blocksToReplicateForRightEye = []
    for imagesetName, oneImagesetBlock in imagesetBlocks.iteritems():
      if 'camera' in oneImagesetBlock:
        if oneImagesetBlock['camera'] in gOaptCameraNameList:
          blocksToReplicateForRightEye.append(
              [imagesetName, oneImagesetBlock])
    if PhactoriDbg():
      myDebugPrint3("blocksToReplicateForRightEye:\n"
        + str(blocksToReplicateForRightEye) + "\n")
    for iterItem in blocksToReplicateForRightEye:
      oneImagesetName = iterItem[0]
      oneImagesetBlock = iterItem[1]
      newImagesetBlock = {}
      for key, value in oneImagesetBlock.iteritems():
        newImagesetBlock[key] = value
      rightCameraName = oneImagesetBlock['camera'] + "Right"
      newImagesetBlock['camera'] = rightCameraName
      imagesetBaseName = str(oneImagesetBlock['image basename'])
      oneImagesetBlock['image basename'] = imagesetBaseName + "left."
      newImagesetBlock['image basename'] = imagesetBaseName + "right."
      if PhactoriDbg():
        myDebugPrint3("new imageset block for right eye:\n"
          + str(newImagesetBlock) + "\n")
      imagesetBlocks[oneImagesetName + "Right"] = newImagesetBlock

  #swap out operations from input deck for testing operations which aren't
  #available through the input deck syntax yet
  global gEnableSubstituteOperationTesting
  if gEnableSubstituteOperationTesting:
    global gSubstituteOperationTestingList
    global gSubstituteOperationTestingMap
    for oneOpName in gSubstituteOperationTestingList:
      if oneOpName in operationBlocks:
        if PhactoriDbg():
          myDebugPrint3("substituting operation named: " + oneOpName + "\n")
        oneOpBlock = operationBlocks[oneOpName]
        subOpBlock = gSubstituteOperationTestingMap[oneOpName]
        if 'input' in oneOpBlock:
          subOpBlock['input'] = oneOpBlock['input']
        operationBlocks[oneOpName] = subOpBlock

  #swap out imagesets from input deck for testing imagesets which aren't
  #available through the input deck syntax yet
  global gEnableSubstituteImagesetTesting
  if gEnableSubstituteImagesetTesting:
    global gSubstituteImagesetTestingList
    global gSubstituteImagesetTestingMap
    for oneIsName in gSubstituteImagesetTestingList:
      if oneIsName in imagesetBlocks:
        oneIsBlock = imagesetBlocks[oneIsName]
        subIsBlock = gSubstituteImagesetTestingMap[oneIsName]
        for key, value in subIsBlock.iteritems():
          oneIsBlock[key] = value

  MakeFiltersFromViewMapOperationsC(gPipeAndViewsState, operationBlocks)
  if PhactoriDbg():
    myDebugPrint3('  gPipeAndViewsState.mOperationBlocks: ' + \
            str(gPipeAndViewsState.mOperationBlocks) + '\n')

  #parse scatter plot blocks

  #special case; check for 'all variables' in a scatter plot and construct all
  #plot json blocks if necessary
  TestForAndConstructAllScatterPlots(scatterplotBlocks)

  #ParseScatterPlotBlocksC(gPipeAndViewsState, scatterplotBlocks)
  #ParseBlocksC(gPipeAndViewsState.mScatterPlotBlocks, scatterplotBlocks,
  #    PhactoriScatterPlotBlock, ParseOneScatterPlotBlockC,
  #    gPipeAndViewsState)
  ParseBlocksC2(gPipeAndViewsState.mScatterPlotBlocks, scatterplotBlocks,
      PhactoriScatterPlotBlock, gPipeAndViewsState)

  #special case; check for 'all variables' in a plot over time and construct all
  #plot json blocks if necessary
  TestForAndConstructAllPlotsOverTime(timeplotBlocks)

  #parse time plot blocks
  ParseBlocksC(gPipeAndViewsState.mPlotOverTimeBlocks, timeplotBlocks,
      PhactoriPlotOverTimeBlock, ParseOnePlotOverTimeBlockC,
      gPipeAndViewsState)

  #create producers for plots
  CreateScatterPlotProducersC(gPipeAndViewsState)
  CreatePlotOverTimeProducersC(gPipeAndViewsState)

  #parse camera blocks
  #ParseCameraBlocksC(gPipeAndViewsState, cameraBlocks)
  ParseBlocksC(gPipeAndViewsState.mCameraBlocks, cameraBlocks,
      PhactoriCameraBlock, ParseOneCameraBlockC,
      gPipeAndViewsState)

  #parse representation blocks
  #ParseRepresentationBlocksC(gPipeAndViewsState, representationBlocks)
  ParseBlocksC(gPipeAndViewsState.mRepresentationBlocks, representationBlocks,
      PhactoriRepresentationBlock, ParseOneRepresentationBlockC,
      gPipeAndViewsState)

  #parse imageset onoff criteria blocks
  ParseBlocksC(gPipeAndViewsState.mImagesetOnOffCriteriaBlocks, criteriaBlocks,
      PhactoriImagesetOnOffFilterCriteria, ParseOneCriteriaBlockC,
      gPipeAndViewsState)

  #parse visual marker blocks
  numMarkerBlocks = ParseBlocksC(
      gPipeAndViewsState.mMarkerBlocks,
      markerBlocks, PhactoriMarkerBlock,
      ParseOneMarkerBlockC, gPipeAndViewsState)

  #parse text annotation blocks
  numTextAnnotationBlocks = ParseBlocksC(
      gPipeAndViewsState.mTextAnnotationBlocks,
      textannotationBlocks, PhactoriTextAnnotationBlock,
      ParseOneTextAnnotationBlockC, gPipeAndViewsState)

  #construct marker paraview items (e.g. sphere or box source and the
  #corresponding representation
  #this is now done the first time the marker is made visible
  #for oneMarkerName, oneMarkerBlock in \
  #        gPipeAndViewsState.mMarkerBlocks.iteritems():
  #  oneMarkerBlock.CreateParaViewItems()

  #parse imageset blocks
  #ParseImagesetBlocksC(gPipeAndViewsState, imagesetBlocks)
  numImagesetBlocks = ParseBlocksC(gPipeAndViewsState.mImagesetBlocks, imagesetBlocks,
      PhactoriImagesetBlock, ParseOneImagesetBlockC,
      gPipeAndViewsState)

  #if there are no imageset blocks, create default and parse parameters in
  #wrappint imagesetBlock json as if they were from an imageset block
  if numImagesetBlocks == 0:
    if PhactoriDbg():
      myDebugPrint3("  no imagesets, creating a default imageset\n")
    imagesetBlocks['defaultimageset'] = {}
    ParseBlocksC(gPipeAndViewsState.mImagesetBlocks, imagesetBlocks,
        PhactoriImagesetBlock, ParseOneImagesetBlockC,
        gPipeAndViewsState)


#    if 'camera' not in ioImagesetBlockJson:
#      #we have to construct and use a default camera, including parsing
#      #commands in the imageset for the camera
#      myDebugPrint3("  ParseOneImagesetBlockC: for imageset " + \
#          ioImagesetBlock.mName + \
#          " there is no camera, so we must add and reference default\n")
#      if 'camera type' not in ioImagesetBlockJson:
#        ioImagesetBlockJson['camera type'] = 'multicamera8'
#      defaultCameraName = ioImagesetBlock.mName + '_default_camera'
#      defaultCameraBlockAndWrapper = {defaultCameraName: ioImagesetBlockJson}
#      ParseBlocksC(ioPipeAndViewsState.mCameraBlocks,
#          defaultCameraBlockAndWrapper,
#          PhactoriCameraBlock,
#          ParseOneCameraBlockC,
#          ioPipeAndViewsState)
#      ioImagesetBlockJson['camera'] = defaultCameraName
#      #myDebugPrint3(  "done adding camera, here it is:\n")
#      #ioPipeAndViewsState.mCameraBlocks[defaultCameraName].PrintSelf()


  #loop through imageset blocks and create a viz for each--i.e. create a
  #ParaView/Catalyst Representation and View for each imageset block, and
  #have it render info from some stage of the above pipeline

  for imagesetName, imagesetInstance in gPipeAndViewsState.mImagesetBlocks.iteritems():
    CreateParaviewItemsForImagesetC(imagesetInstance)

  CreateScatterPlotViewsC(gPipeAndViewsState)
  CreatePlotOverTimeViewsC(gPipeAndViewsState)

  if gPipeAndViewsState.mInteractionEnabled == True:
    SaveJsonSettingsToFile(gPipeAndViewsState, "PhactoriJsonSettings.txt")

def CreateScatterPlotProducersC(inPipeAndViewsState):
  if PhactoriDbg(100):
    myDebugPrint3('CreateScatterPlotProducersC entered\n', 100)
  for oneScatterPlotName, oneScatterPlot in inPipeAndViewsState.mScatterPlotBlocks.iteritems():
    #myDebugPrint3('  creating scatter plot producer nameed ' + oneScatterPlot.mName + '\n')
    CreateOneScatterPlotProducerC(oneScatterPlot)
  if PhactoriDbg(100):
    myDebugPrint3('CreateScatterPlotProducersC returning\n', 100)

def CreatePlotOverTimeProducersC(inPipeAndViewsState):
  if PhactoriDbg(100):
    myDebugPrint3('CreatePlotOverTimeProducersC entered\n', 100)
  for onePlotOverTimeName, onePlotOverTime in inPipeAndViewsState.mPlotOverTimeBlocks.iteritems():
    #myDebugPrint3('  creating plot over time producer nameed ' + onePlotOverTime.mName + '\n')
    CreateOnePlotOverTimeProducerC(onePlotOverTime)
  if PhactoriDbg(100):
    myDebugPrint3('CreatePlotOverTimeProducersC returning\n', 100)

def CreateScatterPlotViewsC(inPipeAndViewsState):

  if PhactoriDbg(100):
    myDebugPrint3('CreateScatterPlotViewsC entered\n', 100)
  #SetPlotView2StartColors(
  #  inBackgroundColor = [1.0, 1.0, 1.0],
  #  inEdgeColor = [0.0, 0.0, 0.5],
  #  inCubeAxesColor = [0.2, 0.2, 0.2],
  #  inDiffuseColor = [0.2, 0.2, 0.2],
  #  inAmbientColor = [0.2, 0.2, 0.2],
  #  inSelectionColor = [0.2, 0.2, 0.2],
  #  inBackfaceDiffuseColor = [0.2, 0.2, 0.2])

  for oneScatterPlotName, oneScatterPlot in inPipeAndViewsState.mScatterPlotBlocks.iteritems():
    #myDebugPrint3('  creating scatter plot view nameed ' + oneScatterPlot.mName + '\n')
    CreateOneScatterPlotViewC(oneScatterPlot)
  if PhactoriDbg(100):
    myDebugPrint3('CreateScatterPlotViewsC returning\n', 100)


def CreatePlotOverTimeViewsC(inPipeAndViewsState):

  if PhactoriDbg(100):
    myDebugPrint3('CreatePlotOverTimeViewsC entered\n', 100)
  #SetPlotView2StartColors(
  #  inBackgroundColor = [1.0, 1.0, 1.0],
  #  inEdgeColor = [0.0, 0.0, 0.5],
  #  inCubeAxesColor = [0.2, 0.2, 0.2],
  #  inDiffuseColor = [0.2, 0.2, 0.2],
  #  inAmbientColor = [0.2, 0.2, 0.2],
  #  inSelectionColor = [0.2, 0.2, 0.2],
  #  inBackfaceDiffuseColor = [0.2, 0.2, 0.2])

  for onePlotOverTimeName, onePlotOverTime in inPipeAndViewsState.mPlotOverTimeBlocks.iteritems():
    #myDebugPrint3('  creating scatter plot view nameed ' + onePlotOverTime.mName + '\n')
    CreateOnePlotOverTimeViewC(onePlotOverTime)
  if PhactoriDbg(100):
    myDebugPrint3('CreatePlotOverTimeViewsC returning\n', 100)


def ColorByVariableComponentOrMagnitudeXX(inParaViewDataRepresentation,
        inPhactoriRepresentation,
        inVariableName, inComponentOrMagnitude, inWhichComponent,
        inVariableArrayType, inColorMapSettings):
  #myDebugPrint3("ColorByVariableComponentOrMagnitudeXX entered\n", 100)
  #myDebugPrint3("  inVariableName: " + inVariableName + "\n")
  #myDebugPrint3("  inComponentOrMagnitude: " + inComponentOrMagnitude + "\n")
  #myDebugPrint3("  inWhichComponent: " + str(inWhichComponent) + "\n")
  #myDebugPrint3("  inVariableArrayType: " + str(inVariableArrayType) + "\n")

  if gParaViewCatalystVersionFlag >= 40300:
    #for paraview 4.3
    if inVariableArrayType == None:
      myDebugPrint3AndException(\
        "ColorByVariableComponentOrMagnitudeXX "
        "expecting inVariableArrayType\n")
    if PhactoriDbg():
      myDebugPrint3("using ColorBy() to set color var " + \
        inVariableName + "  " + inVariableArrayType + " 1\n")
    ColorBy(inParaViewDataRepresentation,
        (inVariableArrayType, inVariableName))
    #deal with vector mode, component index
    pv_4_3_LUT = GetColorTransferFunction(inVariableName)
    if inComponentOrMagnitude == 'Magnitude':
      pv_4_3_LUT.VectorMode = 'Magnitude'
    else:
      pv_4_3_LUT.VectorMode = 'Component'
      pv_4_3_LUT.VectorComponent = inWhichComponent
    #need to deal with color table
    myRGBPoints, myColorSpace, myNanColor = \
        GetColorMapInfoFromColorLegendCollection(inColorMapSettings,
            -1.0, 1.0)
    pv_4_3_LUT.RGBPoints = myRGBPoints
    pv_4_3_LUT.ColorSpace = myColorSpace
    pv_4_3_LUT.NanColor = myNanColor
    return;

  #after this the function is for paraview 4.1

  if inVariableArrayType != None:
      inParaViewDataRepresentation.ColorAttributeType = inVariableArrayType

  myRGBPoints, myColorSpace, myNanColor = \
      GetColorMapInfoFromColorLegendCollection(inColorMapSettings, -1.0, 1.0)

  newLookupTable = None
  if inComponentOrMagnitude == 'Magnitude':
    newLookupTable = GetLookupTableForArray( inVariableName, 3, Discretize=1, RGBPoints=myRGBPoints, UseLogScale=0, VectorComponent=0, NanColor=myNanColor, NumberOfTableValues=1024, ColorSpace=myColorSpace, VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
  elif inComponentOrMagnitude == 'Component':
    newLookupTable = GetLookupTableForArray( inVariableName, 3, Discretize=1, RGBPoints=myRGBPoints, UseLogScale=0, VectorComponent=inWhichComponent, NanColor=myNanColor, NumberOfTableValues=1024, ColorSpace=myColorSpace, VectorMode='Component', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
  else:
    #bad name, should be 'Component' or 'Magnitude'
    return

  new_var_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 1.0, 1.0] )
  inParaViewDataRepresentation.LookupTable = newLookupTable
  try:
    inParaViewDataRepresentation.ScalarOpacityFunction = new_var_PiecewiseFunction
  except:
    if PhactoriDbg():
      myDebugPrint3('no inParaViewDataRepresentation.ScalarOpacityFunction\n');

  inParaViewDataRepresentation.ColorArrayName = inVariableName


def ColorByVariableScalarXX(inParaViewDataRepresentation, inPhactoriRepresentation, inVariableName, inVariableArrayType, inColorMapSettings):
  """Set the scalar variable (by name) which will be used to color the data
     presentation."""
  if gParaViewCatalystVersionFlag >= 40300:
    #for paraview 4.3
    if inVariableArrayType == None:
      myDebugPrint3AndException(\
        "ColorByVariableComponentOrMagnitudeXX "
        "expecting inVariableArrayType\n")
    if PhactoriDbg():
      myDebugPrint3("using ColorBy() to set color var " + \
        inVariableName + "  " + inVariableArrayType + " 3\n")
    ColorBy(inParaViewDataRepresentation,
        (inVariableArrayType, inVariableName))
    #need to deal with color table
    myRGBPoints, myColorSpace, myNanColor = \
        GetColorMapInfoFromColorLegendCollection(inColorMapSettings,
            -1.0, 1.0)
    pv_4_3_LUT = GetColorTransferFunction(inVariableName)
    pv_4_3_LUT.RGBPoints = myRGBPoints
    pv_4_3_LUT.ColorSpace = myColorSpace
    pv_4_3_LUT.NanColor = myNanColor
    return

  #rest of function for paraview 4.1

  inParaViewDataRepresentation.ColorArrayName = inVariableName
  if inVariableArrayType != None:
    inParaViewDataRepresentation.ColorAttributeType = inVariableArrayType

  new_var_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 1.0, 1.0] )

  myRGBPoints, myColorSpace, myNanColor = \
      GetColorMapInfoFromColorLegendCollection(inColorMapSettings,
          0.0, 1647.1430997030882)

  new_var_PVLookupTable = GetLookupTableForArray( inVariableName, 3, Discretize=1, RGBPoints= myRGBPoints, UseLogScale=0, VectorComponent=0, NanColor=myNanColor, NumberOfTableValues=1024, ColorSpace=myColorSpace, VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
  if PhactoriDbg():
    myDebugPrint3("new lookup table (1):\n");
  if PhactoriDbg():
    myDebugPrint3(str(new_var_PVLookupTable) + "\n")
  if PhactoriDbg():
    myDebugPrint3("RGBPoints:\n");
  if PhactoriDbg():
    myDebugPrint3(str(new_var_PVLookupTable.RGBPoints) + "\n")

  if PhactoriDbg():
    myDebugPrint3(str(inParaViewDataRepresentation) + '\n')
  try:
    inParaViewDataRepresentation.ScalarOpacityFunction = new_var_PiecewiseFunction
  except:
    if PhactoriDbg():
      myDebugPrint3('no inParaViewDataRepresentation.ScalarOpacityFunction\n');
  if PhactoriDbg():
    myDebugPrint3('new or reused lookup table: ' + str(new_var_PVLookupTable) + '\n')
  if PhactoriDbg():
    myDebugPrint3('prior inParaViewDataRepresentation.DataRepresentation.LookupTable: ' + str(inParaViewDataRepresentation.LookupTable) + '\n')
  inParaViewDataRepresentation.LookupTable = new_var_PVLookupTable
  if PhactoriDbg():
    myDebugPrint3(str(new_var_PVLookupTable))
  if PhactoriDbg():
    myDebugPrint3(str(new_var_PVLookupTable.RGBPoints))
  if PhactoriDbg():
    myDebugPrint3(str(inParaViewDataRepresentation) + '\n')

global gNumCycleBlockColors
global gCycleBlockColors

gNumCycleBlockColors = 12L
gCycleBlockColors = [[1.00, 0.00, 0.00],
                     [1.00, 1.00, 0.00],
                     [0.00, 1.00, 0.00],
                     [0.00, 1.00, 1.00],
                     [0.00, 0.00, 1.00],
                     [1.00, 0.00, 1.00],
                     [1.00, 0.50, 0.00],
                     [0.50, 1.00, 0.00],
                     [0.00, 1.00, 0.50],
                     [0.00, 0.50, 1.00],
                     [0.50, 0.00, 1.00],
                     [1.00, 0.00, 0.50]]

def ColorByBlockOneBlock(inInputCsData, ioBlockAndLeafBlockCounter, ioBlockColorData, inCountOnlyFlag):
  #myDebugPrint3('ColorByBlockOneBlock entered\n', 100)
  ioBlockAndLeafBlockCounter[1] += 1L
  #myDebugPrint3('  ioBlockAndLeafBlockCounter: ' + str(ioBlockAndLeafBlockCounter) + '\n')
  if inCountOnlyFlag == False:
    global gNumCycleBlockColors
    global gCycleBlockColors
    colorIndex = (ioBlockAndLeafBlockCounter[1] - 1L) % gNumCycleBlockColors
    #myDebugPrint3('  colorIndex: ' + str(colorIndex) + '  color: ' + str(gCycleBlockColors[colorIndex]) + '\n')
    blockIndex = ioBlockAndLeafBlockCounter[0] - 1L
    ioBlockColorData[blockIndex] = gCycleBlockColors[colorIndex]

  #myDebugPrint3('ColorByBlockOneBlock returning\n', 100)

def ColorByBlockRecurse1(inInputCsData,
        ioBlockAndLeafBlockCounter,
        ioBlockColorData,
        inCountOnlyFlag = False,
        inMaximumCountFlag = False,
        inMaximumCount = 10000):
  #myDebugPrint3('ColorByBlockRecurse1 entered\n', 100)

  ioBlockAndLeafBlockCounter[0] += 1L

  if inMaximumCountFlag:
    #myDebugPrint3("  inMaximumCountFlag true, inMaximumCount: " + str(inMaximumCount) + "  block count: " + str(ioBlockAndLeafBlockCounter[0]) + "\n")
    if ioBlockAndLeafBlockCounter[0] >= inMaximumCount:
      #myDebugPrint3("  inMaximumCount hit, returning early\n", 100)
      return

  icsdClassname = inInputCsData.GetClassName()
  if icsdClassname == "vtkMultiBlockDataSet" or \
     icsdClassname == "vtkExodusIIMultiBlockDataSet":
    #myDebugPrint3('recursing: ' + icsdClassname + '\n')
    numBlocks = inInputCsData.GetNumberOfBlocks()
    for ii in range(0, numBlocks):
      oneBlock = inInputCsData.GetBlock(ii)
      if(oneBlock != None):
        ColorByBlockRecurse1(oneBlock, ioBlockAndLeafBlockCounter, ioBlockColorData, inCountOnlyFlag, inMaximumCountFlag, inMaximumCount)
      else:
        ioBlockAndLeafBlockCounter[0] += 1L
      if inMaximumCountFlag:
        #myDebugPrint3("  inMaximumCountFlag (2) true, inMaximumCount: " + str(inMaximumCount) + "  block count: " + str(ioBlockAndLeafBlockCounter[0]) + "\n")
        if ioBlockAndLeafBlockCounter[0] >= inMaximumCount:
          #myDebugPrint3("  inMaximumCount hit (2), returning early\n", 100)
          return
  else:
    ColorByBlockOneBlock(inInputCsData, ioBlockAndLeafBlockCounter, ioBlockColorData, inCountOnlyFlag)

  #myDebugPrint3('  ioBlockAndLeafBlockCounter after this recursion: ' + str(ioBlockAndLeafBlockCounter) + '\n')

  #myDebugPrint3('ColorByBlockRecurse1 returning\n', 100)

def ColorByBlock(inParaViewDataSource, inParaViewDataRepresentation,
                 inBlockCountLimitOverride = False):
  if PhactoriDbg(100):
    myDebugPrint3('ColorByBlock entered\n', 100)

  csdata = inParaViewDataSource.GetClientSideObject().GetOutputDataObject(0)
  blockAndLeafBlockCounter = [0L,0L]

  #test for extreme case when there are too many blocks
  #tooManyBlocksToColor = 8
  tooManyBlocksToColor = 10000
  ColorByBlockRecurse1(csdata, blockAndLeafBlockCounter, None, True, True, tooManyBlocksToColor)
  if blockAndLeafBlockCounter[0] >= tooManyBlocksToColor:
    if inBlockCountLimitOverride:
      if PhactoriDbg():
        myDebugPrint3('  warning:  color by blocks is being used with large number of blocks\n')
    else:
      if PhactoriDbg():
        myDebugPrint3('  too many blocks (>10000) for coloring: color data could cause issues\n')
      return

  blockAndLeafBlockCounter = [0L,0L]

  blockColorData = inParaViewDataRepresentation.BlockColor
  if PhactoriDbg(100):
    myDebugPrint3('  blockColorData is ' + str(blockColorData) + '\n', 100)
  if PhactoriDbg(100):
    myDebugPrint3('  blockColorData type is ' + str(type(blockColorData)) + '\n', 100)

  ColorByBlockRecurse1(csdata, blockAndLeafBlockCounter, inParaViewDataRepresentation.BlockColor, False)

  #inParaViewDataRepresentation.BlockColor = blockColorData

  if PhactoriDbg(100):
    myDebugPrint3('   block color data: ' + str(inParaViewDataRepresentation.BlockColor) + '\n', 100)
  if PhactoriDbg(100):
    myDebugPrint3('   final count: ' + str(blockAndLeafBlockCounter) + '\n', 100)
  if PhactoriDbg(100):
    myDebugPrint3('ColorByBlock returning\n', 100)


def UseDataRangeForColorValues(inPvDataRepresentation, inRepresentation,
        inOperation):

  if inRepresentation.mColorVariableInfo.mVariableIntendedForUseFlag == False:
    #the color by variable is not set, we do nothing
    return

  if PhactoriDbg(100):
    myDebugPrint3("UseDataRangeForColorValues entered\n", 100)

  rep = inPvDataRepresentation

  localColorArrayName = inRepresentation.mColorVariableInfo.mVariableName
  if inRepresentation.mColorVariableInfo.mVariableType == 'node':
    pointsOrCellsType = gPointsString
  else:
    pointsOrCellsType = gCellsString

  #if gParaViewCatalystVersionFlag <= 40100:
  #  localColorArrayName = rep.ColorArrayName
  #  pointsOrCellsType = rep.ColorAttributeType
  #else:
  #  localColorArrayName = rep.ColorArrayName[1]
  #  pointsOrCellsType = rep.ColorArrayName[0]

  if PhactoriDbg():
    myDebugPrint3('  color variable is ' + localColorArrayName + '\n')
    myDebugPrint3('  color attribute type is ' + \
        pointsOrCellsType + '\n')

  if PhactoriDbg():
    myDebugPrint3("  rep:" + str(rep) + "\n")

  #if rep.ColorArrayName == '':
  #  return

  input = rep.Input

  if PhactoriDbg():
    myDebugPrint3("  input:" + str(input) + "\n")
  UpdatePipelineWithCurrentTimeArgument(input)

  #assume point info
  #datainformation = rep.Input.GetPointDataInformation()

  if pointsOrCellsType == gPointsString:
    datainformation = input.GetPointDataInformation()
  elif pointsOrCellsType == gCellsString:
    datainformation = input.GetCellDataInformation()
  else:
    if PhactoriDbg():
      myDebugPrint3("something strange with color attribute type :" + \
        str(pointsOrCellsType) + "\n")
    return;

  datarange = [sys.float_info.max, -sys.float_info.max]
  if gParaViewCatalystVersionFlag <= 40100:
    lut = rep.LookupTable
  else:
    lut = GetColorTransferFunction(rep.ColorArrayName[1])

  if lut == None:
    #at this point in time, all processes should have a look up table, even
    #if they don't have associated data
    myDebugPrint3AndException(\
        "UseDataRangeForColorValues expecting LookupTable\n")

  if inRepresentation.mColorVariableInfo.mVariableIsVectorComponent == True:
    lut.VectorMode = 'Component'
    index = inRepresentation.mColorVariableInfo.mVariableComponent
    lut.VectorComponent = index
    dataarray = datainformation.GetArray(localColorArrayName)
    if dataarray == None:
      #no data array, nothing to be done---screws up mpi?
      #return
      datarange = [sys.float_info.max, -sys.float_info.max]
    else:
      datarange = datainformation.GetArray(localColorArrayName).\
          GetRange(index)
  elif inRepresentation.mColorVariableInfo.\
      mVariableIsVectorMagnitude == True:
    lut.VectorMode = 'Magnitude'
    lut.VectorComponent = 0
    dataarray = datainformation.GetArray(localColorArrayName)
    if dataarray == None:
      #no data array, nothing to be done---screws up mpi?
      #return
      datarange = [sys.float_info.max, -sys.float_info.max]
    else:
      numComponents = dataarray.GetNumberOfComponents()
      if numComponents != 3:
        myDebugPrint3AndException(\
            "UseDataRangeForColorValues expecting 3 component item here\n")
      localMinMaxInfo = FindThisProcessorMinMaxForVar(inOperation,
          inRepresentation.mColorVariableInfo)
      if localMinMaxInfo[2] == False:
        datarange = [sys.float_info.max, -sys.float_info.max]
      else:
        datarange = [localMinMaxInfo[0], localMinMaxInfo[1]]
  else:
    #assume one component, scalar, assume one vector component and it's 0
    dataarray = datainformation.GetArray(localColorArrayName)
    if dataarray == None:
      #no data array, nothing to be done---screws up mpi?
      #return
      if PhactoriDbg():
        myDebugPrint3(\
            "  no data array, using range "
            "[sys.float_info.max, -sys.float_info.max]\n")
      datarange = [sys.float_info.max, -sys.float_info.max]
    else:
      datarange = datainformation.GetArray(localColorArrayName).GetRange(0)
      if PhactoriDbg():
        myDebugPrint3('  datarange: ' + str(datarange) + '\n')

  import vtkParallelCorePython
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()

  localarray = vtk.vtkDoubleArray()
  localarray.SetNumberOfTuples(2)
  # negate so that MPI_MAX gets min instead of doing a MPI_MIN and MPI_MAX
  localarray.SetValue(0, -datarange[0])
  localarray.SetValue(1,  datarange[1])

  globalarray = vtk.vtkDoubleArray()
  globalarray.SetNumberOfTuples(2)
  globalController.AllReduce(localarray, globalarray, 0)
  globalDatarange = [ \
    -globalarray.GetTuple1(0), globalarray.GetTuple1(1), \
    ]
  if PhactoriDbg():
    myDebugPrint3("local datarange:  " + str(datarange) + \
      "\nglobal datarange: " + str(globalDatarange) + "\n")

  minmaxTrk = inRepresentation.mColorRangeMinMaxTracker
  minmaxTrk.ResetThisCb()

  if globalDatarange[0] >= globalDatarange[1]:
    #min was >= max, so we can't do a good color mapping; take steps
    #this means we really had no data; therefore set min and max to 0.0
    #and let this system handle the color mapping as if all data values
    #were 0.0
    globalDatarange[0] = 0.0
    globalDatarange[1] = 0.0
    if PhactoriDbg():
      myDebugPrint3("no valid data for variable; setting min/max to 0.0\n")

  minmaxTrk.MinMaxTestAndSet(globalDatarange[0])
  minmaxTrk.MinMaxTestAndSet(globalDatarange[1])
  minmaxTrk.UpdateMinMaxToUse()
  #myDebugPrint3("representation color range tracking:\n")
  #myDebugPrint3(minmaxTrk.SelfToStr())

  minToUse = minmaxTrk.GetMinToUse()
  maxToUse = minmaxTrk.GetMaxToUse()

  if PhactoriDbg():
    myDebugPrint3("minmaxTrk min: " + str(minToUse) + "   max: " + str(maxToUse) + "\n")

  if maxToUse - minToUse > 0.000000000000000000001:

    #if required, restrict color range to subset of otherwise calculated range
    if inRepresentation.mUseColorSubrange:
        clrDelta = maxToUse - minToUse
        maxToUse = minToUse + inRepresentation.mColorSubrange[1] * clrDelta
        minToUse = minToUse + inRepresentation.mColorSubrange[0] * clrDelta

    SetMinimumMaximumColorValues(inPvDataRepresentation,
        minToUse, maxToUse,
        inOperation, inRepresentation,
        localColorArrayName, pointsOrCellsType)
  else:
    SetMinimumMaximumColorValues(inPvDataRepresentation,
        minToUse - 0.1, maxToUse + 0.1,
        inOperation, inRepresentation, localColorArrayName,
        pointsOrCellsType)

  if PhactoriDbg(100):
    myDebugPrint3("UseDataRangeForColorValues returning\n", 100)

#def UpdateDataRangesForColorValues():
#  myDebugPrint3("UpdateDataRangesForColorValues entered\n", 100)
#  global PhactoriRenderViewInfoList
#  for ii in PhactoriRenderViewInfoList:
#    UseDataRangeForColorValues(ii)
#  myDebugPrint3("UpdateDataRangesForColorValues returning\n", 100)

def helpFindNewRgbPoint(inBaseRgbPts, inRatio):
  #print "inRatio: ", inRatio
  positionIndex = -2

  if inBaseRgbPts[0] >= inRatio:
    positionIndex = 0
  else:
    positionIndex = -1
    for rIndex in range(4, len(inBaseRgbPts), 4):
      #print "rIndex: ", rIndex
      #print "positionIndex: ", positionIndex
      #print "val: ", inBaseRgbPts[rIndex]
      if inRatio <= inBaseRgbPts[rIndex]:
          positionIndex = rIndex - 4
          #print "stopping"
          break
      #bad case, ratio is off end
      if positionIndex ==  -1:
          #print "bad case"
          positionIndex = len(inBaseRgbPts) - 4

  if PhactoriDbg(100):
    myDebugPrint3("position index: " + str(positionIndex) + "\n", 100)

  if positionIndex >= len(inBaseRgbPts):
    myDebugPrint3AndException("helpFindNewRgbPoint: ran off end\n")

  zz1 = inBaseRgbPts[positionIndex + 0]
  rr1 = inBaseRgbPts[positionIndex + 1]
  gg1 = inBaseRgbPts[positionIndex + 2]
  bb1 = inBaseRgbPts[positionIndex + 3]
  zz2 = inBaseRgbPts[positionIndex + 4]
  rr2 = inBaseRgbPts[positionIndex + 5]
  gg2 = inBaseRgbPts[positionIndex + 6]
  bb2 = inBaseRgbPts[positionIndex + 7]
  deltaPt = zz2 - zz1
  ratioPt = (inRatio - zz1)/deltaPt
  #print "before Pt: ", zz1, rr1, gg1, bb1
  #print "after Pt : ", zz2, rr2, gg2, bb2
  #print "inRatio: ", inRatio
  #print "deltaPt: ", deltaPt
  #print "ratioPt: ", ratioPt

  subrangePt1 = [inRatio,
           rr1 + ratioPt * (rr2 - rr1),
           gg1 + ratioPt * (gg2 - gg1),
           bb1 + ratioPt * (bb2 - bb1)]
  if PhactoriDbg(100):
      myDebugPrint3("subrangePt1: " + str(subrangePt1) + "\n")
  return subrangePt1


def CalculateColorMapRGBPointsWithSubranges(inBaseRgbPoints,
        inMinimum, inMaximum, inSubrangeMin, inSubrangeMax, inSubrangeColor):
  """creates the rgb values for color map with highlighed subranges"""

  if PhactoriDbg(100):
    myDebugPrint3("CalculateColorMapRGBPointsWithSubranges entered \n"
      "         min, max: " + str(inMinimum) + ", " + str(inMaximum) + "\n"
      "subrange min, max: " + str(inSubrangeMin) + ", " + \
              str(inSubrangeMax) + "\n"
      "subrange color: " + str(inSubrangeColor) + "\n", 100)

  #find where subrange fits in overallrange
  allRange = inMaximum - inMinimum
  minHlRatio = (inSubrangeMin - inMinimum)/allRange
  maxHlRatio = (inSubrangeMax - inMinimum)/allRange

  #test for case where whole subrange is above or below overall range
  if minHlRatio >= 1.0 or maxHlRatio <= 0.0:
      if PhactoriDbg(100):
        myDebugPrint3("subrange is above or below overall range, returning\n"
          "CalculateColorMapRGBPointsWithSubranges returning \n", 100)
      return inBaseRgbPoints

  if minHlRatio < 0.0:
    minHlRatio = 0.0
  if maxHlRatio > 1.0:
    maxHlRatio = 1.0
  if maxHlRatio - minHlRatio < 0.004:
    #very narrow subrange; use narrow bands but might color odd
    #issue warning?
    ratDelta = maxHlRatio - minHlRatio
    bandStep = ratDelta * 0.25
  else:
    bandStep = 0.001

  #construct the no-subrange color map, with extra items to track the
  #outer range colors while changing the inner subrange colors
  xrgbPoints = inBaseRgbPoints
  #xrgbPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471,
  #              0.5, 0.865,         0.865,          0.865,
  #              1.0, 0.705882352941, 0.0156862745098, 0.149019607843]
  if PhactoriDbg(100):
    myDebugPrint3("inBaseRgbPoints:\n" + str(inBaseRgbPoints) + "\n", 100)

  #find the 4 points for the subrange
  subrangePt1 = helpFindNewRgbPoint(xrgbPoints, minHlRatio)
  subrangePt2 = [minHlRatio + bandStep,
          inSubrangeColor[0], inSubrangeColor[1], inSubrangeColor[2]]
  subrangePt3 = [maxHlRatio - bandStep,
          inSubrangeColor[0], inSubrangeColor[1], inSubrangeColor[2]]
  subrangePt4 = helpFindNewRgbPoint(xrgbPoints, maxHlRatio)

  #construct the map with the subrange
  xrgbIndex = 0
  myRgbPoints = []

  #first do points from original color map which are below subrange
  while (xrgbIndex < len(xrgbPoints)) and xrgbPoints[xrgbIndex] < subrangePt1[0]:
    myRgbPoints.append(xrgbPoints[xrgbIndex+0])
    myRgbPoints.append(xrgbPoints[xrgbIndex+1])
    myRgbPoints.append(xrgbPoints[xrgbIndex+2])
    myRgbPoints.append(xrgbPoints[xrgbIndex+3])
    xrgbIndex += 4

  #add in subrange points, accounting for bottom/top cases
  if subrangePt1[0] > xrgbPoints[0]:
    for nn in subrangePt1:
      myRgbPoints.append(nn)
    for nn in subrangePt2:
      myRgbPoints.append(nn)
  else:
    #special case, bottom of subrange is bottom of whole range
    subrangePt2[0] = subrangePt1[0]
    for nn in subrangePt2:
      myRgbPoints.append(nn)

  numXrgbPoints = int(len(xrgbPoints) / 4)
  if subrangePt4[0] < xrgbPoints[(numXrgbPoints-1) * 4]:
    for nn in subrangePt3:
      myRgbPoints.append(nn)
    for nn in subrangePt4:
      myRgbPoints.append(nn)
  else:
    #special case, top of subrange is top of whole range
    subrangePt3[0] = subrangePt4[0]
    for nn in subrangePt3:
      myRgbPoints.append(nn)

  #add in points from original color map after subrange, skipping any points
  #which were subsumed in subrange
  while subrangePt4[0] >= xrgbPoints[xrgbIndex]:
    xrgbIndex += 4
    if xrgbIndex >= len(xrgbPoints):
        break

  while xrgbIndex < len(xrgbPoints):
    myRgbPoints.append(xrgbPoints[xrgbIndex+0])
    myRgbPoints.append(xrgbPoints[xrgbIndex+1])
    myRgbPoints.append(xrgbPoints[xrgbIndex+2])
    myRgbPoints.append(xrgbPoints[xrgbIndex+3])
    xrgbIndex += 4

  if PhactoriDbg(100):
    myDebugPrint3('subrange: [' + \
      str(inSubrangeMin) + ', ' + \
      str(inSubrangeMax) + ']\n' + \
      'subrange ratio: [' + \
      str(minHlRatio) + ', ' + \
      str(maxHlRatio) + ']\n' + \
      'myRGBPoints: \n' + str(myRgbPoints) + '\n')

  if PhactoriDbg(100):
    myDebugPrint3("CalculateColorMapRGBPointsWithSubranges returning \n", 100)

  return myRgbPoints


global gDefault_LUT_Num_RGBPoints
global gDefault_LUT_RGBPoints
global gDefault_PWF_Num_Points
global gDefault_PWF_Points
gDefault_LUT_Num_RGBPoints = -1
gDefault_LUT_RGBPoints = []
gDefault_PWF_Num_Points = -1
gDefault_PWF_Points = []

def SetMinimumMaximumColorValues(inPvDataRepresentation,
        inMinimum, inMaximum,
        inPhactoriOperation, inPhactoriRepresentation,
        inVariableName, inPointsOrCells):
  "This sets the data values which will correspond to the colors at the " \
  "bottom and top of the color table.  Data values above or below this " \
  "range will be clamped to the color table bottom or top."
  inColorMapSettings = inPhactoriRepresentation.mColorMapSettings
  if PhactoriDbg(100):
    myDebugPrint3('phactori.SetMinimumMaximumColorValues entered, min: ' + \
      str(inMinimum) + ' max: ' + str(inMaximum) + '\n' + \
      'inVariableName: ' + inVariableName + \
      '   inPointsOrCells: ' + inPointsOrCells + '\n', 100)

  if inPhactoriRepresentation.mColorVariableInfo.mVariableTypeNeedsDetection:
    if PhactoriDbg(100):
      myDebugPrint3("  variable type was undetected, returning (1)\n")
    return
  if gParaViewCatalystVersionFlag >= 40300:
    savedActiveSourceZ = GetActiveSource()
    SetActiveSource(inPhactoriOperation.GetPvFilter())
    if PhactoriDbg():
      myDebugPrint3("using ColorBy() to set color var " + \
        inVariableName + "  " + inPointsOrCells + " 4\n")
    ColorBy(inPvDataRepresentation, (inPointsOrCells, inVariableName))

    if inPhactoriRepresentation.mNameOfPresetToUse != None:
      if PhactoriDbg():
          myDebugPrint3("using preset with name: " + \
                  inPhactoriRepresentation.mNameOfPresetToUse + "\n")
      pv_4_3_LUT = GetColorTransferFunction(inVariableName)
      pv_4_3_PWF = GetOpacityTransferFunction(inVariableName)
      #second argument is rescale to data range; need to handle correctly
      pv_4_3_LUT.ApplyPreset(inPhactoriRepresentation.mNameOfPresetToUse, False)
      pv_4_3_PWF.ApplyPreset(inPhactoriRepresentation.mNameOfPresetToUse, False)
      #pv_4_3_LUT.ApplyPreset(inPhactoriRepresentation.mNameOfPresetToUse, True)
      #pv_4_3_PWF.ApplyPreset(inPhactoriRepresentation.mNameOfPresetToUse, True)
    else:
      pv_4_3_LUT = GetColorTransferFunction(inVariableName)
      if pv_4_3_LUT.Discretize != 1:
        pv_4_3_LUT.Discretize = 1
      if pv_4_3_LUT.NumberOfTableValues != 1024:
        pv_4_3_LUT.NumberOfTableValues = 1024

      pv_4_3_PWF = GetOpacityTransferFunction(inVariableName)

      global gDefault_LUT_Num_RGBPoints
      global gDefault_LUT_RGBPoints
      global gDefault_PWF_Num_Points
      global gDefault_PWF_Points
      if gDefault_LUT_Num_RGBPoints < 0:
          for ii in pv_4_3_LUT.RGBPoints:
              gDefault_LUT_RGBPoints.append(ii)
          gDefault_LUT_Num_RGBPoints = len(gDefault_LUT_RGBPoints)
          for ii in pv_4_3_PWF.Points:
              gDefault_PWF_Points.append(ii)
          gDefault_PWF_Num_Points = len(gDefault_PWF_Points)

      if len(pv_4_3_LUT.RGBPoints) != gDefault_LUT_Num_RGBPoints:
          #revert back to default rgb points starting setup
          pv_4_3_LUT.RGBPoints = gDefault_LUT_RGBPoints

      if len(pv_4_3_PWF.Points) != gDefault_PWF_Num_Points:
          #revert back to default points starting setup
          pv_4_3_PWF.Points = gDefault_PWF_Points

      pv_4_3_LUT.ScalarRangeInitialized = 1.0
      pv_4_3_PWF.ScalarRangeInitialized = 1

      pv_4_3_LUT.RescaleTransferFunction(inMinimum, inMaximum)
      pv_4_3_PWF.RescaleTransferFunction(inMinimum, inMaximum)
      #pv_4_3_LUT.LockScalarRange = 1
      if PhactoriDbg(100):
          myDebugPrint3("rgb points: " + str(pv_4_3_LUT.RGBPoints) + "\n" + \
                  "opacity points: " + str(pv_4_3_PWF.Points) + "\n")
      if inPhactoriRepresentation.mUseHighlightSubranges:
        baseRgbPoints, myColorSpace, myNanColor = \
          GetColorMapInfoFromColorLegendCollection(inColorMapSettings,
            0.0, 1.0)
        if PhactoriDbg(100):
          myDebugPrint3(
            "default base rgb points: " + str( baseRgbPoints) + "\n", 100)
        #add each subrange into the base color map one at a time
        for oneSubrange in inPhactoriRepresentation.mHighlightSubranges:
          myRGBPoints = CalculateColorMapRGBPointsWithSubranges(baseRgbPoints,
            inMinimum, inMaximum, oneSubrange[0], oneSubrange[1], oneSubrange[2])
          baseRgbPoints = myRGBPoints

        #adjust the color range from 0 to 1 to the overall min to max
        for ii in range(0, len(myRGBPoints), 4):
          allRange = inMaximum - inMinimum
          myRGBPoints[ii] = inMinimum + myRGBPoints[ii] * allRange
      else:
        myRGBPoints, myColorSpace, myNanColor = \
          GetColorMapInfoFromColorLegendCollection(inColorMapSettings,
            inMinimum, inMaximum)
      pv_4_3_LUT.RGBPoints = myRGBPoints
      pv_4_3_LUT.ColorSpace = myColorSpace
      pv_4_3_LUT.NanColor = myNanColor
      if PhactoriDbg(100):
        myDebugPrint3('phactori.SetMinimumMaximumColorValues returning\n')

    SetActiveSource(savedActiveSourceZ)

    return

  savedActiveSourceZ = GetActiveSource()
  SetActiveSource(inPhactoriOperation.GetPvFilter())

  new_var_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 1.0, 1.0] )

  new_var_PVLookupTable = inPvDataRepresentation.LookupTable
  myRGBPoints, myColorSpace, myNanColor = \
      GetColorMapInfoFromColorLegendCollection(inColorMapSettings,
          inMinimum, inMaximum)
  new_var_PVLookupTable.RGBPoints = myRGBPoints
  new_var_PVLookupTable.ColorSpace = myColorSpace
  if PhactoriDbg():
    myDebugPrint3("nan color: " + str(myNanColor) + "  len: " + str(len(myNanColor)) + "\n")
  new_var_PVLookupTable.NanColor = myNanColor

  if PhactoriDbg():
    myDebugPrint3("new lookup table (2):\n");
    myDebugPrint3(str(new_var_PVLookupTable) + "\n")
    myDebugPrint3("RGBPoints:\n");
    myDebugPrint3(str(new_var_PVLookupTable.RGBPoints) + "\n")

  try:
    inPvDataRepresentation.ScalarOpacityFunction = new_var_PiecewiseFunction
  except:
    if PhactoriDbg():
      myDebugPrint3('no inPvDataRepresentation.ScalarOpacityFunction\n');

  SetActiveSource(savedActiveSourceZ)

  myDebugPrint3('phactori.SetMinimumMaximumColorValues returning (2)\n')

#'Closest Triad'
#'Furthest Triad'
#'Static Triad'
#'Static Edges'
#'Outer Edges'
#def SetCubeAxesMode(inNewMode):
#  global currentPhactoriRenderViewInfo
#  currentPhactoriRenderViewInfo.DataRepresentation1.CubeAxesFlyMode = inNewMode

def ShowCubeAxesXX(inPvRenderView, inOnOrOff, inShowDataCubeAxesInfo = None):
  if PhactoriDbg(100):
    myDebugPrint3("ShowCubeAxesXX entered\n", 100)
  if(inOnOrOff == 'on'):
    if PhactoriDbg():
      myDebugPrint3("  ShowCubeAxesXX: turning on\n")
    #inPvDataRepresentation.CubeAxesVisibility = 1
    SetAxesGridVisibility(inPvRenderView, 1)
    if inShowDataCubeAxesInfo != None:
      inPvRenderView.AxesGrid.ShowEdges = inShowDataCubeAxesInfo.mShowEdges
      inPvRenderView.AxesGrid.ShowTicks = inShowDataCubeAxesInfo.mShowTicks
      axisXInfo = inShowDataCubeAxesInfo.mXAxisInfo
      axisYInfo = inShowDataCubeAxesInfo.mYAxisInfo
      axisZInfo = inShowDataCubeAxesInfo.mZAxisInfo

      ##axis labels as specified
      if axisXInfo.mUseLabelFlag:
        inPvRenderView.AxesGrid.XTitle = axisXInfo.mAxisLabel
        if PhactoriDbg():
          myDebugPrint3("  x axis label set to: " + inPvRenderView.AxesGrid.XTitle + "\n")
      else:
        inPvRenderView.AxesGrid.XTitle = "X Axis"
      if axisYInfo.mUseLabelFlag:
        inPvRenderView.AxesGrid.YTitle = axisYInfo.mAxisLabel
        if PhactoriDbg():
          myDebugPrint3("  y axis label set to: " + inPvRenderView.AxesGrid.YTitle + "\n")
      else:
        inPvRenderView.AxesGrid.YTitle = "Y Axis"
      if axisZInfo.mUseLabelFlag:
        inPvRenderView.AxesGrid.ZTitle = axisZInfo.mAxisLabel
        if PhactoriDbg():
          myDebugPrint3("  z axis label set to: " + inPvRenderView.AxesGrid.ZTitle + "\n")
      else:
        inPvRenderView.AxesGrid.ZTitle = "Z Axis"
  else:
    if PhactoriDbg():
      myDebugPrint3("  ShowCubeAxesXX: turning off\n")
    #inPvDataRepresentation.CubeAxesVisibility = 0
    SetAxesGridVisibility(inPvRenderView, 0)
  if PhactoriDbg(100):
    myDebugPrint3("ShowCubeAxesXX returning\n", 100)


def ShowDataColorLegendXX(inPvView,
        inOnOffSetting, inColorLegendPositionAndSize, inColorSettings,
        inColorLegendRepRef, inPvDataRep):
  """Turns on or off the display of the color bar legend showing the mapping
     between the color and the data value (and the name of the data value.
     (note this is primarily to do the paraview-side work to turn bar on or
      off, on to set up for rendering in the shared view, off to turn off
      rendering in shared view.  On or off state for rendering is stored
      as a flag in the ioPhactoriImagesetBlock instance"""

  if PhactoriDbg(100):
    myDebugPrint3('phactori.ShowDataColorLegendXX entered, setting:' + \
        inOnOffSetting + '\n', 100)
  if(inOnOffSetting == 'on'):
    myVisibility = 1
  else:
    myVisibility = 0
    if inColorLegendRepRef != None:
      if PhactoriDbg(100):
        myDebugPrint3("A inColorLegendRepRef was " + \
          str(inColorLegendRepRef.Visibility) + \
          " now 0: " + str(inColorLegendRepRef) + "\n")
      inColorLegendRepRef.Visibility = 0
    myDebugPrint3(
        'phactori.ShowDataColorLegendXX returning with none rep: ' + \
        inOnOffSetting + '\n', 100)
    return None

  if gParaViewCatalystVersionFlag <= 40100:
    localColorArrayName = \
        inPvDataRep.ColorArrayName
  else:
    localColorArrayName = \
        inPvDataRep.ColorArrayName[1]

  if inColorLegendRepRef != None:
      #myDebugPrint3('  phactori.ShowDataColorLegend using rep reference\n', 100)
      #if inColorLegendRepRef.ColorArrayName == '':
      #  return
      inPvView.OrientationAxesLabelColor = inColorSettings.mTextColor
      if PhactoriDbg(100):
        myDebugPrint3("B inColorLegendRepRef was " + \
          str(inColorLegendRepRef.Visibility) + \
          " now " + str(myVisibility) + ": " + str(inColorLegendRepRef) + "\n")
      inColorLegendRepRef.Visibility = myVisibility
      inColorLegendRepRef.LookupTable = inPvDataRep.LookupTable
      inColorLegendRepRef.Title = localColorArrayName
      #ioPhactoriImagesetBlock.mColorLegendRepRef.Color = \
      #    inColorSettings.mTextColor
      if PhactoriDbg(100):
        myDebugPrint3(
            'phactori.ShowDataColorLegendXX returning with old rep: ' + \
            inOnOffSetting + '\n', 100)
      return inColorLegendRepRef
  #else:
      #myDebugPrint3('  phactori.ShowDataColorLegend have to create rep reference\n', 100)

  if inColorLegendPositionAndSize[0] == 'parameters':
    legendSizeMultiplier = None
    legendSize = inColorLegendPositionAndSize[3]
    legendFontSize = inColorLegendPositionAndSize[4]
  else:
    legendSizeMultiplier = inColorLegendPositionAndSize[1]
    if gParaViewCatalystVersionFlag <= 40100:
      #our own factor to make legends smaller generally
      legendSizeMultiplier *= 0.6
    else:
      #legendSizeMultiplier *= 1.0
      legendSizeMultiplier *= 0.7
    #legendFontSize = int(12.0 * legendSizeMultiplier)
    legendFontSize = int(9.0 * legendSizeMultiplier)

  if gParaViewCatalystVersionFlag <= 40100:
    colorLegendDefaultLongSize = 0.5
    colorLegendDefaultShortSize = 0.13
  else:
    #[0.12, 0.43]
    #[0.85, 0.05]
    colorLegendDefaultLongSize = 0.43
    colorLegendDefaultShortSize = 0.12

  colorLegendAdjustedLongSize = colorLegendDefaultLongSize * legendSizeMultiplier
  colorLegendAdjustedShortSize = colorLegendDefaultShortSize * legendSizeMultiplier
  horizontalLegendSize = [colorLegendAdjustedLongSize, colorLegendAdjustedShortSize]
  verticalLegendSize = [colorLegendAdjustedShortSize, colorLegendAdjustedLongSize]

  xPosForBottomTop = 0.5 - 0.5 * colorLegendAdjustedLongSize
  yPosForLeftRight = 0.5 - 0.5 * colorLegendAdjustedShortSize

  if gParaViewCatalystVersionFlag < 50400:
    if inColorLegendPositionAndSize[0] == 'top':
      legendOrientation = 'Horizontal'
      legendSize = horizontalLegendSize
      legendPosition=[xPosForBottomTop, 0.85]
    elif inColorLegendPositionAndSize[0] == 'bottom':
      legendOrientation = 'Horizontal'
      legendSize = horizontalLegendSize
      legendPosition=[xPosForBottomTop, 0.02]
    elif inColorLegendPositionAndSize[0] == 'left':
      legendOrientation = 'Vertical'
      legendSize = verticalLegendSize
      legendPosition=[0.065, yPosForLeftRight]
    elif inColorLegendPositionAndSize[0] == 'right':
      legendOrientation = 'Vertical'
      legendSize = verticalLegendSize
      legendPosition=[0.9, yPosForLeftRight]
    elif inColorLegendPositionAndSize[0] == 'top left':
      legendOrientation = 'Horizontal'
      legendSize = horizontalLegendSize
      legendPosition=[0.065, 0.85]
    elif inColorLegendPositionAndSize[0] == 'top right':
      legendOrientation = 'Horizontal'
      legendSize = horizontalLegendSize
      legendPosition=[0.7, 0.85]
    elif inColorLegendPositionAndSize[0] == 'bottom left':
      legendOrientation = 'Horizontal'
      legendSize = horizontalLegendSize
      #legendPosition=[0.065, 0.85]
      legendPosition=[0.065, 0.01]
    elif inColorLegendPositionAndSize[0] == 'bottom right':
      legendOrientation = 'Horizontal'
      legendSize = horizontalLegendSize
      #legendPosition=[0.7, 0.05]
      legendPosition=[0.7, 0.01]
    elif inColorLegendPositionAndSize[0] == 'parameters':
      legendOrientation = inColorLegendPositionAndSize[1]
      legendSize = horizontalLegendSize
      legendPosition = inColorLegendPositionAndSize[2]
    else:
      legendOrientation = 'Horizontal'
      legendSize = horizontalLegendSize
      #legendPosition=[xPosForBottomTop, 0.05]
      legendPosition=[xPosForBottomTop, 0.01]

    if PhactoriDbg():
      myDebugPrint3("legend info:\n  legendSizeMultiplier: " + str(legendSizeMultiplier) + "\n" \
        "  legendSize: " + str(legendSize) + "\n" \
        "  legendPos: " + str(legendPosition) + "\n"\
        "  legendOrientation: " + str(legendOrientation) + "\n"\
        "  legendFontSize: " + str(legendFontSize) + "\n")

  else:
    defaultLegendLength = 0.33
    defaultMidPos = 0.5 - 0.5*defaultLegendLength
    #legendFontSize = 16
    #legendSize = 1.0
    #validPositions = ['UpperLeftCorner', 'UpperRightCorner',
    #    'LowerLeftCorner', 'LowerRightCorner',
    #    'UpperCenter', 'LowerCenter']
    legendPosition=[0.0, 0.0]
    if inColorLegendPositionAndSize[0] == 'top':
      legendOrientation = 'Horizontal'
      legendWindowLocation = 'UpperCenter'
    elif inColorLegendPositionAndSize[0] == 'bottom':
      legendOrientation = 'Horizontal'
      legendWindowLocation = 'LowerCenter'
    elif inColorLegendPositionAndSize[0] == 'left':
      legendOrientation = 'Vertical'
      legendPosition=[0.02, defaultMidPos]
      legendWindowLocation = 'AnyLocation'
    elif inColorLegendPositionAndSize[0] == 'right':
      legendOrientation = 'Vertical'
      legendPosition=[0.89, defaultMidPos]
      legendWindowLocation = 'AnyLocation'
    elif inColorLegendPositionAndSize[0] == 'top left':
      legendOrientation = 'Vertical'
      legendWindowLocation = 'UpperLeftCorner'
    elif inColorLegendPositionAndSize[0] == 'top right':
      legendOrientation = 'Vertical'
      legendWindowLocation = 'UpperRightCorner'
    elif inColorLegendPositionAndSize[0] == 'bottom left':
      legendOrientation = 'Vertical'
      legendWindowLocation = 'LowerLeftCorner'
    elif inColorLegendPositionAndSize[0] == 'bottom right':
      legendOrientation = 'Vertical'
      legendWindowLocation = 'LowerRightCorner'
    elif inColorLegendPositionAndSize[0] == 'parameters':
      legendOrientation = inColorLegendPositionAndSize[1]
      legendPosition = inColorLegendPositionAndSize[2]
      legendWindowLocation = 'AnyLocation'
    else:
      legendOrientation = 'Vertical'
      legendWindowLocation = 'LowerRightCorner'

    #newScalarBarWidgetRepresentation = CreateScalarBar( Title=inPvDataRep.ColorArrayName, Position2=[0.13, 0.5], TitleOpacity=1.0, TitleShadow=0, AutomaticLabelFormat=1, TitleFontSize=12, TitleColor=[1.0, 1.0, 1.0], AspectRatio=20.0, NumberOfLabels=5, ComponentTitle='', Resizable=1, TitleFontFamily='Arial', Visibility=myVisibility, LabelFontSize=12, LabelFontFamily='Arial', TitleItalic=0, Selectable=0, LabelItalic=0, Enabled=0, LabelColor=[1.0, 1.0, 1.0], Position=[0.9, 0.31396255850234012], LabelBold=0, UseNonCompositedRenderer=1, LabelOpacity=1.0, TitleBold=0, LabelFormat='%-#6.3g', Orientation='Vertical', LabelShadow=0, LookupTable=inPvDataRep.LookupTable, Repositionable=1 )
  if gParaViewCatalystVersionFlag <= 40100:
    newScalarBarWidgetRepresentation = CreateScalarBar(Title=localColorArrayName,
      Orientation=legendOrientation,
      Position=legendPosition,
      Position2 = legendSize,
      Visibility=myVisibility,
      LookupTable=inPvDataRep.LookupTable,
      LabelFontSize=legendFontSize,
      TitleOpacity=1.0,
      TitleShadow=0,
      AutomaticLabelFormat=1,
      TitleFontSize=legendFontSize,
      TitleColor=inColorSettings.mTextColor,
      AspectRatio=20.0,
      NumberOfLabels=5,
      ComponentTitle='',
      Resizable=1,
      TitleFontFamily='Arial',
      LabelFontFamily='Arial',
      TitleItalic=0,
      Selectable=0,
      LabelItalic=0,
      Enabled=0,
      LabelColor=inColorSettings.mTextColor,
      LabelBold=0,
      UseNonCompositedRenderer=1,
      LabelOpacity=1.0,
      TitleBold=0,
      LabelFormat='%-#6.3g',
      LabelShadow=0,
      Repositionable=1)
  elif gParaViewCatalystVersionFlag < 50400:
    newScalarBarWidgetRepresentation = CreateScalarBar(Title=localColorArrayName,
      Orientation=legendOrientation,
      Position=legendPosition,
      Position2 = legendSize,
      Visibility=myVisibility,
      LookupTable=inPvDataRep.LookupTable,
      LabelFontSize=legendFontSize,
      #TitleOpacity=1.0,
      #TitleShadow=0,
      #AutomaticLabelFormat=1,
      TitleFontSize=legendFontSize,
      TitleColor=inColorSettings.mTextColor,
      AspectRatio=20.0,
      #NumberOfLabels=5,
      ComponentTitle='',
      #Resizable=1,
      #TitleFontFamily='Arial',
      #LabelFontFamily='Arial',
      #TitleItalic=0,
      #Selectable=0,
      #LabelItalic=0,
      #Enabled=0,
      LabelColor=inColorSettings.mTextColor,
      #LabelBold=0,
      #UseNonCompositedRenderer=1,
      #LabelOpacity=1.0,
      #TitleBold=0,
      #LabelFormat='%-#6.3g',
      #LabelShadow=0,
      #Repositionable=1
      )
  else:
    newScalarBarWidgetRepresentation = CreateScalarBar(
        Title=localColorArrayName, ComponentTitle='')
    newScalarBarWidgetRepresentation.Orientation = legendOrientation
    newScalarBarWidgetRepresentation.WindowLocation = legendWindowLocation
    if legendWindowLocation == 'AnyLocation':
      newScalarBarWidgetRepresentation.Position = legendPosition
    if PhactoriDbg():
      nbwr = newScalarBarWidgetRepresentation
      myDebugPrint3("newScalarBarWidgetRepresentation:\n" +\
        str(nbwr) + "\n" +\
        "  Title: " + str(nbwr.Title) + "\n" +\
        "  ComponentTitle: " + str(nbwr.ComponentTitle) + "\n" +\
        "  WindowLocation: " + str(nbwr.WindowLocation) + "\n" +\
        #"  LockPosition: " + str(nbwr.LockPosition) + "\n" +\
        #"  Repositionable: " + str(nbwr.Repositionable) + "\n" +\
        #"  AutoOrient: " + str(nbwr.AutoOrient) + "\n" +\
        "  Position: " + str(nbwr.Position) + "\n" +\
        "  ScalarBarLength: " + str(nbwr.ScalarBarLength) + "\n" +\
        "  ScalarBarThickness: " + str(nbwr.ScalarBarThickness) + "\n" +\
        "  Orientation: " + str(nbwr.Orientation) + "\n" +\
        "  LabelFontSize: " + str(nbwr.LabelFontSize) + "\n" +\
        "  TitleFontSize: " + str(nbwr.TitleFontSize) + "\n" +\
        "  LabelFontFamily: " + str(nbwr.LabelFontFamily) + "\n" +\
        "  TitleFontFamily: " + str(nbwr.TitleFontSize) + "\n")

  inPvView.OrientationAxesLabelColor = inColorSettings.mTextColor
  inPvView.Representations.append(newScalarBarWidgetRepresentation)
  #ioPhactoriImagesetBlock.mColorLegendRepRef = \
  #    newScalarBarWidgetRepresentation
  if PhactoriDbg():
    myDebugPrint3("current lookup table:\n");
  if PhactoriDbg():
    myDebugPrint3(str(inPvDataRep.LookupTable) + '\n')
  if PhactoriDbg():
    myDebugPrint3("RGBPoints:\n");
  if PhactoriDbg():
    if gParaViewCatalystVersionFlag <= 40100:
      myLocalRGBPoints = inPvDataRep.LookupTable.RGBPoints
    else:
      pv_4_3_LUT = GetColorTransferFunction(
        inPvDataRep.ColorArrayName[1])
      myLocalRGBPoints = pv_4_3_LUT.RGBPoints
    myDebugPrint3(str(myLocalRGBPoints) + '\n')
  if PhactoriDbg():
    myDebugPrint3("widget:\n");
  if PhactoriDbg():
    myDebugPrint3(str(newScalarBarWidgetRepresentation) + '\n')

  if PhactoriDbg(100):
    myDebugPrint3('phactori.ShowDataColorLegendXX returning with new rep: ' + \
        inOnOffSetting + '\n', 100)

  return newScalarBarWidgetRepresentation

def ApplyClipPlane(inNormal, inOrigin, inClipPlaneFilterName):
  newClip = Clip(ClipType = "Plane")
  newClip.ClipType.Normal = inNormal
  newClip.ClipType.Origin = inOrigin
  SetActiveSource(newClip);
  AddFilterToFilterMap(inClipPlaneFilterName, newClip)

global gThresholdFilterNameCounter
gThresholdFilterNameCounter = 0

def ThresholdFilter(inVariableName, inType, inRange, inThresholdFilterName = None):
  "Apply a threshold filter.  inVariableName is the variable to use for "
  "thresholding, inType is 'POINTS' or 'CELLS', inRange is the threshold "
  "range, such as [0.5, 1.5] or [-10.0, 10.0]"
  if PhactoriDbg(100):
    myDebugPrint3('phactori.ThresholdFilter entered, setting:' + inVariableName + ' ' + inType + ' ' + str(inRange) + '\n', 100)

  activeSource = GetActiveSource()

  CellDataList = []
  for ii in range(activeSource.CellData.GetNumberOfArrays()):
      CellDataList.append(activeSource.CellData.GetArray(ii).Name)
  if PhactoriDbg():
    myDebugPrint3('before threshold cell data items: ' + str(CellDataList) + '\n');

  myThreshold = Threshold(activeSource)
  myThreshold.ThresholdRange = inRange
  myThreshold.Scalars = [inType, inVariableName]
  SetActiveSource(myThreshold)
  if inThresholdFilterName == None:
    global gThresholdFilterNameCounter
    inThresholdFilterName = "ThresholdFilter" + str(gThresholdFilterNameCounter)
    gThresholdFilterNameCounter += 1
  AddFilterToFilterMap(inThresholdFilterName, myThreshold)

  activeSource = GetActiveSource()
  CellDataList = []
  for ii in range(activeSource.CellData.GetNumberOfArrays()):
      CellDataList.append(activeSource.CellData.GetArray(ii).Name)
  if PhactoriDbg():
    myDebugPrint3('after threshold cell data items: ' + str(CellDataList) + '\n');

  if PhactoriDbg(100):
    myDebugPrint3('phactori.ThresholdFilter returning \n', 100)


#def WarpMeshByDisplacement(inDisplacementVariableName, inWarpDisplacementFilterName = 'WarpDisplacementCalculator'):
#  "Warp the mesh by the displacement variable; in paraview this means we "
#  "add a calculator with adds the displacement value to the coordinates "
#  "of each corresponding mesh element"
#
#  print 'WarpMeshByDisplacement entered'
#  print '  active view A: ' + str(GetActiveView())
#  print '  active source A: ' + str(GetActiveSource())
#
#
#
#  myDebugPrint3('phactori.WarpMeshByDisplacement entered, setting:' + inDisplacementVariableName + '\n')
#
#  PointDataList = []
#
#  FunctionString = 'coords + ' + inDisplacementVariableName
#  DisplacementWarpCalculator = Calculator( guiName="Calculator1", Function=FunctionString, ReplacementValue=0.0, ResultArrayName='Result', ReplaceInvalidResults=1, CoordinateResults=1 )
#  #DisplacementWarpCalculator = Calculator( guiName="Calculator1", Function=FunctionString, ReplacementValue=0.0, ResultArrayName='Result', ReplaceInvalidResults=1, AttributeMode='point_data', CoordinateResults=1 )
#  myDebugPrint3('DisplacementWarpCalculator has ' + str(DisplacementWarpCalculator.PointData.GetNumberOfArrays()) + ' point data items\n');
#  for ii in range(DisplacementWarpCalculator.PointData.GetNumberOfArrays()):
#      PointDataList.append(DisplacementWarpCalculator.PointData.GetArray(ii).Name)
#  myDebugPrint3('DisplacementWarpCalculator point data items: ' + str(PointDataList) + '\n');
#  #tempDataRepresentation = GetDisplayProperties(currentSi.ActiveSource)
#  #tempDataRepresentation.Visibility = 0
#  SetActiveSource(DisplacementWarpCalculator);
#  AddFilterToFilterMap(inWarpDisplacementFilterName, DisplacementWarpCalculator)
#
#  print '  active view B: ' + str(GetActiveView())
#  print '  active source B: ' + str(GetActiveSource())
#  print 'WarpMeshByDisplacement returning'

def TempAddTestFiltersAfterDispl():
    """create the clip plane filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3('TempAddTestFiltersAfterDispl entered\n', 100)
    #info in block class should already be parsed and checked

    #savedActiveSource = GetActiveSource()
    #newCellToPointFilter = CellDatatoPointData()

    ##SetActiveSource(newCellToPointFilter)
    ##SetActiveSource(savedActiveSource)

    #SetActiveSource(newCellToPointFilter)
    #newCellToPointFilter.UpdatePipeline()
    #AddFilterToFilterMap('initialCellToPointFilter', newCellToPointFilter)

    #newContourFilter = Contour(PointMergeMethod="Uniform Binning")
    #myDebugPrint3("  contour filter is: " + str(newContourFilter) + "\n")
    #newContourFilter.PointMergeMethod = "Uniform Binning"
    #newContourFilter.Isosurfaces = [9687]
    #myDebugPrint3("  variable is: " + str('VON_MISES') + "\n")
    #newContourFilter.ContourBy = ['POINTS', 'Global Node Id']

    #myDebugPrint3("  still okay\n")
    #SetActiveSource(newContourFilter)
    #AddFilterToFilterMap('initialContourFilter', newContourFilter)

    #Contour2 = Contour( guiName="Contour2", Isosurfaces=[9882.0], ComputeNormals=1, GenerateTriangles=1, ComputeGradients=0, ComputeScalars=0, ContourBy=['POINTS', 'GlobalNodeId'], PointMergeMethod="Uniform Binning" )
    #Contour2.PointMergeMethod.Numberofpointsperbucket = 8
    #Contour2.PointMergeMethod.Divisions = [50, 50, 50]

    if PhactoriDbg(100):
      myDebugPrint3('TempAddTestFiltersAfterDispl returning\n', 100)


#makes a filter which pulls a component out of a vector so it can be treated
#separately
def ExtractComponentFromVectorVariable(inVariableName, inWhichComponent, inResultName = 'Result', inFilterName = None):

  if inFilterName == None:
    inFilterName = "ExtractComponent_" + str(inWhichComponent) + GetSeparatorString() + inVariableName
  if PhactoriDbg(100):
    myDebugPrint3('phactori.ExtractComponentFromVectorVariable entered, variable name:' + \
      inVariableName + ', ' + inResultName + ', ' + \
      str(inWhichComponent) + '\n', 100)

  #theActiveSource = GetActiveSource()
  #if(theActiveSource == None):
  #  myDebugPrint3('phactori.ColorByVectorVariable no active source, returning:\n', 100)
  #  return

  #debugPrintRepresentationsAndVisibility(currentRvi)
  newFunctionString = ""
  if(inWhichComponent == 0):
    newFunctionString = inVariableName + '_X'
  elif(inWhichComponent == 1):
    newFunctionString = inVariableName + '_Y'
  elif(inWhichComponent == 2):
    newFunctionString = inVariableName + '_Z'
  else:
    if PhactoriDbg():
      myDebugPrint3('phactori.ColorByVectorVariable unexpected component (needs 0, 1, or 2)\n')
    exit()

  #DisplacementWarpCalculator = Calculator( guiName="Calculator1", Function=FunctionString, ReplacementValue=0.0, ResultArrayName='Result', ReplaceInvalidResults=1, AttributeMode='point_data', CoordinateResults=1 )

  calc2Source = Calculator( guiName=inFilterName, Function=newFunctionString, ReplacementValue=0.0, ResultArrayName=inResultName, ReplaceInvalidResults=1)

  if PhactoriDbg():
    myDebugPrint3('Calulator2 has ' + str(calc2Source.PointData.GetNumberOfArrays()) + ' point data items\n');
  PointDataList = []
  for ii in range(calc2Source.PointData.GetNumberOfArrays()):
      PointDataList.append(calc2Source.PointData.GetArray(ii).Name)
  if PhactoriDbg():
    myDebugPrint3('calc2Source point data items: ' + str(PointDataList) + '\n');
  #SetActiveSource(currentSi.ActiveSource);

  #tempDataRepresentation = GetDisplayProperties(currentSi.ActiveSource)
  #tempDataRepresentation.Visibility = 0
  SetActiveSource(calc2Source);
  #currentSi.ActiveSource.ResultArrayName = currentSi.ActiveSource.Function
  #currentSi.ActiveSource = calc2Source;
  AddFilterToFilterMap(inFilterName, calc2Source)
  #myDebugPrint3('calculator function: ' + currentSi.ActiveSource.Function + '\n');
  #originalSi = currentSi

global collectFrameCount
collectFrameCount = 0

def CollectCells2(inClientSideData, inBlockIndex, ioCollectList, ioTotalCellCount, cellVariableName, thresholdValue, thresholdDir):
  if PhactoriDbg():
    myDebugPrint3('  proc ' + str(SmartGetLocalProcessId()) + ' Entering CollectCells2: ' + inClientSideData.GetClassName() + '\n')
  if inClientSideData.GetClassName() == "vtkMultiBlockDataSet":
    numBlocks = inClientSideData.GetNumberOfBlocks()
    if PhactoriDbg():
      myDebugPrint3('  proc ' + str(SmartGetLocalProcessId()) + ' data is vtkMultiBlockDataSet, recursing, numblocks: ' + str(numBlocks) + '\n')
    for ii in range(0, numBlocks):
      oneBlock = inClientSideData.GetBlock(ii)
      if(oneBlock != None):
        CollectCells2(oneBlock, ii, ioCollectList, ioTotalCellCount, cellVariableName, thresholdValue, thresholdDir)
        #if oneBlock.GetClassName() == "vtkExodusIIMultiBlockDataSet":
        #  myDebugPrint3('  proc ' + str(SmartGetLocalProcessId()) + ' vtkExodusIIMultiBlockDataSet breaks loop\n')
        #  break
  elif inClientSideData.GetClassName() == "vtkExodusIIMultiBlockDataSet":
    numBlocks = inClientSideData.GetNumberOfBlocks()
    if PhactoriDbg():
      myDebugPrint3('  proc ' + str(SmartGetLocalProcessId()) + ' data is vtkExodusIIMultiBlockDataSet, not recursing, numblocks: ' + str(numBlocks) + '\n')
    for ii in range(0, numBlocks):
      oneBlock = inClientSideData.GetBlock(ii)
      if(oneBlock != None):
        CollectCells2(oneBlock, ii, ioCollectList, ioTotalCellCount, cellVariableName, thresholdValue, thresholdDir)
  else:
    cellData = inClientSideData.GetCellData()
    tearArray = cellData.GetArray(cellVariableName)
    #idArray = cellData.GetArray('GlobalElementId')
    numTuples = tearArray.GetNumberOfTuples()
    if PhactoriDbg():
      myDebugPrint3('  proc ' + str(SmartGetLocalProcessId()) + ' data is not vtkMultiBlockDataSet, collecting ' + str(numTuples) + '\n')
    global collectFrameCount
    appendCount = 0
    for jj in range(0, numTuples):
      ioTotalCellCount[0] = ioTotalCellCount[0] + 1
      tearVal = tearArray.GetTuple1(jj)
      testResult = False
      if thresholdDir == 1:
        if tearVal > thresholdValue:
          testResult = True
      else:
        if tearVal < thresholdValue:
          testResult = True
      if testResult:
        #elemId = idArray.GetTuple1(jj)
        listItem = [SmartGetLocalProcessId(), inBlockIndex, jj, tearVal, collectFrameCount]
        ioCollectList.append(listItem)
        if(appendCount == 0):
          if PhactoriDbg():
            myDebugPrint3('  proc ' + str(SmartGetLocalProcessId()) + 'listing first few appended (up to 20)\n')
        if(appendCount < 20):
          if PhactoriDbg():
            myDebugPrint3('  proc ' + str(SmartGetLocalProcessId()) + ' appending item: ' + str(listItem) + '\n')
        appendCount += 1
  if PhactoriDbg():
    myDebugPrint3('  proc ' + str(SmartGetLocalProcessId()) + ' Leaving CollectCells2\n')


def CollectCells1(cellVariableName, thresholdValue, thresholdDir):
  if PhactoriDbg():
    myDebugPrint3('Entering CollectCells1  var ' + cellVariableName + '\n')
  global collectFrameCount
  collectFrameCount = collectFrameCount + 1
  xx = GetActiveSource()
  collectedCells = []
  xxdata = xx.GetClientSideObject().GetOutputDataObject(0)
  totalCellCount = [0]
  CollectCells2(xxdata, 0, collectedCells, totalCellCount, cellVariableName, thresholdValue, thresholdDir)
  #myDebugPrint3('  collection result:\n' + str(collectedCells) + '\n')
  if PhactoriDbg():
    myDebugPrint3('  proc ' + str(SmartGetLocalProcessId()) + ' total cells: ' + str(totalCellCount) + ' number of collected cells: ' + str(len(collectedCells)) + '\n')
  if PhactoriDbg():
    myDebugPrint3('  last 20: \n')
  for ii in range(0,20):
    index = len(collectedCells) - 1 - ii
    if index >= 0:
      if PhactoriDbg():
        myDebugPrint3('    ' + str(collectedCells[index]) + ' ')
      if PhactoriDbg():
        myDebugPrint3('\n')
  if PhactoriDbg():
    myDebugPrint3('Leaving CollectCells1\n')
  return collectedCells

def UseReduceToSumArrayOfInts(ioListOfValues):
  if PhactoriDbg(100):
    myDebugPrint3('UseReduceToSumArrayOfInts entered\n', 100)
  import vtkParallelCorePython
  if PhactoriDbg():
    myDebugPrint3('  before reduced: ' + str(ioListOfValues) + '\n');
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()
  localarray = vtk.vtkIntArray()
  numberOfValues = len(ioListOfValues)
  localarray.SetNumberOfTuples(numberOfValues)
  for jj in range(0, numberOfValues):
    localarray.SetValue(jj, ioListOfValues[jj])
  globalarray = vtk.vtkIntArray()
  globalarray.SetNumberOfTuples(numberOfValues)
  globalController.AllReduce(localarray, globalarray, 2)
  for jj in range(0, numberOfValues):
    ioListOfValues[jj] = globalarray.GetTuple1(jj)
  if PhactoriDbg():
    myDebugPrint3('         reduced: ' + str(ioListOfValues) + '\n');
  if PhactoriDbg(100):
    myDebugPrint3('UseReduceToSumArrayOfInts returning\n', 100)

#called by UseReduceToSpreadValues to actually do Reduce
def UseReduceToSpreadPositiveValues(ioListOfValues):
  if PhactoriDbg(100):
    myDebugPrint3('UseReduceToSpreadPositiveValues entered\n', 100)
  import vtkParallelCorePython
  if PhactoriDbg():
    myDebugPrint3('  before reduced: ' + str(ioListOfValues) + '\n');
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()
  localarray = vtk.vtkDoubleArray()
  numberOfValues = len(ioListOfValues)
  localarray.SetNumberOfTuples(numberOfValues)
  for jj in range(0, numberOfValues):
    localarray.SetValue(jj, ioListOfValues[jj])
  globalarray = vtk.vtkDoubleArray()
  globalarray.SetNumberOfTuples(numberOfValues)
  globalController.AllReduce(localarray, globalarray, 0)
  for jj in range(0, numberOfValues):
    ioListOfValues[jj] = globalarray.GetTuple1(jj)
  if PhactoriDbg():
    myDebugPrint3('         reduced: ' + str(ioListOfValues) + '\n');
  if PhactoriDbg(100):
    myDebugPrint3('UseReduceToSpreadPositiveValues returning\n', 100)


#used by GetAndReduceViewControl to read view control values from a file on
#one process, then use Reduce to spread the values to all processes
def UseReduceToSpreadValues(ioListOfValues):
  #assumes values for all processes are 0.0 except one, uses trick to only do max

  if PhactoriDbg(100):
    myDebugPrint3('UseReduceToSpreadValues entered\n', 100)
  if PhactoriDbg():
    myDebugPrint3('  ' + 'proc ' + str(SmartGetLocalProcessId()) + ' start values: ' + str(ioListOfValues) + '\n')
  listOfValuePairs = []
  numberOfValues = len(ioListOfValues)
  for ii in range(0, numberOfValues):
    if ioListOfValues[ii] >= 0.0:
      listOfValuePairs.append(ioListOfValues[ii])
      listOfValuePairs.append(0.0)
    else:
      listOfValuePairs.append(0.0)
      listOfValuePairs.append(-ioListOfValues[ii])
  if PhactoriDbg():
    myDebugPrint3('  ' + 'proc ' + str(SmartGetLocalProcessId()) + ' paired values: ' + str(listOfValuePairs) + '\n')
  UseReduceToSpreadPositiveValues(listOfValuePairs)
  for ii in range(0, numberOfValues):
    value = listOfValuePairs[ii*2]
    if value == 0.0:
      value = listOfValuePairs[ii*2+1]
      if value != 0.0:
        value = -value
    ioListOfValues[ii] = value

  if PhactoriDbg():
    myDebugPrint3('  ' + 'proc ' + str(SmartGetLocalProcessId()) + ' final values: ' + str(ioListOfValues) + '\n')
  if PhactoriDbg(100):
    myDebugPrint3('UseReduceToSpreadValues returning\n', 100)

def GetCurrentSource():
  return GetActiveSource()

def dummyTest1():
  print '####in dummyTest1'

global gPlotView1Pts
gPlotView1Pts = None

global ProgrammableSourceInput
ProgrammableSourceInput = None
global ProgrammableSource2
ProgrammableSource2 = None

class PlotValMinMaxTrkC:
  def __init__(self):
    self.mThisCbMin = 0.0
    self.mThisCbMax = 0.0
    self.mThisCbInitialized = False
    self.mAllMin = 0.0
    self.mAllMax = 0.0
    self.mAllInitialized = False

    #sets minimum and maximum plot ranges, and allows there to not be
    #a minimmum and/or maximum as well
    self.mLowestTop = 0.0
    self.mUseLowestTop = False
    self.mLowestBot = 0.0
    self.mUseLowestBot = False
    self.mHighestTop = 0.0
    self.mUseHighestTop = False
    self.mHighestBot = 0.0
    self.mUseHighestBot = False

    self.mMinToUse = 0.0
    self.mMaxToUse = 0.0

    self.mUseCumulativeRange = False

  def SetFromRestartInfo(self, inJson):
    """given a map (json format), use the info in the map to set the
       tracker state--this reads the info created out in GetRestartInfo"""

    inJsonIsGood = True
    if 'mAll' not in inJson:
      if PhactoriDbg():
        myDebugPrint3(
            'PlotValMinMaxTrkC::SetFromRestartInfo mAll not in inJson\n')
      inJsonIsGood = False
    if 'mLowestTop' not in inJson:
      if PhactoriDbg():
        myDebugPrint3(
            'PlotValMinMaxTrkC::SetFromRestartInfo mLowestTop not in inJson\n')
      inJsonIsGood = False
    if 'mLowestBot' not in inJson:
      if PhactoriDbg():
        myDebugPrint3(
            'PlotValMinMaxTrkC::SetFromRestartInfo mLowestBot not in inJson\n')
      inJsonIsGood = False
    if 'mHighestTop' not in inJson:
      if PhactoriDbg():
        myDebugPrint3(
            'PlotValMinMaxTrkC::SetFromRestartInfo mHighestTop not in inJson\n')
      inJsonIsGood = False
    if 'mHighestBot' not in inJson:
      if PhactoriDbg():
        myDebugPrint3(
            'PlotValMinMaxTrkC::SetFromRestartInfo mHighestBot not in inJson\n')
      inJsonIsGood = False
    if 'mMinMaxToUse' not in inJson:
      if PhactoriDbg():
        myDebugPrint3(
            'PlotValMinMaxTrkC::SetFromRestartInfo mMinMaxToUse not in inJson\n')
      inJsonIsGood = False
    if 'mUseCumulativeRange' not in inJson:
      if PhactoriDbg():
        myDebugPrint3('PlotValMinMaxTrkC::SetFromRestartInfo " + \
            "mUseCumulativeRange not in inJson\n')
      inJsonIsGood = False
    if inJsonIsGood == False:
      return

    #thisCbItem = inJson['mThisCb']
    #self.mThisCbInitialized = thisCbItem[0]
    #self.mThisCbMin = thisCbItem[1]
    #self.mThisCbMax = thisCbItem[2]

    allItem = inJson['mAll']
    self.mAllInitialized = allItem[0]
    self.mAllMin = allItem[1]
    self.mAllMax = allItem[2]

    lowestTopItem = inJson['mLowestTop']
    self.mUseLowestTop = lowestTopItem[0]
    self.mLowestTop = lowestTopItem[1]

    lowestBotItem = inJson['mLowestBot']
    self.mUseLowestBot = lowestBotItem[0]
    self.mLowestBot = lowestBotItem[1]

    highestTopItem = inJson['mHighestTop']
    self.mUseHighestTop = highestTopItem[0]
    self.mHighestTop = highestTopItem[1]

    highestBotItem = inJson['mHighestBot']
    self.mUseHighestBot = highestBotItem[0]
    self.mHighestBot = highestBotItem[1]

    minMaxToUseItem = inJson['mMinMaxToUse']
    self.mMinToUse = minMaxToUseItem[0]
    self.mMaxToUse = minMaxToUseItem[1]

    self.mUseCumulativeRange = inJson["mUseCumulativeRange"]
    if PhactoriDbg():
      myDebugPrint3("PlotValMinMaxTrkC::SetFromRestartInfo state loaded:\n"
        + self.SelfToStr())

  def GetRestartInfo(self):
    """construct, in python map/json info from this PlotValMinMaxTrkC instance
       which contains the information which would be needed to restore the
       instance to the proper state after a simulation restart.
       Return the restart info map/json"""
    newJsonItem = {}
    #newJsonItem["mThisCb"] = [self.mThisCbInitialized,
    #    self.mThisCbMin, self.mThisCbMax]
    newJsonItem["mAll"] = [self.mAllInitialized,
        self.mAllMin, self.mAllMax]
    newJsonItem["mLowestTop"] = [self.mUseLowestTop, self.mLowestTop]
    newJsonItem["mLowestBot"] = [self.mUseLowestBot, self.mLowestBot]
    newJsonItem["mHighestTop"] = [self.mUseHighestTop, self.mHighestTop]
    newJsonItem["mHighestBot"] = [self.mUseHighestBot, self.mHighestBot]
    newJsonItem["mMinMaxToUse"] = [self.mMinToUse, self.mMaxToUse]
    newJsonItem["mUseCumulativeRange"] = self.mUseCumulativeRange
    if PhactoriDbg():
      myDebugPrint3("PlotValMinMaxTrkC::GetRestartInfo state saved:\n"
        + self.SelfToStr())
    return newJsonItem

  def SelfToStr(self):
    retStr = "PlotValMinMaxTrkC info:" +\
      "\n  mUseCumulativeRange: " + str(self.mUseCumulativeRange) +\
      "\n  mThisCbInitialized: " + str(self.mThisCbInitialized) +\
      "\n  mThisCbMin:  " + repr(self.mThisCbMin) +\
      "\n  mThisCbMax:  " + repr(self.mThisCbMax) +\
      "\n  mAllInitialized: " + str(self.mAllInitialized) +\
      "\n  mAllMin:     " + repr(self.mAllMin) +\
      "\n  mAllMax:     " + repr(self.mAllMax) +\
      "\n  mLowestBot:  " + repr(self.mLowestBot) +\
      "\n  mHighestBot: " + repr(self.mHighestBot) +\
      "\n  mLowestTop:  " + repr(self.mLowestTop) +\
      "\n  mHighestTop: " + str(self.mHighestTop) +\
      "\n  mUseLowestBot:  " + str(self.mUseLowestBot) +\
      "\n  mUseHighestBot: " + str(self.mUseHighestBot) +\
      "\n  mUseLowestTop:  " + str(self.mUseLowestTop) +\
      "\n  mUseHighestTop: " + str(self.mUseHighestTop) +\
      "\n  ---in hex--- " +\
      "\n  mThisCbMin:  " + self.mThisCbMin.hex() +\
      "\n  mThisCbMax:  " + self.mThisCbMax.hex() +\
      "\n  mAllMin:     " + self.mAllMin.hex() +\
      "\n  mAllMax:     " + self.mAllMax.hex() +\
      "\n  mLowestBot:  " + self.mLowestBot.hex() +\
      "\n  mHighestBot: " + self.mHighestBot.hex() +\
      "\n  mLowestTop:  " + self.mLowestTop.hex() +\
      "\n  mHighestTop: " + self.mHighestTop.hex() +\
      "\n"
    return retStr

  def PlotValMinMaxTrkCParseJson(self, inJson, inKeyPrefix):
    if PhactoriDbg(100):
      myDebugPrint3("PlotValMinMaxTrkCParseJson entered, key prefix is ->" + inKeyPrefix + "<-\n", 100)
    #self.mParseHadKeys = False
    keyCurrentDataRange = inKeyPrefix + 'use current data range'
    keyCumulativeDataRange = inKeyPrefix + 'use cumulative data range'
    keyRange = inKeyPrefix + 'range'
    keyMinimumRange = inKeyPrefix + 'minimum range'
    keyMaximumRange = inKeyPrefix + 'maximum range'

    if keyCurrentDataRange in inJson:
      self.mUseCumulativeRange = False
    elif keyCumulativeDataRange in inJson:
      self.mUseCumulativeRange = True

    if keyRange in inJson:
      range = inJson[keyRange]
      ff0 = float(range[0])
      ff1 = float(range[1])
      range = [ff0, ff1]
      if PhactoriDbg():
        myDebugPrint3("PlotValMinMaxTrkCParseJson got range " + str(range) + "\n")
      self.mLowestBot = range[0]
      self.mUseLowestBot = True
      self.mHighestBot = range[0]
      self.mUseHighestBot = True
      self.mLowestTop = range[1]
      self.mUseLowestTop = True
      self.mHighestTop = range[1]
      self.mUseHighestTop = True
    else:
      if keyMinimumRange in inJson:
        range = inJson[keyMinimumRange]
        ff0 = float(range[0])
        ff1 = float(range[1])
        range = [ff0, ff1]
        if PhactoriDbg():
          myDebugPrint3("PlotValMinMaxTrkCParseJson got minimum range " + str(range) + "\n")
        self.mHighestBot = range[0]
        self.mUseHighestBot = True
        self.mLowestTop = range[1]
        self.mUseLowestTop = True
      if keyMaximumRange in inJson:
        range = inJson[keyMaximumRange]
        ff0 = float(range[0])
        ff1 = float(range[1])
        range = [ff0, ff1]
        if PhactoriDbg():
          myDebugPrint3("PlotValMinMaxTrkCParseJson got maximum range " + str(range) + "\n")
        self.mLowestBot = range[0]
        self.mUseLowestBot = True
        self.mHighestTop = range[1]
        self.mUseHighestTop = True

    if PhactoriDbg():
      myDebugPrint3(self.SelfToStr())

  def GetMinToUse(self):
    return self.mMinToUse

  def GetMaxToUse(self):
    return self.mMaxToUse

  def ResetThisCb(self):
    """indicate we are starting a new callback search for min and max, so we
       can use the first test value as the initial min and max"""
    self.mThisCbInitialized = False
    self.mThisCbMin = 0.0
    self.mThisCbMax = 0.0

  def MinMaxTestAndSet(self, inTestXyz):
    if self.mAllInitialized:
      if inTestXyz < self.mAllMin:
        self.mAllMin = inTestXyz
      else:
        if inTestXyz > self.mAllMax:
          self.mAllMax = inTestXyz
    else:
      self.mAllMin = inTestXyz
      self.mAllMax = inTestXyz
      self.mAllInitialized = True
    if self.mThisCbInitialized:
      if inTestXyz < self.mThisCbMin:
        self.mThisCbMin = inTestXyz
      else:
        if inTestXyz > self.mThisCbMax:
          self.mThisCbMax = inTestXyz
    else:
      self.mThisCbMin = inTestXyz
      self.mThisCbMax = inTestXyz
      self.mThisCbInitialized = True

  def UpdateMinMaxToUse(self):
    """assumes we have obtained this callback (ThisCb) mins and maxes; uses
       Initial settings and locks to determine what mins and maxes are
       now in force"""
    if self.mUseCumulativeRange:
      localMin = self.mAllMin
    else:
      localMin = self.mThisCbMin
    if self.mUseHighestBot:
      if localMin > self.mHighestBot:
        localMin = self.mHighestBot
    if self.mUseLowestBot:
      if localMin < self.mLowestBot:
        localMin = self.mLowestBot

    if self.mUseCumulativeRange:
      localMax = self.mAllMax
    else:
      localMax = self.mThisCbMax
    if self.mUseHighestTop:
      if localMax > self.mHighestTop:
        localMax = self.mHighestTop
    if self.mUseLowestTop:
      if localMax < self.mLowestTop:
        localMax = self.mLowestTop

    if localMin > localMax:
      if self.mUseLowestBot:
        localMax = localMin
      if self.mUseHighestTop:
        localMin = localMax
      if localMin > localMax:
        if PhactoriDbg():
          myDebugPrint3("UpdateMinMaxToUse: weird case with min/max\n")
        localMin = localMax

    self.mMinToUse = localMin
    self.mMaxToUse = localMax

class PlotXYZMinMaxTrkC:
  def __init__(self):
    self.mXyzTrk = []
    self.mXyzTrk.append(PlotValMinMaxTrkC())
    self.mXyzTrk.append(PlotValMinMaxTrkC())
    self.mXyzTrk.append(PlotValMinMaxTrkC())

    self.currentYScaleFactor = 1.0

    #for scaling axes, use either cumulative data range or current

  def PlotXYZMinMaxTrkCParseJson(self, inJson, inKeyPrefix):
    if PhactoriDbg(100):
      myDebugPrint3("PlotXYZMinMaxTrkCParseJson entered, key prefix is ->" + inKeyPrefix + "<-\n", 100)
    xPrefix = 'x ' + inKeyPrefix
    yPrefix = 'y ' + inKeyPrefix
    zPrefix = 'z ' + inKeyPrefix
    self.mXyzTrk[0].PlotValMinMaxTrkCParseJson(inJson, xPrefix)
    self.mXyzTrk[1].PlotValMinMaxTrkCParseJson(inJson, yPrefix)
    self.mXyzTrk[2].PlotValMinMaxTrkCParseJson(inJson, zPrefix)

  def GetMinToUse(self):
    retVal = [self.mXyzTrk[0].GetMinToUse(), self.mXyzTrk[1].GetMinToUse(), self.mXyzTrk[2].GetMinToUse()]
    return retVal

  def GetMaxToUse(self):
    retVal = [self.mXyzTrk[0].GetMaxToUse(), self.mXyzTrk[1].GetMaxToUse(), self.mXyzTrk[2].GetMaxToUse()]
    return retVal

  def ResetThisCb(self):
    """indicate we are starting a new callback search for min and max, so we
       can use the first test value as the initial min and max"""
    self.mXyzTrk[0].ResetThisCb()
    self.mXyzTrk[1].ResetThisCb()
    self.mXyzTrk[2].ResetThisCb()

  def MinMaxTestAndSet(self, inTestXyz):
    """test a data value and adjust recorded mins and maxes.  Adjust overall
       min/max if necessary, adjust ThisCb min max if necessary, also see
       if this is first call since last ResetThisCb"""
    self.mXyzTrk[0].MinMaxTestAndSet(inTestXyz[0])
    self.mXyzTrk[1].MinMaxTestAndSet(inTestXyz[1])
    self.mXyzTrk[2].MinMaxTestAndSet(inTestXyz[2])

  def UpdateMinMaxToUse(self):
    """assumes we have obtained this callback (ThisCb) mins and maxes; uses
       Initial settings and locks to determine what mins and maxes are
       now in force"""
    self.mXyzTrk[0].UpdateMinMaxToUse()
    self.mXyzTrk[1].UpdateMinMaxToUse()
    self.mXyzTrk[2].UpdateMinMaxToUse()

  def SetFromRestartInfo(self, inJson):
    """given a map (json format), use the info in the map to set the
       state--this reads the info created out in GetRestartInfo"""
    if len(inJson) != 3:
      if PhactoriDbg():
        myDebugPrint3("PlotXYZMinMaxTrkC::SetFromRestartInfo bad inJson\n")
      return
    self.mXyzTrk[0].SetFromRestartInfo(inJson[0])
    self.mXyzTrk[1].SetFromRestartInfo(inJson[1])
    self.mXyzTrk[2].SetFromRestartInfo(inJson[2])

  def GetRestartInfo(self):
    """construct, in python map/json info from this PlotXYZMinMaxTrkC instance
       which contains the information which would be needed to restore the
       instance to the proper state after a simulation restart.
       Basically just calls the GetRestartInfo() for each underlying axis
       Return the restart info map/json"""
    #this is a list of three items, which is an acceptable json type
    newJsonInfo = []
    newJsonInfo.append(self.mXyzTrk[0].GetRestartInfo())
    newJsonInfo.append(self.mXyzTrk[1].GetRestartInfo())
    newJsonInfo.append(self.mXyzTrk[2].GetRestartInfo())
    return newJsonInfo


#class PlotColorInfo:
#  m_BackgroundColor = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]
#  m_EdgeColor = [0.0, 0.0, 0.5000076295109483]
#  m_DiffuseColor = [1.0, 1.0, 1.0]
#  m_AmbientColor = [1.0, 1.0, 1.0]
#  m_SelectionColor = [1.0, 0.0, 1.0]
#  m_BackfaceDiffuseColor = [1.0, 1.0, 1.0]
#  m_CubeAxesColor = [1.0, 1.0, 1.0]
#  def CopyTo(self, outTarget):
#    outTarget.m_BackgroundColor = self.m_BackgroundColor
#    outTarget.m_EdgeColor = self.m_EdgeColor
#    outTarget.m_DiffuseColor = self.m_DiffuseColor
#    outTarget.m_AmbientColor = self.m_AmbientColor
#    outTarget.m_SelectionColor = self.m_SelectionColor
#    outTarget.m_BackfaceDiffuseColor = self.m_BackfaceDiffuseColor
#    outTarget.m_CubeAxesColor = self.m_CubeAxesColor
#  def SetRepresentationAndViewColors(self, ioRenderViewInfo):
#    myDebugPrint3('PlotColorInfo.SetRepresentationColors entered\n', 100)
#    rep1 = ioRenderViewInfo.DataRepresentation1
#    view1 = ioRenderViewInfo.RenderView1
#    view1.Background = self.m_BackgroundColor
#    rep1.EdgeColor = self.m_EdgeColor
#    rep1.DiffuseColor = self.m_DiffuseColor
#    rep1.AmbientColor = self.m_AmbientColor
#    rep1.SelectionColor = self.m_SelectionColor
#    rep1.BackfaceDiffuseColor = self.m_BackfaceDiffuseColor
#    rep1.CubeAxesColor = self.m_CubeAxesColor
#    myDebugPrint3('PlotColorInfo.SetRepresentationColors returning\n', 100)

#global gPlotColorInfoPreset
#gPlotColorInfoPreset = PlotColorInfo()

global gPlotOverTimeMap
gPlotOverTimeMap = {}

global gScatterPlotMap
gScatterPlotMap = {}

#def SetPlotView2StartColors(
#  inBackgroundColor = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052],
#  inEdgeColor = [0.0, 0.0, 0.5000076295109483],
#  inDiffuseColor = [1.0, 1.0, 1.0],
#  inAmbientColor = [1.0, 1.0, 1.0],
#  inSelectionColor = [1.0, 0.0, 1.0,],
#  inBackfaceDiffuseColor = [1.0, 1.0, 1.0],
#  inCubeAxesColor = [1.0, 1.0, 1.0]):
#  global gPlotColorInfoPreset
#  gPlotColorInfoPreset.m_BackgroundColor = inBackgroundColor
#  gPlotColorInfoPreset.m_EdgeColor = inEdgeColor
#  gPlotColorInfoPreset.m_DiffuseColor = inDiffuseColor
#  gPlotColorInfoPreset.m_AmbientColor = inAmbientColor
#  gPlotColorInfoPreset.m_SelectionColor = inSelectionColor
#  gPlotColorInfoPreset.m_BackfaceDiffuseColor = inBackfaceDiffuseColor
#  gPlotColorInfoPreset.m_CubeAxesColor = inBackfaceDiffuseColor

def UpdatePlotView2(inPlotInfo):
  if PhactoriDbg(100):
    myDebugPrint3('UpdatePlotView2 entered (' + inPlotInfo.mName + ')\n', 100)
  if PhactoriDbg():
    myDebugPrint3('  render info(' + str(inPlotInfo.m_PhactoriRenderViewInfo) + ')\n')
  if PhactoriDbg():
    myDebugPrint3('  minmax info: ' + str(inPlotInfo.m_PhactoriRenderViewInfo) + '\n')
  #DataXYMinMax = [0.0, 0.0, 0.0, 0.0]
  UpdatePipelineWithCurrentTimeArgument(
      inPlotInfo.mInputOperation.GetPvFilter())
  UpdatePipelineWithCurrentTimeArgument(inPlotInfo.m_producer)
  SetPlotPointsFromData(inPlotInfo)
  UpdatePlotViewLook(inPlotInfo)
  if PhactoriDbg(100):
    myDebugPrint3('UpdatePlotView2 returning\n', 100)

def UpdateAllScatterPlots():
  if PhactoriDbg(100):
    myDebugPrint3('UpdateAllScatterPlots entered\n', 100)
  global gScatterPlotMap
  for ii in gScatterPlotMap.values():
    UpdatePlotView2(ii)
  if PhactoriDbg(100):
    myDebugPrint3('UpdateAllScatterPlots returning\n', 100)

def UpdateOnePlotOverTime(inPlotInfo):
  if PhactoriDbg(100):
    myDebugPrint3('UpdateOnePlotOverTime entered (' + inPlotInfo.mName + ')\n', 100)
  if PhactoriDbg():
    myDebugPrint3('  render info(' + str(inPlotInfo.m_PhactoriRenderViewInfo) + ')\n')
  if PhactoriDbg():
    myDebugPrint3('  minmax info: ' + str(inPlotInfo.m_xyzMinMaxTrkC) + '\n')
  #DataXYMinMax = [0.0, 0.0, 0.0, 0.0]
  UpdatePipelineWithCurrentTimeArgument(
      inPlotInfo.mInputOperation.GetPvFilter())
  UpdatePipelineWithCurrentTimeArgument(inPlotInfo.m_producer)
  SetPlotOverTimePointsFromData(inPlotInfo)
  UpdatePipelineWithCurrentTimeArgument(
      inPlotInfo.mInputOperation.GetPvFilter())
  UpdatePipelineWithCurrentTimeArgument(inPlotInfo.m_producer)
  UpdatePlotViewLook(inPlotInfo)
  if PhactoriDbg(100):
    myDebugPrint3('UpdateOnePlotOverTime returning\n', 100)

def UpdateAllPlotsOverTime():
  if PhactoriDbg(100):
    myDebugPrint3('UpdateAllPlotsOverTime entered\n', 100)
  global gPlotOverTimeMap
  for ii in gPlotOverTimeMap.values():
    UpdateOnePlotOverTime(ii)
  if PhactoriDbg(100):
    myDebugPrint3('UpdateAllPlotsOverTime returning\n', 100)

def PrintAllPointsA(inPlotInfo):
  ppdata = inPlotInfo.m_producer.GetClientSideObject().GetOutputDataObject(0)
  programmableFilterPoints = ppdata.GetPoints()
  if PhactoriDbg():
    myDebugPrint3('programmableFilterPoints: ' + str(programmableFilterPoints))
  for ii in range(0, programmableFilterPoints.GetNumberOfPoints()):
    onePoint = programmableFilterPoints.GetPoint(ii)
    if PhactoriDbg():
      myDebugPrint3(str(ii) + ': ' + str(onePoint))

#def TryToDetectVariableCellOrElementType(inInputCsData, inVariableName, inInputIsProxy = False, inAllowGlobalVariableDetection = False):
#  if PhactoriDbg():
#    myDebugPrint3('  trying to detect variable type (node/element)\n')
#
#  if inInputIsProxy:
#    pointData = inInputCsData.PointData
#  else:
#    pointData = inInputCsData.GetPointData()
#  if pointData != None:
#    testPointArray = pointData.GetArray(inVariableName)
#    if testPointArray != None:
#      if PhactoriDbg():
#        myDebugPrint3('  type node detected!\n')
#      return [True, 'node']
#
#  if inInputIsProxy:
#    cellData = inInputCsData.CellData
#  else:
#    cellData = inInputCsData.GetCellData()
#  if cellData != None:
#    testCellArray = cellData.GetArray(inVariableName)
#    if testCellArray != None:
#      if PhactoriDbg():
#        myDebugPrint3('  type element detected!\n')
#      return [True, 'element']
#
#  if inAllowGlobalVariableDetection:
#    if inInputIsProxy:
#      fieldData = inInputCsData.FieldData
#    else:
#      fieldData = inInputCsData.GetFieldData()
#    if fieldData != None:
#      testFieldArray = fieldData.GetArray(inVariableName)
#      if testFieldArray != None:
#        if PhactoriDbg():
#          myDebugPrint3('  type global detected!\n')
#        return [True, 'global']
#
#  return [False, '']

def FindThisProcessorMinMaxForVarForOneBlock(inInputCsData, inVariableInfo,
        ioMinMaxInfo):
  if PhactoriDbg(100):
    myDebugPrint3('FindThisProcessorMinMaxForVarForOneBlock entered\n', 100)
  if PhactoriDbg():
    myDebugPrint3(' variable is ' + inVariableInfo.mVariableName + '\n')

  #detect variable type (node/element) if necessary, and save info if detected
  detectResult = inVariableInfo.DetectVariableType(inInputCsData,
      False, False)
  if detectResult == False:
    if PhactoriDbg(100):
      myDebugPrint3('no detection: returning\n', 100)
    return

  if(inVariableInfo.mVariableType == 'node'):
    cellOrPointData = inInputCsData.GetPointData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no point data: returning\n', 100)
      return
  else:
    cellOrPointData = inInputCsData.GetCellData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no cell data: returning\n', 100)
      return

  theData = cellOrPointData.GetArray(inVariableInfo.mVariableName)
  if theData == None:
    if PhactoriDbg(100):
      myDebugPrint3('no theData: returning\n', 100)
    return
  if PhactoriDbg():
    myDebugPrint3('  theData: ' + str(theData) + '\n')
  if PhactoriDbg():
    myDebugPrint3('  theData name: ' + theData.GetName() + '\n')
  numTuplesX = theData.GetNumberOfTuples()

  if numTuplesX <= 0:
    return

  if ioMinMaxInfo[2] == False:
    if inVariableInfo.mVariableIsVectorComponent:
      vecVal = theData.GetTuple3(0)
      vv = vecVal[inVariableInfo.mVariableComponent]
    elif inVariableInfo.mVariableIsVectorMagnitude:
      vecVal = theData.GetTuple3(0)
      aa = vecVal[0]
      bb = vecVal[1]
      cc = vecVal[2]
      vv = math.sqrt(aa*aa + bb*bb + cc*cc)
    else:
      vv = theData.GetTuple1(0)

    ioMinMaxInfo[0] = vv
    ioMinMaxInfo[1] = vv
    ioMinMaxInfo[2] = True

  for ii in range(0, numTuplesX):
    if inVariableInfo.mVariableIsVectorComponent:
      vecVal = theData.GetTuple3(ii)
      vv = vecVal[inVariableInfo.mVariableComponent]
    elif inVariableInfo.mVariableIsVectorMagnitude:
      vecVal = theData.GetTuple3(ii)
      aa = vecVal[0]
      bb = vecVal[1]
      cc = vecVal[2]
      vv = math.sqrt(aa*aa + bb*bb + cc*cc)
    else:
      vv = theData.GetTuple1(ii)

    if(vv < ioMinMaxInfo[0]):
      ioMinMaxInfo[0] = vv
    if(vv > ioMinMaxInfo[1]):
      ioMinMaxInfo[1] = vv

  if PhactoriDbg():
    myDebugPrint3('  min max after this block: ' + str(ioMinMaxInfo) + '\n')
  if PhactoriDbg(100):
    myDebugPrint3('FindThisProcessorMinMaxForVarForOneBlock returning\n', 100)

def FindNodeOrElementIdForMinMaxOneBlock(inInputCsData, ioVariableInfo):
  if PhactoriDbg(100):
    myDebugPrint3('FindNodeOrElementIdForMinMaxOneBlock entered\n'
      ' variable is ' + ioVariableInfo.mVariableName + '\n'
      ' min: ' + str(ioVariableInfo.mStats.mMin) + '\n'
      ' max: ' + str(ioVariableInfo.mStats.mMax) + '\n', 100)

  #detect variable type (node/element) if necessary, and save info if detected
  detectResult = ioVariableInfo.DetectVariableType(inInputCsData,
      False, True)
  if detectResult == False:
    if PhactoriDbg(100):
      myDebugPrint3('no detection: returning\n', 100)
    return

  if(ioVariableInfo.mVariableType == 'node'):
    cellOrPointData = inInputCsData.GetPointData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no point data: returning\n', 100)
      return
  elif(ioVariableInfo.mVariableType == 'element'):
    cellOrPointData = inInputCsData.GetCellData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no cell data: returning\n', 100)
      return
  elif(ioVariableInfo.mVariableType == 'global'):
    cellOrPointData = inInputCsData.GetFieldData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no field data: returning\n', 100)
      return
  else:
      if PhactoriDbg(100):
        myDebugPrint3('no data type, should not be here,: returning\n', 100)
      return

  varData = cellOrPointData.GetArray(ioVariableInfo.mVariableName)
  if varData == None:
    if PhactoriDbg(100):
      myDebugPrint3('no varData: returning\n', 100)
    return
  if PhactoriDbg():
    myDebugPrint3('  y data: ' + str(varData) + '\n')
  if PhactoriDbg():
    myDebugPrint3('  y data name: ' + varData.GetName() + '\n')
  numTuplesX = varData.GetNumberOfTuples()

  if numTuplesX <= 0:
    return

  globalIdArray = None
  if(ioVariableInfo.mVariableType == 'node'):
    globalIdArray = cellOrPointData.GetArray('GlobalNodeId')
    if globalIdArray == None:
      if PhactoriDbg():
        myDebugPrint3('  no GlobalNodeId to find ids in list\n')
  else:
    globalIdArray = cellOrPointData.GetArray('GlobalElementId')
    if globalIdArray == None:
      if PhactoriDbg():
        myDebugPrint3('  no GlobalNodeId to find ids in list\n')

  for ii in range(0, numTuplesX):
    #vv = varData.GetTuple1(ii)
    if ioVariableInfo.mVariableIsVectorComponent:
      vecVal = varData.GetTuple3(ii)
      vv = vecVal[ioVariableInfo.mVariableComponent]
    elif ioVariableInfo.mVariableIsVectorMagnitude:
      vecVal = varData.GetTuple3(ii)
      aa = vecVal[0]
      bb = vecVal[1]
      cc = vecVal[2]
      vv = math.sqrt(aa*aa + bb*bb + cc*cc)
    else:
      vv = varData.GetTuple1(ii)

    if vv == ioVariableInfo.mStats.mMin:
      #if ioVariableInfo.mStats.mLocalFoundMinId == False:
      #  ioVariableInfo.mStats.mMinId = globalIdArray.GetTuple1(ii)
      #  ioVariableInfo.mStats.mLocalFoundMinId = True
      #  if PhactoriDbg(100):
      #    myDebugPrint3(" found min node/element id: " + \
      #      str(ioVariableInfo.mStats.mMinId) + "\n", 100)
      ioVariableInfo.mStats.mLocalMinIdCount += 1
      newGId = globalIdArray.GetTuple1(ii)
      if ioVariableInfo.mStats.mLocalMinIdCount <= 5:
        if PhactoriDbg(100):
          myDebugPrint3(" found min node/element id: " + \
            str(newGId) + "\n"
            " count: " + str(ioVariableInfo.mStats.mLocalMinIdCount) + "\n",
            100)
      if ioVariableInfo.mStats.mLocalFoundMinId == False:
        ioVariableInfo.mStats.mMinId = newGId
        ioVariableInfo.mStats.mLocalFoundMinId = True
      else:
        if newGId > ioVariableInfo.mStats.mMinId:
          ioVariableInfo.mStats.mMinId = newGId

    if vv == ioVariableInfo.mStats.mMax:
      #if ioVariableInfo.mStats.mLocalFoundMaxId == False:
      #  ioVariableInfo.mStats.mMaxId = globalIdArray.GetTuple1(ii)
      #  ioVariableInfo.mStats.mLocalFoundMaxId = True
      #  if PhactoriDbg(100):
      #    myDebugPrint3(" found max node/element id: " + \
      #      str(ioVariableInfo.mStats.mMaxId) + "\n", 100)
      ioVariableInfo.mStats.mLocalMaxIdCount += 1
      newGId = globalIdArray.GetTuple1(ii)
      if ioVariableInfo.mStats.mLocalMaxIdCount <= 5:
        if PhactoriDbg(100):
          myDebugPrint3(" found max node/element id: " + \
            str(newGId) + "\n"
            " count: " + str(ioVariableInfo.mStats.mLocalMaxIdCount) + "\n",
            100)
      if ioVariableInfo.mStats.mLocalFoundMaxId == False:
        ioVariableInfo.mStats.mMaxId = newGId
        ioVariableInfo.mStats.mLocalFoundMaxId = True
      else:
        if newGId > ioVariableInfo.mStats.mMaxId:
          ioVariableInfo.mStats.mMaxId = newGId

  if PhactoriDbg(100):
    myDebugPrint3('FindNodeOrElementIdForMinMaxOneBlock returning\n', 100)

def FindNodeOrElementIdForMinMaxRecurse1(inInputCsData, ioVariableInfo):
  #myDebugPrint3('FindNodeOrElementIdForMinMaxRecurse1 entered\n', 100)

  icsdClassname = inInputCsData.GetClassName()
  if icsdClassname == "vtkMultiBlockDataSet" or \
     icsdClassname == "vtkExodusIIMultiBlockDataSet":
    #myDebugPrint3('recursing: ' + icsdClassname + '\n')
    numBlocks = inInputCsData.GetNumberOfBlocks()
    for ii in range(0, numBlocks):
      oneBlock = inInputCsData.GetBlock(ii)
      if(oneBlock != None):
        FindNodeOrElementIdForMinMaxRecurse1(oneBlock, ioVariableInfo)
  else:
    #myDebugPrint3('finding min/max: ' + icsdClassname + '\n')
    FindNodeOrElementIdForMinMaxOneBlock(inInputCsData, ioVariableInfo)

  #myDebugPrint3('FindNodeOrElementIdForMinMaxRecurse1 returning\n', 100)

def FindNodeOrElementIdForMinMax(inputcsData, inVariableInfo):
  """fairly complex high level mpi-operation routine.  Takes a variable,
     finds the min/max (global) and finds a node id or element id corresponding
     to the node or element for both the min and the max value.  Has some
     intelligence to avoid recaculating and doing mpi communication if the
     information had already been found for this catalyst callback"""

  if PhactoriDbg(100):
    myDebugPrint3("FindNodeOrElementIdForMinMax entered: " + \
      str(inVariableInfo.mVariableName) + "\n", 100)

  #if already been done for this catalyst callback, just return
  if inVariableInfo.mStats.mIdsTestCounter == \
      gPipeAndViewsState.mFrameTagCounter:
    if PhactoriDbg(100):
      myDebugPrint3("already calculated for this callback, returning\n", 100)
    return

  #if variable is global, this shouldn't be called
  if inVariableInfo.mVariableType == "global":
    myDebugPrint3AndException("FindNodeOrElementIdForMinMax:\n"
      "this should not be called with global variable " + \
      inVariableInfo.mVariableName + "\n")

  #first get the min/max for the variable if necessary
  if inVariableInfo.mStats.mStatsTestCounter != \
      gPipeAndViewsState.mFrameTagCounter:
    if PhactoriDbg(100):
      myDebugPrint3("min/max not available, calculating\n", 100)
    DataMinMax = [0.0, 0.0, False]
    DataSumCnt = [0.0, 0]
    FindMinMaxSumCntFromData(inputcsData, inVariableInfo,
        None, DataMinMax, DataSumCnt,
        None, None)
  else:
    if PhactoriDbg(100):
      myDebugPrint3("min/max previously calculated\n", 100)

  #next, get the local node/element id, if any, which contain min and max
  #values
  inVariableInfo.mStats.mLocalFoundMinId = False
  inVariableInfo.mStats.mLocalMinIdCount = 0
  inVariableInfo.mStats.mMinId = -1
  inVariableInfo.mStats.mLocalFoundMaxId = False
  inVariableInfo.mStats.mLocalMaxIdCount = 0
  inVariableInfo.mStats.mMaxId = -1
  FindNodeOrElementIdForMinMaxRecurse1(inputcsData, inVariableInfo)

  #now get every process to agree on which node or element id contains the
  #max and min value; for multiple instances we choose the id with the highest
  #integer value
  idForMinAndidForMax = [inVariableInfo.mStats.mMinId,
      inVariableInfo.mStats.mMaxId]
  UseReduceToGetMaximumIntegers(idForMinAndidForMax)
  inVariableInfo.mStats.mMinId = idForMinAndidForMax[0]
  inVariableInfo.mStats.mMaxId = idForMinAndidForMax[1]

  #mark done for this catalyst callback, so we won't do extra mpi
  #communication if we call this routine again
  inVariableInfo.mStats.mIdsTestCounter = gPipeAndViewsState.mFrameTagCounter

  if PhactoriDbg(100):
    myDebugPrint3(
      "counter:   " + str(inVariableInfo.mStats.mIdsTestCounter) + "\n"
      "min:   " + str(inVariableInfo.mStats.mMin) + "\n"
      "at id: " + str(inVariableInfo.mStats.mMinId) + "\n"
      "max:   " + str(inVariableInfo.mStats.mMax) + "\n"
      "at id: " + str(inVariableInfo.mStats.mMaxId) + "\n"
      "FindNodeOrElementIdForMinMax returning\n")

def FindMinMaxSumCntFromDataOneBlock(inInputCsData, ioVariableInfo,
        ioSecondaryVariableInfo, ioDataMinMax, ioDataSumCnt,
        ioDataForIds, inPlotLineList):
  if PhactoriDbg(100):
    myDebugPrint3('FindMinMaxSumCntFromDataOneBlock entered\n', 100)

  if PhactoriDbg():
    myDebugPrint3(' y axis variable is ' + ioVariableInfo.mVariableName + '\n')

  #detect variable type (node/element) if necessary, and save info if detected
  detectResult = ioVariableInfo.DetectVariableType(inInputCsData,
      False, True)
  if detectResult == False:
    if PhactoriDbg(100):
      myDebugPrint3('no detection: returning\n', 100)
    return

  if ioSecondaryVariableInfo != None:
    ioSecondaryVariableInfo.CopyVariableTypeFrom(ioVariableInfo)

  if(ioVariableInfo.mVariableType == 'node'):
    cellOrPointData = inInputCsData.GetPointData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no point data: returning\n', 100)
      return
  elif(ioVariableInfo.mVariableType == 'element'):
    cellOrPointData = inInputCsData.GetCellData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no cell data: returning\n', 100)
      return
  elif(ioVariableInfo.mVariableType == 'global'):
    cellOrPointData = inInputCsData.GetFieldData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no field data: returning\n', 100)
      return
  else:
      if PhactoriDbg(100):
        myDebugPrint3('no data type, should not be here,: returning\n', 100)
      return

  #print all element names
  #myDebugPrint3('printing array names:\n')
  #numCellArrays = cellOrPointData.GetNumberOfArrays()
  #for ii in range(0, numCellArrays):
  #  myDebugPrint3(str(ii) + ": " + cellOrPointData.GetArray(ii).GetName() + "\n")

  yData = cellOrPointData.GetArray(ioVariableInfo.mVariableName)
  if yData == None:
    if PhactoriDbg(100):
      myDebugPrint3('no yData: returning\n', 100)
    return
  if PhactoriDbg():
    myDebugPrint3('  y data: ' + str(yData) + '\n')
  if PhactoriDbg():
    myDebugPrint3('  y data name: ' + yData.GetName() + '\n')
  numTuplesX = yData.GetNumberOfTuples()

  if numTuplesX <= 0:
    return

  if ioDataMinMax[2] == False:
    #vv = yData.GetTuple1(0)
    if ioVariableInfo.mVariableIsVectorComponent:
      vecVal = yData.GetTuple3(0)
      vv = vecVal[ioVariableInfo.mVariableComponent]
    elif ioVariableInfo.mVariableIsVectorMagnitude:
      vecVal = yData.GetTuple3(0)
      aa = vecVal[0]
      bb = vecVal[1]
      cc = vecVal[2]
      vv = math.sqrt(aa*aa + bb*bb + cc*cc)
    else:
      vv = yData.GetTuple1(0)

    ioDataMinMax[0] = vv
    ioDataMinMax[1] = vv
    ioDataMinMax[2] = True

  for ii in range(0, numTuplesX):
    #vv = yData.GetTuple1(ii)
    if ioVariableInfo.mVariableIsVectorComponent:
      vecVal = yData.GetTuple3(ii)
      vv = vecVal[ioVariableInfo.mVariableComponent]
    elif ioVariableInfo.mVariableIsVectorMagnitude:
      vecVal = yData.GetTuple3(ii)
      aa = vecVal[0]
      bb = vecVal[1]
      cc = vecVal[2]
      vv = math.sqrt(aa*aa + bb*bb + cc*cc)
    else:
      vv = yData.GetTuple1(ii)

    if(vv < ioDataMinMax[0]):
      ioDataMinMax[0] = vv
    if(vv > ioDataMinMax[1]):
      ioDataMinMax[1] = vv
    ioDataSumCnt[0] += vv
    ioDataSumCnt[1] += 1

  if ioDataForIds != None:
    numIds = len(inPlotLineList)
    globalIdArray = None
    if(ioVariableInfo.mVariableType == 'node'):
      globalIdArray = cellOrPointData.GetArray('GlobalNodeId')
      if globalIdArray == None:
        if PhactoriDbg():
          myDebugPrint3('  no GlobalNodeId to find ids in list\n')
    else:
      globalIdArray = cellOrPointData.GetArray('GlobalElementId')
      if globalIdArray == None:
        if PhactoriDbg():
          myDebugPrint3('  no GlobalNodeId to find ids in list\n')
    if globalIdArray != None:
      for ii in range(0, numIds):
        thisIdPlotLine = inPlotLineList[ii]
        idToFind = thisIdPlotLine.m_Id
        if PhactoriDbg():
          myDebugPrint3('    trying to find id: ' + str(idToFind) + '\n')
        if thisIdPlotLine.m_FoundOnThisProcessorLastTime == True:
          if inInputCsData == thisIdPlotLine.m_FoundBlockRef:
            idFromMesh = globalIdArray.GetTuple1(thisIdPlotLine.m_FoundIndex)
            if idFromMesh == idToFind:
              if PhactoriDbg():
                myDebugPrint3('  found id ' + str(idToFind) + ' on this processor again no search\n')
              #vv = yData.GetTuple1(thisIdPlotLine.m_FoundIndex)
              if ioVariableInfo.mVariableIsVectorComponent:
                vecVal = yData.GetTuple3(thisIdPlotLine.m_FoundIndex)
                vv = vecVal[ioVariableInfo.mVariableComponent]
              elif ioVariableInfo.mVariableIsVectorMagnitude:
                vecVal = yData.GetTuple3(thisIdPlotLine.m_FoundIndex)
                aa = vecVal[0]
                bb = vecVal[1]
                cc = vecVal[2]
                vv = math.sqrt(aa*aa + bb*bb + cc*cc)
              else:
                vv = yData.GetTuple1(thisIdPlotLine.m_FoundIndex)
              ioDataForIds[ii] = vv
              #go on to next id
              continue
            else:
              thisIdPlotLine.m_FoundOnThisProcessorLastTime = False
        for jj in range(0, numTuplesX):
          idFromMesh = globalIdArray.GetTuple1(jj)
          if idFromMesh == idToFind:
            if PhactoriDbg():
              myDebugPrint3('  found id ' + str(idToFind) + ' on this processor\n')
            #vv = yData.GetTuple1(jj)
            if ioVariableInfo.mVariableIsVectorComponent:
              vecVal = yData.GetTuple3(jj)
              vv = vecVal[ioVariableInfo.mVariableComponent]
            elif ioVariableInfo.mVariableIsVectorMagnitude:
              vecVal = yData.GetTuple3(jj)
              aa = vecVal[0]
              bb = vecVal[1]
              cc = vecVal[2]
              vv = math.sqrt(aa*aa + bb*bb + cc*cc)
            else:
              vv = yData.GetTuple1(jj)
            ioDataForIds[ii] = vv
            thisIdPlotLine.m_FoundOnThisProcessorLastTime = True
            thisIdPlotLine.m_FoundBlockRef = inInputCsData
            thisIdPlotLine.m_FoundIndex = jj

  if PhactoriDbg():
    myDebugPrint3('  min max after this block: ' + str(ioDataMinMax) + '\n')

  if PhactoriDbg(100):
    myDebugPrint3('FindMinMaxSumCntFromDataOneBlock returning\n', 100)


def SetPlotPointsFromOneBlock(inInputCsData, ioPlotInfo, ioIndex):
  if PhactoriDbg(100):
    myDebugPrint3('SetPlotPointsFromOneBlock entered\n', 100)
  if PhactoriDbg():
    myDebugPrint3(' x axis variable: ' + ioPlotInfo.m_XAxisVariableInfo.mVariableName + \
       '\n y axis variable: ' + ioPlotInfo.m_YAxisVariableInfo.mVariableName + '\n')

  #detect variable type (node/element) if necessary, and save info if detected
  detectResult = ioPlotInfo.m_YAxisVariableInfo.DetectVariableType(
      inInputCsData, False, False)
  if detectResult == False:
    if PhactoriDbg(100):
      myDebugPrint3('no detection: returning\n', 100)
    return

  ioPlotInfo.m_XAxisVariableInfo.CopyVariableTypeFrom(
      ioPlotInfo.m_YAxisVariableInfo)

  if(ioPlotInfo.m_YAxisVariableInfo.mVariableType == 'node'):
    cellOrPointData = inInputCsData.GetPointData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no point data: returning\n', 100)
      return
  else:
    cellOrPointData = inInputCsData.GetCellData()
    if cellOrPointData == None:
      if PhactoriDbg(100):
        myDebugPrint3('no cell data: returning\n', 100)
      return
  #nodeData = inClientSideData.GetNodeData()


  #print all element names
  #myDebugPrint3('printing array names:\n')
  #numCellArrays = cellOrPointData.GetNumberOfArrays()
  #for ii in range(0, numCellArrays):
  #  myDebugPrint3(str(ii) + ": " + cellOrPointData.GetArray(ii).GetName() + "\n")


  xData = cellOrPointData.GetArray(ioPlotInfo.m_XAxisVariableInfo.mVariableName)
  yData = cellOrPointData.GetArray(ioPlotInfo.m_YAxisVariableInfo.mVariableName)
  if xData == None:
    if PhactoriDbg(100):
      myDebugPrint3('no xData: returning\n', 100)
    return
  if yData == None:
    if PhactoriDbg(100):
      myDebugPrint3('no yData: returning\n', 100)
    return
  if PhactoriDbg():
    myDebugPrint3('  x data: ' + str(xData) + '  y data: ' + str(yData) + '\n')
  if PhactoriDbg():
    myDebugPrint3('  x data name: ' + xData.GetName() + '  y data name : ' + yData.GetName() + '\n')
  numTuplesX = xData.GetNumberOfTuples()
  numTuplesY = yData.GetNumberOfTuples()

  if numTuplesX != numTuplesY:
    if PhactoriDbg():
      myDebugPrint3('number of tuples do not match in SetPlotPointsFromOneBlock')
    return

  if 1:
    ptListCount = numTuplesX
    if ptListCount > 100:
      ptListCount = 100
    for jj in range(0, ptListCount):
      if PhactoriDbg():
        myDebugPrint3('  cell ' + str(jj) + ': ' + str(xData.GetTuple1(jj)) + '  data: ' + str(yData.GetTuple1(jj)) + '\n')

  outOutputPoints = ioPlotInfo.m_Points

  #myDebugPrint3('ioIndex is ' + str(ioIndex))
  #myDebugPrint3('outOutputPoints is ' + str(outOutputPoints))
  #myDebugPrint3('outOutputPoints number of points is  ' + str(outOutputPoints.GetNumberOfPoints()))

  localIndex = ioIndex[0]
  for ii in range(0, numTuplesX):

    if ioPlotInfo.m_XAxisVariableInfo.mVariableIsVectorComponent:
      vecVal = xData.GetTuple3(ii)
      xx = vecVal[ioPlotInfo.m_XAxisVariableInfo.mVariableComponent]
    elif ioPlotInfo.m_XAxisVariableInfo.mVariableIsVectorMagnitude:
      vecVal = xData.GetTuple3(ii)
      aa = vecVal[0]
      bb = vecVal[1]
      cc = vecVal[2]
      xx = math.sqrt(aa*aa + bb*bb + cc*cc)
    else:
      xx = xData.GetTuple1(ii)

    if ioPlotInfo.m_YAxisVariableInfo.mVariableIsVectorComponent:
      vecVal = yData.GetTuple3(ii)
      yy = vecVal[ioPlotInfo.m_YAxisVariableInfo.mVariableComponent]
    elif ioPlotInfo.m_YAxisVariableInfo.mVariableIsVectorMagnitude:
      vecVal = yData.GetTuple3(ii)
      aa = vecVal[0]
      bb = vecVal[1]
      cc = vecVal[2]
      yy = math.sqrt(aa*aa + bb*bb + cc*cc)
    else:
      yy = yData.GetTuple1(ii)

    zz = 0.0
    if localIndex < outOutputPoints.GetNumberOfPoints():
      if PhactoriDbg():
        myDebugPrint3(' sppfob changing existing point ' + str(localIndex) +': ' + str(xx) + ', ' + str(yy) + ', ' + str(zz) + '\n')
      outOutputPoints.SetPoint(localIndex, xx, yy, zz)
    else:
      if PhactoriDbg():
        myDebugPrint3(' sppfob inserting new point ' + str(localIndex) +': ' + str(xx) + ', ' + str(yy) + ', ' + str(zz) + '\n')
      outOutputPoints.InsertPoint(localIndex, xx, yy, zz)
    localIndex += 1
  ioIndex[0] = localIndex
  #PrintAllPointsA()
  if PhactoriDbg(100):
    myDebugPrint3('SetPlotPointsFromOneBlock returning\n', 100)

def FindMinMaxSumCntFromDataRecurse1(inInputCsData, ioVariableInfo,
        ioSecondaryVariableInfo, ioDataMinMax, ioDataSumCnt,
        ioDataForIds, inPlotIdLineList):
  #myDebugPrint3('FindMinMaxSumCntFromDataRecurse1 entered\n', 100)

  icsdClassname = inInputCsData.GetClassName()
  if icsdClassname == "vtkMultiBlockDataSet" or \
     icsdClassname == "vtkExodusIIMultiBlockDataSet":
    #myDebugPrint3('recursing: ' + icsdClassname + '\n')
    numBlocks = inInputCsData.GetNumberOfBlocks()
    for ii in range(0, numBlocks):
      oneBlock = inInputCsData.GetBlock(ii)
      if(oneBlock != None):
        FindMinMaxSumCntFromDataRecurse1(oneBlock,
            ioVariableInfo, ioSecondaryVariableInfo,
            ioDataMinMax, ioDataSumCnt, ioDataForIds, inPlotIdLineList)
  else:
    #myDebugPrint3('finding min/max: ' + icsdClassname + '\n')
    FindMinMaxSumCntFromDataOneBlock(inInputCsData,
        ioVariableInfo, ioSecondaryVariableInfo,
        ioDataMinMax, ioDataSumCnt, ioDataForIds, inPlotIdLineList)

  #myDebugPrint3('  min max after this recursion: ' + str(ioDataMinMax) + '\n')

  #myDebugPrint3('FindMinMaxFromDataSumCntRecurse1 returning\n', 100)


def FindThisProcessorMinMaxForVarRecurse1(inInputCsData, inVariableInfo, ioMinMaxInfo):
  #myDebugPrint3('FindThisProcessorMinMaxForVarRecurse1 entered\n', 100)

  icsdClassname = inInputCsData.GetClassName()
  if icsdClassname == "vtkMultiBlockDataSet" or \
     icsdClassname == "vtkExodusIIMultiBlockDataSet":
    #myDebugPrint3('recursing: ' + icsdClassname + '\n')
    numBlocks = inInputCsData.GetNumberOfBlocks()
    for ii in range(0, numBlocks):
      oneBlock = inInputCsData.GetBlock(ii)
      if(oneBlock != None):
        FindThisProcessorMinMaxForVarRecurse1(oneBlock, inVariableInfo, ioMinMaxInfo)
  else:
    FindThisProcessorMinMaxForVarForOneBlock(inInputCsData, inVariableInfo, ioMinMaxInfo)
  #myDebugPrint3('FindThisProcessorMinMaxForVarRecurse1 returning\n', 100)

def SetPlotPointsFromDataRecurse1(inInputCsData, ioPlotInfo, ioIndex):
  #myDebugPrint3('SetPlotPointsFromDataRecurse1 entered\n', 100)

  icsdClassname = inInputCsData.GetClassName()
  if icsdClassname == "vtkMultiBlockDataSet" or \
     icsdClassname == "vtkExodusIIMultiBlockDataSet":
    #myDebugPrint3('recursing: ' + icsdClassname + '\n')
    numBlocks = inInputCsData.GetNumberOfBlocks()
    for ii in range(0, numBlocks):
      oneBlock = inInputCsData.GetBlock(ii)
      if(oneBlock != None):
        SetPlotPointsFromDataRecurse1(oneBlock, ioPlotInfo, ioIndex)
  else:
    #myDebugPrint3('making points: ' + icsdClassname + '\n')
    SetPlotPointsFromOneBlock(inInputCsData, ioPlotInfo, ioIndex)
  #myDebugPrint3('SetPlotPointsFromDataRecurse1 returning\n', 100)

def FindPlotXyzMinMaxPtLst(inPlotPointsList, ioPlotXYZMinMaxTrkC):
  if len(inPlotPointsList) < 1:
    errStr = 'FindPlotXyzMinMaxPtLst needs to have at least one item in list'
    if PhactoriDbg():
      myDebugPrint3(errStr)
    raise Exception(errStr)

  #if ioPlotXYZMinMaxTrkC.mThisDbInitized == False)
  #  myDebugPrint3('no points, returning\n', 100)
  #  return

  firstPointForInitialization = None
  for onePointGroup in inPlotPointsList:
    if onePointGroup.GetNumberOfPoints() > 0:
      firstPointForInitialization = onePointGroup.GetPoint(0)
      break

  if firstPointForInitialization == None:
    if PhactoriDbg(100):
      myDebugPrint3('no points, returning\n', 100)
    return

  ioPlotXYZMinMaxTrkC.ResetThisCb()

  if PhactoriDbg():
    myDebugPrint3('  xyzminmax C incoming: \n' \
      '    x: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[0].mThisCbMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[0].mThisCbMax) + '\n' \
      '    y: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[1].mThisCbMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[1].mThisCbMax) + '\n' \
      '    z: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[2].mThisCbMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[2].mThisCbMax) + '\n' \
      '    xot: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[0].mAllMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[0].mAllMax) + '\n' \
      '    yot: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[1].mAllMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[1].mAllMax) + '\n' \
      '    zot: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[2].mAllMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[2].mAllMax) + '\n' \
    )
  if PhactoriDbg():
    myDebugPrint3(str(ioPlotXYZMinMaxTrkC))

  for onePointGroup in inPlotPointsList:
    numPoints = onePointGroup.GetNumberOfPoints()
    for ii in range(0, numPoints):
      onePoint = onePointGroup.GetPoint(ii)
      ioPlotXYZMinMaxTrkC.MinMaxTestAndSet(onePoint)

  xmin = ioPlotXYZMinMaxTrkC.mXyzTrk[0].mThisCbMin
  ymin = ioPlotXYZMinMaxTrkC.mXyzTrk[1].mThisCbMin
  zmin = ioPlotXYZMinMaxTrkC.mXyzTrk[2].mThisCbMin
  xmax = ioPlotXYZMinMaxTrkC.mXyzTrk[0].mThisCbMax
  ymax = ioPlotXYZMinMaxTrkC.mXyzTrk[1].mThisCbMax
  zmax = ioPlotXYZMinMaxTrkC.mXyzTrk[2].mThisCbMax
  minmaxlist = [xmin, xmax, ymin, ymax, zmin, zmax]
  UseReduceToGetMinMaxPairs(minmaxlist)
  xmin = minmaxlist[0]
  xmax = minmaxlist[1]
  ymin = minmaxlist[2]
  ymax = minmaxlist[3]
  zmin = minmaxlist[4]
  zmax = minmaxlist[5]

  UseReduceToGetMinMaxPairs(minmaxlist)
  ioPlotXYZMinMaxTrkC.MinMaxTestAndSet([xmin, ymin, zmin])
  ioPlotXYZMinMaxTrkC.MinMaxTestAndSet([xmax, ymax, zmax])

  ioPlotXYZMinMaxTrkC.UpdateMinMaxToUse()

  if PhactoriDbg():
    myDebugPrint3('  xyzminmax C outgoing: \n' \
      '    x: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[0].mThisCbMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[0].mThisCbMax) + '\n' \
      '    y: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[1].mThisCbMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[1].mThisCbMax) + '\n' \
      '    z: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[2].mThisCbMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[2].mThisCbMax) + '\n' \
      '    xot: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[0].mAllMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[0].mAllMax) + '\n' \
      '    yot: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[1].mAllMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[1].mAllMax) + '\n' \
      '    zot: ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[2].mAllMin) + ', ' + str(ioPlotXYZMinMaxTrkC.mXyzTrk[2].mAllMax) + '\n' \
    )
  if PhactoriDbg():
    myDebugPrint3(str(ioPlotXYZMinMaxTrkC))


def FindPlotXYZMinMax(inPlotPoints, ioPlotXYZMinMaxTrkC):
  FindPlotXyzMinMaxPtLst([inPlotPoints], ioPlotXYZMinMaxTrkC)

def ScalePlotYForFit(ioPlotPts, ioPlotXYZMinMaxTrkC, inAspectRatio):

  xyzmin = ioPlotXYZMinMaxTrkC.GetMinToUse()
  xmin = xyzmin[0]
  ymin = xyzmin[1]
  zmin = xyzmin[2]
  xyzmax = ioPlotXYZMinMaxTrkC.GetMaxToUse()
  xmax = xyzmax[0]
  ymax = xyzmax[1]
  zmax = xyzmax[2]

  xSpan = xmax - xmin
  ySpan = ymax - ymin
  if PhactoriDbg():
    myDebugPrint3('xSpan: ' + str(xSpan) + '  ySpan: ' + str(ySpan))
  yScaleFactor = 1.0
  if PhactoriDbg():
    myDebugPrint3('inAspectRatio: ' + str(inAspectRatio) + '\n')

  if xSpan > 0.0:
    if ySpan > 0.0:
      yScaleFactor = xSpan / ySpan
      yScaleFactor /= inAspectRatio

  ioPlotXYZMinMaxTrkC.currentYScaleFactor = yScaleFactor
  if PhactoriDbg():
    myDebugPrint3('yScaleFactor: ' + str(yScaleFactor))

  numPoints = ioPlotPts.GetNumberOfPoints()
  for ii in range(0, numPoints):
    onePoint = ioPlotPts.GetPoint(ii)
    xx = onePoint[0]
    yy = onePoint[1]
    zz = onePoint[2]
    ioPlotPts.SetPoint(ii, xx, yy*yScaleFactor, zz)

def FindThisProcessorMinMaxForVar(inOperation, inVariableInfo):

  myDataMinMax = [0.0, 0.0, False]
  pvClientSideData = inOperation.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
  if pvClientSideData == None:
    if PhactoriDbg(100):
      myDebugPrint3('pvClientSideData is None, returning', 100)
    return myDataMinMax

  FindThisProcessorMinMaxForVarRecurse1(pvClientSideData, inVariableInfo, myDataMinMax)

  return myDataMinMax


def FindMinMaxSumCntFromData(inputcsData, inVariableInfo,
        inSecondaryVariableInfo, DataMinMax, DataSumCnt,
        DataForIds, inPlotIdLineList):
  """fairly complex high level mpi-operation routine.  This takes a paraview
   filter client side object, gets the client side data
   (which is presumably a multiblock
   dataset) and finds the min, max, sum, and count of items for a variable.
   It uses MPI (if the item is not global) to share values and get them
   across all processors.  Used for plotting min/max/mean and for evaluating
   boolean criteria (such as whether or not to produce images) for min/max
   /mean/count.  Data return is a bit messy and should be updated--it fills
   in DataMinMax, DataSumCnt, and DataForIds if it is not None"""

  if PhactoriDbg():
    myDebugPrint3("FindMinMaxSumCntFromData entered\n")

  FindMinMaxSumCntFromDataRecurse1(inputcsData, inVariableInfo,
        inSecondaryVariableInfo, DataMinMax, DataSumCnt,
        DataForIds, inPlotIdLineList)

  #if variable is global, it is same on all processors so we don't need to do
  #any mpi communication
  if inVariableInfo.mVariableType != 'global':
    #do mpi to share values around
    UseReduceToSumDoubleAndIntPair(DataSumCnt)
    if DataSumCnt[1] == 0:
      #no data; we shouldn't add anything to plot lines
      if PhactoriDbg(100):
        myDebugPrint3('no data for variable; returning\n', 100)
      return

    #bad situation; if we have no data on this processor we don't want
    #it to mess up other processors min/max
    if DataMinMax[2] == False:
      DataMinMax[0] = sys.float_info.max
      DataMinMax[1] = -sys.float_info.max

    UseReduceToGetMinMaxPairs(DataMinMax)
    if DataForIds != None:
      UseReduceToSpreadValues(DataForIds)

  if PhactoriDbg():
    myDebugPrint3('min/max result: ' + str(DataMinMax) + '\n')
  if PhactoriDbg():
    myDebugPrint3('  sum result: ' + str(DataSumCnt[0]) + ' count result: ' + str(DataSumCnt[1]) + '\n')
  #myDebugPrint3('  average result: ' + str(DataSumCnt[0]/DataSumCnt[1]) + '\n')
  if PhactoriDbg():
    myDebugPrint3('  average result: ' + \
        str(DataSumCnt[0]/float(DataSumCnt[1])) + '\n')
  if PhactoriDbg():
    myDebugPrint3('DataForIds result: ' + str(DataForIds) + '\n')

  inVariableInfo.mStats.mMin = DataMinMax[0]
  inVariableInfo.mStats.mMax = DataMinMax[1]
  inVariableInfo.mStats.mSum = DataSumCnt[0]
  inVariableInfo.mStats.mCount = DataSumCnt[1]
  inVariableInfo.mStats.mStatsTestCounter = \
      gPipeAndViewsState.mFrameTagCounter

  if PhactoriDbg():
    myDebugPrint3("FindMinMaxSumCntFromData returning\n")


def SetPlotOverTimePointsFromData(ioPlotInfo):
  if PhactoriDbg(100):
    myDebugPrint3("SetPlotOverTimePointsFromData entered\n", 100)
  if PhactoriDbg(100):
    myDebugPrint3("variable is: " + str(ioPlotInfo.m_YAxisVariableInfo.mVariableName) + "\n", 100)

  inputcsData = ioPlotInfo.mInputOperation.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
  if inputcsData == None:
    if PhactoriDbg(100):
      myDebugPrint3('inputcsData is None, returning', 100)
    return

  if PhactoriDbg():
    myDebugPrint3('finding min/max for data for this processor at this timestep\n')
  DataMinMax = [0.0, 0.0, False]
  DataSumCnt = [0.0, 0]

  numIds = len(ioPlotInfo.m_IdPlotLineList)
  if numIds > 0:
    DataForIds = []
    for ii in range(0, numIds):
      DataForIds.append(0.0)
  else:
    DataForIds = None

  FindMinMaxSumCntFromData(inputcsData,
      ioPlotInfo.m_YAxisVariableInfo, ioPlotInfo.m_XAxisVariableInfo,
      DataMinMax, DataSumCnt, DataForIds, ioPlotInfo.m_IdPlotLineList)

  if DataSumCnt[1] == 0:
    if PhactoriDbg(100):
      myDebugPrint3('num elements in data is 0, returning', 100)
    return

  global gPipeAndViewsState
  plotPointTime = gPipeAndViewsState.CurrentDatadescription.GetTime()

  #append time to time column
  ioPlotInfo.m_TimeColumn.AppendPlotValue(plotPointTime)

  #append max value to plot
  plotDataVal = DataMinMax[1]
  ioPlotInfo.m_MaxPlotLine.AppendPlotValue(plotDataVal)

  #append min value to plot
  #plotDataVal = -0.1*(DataMinMax[0] + DataMinMax[1])
  plotDataVal = DataMinMax[0]
  ioPlotInfo.m_MinPlotLine.AppendPlotValue(plotDataVal)

  #append mean value to plot
  #ioPlotInfo.m_MeanPlotLine.AppendPlotValue(plotPointTime,
  #    DataSumCnt[0]/DataSumCnt[1])
  plotDataVal = DataSumCnt[0]/float(DataSumCnt[1])
  ioPlotInfo.m_MeanPlotLine.AppendPlotValue(DataSumCnt[0]/float(DataSumCnt[1]))

  if PhactoriDbg():
    myDebugPrint3("ioPlotInfo.m_vtkTable numcols B:" + str(ioPlotInfo.m_vtkTable.GetNumberOfColumns()) + "\n")
    myDebugPrint3("ioPlotInfo.m_vtkTable numrows B:" + str(ioPlotInfo.m_vtkTable.GetNumberOfRows()) + "\n")

  #append node or element values (by id) to plot lines for each id
  numIdValues = len(ioPlotInfo.m_IdPlotLineList)
  for ii in range(0,numIdValues):
   ioPlotInfo.m_IdPlotLineList[ii].m_PlotColumn.AppendPlotValue(DataForIds[ii])

  if PhactoriDbg():
    myDebugPrint3(' done updating ids\n')


def SetPlotPointsFromData(ioPlotInfo):
  if PhactoriDbg(100):
    myDebugPrint3("SetPlotPointsFromData entered\n", 100)

  inputcsData = ioPlotInfo.mInputOperation.GetPvFilter().\
          GetClientSideObject().GetOutputDataObject(0)
  if inputcsData == None:
    if PhactoriDbg(100):
      myDebugPrint3('inputcsData is None, returning', 100)
    return

  plotPointIndex = [0]
  #myDebugPrint3('  inputcsData: ' + str(inputcsData) + '\n')
  #if ioPlotInfo.mVariableType == 'element':
  #  cellOrPointFlag = 0
  #elif ioPlotInfo.mVariableType == 'node':
  #  cellOrPointFlag = 1
  #else:
  #  errStr = 'SetPlotPointsFromData error! can only handle element or node variable types for now'
  #  myDebugPrint3(errStr)
  #  raise Exception(errStr)
  #SetPlotPointsFromDataRecurse1(inputcsData, ioPlotInfo,
  #  cellOrPointFlag, plotPointIndex)
  SetPlotPointsFromDataRecurse1(inputcsData, ioPlotInfo, plotPointIndex)

  plotPts = ioPlotInfo.m_Points
  numPts = plotPts.GetNumberOfPoints()
  #myDebugPrint3(' need to update ids, get number of points: ' + str(numPts) + '\n')

  ioPlotInfo.m_Vertex.GetPointIds().SetNumberOfIds(numPts)

  for i in range(0,numPts):
    #Add the points to the line. The first value indicates
    #the order of the point on the line. The second value
    #is a reference to a point in a vtkPoints object. Depends
    #on the order that Points were added to vtkPoints object.
    #Note that this will not be associated with actual points
    #until it is added to a vtkPolyData object which holds a
    #vtkPoints object.
    #print str(i)
    ioPlotInfo.m_Vertex.GetPointIds().SetId(i, i)

  ioPlotInfo.m_PolyData.Reset()
  ioPlotInfo.m_PolyData.SetPoints(plotPts)
  ioPlotInfo.m_PolyData.InsertNextCell(ioPlotInfo.m_Vertex.GetCellType(), ioPlotInfo.m_Vertex.GetPointIds())

  #ioPlotInfo.m_PolyData.Modified()

  #using ScatterPlot() now
  #ioPlotInfo.m_producer.GetClientSideObject().SetOutput(None)
  #ioPlotInfo.m_producer.GetClientSideObject().SetOutput(ioPlotInfo.m_PolyData)

  #myDebugPrint3(' done updating ids\n')



  FindPlotXYZMinMax(plotPts, ioPlotInfo.m_xyzMinMaxTrkC)
  ScalePlotYForFit(plotPts, ioPlotInfo.m_xyzMinMaxTrkC,
      ioPlotInfo.mImageSettings.GetAspectRatioInsidePixelBorder())

  #find x and y min and max for data (current and over time)

  #myDebugPrint3('before--->')
  #PrintAllPointsA()

  #for ii in range(0,130):
  #  onePoint = plotPts.GetPoint(ii)
  #  xx = onePoint[0]
  #  yy = onePoint[1]
  #  #yy -= 10.0
  #  yy = 50.0
  #  zz = onePoint[2]
  #  plotPts.SetPoint(ii, xx, yy, zz)
  #myDebugPrint3('after--->')
  if PhactoriDbg():
    myDebugPrint3("points in plot after updating\n")
  #PrintAllPointsA(ioPlotInfo)
  if PhactoriDbg(100):
    myDebugPrint3("SetPlotPointsFromData returning\n", 100)

def CreateOnePlotOverTimeProducerC(ioPlotOverTime):
  if PhactoriDbg(100):
    myDebugPrint3('CreateOnePlotOverTimeProducerC entered (' + ioPlotOverTime.mName + ')\n', 100)

  global gPlotOverTimeMap

  #ioPlotOverTime.m_YAxisVariableName = ioinYAxisVariableName
  #ioPlotOverTime.m_CellOrPointOrGlobalFlag = inCellOrPointOrGlobalFlag

  UpdatePipelineWithCurrentTimeArgument(
      ioPlotOverTime.mInputOperation.GetPvFilter())

  savedActiveSource = GetActiveSource()

  if PhactoriDbg():
    myDebugPrint3('before active source: ' + str(savedActiveSource) + "\n")

  SetActiveSource(ioPlotOverTime.mInputOperation.GetPvFilter())

  ioPlotOverTime.m_vtkTable = vtk.vtkTable()

  ioPlotOverTime.m_vtkTable.AddColumn(ioPlotOverTime.m_TimeColumn.m_TableColumn)

  if ioPlotOverTime.mPlotMaximumFlag:
    ioPlotOverTime.m_vtkTable.AddColumn(ioPlotOverTime.m_MaxPlotLine.m_TableColumn)
  if ioPlotOverTime.mPlotMinimumFlag:
    ioPlotOverTime.m_vtkTable.AddColumn(ioPlotOverTime.m_MinPlotLine.m_TableColumn)
  if ioPlotOverTime.mPlotMeanFlag:
    ioPlotOverTime.m_vtkTable.AddColumn(ioPlotOverTime.m_MeanPlotLine.m_TableColumn)

  if PhactoriDbg():
    myDebugPrint3("ioPlotOverTime.m_vtkTable numcols A:" + str(ioPlotOverTime.m_vtkTable.GetNumberOfColumns()) + "\n")
    myDebugPrint3("ioPlotOverTime.m_vtkTable numrows A:" + str(ioPlotOverTime.m_vtkTable.GetNumberOfRows()) + "\n")

  if PhactoriDbg():
    myDebugPrint3("initializing id plot lines (" + str(len(ioPlotOverTime.m_IdPlotLineList)) + ")\n")
  for oneIdPlotLine in ioPlotOverTime.m_IdPlotLineList:
    if PhactoriDbg():
      myDebugPrint3("  doing line for id " + str(oneIdPlotLine.m_Id) + "\n")
    ioPlotOverTime.m_vtkTable.AddColumn(oneIdPlotLine.m_PlotColumn.m_TableColumn)

  producer = PVTrivialProducer()
  #producer.GetClientSideObject().SetOutput(ioPlotOverTime.m_PolyData)

  myDebugPrint3("ioPlotOverTime.m_vtkTable numcols C:" + str(ioPlotOverTime.m_vtkTable.GetNumberOfColumns()) + "\n")
  myDebugPrint3("ioPlotOverTime.m_vtkTable numrows C:" + str(ioPlotOverTime.m_vtkTable.GetNumberOfRows()) + "\n")
  producer.GetClientSideObject().SetOutput(ioPlotOverTime.m_vtkTable)

  ioPlotOverTime.m_producer = producer

  UpdatePipelineWithCurrentTimeArgument(
      ioPlotOverTime.mInputOperation.GetPvFilter())
  UpdatePipelineWithCurrentTimeArgument(ioPlotOverTime.m_producer)
  SetPlotOverTimePointsFromData(ioPlotOverTime)
  UpdatePipelineWithCurrentTimeArgument(ioPlotOverTime.m_producer)
  if PhactoriDbg():
    myDebugPrint3('done filter output 2nd time:')
  if PhactoriDbg():
    myDebugPrint3('After 2nd update pipeline in CreatePlotOverTimeProducer:')
  #PrintAllPointsA(ioPlotOverTime)

  #print '  active source after filter before reset: ' + str(GetActiveSource())
  SetActiveSource(savedActiveSource)
  if PhactoriDbg():
    myDebugPrint3('ioPlotOverTime.m_producer.GetDataInformation().GetBounds():\n' + str(ioPlotOverTime.m_producer.GetDataInformation().GetBounds()) + '\n')

  gPlotOverTimeMap[ioPlotOverTime.mName] = ioPlotOverTime
  if PhactoriDbg():
    myDebugPrint3('gPlotOverTimeMap has ' + str(len(gPlotOverTimeMap)) + ' items\n')
  count = 0
  for ii in gPlotOverTimeMap.values():
    if PhactoriDbg():
      myDebugPrint3(str(count) + ': ' + str(ii) + '\n')
    if PhactoriDbg():
      myDebugPrint3(str(count) + ': ' + ii.mName + '\n')
    count += 1

  AddFilterToFilterMap(ioPlotOverTime.mName, ioPlotOverTime.m_producer)

  if PhactoriDbg(100):
    myDebugPrint3('CreateOnePlotOverTimeProducerC returning\n', 100)


def CreateOneScatterPlotProducerC(ioScatterPlot):
  if PhactoriDbg(100):
    myDebugPrint3('CreateOneScatterPlotProducerC entered (' + ioScatterPlot.mName + ')\n', 100)

  global gScatterPlotMap

  #needs to change to source specified

  UpdatePipelineWithCurrentTimeArgument(
      ioScatterPlot.mInputOperation.GetPvFilter())

  savedActiveSource = GetActiveSource()
  if PhactoriDbg():
    myDebugPrint3('before active source: ' + str(savedActiveSource) + "\n")

  SetActiveSource(ioScatterPlot.mInputOperation.GetPvFilter())

  #create poly data to use
  ioScatterPlot.m_PolyData = vtk.vtkPolyData()
  #ioScatterPlot.m_PolyLine = vtk.vtkPolyLine()
  ioScatterPlot.m_Points = vtk.vtkPoints()
  ioScatterPlot.m_Vertex = vtk.vtkVertex()

  #numPts = 100
  #numPts = 1
  numPts = 0
  if PhactoriDbg():
    myDebugPrint3('numPts: ' + str(numPts) + "\n")

  for ii in range(0, numPts):
    xx = float(ii) / float(numPts)
    yy = (float(ii) * float(ii)) / (float(numPts) * float(numPts))
    zz = 0.0
    ioScatterPlot.m_Points.InsertPoint(ii, xx,yy,zz)
  ioScatterPlot.m_PolyData.SetPoints(ioScatterPlot.m_Points)

  #ioScatterPlot.m_PolyLine.GetPointIds().SetNumberOfIds(numPts)
  ioScatterPlot.m_Vertex.GetPointIds().SetNumberOfIds(numPts)

  for i in range(0,numPts):
    #Add the points to the line. The first value indicates
    #the order of the point on the line. The second value
    #is a reference to a point in a vtkPoints object. Depends
    #on the order that Points were added to vtkPoints object.
    #Note that this will not be associated with actual points
    #until it is added to a vtkPolyData object which holds a
    #vtkPoints object.
    #print str(i)
    #ioScatterPlot.m_PolyLine.GetPointIds().SetId(i, i)
    ioScatterPlot.m_Vertex.GetPointIds().SetId(i, i)

  ioScatterPlot.m_PolyData.Allocate(1, 1)
  ioScatterPlot.m_PolyData.InsertNextCell(ioScatterPlot.m_Vertex.GetCellType(),
      ioScatterPlot.m_Vertex.GetPointIds())

  #producer = PVTrivialProducer()
  #producer.GetClientSideObject().SetOutput(ioScatterPlot.m_PolyData)
  ioScatterPlot.m_MergeBlocks = MergeBlocks(Input=ioScatterPlot.mInputOperation.GetPvFilter())
  producer = ScatterPlot(Input=ioScatterPlot.m_MergeBlocks)

  if PhactoriDbg():
    csdoobj = producer.GetClientSideObject().GetOutputDataObject(0)
    myDebugPrint3("csdoobj: " + str(csdoobj) + "\n")
    myDebugPrint3("csdoobj.GetNumberOfPoints(): " + str(csdoobj.GetNumberOfPoints()) + "\n")
    myDebugPrint3("csdoobj.GetNumberOfCells(): " + str(csdoobj.GetNumberOfCells()) + "\n")

  ioScatterPlot.m_producer = producer

  UpdatePipelineWithCurrentTimeArgument(
      ioScatterPlot.mInputOperation.GetPvFilter())
  ioScatterPlot.ChooseDefaultVariableIfNecessary()
  UpdatePipelineWithCurrentTimeArgument(ioScatterPlot.m_producer)
  SetPlotPointsFromData(ioScatterPlot)
  UpdatePipelineWithCurrentTimeArgument(ioScatterPlot.m_producer)
  if PhactoriDbg():
    myDebugPrint3('done filter output 2nd time:')
  if PhactoriDbg():
    myDebugPrint3('After 2nd update pipeline in CreateScatterPlotProducer:')
  #PrintAllPointsA(ioScatterPlot)

  #print '  active source after filter before reset: ' + str(GetActiveSource())
  SetActiveSource(savedActiveSource)
  if PhactoriDbg():
    myDebugPrint3('ioScatterPlot.m_producer.GetDataInformation().GetBounds():\n' + str(ioScatterPlot.m_producer.GetDataInformation().GetBounds()) + '\n')

  gScatterPlotMap[ioScatterPlot.mName] = ioScatterPlot
  if PhactoriDbg():
    myDebugPrint3('gScatterPlotMap has ' + str(len(gScatterPlotMap)) + ' items\n')
  count = 0
  for ii in gScatterPlotMap.values():
    if PhactoriDbg():
      myDebugPrint3(str(count) + ': ' + str(ii) + '\n')
    if PhactoriDbg():
      myDebugPrint3(str(count) + ': ' + ii.mName + '\n')
    count += 1

  AddFilterToFilterMap(ioScatterPlot.mName, ioScatterPlot.m_producer)

  if PhactoriDbg(100):
    myDebugPrint3('CreateOneScatterPlotProducerC returning\n', 100)


def CreateOneScatterPlotViewC(ioScatterPlot):
  if PhactoriDbg(100):
    myDebugPrint3('CreateOneScatterPlotViewC entered (' + ioScatterPlot.mName + ')\n', 100)
  saveActiveSource = GetActiveSource()
  saveActiveView = GetActiveView()
  global currentPhactoriRenderViewInfo
  saveCurrentInfo = currentPhactoriRenderViewInfo

  SetActiveSource(ioScatterPlot.m_producer)

  #global gDefaultImageSizeX
  #global gDefaultImageSizeY
  #plot_view_setup_1 = { \
  #  'camera_setup': {'look_direction': [0.0, 0.0, -1.0] }, \
  #  'show_color_legend': False, \
  #  'show_data_cube_axes': True, \
  #  'show_orientation_axes': True, \
  #  'image_settings' : { \
  #    'basename': ioScatterPlot.mName, \
  #    'size': [gDefaultImageSizeX, gDefaultImageSizeY] \
  #    }, \
  #  }
  #
  #CreateViewSetFromPhactoriViewMapB(plot_view_setup_1)

  ioScatterPlot.mCamera.mName = 'dummyscatterplotcamera'
  ioScatterPlot.mCamera.mType = 'camera'
  ioScatterPlot.mCamera.mLookAtDistanceType = 'datasize relative'
  ioScatterPlot.mCamera.mLookDirection = [0.0, 0.0, -1.0]
  ioScatterPlot.mCamera.mLookAtDistance = 1.0

  if ioScatterPlot.mImageSettings.mUsingDefaultGeneratedImageBasename:
    imageBasename = ioScatterPlot.mName
    extraImagenameItem = "sctr."
  else:
    imageBasename = ioScatterPlot.mImageSettings.mImageBasename
    extraImagenameItem = ""

  SetUpOneParaViewRepresentationAndViewC(ioScatterPlot.mCamera,
    [0.0, 0.0, -1.0],
    inImagesetInfo = ioScatterPlot,
    inColorSettingsX = ioScatterPlot.mColorSettings,
    inMeshRenderControl = 'Surface With Edges',
    inShowDataCubeAxes = True,
    inShowDataCubeAxesInfo = ioScatterPlot.m_DataCubeAxesInfo,
    inShowOrientationAxes = False,
    inFixedColorRange = None,
    inIsPlotFlag = True,
    inRepresentationFilenameAddon = "",
    inLookDirectionFilenameAddon = extraImagenameItem,
    inPhactoriRepresentation = None)

  global gScatterPlotMap
  #myDebugPrint3('qqxxqq gScatterPlotMap has ' + str(len(gScatterPlotMap)) + ' items\n')
  count = 0
  for ii in gScatterPlotMap.values():
    if PhactoriDbg():
      myDebugPrint3(str(count) + ': ' + str(ii) + '\n')
    if PhactoriDbg():
      myDebugPrint3(str(count) + ': ' + ii.mName + '\n')
    count += 1

  if ioScatterPlot.mName in gScatterPlotMap:
    thePlot = gScatterPlotMap[ioScatterPlot.mName]
    thePlot.m_PhactoriRenderViewInfo = currentPhactoriRenderViewInfo

    #set up colors based .on preset and then set DataRepresentation colors
    #global gPlotColorInfoPreset
    #gPlotColorInfoPreset.CopyTo(thePlot.m_plotColorInfo)
    #thePlot.m_plotColorInfo.SetRepresentationAndViewColors(
    #  thePlot.m_PhactoriRenderViewInfo)
    #thePlot.mColorSettings.SetParaviewRvRepColors(
    #    thePlot.m_PhactoriRenderViewInfo.RenderView1,
    #    thePlot.m_PhactoriRenderViewInfo.DataRepresentation1)

    UpdatePlotViewLook(thePlot)
  else:
    if PhactoriDbg():
      myDebugPrint3('error 3 making scatter plot')

  SetActiveSource(saveActiveSource)
  SetActiveView(saveActiveView)
  currentPhactoriRenderViewInfo = saveCurrentInfo
  if PhactoriDbg(100):
    myDebugPrint3('CreateOneScatterPlotViewC returning\n', 100)

def UseReduceToSumDoubleAndIntPair(ioDoubleAndIntPair):
  import vtkParallelCorePython
  if PhactoriDbg(100):
    myDebugPrint3("UseReduceToSumDoubleAndIntPair entered\n", 100)
  if PhactoriDbg():
    myDebugPrint3("  local double: " + str(ioDoubleAndIntPair[0]) + " local int: " + str(ioDoubleAndIntPair[1]) + "\n")
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()

  localarray = vtk.vtkDoubleArray()
  localarray.SetNumberOfTuples(1)
  globalarray = vtk.vtkDoubleArray()
  globalarray.SetNumberOfTuples(1)
  localarray.SetValue(0, ioDoubleAndIntPair[0])

  #fixed bug here; 0 argument is mpi find max, while 2 argument is mpi sum
  #which is what we want here
  #globalController.AllReduce(localarray, globalarray, 0)
  globalController.AllReduce(localarray, globalarray, 2)
  ioDoubleAndIntPair[0] = globalarray.GetTuple1(0)

  localarray = vtk.vtkIntArray()
  localarray.SetNumberOfTuples(1)
  globalarray = vtk.vtkIntArray()
  globalarray.SetNumberOfTuples(1)
  localarray.SetValue(0, ioDoubleAndIntPair[1])
  globalController.AllReduce(localarray, globalarray, 2)
  ioDoubleAndIntPair[1] = globalarray.GetTuple1(0)

  if PhactoriDbg():
    myDebugPrint3("  reduced double: " + str(ioDoubleAndIntPair[0]) + " reduced int: " + str(ioDoubleAndIntPair[1]) + "\n")

def UseReduceToGetMaximumIntegers(ioListOfMaximumIntegers):
  import vtkParallelCorePython
  if PhactoriDbg(100):
    myDebugPrint3("UseReduceToGetMaximums entered\n", 100)
    myDebugPrint3("  local ioListOfMaximumIntegers: " +
      str(ioListOfMaximumIntegers) + "\n")
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()
  numMaxes = len(ioListOfMaximumIntegers)
  localarray = vtk.vtkIntArray()
  localarray.SetNumberOfTuples(numMaxes)
  for ii in range(0, numMaxes):
    localarray.SetValue(ii, int(ioListOfMaximumIntegers[ii]))
  globalarray = vtk.vtkIntArray()
  globalarray.SetNumberOfTuples(numMaxes)
  globalController.AllReduce(localarray, globalarray, 0)
  for ii in range(0, numMaxes):
    ioListOfMaximumIntegers[ii] = globalarray.GetTuple1(ii)

  if PhactoriDbg():
    myDebugPrint3("  global ioListOfMaximumIntegers: " +
      str(ioListOfMaximumIntegers) + "\n")

def UseReduceToGetMinMaxPairs(ioListOfMinMaxPairs):
  import vtkParallelCorePython
  if PhactoriDbg(100):
    myDebugPrint3("UseReduceToGetMinMaxPairs entered\n", 100)
  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()
  numPairs = len(ioListOfMinMaxPairs) / 2
  localarray = vtk.vtkDoubleArray()
  localarray.SetNumberOfTuples(numPairs*2)
  # we negate mins so we can do a single MPI_MAX reduce
  for ii in range(0, numPairs):
    ndx1 = ii*2
    localarray.SetValue(ndx1, -ioListOfMinMaxPairs[ndx1])
    ndx1 += 1
    localarray.SetValue(ndx1, ioListOfMinMaxPairs[ndx1])
  globalarray = vtk.vtkDoubleArray()
  globalarray.SetNumberOfTuples(numPairs*2)
  globalController.AllReduce(localarray, globalarray, 0)
  if PhactoriDbg():
    myDebugPrint3("  local minmax pairs: " + str(ioListOfMinMaxPairs) + "\n")
  for ii in range(0, numPairs):
    ndx1 = ii*2
    ioListOfMinMaxPairs[ndx1] = -globalarray.GetTuple1(ndx1)
    ndx1 += 1
    ioListOfMinMaxPairs[ndx1] = globalarray.GetTuple1(ndx1)
  if PhactoriDbg():
    myDebugPrint3("  global minmax pairs: " + str(ioListOfMinMaxPairs) + "\n")
  if PhactoriDbg(100):
    myDebugPrint3("UseReduceToGetMinMaxPairs returning\n", 100)

#set the 3d viewpoint to look at the xy scatter plot
def UpdatePlotViewLook(inPlotInfo):
  return

def CreateOnePlotOverTimeViewC(ioPlotOverTime):
  if PhactoriDbg(100):
    myDebugPrint3('CreateOnePlotOverTimeViewC entered (' + ioPlotOverTime.mName + ')\n', 100)
  saveActiveSource = GetActiveSource()
  saveActiveView = GetActiveView()
  global currentPhactoriRenderViewInfo
  saveCurrentInfo = currentPhactoriRenderViewInfo

  SetActiveSource(ioPlotOverTime.m_producer)

  #global gDefaultImageSizeX
  #global gDefaultImageSizeY
  #plot_view_setup_1 = { \
  #  'camera_setup': {'look_direction': [0.0, 0.0, -1.0] }, \
  #  'show_color_legend': False, \
  #  'show_data_cube_axes': True, \
  #  'show_orientation_axes': True, \
  #  'image_settings' : { \
  #    'basename': ioPlotOverTime.mName, \
  #    'size': [gDefaultImageSizeX, gDefaultImageSizeY], \
  #    }, \
  #  }

  #CreateViewSetFromPhactoriViewMapB(plot_view_setup_1)

  ioPlotOverTime.mCamera.mName = 'dummyscatterplotcamera'
  ioPlotOverTime.mCamera.mType = 'camera'
  ioPlotOverTime.mCamera.mLookAtDistanceType = 'datasize relative'
  ioPlotOverTime.mCamera.mLookDirection = [0.0, 0.0, -1.0]
  ioPlotOverTime.mCamera.mLookAtDistance = 1.0

  if ioPlotOverTime.mImageSettings.mUsingDefaultGeneratedImageBasename:
    imageBasename = ioPlotOverTime.mName
    extraImagenameItem = "plot."
  else:
    imageBasename = ioPlotOverTime.mImageSettings.mImageBasename
    extraImagenameItem = ""

  if PhactoriDbg():
    myDebugPrint3("  use specific y axis label: " + str(ioPlotOverTime.m_DataCubeAxesInfo.mYAxisInfo.mUseLabelFlag) + "\n")
  if ioPlotOverTime.m_DataCubeAxesInfo.mYAxisInfo.mUseLabelFlag:
    if PhactoriDbg():
      myDebugPrint3("      specific y axis label: " + str(ioPlotOverTime.m_DataCubeAxesInfo.mYAxisInfo.mAxisLabel) + "\n")

  SetUpOneParaViewRepresentationAndViewC(ioPlotOverTime.mCamera,
    [0.0, 0.0, -1.0],
    inImagesetInfo = ioPlotOverTime,
    inColorSettingsX = ioPlotOverTime.mColorSettings,
    inMeshRenderControl = 'Surface With Edges',
    inShowDataCubeAxes = True,
    inShowDataCubeAxesInfo = ioPlotOverTime.m_DataCubeAxesInfo,
    inShowOrientationAxes = False,
    inFixedColorRange = None,
    inIsPlotFlag = True,
    inRepresentationFilenameAddon = "",
    inLookDirectionFilenameAddon = extraImagenameItem,
    inPhactoriRepresentation = None)


  global gPlotOverTimeMap
  if PhactoriDbg():
    myDebugPrint3('qqxxqq gPlotOverTimeMap has ' + str(len(gPlotOverTimeMap)) + ' items\n')
  count = 0
  for ii in gPlotOverTimeMap.values():
    if PhactoriDbg():
      myDebugPrint3(str(count) + ': ' + str(ii) + '\n')
    if PhactoriDbg():
      myDebugPrint3(str(count) + ': ' + ii.mName + '\n')
    count += 1

  if ioPlotOverTime.mName in gPlotOverTimeMap:
    thePlot = gPlotOverTimeMap[ioPlotOverTime.mName]

    UpdatePlotViewLook(thePlot)
  else:
    if PhactoriDbg():
      myDebugPrint3('error 3 making plot over time')

  SetActiveSource(saveActiveSource)
  SetActiveView(saveActiveView)
  currentPhactoriRenderViewInfo = saveCurrentInfo
  if PhactoriDbg(100):
    myDebugPrint3('CreateOnePlotOverTimeViewC returning\n', 100)

global gPhactoriPipeRootMap
gPhactoriPipeRootMap = {}

def convertJsonUnicodeToStrings(inJson):
    """go through a dict/list/item which was loaded via json.load, and convert
       unicode strings to python strings"""
    #for python 2.6 or earlier, replace
    #  return {convertJsonUnicodeToStrings(key): convertJsonUnicodeToStrings(value) for key, value in inJson.iteritems()}
    #with
    #  return dict([(convert(key), convert(value)) for key, value in input.iteritems()])
    #(not tested, this routine was copied from web answer to question)
    if isinstance(inJson, dict):
        return {convertJsonUnicodeToStrings(key): convertJsonUnicodeToStrings(value) for key, value in inJson.iteritems()}
    elif isinstance(inJson, list):
        return [convertJsonUnicodeToStrings(element) for element in inJson]
    elif isinstance(inJson, unicode):
        return inJson.encode('utf-8')
    else:
        return inJson


def IssueErrorOrWarningThroughSierraIO(datadescription, outString, isErrorFlag):
  if PhactoriDbg():
    myDebugPrint3("IssueErrorOrWarningThroughSierraIO entered")

  if datadescription == None:
    errStr = "IssueErrorOrWarningThroughSierraIO no datadescription\n" + \
            outString
    raise Exception(errStr)
    return

  fd = datadescription.GetUserData()

  if fd == None:
    errStr = "IssueErrorOrWarningThroughSierraIO, no user data\n" + \
            outString
    if PhactoriDbg():
      myDebugPrint3(errStr)
    raise Exception(errStr)
    return

  #numArrays = fd.GetNumberOfArrays()
  #myDebugPrint3("abstract arrays: number of arrays: " + str(numArrays) + "\n")
  #for ii in range(numArrays):
  #  myDebugPrint3("array " + str(ii) + ":\n")
  #  oneArray = fd.GetAbstractArray(ii)
  #  myDebugPrint3("name: " + oneArray.GetName() + "\n")

  errorCodeArray = fd.GetAbstractArray("catalyst_sierra_error_codes")
  errorMessageArray = fd.GetAbstractArray("catalyst_sierra_error_messages")

  if errorCodeArray == None or errorMessageArray == None:
    if PhactoriDbg(999):
      myDebugPrint3("IssueErrorOrWarningThroughSierraIO::\n"
        "catalyst_sierra_error_codes and/or catalyst_sierra_error_messages\n"
        "were missing from UserData.  Item follows:\n" + outString, 999)
  else:
    if isErrorFlag:
      errorCode = 0
    else:
      errorCode = 1
    errorCodeArray.InsertNextValue(errorCode)
    errorMessageArray.InsertNextValue(outString)

  if PhactoriDbg():
    myDebugPrint3("IssueErrorOrWarningThroughSierraIO reached end")

def GetUserDataStringArrayAccountingForBypass(datadescription):
  if datadescription == None:
    if gBypassUserData:
      if PhactoriDbg():
        myDebugPrint3(
          "GetUserDataStringArrayAccountingForBypass returning " \
                  "bypass strings 1\n")
      return gUserDataBypassStrings
    else:
      if PhactoriDbg():
        myDebugPrint3(
          "GetUserDataStringArrayAccountingForBypass returning None 1\n")
      return None

  fd = datadescription.GetUserData()
  if fd == None:
    if gBypassUserData:
      if PhactoriDbg():
        myDebugPrint3(
          "GetUserDataStringArrayAccountingForBypass returning " \
                  "bypass strings 2\n")
      return gUserDataBypassStrings
    else:
      if PhactoriDbg():
        myDebugPrint3(
          "GetUserDataStringArrayAccountingForBypass returning None 2\n")
      return None

  sa = fd.GetAbstractArray(0)
  if sa == None or sa.GetNumberOfValues() < 1:
    if gBypassUserData:
      if PhactoriDbg():
        myDebugPrint3(
          "GetUserDataStringArrayAccountingForBypass returning " \
                  "bypass strings 3\n")
      return gUserDataBypassStrings
    else:
      if PhactoriDbg():
        myDebugPrint3(
          "GetUserDataStringArrayAccountingForBypass returning None 3\n")
      return None

  returnStringArray = ["","","","","","","","",""]

  for ii in range(0, sa.GetNumberOfValues()):
    returnStringArray[ii] = sa.GetValue(ii)

  if gBypassUserData:
    #sub in json string
    if PhactoriDbg():
      myDebugPrint3("gBypassUserData is true, user data is valid, subbing\n")
    returnStringArray[0] = gUserDataBypassStrings[0]
    #also add in missing user values if sa was short
    for jj in range(sa.GetNumberOfValues(), len(gUserDataBypassStrings)):
      returnStringArray[jj] = gUserDataBypassStrings[0]

  return returnStringArray


def GetViewMapCFromUserData(datadescription):
  """given a datadescription, get the json string view map from the user data
     and convert it into a python dict using the json library, and return
     that item.  Also determine the separator character and input deck
     filename"""
  import json

  #the following code forces a UserData item for internal testing
  #myDebugPrint3("GetViewMapCFromUserData entered\n")
  #newFd = vtk.vtkFieldData()
  #newStringArray = vtk.vtkStringArray()
  #xxViewMapCStr = '{ "camera blocks": { }, "representation blocks": { }, \
  #    "operation blocks": { }, "imageset blocks": {}, \
  #    "scatter plot blocks": { }, "plot over time blocks": { } }'
  #newStringArray.InsertNextValue(xxViewMapCStr)
  #x2str = '_'
  #x3str = 'cool_input_deck'
  #newStringArray.InsertNextValue(x2str)
  #newStringArray.InsertNextValue(x3str)
  #newFd.AddArray(newStringArray)

  #datadescription.SetUserData(newFd)

  global gBypassUserData

  sa = GetUserDataStringArrayAccountingForBypass(datadescription)

  testJsonString = sa[0]
  separatorCharacterString = sa[1]
  inputDeckFilename = sa[2]

  if PhactoriDbg():
    myDebugPrint3("gGetJsonViewMapCFromUserData string:\n" + \
        str(testJsonString) + "\n")
  if PhactoriDbg():
    myDebugPrint3("separator string: " + separatorCharacterString + \
        "\ninputDeckFilename: " + inputDeckFilename + "\n")
  #if PhactoriDbg():
    #myDebugPrint3("num strings2: " + str(sa.GetNumberOfValues()))
  if PhactoriDbg():
    myDebugPrint3("json.loads begin ---------------\n");
    myDebugPrint3(testJsonString);
    myDebugPrint3("\njson.loads end ---------------\n");

  if gBypassUserData == True:
    global gBypassUserDataJson
    returnViewMapC = gBypassUserDataJson
    #returnViewMapC = bccolli_controls.catalystSierraInputInJsonFormat
  else:
    returnViewMapC = json.loads(testJsonString)
    returnViewMapC = convertJsonUnicodeToStrings(returnViewMapC)

  SetDefaultImageBasename(inputDeckFilename)

  if PhactoriDbg():
    myDebugPrint3("GetViewMapCFromUserData returning\n")

  return returnViewMapC


def CreatePipelineFromDataDescription(datadescription):
  theSource = GetCurrentSource()

  theCellArrays = None
  theCellArrays = theSource.CellData
  thePointArrays = None
  thePointArrays = theSource.PointData

  theViewMapC = GetViewMapCFromUserData(datadescription)

  CreateViewSetFromPhactoriViewMapC(theViewMapC)

  #SetPlotView2StartColors(
  #  inBackgroundColor = [1.0, 1.0, 1.0],
  #  inEdgeColor = [0.0, 0.0, 0.5],
  #  inCubeAxesColor = [0.2, 0.2, 0.2],
  #  inDiffuseColor = [0.2, 0.2, 0.2],
  #  inAmbientColor = [0.2, 0.2, 0.2],
  #  inSelectionColor = [0.2, 0.2, 0.2],
  #  inBackfaceDiffuseColor = [0.2, 0.2, 0.2])

def parseDbDummyFname(dbDummyFname):
  """given user data item which is dummy database name used for remesh and
     restart purposes, pull out the base name (minus -s0002, -sXXXX, etc.)
     and return the base name and extension (extension in an empty string
     if appropriate)"""
  remeshRestartTagIndex = dbDummyFname.rfind("-s")
  if remeshRestartTagIndex == -1:
    baseId = dbDummyFname
    remeshRestartTag = ''
  else:
    baseId = dbDummyFname[0:remeshRestartTagIndex]
    remeshRestartTag = dbDummyFname[remeshRestartTagIndex:]
  return baseId, remeshRestartTag

def DoUserInteractionWithSimulationPausedIfEnabled():
  if PhactoriDbg():
    myDebugPrint3("DoUserInteractionWithSimulationPausedIfEnabled entered\n")
  global gPipeAndViewsState

  if gPipeAndViewsState.mInteractionRepeatPauseSim == False:
    if PhactoriDbg():
      myDebugPrint3("DoUserInteractionWithSimulationPausedIfEnabled returning 0\n")
    return 0

  updateWasDone = HandleUserInteractionIfEnabled(gPipeAndViewsState)
  if updateWasDone:
    returnVal = 2
  else:
    returnVal = 1
  if PhactoriDbg():
    myDebugPrint3("DoUserInteractionWithSimulationPausedIfEnabled returning " + \
        str(returnVal) + "\n")
  return returnVal

def HandleUserInteractionIfEnabled(ioPipeAndViewsState):
  """if the pipe (or imageset maybe) is in interactive mode, check for the
     trigger which will cause a read-and-update of the state.  Provides
     control to skip interaction when it is not enabled, the default case.
     Allows for a wait for an interaction trigger of N tries at M seconds
     per try, mainly to deal with situations where images are being
     rendered quickly.  Eventually may be altered to help deal with cases
     where we want to keep re-rendering with changed pipeline states
     until we are happy and then let the simulation return to operation"""
  if PhactoriDbg():
    myDebugPrint3("HandleUserInteractionIfEnabled entered\n")
  if ioPipeAndViewsState.mInteractionEnabled == False:
    return False

  if SmartGetLocalProcessId() != 0:
    return False
  #else:
    #do parallel stuff to receive info


  #check for trigger.  Only update if trigger file is incremented
  doUpdate = False
  triggerTestCount = 0

  while doUpdate == False and \
      triggerTestCount < ioPipeAndViewsState.mInteractionTriggerTries:
    if PhactoriDbg():
      myDebugPrint3("testing for trigger: " + str(triggerTestCount) + \
          " of " + str(ioPipeAndViewsState.mInteractionTriggerTries) + "\n")
    triggerTestCount += 1
    try:
      import json
      inFile = open('PhactoriInteractionTrigger.txt', 'rb')
      userTrigger = json.load(inFile)
      inFile.close()
      userTrigger = convertJsonUnicodeToStrings(userTrigger)
    except:
      if PhactoriDbg():
        myDebugPrint3("error dealing with PhactoriInteractionTrigger.txt\n")
      return

    if 'InteractionTriggerCounter' not in userTrigger:
      if PhactoriDbg():
        myDebugPrint3("no InteractionTriggerCounter in \
            PhactoriInteractionTrigger.txt\n")
    else:
      testTriggerCount = userTrigger['InteractionTriggerCounter']
      if testTriggerCount != ioPipeAndViewsState.mInteractionTriggerCounter:
        ioPipeAndViewsState.mInteractionTriggerCounter = testTriggerCount
        doUpdate = True

    if doUpdate == False:
      import time
      time.sleep(ioPipeAndViewsState.mInteractionTriggerSleep)

  if doUpdate:
    if PhactoriDbg():
      myDebugPrint3("update triggered, doing update\n")
    UpdatePipeAndViewsStateFromUserA(ioPipeAndViewsState,
        "PhactoriJsonSettings.txt")
    if 'InteractionState' in userTrigger:
      if userTrigger['InteractionState'] == 'RepeatWithPausedSimulation':
        ioPipeAndViewsState.mInteractionRepeatPauseSim = True
      else:
        ioPipeAndViewsState.mInteractionRepeatPauseSim = False
    else:
      ioPipeAndViewsState.mInteractionRepeatPauseSim = False

  if PhactoriDbg():
    myDebugPrint3("HandleUserInteractionIfEnabled returning\n")
  return doUpdate


class PhactoriCriteriaThreshold:
  """
  maintains a threshold value TT, and tracks when a min/max/mean/sum/count
  of a variable crosses the threshold; also keeps track of how many frames
  it has been since the crossing and if this is within a frame count setting
  and also tracks how many times it has been triggered and no longer triggers
  after a max has been reached
  """

  def __init__(self, inThresholdValue, inFramesAfterTrigger, inMaxTriggers):
    self.mThresholdValue = inThresholdValue
    self.mFramesAfterTrigger = inFramesAfterTrigger
    self.mMaxTriggers = inMaxTriggers
    self.mCrossingsSoFar = 0
    self.mLastTestValue = 0.0
    self.mObserveCallCount = 0
    self.mCountAtLastCrossing = 0
    self.mCountAtLastTrigger = 0
    self.mFramesSinceLastTrigger = 0
    self.mCurrentlyTriggered = False
    if PhactoriDbg(100):
      myDebugPrint3(
        "PhactoriCriteriaThreshold::__init__:\n"
        "self.mThresholdValue: " + str(self.mThresholdValue) + "\n"
        "self.mFramesAfterTrigger: " + str(self.mFramesAfterTrigger) + "\n"
        "self.mMaxTriggers: " + str(self.mMaxTriggers) + "\n", 100)


  def ObservePotentialCrossingValue(self, mTestValue):
    if PhactoriDbg(100):
      if self.mObserveCallCount != 0:
        if PhactoriDbg(100):
          myDebugPrint3("ObservePotentialCrossingValue:\n" +
            str(mTestValue) + " <> " + str(self.mThresholdValue) + " <> " + \
                    str(self.mLastTestValue) + "\n", 100)
      else:
        if PhactoriDbg(100):
          myDebugPrint3("ObservePotentialCrossingValue:\n" +
            str(mTestValue) + " <> " + str(self.mThresholdValue) + " <> " + \
                    str("first call") + "\n", 100)

    if self.mObserveCallCount != 0:
      if((self.mLastTestValue <= self.mThresholdValue) and \
              (mTestValue >= self.mThresholdValue)) or \
        ((self.mLastTestValue >= self.mThresholdValue) and \
                (mTestValue <= self.mThresholdValue)):
        #the test value crosses the threshold from the last test value
        self.mCrossingsSoFar += 1
        self.mCountAtLastCrossing = self.mObserveCallCount
        if PhactoriDbg(100):
            myDebugPrint3("crossing occurred: " + str(self.mCrossingsSoFar) + \
                    ", " + str(self.mCountAtLastCrossing) + "\n")
        if(self.mCrossingsSoFar <= self.mMaxTriggers):
          self.mCountAtLastTrigger = self.mObserveCallCount
          if PhactoriDbg(100):
              myDebugPrint3("self.mCountAtLastTrigger changed to " +
                      str(self.mCountAtLastTrigger) + "\n")
        else:
          if PhactoriDbg(100):
              myDebugPrint3("no trigger due to self.mMaxTriggers")

      self.mCurrentlyTriggered = False
      if self.mCrossingsSoFar != 0:
        if (self.mFramesAfterTrigger <= 0) or \
          (self.mObserveCallCount - self.mCountAtLastTrigger) < \
                self.mFramesAfterTrigger:
          self.mCurrentlyTriggered = True

      if PhactoriDbg(100):
          myDebugPrint3("self.mCurrentlyTriggered now " + \
                  str(self.mCurrentlyTriggered) + "\n")

    self.mLastTestValue = mTestValue
    self.mObserveCallCount += 1
    return

  def CurrentlyTriggered(self):
    return self.mCurrentlyTriggered

class PhactoriImagesetOnOffFilterCriteria:
  """
  this class represents one test used to determine if images are or are
  not shown.  The general idea is that the test is whether a variable
  minimum/maximum/mean is above/below/between a certain value or values.
  The variable can be scalar, vector magnitude, vector component, or
  tensor component, or global, and can be cell/element or point/node
  """

  def __init__(self):
    self.mName = ""
    self.mVariableInfo = PhactoriVariableInfo()
    self.mRange = [0.0, 1.0]
    #min/mean/max
    self.mFunctionType = "maximum"
    self.mInputOperationName = None
    self.mInputOperation = None
    #default number of frames after trigger is 5
    self.mNumberOfFramesAfterTrigger = 5
    self.mHasBeenTriggered = False
    self.mCallbackTagAtFirstTrigger = -1
    self.mThresholdList = []
    self.mMaxTriggersPerThreshold = 1
    self.mFrameTagCounter = -1
    self.mReturnValAtCurrentFrameTagCounter = False

  def TestForTruth(self, ioPipeAndViewsState):
    """returns true if this criteria indicates we should render a frame, false
       otherwise"""

    if PhactoriDbg(100):
      myDebugPrint3(
        "PhactoriImagesetOnOffFilterCriteria::TestForTruth entered", 100)

    if self.mFrameTagCounter == ioPipeAndViewsState.mFrameTagCounter:
      if PhactoriDbg(100):
        myDebugPrint3(
          "additional call, same mFrameTagCounter, just return value", 100)
      return self.mReturnValAtCurrentFrameTagCounter
    self.mFrameTagCounter = ioPipeAndViewsState.mFrameTagCounter

    #if this has been triggered, it remains true no matter what variables do.
    #if there is a limit to the number of frames this trigger is true, it
    #will be false after that many frames no matter what
    #if self.mHasBeenTriggered:
    #  if self.mNumberOfFramesAfterTrigger < 0:
    #    #no frame limit; always true now
    #    return True
    #  else:
    #    countSoFar = ioPipeAndViewsState.mFrameTagCounter - \
    #        self.mCallbackTagAtFirstTrigger
    #    if PhactoriDbg(100):
    #      myDebugPrint3(
    #        "mFrameTagCounter: " + \
    #          str(ioPipeAndViewsState.mFrameTagCounter) + "\n"
    #        "mCallbackTagAtFirstTrigger: " + \
    #          str(self.mCallbackTagAtFirstTrigger) + "\n"
    #        "countSoFar: " + str(countSoFar) + "\n"
    #        "mNumberOfFramesAfterTrigger: " + \
    #          str(self.mNumberOfFramesAfterTrigger) + "\n",
    #        100)
    #    if countSoFar >= self.mNumberOfFramesAfterTrigger:
    #      #reached frame limit
    #      return False
    #    else:
    #      #haven't reached frame limit, still true
    #      return True

    if self.mInputOperation == None:
      if self.mInputOperationName == None:
        self.mInputOperation = ioPipeAndViewsState.mIncomingDefaultOperation
      else:
        if self.mInputOperationName not in ioPipeAndViewsState.mOperationBlocks:
          errStr = 'error! in PhactoriImagesetOnOffFilterCriteria::TestForTruth calls for nonexistent input operation with name ' + str(self.mInputOperationName) + '\n'
          if PhactoriDbg():
            myDebugPrint3(errStr)
          raise Exception(errStr)
        else:
          self.mInputOperation = \
              ioPipeAndViewsState.mOperationBlocks[self.mInputOperationName]

    theParaViewFilter = self.mInputOperation.GetPvFilter()
    UpdatePipelineWithCurrentTimeArgument(theParaViewFilter)

    ###find min/max/avg/sum/count
    inputcsData = theParaViewFilter.GetClientSideObject().\
        GetOutputDataObject(0)
    if inputcsData == None:
      if PhactoriDbg(100):
        myDebugPrint3('inputcsData is None, returning false', 100)
      self.mReturnValAtCurrentFrameTagCounter = False
      return False

    if PhactoriDbg():
      myDebugPrint3('finding min/max data for this processor at this timestep\n')
    DataMinMax = [0.0, 0.0, False]
    DataSumCnt = [0.0, 0]

    FindMinMaxSumCntFromData(inputcsData, self.mVariableInfo, None,
        DataMinMax, DataSumCnt, None, None)

    #just return false if there is no min/max found
    #if (DataMinMax[0] == 0.0) and (DataMinMax[1] == 0.0):
    if (DataSumCnt[1] == 0):
      if PhactoriDbg(100):
        myDebugPrint3('no data value for min/max found, returning', 100)
      if PhactoriDbg(100):
        myDebugPrint3("PhactoriImagesetOnOffFilterCriteria returning false",
            100)
      self.mReturnValAtCurrentFrameTagCounter = False
      return False

    if self.mFunctionType == "maximum":
      mTestValue = DataMinMax[1]
      if PhactoriDbg():
        myDebugPrint3("test is against max, " + str(mTestValue) + "\n")
    elif self.mFunctionType == "minimum":
      mTestValue = DataMinMax[0]
      if PhactoriDbg():
        myDebugPrint3("test is against min, " + str(mTestValue) + "\n")
    elif self.mFunctionType == "mean":
      mTestValue = float(DataSumCnt[1])/float(DataSumCnt[0])
      if PhactoriDbg():
        myDebugPrint3("test is against mean, " + str(mTestValue) + "\n")
    elif self.mFunctionType == "sum":
      mTestValue = DataSumCnt[0]
      if PhactoriDbg():
        myDebugPrint3("test is against sum, " + str(mTestValue) + "\n")
    elif self.mFunctionType == "count":
      mTestValue = DataSumCnt[1]
      if PhactoriDbg():
        myDebugPrint3("test is against count, " + str(mTestValue) + "\n")
    else:
      myDebugPrint3AndException("PhactoriOnOffCriteria::TestForTruth"
      "bad mFunctionType: " + str(self.mFunctionType) + "\n")

    #let each threshold see if this value crossed from the previous value,
    #and update its internal frame count as appropriate
    for oneThrsld in self.mThresholdList:
      oneThrsld.ObservePotentialCrossingValue(mTestValue)

    #if any threshold was just triggered or was triggered and is within
    #the frame count limit, issue a 'True'
    for oneThrsld in self.mThresholdList:
      if oneThrsld.CurrentlyTriggered():
        self.mReturnValAtCurrentFrameTagCounter = True
        return True

    self.mReturnValAtCurrentFrameTagCounter = False
    return False

    #if self.mFunctionType == "maximum":
    #  #possible issue, but probably not; if min/max is 0.0 we assume no
    #  #data existed
    #  #DataMinMax[2] == False: don't use this, it's not mpi shared
    #  if (DataMinMax[0] == 0.0) and (DataMinMax[1] == 0.0):
    #    if PhactoriDbg(100):
    #      myDebugPrint3('no data value for max found, returning', 100)
    #    if PhactoriDbg(
    #          100):
    #      myDebugPrint3("PhactoriImagesetOnOffFilterCriteria returning false",
    #          100)
    #    return False
    #  else:
    #    mTestValue = DataMinMax[1]
    #    if PhactoriDbg():
    #      myDebugPrint3("test is against max, " + str(mTestValue) + "\n")
    #elif self.mFunctionType == "minimum":
    #  #possible issue, but probably not; if min/max is 0.0 we assume no
    #  #data existed
    #  #DataMinMax[2] == False: don't use this, it's not mpi shared
    #  if (DataMinMax[0] == 0.0) and (DataMinMax[1] == 0.0):
    #    if PhactoriDbg(100):
    #      myDebugPrint3('no data value for min found, returning', 100)
    #    if PhactoriDbg(
    #          100):
    #      myDebugPrint3("PhactoriImagesetOnOffFilterCriteria returning false",
    #          100)
    #    return False
    #  else:
    #    mTestValue = DataMinMax[0]
    #    if PhactoriDbg():
    #      myDebugPrint3("test is against min, " + str(mTestValue) + "\n")
    #elif self.mFunctionType == "mean":
    #  if DataSumCnt[1] == 0:
    #    if PhactoriDbg(100):
    #      myDebugPrint3('no data values for mean found, returning false', 100)
    #    return False
    #  else:
    #    mTestValue = float(DataSumCnt[1])/float(DataSumCnt[0])
    #    if PhactoriDbg():
    #      myDebugPrint3("test is against mean, " + str(mTestValue) + "\n")
    #elif self.mFunctionType == "sum":
    #  mTestValue = DataSumCnt[0]
    #  if PhactoriDbg():
    #    myDebugPrint3("test is against sum, " + str(mTestValue) + "\n")
    #elif self.mFunctionType == "count":
    #  mTestValue = DataSumCnt[1]
    #  if PhactoriDbg():
    #    myDebugPrint3("test is against count, " + str(mTestValue) + "\n")

    #if PhactoriDbg():
    #  myDebugPrint3("test range: " + str(self.mRange) + "\n")
    #if mTestValue < self.mRange[0]:
    #  if PhactoriDbg(100):
    #    myDebugPrint3('below range, returning false', 100)
    #  return False
    #if mTestValue > self.mRange[1]:
    #  if PhactoriDbg(100):
    #    myDebugPrint3('above range, returning false', 100)
    #  return False

    #if PhactoriDbg(100):
    #  myDebugPrint3('in range, returning true', 100)

    #if PhactoriDbg(100):
    #  myDebugPrint3("cp 6000")


    #self.mHasBeenTriggered = True
    #self.mCallbackTagAtFirstTrigger = ioPipeAndViewsState.mFrameTagCounter
    #return True

  def ParseParametersFromJson(self, inJson):
    if PhactoriDbg(100):
      myDebugPrint3('PhactoriImagesetOnOffFilterCriteria::'
        'ParseParametersFromJson entered\n', 100)

    foundVariableFlag = self.mVariableInfo.\
        ParseVariableNameAndVectorOrTensorComponent(
            inJson, 'variable ')

    if foundVariableFlag == False:
      errStr = """error!  inJson has no variable info in
          PhactoriImagesetOnOffFilterCriteria.ParseParametersFromJson\n"""
      if PhactoriDbg():
        myDebugPrint3(errStr)
      raise Exception(errStr)

    if 'input' in inJson:
      self.mInputOperationName = inJson['input']
    else:
      noticeStr = \
          'notice!  inJson has no input key, using default pipeline input\n'
      self.mInputOperationName = None
      if PhactoriDbg():
        myDebugPrint3(noticeStr)

    if PhactoriDbg():
      myDebugPrint3("variable name: " + \
          str(self.mVariableInfo.mVariableName) + "\n")

    if 'type' in inJson:
      theType = inJson['type']
      #bug in parsing right now requires different if test
      #if theType == 'triggerthresholds':
      if theType == 'triggerthresholds' or theType == 'onoffcriteriablock':
        if 'threshold variable function' in inJson:
          self.mFunctionType = inJson['threshold variable function']
        if 'frames after trigger' in inJson:
          self.mNumberOfFramesAfterTrigger = \
              inJson['frames after trigger']
        if 'maximum triggers per threshold' in inJson:
          self.mMaxTriggersPerThreshold = \
              inJson['maximum triggers per threshold']
        for oneThreshold in inJson['threshold values']:
          newThrshld = PhactoriCriteriaThreshold(oneThreshold,
                  self.mNumberOfFramesAfterTrigger,
                  self.mMaxTriggersPerThreshold)
          self.mThresholdList.append(newThrshld)

        if PhactoriDbg():
          myDebugPrint3("parsed min/max triggerthresholds\n")
        return

    if 'keep between' in inJson:
      self.mRange = inJson['keep between']
      self.mFunctionType = "maximum"
    elif 'keep below' in inJson:
      rangeMin = -sys.float_info.max
      rangeMax = inJson['keep below']
      self.mRange = [rangeMin, rangeMax]
      self.mFunctionType = "minimum"
    elif 'keep above' in inJson:
      self.mFunctionType = "maximum"
      rangeMin = inJson['keep above']
      rangeMax = sys.float_info.max
      self.mRange = [rangeMin, rangeMax]
    elif 'trigger when maximum value is at or above' in inJson:
      self.mFunctionType = "maximum"
      rangeMin = inJson['trigger when maximum value is at or above']
      rangeMax = sys.float_info.max
      self.mRange = [rangeMin, rangeMax]
    elif 'trigger on maximum value' in inJson:
      self.mFunctionType = "maximum"
      rangeMin = inJson['trigger on maximum value']
      rangeMax = sys.float_info.max
      self.mRange = [rangeMin, rangeMax]
    elif 'trigger on minimum value' in inJson:
      self.mFunctionType = "minimum"
      rangeMin = -sys.float_info.max
      rangeMax = inJson['trigger on minimum value']
      self.mRange = [rangeMin, rangeMax]
    else:
      if PhactoriDbg():
        myDebugPrint3("  no keep between/above/below, using keep above 0.0\n")
      self.mFunctionType = "maximum"
      rangeMin = 0.0
      rangeMax = sys.float_info.max
      self.mRange = [rangeMin, rangeMax]

    if 'number of frames after trigger' in inJson:
      self.mNumberOfFramesAfterTrigger = \
          inJson['number of frames after trigger']

    if 'frames after trigger' in inJson:
      self.mNumberOfFramesAfterTrigger = \
          inJson['frames after trigger']

    if PhactoriDbg():
      myDebugPrint3("range: " + str(self.mRange) + "\n")


class PhactoriImagesetOnOffFilter:
  """
  this class manages the system whereby users can selectively
  choose to render or not render images based on simulation data.  The first
  simple capability allows users to begin rendering after some simulation
  data criteria is met (e.g. maximum of element variable X is above Y), and
  stop rendering after some criteria is met"""

  def __init__(self):
    self.mEnabled = False
    self.mLastStartImagesFlag = False
    self.mStartCriteriaList = []
    #self.mStopCriteriaList = []

  def AddStartCriteria(self, inNewCriteria):
    self.mEnabled = True
    self.mStartCriteriaList.append(inNewCriteria)

  #def AddStopCriteria(self, inNewCriteria):
  #  self.mEnabled = True
  #  self.mStopCriteriaList.append(inNewCriteria)

  def TestDrawImagesThisCallback(self, ioPipeAndViewsState):
    """returns True if images should be rendered, False if they should not"""

    if PhactoriDbg(
          100):
      myDebugPrint3(
          'PhactoriImagesetOnOffFilter::TestDrawImagesThisCallback entered\n',
          100)

    #criteriaIndex = 0
    #for oneCriteria in self.mStartCriteriaList:
    #  criteriaIndex += 1
    #  myDebugPrint3("trying criteria: " + str(criteriaIndex))
    #  if oneCriteria.TestForTruth(ioPipeAndViewsState):
    #    myDebugPrint3("criteria returned true")
    #  else:
    #    myDebugPrint3("criteria returned true")
    #return True

    if self.mEnabled == False:
      if PhactoriDbg():
        myDebugPrint3("not enabled, return True\n")
      return True

    #if images had already been started, keep doing them until a stop
    #criteria is hit
    #if self.mLastStartImagesFlag == True:
    #  if PhactoriDbg():
    #    myDebugPrint3("  returned true last time\n")
    #  startImages = True
    #else:
    if True:
      #if no startup criteria are present, startup is assumed to be true
      if len(self.mStartCriteriaList) == 0:
        if PhactoriDbg():
          myDebugPrint3("  no criteria, default true\n")
        startImages = True
      else:
        if PhactoriDbg():
          myDebugPrint3("  criteria exists, default false\n")
        startImages = False
        for oneCriteria in self.mStartCriteriaList:
          if oneCriteria.TestForTruth(ioPipeAndViewsState):
            if PhactoriDbg():
              myDebugPrint3("  criteria returned true\n")
            startImages = True
            #at least one startup criteria is true
            break

    #if images are going (even if they just started), check the stop criteria
    #to see if they need stopping
    #if startImages:
    #  if len(self.mStopCriteriaList) != 0:
    #    for oneCriteria in self.mStopCriteriaList:
    #      if oneCriteria.TestForTruth(ioPipeAndViewsState):
    #        startImages = False
    #        #at least one stop criteria is true
    #        break

    self.mLastStartImagesFlag = startImages
    if PhactoriDbg(100):
      myDebugPrint3('returning with value: ' + str(startImages) + '\n', 100)
    return startImages

def ParseOneImageStartStopFilterFromViewMap(inOperationParamsJson):
  newFilterCriteria = PhactoriImagesetOnOffFilterCriteria()
  newFilterCriteria.ParseParametersFromJson(inOperationParamsJson)
  global gPipeAndViewsState
  gPipeAndViewsState.mImagesetOnOffFilter.AddStartCriteria(newFilterCriteria)


def WriteOutImagesTest(datadescription, coprocessor):
  global gPipeAndViewsState
  return gPipeAndViewsState.mImagesetOnOffFilter.TestDrawImagesThisCallback(gPipeAndViewsState)


def InitializePerPipeRoot(datadescription, inCoprocessor):
  if PhactoriDbg():
    myDebugPrint3("InitializePerPipeRoot entered\n")

  sa = GetUserDataStringArrayAccountingForBypass(datadescription)

  testJsonString = sa[0]
  separatorCharacterString = sa[1]
  inputDeckFilename = sa[2]
  extraFileString = sa[3]
  sa4 = sa[4]
  sa5 = sa[5]
  dbDummyFname = sa[6]
  if(len(sa) > 7):
    defaultDirectory = sa[7]
  else:
    defaultDirectory = ""
  if(len(sa) > 8):
    catalyst_script_extra_file = sa[8]
  else:
    catalyst_script_extra_file = ""

  outputResultBaseId, remeshRestartTag = parseDbDummyFname(dbDummyFname)
  if PhactoriDbg():
    myDebugPrint3("  extra data list:\n")
  if PhactoriDbg():
    myDebugPrint3("  separatorCharacterString: ->" + separatorCharacterString + "<-\n")
  if PhactoriDbg():
    myDebugPrint3("  inputDeckFilename: ->" + inputDeckFilename + "<-\n")
  if PhactoriDbg():
    myDebugPrint3("  extraFileString: ->" + extraFileString + "<-\n")
  if PhactoriDbg():
    myDebugPrint3("  sa4: ->" + sa4 + "<-\n")
  if PhactoriDbg():
    myDebugPrint3("  sa5: ->" + sa5 + "<-\n")
  if PhactoriDbg():
    myDebugPrint3("  dbDummyFname: ->" + dbDummyFname + "<-\n")
  if PhactoriDbg():
    myDebugPrint3("  outputResultBaseId: ->" + outputResultBaseId + "<-\n")
  if PhactoriDbg():
    myDebugPrint3("  remeshRestartTag: ->" + remeshRestartTag + "<-\n")
  if PhactoriDbg():
    myDebugPrint3("  defaultDirectory: ->" + defaultDirectory + "<-\n")
  if PhactoriDbg():
    myDebugPrint3("  catalyst_script_extra_file: ->" + catalyst_script_extra_file + "<-\n")

  global gPipeAndViewsState

  if outputResultBaseId in gPhactoriPipeRootMap:
    if PhactoriDbg():
      myDebugPrint3("  I've seen this outputResults block id\n")
    oneRoot = gPhactoriPipeRootMap[outputResultBaseId]
    oneRoot.CurrentDatadescription = datadescription
    gPipeAndViewsState = oneRoot
    if remeshRestartTag != oneRoot.mRemeshRestartTag:
      oneRoot.mRemeshRestartTag = remeshRestartTag
      oneRoot.mRemeshCount += 1
    HandleUserInteractionIfEnabled(oneRoot)
  else:
    if PhactoriDbg():
      myDebugPrint3("  first time for this outputResults block id\n")
    SetUpCoProcessor(inCoprocessor)
    newRoot = PhactoriPipeAndViewsState()
    newRoot.mDefaultBasedirectory = defaultDirectory
    newRoot.CurrentDatadescription = datadescription
    newRoot.mOutputResultsBlockId = outputResultBaseId
    newRoot.mCurrentDatabaseDummyFileName = dbDummyFname
    newRoot.mRemeshRestartTag = remeshRestartTag
    newRoot.mSeparatorString = separatorCharacterString
    newRoot.mOutputResultsBlockCountId = len(gPhactoriPipeRootMap)
    newRoot.mBlockIdForRestart = "restart." + str(len(gPhactoriPipeRootMap))
    gPhactoriPipeRootMap[outputResultBaseId] = newRoot
    gPipeAndViewsState = newRoot
    newRoot.mProducer = inCoprocessor.CreateProducer( datadescription, "input" )
    CreatePipelineFromDataDescription(datadescription)
    HandleRestartUpdateForOutputResultsBlock(newRoot)

  gPipeAndViewsState.mCallbackDateTime = datetime.datetime.now()

  if PhactoriDbg():
    myDebugPrint3("InitializePerPipeRoot returning\n")


global gColorLegendCollection
gColorLegendCollection = {
  "Default":{
    "ColorSpace":"Diverging",
    "NanColor":[0.25, 0.0, 0.0],
    "RGBPoints":[
      0.0, 0.23000000000000001, 0.29899999999999999, 0.754,
      1.0, 0.70599999999999996, 0.016, 0.14999999999999999
    ]
  },
  "Cool to Warm":{
    "ColorSpace":"Diverging",
    "NanColor":[0.247058823529, 0.0, 0.0],
    "RGBPoints":[
      0.0, 0.23137254902, 0.298039215686, 0.752941176471,
      0.5, 0.865, 0.865, 0.865,
      1.0, 0.705882352941, 0.0156862745098, 0.149019607843
    ]
  },
  "Blue to Red Rainbow":{
    "ColorSpace":"HSV",
    "NanColor":[0.498039215686, 0.498039215686, 0.498039215686],
    "RGBPoints":[
      0.0, 0.0, 0.0, 1.0,
      1.0, 1.0, 0.0, 0.0
    ]
  },
  "Red to Blue Rainbow":{
    "ColorSpace":"HSV",
    "NanColor":[0.498039215686, 0.498039215686, 0.498039215686],
    "RGBPoints":[
      0.0, 1.0, 0.0, 0.0,
      1.0, 0.0, 0.0, 1.0
    ]
  },
  "Grayscale":{
    "ColorSpace":"RGB",
    "NanColor":[1.0, 0.0, 0.0],
    "RGBPoints":[
      0.0, 0.0, 0.0, 0.0,
      1.0, 1.0, 1.0, 1.0
    ]
  },
  "X Ray":{
    "ColorSpace":"RGB",
    "NanColor":[1.0, 0.0, 0.0],
    "RGBPoints":[
      0.0, 1.0, 1.0, 1.0,
      1.0, 0.0, 0.0, 0.0
    ]
  },
  "Blue to Yellow":{
    "ColorSpace":"RGB",
    "NanColor":[1.0, 0.0, 0.0],
    "RGBPoints":[
      0.0, 0.0392156862745, 0.0392156862745, 0.949019607843,
      1.0, 0.949019607843, 0.949019607843, 0.0392156862745
    ]
  },
  #"Black-Body Radiation":{}
  "Black Body Radiation":{
    "ColorSpace":"RGB",
    "NanColor":[0.0, 0.498039215686, 1.0],
    "RGBPoints":[
      0.0, 0.0, 0.0, 0.0,
      0.4, 0.901960784314, 0.0, 0.0,
      0.8, 0.901960784314, 0.901960784314, 0.0,
      1.0, 1.0, 1.0, 1.0
    ]
  },
  "CIELab Blue to Red":{
    "ColorSpace":"Lab",
    "NanColor":[1.0, 1.0, 0.0],
    "RGBPoints":[
      0.0, 0.0, 0.6, 0.749019607843,
      1.0, 0.76862745098, 0.466666666667, 0.341176470588
    ]
  },
  "Black, Blue and White":{
    "ColorSpace":"RGB",
    "NanColor":[1.0, 1.0, 0.0],
    "RGBPoints":[
      0.0, 0.0, 0.0, 0.0,
      0.333, 0.0, 0.0, 0.501960784314,
      0.666, 0.0, 0.501960784314, 1.0,
      1.0, 1.0, 1.0, 1.0
    ]
  },
  "Black, Orange and White":{
    "ColorSpace":"RGB",
    "NanColor":[1.0, 1.0, 0.0],
    "RGBPoints":[
      0.0, 0.0, 0.0, 0.0,
      0.333, 0.501960784314, 0.0, 0.0,
      0.666, 1.0, 0.501960784314, 0.0,
      1.0, 1.0, 1.0, 1.0
    ]
  },
  "Cold and Hot":{
    "ColorSpace":"RGB",
    "NanColor":[1.0, 1.0, 0.0],
    "RGBPoints":[
      0.0, 0.0, 1.0, 1.0,
      0.45, 0.0, 0.0, 1.0,
      0.5, 0.0, 0.0, 0.501960784314,
      0.55, 1.0, 0.0, 0.0,
      1.0, 1.0, 1.0, 0.0
    ]
  },
  "Rainbow Desaturated":{
    "ColorSpace":"RGB",
    "NanColor":[1.0, 1.0, 0.0],
    "RGBPoints":[
      0.0, 0.278431372549, 0.278431372549, 0.858823529412,
      0.143, 0.0, 0.0, 0.360784313725,
      0.285, 0.0, 1.0, 1.0,
      0.429, 0.0, 0.501960784314, 0.0,
      0.571, 1.0, 1.0, 0.0,
      0.714, 1.0, 0.380392156863, 0.0,
      0.857, 0.419607843137, 0.0, 0.0,
      1.0, 0.878431372549, 0.301960784314, 0.301960784314
    ]
  },
  "Rainbow Blended White":{
    "ColorSpace":"RGB",
    "NanColor":[1.0, 1.0, 0.0],
    "RGBPoints":[
      0.0, 1.0, 1.0, 1.0,
      0.17, 0.0, 0.0, 1.0,
      0.34, 0.0, 1.0, 1.0,
      0.5, 0.0, 1.0, 0.0,
      0.67, 1.0, 1.0, 0.0,
      0.84, 1.0, 0.0, 0.0,
      1.0, 0.878431372549, 0.0, 1.0
    ]
  },
  "Rainbow Blended Grey":{
    "ColorSpace":"RGB",
    "NanColor":[1.0, 1.0, 0.0],
    "RGBPoints":[
      0.0, 0.317647058824, 0.341176470588, 0.43137254902,
      0.17, 0.0, 0.0, 1.0,
      0.34, 0.0, 1.0, 1.0,
      0.5, 0.0, 1.0, 0.0,
      0.67, 1.0, 1.0, 0.0,
      0.84, 1.0, 0.0, 0.0,
      1.0, 0.878431372549, 0.0, 1.0
    ]
  },
  #"Cool to Warm with Top":{
  #  "ColorSpace":"Diverging",
  #  "NanColor":[0.247058823529, 0.0, 0.0],
  #  "RGBPoints":[
  #    0.0, 0.23137254902, 0.298039215686, 0.752941176471,
  #    0.5, 0.865, 0.865, 0.865,
  #    0.98, 0.705882352941, 0.0156862745098, 0.149019607843,
  #    1.0, 0.0, 1.0, 1.0,
  #  ]
  #},
  "Cool to Warm with Top":{
    "ColorSpace":"Diverging",
    #"ColorSpace":"RGB",
    "NanColor":[0.247058823529, 0.0, 0.0],
    "RGBPoints":[
      0.0, 0.23137254902, 0.298039215686, 0.752941176471,
      0.5, 0.865, 0.865, 0.865,
      0.799, 0.769528, 0.355412, 0.435412,
      0.800, 0.0, 1.0, 1.0,
      0.814, 0.0, 1.0, 1.0,
      0.815, 0.769528, 0.355412, 0.435412,
      1.0, 0.705882352941, 0.0156862745098, 0.149019607843,
    ]
  },
}

def GetColorMapInfoFromColorLegendCollection(inColorMapSettings, inMin, inMax):
  global gColorLegendCollection

  inStrId = inColorMapSettings.mColorMapNameId

  #convert underscores (if any) to spaces in key
  inStrId = inStrId.replace('_', ' ')

  if inStrId not in gColorLegendCollection:
    if PhactoriDbg(600):
      myDebugPrint3("GetColorMapInfoFromColorLegendCollection: bad map name:" +\
        str(inStrId) + " using Default\n", 600)
    inStrId = "Default"
  colorLegendItem = gColorLegendCollection[inStrId]
  outRGBPoints = list(colorLegendItem["RGBPoints"])
  outNanColor = list(colorLegendItem["NanColor"])
  outColorSpace = str(colorLegendItem["ColorSpace"])

  #invert color map if requested
  if inColorMapSettings.mInvertColorMap:
    newOutRGBPoints = list()
    #step through the points (4 doubles) backwards and put them on a new
    #list in reverse order, also reversing the position value
    for ii in range(len(outRGBPoints) - 4, -1, -4):
      newOutRGBPoints.append(1.0 - outRGBPoints[ii])
      newOutRGBPoints.append(outRGBPoints[ii+1])
      newOutRGBPoints.append(outRGBPoints[ii+2])
      newOutRGBPoints.append(outRGBPoints[ii+3])
    #myDebugPrint3("uninverted list:")
    #for ii in range(0, len(outRGBPoints), 4):
    #  myDebugPrint3("  " + \
    #      str(outRGBPoints[ii]) + ", " + \
    #      str(outRGBPoints[ii+1]) + ", " + \
    #      str(outRGBPoints[ii+2]) + ", " + \
    #      str(outRGBPoints[ii+3]) + "\n")
    #myDebugPrint3("inverted list:")
    #for ii in range(0, len(newOutRGBPoints), 4):
    #  myDebugPrint3("  " + \
    #      str(newOutRGBPoints[ii]) + ", " + \
    #      str(newOutRGBPoints[ii+1]) + ", " + \
    #      str(newOutRGBPoints[ii+2]) + ", " + \
    #      str(newOutRGBPoints[ii+3]) + "\n")
    outRGBPoints = newOutRGBPoints

  numVals = len(outRGBPoints)
  delta = inMax - inMin
  for ii in range(0, numVals, 4):
    ratio = outRGBPoints[ii]
    outRGBPoints[ii] = inMin + ratio * delta
  return (outRGBPoints, outColorSpace, outNanColor)



myDebugPrint3("gParaViewCatalystVersionFlag is: " + str(gParaViewCatalystVersionFlag) +"\n")


