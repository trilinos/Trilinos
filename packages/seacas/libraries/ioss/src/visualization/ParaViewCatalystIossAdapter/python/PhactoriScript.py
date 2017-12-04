# Copyright(C) 1999-2017 National Technology & Engineering Solutions
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
import phactori

#import tutorial_step_12b as current_tutorial_step

global AllAxesPlusQuadrants1And8
AllAxesPlusQuadrants1And8 = \
[ \
   #[[ -1.0,  0.0,  0.0], [800, 800], "_dir_x-_",     1.0,  'Surface'], \
   #[[ -1.0,  1.0,  1.0], [800, 800], "_dir_xyz+_",   1.0,  'Surface With Edges'], \
   [[ -1.0,  1.0,  1.0], [800, 800], "_adjustable",   1.0,  'Surface With Edges', True, None], \
]

global dirSet2
dirSet2 = \
[ \
   #[[ -1.0,  0.0,  0.0], [800, 800], "_dir_x-_",     1.0,  'Surface'], \
   [[ -1.0,  0.0,  1.0], [800, 800], "_dir_xyz+_",   0.15,  'Surface With Edges', False, [0.0, 0.0, -0.4]], \
]

global dirSet3
dirSet3 = \
[ \
   #[[ -1.0,  0.0,  0.0], [800, 800], "_dir_x-_",     1.0,  'Surface'], \
   [[ -1.0,  1.0,  1.0], [800, 800], "_dir_xyz2+_",   1.0,  'Surface', False, [0.0, 0.0, -0.4]], \
]

#AllAxesPlusQuadrants1And8 = \
#[ \
#   [[ -1.0,  0.0,  0.0], [800, 800], "_dir_x-_",   1.0,  'Surface With Edges'], \
#   [[  1.0,  0.0,  0.0], [800, 800], "_dir_x+_",   1.0,  'Surface With Edges'], \
#   [[  0.0, -1.0,  0.0], [800, 800], "_dir_y-_",   1.0,  'Surface With Edges'], \
#   [[  0.0,  1.0,  0.0], [800, 800], "_dir_y+_",   1.0,  'Surface With Edges'], \
#   [[  0.0,  0.0, -1.0], [800, 800], "_dir_z-_",   1.0,  'Surface With Edges'], \
#   [[  0.0,  0.0,  1.0], [800, 800], "_dir_z+_",   1.0,  'Surface With Edges'], \
#   [[ -1.0, -1.0, -1.0], [800, 800], "_dir_xyz+_", 1.0,  'Surface With Edges'], \
#   [[  1.0,  1.0,  1.0], [800, 800], "_dir_xyz-_", 1.0,  'Surface With Edges'], \
#   [[  1.0, -1.0,  1.0], [800, 800], "_dir_xyz2+_", 0.5,  'Surface'], \
#   [[  1.0,  1.0, -1.0], [800, 800], "_dir_xyz2-_", 0.5,  'Surface'], \
#]

#global AllAxesAndAllDiagonals
#AllAxesAndAllDiagonals =
#[
#   [ -1.0, -1.0, -1.0], 
#   [ -1.0, -1.0,  0.0], 
#   [ -1.0, -1.0,  1.0], 
#   [ -1.0,  0.0, -1.0], 
#   [ -1.0,  0.0,  0.0], 
#   [ -1.0,  0.0,  1.0], 
#   [ -1.0,  1.0, -1.0], 
#   [ -1.0,  1.0,  0.0], 
#   [ -1.0,  1.0,  1.0], 
#
#   [  0.0, -1.0, -1.0], 
#   [  0.0, -1.0,  0.0], 
#   [  0.0, -1.0,  1.0], 
#   [  0.0,  0.0, -1.0], 
#   #[  0.0,  0.0,  0.0], 
#   [  0.0,  0.0,  1.0], 
#   [  0.0,  1.0, -1.0], 
#   [  0.0,  1.0,  0.0], 
#   [  0.0,  1.0,  1.0], 
#
#   [  1.0, -1.0, -1.0], 
#   [  1.0, -1.0,  0.0], 
#   [  1.0, -1.0,  1.0], 
#   [  1.0,  0.0, -1.0], 
#   [  1.0,  0.0,  0.0], 
#   [  1.0,  0.0,  1.0], 
#   [  1.0,  1.0, -1.0], 
#   [  1.0,  1.0,  0.0], 
#   [  1.0,  1.0,  1.0], 
#]

global gMyViewMapB
global gMyDataSetupMapB
gMyViewMapB = None
gMyDataSetupMapB = None

def UseViewMapB_ps(inViewMapB):
  global gMyViewMapB
  print 'entering PhactoriScript.UseViewMapB_ps'
  gMyViewMapB = inViewMapB

def UseDataSetupMapB_ps(inDataSetupMapB):
  print 'entering PhactoriScript.UseDataSetupMapB_ps'
  global gMyDataSetupMapB
  gMyDataSetupMapB = inDataSetupMapB

global gDoNewScatterPlotsC
gDoNewScatterPlotsC = True
#gDoNewScatterPlotsC = False

global gRunTestsA
gRunTestsA = False

global gGetJsonViewMapCFromUserData
gGetJsonViewMapCFromUserData = True

def GetViewMapCFromTestSuite():
  #need to do some locking to avoid parallel issues
  try:
    #ff = open('NextTestToRun.txt', 'rb')
    #label = ff.readline()
    #testDirectory = ff.readline()
    #label = ff.readline()
    #testInputDeck = ff.readline()
    #label = ff.readline()
    #testJsonFile = ff.readline()
    #myLen = len(testJsonFile)
    #if testJsonFile[myLen -1 ] == '\n':
    #  testJsonFile = testJsonFile[0:myLen - 1]
    #phactori.myDebugPrint2("running test ->xzq<-\ntest directory: " + \
    #  testDirectory + "\ntest input deck: " + testInputDeck + \
    #  "\ntest json file: ->" + testJsonFile + "<-\n")
    #ff.close()
    import glob
    jsonFilelist = glob.glob('*.json.correct')
    phactori.myDebugPrint2("about to open json file\n")
    phactori.myDebugPrint2("json file name: " + jsonFilelist[0] + "\n")
    jff = open(jsonFilelist[0], 'rb')
    phactori.myDebugPrint2("about to read json file\n")
    testJsonString = jff.read()
    phactori.myDebugPrint2("test json string: " + testJsonString + "\n")
    import json
    testJson = json.loads(testJsonString)
    testJson = phactori.convertJsonUnicodeToStrings(testJson)

  except:
    phactori.myDebugPrint2("GetViewMapCFromTestSuite error\n")
    exit()

  return testJson

def GetViewMapCFromUserData(datadescription):
  """given a datadescription, get the json string view map from the user data
     and convert it into a python dict using the json library, and return
     that item.  Also determine the separator character and input deck
     filename"""
  import json

  #the following code forces a UserData item for internal testing
  #phactori.myDebugPrint2("GetViewMapCFromUserData entered\n")
  #import paraview.vtk as vtk
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


  fd = datadescription.GetUserData()

  if fd == None:
    phactori.myDebugPrint2("no user data, returning {}\n")
    returnViewMapC = {}
    return returnViewMapC

  sa = fd.GetAbstractArray(0)

  testJsonString = sa.GetValue(0)
  separatorCharacterString = sa.GetValue(1)
  inputDeckFilename = sa.GetValue(2)
  phactori.myDebugPrint2("gGetJsonViewMapCFromUserData string:\n" + \
      str(testJsonString) + "\n")
  phactori.myDebugPrint2("separator string: " + separatorCharacterString + \
      "\ninputDeckFilename: " + inputDeckFilename + "\n")
  extraFileString = sa.GetValue(3) 
  phactori.myDebugPrint2("extraFileString: " + extraFileString + "\n")
  phactori.myDebugPrint2("num strings2: " + str(sa.GetNumberOfValues()))
  phactori.SetSeparatorString(separatorCharacterString)
  phactori.SetExtraFileString(extraFileString)
  returnViewMapC = json.loads(testJsonString)
  returnViewMapC = phactori.convertJsonUnicodeToStrings(returnViewMapC)

  phactori.SetDefaultImageBasename(inputDeckFilename)

  phactori.myDebugPrint2("GetViewMapCFromUserData returning\n")

  return returnViewMapC

def CreatePipeline(datadescription):

  #phactori.TestProgrammableFilter()

  #phactori.ApplyClipPlane([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], "clip1")
  #phactori.WarpMeshByDisplacement('DISPL_')
  phactori.TempAddTestFiltersAfterDispl()
  #phactori.ThresholdFilter('STATUS', 'CELLS', [0.7, 100.0], 'DeathStatusFilter')

  #phactori.SetPlotView2StartXYZMinMax(4.0, 131.0, -600000.0, 600000.0, 0.0, 0.0)
  #phactori.SetPlotView2StartXYZMinMax(0.0, 1000.0, -1000.0, 1000.0, 0.0, 0.0)
  #phactori.SetPlotView2StartXYZMinMax(0.0, 1.0, 0.0, 0.01, 0.0, 1.0)

  theSource = phactori.GetCurrentSource()

  theCellArrays = None
  theCellArrays = theSource.CellData
  thePointArrays = None
  thePointArrays = theSource.PointData

  theViewMapC = {
     'camera blocks': {
       'overviewcamera8' : {
         'camera type': 'multicamera8',
         'zoom factor': 1.0,
         'look at point relative to data size': [0.0, 0.0, 0.0]
       },
       'zoom5camera8' : {
         'camera type': 'multicamera8',
         'zoom factor': 0.2,
         'look at point': [0.0, 0.0, 0.0]
       },
       'node234camera8' : {
         'camera type': 'multicamera8',
         'zoom factor': 1.0,
         'look at node': 234
       },
       'element222camera8' : {
         'camera type': 'multicamera8',
         'zoom factor': 1.0,
         'look at element': 222
       },
       'thresh1camera8' : {
         'camera type': 'multicamera8',
         'zoom factor': 2.0,
         'look at point': [0.0, 0.0, 0.0]
       }
     },
     'representation blocks': {
       'rep1': {
         'color by variable': 'VON_MISES',
         #'variable type': 'element',
         'show element': 'surfaces and edges'
       },
       'rep2': {
         'color by variable': 'VON_MISES',
         #'variable type': 'element',
         'show element': 'surfaces'
       },
       'rep3dvecmag': {
         'color by variable 3d vector magnitude': 'DISPL',
         #'variable type': 'node',
         'show element': 'surfaces and edges'
       },
       'rep3dveccompz': {
         'color by variable 3d vector component': 'DISPL_z',
         #'variable type': 'node',
         'show element': 'surfaces and edges'
       },
       'rep3dveccompy': {
         'color by variable 3d vector component': 'DISPL_y',
         #'variable type': 'node',
         'show element': 'surfaces and edges',
         'color map value range': [-0.03, 0.03]
       },
       'rep6dtencompyy': {
         'color by variable 6d tensor component': 'STRESS_yy',
         #'variable type': 'element',
         'show element': 'surfaces and edges'
       }
     },
     #phactori.ThresholdFilter('STATUS', 'CELLS', [0.7, 100.0], 'DeathStatusFilter')
     'operation blocks': {
       'death_status_filter_1': {
         'type': 'threshold',
         'variable': 'STATUS',
         'range': [0.7, 100.0],
         'keep inside range': True,
         'variable type': 'element'
       },
       'clip_plane_1': {
         'type': 'clip',
         'point on plane': [0.0, 0.0, 0.0],
         'plane normal': [1.0, 1.0, 0.0],
         'crinkle cut': True
       },
       'clip_plane_2': {
         'type': 'clip',
         'point on plane': [0.0, 0.0, 0.0],
         'plane normal': [1.0, 1.0, 0.0],
         'crinkle cut': True,
         'input operation': 'death_status_filter_1'
       },
       #'vonmises_threshold_1': {
       #  'type': 'threshold',
       #  'variable': 'VON_MISES',
       #  'range': [50000.0, 5000000.0],
       #  #'range': [-1.0, 50000.0],
       #  'keep inside range': True,
       #  'variable type': 'element'
       #}
     },
     'imageset blocks': {
       #'defaultimageset': {
       #},
       #'defaultimageset2': {
       #  'look at point': [0.0, 0.0, 0.0],
       #  'color by variable': 'VON_MISES',
       #  #'variable type': 'element',
       #  'zoom factor': 0.2,
       #},
       #'defaultimageset3': {
       #  'look at point': [0.0, 0.0, 0.0],
       #  'color by variable': 'DISPL_',
       #  #'variable type': 'element',
       #  'zoom factor': 0.2,
       #},
       #'defaultzoom5': {
       #  'zoom factor': 0.2,
       #  'look at point': [0.0, 0.0, 0.0],
       #  'image basename': 'default_zoom5_',
       #  'representation': 'rep2',
       #},
       #'overviewimageset': {
       #  'representation': 'rep2',
       #  'camera': 'overviewcamera8',
       #  'image basename': 'overview_'
       #},
       #'zoom5imageset': {
       #  'representation': 'rep1',
       #  'camera': 'zoom5camera8',
       #  'image basename': 'zoom5_'
       #},
       #'zoom5imageset2': {
       #  'representation': 'rep1',
       #  'camera': 'zoom5camera8',
       #  'image basename': 'zoom5_nodeath_',
       #  'operation': 'death_status_filter_1'
       #},
       #'clip_no_thresh_overview': {
       #  'representation': 'rep1',
       #  'camera': 'overviewcamera8',
       #  'image basename': 'clip_no_thresh_overview_',
       #  'operation': 'clip_plane_1'
       #},
       #'clip_no_thresh': {
       #  'representation': 'rep1',
       #  'camera': 'zoom5camera8',
       #  'image basename': 'clip_no_thresh_zoom5_',
       #  'operation': 'clip_plane_1'
       #},
       #'clip_with_thresh': {
       #  'representation': 'rep1',
       #  'camera': 'zoom5camera8',
       #  'image basename': 'clip_with_thresh_zoom5_',
       #  'operation': 'clip_plane_2'
       #},
       #'zoom5vecmag': {
       #  'representation': 'rep3dvecmag',
       #  'camera': 'zoom5camera8',
       #  'image basename': 'zoom5vecmag_'
       #},
       #'zoom5veccompz': {
       #  'representation': 'rep3dveccompz',
       #  'camera': 'zoom5camera8',
       #  'image basename': 'zoom5veccompz_'
       #},
       'zoom5veccompy': {
         'representation': 'rep3dveccompy',
         'camera': 'zoom5camera8',
         'image basename': 'zoom5veccompy_'
       },
       #'zoom5tenyy': {
       #  'representation': 'rep6dtencompyy',
       #  'camera': 'zoom5camera8',
       #  'image basename': 'zoom5tenyy_'
       #},
       #'node234': {
       #  'representation': 'rep1',
       #  'camera': 'node234camera8',
       #  'image basename': 'node234_'
       #},
       #'element222': {
       #  'representation': 'rep1',
       #  'camera': 'element222camera8',
       #  'image basename': 'element222_'
       #},
       #'thresh1set': {
       #  'representation': 'rep1',
       #  'camera': 'thresh1camera8',
       #  'image basename': 'thresh1_',
       #  'operation': 'vonmises_threshold_1'
       #}
     },
     'scatter plot blocks': {
       'VonMisesVsElementId': {
         'x axis variable': 'GlobalElementId',
         'y axis variable': 'VON_MISES',
         'image basename': 'vm_vs_id_',
         'variable type': 'element'
       },
       'displacementVsNodeId': {
         'x axis variable': 'GlobalNodeId',
         'y axis variable': 'DISPL_',
         'image basename': 'displ_vs_id_',
         'variable type': 'node'
       }
     },
     'plot over time blocks': {
       'VonMisesVsTime': {
         'y axis variable': 'VON_MISES',
         'image basename': 'vm_over_time_',
         'variable type': 'element'
       },
       'DisplacementVsTime': {
         'y axis variable': 'DISPL_',
         'image basename': 'displacement_over_time_',
         'variable type': 'node'
       }
     }
    }

  myLookAtPoint = [0.0, 0.0, 0.0]
  theViewMapC = {
     'camera blocks': {
     },
     'representation blocks': {
     },
     'operation blocks': {
       #'element_to_node': {
       #  'type': 'element data to node data',
       #},
       #'contour1': {
       #  'type': 'contour',
       #  #'input operation': 'element_to_node',
       #  'variable': 'VON_MISES',
       #  'contour value': 50000.0,
       #  'variable type': 'node',
       #},
     },
     'imageset blocks': {
       'zoom1': {
         'color by scalar': 'VON_MISES',
         'image digit count': 5,
         'look at point absolute': [0.0, 0.0, 0.0],
         'zoom factor': 1.0,
         'image basename': 'zoom_vm_'
         #'operation': 'element_to_node'
       },
       'zoom2': {
         'color by vector component': 'VEL_Y',
         'image digit count': 5,
         'look at point absolute': [0.0, 0.0, 0.0],
         'zoom factor': 1.0,
         'image basename': 'zoom_vel_y_'
         #'operation': 'element_to_node'
       },
       #'zcontour_imageset': {
       #  'color by variable': 'VON_MISES',
       #  'operation': 'contour1',
       #  'image basename': 'contour_',
       #},
     },
     'scatter plot blocks': {
       'VonMisesPointVsNodeId': {
         'x axis variable scalar': 'GlobalElementId',
         'y axis variable scalar': 'VON_MISES',
         'variable type': 'element',
         #'image basename': 'vm_vs_id_',
         #'variable type': 'node',
         #'operation': 'element_to_node'
       },
       #'SurfVonMisesPointVsNodeId': {
       #  'x axis variable': 'GlobalNodeId',
       #  'y axis variable': 'VON_MISES',
       #  'image basename': 'surf_vm_vs_id_',
       #  'variable type': 'node',
       #  #'operation': 'element_to_node'
       #  'operation': 'contour1'
       #},
     },
     'plot over time blocks': {
     }
    }

  global gRunTestsA
  if gRunTestsA:
    theViewMapC = GetViewMapCFromTestSuite()

  global gGetJsonViewMapCFromUserData
  if gGetJsonViewMapCFromUserData == True:
      theViewMapC = GetViewMapCFromUserData(datadescription)

  phactori.CreateViewSetFromPhactoriViewMapC(theViewMapC)

  global gMyViewMapB
  if gMyViewMapB != None:
    phactori.CreateViewSetFromPhactoriViewMapB(gMyViewMapB)

  #view_setup_2 = { 'color_variable': {'name': 'VON_MISES' } }
  #phactori.SetPlotView2StartColors(
  #  inBackgroundColor = [0.0, 0.0, 0.0],
  #  inEdgeColor = [1.0, 1.0, 1.0])

  phactori.SetPlotView2StartColors(
    inBackgroundColor = [1.0, 1.0, 1.0],
    inEdgeColor = [0.0, 0.0, 0.5],
    inCubeAxesColor = [0.2, 0.2, 0.2],
    inDiffuseColor = [0.2, 0.2, 0.2],
    inAmbientColor = [0.2, 0.2, 0.2],
    inSelectionColor = [0.2, 0.2, 0.2],
    inBackfaceDiffuseColor = [0.2, 0.2, 0.2])

  #test1 = phactori.CollectCells1('TEAR_DOUBLE', 0.01, 1)
  #test2 = phactori.CollectCells1('STATUS', 0.8, -1)
  #compareTearDeath(test1, test2)
  #GetAndReduceViewControl()

global gSavedViewControlDataList
gSavedViewControlDataList = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
global gViewControlTriggerTime
gViewControlTriggerTime = 0.0


import os

def GetViewControlFromTextFile(outDataList):

  global gSavedViewControlDataList
  global gViewControlTriggerTime

  triggerTime = os.path.getmtime("CatalystViewControl_trigger.txt")
  if triggerTime <= gViewControlTriggerTime:
    for ii in range(0, 12):
      outDataList[ii] = gSavedViewControlDataList[ii]
    return

  gViewControlTriggerTime = triggerTime

  ff = open('CatalystViewControl.txt', 'rb')
  fileLines = ff.readlines()
  xx = float(fileLines[1])
  yy = float(fileLines[3])
  zz = float(fileLines[5])
  zoom = float(fileLines[7])
  fx = float(fileLines[9])
  fy = float(fileLines[11])
  fz = float(fileLines[13])
  showColorKey = float(fileLines[15])
  showAxes = float(fileLines[17])
  showOrientation = float(fileLines[19])
  colorVariableCellOrPoint = float(fileLines[21])
  colorVariableName = str(fileLines[23])
  endIndex = len(colorVariableName)-1
  if colorVariableName[endIndex] == '\n':
    #newEndIndex = len(colorVariableName)-1
    #colorVariableName = colorVariableName[0..newEndIndex]
    newColorVar = ""
    for ii in range(0,endIndex):
      newColorVar = newColorVar + colorVariableName[ii]
    colorVariableName = newColorVar
    
  ff.close()

  #theSource = GetActiveSource()
  theSource = phactori.GetCurrentSource()

  theArrays = None
  if(colorVariableCellOrPoint == 1):
    phactori.myDebugPrint("GetViewControlFromTextFile doing cell array\n")
    theArrays = theSource.CellData
  else:
    phactori.myDebugPrint("GetViewControlFromTextFile doing point array\n")
    theArrays = theSource.PointData

  colorVariableIndex = 0
  phactori.myDebugPrint("GetViewControlFromTextFile finding ->" + colorVariableName + "<-\n")
  for ii in range(theArrays.GetNumberOfArrays()):
      phactori.myDebugPrint("GetViewControlFromTextFile trying ->" + theArrays.GetArray(ii).Name + "<-\n")
      if theArrays.GetArray(ii).Name == colorVariableName:
        phactori.myDebugPrint("match! " + str(ii) + "\n");
        colorVariableIndex = ii
      else:
        phactori.myDebugPrint("no match\n");

  phactori.myDebugPrint('PhactoriScript.GetViewControlFromTextFile   [' + \
    str(xx) + ', ' + str(yy) + ', ' + str(zz) + '] ' + str(zoom) +
    ' [' + str(fx) + ', ' + str(fy) + ',' + str(fz) + 
    ']\nshowColorKey: ' + str(showColorKey) + ' showAxes: ' + str(showAxes) + \
    ' showOrientation: ' + str(showOrientation) + \
    '\ncolorVariableCellOrPoint: ' + str(colorVariableCellOrPoint) + \
    ' colorVariableIndex: ' + str(colorVariableIndex) + '\n')
  outDataList[0] = xx
  outDataList[1] = yy
  outDataList[2] = zz
  outDataList[3] = zoom
  outDataList[4] = fx
  outDataList[5] = fy
  outDataList[6] = fz
  outDataList[7] = showColorKey
  outDataList[8] = showAxes
  outDataList[9] = showOrientation
  outDataList[10] = colorVariableCellOrPoint
  outDataList[11] = colorVariableIndex
  for ii in range(0, 12):
    gSavedViewControlDataList[ii] = outDataList[ii]

def GetAndReduceViewControl():
  phactori.myDebugPrint('GetAndReduceViewControl entered\n')
  myProcId = phactori.SmartGetLocalProcessId()
  dataList = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  if(myProcId == 0):
    phactori.myDebugPrint('  I am process 0, reading file\n')
    GetViewControlFromTextFile(dataList)
  else:
    phactori.myDebugPrint('  I am not process 0, reading file\n')

  phactori.UseReduceToSpreadValues(dataList)

  xx = dataList[0]
  yy = dataList[1]
  zz = dataList[2]
  zoom = dataList[3]
  fx = dataList[4]
  fy = dataList[5]
  fz = dataList[6]
  showColorKey =             dataList[7]
  showAxes =                 dataList[8]
  showOrientation =          dataList[9]
  colorVariableCellOrPoint = dataList[10]
  colorVariableIndex =       dataList[11]

  #theSource = GetActiveSource()
  theSource = phactori.GetCurrentSource()

  theArrays = None
  if(colorVariableCellOrPoint == 1):
    theArrays = theSource.CellData
  else:
    theArrays = theSource.PointData

  phactori.SetCameraLookAtPointAndLookDirection(\
    inLookDirection = [xx, yy, zz], \
    inFocalPoint = [fx, fy, fz], \
    inEyePositionFactor = zoom)

  phactori.ShowDataColorLegend('off')
  phactori.myDebugPrint("  GetAndReduceViewControl color variable index " + \
    str(colorVariableIndex) + " name: " + \
    theArrays[colorVariableIndex].Name + "\n")
  phactori.ColorByVariable(theArrays[colorVariableIndex].Name)
  if showColorKey == 1.0:
    phactori.ShowDataColorLegend('on')
  else:
    phactori.ShowDataColorLegend('off')
  if showAxes == 1.0:
    phactori.ShowCubeAxes('on')
  else:
    phactori.ShowCubeAxes('off')
  if showOrientation == 1.0:
    phactori.SetOrientationAxesVisibility('on')
  else:
    phactori.SetOrientationAxesVisibility('off')

  phactori.myDebugPrint('GetAndReduceViewControl returning\n')


global tearListPersistent
tearListPersistent = []
global deathListPersistent
deathListPersistent = []

def mergeCipCompare(item1, item2):
  if item1[1] > item2[1]:
   return 1
  if item1[1] < item2[1]:
   return -1
  if item1[2] > item2[2]:
   return 1
  if item1[2] < item2[2]:
   return -1
  return 0

def mergeCurrentIntoPersistent(pList, cList):
  phactori.myDebugPrint2('PhactoriScript.mergeCurrentIntoPersistent entered\n')
  pListLen = len(pList)
  cListLen = len(cList)
  phactori.myDebugPrint2('  pList has ' + str(pListLen) + ' elements, cList has ' + str(cListLen) + ' elements\n')
  pIndex = 0
  cIndex = 0
  mergeList = []
  while (pIndex < pListLen) or (cIndex < cListLen):
    if(pIndex >= pListLen):
      mergeList.append(cList[cIndex])
      #phactori.myDebugPrint2('  pIndex ' + str(pIndex) + ' cIndex ' + str(cIndex) + '\n')
      #phactori.myDebugPrint2('  (p ended) [cIndex] ' + str(cList[cIndex]) + '\n')
      cIndex += 1
      continue
    if(cIndex >= cListLen):
      mergeList.append(pList[pIndex])
      #phactori.myDebugPrint2('  pIndex ' + str(pIndex) + ' cIndex ' + str(cIndex) + '\n')
      #phactori.myDebugPrint2('  (c ended) [pIndex] ' + str(cList[pIndex]) + '\n')
      pIndex += 1
      continue
    compareResult = mergeCipCompare(cList[cIndex], pList[pIndex])
    #phactori.myDebugPrint2('  pIndex ' + str(pIndex) + ' cIndex ' + str(cIndex) + '\n')
    #phactori.myDebugPrint2('  [pIndex] ' + str(pList[pIndex]) + ' [cIndex] ' + str(cList[cIndex]) + '\n')
    if compareResult == 1:
      #phactori.myDebugPrint2('  item from pList is less\n')
      mergeList.append(pList[pIndex])
      pIndex += 1
    if compareResult == -1:
      #phactori.myDebugPrint2('  item from cList is less\n')
      mergeList.append(cList[cIndex])
      cIndex += 1
    if compareResult == 0:
      #phactori.myDebugPrint2('  items from cList and pList are same element\n')
      mergeList.append(pList[pIndex])
      pIndex += 1
      cIndex += 1
  phactori.myDebugPrint2('  merged list has ' + str(len(mergeList)) + ' elements\n')
  phactori.myDebugPrint2('PhactoriScript.mergeCurrentIntoPersistent returning\n')
  return mergeList

def compareTearDeath(tList, dList):
  phactori.myDebugPrint2('PhactoriScript.compareTearDeath entered\n')
  phactori.myDebugPrint2('PhactoriScript.compareTearDeath returning\n')
   
def PerFrameUpdate():
  phactori.myDebugPrint('PhactoriScript.PerFrameUpdate entered\n')
  #GetAndReduceViewControl()
  global tearListPersistent
  global deathListPersistent

  global gDoNewScatterPlotsC
  if gDoNewScatterPlotsC:
    phactori.myDebugPrint2('doing phactori.UpdateAllScatterPlots (C)\n')
    phactori.UpdateAllScatterPlots()
    phactori.myDebugPrint2('doing phactori.UpdateAllPlotsOverTime\n')
    phactori.UpdateAllPlotsOverTime()
    phactori.myDebugPrint2('did plot updates\n')


  #currentFrameTearList = phactori.CollectCells1('TEAR_DOUBLE', 0.01, 1)
  #tearListPersistent = mergeCurrentIntoPersistent(tearListPersistent, currentFrameTearList)

  #currentFrameDeathList = phactori.CollectCells1('STATUS', 0.8, -1)
  #deathListPersistent = mergeCurrentIntoPersistent(deathListPersistent, currentFrameDeathList)

  #compareTearDeath(tearListPersistent, deathListPersistent)

  phactori.myDebugPrint('PhactoriScript.PerFrameUpdate exiting\n')




