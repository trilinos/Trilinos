#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
##paraview.simple._DisableFirstRenderCameraReset()

import os

def strnum(i):
  num = ""
  if (i  < 10):
    num = "000" + str(i)
  elif (i < 100):
    num = "00" + str(i)
  elif (i < 1000):
    num = "0" + str(i)
  else:
    num = str(i)
  return num

# create a new 'ExodusIIReader'
fileName = os.getcwd() + '/bridgeFalls.e'
bridgeFallse = ExodusIIReader(FileName=[fileName])
bridgeFallse.ElementVariables = []
bridgeFallse.PointVariables = []
bridgeFallse.GlobalVariables = []
bridgeFallse.SideSetArrayStatus = []

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on bridgeFallse
bridgeFallse.ElementVariables = ['death', 'stress', 'processor_id']
bridgeFallse.PointVariables = ['displ_', 'points']
bridgeFallse.GlobalVariables = ['external_energy', 'internal_energy', 'kinetic_energy', 'momentum_', 'timestep']
bridgeFallse.ElementBlocks = ['block_1', 'block_2', 'block_3', 'block_4', 'block_5', 'block_6']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [640, 480]

# show data in view
bridgeFallseDisplay = Show(bridgeFallse, renderView1)
# trace defaults for the display properties.
bridgeFallseDisplay.Representation = 'Surface'
bridgeFallseDisplay.ColorArrayName = [None, '']
bridgeFallseDisplay.OSPRayScaleArray = 'GlobalNodeId'
bridgeFallseDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
bridgeFallseDisplay.SelectOrientationVectors = 'GlobalNodeId'
bridgeFallseDisplay.ScaleFactor = 20.0
bridgeFallseDisplay.SelectScaleArray = 'GlobalNodeId'
bridgeFallseDisplay.GlyphType = 'Arrow'
bridgeFallseDisplay.GlyphTableIndexArray = 'GlobalNodeId'
bridgeFallseDisplay.DataAxesGrid = 'GridAxesRepresentation'
bridgeFallseDisplay.PolarAxes = 'PolarAxesRepresentation'
bridgeFallseDisplay.ScalarOpacityUnitDistance = 12.640189659732366
bridgeFallseDisplay.GaussianRadius = 10.0
bridgeFallseDisplay.SetScaleArray = ['POINTS', 'GlobalNodeId']
bridgeFallseDisplay.ScaleTransferFunction = 'PiecewiseFunction'
bridgeFallseDisplay.OpacityArray = ['POINTS', 'GlobalNodeId']
bridgeFallseDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(bridgeFallseDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
bridgeFallseDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# create a new 'Threshold'
threshold1 = Threshold(Input=bridgeFallse)
threshold1.Scalars = ['POINTS', 'GlobalNodeId']
threshold1.ThresholdRange = [1.0, 25728.0]

# Properties modified on threshold1
threshold1.Scalars = ['CELLS', 'death']
threshold1.ThresholdRange = [0.1, 25728.0]

# show data in view
threshold1Display = Show(threshold1, renderView1)
# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = [None, '']
threshold1Display.OSPRayScaleArray = 'GlobalNodeId'
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'GlobalNodeId'
threshold1Display.ScaleFactor = 20.0
threshold1Display.SelectScaleArray = 'GlobalNodeId'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'GlobalNodeId'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityUnitDistance = 12.640189659732366
threshold1Display.GaussianRadius = 10.0
threshold1Display.SetScaleArray = ['POINTS', 'GlobalNodeId']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = ['POINTS', 'GlobalNodeId']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(bridgeFallse, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(threshold1Display, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

animationScene1.GoToLast()

# hide color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, False)

animationScene1.GoToFirst()

# current camera placement for renderView1
renderView1.CameraPosition = [0, 0, 300.0]
renderView1.CameraFocalPoint = [0, 0, 0]
renderView1.CameraParallelScale = 100.0

numsteps = len(bridgeFallse.TimestepValues)
for i in range(numsteps):
  picName = 'my_pic_' + strnum(i+1) + '.png'
  # save screenshot
  SaveScreenshot(picName, renderView1, ImageResolution=[640, 480],
    OverrideColorPalette='PrintBackground')
  animationScene1.GoToNext()
