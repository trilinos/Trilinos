# script-version: 2.0
# Catalyst state generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [910, 528]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.013602, 0.0325, 0.005000000000000002]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.06420701247294598, -0.03333954738404109, 0.1323247768943745]
renderView1.CameraFocalPoint = [0.013602, 0.0325, 0.005000000000000002]
renderView1.CameraViewUp = [0.4116987268415183, 0.8658965574555053, 0.28412551821999]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.03934331460362739
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(910, 528)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'CGNS Series Reader'
#aero_blunt_wedge_test3cgns = CGNSSeriesReader(registrationName='input', FileNames=['C:\\Users\\jamauld\\sparc_work\\data_2021Mar03\\aero_blunt_wedge_test3.cgns'])
aero_blunt_wedge_test3cgns = Wavelet(registrationName='input')
aero_blunt_wedge_test3cgns.Bases = ['Base']
aero_blunt_wedge_test3cgns.Families = ['bc-1-wall', 'bc-2-inflow', 'bc-3-outflow', 'bc-4-symm-y', 'bc-5-symm-z1', 'bc-6-symm-z2']
aero_blunt_wedge_test3cgns.CellArrayStatus = ['density', 'pressure', 'temperature', 'velocity_x', 'velocity_y', 'velocity_z']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from aero_blunt_wedge_test3cgns
aero_blunt_wedge_test3cgnsDisplay = Show(aero_blunt_wedge_test3cgns, renderView1, 'StructuredGridRepresentation')

# get color transfer function/color map for 'pressure'
pressureLUT = GetColorTransferFunction('pressure')
pressureLUT.RGBPoints = [159147.93283308673, 0.231373, 0.298039, 0.752941, 660581.9293834264, 0.865003, 0.865003, 0.865003, 1162015.9259337662, 0.705882, 0.0156863, 0.14902]
pressureLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'pressure'
pressurePWF = GetOpacityTransferFunction('pressure')
pressurePWF.Points = [159147.93283308673, 0.0, 0.5, 0.0, 1162015.9259337662, 1.0, 0.5, 0.0]
pressurePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
aero_blunt_wedge_test3cgnsDisplay.Representation = 'Surface With Edges'
aero_blunt_wedge_test3cgnsDisplay.ColorArrayName = ['CELLS', 'pressure']
aero_blunt_wedge_test3cgnsDisplay.LookupTable = pressureLUT
aero_blunt_wedge_test3cgnsDisplay.SelectTCoordArray = 'None'
aero_blunt_wedge_test3cgnsDisplay.SelectNormalArray = 'None'
aero_blunt_wedge_test3cgnsDisplay.SelectTangentArray = 'None'
aero_blunt_wedge_test3cgnsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
aero_blunt_wedge_test3cgnsDisplay.SelectOrientationVectors = 'None'
aero_blunt_wedge_test3cgnsDisplay.ScaleFactor = 0.006500000000000001
aero_blunt_wedge_test3cgnsDisplay.SelectScaleArray = 'None'
aero_blunt_wedge_test3cgnsDisplay.GlyphType = 'Arrow'
aero_blunt_wedge_test3cgnsDisplay.GlyphTableIndexArray = 'None'
aero_blunt_wedge_test3cgnsDisplay.GaussianRadius = 0.00032500000000000004
aero_blunt_wedge_test3cgnsDisplay.SetScaleArray = [None, '']
aero_blunt_wedge_test3cgnsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
aero_blunt_wedge_test3cgnsDisplay.OpacityArray = [None, '']
aero_blunt_wedge_test3cgnsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
aero_blunt_wedge_test3cgnsDisplay.DataAxesGrid = 'GridAxesRepresentation'
aero_blunt_wedge_test3cgnsDisplay.PolarAxes = 'PolarAxesRepresentation'
aero_blunt_wedge_test3cgnsDisplay.ScalarOpacityFunction = pressurePWF
aero_blunt_wedge_test3cgnsDisplay.ScalarOpacityUnitDistance = 0.031226809494856976

# setup the color legend parameters for each legend in this view

# get color legend/bar for pressureLUT in view renderView1
pressureLUTColorBar = GetScalarBar(pressureLUT, renderView1)
pressureLUTColorBar.Title = 'pressure'
pressureLUTColorBar.ComponentTitle = ''

# set color bar visibility
pressureLUTColorBar.Visibility = 1

# show color legend
aero_blunt_wedge_test3cgnsDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
#pNG1.Writer.FileName = 'RenderView1_%.6ts%cm.png'
pNG1.Writer.FileName = 'test8_%.6ts%cm.png'
pNG1.Writer.ImageResolution = [910, 528]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.ExtractsOutputDirectory = '.'
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
