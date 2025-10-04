#!/usr/bin/env python
import sys
import os
import vtk
import random
from subprocess import Popen, PIPE, STDOUT
from vtk import *

############################################
# read in nodes file
def read_nodecoords_from_file(filename):
  global_nodecoords = []
  for l in file(filename):
    line = l.strip()
    if len(line)==1:
      dimension = int(line)
    else:
      line = line.split(" ")
      x = float(line[0]) #+ 0.0001 * random.random()
      y = float(line[1]) #+ 0.0001 * random.random()
      if len(line) == 3:
	z = float(line[2]) #+ 0.0001 * random.random()
	global_nodecoords.append((x,y,z))
      else:
	global_nodecoords.append((x,y))
  return global_nodecoords,dimension

def read_finelevel_nodecoords_from_file(filename):
  global_nodecoords = []
  dimension = 2 # we only support 2D aggregates here...
  for l in file(filename):
    line = l.strip()
    if line.find("#") == 0:
      continue
    else:
      while '  ' in line:
	line = line.replace('  ', ' ')
      line = line.split(" ")
      x = float(line[2]) #+ 0.0001 * random.random()
      y = float(line[3]) #+ 0.0001 * random.random()
      global_nodecoords.append((x,y))
  return global_nodecoords,dimension

# dimension: problem dimension, i.e. either 2 or 3
# agg_globalnodeidx: global node ids of current aggregate
# global_nodecoords: coordinates of nodes
# aggpolygons: vtk append object for aggregates
# aggid: aggregate id (global)
# aggid2nodes: map aggid -> nodeid  (for visualization)
# aggid2procs: map aggid -> proc id (for visualization)
def prepareDelaunayData3d(dimension, agg_globalnodeidx, global_nodecoords, aggpolygons, aggid, aggid2nodes, aggid2procs):

  local_nodeidx2global_nodeidx = {}
  no_of_aggnodes = len(agg_globalnodeidx)
  no_aggs = len(aggid2nodes)

  Points   = vtk.vtkPoints()
  Vertices = vtk.vtkCellArray()

  for i in range(0,len(agg_globalnodeidx)):
    id = -1
    local_nodeidx2global_nodeidx[i] = agg_globalnodeidx[i]
    nodecoords = global_nodecoords[int(agg_globalnodeidx[i])]
    id = Points.InsertNextPoint(nodecoords[0]+ 0.0001 * random.random(),nodecoords[1]+ 0.0001 * random.random(),nodecoords[2]+ 0.0001 * random.random())
    Vertices.InsertNextCell(1)
    Vertices.InsertCellPoint(id)

  polydata2 = vtk.vtkPolyData()
  polydata2.SetPoints(Points)
  polydata2.Modified()
  polydata2.Update()

  delaunay = vtk.vtkDelaunay3D()
  delaunay.SetInput(polydata2)
  delaunay.Update()

  # create surfaceFilter
  surfaceFilter = vtk.vtkDataSetSurfaceFilter()
  surfaceFilter.SetInputConnection(delaunay.GetOutputPort())
  surfaceFilter.Update()

  pt_polydata = surfaceFilter.GetOutput()

  lookupTable = vtk.vtkLookupTable()
  lookupTable.SetNumberOfTableValues(no_aggs)
  lookupTable.Build()

  Ids = vtk.vtkUnsignedIntArray()
  Ids.SetNumberOfComponents(1)
  Ids.SetName("Ids")
  for i in range(0,Points.GetNumberOfPoints()):
    Ids.InsertNextTuple1(int(aggid))
  Ids.SetLookupTable(lookupTable)

  Procs = vtk.vtkUnsignedCharArray()
  Procs.SetNumberOfComponents(1)
  Procs.SetName("proc")
  for i in range(0,Points.GetNumberOfPoints()):
    Procs.InsertNextTuple1(aggid2procs[aggid])

  polydata3 = vtk.vtkPolyData()
  polydata3 = surfaceFilter.GetOutput()
  polydata3.GetPointData().SetScalars(Ids)
  polydata3.GetPointData().AddArray(Procs)

  polydata4 = vtk.vtkPolyData()
  polydata4.SetPoints(Points)
  polydata4.SetVerts(Vertices)
  polydata4.GetPointData().SetScalars(Ids)
  polydata4.GetPointData().AddArray(Procs)

  #datamapper = vtk.vtkDataSetMapper()
  #datamapper.SetInputConnection(delaunay.GetOutputPort())
  #datamapper.SetInput(polydata3)

  #actor = vtk.vtkActor()
  #actor.SetMapper(datamapper)

  #renderer = vtk.vtkRenderer()
  #renderWindow = vtk.vtkRenderWindow()
  #renderWindow.AddRenderer(renderer)
  #renderWindowInteractor = vtk.vtkRenderWindowInteractor()
  #renderWindowInteractor.SetRenderWindow(renderWindow)
  #renderer.AddActor(actor)
  #renderWindow.Render()
  #renderWindowInteractor.Start()

  #print polydata.GetVertices()

  aggpolygons.AddInput(polydata3)
  aggpolygons.AddInput(polydata4)

# dimension: problem dimension, i.e. either 2 or 3
# agg_globalnodeidx: global node ids of current aggregate
# global_nodecoords: coordinates of nodes
# aggpolygons: vtk append object for aggregates
# aggid: aggregate id (global)
# aggid2nodes: map aggid -> nodeid  (for visualization)
# aggid2procs: map aggid -> proc id (for visualization)
def prepareDelaunayData(dimension, agg_globalnodeidx, global_nodecoords, aggpolygons, aggid, aggid2nodes, aggid2procs):
  local_nodeidx2global_nodeidx = {}
  no_of_aggnodes = len(agg_globalnodeidx)
  dim = len(global_nodecoords[0])

  no_aggs = len(aggid2nodes)

  Points   = vtk.vtkPoints()
  Vertices = vtk.vtkCellArray()

  for i in range(0,len(agg_globalnodeidx)):
    local_nodeidx2global_nodeidx[i] = agg_globalnodeidx[i]
    nodecoords = global_nodecoords[int(agg_globalnodeidx[i])]
    if dimension==2:
      id = Points.InsertNextPoint(nodecoords[0],nodecoords[1],0.0)
    elif dimension==3:
      id = Points.InsertNextPoint(nodecoords[0]+ 0.001 * random.random(),nodecoords[1]+ 0.001 * random.random(),nodecoords[2]+ 0.001 * random.random())
    Vertices.InsertNextCell(1)
    Vertices.InsertCellPoint(id)

  # create polygon for current aggregate
  polydata = vtk.vtkPolyData()
  polydata.SetPoints(Points)
  polydata.SetVerts(Vertices)
  polydata.Modified()
  polydata.Update()

  polydata2 = vtk.vtkPolyData()
  if Points.GetNumberOfPoints()>2: # todo: avoid error messages + add support for lines/surfaces
    # create delaunay object
    if dimension==2:
      delaunay = vtk.vtkDelaunay2D()
    elif dimension==3:
      delaunay = vtk.vtkDelaunay3D()
      #delaunay.SetAlpha(0.1)
    delaunay.SetInput(polydata)
    delaunay.Update()

    # create surfaceFilter
    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    surfaceFilter.SetInputConnection(delaunay.GetOutputPort())
    surfaceFilter.Update()

    polydata2 = surfaceFilter.GetOutput()

  lookupTable = vtk.vtkLookupTable()
  lookupTable.SetNumberOfTableValues(no_aggs)
  lookupTable.Build()

  Ids = vtk.vtkUnsignedIntArray()
  Ids.SetNumberOfComponents(1)
  Ids.SetName("Ids")
  for i in range(0,Points.GetNumberOfPoints()):
    Ids.InsertNextTuple1(int(aggid))
  Ids.SetLookupTable(lookupTable)

  Procs = vtk.vtkUnsignedCharArray()
  Procs.SetNumberOfComponents(1)
  Procs.SetName("proc")
  for i in range(0,Points.GetNumberOfPoints()):
    Procs.InsertNextTuple1(aggid2procs[aggid])


  polydata2.SetPoints(Points)
  polydata2.SetVerts(Vertices)
  polydata2.GetPointData().SetScalars(Ids)
  polydata2.GetPointData().AddArray(Procs)
  polydata2.Modified()
  polydata2.Update()

  aggpolygons.AddInput(polydata2)


################################################################################
# READ IN AGGREGATES
# routines that read aggregation information from files

def checkAggregateLine(line):
	if line.find("Agg ") == 0:
		return 1,line
	else:
		return 0,line

def read_aggregates_from_file(filename,procid):
  	aggid2nodes = {}
  	aggid2procs = {}
	for l in file(filename):
		line = l.strip()

		# filter out only Agg lines
		ret,line = checkAggregateLine(line)
		if ret == 0:
			continue

		# line now contains a list of all tokens in that line
		line = line.split(": ")

		# extract aggid and proc number
		agginfo = line[0]
		agginfo = agginfo.split(" ")

		aggid = agginfo[1]
		procid = agginfo[3]

		# handle node ids for aggregate
		aggnodeids = line[1]
		aggnodeids = aggnodeids.split(" ")

		# fill in data variables
		aggid2nodes[aggid] = aggnodeids
		aggid2procs[aggid] = int(procid)

	return aggid2nodes,aggid2procs

################################################################################
# read in aggregation info from file
# input: filename_prototype string with prototype for filename, e.g. aggs_level%LEVEL_proc%PROCS.outlines
#        the variables %LEVEL and %PROCS are replaced by the corresponding values
# input: procs: number of processors (4 means that information from processors 0..3 is expected)
# input: level: level number
def readin_aggregates(filename_prototype,procs,level):
  aggid2nodes = {}
  aggid2procs = {}
  for proc in range(0,procs):
    #filename = "aggs_level" + str(level) + "_proc" + str(proc) + ".out"
    filename = filename_prototype
    filename = filename.replace("%LEVEL",str(level))
    filename = filename.replace("%PROC",str(proc))
    print "process ", filename
    if os.path.exists(filename):
      [aggid2nodesfromproc,aggid2procsfromproc] = read_aggregates_from_file(filename,proc)
      aggid2nodes.update(aggid2nodesfromproc)
      aggid2procs.update(aggid2procsfromproc)

  return aggid2nodes,aggid2procs

################################################################################
# HELPER ROUTINES
# for generating next level information (nodesX.txt)

################################################################################
# get_agg_coords (helper function for get_rootnodes)
# input: list of nodes (coordinates)
# input: map aggid2gids: map of local agg ids -> GIDs
# input: aggid: global aggregate id
# output: returns set of node coordinates for aggregate with global aggid
def get_agg_coords(nodes,aggid2nodes,aggid):
  agg_nodes = aggid2nodes[aggid]
  nodeset = []
  for node in range(0,len(agg_nodes)):
    nodeset.append(nodes[int(agg_nodes[node])])
  return nodeset

################################################################################
# get rootnodes
# input: aggs = map: aggid -> list of nodeids in this agg
# output: list of rootnodes
# note: we calculate the "midpoint" of each aggregate
# TODO extend me for 3d!
def get_rootnodes(aggid2nodes,nodes):
	dim = 2
	if len(nodes[0]) == 3:
		dim = 3

	rootnodes = []
	for i in range(0,len(aggid2nodes.keys())):
		rootnodes.append((0,0))
	for k in aggid2nodes.keys():

		nodecoords = get_agg_coords(nodes,aggid2nodes,k)

		x = 0.0
		y = 0.0
		z = 0.0
		for m in nodecoords:

			x = x + m[0]
			y = y + m[1]
			if dim==3:
				z = z + m[2]
		x = x/len(aggid2nodes[k])
		y = y/len(aggid2nodes[k])
		if dim == 3:
			z = z/len(aggid2nodes[k])

		if dim == 2:
			rootnodes[int(k)] = (x,y)
		elif dim == 3:
			rootnodes[int(k)] = (x,y,z)
		else: print "error: dim must be 2 or 3 but it is " + str(dim)

	return rootnodes

# write nodes file
# input: filename: filename for nodes file (should follow nodeX.txt style)
# input: aggid2nodes map for aggid to list of global nodeidx
# input: nodes list of node coordinates
# input: dimension (2 or 3)
def write_nodes_file(filename,aggid2nodes,nodes,dimension):

  # calculate root nodes (works only for 2d)
  rootnodes = get_rootnodes(aggid2nodes,nodes)

  # write nodes file
  f = open(filename,"w")
  f.write(str(dimension))
  f.write("\r\n")
  for i in range(len(rootnodes)):
    rootnode = rootnodes[i]
    f.write(str(rootnode[0]))
    f.write(" ")
    f.write(str(rootnode[1]))
    if len(rootnode)==3:
      f.write(" ")
      f.write(str(rootnode[2]))
    f.write("\r\n")
  f.close()
  print "node file " + filename + " generated: OK"

################################################################################
# check if all files exist to proceed with next level
# we need a nodesX.txt file for the node coordinates
# and all aggregation information files (from the AggregationExportFactory)
# input: nextlevel: id for next level
# procs: number of procs
# file_prototype: prototype for filename of aggregation information
def check_files_for_next_level(nextlevel,procs,file_prototype):

  if nextlevel==0:
    if os.path.isfile("example.txt") == False:
      return False
  else:
    # check if coarse level node coordinates are available
    if os.path.isfile("nodes"+str(nextlevel)+".txt") == False:
      return False

  #for p in range(0,procs):
  for p in range(0,1): # check only processor one
    filename = file_prototype
    filename = filename.replace("%LEVEL",str(nextlevel))
    filename = filename.replace("%PROC",str(p))
    if os.path.isfile(filename) == False:
      return False

  return True

###########
# MAIN routine
def main(argv=None):
  dimension = 2
  numprocs = 2
  level = 0  # startlevel

  no_multigridlevels = 0

  # check how many processors generated aggregation output
  #while check_files_for_next_level(0,numprocs, "aggs_level%LEVEL_proc%PROC.out") == True:
  #  numprocs = numprocs + 1
  #numprocs = numprocs - 1
  #print "Aggregtaion information for " + str(numprocs) + " processors found"

  # process all multigrid levels
  while check_files_for_next_level(level,numprocs,"aggs_level%LEVEL_proc%PROC.out"):
    global_nodecoords = []

    print "Level " + str(level)

    if level==0:  # read in coordinates (finest level
      global_nodecoords,dimension = read_finelevel_nodecoords_from_file("example.txt")
    else:
      global_nodecoords,dimension = read_nodecoords_from_file("nodes"+str(level)+".txt")

    # read aggregates
    aggid2nodes, aggid2procs = readin_aggregates("aggs_level%LEVEL_proc%PROC.out",numprocs,level)

    # vtk polygon for output
    aggpolygons = vtk.vtkAppendPolyData()

    # collect all aggregates
    for aggid,agg_nodes in aggid2nodes.iteritems():
      # build an aggregate
      if dimension==2:
	prepareDelaunayData(dimension, agg_nodes, global_nodecoords, aggpolygons, aggid, aggid2nodes, aggid2procs)
      else:
	prepareDelaunayData3d(dimension, agg_nodes, global_nodecoords, aggpolygons, aggid, aggid2nodes, aggid2procs)

    #aggpolygons.GetOutput().GetPointData().SetVectors(vtkDisplacementVector)
    #aggpolygons.Update()

    writer = vtk.vtkXMLPolyDataWriter()
    fname = "aggs"+str(level)+".vtp"
    writer.SetFileName(fname)
    writer.SetInput(aggpolygons.GetOutput())
    writer.Write()


    write_nodes_file("nodes"+str(level+1)+".txt",aggid2nodes,global_nodecoords,dimension)

    # increment number of multigrid levels that have been found in the files
    if no_multigridlevels < level:
      no_multigridlevels = level

    print "VTK Export for level " + str(level) + " finished...\r\n"

    level = level + 1


if __name__ == "__main__":
  sys.exit(main())
