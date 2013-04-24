#!/usr/bin/python
import sys
from subprocess import Popen, PIPE, STDOUT
from os import path

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
      x = float(line[0])
      y = float(line[1])
      if len(line) == 3:
	z = float(line[2])
	global_nodecoords.append((x,y,z))
      else:
	global_nodecoords.append((x,y))
  return global_nodecoords,dimension

#############################################
# create qconvex input file for aggregate
def create_qconvex_input(agg_globalnodeidx, global_nodecoords):

  local_nodeidx2global_nodeidx = {}

  no_of_aggnodes = len(agg_globalnodeidx)

  f = open("local_agg.inp","w")
  f.write(str(len(global_nodecoords[0])))	# write dimension
  f.write("\r\n")
  f.write(str(no_of_aggnodes)) # no of nodes for aggregate
  f.write("\r\n")

  for i in range(0,len(agg_globalnodeidx)):
    # store information in map local_nodeidx 2 global_nodeidx
    local_nodeidx2global_nodeidx[i] = agg_globalnodeidx[i]
    nodecoords = global_nodecoords[int(agg_globalnodeidx[i])]
    f.write(str(nodecoords[0]))
    f.write(" ")
    f.write(str(nodecoords[1]))
    if len(nodecoords)==3:
      f.write(" ")
      f.write(str(nodecoords[2]))
    f.write("\r\n")
  f.close()
  return local_nodeidx2global_nodeidx

##############################################
# run qconvex for some agg input
def run_qconvex_with_input(inputfile,local_nodeidx2global_nodeidx):
  f = file("local_agg.inp","r") # std inputfile is "local_agg.inp"
  p = Popen(['qconvex','o'], stdout=PIPE, stdin=f, stderr=STDOUT)
  qconvex_stdout = p.communicate()[0]
  f.close()
  outlines = qconvex_stdout.splitlines()

  # print qconvex_stdout
  if "ERRONEOUS FACET:" in qconvex_stdout:
    print "Error in qhull"
    lfacets = []
    lfacetslength = 0
    globalnodeidx_for_polygon = []
    return lfacets, lfacetslength, globalnodeidx_for_polygon

  if "qhull input error: not enough points" in qconvex_stdout:
    print "Error in qhull: not enough input points"
    lfacets = []
    lfacetslength = 0
    globalnodeidx_for_polygon = []
    return lfacets, lfacetslength, globalnodeidx_for_polygon

  if "qhull input error: input is less than 3-dimensional" in qconvex_stdout:
    print "Error in qhull: input is less than 3-dimensional"
    lfacets = []
    lfacetslength = 0
    globalnodeidx_for_polygon = []
    return lfacets, lfacetslength, globalnodeidx_for_polygon

  if "qhull input error: dimension 1 (first number) should be at least 2" in qconvex_stdout:
    print "Error in qhull: dimension is 1"
    lfacets = []
    lfacetslength = 0
    globalnodeidx_for_polygon = []
    return lfacets, lfacetslength, globalnodeidx_for_polygon

  # interpret output from qconvex
  noflines  = len(outlines) # number of lines of qconvex output
  dimension = outlines[0]   # line 0
  line1 = outlines[1].split(" ")

  number_of_points = int(line1[0])
  number_of_facets  = int(line1[1])

  lfacets = []    # list of all facets with all nodeidx stored
  for i in range(2+number_of_points,2+number_of_points+number_of_facets):
    # split face i
    linei = outlines[i].split(" ")
    no_nodes_in_facet_i = int(linei[0])

    gnodeidx = [] # list of local node ids for current face

    for j in range(1,no_nodes_in_facet_i+1):
      lnodeidxi = linei[j]
      gnodeidxi = local_nodeidx2global_nodeidx[int(lnodeidxi)]
      gnodeidx.append(gnodeidxi) # append global node idx to current faces of agg

    lfacets.append(gnodeidx) # add face to list of local faces of current aggregate

  # calculate some additional data for vtk file format
  lfacetslength = 0
  for i in range(0,len(lfacets)):
    lfacetslength = lfacetslength + len(lfacets[i]) + 1

  # check if 2d problem and all facets have length 2 (=lines)
  bIs2d = 1
  for i in range(0,len(lfacets)):
    if len(lfacets[i])==3:
      bIs2d = 0
      break

  # for 2d problems form a polygon
  globalnodeidx_for_polygon = []
  if bIs2d == 1:
    lfacet = lfacets[0]
    globalnodeidx_for_polygon.append(lfacet[0])
    globalnodeidx_for_polygon.append(lfacet[1])
    for i in range(1,len(lfacets)): # = number of missing nodes
      # look for next facet
      for fac in lfacets:
	if fac[0] == lfacet[1]:
	  globalnodeidx_for_polygon.append(fac[1]) # found facet
	  lfacet = fac
	  break

  # return list of faces with global node ids
  # and length of all face entries for this aggregate (for vtk output)
  return lfacets, lfacetslength, globalnodeidx_for_polygon

################################################################################
# CREATE VTK files
# routines that write aggregation information to VTK files

################################################################################
# create vtk file from global list of all nodes
# a list of faces of all aggregates
# and the overall facelength
# dimension is 2 or 3
# agg_polygons is a list of all aggregate polygons (only 2d)
# input: filename_prototype: string with prototype for output filename. %LEVEL is replaced by level
# input: level: integer with level index
# input: global_nodecoords: node coordinates of all nodes
# input: agg_facets: list of facets for all aggregates (put together from qconvex output)
# input: agg_facetlengths: overall length of facets information
# input: dimension (2 or 3) if dimension is 2 facets are lines, then we build the aggregate polygon
# input agg_polygons: list aggregate polygons (each aggregate = one polygon), empty if dimension == 3
def create_vtk_file(filename_prototype, level, global_nodecoords, agg_facets, agg_facetlengths, dimension, agg_polygons):
  # calculate number of all faces
  number_of_facets = 0
  for i in range(0, len(agg_facets)):
    number_of_facets = number_of_facets + len(agg_facets[i])

  # generate filename
  filename = filename_prototype
  filename = filename.replace("%LEVEL", str(level))

  # generate vtk file
  f = open(filename,"w")
  f.write("# vtk DataFile Version 1.0\r\n"); # print header
  f.write("Test example\r\n");
  f.write("ASCII\r\n\r\n");
  f.write("DATASET POLYDATA\r\n")
  f.write("POINTS ")                          # print point coords (->global?!)
  f.write(str(len(global_nodecoords)))
  f.write(" float\r\n")
  for i in range(0,len(global_nodecoords),1): # print point coordinates in list
    node_coordsi = global_nodecoords[i]
    f.write(str(node_coordsi[0]))
    f.write(" ")
    f.write(str(node_coordsi[1]))
    f.write(" ")
    if len(node_coordsi) == 3:
      f.write(str(node_coordsi[2]))
    else:
      f.write("0.0")
    f.write("\r\n")
  f.write("VERTICES 1 ")            # print vertices
  f.write(str(len(global_nodecoords)+1))
  f.write("\r\n")
  f.write(str(len(global_nodecoords)))
  for i in range(0,len(global_nodecoords),1): # print point coordinates in list
    f.write(" ")
    f.write(str(i))
  f.write("\r\n")
  if dimension == 2:
    f.write("LINES ")
  else:
    f.write("POLYGONS ")          # print faces
  f.write(str(number_of_facets))
  f.write(" ")
  f.write(str(agg_facetlengths))
  f.write("\r\n")
  # loop over all aggregates
  for agg in range(0, len(agg_facets)):
    cur_agg_facets = agg_facets[agg]
    for i in range(0,len(cur_agg_facets)): # print out all faces of current aggregate
      agg_facet_i = cur_agg_facets[i]
      f.write(str(len(agg_facet_i)))
      f.write(" ")
      for j in range(0,len(agg_facet_i)):
	f.write(str(agg_facet_i[j]))
	f.write(" ")
      f.write("\r\n")
  # for 2d problems write polygon
  if dimension == 2:
    poly_length = 0
    for agg in agg_polygons:
      poly_length = poly_length + len(agg) + 1
    f.write("POLYGONS ")
    f.write(str(len(agg_polygons)))
    f.write(" ")
    f.write(str(poly_length))
    f.write("\r\n")
    for agg in range(0, len(agg_polygons)):
      cur_agg_poly = agg_polygons[agg]
      f.write(str(len(cur_agg_poly)))
      for i in range(0,len(cur_agg_poly)):
	f.write(" ")
	f.write(str(cur_agg_poly[i]))
      f.write("\r\n")
  f.close();
  print "VTK file " + filename + " generated: OK"

def add_nodeinformation_to_vtk(filename_prototype, level, global_nodecoords,aggid2nodes,aggid2procs):
  nodeid2aggid = {} # empty map nodeid2aggid
  nodeid2proc = {}  # empty map nodeid2proc
  for aggid,nodes in aggid2nodes.iteritems():
    for n in nodes:
      nodeid2aggid[n] = aggid
      nodeid2proc[n]  = aggid2procs[aggid]

  # generate filename
  filename = filename_prototype
  filename = filename.replace("%LEVEL", str(level))

  f = open(filename,"a")
  f.write("POINT_DATA ")
  f.write(str(len(global_nodecoords)))
  f.write("\r\nSCALARS aggid float\r\nLOOKUP_TABLE t\r\n")
  for i in range(len(global_nodecoords)):
    if str(i) in nodeid2aggid:
      f.write(str(nodeid2aggid[str(i)]))
      f.write("\r\n")
    else:
      f.write("-1.0\r\n")
  f.write("SCALARS proc float\r\nLOOKUP_TABLE t2\r\n")
  for i in range(len(global_nodecoords)):
    if str(i) in nodeid2proc:
      f.write(str(nodeid2proc[str(i)][0]))
      f.write("\r\n")
    else:
      f.write("-1.0\r\n") # error
  f.close()
  print "VTK file " + filename + " filled with information: OK"

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
		aggid2procs[aggid] = [int(procid)]

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
def check_files_for_next_level(nextlevel, file_prototype):
  p = 0
  v = 1
  while v == 1:
    filename = file_prototype
    filename = filename.replace("%LEVEL", str(nextlevel))
    filename = filename.replace("%PROC",  str(p))
    if path.isfile(filename) == False:
      break
    p = p+1

  return p

###########
# MAIN routine
def main(argv=None):

  numprocs = 1
  level = 0  # startlevel

  # check if coords for finest level are available
  if path.isfile("nodes0.txt") == False:
    print "No \"nodes0.txt\" file found, exiting..."
    return

  v = 1
  while v == 1:
    numprocs = check_files_for_next_level(level,"aggs_level%LEVEL_proc%PROC.out")
    if numprocs <= 0:
      print "No aggregation information for level " + str(level) + " found, exiting..."
      break

    print "Aggregation information for " + str(numprocs) + " processors found"

    global_nodecoords,dimension = read_nodecoords_from_file("nodes"+str(level)+".txt")

    agg_faces = []  # empty list of all faces of all aggregates
    agg_facelengths = 0
    agg_polygons = [] # empty list with polygons for all aggregates (only 2d)

    aggid2nodes, aggid2procs = readin_aggregates("aggs_level%LEVEL_proc%PROC.out",numprocs,level)

    # collect all aggregates
    for aggid,agg_nodes in aggid2nodes.iteritems():
      # build an aggregate
      local_nodeidx2global_nodeidx = create_qconvex_input(agg_nodes, global_nodecoords)
      lfaces,lfaceslength,polygon = run_qconvex_with_input("local_agg.inp",local_nodeidx2global_nodeidx)
      agg_faces.append(lfaces)
      agg_facelengths = agg_facelengths + lfaceslength
      agg_polygons.append(polygon)

    # create vtk file
    create_vtk_file("output%LEVEL.vtk",level,global_nodecoords, agg_faces, agg_facelengths, dimension, agg_polygons)
    add_nodeinformation_to_vtk("output%LEVEL.vtk",level,global_nodecoords,aggid2nodes,aggid2procs)

    write_nodes_file("nodes"+str(level+1)+".txt",aggid2nodes,global_nodecoords,dimension)

    print "VTK Export for level " + str(level) + " finished...\r\n"

    level = level + 1

  #print check_files_for_next_level(level+1,numprocs,"aggs_level%LEVEL_proc%PROC.out")

if __name__ == "__main__":
  sys.exit(main())
