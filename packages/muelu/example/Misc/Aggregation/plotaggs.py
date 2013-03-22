#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import sys
import getopt
from Numeric import *
import numpy as np
import matplotlib as pltlib
import matplotlib.pyplot as plt
import convhull as ch
from matplotlib.path import Path
import matplotlib.patches as patches
import random

#############################################################################
# check, if line is a NODE line or not
def checkline(line):
	if line.find("NODE ") == 0:
		return 1,line
	else:
		return 0,line

################################################################################
# read in node file (from BACI dat files)
def read_nodes(filename):

	nodes = []
	for l in file(filename):
		line = l.strip()

		# filter out only NODE lines
		ret,line = checkline(line)
		if ret == 0:
			continue

		# line now contains a list of all tokens in that line
		line = line.split("\t")


		x = float(line[2])
		y = float(line[3])
		nodes.append((x,y))
	return nodes,len(nodes)

################################################################################
# helper function: uniquify list (non-order-preserving
def uniquify_list(seq):
    # Not order preserving
    keys = {}
    for e in seq:
        keys[e] = 1
    return keys.keys()

################################################################################
# add an aggregate patch to current axis ax
def add_aggregate_patch(ax,aggvertices,owningprocs):
	procs = uniquify_list(owningprocs)

	# translate tuple -> list
	curagg = []
	for node in range(0,len(aggvertices)):
		curagg.append(aggvertices[node])

	# prepare aggregate
	codes = [Path.MOVETO]
	for node in range(0,len(curagg)):
		codes.append(Path.LINETO)
	curagg.append(curagg[0])	# close aggregate (repeat first node)

	# create path from curagg
	path = Path(curagg, codes)

	# create patch for aggregate
	if len(procs) > 1:
	  fcol = (0.9,0.9,0.9,1)
	elif procs[0]==0:
	  fcol = (random.uniform(0.0,0.3),random.uniform(0.7,1.0),random.uniform(0.0,0.3))
	elif procs[0]==1:
	  fcol = (random.uniform(0.7,1.0),random.uniform(0.0,0.3),random.uniform(0.0,0.3))
	elif procs[0]==2:
	  fcol = (random.uniform(0.9,1.0),random.uniform(0.7,1.0),random.uniform(0.0,0.3))
	elif procs[0]==3:
	  fcol = (random.uniform(0.0,0.1),random.uniform(0.0,0.3),random.uniform(0.7,1.0))
	patch = patches.PathPatch(path, facecolor=fcol)
	ax.add_patch(patch)

################################################################################
# get rootnodes
# input: aggs = map: aggid -> list of nodeids in this agg
# output: list of rootnodes
# note: we calculate the "midpoint" of each aggregate
def get_rootnodes(aggid2gids,nodes):
	rootnodes = []
	for k in aggid2gids.keys():

		nodecoords = get_agg_coords(nodes,aggid2gids,k)


		firstpoint = nodecoords[0]
		x = 0.0
		y = 0.0
		for m in nodecoords:

			x = x + m[0]
			y = y + m[1]
		x = x/len(aggid2gids[k])
		y = y/len(aggid2gids[k])
		rootnodes.append((x,y))

	return rootnodes


################################################################################
# read in aggregation info from file
# input/output: gid2aggid: map for global row id -> global aggregation id
# input/output: gid2procid: map for global row id -> processor id
# input: procid: processor id (=value for gid2procid)
def read_aggregation_info(filename,gid2aggid,gid2procid,procid):
	lid2aggid = {}
	lid = 0
	for l in file(filename):
		line = l.strip()
		ret = line.partition(" ")
		mypid = ret[0]
		if mypid.find("MyPID")==0:	# skip first row
			continue
		mypid = int(mypid)
		line = ret[2]
		line = line.strip()
		ret = line.partition(" ")
		gid = int(ret[0])
		aggid = int(ret[2].strip())
		gid2aggid[gid]=aggid
		gid2procid[gid]=procid
		lid2aggid[lid]=aggid
		lid = lid + 1
	return lid2aggid, gid2aggid

################################################################################
# read in aggregation info from file
def readin_aggregates(procs,level):

  gid2aggid = {}
  gid2procid = {}

  for proc in range(0,procs):
    filename = "aggs_level" + str(level) + "_proc" + str(proc) + ".out"
    print "process ", filename
    [proclid2aggid,procgid2aggid] = read_aggregation_info(filename,gid2aggid,gid2procid,proc)
  return gid2aggid, gid2procid

################################################################################
# fill aggs with nodes
# input: gid2aggid: map from global row id to corresponding (global) aggregate id
# input: gid2procid: map from global row id to owning processors of corresponding agg
# output: aggid2gids = map of (global) aggregate id -> tuple of gids
def fill_aggs_with_gids(gid2aggid,gid2procid):

	# number of gids
	numgids = len(gid2aggid)

	################# setup aggs
	# aggid2gids is a map: globalaggid -> list of gids
	aggid2gids = {}
	aggid2procs = {}  # map: global agg id -> list of owning procs
	# loop over all nodes by nodeid
	for gid in range(0,numgids):
		globalaggid = gid2aggid[gid]

		# check if aggregate is already existing
		if globalaggid in aggid2gids:
			gidsinagg = aggid2gids[globalaggid]
			gidsinagg.append(gid)
			aggid2gids[globalaggid] = gidsinagg
			procsinagg = aggid2procs[globalaggid]
			procsinagg.append(gid2procid[gid])
			aggid2procs[globalaggid] = procsinagg
		else:
			# add new aggregate
			gidsinagg = [gid]
			aggid2gids[globalaggid] = gidsinagg
			procsinagg = [gid2procid[gid]]
			aggid2procs[globalaggid] = procsinagg
	return aggid2gids, aggid2procs

################################################################################
# get_agg_coords
# input: list of nodes (coordinates)
# input: map aggid2gids: map of local agg ids -> GIDs
# input: aggid: global aggregate id
# output: returns set of node coordinates for aggregate with global aggid
def get_agg_coords(nodes,aggid2gids,aggid):
  agg_gids = aggid2gids[aggid]
  nodeset = []
  for gid in range(0,len(agg_gids)):
    nodeset.append(nodes[agg_gids[gid]])
  return nodeset

################################################################################
# plot_aggregates
# input: ax = axis for plotting
# input: level = int with current level (0...)
# input: nodes = list of nodes (only used for finest level)
# input: numprocs = number of processors
# output: return rootnodes for next level
def plot_aggregates_level(ax,level,nodes,numprocs):

	print "plot aggregates level " + str(level)

	# read in aggs for current level
	[gid2aggid,gid2procid] = readin_aggregates(numprocs,level)

	# calculate new aggregates
	[aggid2gids,aggid2procs] = fill_aggs_with_gids(gid2aggid,gid2procid)

	# determine root node vertices from old aggregates
	#if level==0:
	#  rootnodes = nodes
	#else:
	rootnodes = get_rootnodes(aggid2gids, nodes)
	print "number of rootnodes " + str(len(rootnodes))

	# plot aggregates on current level
	for k in aggid2gids.keys():

	  nodecoords = get_agg_coords(nodes,aggid2gids,k)
	  owningprocs = aggid2procs[k]

	  if len(nodecoords) > 5:
	    arraydata = np.transpose(np.array(nodecoords))
	    convhulldata = ch.convex_hull(arraydata)
	    add_aggregate_patch(ax,convhulldata,owningprocs)
	  elif len(nodecoords) > 2:
	    pt1 = nodecoords[0]
	    pt2 = nodecoords[1]
	    pt3 = nodecoords[2]
	    nodecoords.append(((pt1[0]+pt2[0])/2,(pt1[1]+pt2[1])/2));
	    nodecoords.append(((pt1[0]+pt3[0])/2,(pt1[1]+pt3[1])/2));
	    nodecoords.append(((pt2[0]+pt3[0])/2,(pt2[1]+pt3[1])/2));
	    arraydata = np.transpose(np.array(nodecoords))
	    convhulldata = ch.convex_hull(arraydata)
	    add_aggregate_patch(ax,convhulldata,owningprocs)
	  elif len(nodecoords) == 2:
	    print "2 point aggregate"
	    plt.plot(nodecoords[0][0],nodecoords[0][1],'ro')
	    plt.plot(nodecoords[1][0],nodecoords[1][1],'bo')
	  else:
	    print "1 point aggregate"
	    plt.plot(nodecoords[0][0],nodecoords[0][1],'ro')

	return rootnodes

################################################################################
# plot aggregates
# input: number of levels
def plot_aggregates(nlevels,procs):
	print "plot aggregates"

	gid2aggid = {}
	curlevel = 0 # finest level

	# read in nodes
	nodes, numnodes = read_nodes("nodes" + str(curlevel) + ".txt")


	############# create new figure
	fig = plt.figure()

	width=nlevels
	coarseLevelnodes = nodes
	for i in range(0,width):

		ax = fig.add_subplot(1,width,i+1,aspect='equal')

		coarseLevelnodes = plot_aggregates_level(ax,i,coarseLevelnodes,procs)
	plt.show()


################################################################################
# MAIN routine
def main(argv=None):

	# plot aggregates
	# input: number of levels (e.g. 3)
	# input: number of processors
	plot_aggregates(3,2)

if __name__ == "__main__":
	sys.exit(main())



