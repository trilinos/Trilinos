/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef _DR_CONST_H
#define _DR_CONST_H

#include <stdio.h>
#include <stdlib.h>
#include "zoltan.h"

/*****************************************************************************
 *  Definitions for the LB library driver program.
 *****************************************************************************/

#define DRIVER_NAME "zdrive"
#define VER_STR "1.0"


/* If it doesn't get defined in stdio.h then use this as a default */
#ifndef FILENAME_MAX
#define FILENAME_MAX    1024
#endif

#define MAX_NP_ELEM	27 /* max nodes per element */
#define MAX_DIM		 3 /* max number of coordinate dimensions */
#define MAX_CPU_WGTS	10 /* max number of cpu weights */

/*
 * Structure used to describe an element. Each processor will
 * allocate an array of these structures.
 */
struct Element_Description
{
  int      border;	/* set to 1 if this element is a border element */
  int      globalID;	/* Global ID of this element, the local ID is the
			   position in the array of elements              */
  int      elem_blk;    /* element block number which this element is in */
  float    cpu_wgt[MAX_CPU_WGTS]; /* the computational weight(s) associated with the elem */
  float    mem_wgt;	/* the memory weight associated with the elem */
  float  **coord;	/* array for the coordinates of the element.
                             for Nemesis meshes, nodal coordinates are stored;
                             for Chaco graphs with geometry, one set of coords
                                 is stored. */
  int     *connect;	/* list of nodes that make up this element, the node
			   numbers in this list are global and not local    */
  int     *adj;		/* list of adjacent elements .
                           For Nemesis input, the list is ordered by
                           side number, to encode side-number info needed to
                           rebuild communication maps.  Value -1 represents 
                           sides with no neighboring element (e.g., along mesh
                           boundaries).  Chaco doesn't have "sides," so the 
                           ordering is irrelevent for Chaco input. */
  int     *adj_proc;	/* list of processors for adjacent elements */
  float   *edge_wgt;	/* edge weights for adjacent elements */
  int      nadj;	/* number of entries in adj */
  int      adj_len;	/* allocated length of adj/adj_proc/edge_wgt arrays */
};
typedef struct Element_Description ELEM_INFO;
typedef struct Element_Description *ELEM_INFO_PTR;

/*
 * structure for general mesh information
 */
/* Structure used to store information about the mesh */
struct Mesh_Description
{
  int     num_nodes;		/* number of nodes on this processor         */
  int     num_elems;		/* number of elements on this processor      */
  int     num_dims;		/* number of dimensions for the mesh         */
  int     num_el_blks;		/* number of element blocks in the mesh      */
  int     num_node_sets;	/* number of node sets in the mesh           */
  int     num_side_sets;	/* number of side sets in the mesh           */
  char  **eb_names;		/* element block element names               */
  int    *eb_etypes;            /* element block element types               */
  int    *eb_ids;		/* element block ids                         */
  int    *eb_cnts;		/* number of elements in each element block  */
  int    *eb_nnodes;		/* number of nodes per element in each
				   element block                           
                                      for Nemesis meshes, this value depends
                                          on element type;
                                      for Chaco graphs, only one "node" per
                                          element.                           */
  int    *eb_nattrs;		/* number of attributes per element in each
				   element block                             */
  int     necmap;               /* number of elemental communication maps.   */
  int    *ecmap_id;             /* IDs of each elemental communication map.  */
  int    *ecmap_cnt;            /* number of elements in each elemental
                                   communication map.                        */
  int    *ecmap_elemids;        /* element ids of elements for all elemental
                                   communication maps. (local numbering)     */
  int    *ecmap_sideids;        /* side ids of elements for all elemental 
                                   communication maps.                       */
  int    *ecmap_neighids;       /* elements ids of neighboring elements 
                                   for all elemental communication maps. 
                                   (global numbering)                        */
  int     elem_array_len;	/* length that the ELEM_INFO array is
				   allocated for. Need to know this when array
				   is not completely filled during migration */
  ELEM_INFO_PTR elements;       /* array of elements that are in the mesh.   */

};
typedef struct Mesh_Description  MESH_INFO;
typedef struct Mesh_Description *MESH_INFO_PTR;

/* Structure for the problem description. */
typedef char Parameter_Pair[2][128]; /* typedef for parameter strings. 
                                        Parameters are specified as pairs
                                        of strings:
                                          param_str = value_str              */

struct Problem_Description
{
  char   method[32];                 /* this is the method string that will
                                        be passed unchanged to Zoltan        */
  int num_params;                    /* number of parameters read.           */
  Parameter_Pair *params;            /* parameter array to be passed to
                                        Zoltan.  Parameters are specified as
                                        pairs of strings:
                                          param_str = value_str              */
};
typedef struct Problem_Description  PROB_INFO;
typedef struct Problem_Description *PROB_INFO_PTR;

extern int Debug_Driver;
extern int DDirectory_Test;
extern int Gnuplot_Output;
extern int Number_Iterations;
extern int Debug_Chaco_Input;
extern int Chaco_In_Assign_Inv;

#define DEBUG_TRACE_START(proc,yo) \
  if (((proc) == 0 && Debug_Driver > 1) || (Debug_Driver > 2))  \
    printf("%d DRIVER ENTERING %s\n", (proc), yo);
#define DEBUG_TRACE_END(proc,yo) \
  if (((proc) == 0 && Debug_Driver > 1) || (Debug_Driver > 2))  \
    printf("%d DRIVER LEAVING %s\n", (proc), yo);
#define DEBUG_TRACE_DETAIL(proc,yo,str) \
  if (Debug_Driver > 2) \
    printf("%d DRIVER %s: %s\n", proc,yo, str);

#endif /* _DR_CONST_H */
