
/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_dr_const_id = "$Id$";
#endif

#ifndef _DR_CONST_H
#define _DR_CONST_H

#include "lbi_const.h"

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
#define MAX_DIM		3  /* max number of dimensions */

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
  float    cpu_wgt;	/* the computational weight associated with the elem */
  float    mem_wgt;	/* the memory weight associated with the elem */
  float  **coord;	/* array for the coordinates of the element.
                             for Nemesis meshes, nodal coordinates are stored;
                             for Chaco graphs with geometry, one set of coords
                                 is stored. */
  int     *connect;	/* list of nodes that make up this element, the node
			   numbers in this list are global and not local    */
  int     *adj;		/* list of adjacent elements */
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
  int     elem_array_len;	/* length that the ELEM_INFO array is
				   allocated for. Need to know this when array
				   is not completely filled during migration */
};
typedef struct Mesh_Description  MESH_INFO;
typedef struct Mesh_Description *MESH_INFO_PTR;

/*
 * global struct for mesh description
 * The Zoltan callback functions need both the element information struct
 * array and the mesh information struct. It is a lot easier to just pass
 * the element struct array as a data pointer, and have the mesh information
 * as a global variable.
 */
extern MESH_INFO Mesh;


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
  int    gen_graph;                  /* set if method being used needs graph */
  int    read_coord;                 /* set if method being used needs geom  */
};
typedef struct Problem_Description  PROB_INFO;
typedef struct Problem_Description *PROB_INFO_PTR;

extern void print_input_info(FILE *, int, PROB_INFO_PTR);

#endif /* _DR_CONST_H */
