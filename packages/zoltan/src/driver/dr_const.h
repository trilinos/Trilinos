/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifndef _DR_CONST_H
#define _DR_CONST_H

#include "zoltan.h"
#include <stdio.h>
#include <stdlib.h>

#define MIN(A,B)                (((A) < (B)) ? (A) : (B))

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/********  Trilinos Build Environment *******/
/* This block should only be executed for an Autotools build. */
#ifndef TRILINOS_NO_CONFIG_H

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for 
 * each package and need to
 * be undef'd here to avoid warnings when this file is included 
 * from another package.
 * KL 11/25/02
 */
#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

/* This file passes values from configure to the source code. */
#include "Zoltan_config.h"

#ifdef HAVE_PARMETIS
#define ZOLTAN_PARMETIS
#endif

#ifdef HAVE_METIS
#define ZOLTAN_METIS
#endif

#ifdef HAVE_SCOTCH
#define ZOLTAN_SCOTCH
#endif

#ifdef HAVE_OVIS
#define ZOLTAN_OVIS
#endif

#ifdef HAVE_GZIP
#define ZOLTAN_GZIP
#endif

#ifdef HAVE_PATOH
#define ZOLTAN_PATOH
#endif

#ifdef HAVE_PARKWAY
#define ZOLTAN_PARKWAY
#endif

#ifdef HAVE_NEMESIS_EXODUS
#define ZOLTAN_NEMESIS
#endif

#ifdef HAVE_PURIFY
#define ZOLTAN_PURIFY
#define strcasecmp Zoltan_strcasecmp
#define strncasecmp Zoltan_strncasecmp
#endif

#endif /* TRILINOS_NO_CONFIG_H */
/*****************************************************************************
 *  Definitions for the LB library driver program.
 *****************************************************************************/

#define DRIVER_NAME "zdrive"
#define VER_STR "1.0"


/* If it doesn't get defined in stdio.h then use this as a default */
#ifndef FILENAME_MAX
#define FILENAME_MAX    1024
#endif

/* Maximum defined by the NIST standard */
#define MATRIX_MARKET_MAX_LINE  1024

#define MAX_NP_ELEM	27 /* max nodes per element */
#define MAX_DIM		 3 /* max number of coordinate dimensions */
#define MAX_CPU_WGTS	10 /* max number of cpu weights */

enum LISTS {  /* NULL lists to pass to Zoltan_Migrate */
  NONE = 0,
  IMPORT_LISTS,
  EXPORT_LISTS
};

enum DATA_TYPE {
  MESH = 0,
  ZOLTAN_GRAPH,
  HYPERGRAPH
};

enum PARTITIONING_TYPE {
  HYPERGRAPH_PARTITIONING= 0,
  GRAPH_PARTITIONING,
  OBJECT_PARTITIONING
};

/*
 * Structure used to describe an element. Each processor will
 * allocate an array of these structures.
 */
struct Element_Description
{
  int      border;	/* set to 1 if this element is a border element */
  ZOLTAN_ID_TYPE globalID;	/* Global ID of this element, the local ID is the
			   position in the array of elements              */
  int      elem_blk;    /* element block number which this element is in */
  int      my_part;     /* Partition to which the element is assigned; 
                           default is processor number. */
  int      fixed_part;  /* Partition to which the element should be fixed;
                           this may not be equal to my_part before partitioning,
                           but a correct partitioning should assign the 
                           element to fixed_part.  fixed_part = -1 if the
                           element is not to be assigned to a fixed partition.*/
  int      perm_value;  /* Value for this element in permutation vector 
                           generated by Zoltan_Order.   
                           Default is -1 (no ordering done). */
  int      invperm_value;/* Value for this element in inverse permutation vector
                           generated by Zoltan_Order.   
                           Default is -1 (no ordering done). */
  float    cpu_wgt[MAX_CPU_WGTS]; /* the computational weight(s) associated with the elem */
  float    mem_wgt;	/* the memory weight associated with the elem */
  double   avg_coord[3];/* average coordinates (centroid) for the element */
  float  **coord;	/* array for the coordinates of the element.
                             for Nemesis meshes, nodal coordinates are stored;
                             for Chaco graphs with geometry, one set of coords
                                 is stored. */
  ZOLTAN_ID_TYPE *connect;	/* list of nodes that make up this element, the node
			   numbers in this list are global and not local    */
  ZOLTAN_ID_TYPE *adj;	/* list of adjacent elements .
                           For Nemesis input, the list is ordered by
                           side number, to encode side-number info needed to
                           rebuild communication maps.  Value -1 represents 
                           sides with no neighboring element (e.g., along mesh
                           boundaries).  Chaco doesn't have "sides," so the 
                           ordering is irrelevent for Chaco input. */
  int     *adj_proc;	/* list of processors for adjacent elements */
  int     *adj_blank;   /* NULL if not blanking, else 1/0 for blanked/not */
  float   *edge_wgt;	/* edge weights for adjacent elements */
  int      nadj;	/* number of entries in adj */
  int      adj_len;	/* allocated length of adj/adj_proc/edge_wgt arrays */
};
typedef struct Element_Description ELEM_INFO;
typedef struct Element_Description *ELEM_INFO_PTR;

/*
 * structure for general mesh information
 */
/*pStructure used to store information about the mesh */
struct Mesh_Description
{
  struct Zoltan_DD_Struct *dd;  /* Only used C and Fortran test driver */

  enum DATA_TYPE data_type;     /* Type of data stored in this data structure,
                                   based on input file type.
                                   Valid types are MESH, GRAPH, or HYPERGRAPH.*/
  int     vwgt_dim;             /* number of weights per element.            */
  int     ewgt_dim;             /* number of weights per graph edge.         */
  int     num_nodes;		/* number of nodes on this processor.        */
  int     num_elems;		/* number of elements on this processor.     */
  int     num_dims;		/* number of dimensions for the mesh         */
  int     num_el_blks;		/* number of element blocks in the mesh      */
  int     num_node_sets;	/* number of node sets in the mesh           */
  int     num_side_sets;	/* number of side sets in the mesh           */
  char  **eb_names;		/* element block element names               */
  int    *eb_etypes;            /* element block element types               */
  int    *eb_ids;		/* element block ids                         */
  ZOLTAN_ID_TYPE *eb_cnts;	/* number of elements in each element block  */
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
  ZOLTAN_ID_TYPE *ecmap_neighids;    /* elements ids of neighboring elements 
                                   for all elemental communication maps. 
                                   (global numbering)                        */
  int     elem_array_len;	/* length that the ELEM_INFO array is
				   allocated for. Need to know this when array
				   is not completely filled during migration */
  ELEM_INFO_PTR elements;       /* array of elements that are in the mesh.   */
  int     *blank;               /* 1 if my element is blanked, 0 if not      */
  int     blank_count;          /* number of my elements that are blanked    */
  ZOLTAN_ID_TYPE global_blank_count;   /* number of all elements that are blanked   */
  ZOLTAN_ID_TYPE gnhedges;             /* for hypergraphs, the number of global
                                   hyperedges.*/
  int     hewgt_dim;            /* for hypergraphs, the number of weights per
                                   hyperedge.                                */

                                /* The following fields usually supply the
                                   pins in the hyperedges (rows).  But if 
                                  "format" indicates columns instead of rows, 
                                   we store vertices and their pins here. */
  int     format;               /* rows (edges) or columns (vertices) */
  int     nhedges;              /* # local edges (if cols: # pin vertices) */
  ZOLTAN_ID_TYPE *hgid;        /* Global number for edges (or pin vertices),
                                   derived implicitly from order pins
                                   are read from file. Numbering is 0-based. */
  int    *hindex;               /* for hypergraphs, an entry for each 
                                   edge (or vertex), giving the starting index
                                   into hvertex for hyperedge (or vertex).*/ 
  ZOLTAN_ID_TYPE *hvertex;      /* row storage: global number for each pin
                                   vertex, col storage: global number for
                                   each pin hyperedge                     */
  int    *hvertex_proc;         /* row storage: array listing the processor 
                                   owning vertices in hvertex, or NULL if
                                   don't care.  col storage: NULL */

  int    heNumWgts;             /* number of edges for which we have weights */
  ZOLTAN_ID_TYPE *heWgtId;              /* global edge ID of the heNumWgts edges,
                                    if NULL it's the same as hgid            */
  float  *hewgts;               /* for hypergraphs, an array of hyperedge
                                   weights; size = hewgt_dim * heNumWgts;  */
  ZOLTAN_ID_TYPE  visible_nvtx;          /* #vertices to use, may be < num_elems */
  int    proc;        /* my rank, want to know if adj elements are on my proc */
};
typedef struct Mesh_Description  MESH_INFO;
typedef struct Mesh_Description *MESH_INFO_PTR;

/* Structure for the problem description. */
struct Parameter_Description {
  char Name[128];     /* parameter name */
  char Val[128];      /* parameter value */
  int  Index;         /* index for vector params */
};
typedef struct Parameter_Description Parameter_Pair;

struct Problem_Description
{
  char   method[32];                 /* this is the method string that will
                                        be passed unchanged to Zoltan        */
  int num_params;                    /* number of parameters read.           */
  Parameter_Pair *params;            /* parameter array to be passed to
                                        Zoltan.  Parameters are specified as
                                        pairs of strings:
                                          param_str = value_str              */
  char zoltanParams_file[FILENAME_MAX]; /* file name to get more
				       Zoltan parameters from separate
				       file (for hier support) */
};
typedef struct Problem_Description  PROB_INFO;
typedef struct Problem_Description *PROB_INFO_PTR;

/* Structure for driver flags for various test options. */
struct Test_Flags {
  int DDirectory;           /* Exercises data directories */
  int Local_Parts;          /* Sets NUM_LOCAL_PARTS parameter in various
                               ways. */
  int Fixed_Objects;        /* Registers functions for assigning fixed
                               objects; sets fixed_part within elements in
                               various way; sets RETURN_LISTS to 
                               "PARTITION ASSIGNMENTS" */
  int Drops;                /* Exercises point- and box-assign. */
  int RCB_Box;              /* Exercises Zoltan_RCB_Box. */
  int Multi_Callbacks;      /* Exercises list-based callback functions. */
  int Graph_Callbacks;      /* Register and test graph callbacks */
  int Hypergraph_Callbacks; /* Register and test hypergraph callbacks */
  int No_Global_Objects;    /* Test case where there are no objects on any process */
  int Gen_Files;            /* Exercise output file generation. */
  int Null_Lists;           /* Exercises null import or export lists to
                               Zoltan_Migrate. */
  float Dynamic_Weights;    /* Perturb weights between iterations. */
  float Dynamic_Graph;      /* Graph pertubation between iterations. */
  int Vtx_Inc;              /* Increment #vertices for each iteration. */
};

/* Structure for output flags for various types of output. */
struct Output_Flags {
  int Text;
  int Gnuplot;
  int Nemesis;
  int Plot_Partition;
  int Mesh_Info_File;
};

/* Global variables for driver */
extern int Debug_Driver;
extern int Debug_Chaco_Input;
extern int Number_Iterations;
extern int Driver_Action;
extern int Chaco_In_Assign_Inv;
extern struct Test_Flags Test;
extern struct Output_Flags Output;

extern double Total_Partition_Time;

#define DEBUG_TRACE_START(proc,yo) \
  if (((proc) == 0 && Debug_Driver > 1) || (Debug_Driver > 2))  \
    printf("%d DRIVER ENTERING %s\n", (proc), yo);
#define DEBUG_TRACE_END(proc,yo) \
  if (((proc) == 0 && Debug_Driver > 1) || (Debug_Driver > 2))  \
    printf("%d DRIVER LEAVING %s\n", (proc), yo);
#define DEBUG_TRACE_DETAIL(proc,yo,str) \
  if (Debug_Driver > 2) \
    printf("%d DRIVER %s: %s\n", proc,yo, str);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
#endif /* _DR_CONST_H */
