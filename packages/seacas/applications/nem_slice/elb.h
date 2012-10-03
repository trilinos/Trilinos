/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef _EXOIILB_CONST_H_
#define _EXOIILB_CONST_H_

#include <string>
#include <vector>
#include <stdio.h>
#include <exodusII.h>
#include "elb_elem.h"

#define ELB_VERSION	"4.04"
#define UTIL_NAME	"nem_slice"
#define ELB_FALSE	0
#define ELB_TRUE	1

#define TOPTR(x) (x.empty() ? NULL : &x[0])

/* Macro for maximum value */
#ifndef MAX
#define MAX(x,y)  ((x > y) ? x : y)
#endif

/*
 * Constants for memory allocation of graph structures. The smaller these
 * values, the more efficient the code will be, memory wise. Larger values
 * will likely speed execution and prevent swap thrashing.
 */
#define SURND_ALLOC	8
#define ADJ_ALLOC	8
#define MEM_CHUNK_SIZE	16	/* Value MUST be >= 2 */
#define MEM_GROWTH      1.5

#define MAX_INP_LINE 10240

template <typename INT>
void vec_free(std::vector<INT> &V)
{ std::vector<INT>().swap(V);}

/* Prototype for timing function */
extern double get_time(void);

/* Structure used for the description of the machine for which the
 * load balance is to be constructed. */
struct Machine_Description
{
  int type;
  int num_dims;
  int dim[3];
  int num_boxes;     /* added for cluster type machines */
  int procs_per_box; /* added for cluster type machines, if there is only
                        one box, then this is the same as num_procs */
  int num_procs;

  Machine_Description() :
    type(-1), num_dims(-1), num_boxes(-1), procs_per_box(-1)
  { dim[0] = dim[1] = dim[2] = -1;}
};

/* Structure used for the description of what type of load balance is
 * to be performed. */
template <typename INT>
struct LB_Description
{
  int      type;
  int      refine;
  int      num_sects;
  int      cnctd_dom;
  int      outfile;
  std::string file;

  /* Calculated quantities */
  int     *vertex2proc;

  /* Nodal */
  std::vector<std::vector<INT> > int_nodes;
  std::vector<std::vector<INT> > bor_nodes;
  std::vector<std::vector<INT> > ext_nodes;
  std::vector<std::vector<INT> > ext_procs;

  /* Elemental */
  std::vector<std::vector<std::vector<INT> > > born_procs;
  std::vector<std::vector<INT> > int_elems;
  std::vector<std::vector<INT> > bor_elems;
  std::vector<std::vector<INT> > e_cmap_elems;
  std::vector<std::vector<INT> > e_cmap_sides;
  std::vector<std::vector<INT> > e_cmap_procs;
  std::vector<std::vector<INT> > e_cmap_neigh;

  LB_Description() :
    type(-1), refine(-1), num_sects(-1), cnctd_dom(-1), outfile(-1)
  {}
};

/* Structure for the problem description. */
struct Problem_Description
{
  int   type;
  int   read_coords;
  int   coarse_flag;
  int   alloc_graph;
  size_t num_vertices;
  int   vis_out;
  int   skip_checks;      /* put in to skip some error checks for some meshes  */
  int   face_adj;         /* true if using face definition of adjacencies      */
  int   partial_adj;      /* true if allowing partial (3/4) of nodes to */
                          /* determine adjancencies */
  int   global_mech;      /* true if looking for mechanisms in original mesh   */
  int   local_mech;       /* true if looking for mechanisms in subdivided mesh */
  int   find_cnt_domains; /* true if finding number of connected domains in a graph */
  int   mech_add_procs;   /* adds processors in cases of mechanisms       */
  int   dsd_add_procs;    /* adds processors in cases of disconnected subdomains */
  int   no_sph;
  char *groups;
  int  *group_no;
  int   num_groups;
  int   int64db;         /* integer types for output mesh database */
  int   int64api;        /* integer types for exodus api calls */

  Problem_Description() :
    type(-1), read_coords(-1), coarse_flag(-1), alloc_graph(-1), num_vertices(0), vis_out(-1),
    skip_checks(-1), face_adj(-1), partial_adj(0), global_mech(-1), local_mech(-1),
    find_cnt_domains(-1), mech_add_procs(-1), dsd_add_procs(-1), no_sph(-1),
    groups(NULL), group_no(NULL), num_groups(-1), int64db(0), int64api(0)
  {}
};

/* Structure for parameters needed for the Eigensolver in Chaco */
struct Solver_Description
{
  double tolerance;
  int    rqi_flag;
  int    vmax;
  
  Solver_Description() :
    tolerance(-1.0), rqi_flag(-1), vmax(-1)
  {}
};

/* Structure used to store information about the weighting scheme, if
 * any, that is to be used. */
template <typename INT>
struct Weight_Description
{
  int   type;  /* See weight type below for possible types */
  int   ow_read; /* 1 if element block settings overwrite exodus file read */

  std::string exo_filename;
  std::string exo_varname;
  
  int   exo_tindx;
  int   exo_vindx;

  /* Variable parameters */
  int    nvals;

  /* vectors to hold element block weights */
  std::vector<INT> elemblk; /* Id of element block */
  std::vector<INT> elemblk_wgt; /* Weight of that element block */

  /* vector to indicate if weight value has already been overwritten */
  std::vector<INT> ow;

  std::vector<int> vertices;
  std::vector<float> edges;

  Weight_Description<INT>() :
    type(-1), ow_read(0), exo_tindx(-1), exo_vindx(-1)
  {}
};

/* Structure used to store information about the FEM mesh */
template <typename INT>
struct Mesh_Description
{
  size_t     num_nodes;
  size_t     num_elems;
  size_t     num_dims;
  size_t     num_el_blks;
  INT       *eb_cnts;
  size_t     num_node_sets;
  size_t     num_side_sets;
  size_t     max_np_elem;
  size_t     ns_list_len;
  char    title[MAX_LINE_LENGTH+1];
  float  *coords;
  E_Type *elem_type;
  INT   **connect;
  
  Mesh_Description() :
    num_nodes(0), num_elems(0), num_dims(0), num_el_blks(0),
    eb_cnts(NULL), num_node_sets(0), num_side_sets(0),
    max_np_elem(0), ns_list_len(0), coords(NULL), elem_type(NULL),
    connect(NULL)
  {}
};

/* Structure for handling meshes with spheres */
struct Sphere_Info
{
  size_t num;
  int  *adjust;
  int  *begin;
  int  *end;

  Sphere_Info() :
    num(0), adjust(NULL), begin(NULL), end(NULL)
  {}
};

/* Structure used to store various information about the graph */
template <typename INT>
struct Graph_Description
{
  size_t nadj;
  size_t max_nsur;
  std::vector<INT> adj;
  std::vector<INT> start;
  std::vector<std::vector<INT> > sur_elem;
  Graph_Description<INT>() :
    nadj(0),  max_nsur(0)
  {}
};

/* Various constants */
#define NODAL 0
#define ELEMENTAL 1

#define UTIL_NAME "nem_slice"

/* Load balance types */
#define MULTIKL		0
#define SPECTRAL	1
#define INERTIAL	2
#define LINEAR 		3
#define RANDOM 		4
#define SCATTERED 	5
#define INFILE		6
#define KL_REFINE 	7
#define NO_REFINE 	8
#define NUM_SECTS	9
#define CNCT_DOM	10
#define OUTFILE		11
#define ZPINCH          12
#define BRICK           13
#define ZOLTAN_RCB      14
#define ZOLTAN_RIB      15
#define ZOLTAN_HSFC     16

/* Machine types */
#define MESH 		0
#define HCUBE		1
#define HYPERCUBE	2
#define CLUSTER		3

/* Solver options */
#define TOLER 		0
#define USE_RQI 	1
#define VMAX		2

/* ISSUES options */

#define LOCAL_ISSUES 0
#define GLOBAL_ISSUES 1

/* Weighting options */
/*
 * NOTE: the position of NO_WEIGHT, READ_EXO, EL_BLK, and EWGT_ON
 * should not be changed. These are the possible values for the
 * "type" variable in the Weight struct. They need to b 0, 1, 2, & 4
 * to allow bit masking for the type. The other variables are not
 * currently used in the type, but are needed since they appear
 * on the command line.
 */
#define NO_WEIGHT	0
#define READ_EXO	1
#define EL_BLK		2
#define VAR_INDX	3
#define EDGE_WGT	4
#define TIME_INDX	5
#define VAR_NAME	6

#endif /* _EXOIILB_CONST_H_ */
