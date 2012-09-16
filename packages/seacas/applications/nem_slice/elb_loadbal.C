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


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Functions contained in this file:
 *	generate_loadbal()
 *	generate_maps()
 *	nodal_dist()
 *	elemental_dist()
 *	ilog2i()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <assert.h>                     // for assert
#include <float.h>                      // for FLT_MAX
#include <limits.h>                     // for INT_MAX
#include <math.h>           /* Needed for ZPINCH_assign */
#include <stdio.h>                      // for printf, NULL, fprintf, etc
#include <stdlib.h>                     // for malloc, free, realloc, exit, etc
#include <string.h>                     // for strcat, strcpy

#include "chaco.h"                      // for input_assign, interface
#include "elb.h"                  // for LB_Description<INT>, etc
#include "elb_elem.h"             // for E_Type, get_elem_info, etc
#include "elb_err.h"              // for Gen_Error
#include "elb_graph.h"            // for generate_graph
#include "elb_groups.h"           // for get_group_info
#include "elb_loadbal.h"
#include "elb_util.h"             // for find_inter, etc

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#define M_PI_2 (M_PI/2.0)
#endif

#ifdef USE_ZOLTAN
#include <mpi.h>                        // for MPI_Finalize, etc
#include "zoltan.h"                     // for Zoltan_Set_Param, etc
#include "zoltan_types.h"               // for ZOLTAN_ID_PTR*, ZOLTAN_OK, etc
#endif

int is_hex(E_Type etype)
{
  if (etype == HEX8 || etype == HEX27 || etype == HEX20 || etype == HEXSHELL)
    return 1;
  else
    return 0;
}

int is_tet(E_Type etype)
{
  if (etype == TET4 || etype == TET10 || etype == TET8)
    return 1;
  else
    return 0;
}

int is_3d_element(E_Type etype)
{
  return (is_hex(etype) ||
	  is_tet(etype) ||
	  etype == WEDGE6 ||
	  etype == WEDGE15 ||
	  etype == WEDGE16 ||
	  etype == PYRAMID5 ||
	  etype == PYRAMID13 );
}

int ilog2i (size_t n)
{
  size_t i = 0;
  unsigned int n1 = n;
  while(n1 >>= 1) ++i;

  if((size_t)1 <<i != n)
    return(-1);
  else
    return(i);
}

namespace {
  template <typename INT>
  int nodal_dist(LB_Description<INT>*, Machine_Description*, Mesh_Description<INT>*, Graph_Description<INT>*);

  template <typename INT>
  int elemental_dist(LB_Description<INT>*, Machine_Description*, Mesh_Description<INT>*, Graph_Description<INT>*,
		     Problem_Description*);

  /* ZPINCH partitioning interface */
  int ZPINCH_assign(Machine_Description*, int, float *, float *, float *, 
			   int *);
  
  /* BRICK partitioning interface */
  int BRICK_assign(Machine_Description*, int, float *, float *, float *, 
			  int *);
  
#ifdef USE_ZOLTAN
  /* ZOLTAN_RCB partitioning interface */
  int ZOLTAN_assign(const char *, int, size_t, int *, float *, float *, float *,
			   int *, int, char **);
#endif
  void BALANCE_STATS(Machine_Description*, int *, size_t, int *);

  template <typename INT>
  int identify_mechanisms(Machine_Description* machine,
			  Problem_Description* problem,
			  Mesh_Description<INT>* mesh,
			  LB_Description<INT>* lb,
			  Graph_Description<INT>* graph,
			  int check_type);

  int extract_connected_lists( int nrow, const int* columns,  
			       const int* rows, int* list, 
			       int** list_ptr );
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function generate_loadbal() begins
 *----------------------------------------------------------------------------
 * This function calls Chaco to load balance the given problem.
 *****************************************************************************/
template int generate_loadbal(Machine_Description* machine, Problem_Description* problem,
			      Mesh_Description<int>* mesh, LB_Description<int>* lb,
			      Solver_Description* solve, Graph_Description<int>* graph,
			      Weight_Description<int>* weight, Sphere_Info* sphere,
			      int argc, char *argv[]);
template int generate_loadbal(Machine_Description* machine, Problem_Description* problem,
			      Mesh_Description<int64_t>* mesh, LB_Description<int64_t>* lb,
			      Solver_Description* solve, Graph_Description<int64_t>* graph,
			      Weight_Description<int64_t>* weight, Sphere_Info* sphere,
			      int argc, char *argv[]);


template <typename INT>
int generate_loadbal(Machine_Description* machine,
                     Problem_Description* problem,
                     Mesh_Description<INT>* mesh,
                     LB_Description<INT>* lb,
                     Solver_Description* solve,
                     Graph_Description<INT>* graph,
                     Weight_Description<INT>* weight,
                     Sphere_Info* sphere,
                     int argc,
                     char *argv[])
{
  const char   *assignfile = NULL;

  int     flag=0, arch=0, refine=0, num_level=0, totalproc=0, glob_method=0;
  int     start_proc = 0;

  size_t *tmp_start=NULL;
  INT    *tmp_adj=NULL;
  int    *tmp_vwgts=NULL;
  size_t tmp_nv;
  int    *nprocg=NULL, *nelemg=NULL, *nadjg=NULL;
  size_t  max_vtx, max_adj;
  int    group;
  INT    *elem_map=NULL;
  int     dim[1], tmpdim[3], tmp_arch, tmp_lev;
  int    *tmp_v2p=NULL;

  float  *x_ptr=NULL, *y_ptr=NULL, *z_ptr=NULL;
  float  *x_node_ptr=NULL, *y_node_ptr=NULL, *z_node_ptr=NULL;
  float  *x_elem_ptr=NULL, *y_elem_ptr=NULL, *z_elem_ptr=NULL;
  float  *tmp_x=NULL, *tmp_y=NULL, *tmp_z=NULL;
  float  *tmp_ewgts=NULL;

  long    seed=1;
  double *goal=NULL;
  double  time1, time2;
  FILE   *fp=NULL;

  /* Variables used in Chaco */
  extern int FREE_GRAPH;
  extern int CONNECTED_DOMAINS;
  extern int OUTPUT_ASSIGN;
  /*-----------------------------Execution Begins------------------------------*/
  
  tmpdim[0] = 0;
  tmpdim[1] = 0;
  tmpdim[2] = 0;

  /* if user requests, check for mechanisms in the original mesh  before
     working on the loadbalance
  */

  if(problem->type == ELEMENTAL && problem->global_mech == 1 && 
     problem->alloc_graph == ELB_TRUE)
    {
      printf("\n==============Looking For Global Issues=================\n");
      identify_mechanisms(machine, problem, mesh, lb, graph, GLOBAL_ISSUES);
      printf("============================================================\n");
      printf("\n");
    }


  /* Allocate the graph structure as Chaco expects it */
  if(problem->alloc_graph == ELB_TRUE) {
    /*
     * Increment the graph vertices to start a "1" instead of "0",
     * as expected by Chaco.
     */
    for(size_t cnt=0; cnt < graph->nadj; cnt++)
      graph->adj[cnt]++;
  }

  if(problem->read_coords == ELB_TRUE)
    {
      switch(mesh->num_dims) {
      case 3:
	x_node_ptr = (mesh->coords);
	y_node_ptr = (mesh->coords) + (mesh->num_nodes);
	z_node_ptr = (mesh->coords) + 2*(mesh->num_nodes);
	break;

      case 2:
	x_node_ptr = (mesh->coords);
	y_node_ptr = (mesh->coords) + (mesh->num_nodes);
	z_node_ptr = (float*)calloc(mesh->num_nodes, sizeof(float));
	break;

      case 1:
	x_node_ptr = (mesh->coords);
	y_node_ptr = (float*)calloc(mesh->num_nodes, sizeof(float));
	z_node_ptr = (float*)calloc(mesh->num_nodes, sizeof(float));
	break;
      }
    }
  else
    x_node_ptr = y_node_ptr = z_node_ptr = NULL;

  /* now set the pointers that are being sent to Chaco */
  x_ptr = x_node_ptr;
  y_ptr = y_node_ptr;
  z_ptr = z_node_ptr;

  /*
   * For an elemental decomposition using inertial, ZPINCH, BRICK
   * or ZOLTAN geometric load balancing,
   * I need to come up with coordinates for the elements.
   */
  if ((problem->type == ELEMENTAL) && 
      (lb->type == INERTIAL || lb->type == ZPINCH || lb->type == BRICK ||
       lb->type == ZOLTAN_RCB || lb->type == ZOLTAN_RIB || 
       lb->type == ZOLTAN_HSFC))
    {
      /* just check to make sure that there are vertices to get coords for */
      if (problem->num_vertices > 0) {
	/* allocate memory for element coordinates */
	x_elem_ptr = (float*)malloc(problem->num_vertices * sizeof(float));
	y_elem_ptr = (float*)malloc(problem->num_vertices * sizeof(float));
	z_elem_ptr = (float*)malloc(problem->num_vertices * sizeof(float));
	if (!(x_elem_ptr) || !(y_elem_ptr) || !(z_elem_ptr))
	  {
	    Gen_Error(0, "fatal: insufficient memory");
	    return 0;
	  }

	size_t cnt = 0;
	for (size_t ecnt=0; ecnt < mesh->num_elems; ecnt++)
	  {

	    if(mesh->elem_type[ecnt] != SPHERE || 
	       (mesh->elem_type[ecnt] == SPHERE && problem->no_sph == 1))
	      {
		/*
		 * for our purposes, the coordinate of the element will
		 * be the average of the coordinates of the nodes that make
		 * up that element
		 */
		x_elem_ptr[cnt] = 0.0;
		y_elem_ptr[cnt] = 0.0;
		z_elem_ptr[cnt] = 0.0;
		size_t nnodes = get_elem_info(NNODES, mesh->elem_type[cnt]);
		for (size_t ncnt=0; ncnt < nnodes; ncnt++)
		  {
		    x_elem_ptr[cnt] += x_ptr[(mesh->connect[ecnt][ncnt])];
		    y_elem_ptr[cnt] += y_ptr[(mesh->connect[ecnt][ncnt])];
		    z_elem_ptr[cnt] += z_ptr[(mesh->connect[ecnt][ncnt])];
		  }
		x_elem_ptr[cnt] = x_elem_ptr[cnt] / nnodes;
		y_elem_ptr[cnt] = y_elem_ptr[cnt] / nnodes;
		z_elem_ptr[cnt] = z_elem_ptr[cnt] / nnodes;
		cnt++;
	      }

	  } /* End "for (cnt=0; cnt < mesh->num_elem; cnt++)" */

	/* and use different pointers for Chaco */
	x_ptr = x_elem_ptr;
	y_ptr = y_elem_ptr;
	z_ptr = z_elem_ptr;

      } /* End "if (problem->num_vertices > 0)" */
    } /* End "if ((problem->type == ELEMENTAL) && 
	 (lb->type==INERTIAL||ZPINCH||BRICK||ZOLTAN))"*/


  /* Allocate memory for the vertex to processor vector */
  if(problem->type == ELEMENTAL)
    lb->vertex2proc = (int*)malloc(((problem->num_vertices)+sphere->num) *
				   sizeof(int));
  else
    lb->vertex2proc = (int*)malloc(problem->num_vertices * sizeof(int));

  if(!(lb->vertex2proc)) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }

  if(machine->type == HCUBE)
    {
      arch = 0;
      num_level = ilog2i((unsigned int) machine->procs_per_box);
    }
  if(machine->type == MESH)
    {
      arch = machine->num_dims;
      num_level = 0;
    }

  if(lb->refine == KL_REFINE)
    refine = 1;
  else
    refine = 2;

  if(lb->type == INFILE) {
    assignfile = lb->file.c_str();
    fp = fopen(assignfile, "r");
    if (fp == NULL) {
      Gen_Error(0, "fatal: could not open assignment file");
      return 0;
    }
  }

  if(lb->outfile) {
    assignfile = lb->file.c_str();
    OUTPUT_ASSIGN = 1;
  }

  switch(lb->type)
    {
    case MULTIKL:
      glob_method = 1;
      break;
    case SPECTRAL:
      glob_method = 2;
      break;
    case INERTIAL:
      glob_method = 3;
      break;
    case ZPINCH:
    case BRICK:
    case ZOLTAN_RCB:
    case ZOLTAN_RIB:
    case ZOLTAN_HSFC:
      glob_method = 999;  /* Chaco methods don't apply to ZPINCH, BRICK 
			     ZOLTAN_RCB, ZOLTAN_RIB, ZOLTAN_HSFC */
      break;
    case LINEAR:
      glob_method = 4;
      break;
    case RANDOM:
      glob_method = 5;
      break;
    case SCATTERED:
      glob_method = 6;
      break;
    case INFILE:
      glob_method = 7;
      break;
    }

  /* check if Chaco is supposed to make sure that the domains are connected */
  if (lb->cnctd_dom) CONNECTED_DOMAINS = 1;

  /*
   * check if the data needs to be sent to Chaco in groups, and
   * do some preliminary calculations and memory allocation
   *
   * Group or box designations for elements are stored in the
   * vertex2proc array as negative numbers. So, group 1 will
   * be designated as -1, etc.
   *
   * by default nloops = 1
   */
  size_t nloops = 1;
  if (problem->num_groups > 1) {
    /* allocate space to hold number of procs and elem per group */
    nprocg = (int *) (int*)malloc(2 * problem->num_groups * sizeof(int));
    if (!nprocg) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
    nelemg = nprocg + problem->num_groups;


    if (!get_group_info(machine, problem, mesh, graph, lb->vertex2proc, nprocg,
                        nelemg, &max_vtx, &max_adj)) {
      Gen_Error(0, "fatal: Error obtaining group information.");
      return 0;
    }

    nloops = problem->num_groups;
  }

  if (machine->num_boxes > 1) {
    /* allocate space to hold number of procs and elem per box */
    nprocg = (int *) malloc(3 * machine->num_boxes * sizeof(int));
    if (!nprocg) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
    nelemg = nprocg + machine->num_boxes;
    nadjg = nelemg + machine->num_boxes;

    /*
     * Call Chaco to get the initial breakdown of the verticies
     * onto the boxes. The vertex2proc array can then be used
     * to assign elements to groups that can be used in the
     * Chaco calls below.
     */
    tmp_arch = 1;
    tmp_lev = 0;
    dim[0] = machine->num_boxes;
    FREE_GRAPH = 0; /* Don't have Chaco to free the adjacency */

    printf("=======================Call Chaco===========================\n");
    time1 = get_time();
    if (lb->type == INFILE)
      flag = input_assign(fp, (char*)assignfile, problem->num_vertices,
                          lb->vertex2proc);
    if (lb->type == ZPINCH || lb->type == BRICK || 
        lb->type == ZOLTAN_RCB || lb->type == ZOLTAN_RIB || 
        lb->type == ZOLTAN_HSFC) {
      fprintf(stderr, "KDD -- ZPINCH, BRICK, ZOLTAN_RCB, ZOLTAN_RIB, and "
	      "ZOLTAN_HSFC not supported with num_boxes > 1.\n");
      fprintf(stderr, "KDD -- Contact Karen Devine, kddevin@sandia.gov.\n");
      exit(-1);
    }
    else
      flag = interface(problem->num_vertices, (size_t*)TOPTR(graph->start), (int*)TOPTR(graph->adj),
                       TOPTR(weight->vertices), TOPTR(weight->edges), x_ptr, y_ptr, z_ptr,
                       (char*)assignfile, (char *)NULL, lb->vertex2proc, tmp_arch,
                       tmp_lev, dim, goal, glob_method, refine,
                       solve->rqi_flag, solve->vmax, lb->num_sects,
                       solve->tolerance, seed);

    time2 = get_time();
    printf("============================================================\n");
    printf("Time in Chaco: %fs\n", time2-time1);

    if(flag != 0)
      {
	Gen_Error(0, "fatal: Chaco returned an error");
	return 0;
      }

    nloops = machine->num_boxes;
    for (size_t iloop = 0; iloop < nloops; iloop++) {
      nelemg[iloop] = nadjg[iloop] = 0;
      nprocg[iloop] = machine->procs_per_box;
    }

    for (size_t ecnt = 0; ecnt < problem->num_vertices; ecnt++) {
      nelemg[lb->vertex2proc[ecnt]]++;
      if (problem->alloc_graph == ELB_TRUE)
        nadjg[lb->vertex2proc[ecnt]] += graph->start[ecnt+1]
	  - graph->start[ecnt];

      /*
       * use negative numbers to specify the groups, in order to
       * avoid having a problem with group 0, add 1 to the group
       * number
       */
      lb->vertex2proc[ecnt] = -(lb->vertex2proc[ecnt] + 1);
    }

    max_vtx = 0;
    max_adj = 0;
    for (size_t iloop = 0; iloop < nloops; iloop++) {
      if (nelemg[iloop] > (int)max_vtx) max_vtx = nelemg[iloop];
      if (problem->alloc_graph == ELB_TRUE)
        if (nadjg[iloop] > (int)max_adj) max_adj = nadjg[iloop];
    }
  }

  if (nloops > 1) {

    /*
     * now allocate temporary arrays to hold information
     * to pass into Chaco
     */
    if (problem->alloc_graph == ELB_TRUE) {
      tmp_start = (size_t *) malloc((max_vtx + 1) * sizeof(size_t));
      tmp_adj = (INT *) malloc(max_adj * sizeof(INT));
      if (!tmp_start || !tmp_adj) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
      /* and need a group map for the elements */
      elem_map = (INT *) malloc(problem->num_vertices * sizeof(INT));
      if (!elem_map) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }
    if (!weight->vertices.empty()) {
      tmp_vwgts = (int *) malloc(max_adj * sizeof(int));
      if (!tmp_vwgts) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }
    if (!weight->edges.empty()) {
      tmp_ewgts = (float *) malloc(max_adj * sizeof(float));
      if (!tmp_ewgts) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }
    if (x_ptr != NULL) {
      tmp_x = (float *) malloc(max_vtx * sizeof(float));
      if (!tmp_x) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }
    if (y_ptr != NULL) {
      tmp_y = (float *) malloc(max_vtx * sizeof(float));
      if (!tmp_y) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }
    if (z_ptr != NULL) {
      tmp_z = (float *) malloc(max_vtx * sizeof(float));
      if (!tmp_z) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }
    tmp_v2p = (int *) malloc(max_vtx * sizeof(int));
    if (!tmp_v2p) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }

    start_proc = 0; /* counter to keep track of proccessors for each group */
  }

  /*
   * need to check to make sure that the mesh is not made up
   * entirely of SPHERE elements. If it is, don't call Chaco.
   * The spheres are linearly distributed after chaco is called.
   */

  if (sphere->num < mesh->num_elems) {

    /* now loop over the number of groups */
    for (size_t iloop = 0; iloop < nloops; iloop++) {

      /* group #'s 1 based */
      group = iloop + 1;
   
      /* Call chaco to generate the load balance */

      if (nloops > 1) {
        if (nelemg[iloop] <= 0) continue;

        tmp_nv = nelemg[iloop];


        if (problem->alloc_graph == ELB_TRUE) {
          /*
           * Build an element map that is local to this group. This
           * is necessary so that the adjacency list for each element in
           * this group can is correctly mapped. Elements that are in this
           * group are mapped 1:n and elements that are not in the group
           * are assigned -1.
           */
          INT elemp = 1;
          for (size_t cnt = 0; cnt < problem->num_vertices; cnt++) {
            if (group == -lb->vertex2proc[cnt])
              elem_map[cnt] = elemp++;
            else
              elem_map[cnt] = -1;
          }

          tmp_start[0] = 0;
        }

        size_t adjp = 0;
	size_t elemp = 0;

        /* fill the temporary vectors */
        for (size_t ecnt = 0; ecnt < problem->num_vertices; ecnt++) {
          if (group == -lb->vertex2proc[ecnt]) {
            if (problem->alloc_graph == ELB_TRUE) {
              for (int cnt = graph->start[ecnt]; cnt < graph->start[ecnt+1]; cnt++) {
                if (elem_map[graph->adj[cnt]-1] > 0) {
                  tmp_adj[adjp] = elem_map[graph->adj[cnt]-1];
                  if (!weight->edges.empty())
                    tmp_ewgts[adjp] = weight->edges[cnt];
                  adjp++;
                }
              }
              tmp_start[elemp+1] = adjp;
            }
            if (!weight->vertices.empty())
              tmp_vwgts[elemp] = weight->vertices[ecnt];

            if (x_ptr)
              tmp_x[elemp] = x_ptr[ecnt];
            if (y_ptr)
              tmp_y[elemp] = y_ptr[ecnt];
            if (z_ptr)
              tmp_z[elemp] = z_ptr[ecnt];

            elemp++; /* now increment the element # for the group */
          }
        }

        /*
         * set number of procs for this group
         * If this is a cluster machine, then machine-dim, arch, and
         * num_level are already set for each box.
         */
        if (problem->num_groups > 1) {
          if (machine->type == MESH) {

	    /*
	     * mesh and groups are in conflict.  the only way to resolve
	     * this is to assume a 1-d mesh. For example if group 1 requires
	     * 7 processors, it is impossible to apply this to a 2d mesh.
	     */

            tmpdim[0] = nprocg[iloop];
            tmpdim[1] = 1;
            tmpdim[2] = 1;
            totalproc = nprocg[iloop];
          }
          else
            num_level = nprocg[iloop];
	  totalproc = nprocg[iloop];
        }
        FREE_GRAPH = 0; /* Don't have Chaco to free the adjacency */
      }
      else {
        tmp_nv    = problem->num_vertices;
        tmp_start = &graph->start[0];
        tmp_adj   = &graph->adj[0];
        tmp_vwgts = &weight->vertices[0];
        tmp_ewgts = &weight->edges[0];
        tmp_x     = x_ptr;
        tmp_y     = y_ptr;
        tmp_z     = z_ptr;
        tmp_v2p   = lb->vertex2proc;

	if (problem->local_mech == 1)
	  FREE_GRAPH = 0;
	else
	  FREE_GRAPH = 0;/* Have Chaco to free the adjacency -- not anymore, using vectors*/	  

        for(int cnt = 0; cnt < machine->num_dims; cnt++) 
          tmpdim[cnt] = machine->dim[cnt];
        if (machine->type == MESH) {
          totalproc = tmpdim[0];
          if(tmpdim[1] != 0) totalproc *= tmpdim[1];
          if(tmpdim[2] != 0) totalproc *= tmpdim[2];
        }
        else
          totalproc = num_level;
      }


      if (lb->type == INFILE)
        flag = input_assign(fp, (char*)assignfile, tmp_nv, tmp_v2p);
      else if (lb->type == ZPINCH) {
        flag = ZPINCH_assign(machine, tmp_nv, tmp_x, tmp_y, tmp_z, tmp_v2p);
        BALANCE_STATS(machine, NULL, tmp_nv, tmp_v2p);
      }
      else if (lb->type == BRICK) {
        flag = BRICK_assign(machine, tmp_nv, tmp_x, tmp_y, tmp_z,tmp_v2p);
        BALANCE_STATS(machine, NULL, tmp_nv, tmp_v2p); 
      }
#ifdef USE_ZOLTAN
      else if (lb->type == ZOLTAN_RCB) {
        flag = ZOLTAN_assign("RCB", totalproc, tmp_nv, tmp_vwgts,
                             tmp_x, tmp_y, tmp_z, tmp_v2p, argc, argv);
        BALANCE_STATS(machine, tmp_vwgts, tmp_nv, tmp_v2p);
      }
      else if (lb->type == ZOLTAN_RIB) {
        flag = ZOLTAN_assign("RIB", totalproc, tmp_nv, tmp_vwgts,
                             tmp_x, tmp_y, tmp_z, tmp_v2p, argc, argv);
        BALANCE_STATS(machine, tmp_vwgts, tmp_nv, tmp_v2p);
      }
      else if (lb->type == ZOLTAN_HSFC) {
        flag = ZOLTAN_assign("HSFC", totalproc, tmp_nv, tmp_vwgts,
                             tmp_x, tmp_y, tmp_z, tmp_v2p, argc, argv);
        BALANCE_STATS(machine, tmp_vwgts, tmp_nv, tmp_v2p);
      }
#endif
      else if (lb->type == LINEAR) {
	fprintf(stderr, "Using internal linear decomposition\n");
	/* assign the elements to the processors linearly */
	size_t cnt  = mesh->num_elems / machine->num_procs;
  
	/* The left overs */
	int extra = (mesh->num_elems) % (machine->num_procs);

	size_t end = 0;
	for(int iproc=0; iproc < machine->num_procs; iproc++) {
	  size_t start = end;
	  end += cnt;
	  if (iproc < extra)
	    end++;

	  for(size_t ecnt=start; ecnt < end; ecnt++) {
	    lb->vertex2proc[ecnt] = iproc;
	  }
	}
	FREE_GRAPH = 0; /* Chaco not called, */

	for(size_t i=0; i < mesh->num_elems; i++) {
	  assert(lb->vertex2proc[i] < machine->num_procs);
	}
      }
      else if (lb->type == SCATTERED) {
	// Bad decomposition, useful for maximizing communication for testing algorithms...
	/* each element 'e' assigned to processor 'e' % 'num_proc' */
	
	for(size_t e=0; e < mesh->num_elems; e+=machine->num_procs) {
	  size_t iproc = 0;
	  size_t elem  = e;
	  size_t end   = elem + machine->num_procs;
	  if (end > mesh->num_elems)
	    end = mesh->num_elems;
	  while (elem < end) {
	    lb->vertex2proc[elem++] = iproc++;
	  }
	}
	FREE_GRAPH = 0; /* Chaco not called, */
      }
      else {
        printf("===================Call Chaco===========================\n");
        time1 = get_time();
        flag = interface(tmp_nv, (size_t*)tmp_start, (int*)tmp_adj,
                         tmp_vwgts, tmp_ewgts, tmp_x, tmp_y, tmp_z,
                         (char*)assignfile, (char *)NULL, tmp_v2p, arch,
                         num_level, tmpdim, goal, glob_method, refine,
                         solve->rqi_flag, solve->vmax, lb->num_sects,
                         solve->tolerance, seed);
        time2 = get_time();
        printf("========================================================\n");
        printf("Time in Chaco: %fs\n", time2-time1);
      }
   
      if(flag != 0)
	{
	  Gen_Error(0, "fatal: Partitioner returned an error");
	  return 0;
	}

      if (nloops > 1) {
        for (int ecnt = 0; ecnt < nelemg[iloop]; ecnt++)
          tmp_v2p[ecnt] += start_proc;
        start_proc += nprocg[iloop];
        /* copy the assignment data back into assignment */
        size_t elemp = 0;
        for (size_t ecnt = 0; ecnt < problem->num_vertices; ecnt++)
          if (-lb->vertex2proc[ecnt] == (int)group)
            lb->vertex2proc[ecnt] = tmp_v2p[elemp++];
      }

    } /* End: "for (iloop = 0; iloop < nloops; iloop++)" */
  } /* End: "if (sphere->num < mesh->num_elems)" */

  /* Free up coordinates if used */
  if(problem->read_coords == ELB_TRUE)
    {
      switch(mesh->num_dims)
	{
	case 2:
	  free(z_node_ptr);
	  break;
	case 1:
	  free(y_node_ptr);
	  free(z_node_ptr);
	}
    }

  /* free up element coordinate memory, if used */
  if (problem->type == ELEMENTAL)
    {
      if(lb->type == INERTIAL || lb->type == ZPINCH || lb->type == BRICK ||
	 lb->type == ZOLTAN_RCB || lb->type == ZOLTAN_RIB || 
	 lb->type == ZOLTAN_HSFC)
	{
	  if (x_elem_ptr) free(x_elem_ptr);
	  if (y_elem_ptr) free(y_elem_ptr);
	  if (z_elem_ptr) free(z_elem_ptr);
	}
    }

  /* free up memory for groups */
  if (nloops > 1) {
    if (problem->alloc_graph == ELB_TRUE) {
      free(tmp_start);
      free(tmp_adj);
      free (elem_map);
    }
    free (nprocg);
    free (tmp_v2p);
    if (tmp_vwgts) free (tmp_vwgts);
    if (tmp_ewgts) free (tmp_ewgts);
    if (tmp_x)     free (tmp_x);
    if (tmp_y)     free (tmp_y);
    if (tmp_z)     free (tmp_z);
    free (problem->group_no);
    free (mesh->eb_cnts);
    /* since Chaco didn't free the graph, need to do it here */
    vec_free(graph->start);
    vec_free(graph->adj);
    vec_free(weight->vertices);
    vec_free(weight->edges);
  }

  /*
   * If this is an elemental load balance and there are spheres present
   * then adjust lb->vertex2proc accordingly. The spheres are then
   * distributed linearly to processors.
   */
  if(problem->type == ELEMENTAL) {
    if(sphere->num > 0) {
      // If type == LINEAR || SCATTERED, then the spheres were handled above (unless the model is all spheres...)
      if ((lb->type != LINEAR && lb->type != SCATTERED) || sphere->num == mesh->num_elems) {

	if (sphere->num != mesh->num_elems) {
	  // If all spheres, don't need to open up space...
	  for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
	    /*
	     * First generate "holes" in the vertex2proc vector where the
	     * sphere assignments will be.
	     */
	    if(mesh->elem_type[ecnt] == SPHERE) {
	      for(size_t cnt=(problem->num_vertices); cnt > ecnt; cnt--)
		lb->vertex2proc[cnt] = lb->vertex2proc[cnt-1];
	      
	      lb->vertex2proc[ecnt] = -1;
	      (problem->num_vertices)++;
	    }
	  }
	}

	/* Now assign the spheres to the processors linearly */
	size_t cnt  = sphere->num / machine->num_procs;
  
	/* The left overs */
	flag = (sphere->num) % (machine->num_procs);
  
	size_t iproc = 0;
	size_t cnt2  = 0;
  
	for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
	  if(mesh->elem_type[ecnt] == SPHERE) {
	    lb->vertex2proc[ecnt] = iproc;
	    cnt2++;
  
	    /*
	     * If the processor ID is lower than the remainder then that
	     * processor gets an extra sphere
	     */
	    if((iproc+1) <= (size_t)flag) {
	      if(cnt2 >= (cnt+1)) {
		iproc++;
		cnt2 = 0;
	      }
	    }
	    else {
	      if(cnt2 >= cnt) {
		iproc++;
		cnt2 = 0;
	      }
	    }
	  }
	}
      }
    } /* End "if(sphere->num > 0)" */
  } /* End "if(lb->type == ELEMENTAL)" */

  if(problem->type == ELEMENTAL) {
    if(problem->local_mech == 1) {
      printf("\n==============Looking For Local Issues==================\n");

      if(problem->face_adj == 1 && problem->alloc_graph == ELB_TRUE) {
	identify_mechanisms(machine, problem, mesh, lb, graph, LOCAL_ISSUES);
      }
      else {
	/* need to free the bigger graph and create a face adjacent graph */
	int tmp_alloc_graph, tmp_adjacency;

	vec_free(graph->start);
	vec_free(graph->adj);
	vec_free(weight->vertices);
	vec_free(weight->edges);
	
	if(!graph->sur_elem.empty()) {
	  for(size_t cnt=0; cnt < mesh->num_nodes; cnt++)
	    vec_free(graph->sur_elem[cnt]);
	  vec_free(graph->sur_elem);
	}

	tmp_alloc_graph = problem->alloc_graph;
	tmp_adjacency   = problem->face_adj;

	problem->alloc_graph = ELB_TRUE;
	problem->face_adj    = 1;
      
	generate_graph(problem, mesh, graph, weight, sphere);

	/*
	 * adjacancy graph sent to identify_mechanisms must be 1 based
	 *
	 */

	for(size_t cnt = 0; cnt < graph->nadj; cnt++)
	  graph->adj[cnt]++;

	identify_mechanisms(machine, problem, mesh, lb, graph, LOCAL_ISSUES);

	problem->alloc_graph = tmp_alloc_graph;
	problem->face_adj    = tmp_adjacency;
      }

      printf("============================================================\n");
      printf("\n");
    }
  }

  /* since Chaco didn't free the graph, need to do it here */
  if (FREE_GRAPH == 0) {
    vec_free(graph->start);
    vec_free(graph->adj);
    vec_free(weight->vertices);
    vec_free(weight->edges);
  }
  return 1;

} /*---------------------- End generate_loadbal() ---------------------------*/


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function assign_fem() begins:
 *----------------------------------------------------------------------------
 * This function takes the load balance information generated by
 * generate_loadbal() and assigns the FEM quantities to processors based on
 * that load balance.
 *****************************************************************************/
template int generate_maps(Machine_Description*    machine,
			   Problem_Description*  problem,
			   Mesh_Description<int>*  mesh,
			   LB_Description<int>*    lb,
			   Graph_Description<int>* graph);

template int generate_maps(Machine_Description*    machine,
			   Problem_Description*  problem,
			   Mesh_Description<int64_t>*  mesh,
			   LB_Description<int64_t>*    lb,
			   Graph_Description<int64_t>* graph);


template <typename INT>
int generate_maps(
  Machine_Description*    machine,
  Problem_Description*  problem,
  Mesh_Description<INT>*  mesh,
  LB_Description<INT>*    lb,
  Graph_Description<INT>* graph)
{

/*-----------------------------Execution Begins------------------------------*/

  /* Generate the map for a nodal load balance */
  if(problem->type == NODAL)
  {
    /*
     * Generate the nodal and elemental distribution for a nodal
     * decomposition.
     */
    if(!nodal_dist(lb, machine, mesh, graph))
    {
      Gen_Error(0, "fatal: unable to find nodal distribution");
      return 0;
    }

  }
  else
  {
    /*
     * Generate the nodal and elemental distribution for an elemental
     * decomposition.
     */
    if(!elemental_dist(lb, machine, mesh, graph, problem))
    {
      Gen_Error(0, "fatal: unable to find elemental distribution");
      return 0;
    }
  }

  return 1;

} /*----------------------End generate_maps()----------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
namespace {
  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************
   *----------------------------------------------------------------------------
   * This function looks for 3D elements that are connected to other elements
   * by only a node or edge resulting in a boundary face.
   *****************************************************************************/
  template <typename INT>
  int identify_mechanisms(Machine_Description*    machine,
			  Problem_Description*  problem,
			  Mesh_Description<INT>*  mesh,
			  LB_Description<INT>*    lb,
			  Graph_Description<INT>* graph,
			  int            check_type)
  {
    INT   *pt_list;
    size_t nelem;
    INT   *hold_elem;
    size_t nhold;
    size_t count;
    size_t num_found = 0;

    int   *problems;
    int nhold2, nsides2;

    int *list_ptr;
    int end;

    /*
     * look for discontinuities in the graph
     *
     */

    if(problem->find_cnt_domains == 1) {
      if(check_type == GLOBAL_ISSUES) {
	size_t nrow = mesh->num_elems;
	int *rows = (int *) malloc((nrow+1)*sizeof(int));
	int *list = (int *) malloc(nrow*sizeof(int)) ;

	if (list == NULL || rows == NULL) {
	  Gen_Error(0, "fatal: insufficient memory");
	  return 0;
	}


	rows[0] = 1;
	for(size_t ecnt=1; ecnt < mesh->num_elems; ecnt++) {
	  ssize_t distance = graph->start[ecnt] - graph->start[ecnt-1] + 1;
	  rows[ecnt] = rows[ecnt-1] + distance;
	}
	rows[nrow] = graph->nadj + mesh->num_elems + 1;
	size_t nedges = rows[nrow] - 1;
	int *columns  = (int *) malloc(nedges*sizeof(int));

	if (columns == NULL) {
	  Gen_Error(0, "fatal: insufficient memory");
	  return 0;
	}

	size_t ki = 0;
	size_t kf = 0;
	for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
	  columns[kf++] = ecnt + 1;
	  size_t distance = rows[ecnt+1] - rows[ecnt] - 1;
	  for(size_t i = 0; i < distance; i++) {
	    columns[kf++] = graph->adj[ki++] + 1;
	  }
	}

	size_t components = extract_connected_lists(nrow, columns, rows, list, &list_ptr);

	if (components) {
    
	  printf("There are %lu connected components.\n",components);
	  for(size_t i=0; i <components; i++){
	    ki = (list_ptr)[i];
	    kf = (list_ptr)[i+1]-1;
	    size_t distance = kf - ki + 1;
	    printf("Connection %lu #elements %lu\n",i+1, distance);
	  }
	}

	free((char *)list_ptr);
	free(columns);
	free(rows);
	free(list);
      }

      if(check_type == LOCAL_ISSUES) {

	int *proc_cnt     = (int*)malloc(machine->num_procs * sizeof(int));
	INT *local_number = (INT*)malloc(mesh->num_elems * sizeof(INT));

	if((!(proc_cnt) || !(local_number))) {
	  Gen_Error(0, "fatal: insufficient memory");
	  return 0;
	}

	/*
	 * after running through chaco, the adjacency graph entries have
	 * been changed to 1 based
	 */
	for(int pcnt=0; pcnt < machine->num_procs; pcnt++) proc_cnt[pcnt] = 0;

	for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
	  int proc = lb->vertex2proc[ecnt];
	  proc_cnt[proc]++;
	  local_number[ecnt] = proc_cnt[proc];
	}

	int tmp_procs = machine->num_procs;
	for(int pcnt=0; pcnt < tmp_procs; pcnt++) {
	  if(proc_cnt[pcnt]) {
	    int nrow         = proc_cnt[pcnt];
	    int *rows         = (int *) malloc((nrow+1)*sizeof(int));
	    int *global_index = (int *) malloc((nrow)*sizeof(int));
	    int *list         = (int *) malloc(nrow*sizeof(int)) ;

	    if (list == NULL || rows == NULL || global_index == NULL) {
	      Gen_Error(0, "fatal: insufficient memory");
	      return 0;
	    }


	    rows[0] = 1;
	    size_t cnt = 0;
	    size_t tcnt = 0;
	    for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
	      int proc = lb->vertex2proc[ecnt];
	      assert(proc < machine->num_procs);
	      if(proc == pcnt) {
		if(ecnt < (mesh->num_elems -1)) {
		  end = graph->start[ecnt+1];
		}
		else {
		  end = graph->nadj;;
		}
		size_t distance = 1;
		for(int i = graph->start[ecnt]; i < end; i++) {
		  int proc2 = lb->vertex2proc[graph->adj[i]-1];
		  assert(proc2 < machine->num_procs);
		  if(proc2 == proc) distance++;
		}
		cnt++;
		rows[cnt] = rows[cnt-1] + distance;
		tcnt += distance;
	      }
	    }
	    rows[nrow] = tcnt + 1;
	    size_t nedges = rows[nrow] - 1;
	    int *columns  = (int *) malloc(nedges*sizeof(int));

	    if (columns == NULL) {
	      Gen_Error(0, "fatal: insufficient memory");
	      return 0;
	    }

	    {
	      size_t kf = 0;
	      size_t ki = 0;
	      for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
		int proc = lb->vertex2proc[ecnt];
		assert(proc < machine->num_procs);
		if(proc == pcnt) {
		  global_index[ki++] = ecnt;
		  columns[kf++]     = local_number[ecnt];
		  if(ecnt < (mesh->num_elems -1)) {
		    end = graph->start[ecnt+1];
		  }
		  else {
		    end = graph->nadj;
		  }
		  for(int i = graph->start[ecnt]; i < end; i++) {
		    int proc2 = lb->vertex2proc[graph->adj[i]-1];
		    assert(proc2 < machine->num_procs);
		    if(proc2 == proc) columns[kf++] = local_number[graph->adj[i]-1];
		  }
		}
	      }
	    }
	    
	    int components = extract_connected_lists(nrow, columns, rows, list, &list_ptr);

	    if (components > 0) {
    
	      printf("For Processor %d there are %d connected components.\n",pcnt, components);
	      for(int i=0; i <components; i++){
		int ki = (list_ptr)[i];
		int kf = (list_ptr)[i+1]-1;
		int distance = kf - ki + 1;
		printf("Connection %d #elements %d\n",i+1, distance);
	      }
	      if(problem->dsd_add_procs == 1) {
		for(int i=1; i <components; i++){
		  int ki = (list_ptr)[i];
		  int kf = (list_ptr)[i+1]-1;
		  for(int k=ki; k <=kf; k++){
		    lb->vertex2proc[global_index[list[k]-1]] = machine->num_procs;
		  }
		  machine->num_procs++;
		}
	      }
	    }

	    free((char *)list_ptr);
	    free(columns);
	    free(rows);
	    free(list);
	    free(global_index);
	  }
	}
	free(proc_cnt);
	free(local_number);

	if(tmp_procs != machine->num_procs) {
	  printf("\n!!! Processor count increased to %d processors\n", machine->num_procs);
	  printf("!!! in order to make connected subdomains\n\n");
	} 
      }
    }


    /*
     * determine local element numbers for diagnostic output
     *
     */

    if(problem->global_mech == 1 || problem->local_mech == 1) {

      pt_list = (INT*)malloc(graph->max_nsur * sizeof(INT));
      hold_elem = (INT*)malloc(graph->max_nsur * sizeof(INT));
      problems  = (int*)malloc(mesh->num_nodes * sizeof(int));

      if(!(pt_list) || !(hold_elem) || !(problems))
	{
	  Gen_Error(0, "fatal: insufficient memory");
	  return 0;
	}

      int *proc_cnt     = (int*)malloc(machine->num_procs * sizeof(int));
      INT *local_number = (INT*)malloc(mesh->num_elems * sizeof(INT));

      if((!(proc_cnt) || !(local_number))) {
	  Gen_Error(0, "fatal: insufficient memory");
	  return 0;
	}

      if(check_type == LOCAL_ISSUES) {
	for(int pcnt=0; pcnt < machine->num_procs; pcnt++)
	  proc_cnt[pcnt] = 0;
	for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
	  int proc = lb->vertex2proc[ecnt];
	  assert(proc < machine->num_procs);
	  proc_cnt[proc]++;
	  local_number[ecnt] = proc_cnt[proc];
	}
      }


      for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
	E_Type etype    = mesh->elem_type[ecnt];

	/* need to check for volume elements */
	if (etype == HEX8     || etype == HEXSHELL || etype == HEX20 ||
	    etype == TET4     || etype == TET10     || etype == WEDGE6 ||
	    etype == WEDGE15   || etype == WEDGE16   || etype == PYRAMID5 || 
	    etype == PYRAMID13 || etype == TET8 ) {

	  int nnodes = get_elem_info(NNODES, mesh->elem_type[ecnt]);

	  for(int ncnt=0; ncnt < nnodes; ncnt++) {
	    size_t node = mesh->connect[ecnt][ncnt];
	    problems[node] = 0;
	  }

	  int nsides   = get_elem_info(NSIDES, etype);
	  int proc = 0;
	  if(check_type == LOCAL_ISSUES)
	    proc = lb->vertex2proc[ecnt];
	  else
	    proc = 0;
	  assert(proc < machine->num_procs);

	  /* check each side of this element */

	  INT side_nodes[MAX_SIDE_NODES];
	  INT side_nodes2[MAX_SIDE_NODES];
	  for (int nscnt = 0; nscnt < nsides; nscnt++) {

	    /* get the list of nodes on this side set */
	    int side_cnt = ss_to_node_list(etype, mesh->connect[ecnt], (nscnt+1),
					   side_nodes);

	    for(int ncnt = 0; ncnt < side_cnt; ncnt++) {

	      size_t node  = side_nodes[ncnt];
	      nhold = graph->sur_elem[node].size();

	      /* 
	       * look for the following cases
	       * 1) bar elements connected to a volume element 
	       * 2) shell element connected to a edge or point of a volume
	       *     element and both its faces are boundaries
	       *
	       */
	      for (size_t pcnt = 0; pcnt < nhold; pcnt++) {

		  size_t el2     = graph->sur_elem[node][pcnt];
		  E_Type etype2  = mesh->elem_type[el2];

		  int proc2 = 0;
		  if(check_type == LOCAL_ISSUES) proc2 = lb->vertex2proc[el2];
		  assert(proc2 < machine->num_procs);

		  if(ecnt != el2 && proc == proc2) {
		    if(etype2 == BAR2 || etype2 == BAR3 || etype2 == SHELL2 || etype2 == SHELL3) {
		      problems[node] = el2+1;
		    }
		    else if(etype2 == SHELL4 || etype2 == SHELL8 ||
			    etype2 == TSHELL3 || etype2 == TSHELL6) { 
		      /*
		       * look for an element connect to one of the shells faces (not edges)
		       * one can look at elements connected to 3 of the nodes - cannot use
		       * diagonals due to triangular shells
		       */

		      nsides2   = get_elem_info(NSIDES, etype2);

		      count = 0;
		      for (int cnt = 0; cnt < nsides2; cnt++) {
                 
			ss_to_node_list(etype2, mesh->connect[el2], (cnt+1),
					side_nodes2);
                
			nhold2 = find_inter(TOPTR(graph->sur_elem[side_nodes2[0]]), 
					    TOPTR(graph->sur_elem[side_nodes2[1]]),
					    graph->sur_elem[side_nodes2[0]].size(),
					    graph->sur_elem[side_nodes2[1]].size(),
					    pt_list);
        
			for(int i = 0; i < nhold2; i++) 
			  hold_elem[i] = graph->sur_elem[side_nodes2[0]][pt_list[i]];

			nelem = find_inter(hold_elem,
					   TOPTR(graph->sur_elem[side_nodes2[2]]),
					   nhold2,
					   graph->sur_elem[side_nodes2[2]].size(),
					   pt_list);

			if(nelem >= 1) {
			  count++;
			  break;
			}
		      }
      
		      /* at this point, a shell was connnect to this volume element.
		       * if count, then the shell has a element connect to one of
		       * its faces and thus this is not a mechanism.  If !count,
		       * then this is a mechanism.
		       */

		      if(!count) problems[node] = el2+1; 
		    }
		  }
		}
	    }
	  }


	  for(int ncnt=0; ncnt < nnodes; ncnt++) {
	    size_t node = mesh->connect[ecnt][ncnt];
	    if(problems[node]) {
	      size_t el2    = problems[node]-1;
	      E_Type etype2 = mesh->elem_type[el2];

	      if(check_type == LOCAL_ISSUES) 
		{
		  printf("WARNING: On Processor %d Local Element %lu (%s) has a mechanism through Global Node %lu with Local Element %lu (%s)\n", 
			 proc, 
			 (size_t)local_number[ecnt], 
			 elem_name_from_enum(etype), 
			 (size_t)node, 
			 (size_t)local_number[el2], 
			 elem_name_from_enum(etype2)); 
		  if(problem->mech_add_procs == 1) lb->vertex2proc[el2] = machine->num_procs;
		}
	      else 
		{
		  printf("WARNING: Element %lu (%s) has a mechanism through Node %lu with Element %lu (%s)\n", 
			 (size_t)ecnt+1, elem_name_from_enum(etype), (size_t)node, (size_t)el2+1, elem_name_from_enum(etype2)); 
		}
	      num_found++;
	    }
	  }
	}
      }

      free(pt_list);
      free(hold_elem);
      free(problems);
      free(proc_cnt);
      free(local_number);

      if(num_found) {
	printf("Total mechanisms found = %lu\n", num_found);
	if(check_type == LOCAL_ISSUES) {
	  if(problem->mech_add_procs == 1) {
	    machine->num_procs++;
	    printf("\n!!! Processor count increased to %d processors\n", machine->num_procs);
	    printf("!!! to move mechanisms to another processor\n\n");
	  }
	}
      }
      else
	printf("NO mechanisms found\n");
    }

    return 1;
  }

  template <typename INT>
  int nodal_dist(LB_Description<INT>* lb,
		 Machine_Description* machine,
		 Mesh_Description<INT>* mesh,
		 Graph_Description<INT>* graph)
  {
    size_t ncnt, i;
    int    proc, proc_n;
    int    internal, flag, elem, node, nnodes;
    E_Type etype;
    double time1, time2;

    /*-----------------------------Execution Begins------------------------------*/
    /* Allocate memory */
    time1 = get_time();
    lb->int_nodes.resize(machine->num_procs);
    lb->bor_nodes.resize(machine->num_procs);
    lb->ext_nodes.resize(machine->num_procs);
    lb->int_elems.resize(machine->num_procs);
    lb->ext_procs.resize(machine->num_procs);

    time2 = get_time();
    printf("Allocation time: %fs\n", time2-time1);

    /* Find the internal, border and external nodes */
    time1 = get_time();
    for(ncnt=0; ncnt < mesh->num_nodes; ncnt++) {
      proc = lb->vertex2proc[ncnt];
      assert(proc < machine->num_procs);
      internal = 1;
      flag     = 0;
      for(size_t ecnt=0; ecnt < graph->sur_elem[ncnt].size(); ecnt++) {
	elem   = graph->sur_elem[ncnt][ecnt];
	etype  = mesh->elem_type[elem];
	nnodes = get_elem_info(NNODES, etype);
	for(i=0; i < (size_t)nnodes; i++) {
	  proc_n = lb->vertex2proc[mesh->connect[elem][i]];
	  assert(proc_n < machine->num_procs);
	  if(proc_n != proc) {
	    /* "ncnt" is a border node and is an external node to proc_n */
	    internal = 0;
	    if(!flag) {
	      flag = 1;
	      lb->bor_nodes[proc].push_back(ncnt);
	    }

	    /*
	    ** to make sure that this node has not already been put
	    ** in the external node list for proc_n I need to check
	    ** only the last element in the current list
	    */
	    if((!lb->ext_nodes[proc_n].empty()) ||
	       (lb->ext_nodes[proc_n][lb->ext_nodes[proc_n].size()-1] != (int)ncnt)) {
	      lb->ext_nodes[proc_n].push_back(ncnt);
	      lb->ext_procs[proc_n].push_back(proc);
	    }
	  }
	} /* End "for(i=0; i < nnodes; i++)" */
      } /* End "for(ecnt=0; ecnt < graph->nsur_elem[ncnt]; ecnt++)" */

      if(internal) {
	/* "ncnt" is an internal node */
	lb->int_nodes[proc].push_back(ncnt);
      }
    } /* End "for(ncnt=0; ncnt < mesh->num_nodes; ncnt++)" */
    time2 = get_time();
    printf("Time for nodal categorization: %fs\n", time2-time1);

    /* Find the internal elements */
    time1 = get_time();
    for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
      etype  = mesh->elem_type[ecnt];
      nnodes = get_elem_info(NNODES, etype);
      for(ncnt=0; ncnt < (size_t)nnodes; ncnt++) {
	  node = mesh->connect[ecnt][ncnt];
	  proc = lb->vertex2proc[node];
	  assert(proc < machine->num_procs);
	  /*
	  ** since the outer loop is on the elements, I don't need to
	  ** search over the entire list to find out if this element is
	  ** already in it. If the element is in the processors list,
	  ** then it must be the last element.
	  */
	  if((lb->int_elems[proc].empty()) ||
	     (lb->int_elems[proc][lb->int_elems[proc].size()-1] != (INT)ecnt)) {
	    lb->int_elems[proc].push_back(ecnt);
	  }
	}
    }
    time2 = get_time();
    printf("Elemental categorization: %fs\n", time2-time1);
    return 1;
  } /*-----------------------------End nodal_dist()----------------------------*/

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/
  template <typename INT>
  int elemental_dist(LB_Description<INT>* lb,
		     Machine_Description* machine,
		     Mesh_Description<INT>* mesh,
		     Graph_Description<INT>* graph,
		     Problem_Description* problem)
  {
    int    sid;
    int    hflag1, hflag2, tflag1, tflag2;
    int    dflag;

    INT   *pt_list, nelem;
    INT   *hold_elem;
    INT    side_nodes[MAX_SIDE_NODES], mirror_nodes[MAX_SIDE_NODES];
    int    side_cnt;

    double time1, time2;

    char   cmesg[256], tmpstr[80];

    /*-----------------------------Execution Begins------------------------------*/

    /* Allocate memory */
    lb->int_nodes.resize(machine->num_procs);
    lb->bor_nodes.resize(machine->num_procs);
    lb->ext_nodes.resize(machine->num_procs); /* Not used in elemental dist */

    lb->int_elems.resize(machine->num_procs);
    lb->bor_elems.resize(machine->num_procs);

    lb->ext_procs.resize(machine->num_procs);
    lb->born_procs.resize(machine->num_procs);

    lb->e_cmap_elems.resize(machine->num_procs);
    lb->e_cmap_sides.resize(machine->num_procs);
    lb->e_cmap_procs.resize(machine->num_procs);
    lb->e_cmap_neigh.resize(machine->num_procs);

    /* allocate space to hold info about surounding elements */
    pt_list   = (INT*)malloc(graph->max_nsur * sizeof(INT));
    hold_elem = (INT*)malloc(graph->max_nsur * sizeof(INT));
    if(!(pt_list) || !(hold_elem)) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }

    /* Find the internal and border elements */
    time1 = get_time();

    for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
      int proc = lb->vertex2proc[ecnt];
      assert(proc < machine->num_procs);
      bool internal = true;
      int flag      = 0;
      E_Type etype  = mesh->elem_type[ecnt];
      int dim1      = get_elem_info(NDIM, etype);

      /* need to check for hex's or tet's */
      hflag1 = is_hex(etype);

      /* a TET10 cannot connect to a HEX */
      tflag1 = is_tet(etype);

      int nsides   = get_elem_info(NSIDES, etype);

      /* check each side of this element */
      for (int nscnt = 0; nscnt < nsides; nscnt++) {

	/* get the list of nodes on this element side */
	side_cnt = ss_to_node_list(etype, mesh->connect[ecnt], (nscnt+1),
				   side_nodes);

	/*
	 * now determine how many side set nodes are needed to 
	 * determine if there is an element connected to this side.
	 *
	 * 2-D - need two nodes, so find one intersection
	 * 3-D - need three nodes, so find two intersections
	 * NOTE: must check to make sure that this number is not
	 *       larger than the number of nodes on the sides (ie - SHELL).
	 */

	int nnodes = mesh->num_dims;
	if (side_cnt < nnodes)   nnodes = side_cnt;
	nnodes--;    /* decrement to find the number of intersections needed */


	nelem = 0;   /* reset this in case no intersections are needed */

	/*
	 * need to handle hex's differently because of
	 * the tet/hex combination
	 */

	if (!hflag1) { /* Not a hex */

	  /* ignore degenerate bars */

	  if (!((etype == BAR2 || etype == SHELL2) && side_nodes[0] == side_nodes[1])) {

	    size_t nhold = graph->sur_elem[side_nodes[0]].size();
	    for (size_t ncnt = 0; ncnt < nhold; ncnt++)
	      hold_elem[ncnt] = graph->sur_elem[side_nodes[0]][ncnt];

	    for (int ncnt = 0; ncnt < nnodes; ncnt++) {
	      /* Find elements connnected to both node '0' and node 'ncnt+1' */
	      nelem = find_inter(hold_elem,
				 TOPTR(graph->sur_elem[side_nodes[(ncnt+1)]]),
				 nhold, graph->sur_elem[side_nodes[(ncnt+1)]].size(),
				 pt_list);

	      if (nelem < 2)
		break;
	      else {
		nhold = nelem;
		for (int ncnt2 = 0; ncnt2 < nelem; ncnt2++)
		  hold_elem[ncnt2] = hold_elem[pt_list[ncnt2]];
	      }
	    }
	  }
	  else {
	    printf("WARNING: Element = %lu is a DEGENERATE BAR\n", ecnt+1);
	  }
	}
	else { /* Is a hex */

	  /*
	   * Hex faces are fairly complicated now. There are two
	   * exceptions to the standard case:
	   *   1. it is a degenerate hex (mimics a wedge). this has
	   *      two special faces, the first is a triangle, and the
	   *      second is a 2d line
	   *   2. two tets are connected to this hex face 
	   */

	  /* first need to check for a degenerate element */
	  dflag = 0;
	  if (side_nodes[0] == side_nodes[1] || side_nodes[0] == side_nodes[3])
	    dflag++;
	  if (side_nodes[2] == side_nodes[1] || side_nodes[2] == side_nodes[3])
	    dflag++;

	  /*
	   * if both flags are set, then this face is the 2d line,
	   * and should be ignored with respect to elemental
	   * communication maps
	   */
	  if (dflag == 2) {
	    nelem = 1;

	  } else {
	    /*
	     * In order to check for two tets connected to this face,
	     * check the intersection of opposite corners of this face.
	     * Both tets should show up in the intersection of one of the
	     * sets of opposite corners (nothing should show up in the
	     * other).
	     */

	    /*
	     * Initial check is side nodes 0 and 2 which are
	     * diagonally opposite
	     */
	    int inode = 0;
	    int node  = 2; 
	    size_t nhold = 0;
	    for (int ncnt = 0; ncnt < nnodes; ncnt++) {
	      /* Find elements connnected to both node 'inode' and node 'node' */
	      nelem = find_inter(TOPTR(graph->sur_elem[side_nodes[inode]]),
				 TOPTR(graph->sur_elem[side_nodes[node]]),
				 graph->sur_elem[side_nodes[inode]].size(),
				 graph->sur_elem[side_nodes[node]].size(),
				 pt_list);

	      if (nelem > 1) {
		if (ncnt == 0) {
		  nhold = nelem;
		  for (int ncnt2 = 0; ncnt2 < nelem; ncnt2++)
		    hold_elem[ncnt2] = graph->sur_elem[side_nodes[inode]][pt_list[ncnt2]];

		  if (dflag) {
		    /*
		     * in this case, need to get an intersection with
		     * another (unique) point since nodes 0 and 2
		     * may represent an edge and not the diagonal
		     */
		    if (side_nodes[1] != side_nodes[0] &&
			side_nodes[1] != side_nodes[2])
		      node = 1;
		    else
		      node = 3;
		  } else {
		    /*
		     * in the non-degenerate case, if an element is connected
		     * to two opposite nodes, then it must share a face.
		     */
		    break;
		  }
		  
		} else {
		  if (!dflag) {
		    fprintf(stderr, "Possible corrupted mesh detected at element %lu, strange connectivity.\n", ecnt);
		  } 
		  /* This is the second or later time through this
		     loop and each time through, there have been two
		     or more elements that are connected to 'node'
		     (which changes) and 'inode'.  We want to make
		     sure that the elements matched this time through
		     were also in the list the first time through so
		     that the elements contain nodes 0 1 2 of the face
		     and not just 0 1 and 0 2...
		     So, this time, only put an element in the list if
		     it was in the list before.
		  */
		  for (size_t ncnt2 = 0; ncnt2 < nhold; ncnt2++)
		    hold_elem[ncnt2] = -hold_elem[ncnt2];
		  for (int ncnt3 = 0; ncnt3 < nelem; ncnt3++) {
		    for (size_t ncnt2 = 0; ncnt2 < nhold; ncnt2++) {
		      if (-hold_elem[ncnt2] == graph->sur_elem[side_nodes[inode]][pt_list[ncnt3]]) {
			hold_elem[ncnt2] = graph->sur_elem[side_nodes[inode]][pt_list[ncnt3]];
			break;
		      }
		    }
		  }
		  /* Now, go through list and cull out element < 0 */
		  size_t ncnt3 = 0;
		  for (size_t ncnt2 = 0; ncnt2 < nhold; ncnt2++) {
		    if (hold_elem[ncnt2] >= 0) {
		      hold_elem[ncnt3] = hold_elem[ncnt2];
		      ncnt3++;
		    }
		  }
		  nelem = ncnt3;
		}
	      }
	      else { /* nelem == 1 or 0 */
		if (!dflag) {
		  nhold = graph->sur_elem[side_nodes[1]].size();
		  for (size_t ncnt2 = 0; ncnt2 < nhold; ncnt2++)
		    hold_elem[ncnt2] = graph->sur_elem[side_nodes[1]][ncnt2];
		}
		inode = 1;
		node = 3; /* The node diagonally opposite node 1 */
	      }
	    }
	  }
	} /* "if (!hflag1)" */
	
	/*
	 * if there is an element on this side of ecnt, then there
	 * will be at least two elements in the intersection (one
	 * will be ecnt)
	 */
	if (nelem > 1) {

	  /*
	   * now go through and check each element in the list to see
	   * if it on a different processor than ecnt.  Don't need to
	   * worry about ecnt (which is in the list) since it is on
	   * the same processor as itself.  Note that due to filtering
	   * done above, we are guaranteed to either have an element
	   * on a different processor or elem==ecnt.
	   */
	  for (int ncnt = 0; ncnt < nelem; ncnt++) {

	    INT elem   = hold_elem[ncnt];
	    int proc2  = lb->vertex2proc[elem];
	    assert(proc2 < machine->num_procs);

	    if (proc != proc2) {

	      E_Type etype2 = mesh->elem_type[elem];

	      int dim2 = get_elem_info(NDIM, etype2);

	      int diff = abs(dim1 - dim2);  

	      /* 
	       * hex's to shells - ok
	       * shells to bar - ok
	       * hex to bar - BAD since a BAR will see a HEX but a HEX will not
	       *              see a BAR
	       */

	      if(diff < 2) { 

		/* need to check for hex's */
		hflag2 = is_hex(etype2);
		tflag2 = is_tet(etype2);

		/* check here for tet/hex combinations */
		if ((tflag1 && hflag2) || (hflag1 && tflag2)) {
		  /*
		   * have to call a special function to get the side id
		   * in these cases. In both cases, the number of side
		   * nodes for the element will not be consistent with
		   * side_cnt, and:
		   *
		   * TET/HEX - side_nodes only contains three of the
		   *           the side nodes of the hex.
		   *
		   * HEX/TET - Have to check that this tet shares a side
		   *           with the hex.
		   */
		  sid = get_side_id_hex_tet(mesh->elem_type[elem],
					    mesh->connect[elem],
					    side_cnt, side_nodes);
		} else {
		  /*
		   * get the side id of elem. Make sure that ecnt is
		   * trying to communicate to a valid side of elem
		   */
		  side_cnt = get_ss_mirror(etype, side_nodes, (nscnt+1),
					   mirror_nodes);

		  /*
		   * small kludge to handle 6 node faces butted up against
		   * 4 node faces
		   */
		  if (etype == HEXSHELL && side_cnt == 6) side_cnt = 4;

		  /*
		   * in order to get the correct side order for elem,
		   * get the mirror of the side of ecnt
		   */
		  sid = get_side_id(mesh->elem_type[elem], mesh->connect[elem],
				    side_cnt, mirror_nodes, problem->skip_checks,
				    problem->partial_adj);
		  
		}

		if (sid > 0) {
		  /* Element is a border element */
		  internal = false;
		  if(!flag) {
		    flag = 1;
		    lb->bor_elems[proc].push_back(ecnt);
		  }

		  /* now put ecnt into proc2's communications map */
		  lb->e_cmap_elems[proc2].push_back(elem);
		  lb->e_cmap_sides[proc2].push_back(sid);
		  lb->e_cmap_procs[proc2].push_back(proc);
		  lb->e_cmap_neigh[proc2].push_back(ecnt);

		} else if ((sid < 0) && (!problem->skip_checks)) {
		  /*
		   * too many errors with bad meshes, print out
		   * more information here for diagnostics
		   */
		  sprintf(cmesg,
			  "Error returned while getting side id for communication map.");
		  Gen_Error(0, cmesg);
		  sprintf(cmesg, "Element 1: %lu", (ecnt+1));
		  Gen_Error(0, cmesg);
		  nnodes = get_elem_info(NNODES, etype);
		  strcpy(cmesg, "connect table:");
		  for (int i = 0; i < nnodes; i++) {
		    sprintf(tmpstr, " %lu", (size_t)(mesh->connect[ecnt][i]+1));
		    strcat(cmesg, tmpstr);
		  }
		  Gen_Error(0, cmesg);
		  sprintf(cmesg, "side id: %lu", (size_t)(nscnt+1));
		  Gen_Error(0, cmesg);
		  strcpy(cmesg, "side nodes:");
		  for (int i = 0; i < side_cnt; i++) {
		    sprintf(tmpstr, " %lu", (size_t)(side_nodes[i]+1));
		    strcat(cmesg, tmpstr);
		  }
		  Gen_Error(0, cmesg);
		  sprintf(cmesg, "Element 2: %lu", (size_t)(elem+1));
		  Gen_Error(0, cmesg);
		  nnodes = get_elem_info(NNODES, etype2);
		  strcpy(cmesg, "connect table:");
		  for (int i = 0; i < nnodes; i++) {
		    sprintf(tmpstr, " %lu", (size_t)(mesh->connect[elem][i]+1));
		    strcat(cmesg, tmpstr);
		  }
		  Gen_Error(0, cmesg);
		  return 0; /* and get out of here */

		} /* End "if sid < 0 && !problem>skip_checks" */
	      } /* End "if (sid > 0)" */
	    } /* End "if (proc != proc2)" */
	  } /* End "for (ncnt = 0; ncnt < nelem; ncnt++)" */
	} /* End "if (nelem > 1)" */
      } /* End "for (nscnt = 0; nscnt < nsides; nscnt++)" */

      if(internal) {
	lb->int_elems[proc].push_back(ecnt);
      }
    } /* End "for(ecnt=0; ecnt < mesh->num_elems; ecnt++)" */

    /* free up memory */
    free (pt_list);
    free (hold_elem);

    time2 = get_time();
    printf("Time for elemental categorization: %fs\n", time2-time1);

    /* Find the internal and border nodes */

    time1 = get_time();
    for(size_t ncnt=0; ncnt < mesh->num_nodes; ncnt++) {
      bool internal = true;
      int  flag = 0;
      int proc = 0;

      /* If a node is not connected to any elements (graph->nsur_elem[ncnt] == 0),
	 then it will be assigned to processor 0 and treated as internal.
      */

      if(!graph->sur_elem[ncnt].empty()) {
	size_t elem = graph->sur_elem[ncnt][0];

	proc = lb->vertex2proc[elem];
	assert(proc < machine->num_procs);
	for(size_t ecnt=1; ecnt < graph->sur_elem[ncnt].size(); ecnt++) {
	  int proc2 = lb->vertex2proc[graph->sur_elem[ncnt][ecnt]];
	  assert(proc2 < machine->num_procs);
	  /* check if the processor for any two surrounding elems differ */
	  if (proc != proc2) {
	    /* ncnt is a border node  of proc */
	    internal = false;
	    /* first, I have to deal with node being border for proc */
	    if (!flag) {
	      flag = 1;                 /* only want to do this once */
	      lb->bor_nodes[proc].push_back(ncnt);
	    }

	    /*
	     * now I have to put ncnt in the border list for proc2
	     * I need to check to make sure that this node has not
	     * already been added to this list. If it has, then it
	     * is in the last position in the array
	     */
	    if ((lb->bor_nodes[proc2].empty()) ||
		((INT)ncnt != lb->bor_nodes[proc2][lb->bor_nodes[proc2].size()-1])) {
	      lb->bor_nodes[proc2].push_back(ncnt);
	    }
	  } /* if (proc != lb->vertex2proc[graph->sur_elem[ncnt][ecnt]]) */
	}   /* for(ecnt=1; ecnt < graph->nsur_elem[ncnt]; ecnt++) */
      } /* if(graph->nsur_elem[ncnt]) */

      if (internal) {
	/*
	 * NOTE: if all of the processors above were the same, then
	 * the one held in proc is the correct one
	 */
	lb->int_nodes[proc].push_back(ncnt);
      }
    }  /* for(ncnt=0; ncnt < machine->num_nodes; ncnt++) */

    time2 = get_time();
    printf("Nodal categorization: %fs\n", time2-time1);

    /* Allocate memory for the border node processor IDs */
    for(int proc=0; proc < machine->num_procs; proc++) {
      if (!lb->bor_nodes[proc].empty()) {
	lb->born_procs[proc].resize(lb->bor_nodes[proc].size());
      }
    }

    /* Now find the processor(s) associated with each border node */
    time1 = get_time();
    for(int pcnt=0; pcnt < machine->num_procs; pcnt++) {
      for(size_t ncnt=0; ncnt < lb->bor_nodes[pcnt].size(); ncnt++) {
	size_t node = lb->bor_nodes[pcnt][ncnt];

	for(size_t ecnt=0; ecnt < graph->sur_elem[node].size(); ecnt++) {
	  size_t elem = graph->sur_elem[node][ecnt];
	  INT proc = lb->vertex2proc[elem];
	  assert(proc < machine->num_procs);
	  if(proc != pcnt) {
	    if(in_list(proc, lb->born_procs[pcnt][ncnt]) < 0) {
	      lb->born_procs[pcnt][ncnt].push_back(proc);
	    }
	  } /* End "if(proc != pcnt)" */
	} /* End "for(ecnt=0; ecnt < graph->nsur_elems[node]; ecnt++)" */
      } /* End "for(ncnt=0; ncnt < lb->num_bor_nodes[pcnt]; ncnt++)" */
    } /* End "for(pcnt=0; pcnt < machine->num_procs; pcnt++)" */

    time2 = get_time();
    printf("Find procs for border nodes: %fs\n", time2-time1);

    /* Order the element communication maps by processor */
    time1 = get_time();
    for(int pcnt=0; pcnt < machine->num_procs; pcnt++) {
      /* Note that this sort is multi-key */
      qsort4(&lb->e_cmap_procs[pcnt][0],   /* 1st key */
	     &lb->e_cmap_elems[pcnt][0],   /* 2nd key */
	     &lb->e_cmap_neigh[pcnt][0],   /* 3rd key */
	     &lb->e_cmap_sides[pcnt][0],   /* 4th key */
	     lb->e_cmap_procs[pcnt].size());   /* Size */
    }
    /*
     * At this point, each processors arrays are sorted on three keys:
     * [processor, element, neighbor]
     */

    time2 = get_time();
    printf("Order elem cmaps: %fs\n", time2-time1);

    /*
     * Now order the elemental communication maps so that they are
     * consistent between processors.
     */
    time1 = get_time();
    for(INT pcnt=1; pcnt < machine->num_procs; pcnt++) {
      int save_fv1 = 0;
      ssize_t size = lb->e_cmap_procs[pcnt].size();       /* Define shortcuts size and procs */
      std::vector<INT> procs = lb->e_cmap_procs[pcnt];

      INT fv1 = -1;
      INT lv1 = -1;
      
      for(int pcnt2=0; pcnt2 < pcnt; pcnt2++) {

	/*
	 * Find the first and last entries for processor "pcnt2" in
	 * the list of processor "pcnt".
	 */

	/* lb->c_cmap_procs[pcnt] is sorted based on processor.
	 * From point 'save_fv1' search for 'pcnt2'
	 * If not found, value is -1; else search for !pcnt2 from
	 * that point forward.
	 */
	ssize_t i = save_fv1;
	while (i < size && procs[i] < pcnt2)
	  i++;
	if (i >= size || procs[i] != pcnt2) {
	  fv1 = -1;
	  lv1 = -1;
	}
	else {
	  fv1 = i;
	  assert(procs[i] == pcnt2);
	  for (lv1 = fv1; lv1 < size; lv1++) {
	    if (procs[lv1] != pcnt2) {
	      lv1 = lv1 - 1;
	      assert(procs[lv1] == pcnt2);
	      break;
	    }
	  }
	}
	    
	if (lv1 >= size)
	  lv1 = size-1;
	  
	if (lv1 != -1)
	  save_fv1 = lv1+1;

#if 0
	/* Old method -- can use for verification by uncommenting this if block  */
	{
	  int tst_fv1, tst_lv1;
	  find_first_last(pcnt2, size, procs, &tst_fv1, &tst_lv1);
	  assert(tst_fv1 == fv1);
	  assert(tst_lv1 == lv1);
	}
#endif
	  
	if(fv1 >= 0) {
	  /* Sort based on neighbor element */
	  sort3(lv1-fv1+1,
		(&lb->e_cmap_neigh[pcnt][fv1]),
		(&lb->e_cmap_elems[pcnt][fv1]),
		(&lb->e_cmap_sides[pcnt][fv1]));
	  /*
	   * Find the first and last entries for processor "pcnt" in
	   * the list of processor "pcnt2".
	   */
	  INT fv2 = -1;
	  INT lv2 = -1;
	  find_first_last(pcnt, lb->e_cmap_procs[pcnt2].size(), &lb->e_cmap_procs[pcnt2][0],
			  &fv2, &lv2);
#if 1
	  if (lv2-fv2 != lv1-fv1) {
	    fprintf(stderr, "%lu: %lu to %lu\n", (size_t)pcnt2, (size_t)fv1, (size_t)lv1);
	    for (i=fv1; i <= lv1; i++)
	      fprintf(stderr, "%lu: %lu\t%lu\t%lu\t%lu\n", (size_t)i,
		      (size_t)lb->e_cmap_elems[pcnt][i],
		      (size_t)lb->e_cmap_neigh[pcnt][i],
		      (size_t)lb->e_cmap_procs[pcnt][i],
		      (size_t)lb->e_cmap_sides[pcnt][i]);
	    fprintf(stderr, "%lu: %lu to %lu\n", (size_t)pcnt, (size_t)fv2, (size_t)lv2);
	    for (i=fv2; i <= lv2; i++)
	      fprintf(stderr, "%lu: %lu\t%lu\t%lu\t%lu\n", (size_t)i,
		      (size_t)lb->e_cmap_elems[pcnt2][i],
		      (size_t)lb->e_cmap_neigh[pcnt2][i],
		      (size_t)lb->e_cmap_procs[pcnt2][i],
		      (size_t)lb->e_cmap_sides[pcnt2][i]);
	  }
#endif
	  assert(lv2-fv2 == lv1-fv1);
		
	  /* Sort based on element -- This will then match order of
	   * the fv1->lv1 arrays.
	   */
	  sort3(lv2-fv2+1,
		(&lb->e_cmap_elems[pcnt2][fv2]),
		(&lb->e_cmap_neigh[pcnt2][fv2]),
		(&lb->e_cmap_sides[pcnt2][fv2]));
	  
	} /* End "if(fv1 >= 0)" */
      } /* End "for(pcnt2=0; pcnt2 < pcnt; pcnt2++)" */
    } /* End "for(pcnt=0; pcnt < machine->num_procs; pcnt++)" */

    time2 = get_time();
    printf("Make cmaps consistent: %fs\n", time2-time1);

    return 1;
  } /*--------------------------End elemental_dist()---------------------------*/


  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/
  /*
   * Determine a graph's connected components
   * Input:
   *   nrow = number of elements in the problem
   *   rows = the pointer into columns. 1 based and of dimension nrow+1
   *          The last value=1+dimension of columns.
   *   columns = the element to element connectivity.  Also 1 based
   *          and includes the self connectivity.
   * Output:
   *   return value = number of pieces
   *   list = all the elements grouped by connectivity, allocated
   *          before calling extract_connected_lists.
   *   list_ptr = pointers to the start of each connected component in
   *              list[]. The dimension of *list_ptr is the number of 
   *              connected components.  Memory is allocated within the
   *              routine. Memory must be freed using free(). List_ptr
   *              is zero based.
   * Example:
   * 
   *        1  2  3  4  5
   *     1  *     *  *
   *     2     *     
   *     3  *     *
   *     4  *        *
   *     5              *
   * There are three connected components: { 1 3 4 }, { 2 } and { 5 }
   * 
   * Input:
   *   nrow=5
   *   rows={   1     4 5   7   9 10}
   *   columns={1 3 4 2 1 3 1 4 5}
   * 
   * Output:
   *   returns 3
   *   list=    { 1 3 4 2 5 }
   *   list_ptr={ 0     3 4 }
   */

  int extract_connected_lists( int nrow, const int* columns,  
			       const int* rows, int* list, 
			       int** list_ptr )
  {
    int root, nordered, ni, nf, nni, nnf, i, ki, kf, k, j,
      components, too_many_components = 0;
    if (nrow == 0) return(0);

    /* Reasonable guess for number of components */
    int pieces = 10 +  (nrow/10);
    *list_ptr = (int *) malloc(pieces*sizeof(int));
    int *mask = (int *) malloc(nrow*sizeof(int));
    if (mask == NULL) {
      fprintf(stderr,"Memory exhausted in extract_connected_lists\n");
      return(0);
    }
    root = 1;
    for(i=0; i<nrow; i++)
      mask[i] = 1;
	
    components                  = 1;
    nordered                    = 1;
    (*list_ptr)[ components - 1 ] = nordered-1;
    list[ nordered-1 ]          = root;
    mask[root-1]                = 0;
    ni = 1;
    nf = 1;
    while( nordered < nrow ){
      if( nf == ni - 1 ){
	++components;
	if( components < pieces ){
	  (*list_ptr)[ components - 1 ] = nordered;
	}
	else{
	  too_many_components = 1;
	}
	for(i=0;i<nrow;i++){
	  if( mask[i] == 1 ){
	    ++nordered;
	    list[ nordered-1 ] = i+1;
	    mask[i]            = 0;
	    nf = ni;
	    break;
	  }
	}
      }
      nni = nf  + 1;
      nnf = nni - 1;
      for(i=ni; i <= nf; i++){
	ki = rows[ list[i-1] -1 ] - 1;
	kf = rows[ list[i-1]    ] - 2;
	for(k=ki; k <= kf;k++){
	  j = columns[k];
	  if( mask[j-1] == 1 ){
	    ++nnf;
	    ++nordered;
	    mask[j-1] = 0;
	    list[ nordered - 1 ] = j;
	  }
	}
      }
      ni = nni;
      nf = nnf;
    }
    if( too_many_components == 0 ){
      (*list_ptr)[ components ] = nordered;
      free((char *)mask);
      return(components);
    }

    /* Start over with *list_ptr correctly allocated */
    free((char *)list_ptr);
    *list_ptr = (int *) malloc(components*sizeof(int));
    if ( *list_ptr == NULL) {
      fprintf(stderr,"Memory exhausted in extract_connected_lists\n");
      return(0);
    }
    for(i=0; i<nrow; i++)
      mask[i] = 1;

    components                  = 1;
    nordered                    = 1;
    (*list_ptr)[ components - 1 ] = nordered-1;
    list[ nordered-1 ]          = root;
    mask[root-1]                = 0;
    ni = 1;
    nf = 1;
    while( nordered < nrow ){
      if( nf == ni - 1 ){
	++components;
	(*list_ptr)[ components - 1 ] = nordered;
	for(i=0;i<nrow;i++){
	  if( mask[i] == 1 ){
	    ++nordered;
	    list[ nordered-1 ] = i+1;
	    mask[i]            = 0;
	    nf = ni;
	    break;
	  }
	}
      }
      nni = nf  + 1;
      nnf = nni - 1;
      for(i=ni; i <= nf; i++){
	ki = rows[ list[i-1] -1 ] - 1;
	kf = rows[ list[i-1]    ] - 2;
	for(k=ki; k <= kf;k++){
	  j = columns[k];
	  if( mask[j-1] == 1 ){
	    ++nnf;
	    ++nordered;
	    mask[j-1] = 0;
	    list[ nordered - 1 ] = j;
	  }
	}
      }
      ni = nni;
      nf = nnf;
    }
    (*list_ptr)[ components ] = nordered;
    free((char *)mask);
    return(components);
  }

  /*****************************************************************************/
  int ZPINCH_assign(
		    Machine_Description* machine,  /* Machine MESH = nwedge * nslice  */
		    int ndot,             /* Length of x, y, z, and part (== # of elements) */
		    float *x,             /* x-coordinates */
		    float *y,             /* y-coordinates */
		    float *z,             /* z-coordinates */
		    int *part           /* Output:  partition assignments for each element */
		    )
  {
    /* Routine to generate a partition of a cylinder. 
     * Assumptions:
     *   Height of cylinder is aligned with z-axis.
     *   Face of cylinder is centered at (x,y) = (0,0).
     * The x,y-plane is partitioned into wedge-shaped slices 
     * (just like cutting a pie!).
     * The z-directon is partitioned into cylindrical slices
     * (just like cutting a hot dog!).
     *
     * Decomposition requested by Chris Garasi for Z-Pinch simulations,
     * NNSA Level-1 milestone for Petaflop Computer justification.
     * Decomposition implemented by Karen Devine, kddevin@sandia.gov.
     * March 18, 2004
     */

    int i;
    int nslice;                /* # of slices in the z-direction */
    int nwedge;                /* # of wedge-shaped slices in the x,y-plane */
    int slice;                 /* the slice to which an element is assigned */
    int wedge;                 /* the wedge to which an element is assigned */
    float zmin =  FLT_MAX;     /* Minimum and maximum z-coordinate values */
    float zmax = -FLT_MAX;
    double dz;                 /* zmax - zmin */
    double theta;              /* angle of (x,y) (polar coordinates) */
    double *wedge_max_theta;   /* max angle in each wedge. */
    double *slice_max_z;       /* max z value in each slice. */
    double epsilon = 1e-07;    /* tolerance that allows a point to be in wedge */

    if (ndot > 0 && (x == NULL || y == NULL || z == NULL || part == NULL)) {
      fprintf(stderr, "KDD -- Bad input to ZPINCH_assign.\n");
      fprintf(stderr, "KDD -- Contact Karen Devine, kddevin@sandia.gov.\n");
      exit(-1);
    }

    if (machine->type != MESH) {
      fprintf(stderr, "KDD -- Machine must be a MESH "
	      "with # wedges * # slices processors.\n");
      fprintf(stderr, "KDD -- Use nem_slice argument -m mesh=AxB, \n");
      fprintf(stderr, "KDD -- where A = # wedges in x,y-plane and \n");
      fprintf(stderr, "KDD -- B = # slices along z-axis.\n");
      exit(-1);
    }
    nwedge = machine->dim[0];
    nslice = machine->dim[1];
    printf("ZPINCH:  Computing\n"
	   "   %d slices in the z-direction\n"
	   "   %d wedges in the x,y-plane\n"
	   "   %d partitions total\n", nslice, nwedge, nslice * nwedge);

    /* Compute the maximum values of z */
    for (i = 0; i < ndot; i++) {
      if (z[i] > zmax) zmax = z[i];
      if (z[i] < zmin) zmin = z[i];
    }
    dz = zmax - zmin;

    /* Compute maximum z for each slice, using uniform partition of zmin - zmax */
    slice_max_z = (double*)malloc(nslice * sizeof(double));
    for (i = 0; i < nslice; i++)
      slice_max_z[i] = zmin + (double) (i+1) * dz / (double) nslice;

    /* Compute maximum angle for each wedge, using uniform partition of 2*M_PI */
    wedge_max_theta = (double*)malloc(nwedge * sizeof(double));
    for (i = 0; i < nwedge; i++)
      wedge_max_theta[i] = (double) (i+1) * (2. * M_PI) / (double) nwedge;


    /* Compute the partition assignment for each set of coordinates */
    for (i = 0; i < ndot; i++) {

      /* Compute the z slice that the element is in. */
      if (dz > 0.) {
	slice = (int) (nslice * (z[i] - zmin) / dz);
	if (slice == nslice) slice--;   /* Handles z[i] == zmax correctly */

	/* Move dots within epsilon of upper end of slice into next slice */
	/* This step reduces jagged edges due to roundoff in coordinate values */
	if (slice != nslice-1 && z[i] > (slice_max_z[slice] - epsilon)) slice++;
      }
      else  /* 2D problem */
	slice = 0;

      /* Compute polar coordinate theta in x,y-plane for the element. */
      if (x[i] == 0.) {
	if (y[i] >= 0.) theta = M_PI_2;
	else            theta = 3. * M_PI_2;
      }
      else {
	theta = atan(y[i] / x[i]);   /* In range -M_PI_2 to M_PI_2 */

	/* Convert to range 0 to 2*M_PI */
	if (x[i] < 0.)       theta += M_PI;
	else if (y[i] < 0.)  theta += 2 * M_PI;
      }

      /* Compute the wedge that the element is in. */
      wedge = (int) (nwedge * theta / (2 * M_PI));
      if (wedge == nwedge) wedge--;   /* Handles theta == 2*M_PI correctly */

      /* Move dots within epsilon of upper angle of wedge into next wedge */
      /* This step reduces jagged edges due to roundoff in coordinate values */
      if (theta > wedge_max_theta[wedge] - epsilon)
	wedge = (wedge + 1) % nwedge;

      /* Compute the part that the element is in. */
      part[i] = (int) (slice * nwedge + wedge);

    }

    free(wedge_max_theta);
    free(slice_max_z);

    return 0;
  }

  /*****************************************************************************/
  void BRICK_slices(
		    int nslices_d,      /* # of subdomains in this dimension */
		    int ndot,           /* # of dots */
		    float *d,           /* Array of ndot coordinates in this dimension */
		    float *dmin,        /* Output:  Smallest value in d[] */
		    float *dmax,        /* Output:  Largest value in d[] */
		    double *delta,      /* Output:  dmax - dmin */
		    double **slices_d   /* Output:  maximum d for each slice in dimension using
					   uniform partition of dmax - dmin */
		    )
  {
    /* Compute the min, max, delta and slices values for a single dimension. */
    /* KDDKDD Note:  This routine could also be used by ZPINCH in the z-direction,
     * KDDKDD Note:  but I don't want to tamper with code that already works. */
    int i;

    *dmin =  FLT_MAX;     /* Minimum and maximum coordinate values */
    *dmax = -FLT_MAX;

    /* Compute the minimum and maximum coordinate values */
    for (i = 0; i < ndot; i++) {
      if (d[i] > *dmax) *dmax = d[i];
      if (d[i] < *dmin) *dmin = d[i];
    }
    *delta = *dmax - *dmin;

    /* Compute maximum coordinate value for each slice, 
       using uniform partition of dmax - dmin */
    *slices_d = (double*)malloc(nslices_d * sizeof(double));
    for (i = 0; i < nslices_d; i++)
      (*slices_d)[i] = *dmin + (double) (i+1) * *delta / (double) nslices_d;
  }

  /*****************************************************************************/
  int BRICK_which_slice(
			int nslices_d,     /* # of subdomains in this dimension */
			float d,           /* Coordinate of dot in this dimension */
			float dmin,        /* Smallest value in d[] */
			double delta,      /* dmax - dmin */
			double *slices_d   /* Maximum d for each slice in dimension using
					      uniform partition of dmax - dmin */
			)
  {
    /* Function returning in which slice a coordinate d lies. */
    /* KDDKDD Note:  This routine could also be used by ZPINCH in the z-direction,
     * KDDKDD Note:  but I don't want to tamper with code that already works. */

    int d_slice;  /* Return value:  The slice to which the coordinate is assigned */
    double epsilon = 5e-06;    /* tolerance that allows a point to be in subdomain*/

    if (delta > 0.) {
      d_slice = (int) (nslices_d * (d - dmin) / delta);
      if (d_slice == nslices_d) d_slice--;   /* Handles d == dmax correctly */

      /* Move dots within epsilon of upper end of slice into next slice */
      /* This step reduces jagged edges due to roundoff in coordinate values */
      if (d_slice != nslices_d-1 && d > (slices_d[d_slice] - epsilon)) d_slice++;
    }
    else  /* not a 3D problem */
      d_slice = 0;

    return d_slice;
  }

  /*****************************************************************************/
  int BRICK_assign(
		   Machine_Description* machine,  /* Machine MESH = nx * ny * nz */
		   int ndot,             /* Length of x, y, z, and part (== # of elements) */
		   float *x,             /* x-coordinates */
		   float *y,             /* y-coordinates */
		   float *z,             /* z-coordinates */
		   int *part           /* Output:  partition assignments for each element */
		   )
  {
    /* Routine to generate a partition of an axis-aligned hexahedral domain.
     * Assumptions:
     *   Hexahedral domain is axis-aligned.
     *   For good balance, the number of elements in each dimension should
     *   be divisible by the number of processors used in that dimension.
     * Each dimension is partitioned into equally sized regions, generating
     * subdomains that are axis-aligned hexahedra.
     *
     * Decomposition requested by Tom Brunner and Chris Garasi for 
     * NNSA Level-1 milestone for Petaflop Computer justification.
     * Decomposition implemented by Karen Devine, kddevin@sandia.gov.
     * June 10, 2005
     */


    if (ndot > 0 && (x == NULL || y == NULL || z == NULL || part == NULL)) {
      fprintf(stderr, "KDD -- Bad input to BRICK_assign.\n");
      fprintf(stderr, "KDD -- Contact Karen Devine, kddevin@sandia.gov.\n");
      exit(-1);
    }

    if (machine->type != MESH) {
      fprintf(stderr, "KDD -- Machine must be a MESH "
	      "with nx * ny * nz processors.\n");
      fprintf(stderr, "KDD -- Use nem_slice argument -m mesh=AxBxC, \n");
      fprintf(stderr, "KDD -- where A = nx, B = ny, and "
	      "C = nz\n");
      exit(-1);
    }
    int nx = machine->dim[0];
    int ny = machine->dim[1];
    int nz = machine->dim[2];

    printf("BRICK:  Computing\n"
	   "   %d subdomains in the x direction\n"
	   "   %d subdomains in the y direction\n"
	   "   %d subdomains in the z direction\n"
	   "   %d partitions total\n", 
	   nx, ny, nz, nx * ny * nz);

    float xmin, xmax;          /* Minimum and maximum x-coordinate values */
    float ymin, ymax;          /* Minimum and maximum y-coordinate values */
    float zmin, zmax;          /* Minimum and maximum z-coordinate values */
    double dx;                 /* xmax - xmin */
    double dy;                 /* ymax - ymin */
    double dz;                 /* zmax - zmin */
    double *slices_x=NULL;     /* max x value in each slice. */
    double *slices_y=NULL;     /* max y value in each slice. */
    double *slices_z=NULL;     /* max z value in each slice. */

    /* Compute bounds of subdomains in each dimension */
    BRICK_slices(nx, ndot, x, &xmin, &xmax, &dx, &slices_x);
    BRICK_slices(ny, ndot, y, &ymin, &ymax, &dy, &slices_y);
    BRICK_slices(nz, ndot, z, &zmin, &zmax, &dz, &slices_z);

    /* Compute the partition assignment for each set of coordinates */
    for (int i = 0; i < ndot; i++) {

      /* Compute the slice in each dimension that the element is in */
      int x_slice = BRICK_which_slice(nx, x[i], xmin, dx, slices_x);
      int y_slice = BRICK_which_slice(ny, y[i], ymin, dy, slices_y);
      int z_slice = BRICK_which_slice(nz, z[i], zmin, dz, slices_z);

      /* Compute the part that the element is in. */
      part[i] = (int) (z_slice * (nx * ny) + y_slice * nx + x_slice);
    }

    free(slices_x);
    free(slices_y);
    free(slices_z);

    return 0;
  }

#ifdef USE_ZOLTAN
  /*****************************************************************************/
  /***** Global data structure used by Zoltan callbacks.                   *****/
  /***** Could implement Zoltan callbacks without global data structure,   *****/
  /***** but using the global data structure makes implementation quick.   *****/

  static struct {
    size_t ndot;             /* Length of x, y, z, and part (== # of elements) */
    int *vwgt;            /* vertex weights */
    float *x;             /* x-coordinates */
    float *y;             /* y-coordinates */
    float *z;             /* z-coordinates */
  } Zoltan_Data;

  /*****************************************************************************/
  /***** ZOLTAN CALLBACK FUNCTIONS *****/

  int zoltan_num_dim(void *data, int *ierr)
  {
    /* Return dimensionality of coordinate data.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
    *ierr = ZOLTAN_OK;
    if (Zoltan_Data.z != NULL) return 3;
    if (Zoltan_Data.y != NULL) return 2;
    return 1;
  }

  int zoltan_num_obj(void *data, int *ierr)
  {
    /* Return number of objects.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
    *ierr = ZOLTAN_OK;
    return Zoltan_Data.ndot;
  }

  void zoltan_obj_list(void *data, int ngid_ent, int nlid_ent,
		       ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
		       int wdim, float *wgts, int *ierr)
  {
    /* Return list of object IDs.  
     * Return only global IDs; don't need local IDs since running in serial.
     * gids are array indices for coordinate and vwgts arrays.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
    for (size_t i = 0; i < Zoltan_Data.ndot; i++) {
      gids[i] = i;
      if (wdim) wgts[i] = (float) Zoltan_Data.vwgt[i];
    }

    *ierr = ZOLTAN_OK;
    return;
  }

  void zoltan_geom(void *data, int ngid_ent, int nlid_ent, int nobj,
		   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
		   int ndim, double *geom, int *ierr)
  {
    /* Return coordinates for objects.
     * gids are array indices for coordinate arrays.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
  
    for (int i = 0; i < nobj; i++) {
      size_t j = gids[i];
      geom[i*ndim] = Zoltan_Data.x[j];
      if (ndim > 1) geom[i*ndim+1] = Zoltan_Data.y[j];
      if (ndim > 2) geom[i*ndim+2] = Zoltan_Data.z[j];
    }
     
    *ierr = ZOLTAN_OK;
    return;
  }
      
  /*****************************************************************************/
  int ZOLTAN_assign(
		    const char *method,         /* Zoltan LB_METHOD to use. */
		    int totalproc,        /* # of processors for which to partition */
		    size_t ndot,             /* Length of x, y, z, and part (== # of elements) */
		    int *vwgt,            /* vertex weights */
		    float *x,             /* x-coordinates */
		    float *y,             /* y-coordinates */
		    float *z,             /* z-coordinates */
		    int *part,          /* Output:  partition assignments for each element */
		    int argc,             /* Fields needed by MPI_Init */
		    char *argv[]          /* Fields needed by MPI_Init */
		    )
  {
    /* Function to allow Zoltan to compute decomposition using RCB.
     * Assuming running Zoltan in serial (as nem_slice is serial).
     * Return PARTITION_ASSIGNMENTS from Zoltan_LB_Partition; they should
     * match what is needed in part array above.
     */
    struct Zoltan_Struct *zz;
    int zngid_ent, znlid_ent;        /* Useful output from Zoltan_LB_Partition */
    int znobj;
    ZOLTAN_ID_PTR zgids, zlids;      /* Useful output from Zoltan_LB_Partition */
    int *zprocs, *zparts;            /* Useful output from Zoltan_LB_Partition */
    ZOLTAN_ID_PTR dummy1, dummy2;    /* Empty output from Zoltan_LB_Partition */
    int dummy0, *dummy3, *dummy4;    /* Empty output from Zoltan_LB_Partition */
    int ierr = ZOLTAN_OK;
    float ver;
    char str[10];
    int changes;

    /* Copy mesh data and pointers into structure accessible from callback fns. */
    Zoltan_Data.ndot = ndot;
    Zoltan_Data.vwgt = vwgt;
    Zoltan_Data.x = x;
    Zoltan_Data.y = y;
    Zoltan_Data.z = z;

    /* Initialize Zoltan */
    Zoltan_Initialize(argc, argv, &ver);
    zz = Zoltan_Create(MPI_COMM_WORLD);
    if (ierr) {
      fprintf(stderr, "Error returned from Zoltan_Create (%s:%d)\n",
	      __FILE__, __LINE__);
      goto End;
    }

    /* Register Callback functions */
    /* Using global Zoltan_Data; could register it here instead as data field. */
    Zoltan_Set_Fn(zz, ZOLTAN_NUM_GEOM_FN_TYPE,
		  (ZOLTAN_VOID_FN *) zoltan_num_dim, NULL);
    Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE,
		  (ZOLTAN_VOID_FN *) zoltan_num_obj, NULL);
    Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE,
		  (ZOLTAN_VOID_FN *) zoltan_obj_list, NULL);
    Zoltan_Set_Fn(zz, ZOLTAN_GEOM_MULTI_FN_TYPE,
		  (ZOLTAN_VOID_FN *) zoltan_geom, NULL);

    /* Set parameters for Zoltan */
    sprintf(str, "%d", totalproc);
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTITIONS", str);
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", method);
    Zoltan_Set_Param(zz, "REMAP", "0");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTITION_ASSIGNMENTS");
    if (vwgt) Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");

    /* Call partitioner */
    printf("Using Zoltan version %f, method %s\n", ver, method);
    ierr = Zoltan_LB_Partition(zz, &changes, &zngid_ent, &znlid_ent, 
			       &dummy0, &dummy1, &dummy2, &dummy3, &dummy4,
			       &znobj, &zgids, &zlids, &zprocs, &zparts);
    if (ierr) {
      fprintf(stderr, "Error returned from Zoltan_LB_Partition (%s:%d)\n",
	      __FILE__, __LINE__);
      goto End;
    }

    /* Sanity check */
    if (ndot != (size_t)znobj) {
      fprintf(stderr, "Sanity check failed; ndot %lu != znobj %lu.\n", 
	      (size_t)ndot, (size_t)znobj);
      goto End;
    }

    /* Convert data types from int to int. */
    for (size_t i = 0; i < ndot; i++)
      part[zgids[i]] = zparts[i];

  End:
    /* Clean up */
    Zoltan_LB_Free_Part(&zgids, &zlids, &zprocs, &zparts);
    Zoltan_Destroy(&zz);
    MPI_Finalize();
    if (ierr) exit(-1);
    return 0;
  }
#endif

  /*****************************************************************************/
  void BALANCE_STATS(
		     Machine_Description* machine,  /* Machine MESH = nwedge * nslice  */
		     int *wgt,            /* element weights; can be NULL if no weights  */
		     size_t ndot,             /* Length of x, y, z, and part (== # of elements) */
		     int *part           /* Partition assignments for each element */
		     )
  {
    /* Routine to print some info about the ZPINCH, BRICK or ZOLTAN_RCB 
       decompositions */
    int npart = machine->dim[0] * machine->dim[1] * machine->dim[2];
    int* cntwgt=NULL;

    int minwgt, maxwgt, sumwgt;
    
    int *cnts       = (int*)calloc(npart, sizeof(int));
    if (wgt) cntwgt = (int*)calloc(npart, sizeof(int));
    for (size_t i = 0; i < ndot; i++) {
      cnts[part[i]]++;
      if (wgt) cntwgt[part[i]] += wgt[i];
    }

    int max = 0;
    int min = ndot;
    size_t sum = 0;
    if (wgt) {
      maxwgt = 0;
      minwgt = INT_MAX;
      sumwgt = 0;
    }

    for (int i = 0; i < npart; i++) {
      if (cnts[i] > max) max = cnts[i];
      if (cnts[i] < min) min = cnts[i];
      sum += cnts[i];
      if (wgt) {
	if (cntwgt[i] > maxwgt) maxwgt = cntwgt[i];
	if (cntwgt[i] < minwgt) minwgt = cntwgt[i];
	sumwgt += cntwgt[i];
      }
      if (cnts[i] == 0) printf("ZERO on %d\n", i);
    }

    printf("CNT STATS:  min = %d  max = %d  avg = %f\n", min, max,
	   (float)sum/(float)npart);
    if (wgt)
      printf("WGT STATS:  min = %d  max = %d  avg = %f\n", minwgt, maxwgt,
	     (float)sumwgt/(float)npart);
    free(cnts);
    if (wgt) free(cntwgt);
  }
}
