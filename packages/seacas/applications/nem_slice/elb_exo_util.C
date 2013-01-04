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

#include <exodusII.h>                   // for ex_close, ex_inquire, etc

#include <vector>
#include <string>
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for printf, NULL, sprintf
#include <stdlib.h>                     // for malloc, free, calloc
#include <string.h>                     // for strlen
#include <algorithm>
#include <assert.h>

#include "elb_exo.h"
#include "elb.h"                  // for Weight_Description<INT>, etc
#include "elb_elem.h"             // for get_elem_type, E_Type, etc
#include "elb_err.h"              // for Gen_Error, MAX_ERR_MSG
#include "elb_groups.h"           // for parse_groups
#include "elb_util.h"             // for in_list, roundfloat

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function read_exo_weights() begins:
 *----------------------------------------------------------------------------
 * This function reads the nodal or elemental values from an ExodusII file
 * which will be used by Chaco for weighting of the graph.
 *****************************************************************************/
template int read_exo_weights(Problem_Description* prob, Weight_Description<int>* weight);
template int read_exo_weights(Problem_Description* prob, Weight_Description<int64_t>* weight);

template <typename INT>
int read_exo_weights(Problem_Description* prob, Weight_Description<INT>* weight)
{
  int    exoid, cpu_ws=0, io_ws=0;
  int    neblks;
  float  version, minval = 1.0f;
  char elem_type[MAX_STR_LENGTH+1];
  char ctemp[1024];
/*---------------------------Execution Begins--------------------------------*/

  /* Open the ExodusII file containing the weights */
  int mode = EX_READ | prob->int64api;
  if((exoid=ex_open(weight->exo_filename.c_str(), mode, &cpu_ws, &io_ws,
                    &version)) < 0)
  {
    sprintf(ctemp, "fatal: could not open ExodusII file %s",
            weight->exo_filename.c_str());
    Gen_Error(0, ctemp);
    return 0;
  }

  std::vector<float> values(weight->nvals);
  if(prob->type == NODAL)
  {
    size_t tmp_nodes = ex_inquire_int(exoid, EX_INQ_NODES);
    /* check to make sure the sizes agree */
    if ((size_t)weight->nvals != tmp_nodes) {
      Gen_Error(0, "fatal: different number of nodes in mesh and weight files");
      ex_close(exoid);
      return 0;
    }

    weight->ow.resize(weight->nvals);
    /* Read in the nodal values */
    if(ex_get_nodal_var(exoid, weight->exo_tindx, weight->exo_vindx,
                        weight->nvals, TOPTR(values)) < 0)
    {
      Gen_Error(0, "fatal: unable to read nodal values");
      ex_close(exoid);
      return 0;
    }
  }
  else
  {
    size_t tmp_elem = ex_inquire_int(exoid, EX_INQ_ELEM);
    /* check to make sure the sizes agree */
    if ((size_t)weight->nvals != tmp_elem) {
      Gen_Error(0, "fatal: different number of elems in mesh and weight files");
      ex_close(exoid);
      return 0;
    }

    /* Get the number of element blocks */
    neblks = ex_inquire_int(exoid, EX_INQ_ELEM_BLK);
    std::vector<INT> eblk_ids(neblks);
    std::vector<INT> eblk_ecnts(neblks);

    if(ex_get_ids(exoid, EX_ELEM_BLOCK, &eblk_ids[0]) < 0)
    {
      Gen_Error(0, "fatal: unable to get element block IDs");
      ex_close(exoid);
      return 0;
    }

    /* Get the count of elements in each element block */
    for(int cnt=0; cnt < neblks; cnt++) {
      INT dum1, dum2;
      if(ex_get_elem_block(exoid, eblk_ids[cnt], elem_type,
                           &(eblk_ecnts[cnt]), &dum1, &dum2) < 0)
      {
        Gen_Error(0, "fatal: unable to get element block");
        ex_close(exoid);
        return 0;
      }
    }

    /* Get the element variables */
    size_t offset = 0;
    for(int cnt=0; cnt < neblks; cnt++)
    {
      if(ex_get_elem_var(exoid, weight->exo_tindx, weight->exo_vindx,
                         eblk_ids[cnt], eblk_ecnts[cnt], &(values[offset])) < 0)
      {
        Gen_Error(0, "fatal: unable to get element variable");
        ex_close(exoid);
        return 0;
      }
      offset += eblk_ecnts[cnt];
    }
  }

  /* Close the ExodusII weighting file */
  if(ex_close(exoid) < 0)
  {
    sprintf(ctemp, "warning: failed to close ExodusII file %s",
            weight->exo_filename.c_str());
    Gen_Error(0, ctemp);
  }

  /* now I need to translate the values to positive integers */

  /* first find the minimum value */
  minval = *std::min_element(values.begin(), values.end());

  /* now translate the values to be greater than 1 and convert to ints */
  for (int cnt=0; cnt < weight->nvals; cnt++) {
    values[cnt] += 1.0 - minval;
    weight->vertices[cnt] = roundfloat(values[cnt]);
  }
  return 1;
} /*------------------------End read_exo_weights()----------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function read_mesh_params() begins:
 *----------------------------------------------------------------------------
 * This function reads in information about the finite element mesh.
 *****************************************************************************/
template int read_mesh_params(const std::string &exo_file, Problem_Description* problem,
			      Mesh_Description<int>* mesh, Sphere_Info* sphere);
template int read_mesh_params(const std::string &exo_file, Problem_Description* problem,
			      Mesh_Description<int64_t>* mesh, Sphere_Info* sphere);


template <typename INT>
int read_mesh_params(const std::string &exo_file,
                     Problem_Description* problem,
                     Mesh_Description<INT>* mesh,
                     Sphere_Info* sphere
  )
{
  int    exoid, cpu_ws=0, io_ws=0;
  float  version;
  char   elem_type[MAX_STR_LENGTH+1];
/*---------------------------Execution Begins--------------------------------*/

  /* Open the ExodusII geometry file */
  int mode = EX_READ | problem->int64api;
  if((exoid=ex_open(exo_file.c_str(), mode, &cpu_ws, &io_ws, &version)) < 0)
  {
    Gen_Error(0, "fatal: unable to open ExodusII file for mesh params");
    return 0;
  }

  /* Get the init info */
  ex_init_params exo;
  if(ex_get_init_ext(exoid, &exo))
  {
    Gen_Error(0, "fatal: unable to get init info from ExodusII file");
    ex_close(exoid);
    return 0;
  }
  strcpy(mesh->title, exo.title);
  mesh->num_dims      = exo.num_dim;
  mesh->num_nodes     = exo.num_nodes;
  mesh->num_elems     = exo.num_elem;
  mesh->num_el_blks   = exo.num_elem_blk;
  mesh->num_node_sets = exo.num_node_sets;
  mesh->num_side_sets = exo.num_side_sets;

  /* Get the length of the concatenated node set node list */
  if(mesh->num_node_sets > 0)
  {
    mesh->ns_list_len = ex_inquire_int(exoid, EX_INQ_NS_NODE_LEN);
  }
  else
    mesh->ns_list_len = 0;

  /* Allocate and initialize memory for the sphere adjustment */
  sphere->adjust = (int*)malloc(sizeof(int)*3*(mesh->num_el_blks));
  if(!(sphere->adjust)) {
    Gen_Error(0, "fatal: insufficient memory");
    ex_close(exoid);
    return 0;
  }
  else {
    sphere->begin = sphere->adjust + mesh->num_el_blks;
    sphere->end   = sphere->begin  + mesh->num_el_blks;
    for(size_t cnt=0; cnt < mesh->num_el_blks; cnt++) {
      sphere->adjust[cnt] = 0;
      sphere->begin[cnt]  = 0;
      sphere->end[cnt]    = 0;
    }
  }

  std::vector<INT> el_blk_ids(mesh->num_el_blks);

  /* Read the element block IDs */
  if(ex_get_elem_blk_ids(exoid, &el_blk_ids[0]) < 0) {
    Gen_Error(0, "fatal: unable to get element block IDs");
    ex_close(exoid);
    return 0;
  }

  /* Determine the maximum number of nodes per element */
  mesh->max_np_elem = 0;
  for(size_t cnt=0; cnt < mesh->num_el_blks; cnt++) {
    INT num_elems, idum;
    INT nodes_in_elem;
    
    if(ex_get_elem_block(exoid, el_blk_ids[cnt], elem_type,
                         &num_elems, &nodes_in_elem, &idum) < 0)
    {
      Gen_Error(0, "fatal: unable to get element block");
      ex_close(exoid);
      return 0;
    }

    if(cnt == 0)
      sphere->end[0] = num_elems;

    if(get_elem_type(elem_type, nodes_in_elem, mesh->num_dims) == SPHERE && problem->no_sph != 1) {
      sphere->num  += num_elems;
      sphere->adjust[cnt] = 0;
    }
    else
      sphere->adjust[cnt] = sphere->num;

    if(cnt != 0) {
      sphere->begin[cnt] = sphere->end[cnt-1];
      sphere->end[cnt]   = sphere->begin[cnt] + num_elems;
    }

    mesh->max_np_elem = MAX(mesh->max_np_elem, (size_t)nodes_in_elem);
  }

  /* Close the ExodusII file */
  if(ex_close(exoid) < 0)
    Gen_Error(1, "warning: unable to close ExodusII file");

  printf("ExodusII mesh information\n");
  if(strlen(mesh->title) > 0)
    printf("\ttitle: %s\n", mesh->title);
  printf("\tgeometry dimension: %lu\n", mesh->num_dims);
  printf("\tnumber of nodes: %lu\tnumber of elements: %lu\n", mesh->num_nodes,
         mesh->num_elems);
  printf("\tnumber of element blocks: %lu\n", mesh->num_el_blks);
  printf("\tnumber of node sets: %lu\tnumber of side sets: %lu\n",
         mesh->num_node_sets, mesh->num_side_sets);

  return 1;

} /*--------------------------End read_mesh_params()-------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function read_mesh_params() begins:
 *----------------------------------------------------------------------------
 * This function reads in the finite element mesh.
 *****************************************************************************/
template int read_mesh(const std::string &exo_file, Problem_Description* problem,
		       Mesh_Description<int>* mesh, Weight_Description<int>* weight);
template int read_mesh(const std::string &exo_file, Problem_Description* problem,
		       Mesh_Description<int64_t>* mesh, Weight_Description<int64_t>* weight);


template <typename INT>
int read_mesh(const std::string &exo_file,
              Problem_Description* problem,
              Mesh_Description<INT>* mesh,
              Weight_Description<INT>* weight
	      )
{
  float  version, *xptr, *yptr, *zptr;
  char   elem_type[MAX_STR_LENGTH+1];
  E_Type blk_elem_type;
  
  /*---------------------------Execution Begins--------------------------------*/

  /* Open the ExodusII file */
  int exoid, cpu_ws=0, io_ws=0;
  int mode = EX_READ | problem->int64api;
  if((exoid=ex_open(exo_file.c_str(), mode, &cpu_ws, &io_ws, &version)) < 0)
    {
      Gen_Error(0, "fatal: unable to open ExodusII mesh file");
      return 0;
    }

  /* Read the coordinates, if desired */
  xptr = yptr = zptr = NULL;
  if(problem->read_coords == ELB_TRUE)
    {
      switch(mesh->num_dims)
	{
	case 3:
	  zptr = (mesh->coords)+2*(mesh->num_nodes);
	  /* FALLTHRU */
	case 2:
	  yptr = (mesh->coords)+(mesh->num_nodes);
	  /* FALLTHRU */
	case 1:
	  xptr = mesh->coords;
	}

      if(ex_get_coord(exoid, xptr, yptr, zptr) < 0)
	{
	  Gen_Error(0, "fatal: unable to read coordinate values for mesh");
	  return 0;
	}

    } /* End "if(problem->read_coords == ELB_TRUE)" */

  /* Read the element block IDs */
  std::vector<INT> el_blk_ids(mesh->num_el_blks);
  std::vector<INT> el_blk_cnts(mesh->num_el_blks);

  if(ex_get_elem_blk_ids(exoid, &el_blk_ids[0]) < 0)
    {
      Gen_Error(0, "fatal: unable to read element block IDs");
      return 0;
    }

  /* Read the element connectivity */
  size_t gelem_cnt=0;
  for(size_t cnt=0; cnt < mesh->num_el_blks; cnt++) {
    INT nodes_per_elem, num_attr;
    if(ex_get_elem_block(exoid, el_blk_ids[cnt], elem_type,
                         &(el_blk_cnts[cnt]), &nodes_per_elem,
                         &num_attr) < 0)
      {
	Gen_Error(0, "fatal: unable to read element block");
	return 0;
      }

    blk_elem_type = get_elem_type(elem_type, nodes_per_elem, mesh->num_dims);

    INT *blk_connect = (INT*)malloc(sizeof(INT)*el_blk_cnts[cnt]*nodes_per_elem);
    if(!blk_connect)
      {
	Gen_Error(0, "fatal: insufficient memory");
	return 0;
      }

    /* Get the connectivity for this element block */
    if(ex_get_elem_conn(exoid, el_blk_ids[cnt], blk_connect) < 0)
      {
	Gen_Error(0, "fatal: failed to get element connectivity");
	return 0;
      }

    /* find out if this element block is weighted */
    int wgt = -1;
    if (weight->type & EL_BLK)
      wgt = in_list(el_blk_ids[cnt], weight->elemblk);
    
    /* Fill the 2D global connectivity array */
    if (((problem->type == ELEMENTAL) && (weight->type & EL_BLK)) ||
        ((problem->type == NODAL) && (weight->type & EL_BLK))) {
      
      for(size_t cnt2=0; cnt2 < el_blk_cnts[cnt]; cnt2++) {
	mesh->elem_type[gelem_cnt] = blk_elem_type;
      
	/* while going through the blocks, take care of the weighting */
	if ((problem->type == ELEMENTAL) && (weight->type & EL_BLK)) {
	  /* is this block weighted */
	  if (wgt >= 0) {
	    /* check if there is a read value */
	    if (weight->vertices[gelem_cnt] >= 1) {
	      /* and if it should be overwritten */
	      if (weight->ow_read)
		weight->vertices[gelem_cnt] = weight->elemblk_wgt[wgt];
	    }
	    else
	      weight->vertices[gelem_cnt] = weight->elemblk_wgt[wgt];
	  }
	  else {
	    /* now check if this weight has been initialized */
	    if (weight->vertices[gelem_cnt] < 1)
	      weight->vertices[gelem_cnt] = 1;
	  }
	}

	for(size_t cnt3=0; cnt3 < nodes_per_elem; cnt3++) {
	  INT node = blk_connect[cnt3 + cnt2*nodes_per_elem] - 1;
	  assert(node >= 0);
	  mesh->connect[gelem_cnt][cnt3] = node;

	  /* deal with the weighting if necessary */
	  if ((problem->type == NODAL) && (weight->type & EL_BLK)) {
	    /* is this block weighted */
	    if (wgt >= 0) {
	      /* check if I read an exodus file */
	      if (weight->type & READ_EXO) {
		/* check if it can be overwritten */
		if (weight->ow_read) {
		  /* check if it has been overwritten already */
		  if (weight->ow[node]) {
		    weight->vertices[node] =
		      MAX(weight->vertices[node], weight->elemblk_wgt[wgt]);
		  }
		  else {
		    weight->vertices[node] = weight->elemblk_wgt[wgt];
		    weight->ow[node] = 1;   /* read value has been overwritten */
		  }
		}
	      }
	      else {
		weight->vertices[node] =
		  MAX(weight->vertices[node], weight->elemblk_wgt[wgt]);
	      }
	    }
	    else {
	      /* now check if this weight has been initialized */
	      if (weight->vertices[node] < 1)
		weight->vertices[node] = 1;
	    }
	  }
	}
	gelem_cnt++;
      }
    } else {
      // No weights...
      for(size_t cnt2=0; cnt2 < el_blk_cnts[cnt]; cnt2++) {
	mesh->elem_type[gelem_cnt] = blk_elem_type;

	for(size_t cnt3=0; cnt3 < nodes_per_elem; cnt3++) {
	  INT node = blk_connect[cnt2*nodes_per_elem + cnt3] - 1;
	  assert(node >= 0);
	  mesh->connect[gelem_cnt][cnt3] = node;
	}

	gelem_cnt++;
      }
    }
    /* Free up memory */
    free(blk_connect);

  } /* End "for(cnt=0; cnt < mesh->num_el_blks; cnt++)" */

  /* if there is a group designator, then parse it here */
  if (problem->groups != NULL) {
    if (!parse_groups(&el_blk_ids[0], &el_blk_cnts[0], mesh, problem)) {
      Gen_Error(0, "fatal: unable to parse group designator");
      ex_close(exoid);
      return 0;
    }
  }
  else problem->num_groups = 1; /* there is always one group */

  /* Close the ExodusII file */
  if(ex_close(exoid) < 0)
    Gen_Error(0, "warning: failed to close ExodusII mesh file");

  return 1;

} /*---------------------------End read_mesh()-------------------------------*/


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function init_weight_struct() begins:
 *----------------------------------------------------------------------------
 * This function initializes the weight structure given the current mesh.
 *****************************************************************************/
template int init_weight_struct(Problem_Description* problem, Mesh_Description<int>* mesh,
				Weight_Description<int>* weight);
template int init_weight_struct(Problem_Description* problem, Mesh_Description<int64_t>* mesh,
				Weight_Description<int64_t>* weight);


template <typename INT>
int init_weight_struct(Problem_Description* problem,
                       Mesh_Description<INT>* mesh,
                       Weight_Description<INT>* weight)
{
  if(problem->type == NODAL) weight->nvals = mesh->num_nodes;
  else                       weight->nvals = mesh->num_elems;

  /* Allocate memory */
  weight->vertices.resize(weight->nvals);
  return 1;
} /*-----------------------End init_weight_struct()--------------------------*/
