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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <exodusII.h>

#include "elb_const.h"
#include "elb_exo_const.h"
#include "elb_err_const.h"
#include "elb_elem_const.h"
#include "elb_util_const.h"
#include "elb_groups_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function read_exo_weights() begins:
 *----------------------------------------------------------------------------
 * This function reads the nodal or elemental values from an ExodusII file
 * which will be used by Chaco for weighting of the graph.
 *****************************************************************************/
int read_exo_weights(PROB_INFO_PTR prob,
                     WEIGHT_INFO_PTR weight
  )
{
  int    exoid, cpu_ws=0, io_ws=0;
  size_t cnt, offset;
  int    dum1, dum2;
  int    neblks, *eblk_ids, *eblk_ecnts;
  float  version, *values, minval = 1.0f;
  char   ctemp[MAX_ERR_MSG+1], elem_type[MAX_STR_LENGTH+1];
/*---------------------------Execution Begins--------------------------------*/

  /* Open the ExodusII file containing the weights */
  if((exoid=ex_open(weight->exo_filename, EX_READ, &cpu_ws, &io_ws,
                    &version)) < 0)
  {
    sprintf(ctemp, "fatal: could not open ExodusII file %s",
            weight->exo_filename);
    Gen_Error(0, ctemp);
    return 0;
  }

  if(prob->type == NODAL)
  {
    if(ex_inquire(exoid, EX_INQ_NODES, &dum1, NULL, NULL) < 0)
    {
      Gen_Error(0, "fatal: unable to read number of nodes");
      ex_close(exoid);
      return 0;
    }

    /* check to make sure the sizes agree */
    if (weight->nvals != dum1)
    {
      Gen_Error(0, "fatal: different number of nodes in mesh and weight files");
      ex_close(exoid);
      return 0;
    }

    /* Allocate memory */
    values = malloc(sizeof(float)*(weight->nvals));
    /* also allocate the overwritten array, and set all values to 0 */
    weight->ow = calloc(weight->nvals, sizeof(float));

    if(!values)
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exoid);
      return 0;
    }

    /* Read in the nodal values */
    if(ex_get_nodal_var(exoid, weight->exo_tindx, weight->exo_vindx,
                        weight->nvals, values) < 0)
    {
      Gen_Error(0, "fatal: unable to read nodal values");
      ex_close(exoid);
      return 0;
    }
  }
  else
  {
    if(ex_inquire(exoid, EX_INQ_ELEM, &dum1, NULL, NULL) < 0)
    {
      Gen_Error(0, "fatal: unable to read number of elements");
      ex_close(exoid);
      return 0;
    }

    /* check to make sure the sizes agree */
    if (weight->nvals != dum1)
    {
      Gen_Error(0, "fatal: different number of elems in mesh and weight files");
      ex_close(exoid);
      return 0;
    }

    /* Allocate memory */
    values = malloc(sizeof(float)*(weight->nvals));
    if(!values)
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exoid);
      return 0;
    }

    /* Get the number of element blocks */
    if(ex_inquire(exoid, EX_INQ_ELEM_BLK, &neblks, NULL, NULL) < 0)
    {
      sprintf(ctemp, "fatal: unable to get element block count in %s",
              weight->exo_filename);
      Gen_Error(0, ctemp);
      ex_close(exoid);
      return 0;
    }

    /* Allocate memory for element block IDs */
    eblk_ids = malloc(sizeof(int)*neblks);
    eblk_ecnts = malloc(sizeof(int)*neblks);
    if(!eblk_ids || !eblk_ecnts)
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exoid);
      return 0;
    }

    if(ex_get_elem_blk_ids(exoid, eblk_ids) < 0)
    {
      Gen_Error(0, "fatal: unable to get element block IDs");
      ex_close(exoid);
      return 0;
    }

    /* Get the count of elements in each element block */
    for(cnt=0; cnt < neblks; cnt++)
    {
      if(ex_get_elem_block(exoid, eblk_ids[cnt], elem_type,
                           &(eblk_ecnts[cnt]), &dum1, &dum2) < 0)
      {
        Gen_Error(0, "fatal: unable to get element block");
        ex_close(exoid);
        return 0;
      }
    }

    /* Get the element variables */
    offset = 0;
    for(cnt=0; cnt < neblks; cnt++)
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

    /* free scratch space */
    free (eblk_ids);
    free (eblk_ecnts);

  }

  /* Close the ExodusII weighting file */
  if(ex_close(exoid) < 0)
  {
    sprintf(ctemp, "warning: failed to close ExodusII file %s",
            weight->exo_filename);
    Gen_Error(0, ctemp);
  }

  /* now I need to translate the values to positive integers */

  /* first find the minimum value */
  for (cnt=0; cnt < weight->nvals; cnt++)
    if (values[cnt] < minval) minval = values[cnt];

  /* now translate the values to be greater than 1 and convert to ints */
  for (cnt=0; cnt < weight->nvals; cnt++) {
    values[cnt] += 1.0 - minval;
    weight->vertices[cnt] = roundfloat(values[cnt]);
  }

  free (values);

  return 1;

} /*------------------------End read_exo_weights()----------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function read_mesh_params() begins:
 *----------------------------------------------------------------------------
 * This function reads in information about the finite element mesh.
 *****************************************************************************/
int read_mesh_params(char exo_file[],
                     PROB_INFO_PTR problem,
                     MESH_INFO_PTR mesh,
                     SPHERE_INFO_PTR sphere
  )
{
  int    exoid, idum, num_elems, cnt, nodes_in_elem, cpu_ws=0, io_ws=0;
  int   *el_blk_ids;
  float  version;
  char   elem_type[MAX_STR_LENGTH+1];
/*---------------------------Execution Begins--------------------------------*/

  /* Open the ExodusII geometry file */
  if((exoid=ex_open(exo_file, EX_READ, &cpu_ws, &io_ws, &version)) < 0)
  {
    Gen_Error(0, "fatal: unable to open ExodusII file for mesh params");
    return 0;
  }

  /* Get the init info */
  if(ex_get_init(exoid, mesh->title, &(mesh->num_dims), &(mesh->num_nodes),
                 &(mesh->num_elems), &(mesh->num_el_blks),
                 &(mesh->num_node_sets), &(mesh->num_side_sets)) < 0)
  {
    Gen_Error(0, "fatal: unable to get init info from ExodusII file");
    ex_close(exoid);
    return 0;
  }

  /* Get the length of the concatenated node set node list */
  if(mesh->num_node_sets > 0)
  {
    if(ex_inquire(exoid, EX_INQ_NS_NODE_LEN, &(mesh->ns_list_len),
                  NULL, NULL) < 0)
    {
      Gen_Error(0, "fatal: unable to obtain node set node list length");
      ex_close(exoid);
      return 0;
    }
  }
  else
    mesh->ns_list_len = 0;

  /* Allocate and initialize memory for the sphere adjustment */
  sphere->adjust = malloc(sizeof(int)*3*(mesh->num_el_blks));
  if(!(sphere->adjust)) {
    Gen_Error(0, "fatal: insufficient memory");
    ex_close(exoid);
    return 0;
  }
  else {
    sphere->begin = sphere->adjust + mesh->num_el_blks;
    sphere->end   = sphere->begin  + mesh->num_el_blks;
    for(cnt=0; cnt < mesh->num_el_blks; cnt++) {
      sphere->adjust[cnt] = 0;
      sphere->begin[cnt]  = 0;
      sphere->end[cnt]    = 0;
    }
  }

  /* Allocate temporary memory for element block IDs */
  el_blk_ids = malloc(sizeof(int)*(mesh->num_el_blks));
  if(!el_blk_ids)
  {
    Gen_Error(0, "fatal: insufficient memory");
    ex_close(exoid);
    return 0;
  }

  /* Read the element block IDs */
  if(ex_get_elem_blk_ids(exoid, el_blk_ids) < 0)
  {
    Gen_Error(0, "fatal: unable to get element block IDs");
    ex_close(exoid);
    return 0;
  }

  /* Determine the maximum number of nodes per element */
  mesh->max_np_elem = 0;
  for(cnt=0; cnt < mesh->num_el_blks; cnt++)
  {
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

    mesh->max_np_elem = MAX(mesh->max_np_elem, nodes_in_elem);
  }

  /* Free scratch memory */
  free(el_blk_ids);

  /* Close the ExodusII file */
  if(ex_close(exoid) < 0)
    Gen_Error(1, "warning: unable to close ExodusII file");

  printf("ExodusII mesh information\n");
  if(strlen(mesh->title) > 0)
    printf("\ttitle: %s\n", mesh->title);
  printf("\tgeometry dimension: %d\n", mesh->num_dims);
  printf("\tnumber of nodes: %d\tnumber of elements: %d\n", mesh->num_nodes,
         mesh->num_elems);
  printf("\tnumber of element blocks: %d\n", mesh->num_el_blks);
  printf("\tnumber of node sets: %d\tnumber of side sets: %d\n",
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
int read_mesh(char exo_file[],
              PROB_INFO_PTR problem,
              MESH_INFO_PTR mesh,
              WEIGHT_INFO_PTR weight
  )
{
  int    exoid, cpu_ws=0, io_ws=0, gelem_cnt=0;
  size_t cnt, cnt2, cnt3;
  int    nodes_per_elem, num_attr;
  int   *el_blk_ids, *el_blk_cnts, *blk_connect;
  int    wgt;
  float  version, *xptr, *yptr, *zptr;
  char   elem_type[MAX_STR_LENGTH+1];
  E_Type blk_elem_type;
  
   /*---------------------------Execution Begins--------------------------------*/

  /* Open the ExodusII file */
  if((exoid=ex_open(exo_file, EX_READ, &cpu_ws, &io_ws, &version)) < 0)
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
  el_blk_ids = malloc(sizeof(int)*2*(mesh->num_el_blks));
  if(!el_blk_ids)
  {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  el_blk_cnts = el_blk_ids + mesh->num_el_blks;

  if(ex_get_elem_blk_ids(exoid, el_blk_ids) < 0)
  {
    Gen_Error(0, "fatal: unable to read element block IDs");
    return 0;
  }

  /* Read the element connectivity */
  for(cnt=0; cnt < mesh->num_el_blks; cnt++)
  {
    if(ex_get_elem_block(exoid, el_blk_ids[cnt], elem_type,
                         &(el_blk_cnts[cnt]), &nodes_per_elem,
                         &num_attr) < 0)
    {
      Gen_Error(0, "fatal: unable to read element block");
      return 0;
    }

    blk_elem_type = get_elem_type(elem_type, nodes_per_elem, mesh->num_dims);

    blk_connect = malloc(sizeof(int)*el_blk_cnts[cnt]*nodes_per_elem);
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
    wgt = -1;
    if (weight->type & EL_BLK)
      wgt = in_list(el_blk_ids[cnt], weight->num_ebw, weight->elemblk);

    /* Fill the 2D global connectivity array */
    for(cnt2=0; cnt2 < el_blk_cnts[cnt]; cnt2++)
    {
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

      for(cnt3=0; cnt3 < nodes_per_elem; cnt3++)
      {
        int node;
        node = blk_connect[cnt3 + cnt2*nodes_per_elem] - 1;

	if (node < 0) {
	  char ctemp[MAX_ERR_MSG+1];
	  sprintf(ctemp, "fatal: Invalid connectivity (%d) for block %d, element %d, local node %d.",
		  node, el_blk_ids[cnt], (int)cnt2+1, (int)cnt3+1);
	  Gen_Error(0, ctemp);
	  return 0;
	}
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
        mesh->connect[gelem_cnt][cnt3] = node;
      }

      gelem_cnt++;
    }

    /* Free up memory */
    free(blk_connect);

  } /* End "for(cnt=0; cnt < mesh->num_el_blks; cnt++)" */

  /* if there is a group designator, then parse it here */
  if (problem->groups != NULL) {
    if (!parse_groups(el_blk_ids, el_blk_cnts, mesh, problem)) {
      Gen_Error(0, "fatal: unable to parse group designator");
      ex_close(exoid);
      return 0;
    }
  }
  else problem->num_groups = 1; /* there is always one group */

  /* Close the ExodusII file */
  if(ex_close(exoid) < 0)
    Gen_Error(0, "warning: failed to close ExodusII mesh file");

  /* Free up unused variables */
  free(el_blk_ids);

  return 1;

} /*---------------------------End read_mesh()-------------------------------*/


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function init_weight_struct() begins:
 *----------------------------------------------------------------------------
 * This function initializes the weight structure given the current mesh.
 *****************************************************************************/
int init_weight_struct(PROB_INFO_PTR problem,
                       MESH_INFO_PTR mesh,
                       WEIGHT_INFO_PTR weight
  )
{
/*---------------------------Execution Begins--------------------------------*/

  if(problem->type == NODAL) weight->nvals = mesh->num_nodes;
  else                       weight->nvals = mesh->num_elems;

  /* Allocate memory */
  weight->vertices = calloc((weight->nvals), sizeof(int));
  if(!(weight->vertices))
  {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  return 1;

} /*-----------------------End init_weight_struct()--------------------------*/
