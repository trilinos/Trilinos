/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "dr_const.h"
#include "dr_elem_util_const.h"
#include "dr_util_const.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Utility functions for element initialization, etc.
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void initialize_element(ELEM_INFO *elem)
{
  int i; /* loop counter */

/*
 * Initializes all fields of an element.
 */
  elem->globalID = -1;
  elem->border = 0;
  elem->elem_blk = -1;
  elem->my_part = -1;
  elem->perm_value = -1;
  elem->invperm_value = -1;
  for (i=0; i<MAX_CPU_WGTS; i++)
    elem->cpu_wgt[i] = 0;
  elem->mem_wgt = 0;
  elem->nadj = 0;
  elem->adj_len = 0;
  elem->avg_coord[0] = elem->avg_coord[1] = elem->avg_coord[2] = 0.;
  elem->coord = NULL;
  elem->connect = NULL;
  elem->adj = NULL;
  elem->adj_proc = NULL;
  elem->edge_wgt = NULL;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void free_element_arrays(ELEM_INFO *elem, MESH_INFO_PTR mesh)
{
/*
 * Frees all memory malloc'ed for an individual element.
 */
int j;

  if (elem->coord != NULL) {
    for (j = 0; j < mesh->eb_nnodes[elem->elem_blk]; j++)
      safe_free((void **) &(elem->coord[j]));
    safe_free((void **) &(elem->coord));
  }
  safe_free((void **) &(elem->connect));
  safe_free((void **) &(elem->adj));
  safe_free((void **) &(elem->adj_proc));
  safe_free((void **) &(elem->edge_wgt));
  elem->avg_coord[0] = elem->avg_coord[1] = elem->avg_coord[2] = 0.;
  elem->globalID = -1;
  elem->border = 0;
  elem->my_part = -1;
  elem->perm_value = -1;
  elem->invperm_value = -1;
  elem->nadj = 0;
  elem->adj_len = 0;
  elem->elem_blk = -1;
  for (j=0; j<MAX_CPU_WGTS; j++)
    elem->cpu_wgt[j] = 0;
  elem->mem_wgt = 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void free_mesh_arrays(MESH_INFO_PTR mesh)
{
int i;

  for (i = 0; i < mesh->elem_array_len; i++) 
    free_element_arrays(&(mesh->elements[i]), mesh);
  safe_free((void **) &(mesh->elements));

  for (i = 0; i < mesh->num_el_blks; i++) 
    safe_free((void **) &(mesh->eb_names[i]));

  safe_free((void **) &(mesh->eb_names));
  safe_free((void **) &(mesh->eb_etypes));
  safe_free((void **) &(mesh->ecmap_id));
  safe_free((void **) &(mesh->ecmap_cnt));
  safe_free((void **) &(mesh->ecmap_elemids));
  safe_free((void **) &(mesh->ecmap_sideids));
  safe_free((void **) &(mesh->ecmap_neighids));
  safe_free((void **) &(mesh->hgid));
  safe_free((void **) &(mesh->hindex));
  safe_free((void **) &(mesh->hvertex));
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
