// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include "dr_const.h"
#include "dr_elem_util_const.h"
#include "dr_util_const.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Utility functions for element initialization, etc.
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void initialize_element(ELEM_INFO *elem)
{
  int i; /* loop counter */

/*
 * Initializes all fields of an element.
 */
  elem->globalID = ZOLTAN_ID_INVALID;
  elem->border = 0;
  elem->elem_blk = 0;
  elem->my_part = -1;
  elem->fixed_part = -1;
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
  elem->adj_blank = NULL;
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
      safe_free((void **)(void *) &(elem->coord[j]));
    safe_free((void **)(void *) &(elem->coord));
  }
  safe_free((void **)(void *) &(elem->connect));
  safe_free((void **)(void *) &(elem->adj));
  safe_free((void **)(void *) &(elem->adj_proc));
  safe_free((void **)(void *) &(elem->adj_blank));
  safe_free((void **)(void *) &(elem->edge_wgt));
  elem->avg_coord[0] = elem->avg_coord[1] = elem->avg_coord[2] = 0.;
  elem->globalID = ZOLTAN_ID_INVALID;
  elem->border = 0;
  elem->my_part = -1;
  elem->fixed_part = -1;
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
  safe_free((void **)(void *) &(mesh->elements));
  safe_free((void **)(void *) &(mesh->blank));

  for (i = 0; i < mesh->num_el_blks; i++) 
    safe_free((void **)(void *) &(mesh->eb_names[i]));

  safe_free((void **)(void *) &(mesh->eb_names));
  safe_free((void **)(void *) &(mesh->eb_etypes));
  safe_free((void **)(void *) &(mesh->eb_cnts));
  safe_free((void **)(void *) &(mesh->ecmap_id));
  safe_free((void **)(void *) &(mesh->ecmap_cnt));
  safe_free((void **)(void *) &(mesh->ecmap_elemids));
  safe_free((void **)(void *) &(mesh->ecmap_sideids));
  safe_free((void **)(void *) &(mesh->ecmap_neighids));
  safe_free((void **)(void *) &(mesh->hgid));
  safe_free((void **)(void *) &(mesh->hindex));
  safe_free((void **)(void *) &(mesh->hvertex));
  safe_free((void **)(void *) &(mesh->hvertex_proc));
  safe_free((void **)(void *) &(mesh->heWgtId));
  safe_free((void **)(void *) &(mesh->hewgts));
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
