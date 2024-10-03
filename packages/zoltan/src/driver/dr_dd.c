// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include "dr_const.h"
#include "dr_util_const.h"
#include "dr_err_const.h"
#include "dr_dd.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/****************************************************************************/
int build_elem_dd(MESH_INFO_PTR mesh) 
{
/* Create a distributed directory of the elements so we can track their
 * processor assignment after migrations.
 */
int maxelems;

  MPI_Allreduce(&(mesh->num_elems), &maxelems, 1, MPI_INT, MPI_MAX,
                zoltan_get_global_comm());
  if (Zoltan_DD_Create(&(mesh->dd), zoltan_get_global_comm(), 1, 0, 0, maxelems, 0) != 0){
    Gen_Error(0, "fatal:  NULL returned from Zoltan_DD_Create()\n");
    return 0;
  }

  return update_elem_dd(mesh);
}

/****************************************************************************/
int update_elem_dd(MESH_INFO_PTR mesh)
{
/*  Update a distributed directory of the elements with processor
 *  assignments initially and after migration.
 */
ELEM_INFO_PTR current_elem;
ZOLTAN_ID_PTR gids; 
int *parts;
int i, j;

  gids = (ZOLTAN_ID_PTR) malloc(mesh->num_elems * sizeof(ZOLTAN_ID_TYPE));
  parts = (int *) malloc(mesh->num_elems * sizeof(int));

  for (j = 0, i = 0; i < mesh->elem_array_len; i++) {
    current_elem = &(mesh->elements[i]);
    if (current_elem->globalID != ZOLTAN_ID_INVALID) {
      gids[j] = (ZOLTAN_ID_TYPE)current_elem->globalID;
      parts[j] = current_elem->my_part;
      j++;
    }
  }

  if (Zoltan_DD_Update(mesh->dd, gids, NULL,NULL, parts, mesh->num_elems)!=0) {
    Gen_Error(0, "fatal:  NULL returned from Zoltan_DD_Update()\n");
    return 0;
  }

  safe_free((void **)(void *) &gids);
  safe_free((void **)(void *) &parts);
  return 1;
}


/****************************************************************************/
int update_hvertex_proc(MESH_INFO_PTR mesh)
{
  int npins;

  npins = mesh->hindex[mesh->nhedges];  

  if (Zoltan_DD_Find(mesh->dd, mesh->hvertex, NULL, NULL, NULL, npins, mesh->hvertex_proc) != 0) {

    Gen_Error(0, "fatal:  NULL returned from Zoltan_DD_Find()\n");
    return 0;
  }

  return 1;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
