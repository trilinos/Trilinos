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
                MPI_COMM_WORLD);
  if (Zoltan_DD_Create(&(mesh->dd), MPI_COMM_WORLD, 1, 0, 0, maxelems, 0) != 0){
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
