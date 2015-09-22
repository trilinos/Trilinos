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

#include "zoltan_dd_cpp.h"

static Zoltan_DD *dd = NULL;

/****************************************************************************/
int build_elem_dd(MESH_INFO_PTR mesh) 
{
  destroy_elem_dd();

  dd = new Zoltan_DD(MPI_COMM_WORLD, 1, 0, 0, 0, 0);

  return update_elem_dd(mesh);
}
/****************************************************************************/
void destroy_elem_dd()
{
  if (dd)
    {
    delete dd;
    dd = NULL;
    }
}
/****************************************************************************/
int update_elem_dd(MESH_INFO_PTR mesh)
{
/*  Update a distributed directory of the elements with processor
 *  assignments initially and after migration.
 */

  ZOLTAN_ID_PTR gids = new ZOLTAN_ID_TYPE [mesh->num_elems];
  int *parts = new int [mesh->num_elems];

  int i, j;

  for (j = 0, i = 0; i < mesh->elem_array_len; i++) {
    ELEM_INFO_PTR current_elem = &(mesh->elements[i]);
    if (current_elem->globalID != ZOLTAN_ID_INVALID) {
      gids[j] = current_elem->globalID;
      parts[j] = current_elem->my_part;
      j++;
    }
  }

  int rc = 1;

  if (dd->Update(gids, NULL,NULL, parts, mesh->num_elems)!=0) {
    Gen_Error(0, "fatal:  NULL returned from Zoltan_DD::Update()\n");
    rc = 0;
  }

  delete [] gids;
  delete [] parts;
  
  return rc;
}

/****************************************************************************/
int update_hvertex_proc(MESH_INFO_PTR mesh)
{
  if (dd->Find(mesh->hvertex, NULL, NULL, NULL, 
                     mesh->hindex[mesh->nhedges], mesh->hvertex_proc) != 0) {
    Gen_Error(0, "fatal:  NULL returned from Zoltan_DD::Find()\n");
    return 0;
  }

  return 1;
}
