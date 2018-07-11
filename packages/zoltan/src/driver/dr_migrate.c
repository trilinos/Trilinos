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

/*--------------------------------------------------------------------------*/
/* Purpose: Call Zoltan to migrate elements.                                */
/*          Contains all of the callback functions that Zoltan needs        */
/*          for the migration.                                              */
/*                                                                          */
/* General migration strategy:                                              */
/*  1. In migrate_pre_process, reset all adjacency info for local elems,    */
/*     using the local ID for adj. elems that are or will be (after         */
/*     migration) on this processor, and the global ID for adj elems that   */
/*     are not or will not be (after migration)  on this processor.         */
/*  2. When exporting elements, convert all the export elems' adjacencies'  */
/*     local IDs to global IDs.                                             */ 
/*  3. When importing elements, convert import elems' adjacencies that are  */
/*     local elements to local ids.                                         */
/*                                                                          */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
/*             Karen D. Devine (9226)                                       */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Revision History:                                                        */
/*    10 May 1999:       Date of creation.                                  */
/*--------------------------------------------------------------------------*/

#include <mpi.h>
#include <stdlib.h>

#include "dr_const.h"
#include "dr_externs.h"
#include "dr_err_const.h"
#include "dr_loadbal_const.h"
#include "dr_par_util_const.h"
#include "dr_util_const.h"
#include "dr_output_const.h"
#include "dr_elem_util_const.h"
#include "dr_maps_const.h"
#include "dr_dd.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

extern int my_rank;

/*
 *  PROTOTYPES for load-balancer interface functions.
 */
ZOLTAN_PRE_MIGRATE_PP_FN migrate_pre_process;
ZOLTAN_POST_MIGRATE_PP_FN migrate_post_process;

/* Object-based migration callbacks; only one of these or the list-based
 * callbacks are actually needed. */
ZOLTAN_PACK_OBJ_FN migrate_pack_elem;
ZOLTAN_UNPACK_OBJ_FN migrate_unpack_elem;

/* List-based migration callbacks; only one of these or the object-based
 * callbacks are actually needed. */
ZOLTAN_PACK_OBJ_MULTI_FN migrate_pack_elem_multi;
ZOLTAN_UNPACK_OBJ_MULTI_FN migrate_unpack_elem_multi;

/*****************************************************************************/
/*
 *  Static global variables to help with migration.
 */
struct New_Elem_Hash_Node{
  ZOLTAN_ID_TYPE globalID;
  int localID;
  int next;
};

static int *New_Elem_Hash_Table = NULL;
static struct New_Elem_Hash_Node *New_Elem_Hash_Nodes = NULL;
                                      /* Hash table containing globalIDs and
                                         localIDs of elements in the new
                                         decomposition; used for quick
                                         globalID -> localID lookup. */

static ZOLTAN_ID_TYPE *New_Elem_Index = NULL; /* Array containing globalIDs of 
                                         elements in the new decomposition,
                                         ordered in the same order as the
                                         elements array.
                                         Built in migrate_pre_process; used
                                         in migrate_pre_process to adjust
                                         element adjacencies; used in 
                                         migrate_unpack_elem to store 
                                         imported elements.                  */
static int New_Elem_Index_Size = 0;   /* Number of integers allocated in
                                         New_Elem_Index.                     */
static int Use_Edge_Wgts = 0;         /* Flag indicating whether elements
                                         store edge weights.                 */
/*static int Vertex_Blanking = 0;        We're dynamically altering the graph
                                         in each iteration by blanking portions
                                         of it, so we must migrate flags 
                                         indicating whether adjacent vertices
                                         of migrated elements were blanked on
                                         the originating process.             */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int migrate_elements(
  int Proc,
  MESH_INFO_PTR mesh,
  struct Zoltan_Struct *zz,
  int num_gid_entries, 
  int num_lid_entries,
  int num_imp,
  ZOLTAN_ID_PTR imp_gids,
  ZOLTAN_ID_PTR imp_lids,
  int *imp_procs,
  int *imp_to_part,
  int num_exp,
  ZOLTAN_ID_PTR exp_gids,
  ZOLTAN_ID_PTR exp_lids,
  int *exp_procs,
  int *exp_to_part)
{
/* Local declarations. */
char *yo = "migrate_elements";

/***************************** BEGIN EXECUTION ******************************/
  DEBUG_TRACE_START(Proc, yo);

  /*
   * register migration functions
   */
  if (!Test.Null_Lists) {
    /* If not passing NULL lists, let Help_Migrate call the
     * pre-processing and post-processing routines.
     */
    if (Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, 
                      (void (*)()) migrate_pre_process,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }

    if (Zoltan_Set_Fn(zz, ZOLTAN_POST_MIGRATE_PP_FN_TYPE, 
                      (void (*)()) migrate_post_process,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }
  }

  if (Test.Multi_Callbacks) {
    if (Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE, 
                      (void (*)()) migrate_elem_size_multi,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }

    if (Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_MULTI_FN_TYPE,
                      (void (*)()) migrate_pack_elem_multi,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }
  
    if (Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE,
                      (void (*)()) migrate_unpack_elem_multi,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }
  }
  else {
    if (Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_FN_TYPE,
                      (void (*)()) migrate_elem_size,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }

    if (Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_FN_TYPE, 
                      (void (*)()) migrate_pack_elem,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }

    if (Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_FN_TYPE,
                      (void (*)()) migrate_unpack_elem,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }
  }


  if (Test.Null_Lists == NO_NULL_LISTS) {
    if (Zoltan_Migrate(zz, num_imp, imp_gids, imp_lids, imp_procs, imp_to_part,
                           num_exp, exp_gids, exp_lids, exp_procs, exp_to_part)
      == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Migrate()\n");
      return 0;
    }
  }
  else {
    /* Call Zoltan_Help_Migrate with empty import lists. */
    /* Have to "manually" call migrate_pre_process and migrate_post_process. */
    int ierr = 0;
    migrate_pre_process((void *) mesh, 1, 1,
                        num_imp, imp_gids, imp_lids, imp_procs, imp_to_part,
                        num_exp, exp_gids, exp_lids, exp_procs, exp_to_part,
                        &ierr);
    if (Test.Null_Lists == IMPORT_LISTS) {
      if (Zoltan_Migrate(zz,
                         -1, NULL, NULL, NULL, NULL,
                         num_exp, exp_gids, exp_lids, exp_procs, exp_to_part)
          == ZOLTAN_FATAL) {
        Gen_Error(0, "fatal:  error returned from Zoltan_Migrate()\n");
        return 0;
      }
    }
    else {
      if (Zoltan_Migrate(zz,
                         num_imp, imp_gids, imp_lids, imp_procs, imp_to_part,
                         -1, NULL, NULL, NULL, NULL)
          == ZOLTAN_FATAL) {
        Gen_Error(0, "fatal:  error returned from Zoltan_Migrate()\n");
        return 0;
      }
    }
    migrate_post_process((void *) mesh, 1, 1,  
                         num_imp, imp_gids, imp_lids, imp_procs, imp_to_part,
                         num_exp, exp_gids, exp_lids, exp_procs, exp_to_part,
                         &ierr);
  }

  DEBUG_TRACE_END(Proc, yo);
  return 1;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
extern unsigned int Zoltan_Hash(ZOLTAN_ID_PTR, int, unsigned int);

void insert_in_hash(
  ZOLTAN_ID_TYPE globalID,
  int localID
)
{
int j;

  New_Elem_Hash_Nodes[localID].globalID = globalID;
  New_Elem_Hash_Nodes[localID].localID = localID;
  j = Zoltan_Hash(&globalID,1,New_Elem_Index_Size);
  New_Elem_Hash_Nodes[localID].next = New_Elem_Hash_Table[j];
  New_Elem_Hash_Table[j] = localID;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int find_in_hash(
  ZOLTAN_ID_TYPE globalID
) 
{
int idx;

  idx = Zoltan_Hash(&globalID, 1, New_Elem_Index_Size);
  idx = New_Elem_Hash_Table[idx];
  while (idx != -1 && New_Elem_Hash_Nodes[idx].globalID != globalID) {
    idx = New_Elem_Hash_Nodes[idx].next;
  }

  return idx;
}

void remove_from_hash(
  ZOLTAN_ID_TYPE globalID
)
{
int idx, hidx, prev;

  hidx = Zoltan_Hash(&globalID, 1, New_Elem_Index_Size);
  idx = New_Elem_Hash_Table[hidx];
  prev = -1;
  while (idx != -1 && New_Elem_Hash_Nodes[idx].globalID != globalID) {
    prev = idx;
    idx = New_Elem_Hash_Nodes[idx].next;
  }
  if (prev == -1) 
    New_Elem_Hash_Table[hidx] = New_Elem_Hash_Nodes[idx].next;
  else
    New_Elem_Hash_Nodes[prev].next = New_Elem_Hash_Nodes[idx].next;

  New_Elem_Hash_Nodes[idx].globalID = ZOLTAN_ID_INVALID;
  New_Elem_Hash_Nodes[idx].localID = -1;
  New_Elem_Hash_Nodes[idx].next = -1;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void migrate_pre_process(void *data, int num_gid_entries, int num_lid_entries, 
                         int num_import, 
                         ZOLTAN_ID_PTR import_global_ids,
                         ZOLTAN_ID_PTR import_local_ids, int *import_procs,
                         int *import_to_part,
                         int num_export, ZOLTAN_ID_PTR export_global_ids,
                         ZOLTAN_ID_PTR export_local_ids, int *export_procs,
                         int *export_to_part,
                         int *ierr)
{
int i, j, k, idx, maxlen, proc, offset;
int *proc_ids = NULL;   /* Temp array of processor assignments for elements.*/
char *change = NULL;    /* Temp array indicating whether local element's adj 
                           list must be updated due to a nbor's migration.  */
int new_proc;           /* New processor assignment for nbor element.       */
int exp_elem;           /* index of an element being exported */
int bor_elem;           /* index of an element along the processor border */
int *send_vec = NULL, *recv_vec = NULL;  /* Communication vecs. */
MESH_INFO_PTR mesh;
ELEM_INFO_PTR elements;
int lid = num_lid_entries-1;
int gid = num_gid_entries-1;
char msg[256];

  *ierr = ZOLTAN_OK;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  mesh = (MESH_INFO_PTR) data;
  elements = mesh->elements;

  for (i=0; i < mesh->num_elems; i++){
    /* don't migrate a pointer created on this process */
    safe_free((void **)(void *)&(elements[i].adj_blank));
  }

  /*
   *  Set some flags. Assume if true for one element, true for all elements.
   *  Note that some procs may have no elements. 
   */

  if (elements[0].edge_wgt != NULL)
    k = 1;
  else
    k = 0;
  /* Make sure all procs have the same value */
  MPI_Allreduce(&k, &Use_Edge_Wgts, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* NOT IMPLEMENTED: blanking information is not sent along.  Subsequent
     lb_eval may be incorrect, since imported elements may have blanked
     adjacencies.

  if (mesh->blank_count > 0)
    k = 1;
  else
    k = 0;
 
  MPI_Allreduce(&k, &Vertex_Blanking, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  */

  /*
   *  For all elements, update adjacent elements' processor information.
   *  That way, when perform migration, will be migrating updated adjacency
   *  information.  
   */
  
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  /*
   *  Build New_Elem_Index array and list of processor assignments.
   */

  New_Elem_Index_Size = mesh->num_elems + num_import - num_export;

  if (mesh->elem_array_len > New_Elem_Index_Size) 
    New_Elem_Index_Size = mesh->elem_array_len;

  New_Elem_Index = (ZOLTAN_ID_TYPE *) malloc(New_Elem_Index_Size * sizeof(ZOLTAN_ID_TYPE));
  New_Elem_Hash_Table = (int *) malloc(New_Elem_Index_Size * sizeof(int));
  New_Elem_Hash_Nodes = (struct New_Elem_Hash_Node *) 
           malloc(New_Elem_Index_Size * sizeof(struct New_Elem_Hash_Node));

  if (New_Elem_Index == NULL || 
      New_Elem_Hash_Table == NULL || New_Elem_Hash_Nodes == NULL) {
    Gen_Error(0, "fatal: insufficient memory");
    *ierr = ZOLTAN_MEMERR;
    return;
  }

  for (i = 0; i < New_Elem_Index_Size; i++) 
    New_Elem_Hash_Table[i] = -1;
  for (i = 0; i < New_Elem_Index_Size; i++) {
    New_Elem_Hash_Nodes[i].globalID = ZOLTAN_ID_INVALID;
    New_Elem_Hash_Nodes[i].localID = -1;
    New_Elem_Hash_Nodes[i].next = -1;
  }

  if (mesh->num_elems > 0) {

    proc_ids = (int *)  malloc(mesh->num_elems * sizeof(int));
    change   = (char *) malloc(mesh->num_elems * sizeof(char));

    if (New_Elem_Index == NULL || proc_ids == NULL || change == NULL ||
        New_Elem_Hash_Table == NULL || New_Elem_Hash_Nodes == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      *ierr = ZOLTAN_MEMERR;
      return;
    }

    for (i = 0; i < mesh->num_elems; i++) {
      New_Elem_Index[i] = elements[i].globalID;
      insert_in_hash(elements[i].globalID, i);
      proc_ids[i] = proc;
      change[i] = 0;
    }
  }

  for (i = mesh->num_elems; i < New_Elem_Index_Size; i++) {
    New_Elem_Index[i] = ZOLTAN_ID_INVALID;
  }

  for (i = 0; i < num_export; i++) {
    if (num_lid_entries)
      exp_elem = export_local_ids[lid+i*num_lid_entries];
    else  /* testing num_lid_entries == 0 */
      search_by_global_id(mesh, export_global_ids[gid+i*num_gid_entries], 
                          &exp_elem);

    if (export_procs[i] != proc) {
      /* Export is moving to a new processor */
      New_Elem_Index[exp_elem] = ZOLTAN_ID_INVALID;
      remove_from_hash(export_global_ids[gid+i*num_gid_entries]);
      proc_ids[exp_elem] = export_procs[i];
    }
  }

  j = 0;
  for (i = 0; i < num_import; i++) {
    if (import_procs[i] != proc) {
      /* Import is moving from a new processor, not just from a new partition */
      /* search for first free location */
      for ( ; j < New_Elem_Index_Size; j++) 
        if (New_Elem_Index[j] == ZOLTAN_ID_INVALID) break;

      New_Elem_Index[j] = import_global_ids[gid+i*num_gid_entries];
      insert_in_hash(import_global_ids[gid+i*num_gid_entries], j);
    }
  }

  /* 
   * Update local information 
   */

  /* Set change flag for elements whose adjacent elements are being exported */

  for (i = 0; i < num_export; i++) {

    if (num_lid_entries)
      exp_elem = export_local_ids[lid+i*num_lid_entries];
    else  /* testing num_lid_entries == 0 */
      search_by_global_id(mesh, export_global_ids[gid+i*num_gid_entries], 
                          &exp_elem);

    elements[exp_elem].my_part = export_to_part[i];

    if (export_procs[i] == proc) 
      continue;  /* No adjacency changes needed if export is changing
                    only partition, not processor. */

    for (j = 0; j < elements[exp_elem].adj_len; j++) {

      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (elements[exp_elem].adj[j] == ZOLTAN_ID_INVALID) continue;

      /* Set change flag for adjacent local elements. */
      if (elements[exp_elem].adj_proc[j] == proc) {
        change[elements[exp_elem].adj[j]] = 1;
      }
    }
  }

  /* Change adjacency information in marked elements */
  for (i = 0; i < mesh->num_elems; i++) {
    if (change[i] == 0) continue;

    /* loop over marked element's adjacencies; look for ones that are moving */
    for (j = 0; j < elements[i].adj_len; j++) {

      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (elements[i].adj[j] == ZOLTAN_ID_INVALID) continue;

      if (elements[i].adj_proc[j] == proc) {
        /* adjacent element is local; check whether it is moving. */
        if ((new_proc = proc_ids[elements[i].adj[j]]) != proc) {
          /* Adjacent element is being exported; update this adjacency entry */
          elements[i].adj[j] = elements[elements[i].adj[j]].globalID;
          elements[i].adj_proc[j] = new_proc;
        }
      }
    }
  }
  safe_free((void **)(void *) &change);

  /*
   * Update off-processor information 
   */

  maxlen = 0;
  for (i = 0; i < mesh->necmap; i++) 
    maxlen += mesh->ecmap_cnt[i];

  if (maxlen > 0) {
    send_vec = (int *) malloc(maxlen * sizeof(int));
    if (send_vec == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      *ierr = ZOLTAN_MEMERR;
      return;
    }

    /* Load send vector */

    for (i = 0; i < maxlen; i++)
      send_vec[i] = proc_ids[mesh->ecmap_elemids[i]];
  }

  safe_free((void **)(void *) &proc_ids);

  if (maxlen > 0)
    recv_vec = (int *) malloc(maxlen * sizeof(int));

  /*  Perform boundary exchange */

  boundary_exchange(mesh, 1, send_vec, recv_vec);
  
  /* Unload receive vector */

  offset = 0;
  for (i = 0; i < mesh->necmap; i++) {
    for (j = 0; j < mesh->ecmap_cnt[i]; j++, offset++) {
      if (recv_vec[offset] == mesh->ecmap_id[i]) {
        /* off-processor element is not changing processors.  */
        /* no changes are needed in the local data structure. */
        continue;
      }
      /* Change processor assignment in local element's adjacency list */
      bor_elem = mesh->ecmap_elemids[offset];
      for (k = 0; k < elements[bor_elem].adj_len; k++) {

        /* Skip NULL adjacencies (sides that are not adj to another elem). */
        if (elements[bor_elem].adj[k] == ZOLTAN_ID_INVALID) continue;

        if (elements[bor_elem].adj[k] == mesh->ecmap_neighids[offset] &&
            elements[bor_elem].adj_proc[k] == mesh->ecmap_id[i]) {
          elements[bor_elem].adj_proc[k] = recv_vec[offset];
          if (recv_vec[offset] == proc) {
            /* element is moving to this processor; */
            /* convert adj from global to local ID. */
            idx = find_in_hash(mesh->ecmap_neighids[offset]);
            if (idx >= 0) 
              idx = New_Elem_Hash_Nodes[idx].localID;
            else {
              sprintf(msg, "fatal: unable to locate element " ZOLTAN_ID_SPEC " in "
                           "New_Elem_Index", mesh->ecmap_neighids[offset]);
              Gen_Error(0, msg);
              *ierr = ZOLTAN_FATAL;
              return;
            }
            elements[bor_elem].adj[k] = (ZOLTAN_ID_TYPE)idx;
          }
          break;  /* from k loop */
        }
      }
    }
  }

  safe_free((void **)(void *) &recv_vec);
  safe_free((void **)(void *) &send_vec);

  /*
   * Allocate space (if needed) for the new element data.
   */

  if (mesh->elem_array_len < New_Elem_Index_Size) {
    mesh->elem_array_len = New_Elem_Index_Size;
    mesh->elements = (ELEM_INFO_PTR) realloc (mesh->elements,
                                     mesh->elem_array_len * sizeof(ELEM_INFO));
    if (mesh->elements == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      return;
    }

    /* initialize the new spots */
    for (i = mesh->num_elems; i < mesh->elem_array_len; i++)
      initialize_element(&(mesh->elements[i]));
  }
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void migrate_post_process(void *data, int num_gid_entries, int num_lid_entries,
                          int num_import, 
                          ZOLTAN_ID_PTR import_global_ids,
                          ZOLTAN_ID_PTR import_local_ids, int *import_procs,
                          int *import_to_part,
                          int num_export, ZOLTAN_ID_PTR export_global_ids,
                          ZOLTAN_ID_PTR export_local_ids, int *export_procs,
                          int *export_to_part,
                          int *ierr)
{
MESH_INFO_PTR mesh;
ELEM_INFO *elements;
int proc, num_proc;
int i, j, k, last;
ZOLTAN_ID_TYPE adj_elem;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  mesh = (MESH_INFO_PTR) data;
  elements = mesh->elements;


  MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

  /* compact elements array, as the application expects the array to be dense */
  for (i = 0; i < New_Elem_Index_Size; i++) {
    if (New_Elem_Index[i] != ZOLTAN_ID_INVALID) continue;

    /* Don't want to shift all elements down one position to fill the  */
    /* blank spot -- too much work to adjust adjacencies!  So find the */
    /* last element in the array and move it to the blank spot.        */

    for (last = New_Elem_Index_Size-1; last >= 0; last--)
      if (New_Elem_Index[last] != ZOLTAN_ID_INVALID) break;

    /* If (last < i), array is already dense; i is just in some blank spots  */
    /* at the end of the array.  Quit the compacting.                     */
    if (last < i) break;

    /* Copy elements[last] to elements[i]. */
    elements[i] = elements[last];

    /* Adjust adjacencies for local elements.  Off-processor adjacencies */
    /* don't matter here.                                                */

    for (j = 0; j < elements[i].adj_len; j++) {

      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (elements[i].adj[j] == ZOLTAN_ID_INVALID) continue;

      adj_elem = elements[i].adj[j];

      /* See whether adjacent element is local; if so, adjust its entry */
      /* for local element i.                                           */
      if (elements[i].adj_proc[j] == proc) {
        for (k = 0; k < elements[adj_elem].adj_len; k++) {
          if (elements[adj_elem].adj[k] == (ZOLTAN_ID_TYPE)last &&
              elements[adj_elem].adj_proc[k] == proc) {
            /* found adjacency entry for element last; change it to i */
            elements[adj_elem].adj[k] = (ZOLTAN_ID_TYPE)i;
            break;
          }
        }
      }
    }

    /* Update New_Elem_Index */
    New_Elem_Index[i] = New_Elem_Index[last];
    New_Elem_Index[last] = ZOLTAN_ID_INVALID;

    /* clear elements[last] */
    elements[last].globalID = ZOLTAN_ID_INVALID;
    elements[last].border = 0;
    elements[last].my_part = -1;
    elements[last].perm_value = -1;
    elements[last].invperm_value = -1;
    elements[last].nadj = 0;
    elements[last].adj_len = 0;
    elements[last].elem_blk = -1;
    for (k=0; k<MAX_CPU_WGTS; k++)
      elements[last].cpu_wgt[k] = 0;
    elements[last].mem_wgt = 0;
    elements[last].avg_coord[0] = elements[last].avg_coord[1] 
                                = elements[last].avg_coord[2] = 0.;
    elements[last].coord = NULL;
    elements[last].connect = NULL;
    elements[last].adj = NULL;
    elements[last].adj_proc = NULL;
    elements[last].edge_wgt = NULL;
  }

  if (New_Elem_Index != NULL) safe_free((void **)(void *) &New_Elem_Index);
  if (New_Elem_Hash_Table != NULL) safe_free((void **)(void *) &New_Elem_Hash_Table);
  if (New_Elem_Hash_Nodes != NULL) safe_free((void **)(void *) &New_Elem_Hash_Nodes);
  New_Elem_Index_Size = 0;

  if (!build_elem_comm_maps(proc, mesh)) {
    Gen_Error(0, "Fatal: error rebuilding elem comm maps");
  }

  if (mesh->data_type == ZOLTAN_HYPERGRAPH && !update_elem_dd(mesh)) {
    Gen_Error(0, "Fatal: error updating element dd");
  }

  if (mesh->data_type == ZOLTAN_HYPERGRAPH && mesh->hvertex_proc &&
      !update_hvertex_proc(mesh)) {
    Gen_Error(0, "Fatal: error updating hyperedges");
  }

  if (Debug_Driver > 3)
    print_distributed_mesh(proc, num_proc, mesh);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int migrate_elem_size(void *data, int num_gid_entries, int num_lid_entries,
    ZOLTAN_ID_PTR elem_gid, ZOLTAN_ID_PTR elem_lid, int *ierr)
/*
 * Function to return size of element information for a single element.
 */
{
int size;
int num_nodes;
MESH_INFO_PTR mesh;
ELEM_INFO *current_elem;
int gid = num_gid_entries-1;
int lid = num_lid_entries-1;
int idx;

  *ierr = ZOLTAN_OK;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  mesh = (MESH_INFO_PTR) data;
  current_elem = (num_lid_entries 
                   ? &(mesh->elements[elem_lid[lid]])
                   : search_by_global_id(mesh, elem_gid[gid], &idx));
  num_nodes = mesh->eb_nnodes[current_elem->elem_blk];

  /*
   * Compute an upper bound of the size of one element's data.
   *   Some values are hard-coded to an upper bound because we want
   *   to get the same test answers on 32-bit arches that we get on
   *   64-bit arches, and we want the same answer whether we have
   *   32-bit ZOLTAN_ID_TYPEs or 64-bit ZOLTAN_ID_TYPEs.
   */

  /* Using 200 instead of sizeof(ELEM_INFO) */

  if (sizeof(ELEM_INFO) > 200){
    fprintf(stderr,"Re-code migrate_elem_size\n");
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  size = 200;
 
  /* Add space to correct alignment so casts work in (un)packing. */
  size = Zoltan_Align(size);

  /* Add space for connect table.  
   * Using "8" instead of  sizeof(ZOLTAN_ID_TYPE). */
   
  if (mesh->num_dims > 0)
    size += num_nodes * 8;

  /* Add space for adjacency info (elements[].adj)
   * Using "8" instead of  sizeof(ZOLTAN_ID_TYPE). */
   
  size += current_elem->adj_len * 8;

  /* Add space for adjacency info (elements[].adj_proc). */
  size += current_elem->adj_len * sizeof(int);

  /* Assume if one element has edge wgts, all elements have edge wgts. */
  if (Use_Edge_Wgts) {
    /* Add space to correct alignment so casts work in (un)packing. */
    size = Zoltan_Align(size);
    size += current_elem->adj_len * sizeof(float);
  }

  /* Add space for coordinate info */
  size = Zoltan_Align(size);
  size += num_nodes * mesh->num_dims * sizeof(float);
  
  /* For dynamic weights test, multiply size by vertex weight. */
  /* This simulates mesh refinement. */
  if (Test.Dynamic_Weights){
    size *= ((current_elem->cpu_wgt[0] > 1.0) ? current_elem->cpu_wgt[0] : 1.0);
  }

  return (size);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void migrate_pack_elem(void *data, int num_gid_entries, int num_lid_entries,
                       ZOLTAN_ID_PTR elem_gid, ZOLTAN_ID_PTR elem_lid,
                       int mig_part, int elem_data_size, char *buf, int *ierr)
{
  MESH_INFO_PTR mesh;
  ELEM_INFO *elem, *elem_mig;
  ELEM_INFO *current_elem;
  int *buf_int;
  ZOLTAN_ID_TYPE *buf_id_type;
  float *buf_float;
  int size;
  int i, j, idx;
  int proc;
  int num_nodes;
  int gid = num_gid_entries-1;
  int lid = num_lid_entries-1;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  mesh = (MESH_INFO_PTR) data;
  elem = mesh->elements;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  current_elem = (num_lid_entries 
                   ? &(elem[elem_lid[lid]])
                   : search_by_global_id(mesh, elem_gid[gid], &idx));
  num_nodes = mesh->eb_nnodes[current_elem->elem_blk];

  elem_mig = (ELEM_INFO *) buf; /* this is the element struct to be migrated */

  /*
   * copy the ELEM_INFO structure
   */
  *elem_mig = *current_elem;
  size = sizeof(ELEM_INFO);

  /*
   * copy the allocated integer fields for this element.
   */

  /* Pad the buffer so the following casts will work.  */
  
  size = Zoltan_Align(size);

  buf_id_type = (ZOLTAN_ID_TYPE *) (buf + size);

  /* copy the connect table */
  if (mesh->num_dims > 0) {
    for (i = 0; i < num_nodes; i++) {
      *buf_id_type = current_elem->connect[i];
      buf_id_type++;
    }
    size += num_nodes * sizeof(ZOLTAN_ID_TYPE);
  }

  /* copy the adjacency info - all global IDs, then all processes owning global IDs */

  for (i =  0; i < current_elem->adj_len; i++) {
    if (current_elem->adj[i] != ZOLTAN_ID_INVALID && current_elem->adj_proc[i] == proc) 
      *buf_id_type++ = New_Elem_Index[current_elem->adj[i]];
    else
      *buf_id_type++ = current_elem->adj[i];
  }
  size += current_elem->adj_len * sizeof(ZOLTAN_ID_TYPE);

  buf_int = (int *) (buf + size);

  for (i =  0; i < current_elem->adj_len; i++) {
    *buf_int++ = current_elem->adj_proc[i];
  }
  size += current_elem->adj_len * sizeof(int);

  /*
   * copy the allocated float fields for this element.
   */

  /* copy the edge_wgt data */
  if (Use_Edge_Wgts) {

    /* Pad the buffer so the following casts will work.  */
    size = Zoltan_Align(size);
    buf_float = (float *) (buf + size);

    for (i = 0; i < current_elem->adj_len; i++) {
      *buf_float = current_elem->edge_wgt[i];
      buf_float++;
    }
    size += current_elem->adj_len * sizeof(float);
  }

  /* Pad the buffer so the following casts will work.  */
  size = Zoltan_Align(size);
  buf_float = (float *) (buf + size);

  /* copy coordinate data */
  for (i = 0; i < num_nodes; i++) {
    for (j = 0; j < mesh->num_dims; j++) {
      *buf_float = current_elem->coord[i][j];
      buf_float++;
    }
  }
  size += num_nodes * mesh->num_dims * sizeof(float);

  /*
   * need to update the Mesh struct to reflect this element
   * being gone
   */
  mesh->num_elems--;
  mesh->eb_cnts[current_elem->elem_blk]--;

  /*
   * need to remove this entry from this procs list of elements
   * do so by setting the globalID to ZOLTAN_ID_INVALID. 
   */
  current_elem->globalID = ZOLTAN_ID_INVALID;
  free_element_arrays(current_elem, mesh);

  /*
   * NOTE: it is not worth the effort to determine the change in the
   * number of nodes on this processor until all of the migration is
   * completed.
   */
  if (size > elem_data_size) 
    *ierr = ZOLTAN_WARN;
  else
    *ierr = ZOLTAN_OK;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void migrate_unpack_elem(void *data, int num_gid_entries, ZOLTAN_ID_PTR elem_gid, 
                         int elem_data_size, char *buf, int *ierr)
{
  MESH_INFO_PTR mesh;
  ELEM_INFO *elem, *elem_mig;
  ELEM_INFO *current_elem;
  ZOLTAN_ID_TYPE *buf_id_type;
  int *buf_int;
  float *buf_float;
  int size, num_nodes;
  int i, j, idx;
  int proc;
  int gid = num_gid_entries-1;
  size_t adjgids_len, adjprocs_len;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  mesh = (MESH_INFO_PTR) data;
  elem = mesh->elements;
  elem_mig = (ELEM_INFO *) buf;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  idx = find_in_hash((int)elem_gid[gid]);
  if (idx >= 0) 
    idx = New_Elem_Hash_Nodes[idx].localID;
  else {
    Gen_Error(0, "fatal: Unable to locate position for element");
    *ierr = ZOLTAN_FATAL;
    return;
  }

  current_elem = &(elem[idx]);
  /* now put the migrated information into the array */
  *current_elem = *elem_mig;
  num_nodes = mesh->eb_nnodes[current_elem->elem_blk];

  size = sizeof(ELEM_INFO);

  /*
   * copy the allocated integer fields for this element.
   */

  /* Pad the buffer so the following casts will work.  */
  size = Zoltan_Align(size);
  buf_id_type = (ZOLTAN_ID_TYPE *) (buf + size);

  /* copy the connect table */
  if (mesh->num_dims > 0) {
    current_elem->connect = (ZOLTAN_ID_TYPE *) malloc(num_nodes * sizeof(ZOLTAN_ID_TYPE));
    if (current_elem->connect == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      *ierr = ZOLTAN_MEMERR;
      return;
    }
    for (i = 0; i < num_nodes; i++) {
      current_elem->connect[i] = *buf_id_type++;
    }
    size += num_nodes * sizeof(ZOLTAN_ID_TYPE);
  }

  /* copy the adjacency info */
  /* globalIDs are received; convert to local IDs when adj elem is local */

  if (current_elem->adj_len > 0) {

    adjgids_len = current_elem->adj_len * sizeof(ZOLTAN_ID_TYPE);
    adjprocs_len = current_elem->adj_len * sizeof(int);

    current_elem->adj      = (ZOLTAN_ID_TYPE *)malloc(adjgids_len);
    current_elem->adj_proc = (int *)malloc(adjprocs_len);

    if (current_elem->adj == NULL || current_elem->adj_proc == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      *ierr = ZOLTAN_MEMERR;
      return;
    }

    buf_id_type = (ZOLTAN_ID_TYPE *) (buf + size);
    buf_int    = (int *) (buf + size + adjgids_len);

    for (i =  0; i < current_elem->adj_len; i++) {
      current_elem->adj[i] =      *buf_id_type++;
      current_elem->adj_proc[i] = *buf_int++;

      if (current_elem->adj[i] != ZOLTAN_ID_INVALID && current_elem->adj_proc[i] == proc) {
        int nidx = find_in_hash(current_elem->adj[i]);
        if (nidx >= 0) 
          current_elem->adj[i] = (ZOLTAN_ID_TYPE)New_Elem_Hash_Nodes[nidx].localID;
        else {
          Gen_Error(0, "fatal: Unable to locate position for neighbor");
          *ierr = ZOLTAN_FATAL;
          return;
        }     
      }
    }
    size += (adjgids_len + adjprocs_len);

    /* copy the edge_wgt data */
    if (Use_Edge_Wgts) {

      /* Pad the buffer so the following casts will work.  */
      size = Zoltan_Align(size);
      buf_float = (float *) (buf + size);

      current_elem->edge_wgt = (float *) malloc(current_elem->adj_len 
                                              * sizeof(float));
      if (current_elem->edge_wgt == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        *ierr = ZOLTAN_MEMERR;
        return;
      }
      for (i = 0; i < current_elem->adj_len; i++) {
        current_elem->edge_wgt[i] = *buf_float;
        buf_float++;
      }
      size += current_elem->adj_len * sizeof(float);
    }
  }

  /* copy coordinate data */
  if (num_nodes > 0) {

    /* Pad the buffer so the following casts will work.  */
    size = Zoltan_Align(size);
    buf_float = (float *) (buf + size);

    current_elem->coord = (float **) malloc(num_nodes * sizeof(float *));
    if (current_elem->coord == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      *ierr = ZOLTAN_MEMERR;
      return;
    }
    for (i = 0; i < num_nodes; i++) {
      current_elem->coord[i] = (float *) malloc(mesh->num_dims * sizeof(float));
      if (current_elem->coord[i] == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        *ierr = ZOLTAN_MEMERR;
        return;
      }
      for (j = 0; j < mesh->num_dims; j++) {
        current_elem->coord[i][j] = *buf_float;
        buf_float++;
      }
    }
    size += num_nodes * mesh->num_dims * sizeof(float);
  }


  /* and update the Mesh struct */
  mesh->num_elems++;
  mesh->eb_cnts[current_elem->elem_blk]++;

  if (size > elem_data_size) 
    *ierr = ZOLTAN_WARN;
  else
    *ierr = ZOLTAN_OK;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void migrate_elem_size_multi(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_ids,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int *num_bytes,
  int *ierr
)
{
int i;

  *ierr = ZOLTAN_OK;
  for (i = 0; i < num_ids; i++) {
    num_bytes[i] = migrate_elem_size(data, num_gid_entries, num_lid_entries,
                  &(global_ids[i*num_gid_entries]),
                  (num_lid_entries!=0 ? &(local_ids[i*num_lid_entries]) : NULL),
                  ierr);
    if (*ierr != ZOLTAN_OK)
      return;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void migrate_pack_elem_multi(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_ids,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int *dest_proc,
  int *size,
  int *index,
  char *buffer,
  int *ierr
)
{
int i;

  *ierr = ZOLTAN_OK;
  for (i = 0; i < num_ids; i++) {
    migrate_pack_elem(data, num_gid_entries, num_lid_entries,
            &(global_ids[i*num_gid_entries]),
            (num_lid_entries!=0 ? &(local_ids[i*num_lid_entries]) : NULL),
            dest_proc[i], size[i], &(buffer[index[i]]), ierr);
    if (*ierr != ZOLTAN_OK)
      return;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void migrate_unpack_elem_multi(
  void *data,
  int num_gid_entries,
  int num_ids,
  ZOLTAN_ID_PTR global_ids,
  int *size,
  int *index,
  char *buffer,
  int *ierr
)
{
int i;

  *ierr = ZOLTAN_OK;
  for (i = 0; i < num_ids; i++) {
    migrate_unpack_elem(data, num_gid_entries,
                        &(global_ids[i*num_gid_entries]),
                        size[i], &(buffer[index[i]]), ierr);
    if (*ierr != ZOLTAN_OK)
      return;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
