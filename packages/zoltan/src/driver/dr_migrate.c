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

/*define NULL_IMPORT_LISTS_TO_HELP_MIGRATE*/
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

#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "dr_const.h"
#include "dr_err_const.h"
#include "dr_loadbal_const.h"
#include "dr_par_util_const.h"
#include "dr_util_const.h"
#include "dr_output_const.h"
#include "dr_elem_util_const.h"
#include "dr_maps_const.h"

/*
 *  PROTOTYPES for load-balancer interface functions.
 */
ZOLTAN_PRE_MIGRATE_FN migrate_pre_process;
ZOLTAN_POST_MIGRATE_FN migrate_post_process;
ZOLTAN_OBJ_SIZE_FN migrate_elem_size;
ZOLTAN_PACK_OBJ_FN migrate_pack_elem;
ZOLTAN_UNPACK_OBJ_FN migrate_unpack_elem;

/*****************************************************************************/
/*
 *  Static global variables to help with migration.
 */
static int *New_Elem_Index = NULL;    /* Array containing globalIDs of 
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
  int num_exp,
  ZOLTAN_ID_PTR exp_gids,
  ZOLTAN_ID_PTR exp_lids,
  int *exp_procs)
{
/* Local declarations. */
char *yo = "migrate_elements";

/***************************** BEGIN EXECUTION ******************************/
  DEBUG_TRACE_START(Proc, yo);

  /*
   * register migration functions
   */
#ifndef NULL_IMPORT_LISTS_TO_HELP_MIGRATE
  if (Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_FN_TYPE, (void (*)()) migrate_pre_process,
                (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_POST_MIGRATE_FN_TYPE, (void (*)()) migrate_post_process,
                (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }
#endif

  if (Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_FN_TYPE, (void (*)()) migrate_elem_size,
               (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_FN_TYPE, (void (*)()) migrate_pack_elem,
                (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_FN_TYPE, (void (*)()) migrate_unpack_elem,
                (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

#ifndef NULL_IMPORT_LISTS_TO_HELP_MIGRATE
  if (Zoltan_Help_Migrate(zz, 
                      num_imp, imp_gids, imp_lids, imp_procs,
                      num_exp, exp_gids, exp_lids, exp_procs) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Help_Migrate()\n");
    return 0;
  }
#else
  {
  /* Call Zoltan_Help_Migrate with empty import lists. */
  /* Have to "manually" call migrate_pre_process and migrate_post_process. */
  int ierr = 0;
  migrate_pre_process((void *) mesh, 1, 1,
                      num_imp, imp_gids, imp_lids, imp_procs,
                      num_exp, exp_gids, exp_lids, exp_procs, &ierr);
  if (Zoltan_Help_Migrate(zz, 
                      -1, NULL, NULL, NULL,
                      num_exp, exp_gids, exp_lids, exp_procs) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Help_Migrate()\n");
    return 0;
  }
  migrate_post_process((void *) mesh, 1, 1,  
                       num_imp, imp_gids, imp_lids, imp_procs,
                       num_exp, exp_gids, exp_lids, exp_procs, &ierr);
  }
#endif


  DEBUG_TRACE_END(Proc, yo);
  return 1;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void migrate_pre_process(void *data, int num_gid_entries, int num_lid_entries, 
                               int num_import, 
                               ZOLTAN_ID_PTR import_global_ids,
                               ZOLTAN_ID_PTR import_local_ids, int *import_procs,
                               int num_export, ZOLTAN_ID_PTR export_global_ids,
                               ZOLTAN_ID_PTR export_local_ids, int *export_procs,
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
  New_Elem_Index = (int *) malloc(New_Elem_Index_Size * sizeof(int));

  if (mesh->num_elems > 0) {

    proc_ids = (int *)  malloc(mesh->num_elems * sizeof(int));
    change   = (char *) malloc(mesh->num_elems * sizeof(char));

    if (New_Elem_Index == NULL || proc_ids == NULL || change == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      *ierr = ZOLTAN_MEMERR;
      return;
    }

    for (i = 0; i < mesh->num_elems; i++) {
      New_Elem_Index[i] = elements[i].globalID;
      proc_ids[i] = proc;
      change[i] = 0;
    }
  }

  for (i = mesh->num_elems; i < New_Elem_Index_Size; i++) {
    New_Elem_Index[i] = -1;
  }

  for (i = 0; i < num_export; i++) {
    if (num_lid_entries)
      exp_elem = export_local_ids[lid+i*num_lid_entries];
    else  /* testing num_lid_entries == 0 */
      search_by_global_id(mesh, export_global_ids[gid+i*num_gid_entries], 
                          &exp_elem);

    New_Elem_Index[exp_elem] = -1;
    proc_ids[exp_elem] = export_procs[i];
  }

  for (i = 0; i < num_import; i++) {
    /* search for first free location */
    for (j = 0; j < New_Elem_Index_Size; j++) 
      if (New_Elem_Index[j] == -1) break;

    New_Elem_Index[j] = import_global_ids[gid+i*num_gid_entries];
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

    for (j = 0; j < elements[exp_elem].adj_len; j++) {

      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (elements[exp_elem].adj[j] == -1) continue;

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
      if (elements[i].adj[j] == -1) continue;

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
  safe_free((void **) &change);

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

  safe_free((void **) &proc_ids);

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
        if (elements[bor_elem].adj[k] == -1) continue;

        if (elements[bor_elem].adj[k] == mesh->ecmap_neighids[offset] &&
            elements[bor_elem].adj_proc[k] == mesh->ecmap_id[i]) {
          elements[bor_elem].adj_proc[k] = recv_vec[offset];
          if (recv_vec[offset] == proc) {
            /* element is moving to this processor; */
            /* convert adj from global to local ID. */
            if ((idx = in_list(mesh->ecmap_neighids[offset],New_Elem_Index_Size,
                              New_Elem_Index)) == -1) {
              sprintf(msg, "fatal: unable to locate element %d in "
                           "New_Elem_Index", mesh->ecmap_neighids[offset]);
              Gen_Error(0, msg);
              *ierr = ZOLTAN_FATAL;
              return;
            }
            elements[bor_elem].adj[k] = idx;
          }
          break;  /* from k loop */
        }
      }
    }
  }

  safe_free((void **) &recv_vec);
  safe_free((void **) &send_vec);

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
                               int num_export, ZOLTAN_ID_PTR export_global_ids,
                               ZOLTAN_ID_PTR export_local_ids, int *export_procs,
                               int *ierr)
{
MESH_INFO_PTR mesh;
ELEM_INFO *elements;
int proc, num_proc;
int i, j, k, last;
int adj_elem;

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
    if (New_Elem_Index[i] != -1) continue;

    /* Don't want to shift all elements down one position to fill the  */
    /* blank spot -- too much work to adjust adjacencies!  So find the */
    /* last element in the array and move it to the blank spot.        */

    for (last = New_Elem_Index_Size-1; last >= 0; last--)
      if (New_Elem_Index[last] != -1) break;

    /* If (last < i), array is already dense; i is just in some blank spots  */
    /* at the end of the array.  Quit the compacting.                     */
    if (last < i) break;

    /* Copy elements[last] to elements[i]. */
    elements[i] = elements[last];

    /* Adjust adjacencies for local elements.  Off-processor adjacencies */
    /* don't matter here.                                                */

    for (j = 0; j < elements[i].adj_len; j++) {

      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (elements[i].adj[j] == -1) continue;

      adj_elem = elements[i].adj[j];

      /* See whether adjacent element is local; if so, adjust its entry */
      /* for local element i.                                           */
      if (elements[i].adj_proc[j] == proc) {
        for (k = 0; k < elements[adj_elem].adj_len; k++) {
          if (elements[adj_elem].adj[k] == last &&
              elements[adj_elem].adj_proc[k] == proc) {
            /* found adjacency entry for element last; change it to i */
            elements[adj_elem].adj[k] = i;
            break;
          }
        }
      }
    }

    /* Update New_Elem_Index */
    New_Elem_Index[i] = New_Elem_Index[last];
    New_Elem_Index[last] = -1;

    /* clear elements[last] */
    elements[last].globalID = -1;
    elements[last].border = 0;
    elements[last].nadj = 0;
    elements[last].adj_len = 0;
    elements[last].elem_blk = -1;
    for (k=0; k<MAX_CPU_WGTS; k++)
      elements[last].cpu_wgt[k] = 0;
    elements[last].mem_wgt = 0;
    elements[last].coord = NULL;
    elements[last].connect = NULL;
    elements[last].adj = NULL;
    elements[last].adj_proc = NULL;
    elements[last].edge_wgt = NULL;
  }

  if (New_Elem_Index != NULL) safe_free((void **) &New_Elem_Index);
  New_Elem_Index_Size = 0;

  if (!build_elem_comm_maps(proc, mesh)) {
    Gen_Error(0, "Fatal: error rebuilding elem comm maps");
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
   * Compute size of one element's data.
   */

  size = sizeof(ELEM_INFO);
 
  /* Add space to correct alignment so casts work in (un)packing. */
  size = Zoltan_Align(size);

  /* Add space for connect table. */
  if (mesh->num_dims > 0)
    size += num_nodes * sizeof(int);

  /* Add space for adjacency info (elements[].adj and elements[].adj_proc). */
  size += current_elem->adj_len * 2 * sizeof(int);

  /* Assume if one element has edge wgts, all elements have edge wgts. */
  if (Use_Edge_Wgts) {
    /* Add space to correct alignment so casts work in (un)packing. */
    size = Zoltan_Align(size);
    size += current_elem->adj_len * sizeof(float);
  }

  /* Add space for coordinate info */
  size = Zoltan_Align(size);
  size += num_nodes * mesh->num_dims * sizeof(float);
  
  return (size);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void migrate_pack_elem(void *data, int num_gid_entries, int num_lid_entries,
                       ZOLTAN_ID_PTR elem_gid, ZOLTAN_ID_PTR elem_lid,
                       int mig_proc, int elem_data_size, char *buf, int *ierr)
{
  MESH_INFO_PTR mesh;
  ELEM_INFO *elem, *elem_mig;
  ELEM_INFO *current_elem;
  int *buf_int;
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

  buf_int = (int *) (buf + size);

  /* copy the connect table */
  if (mesh->num_dims > 0) {
    for (i = 0; i < num_nodes; i++) {
      *buf_int = current_elem->connect[i];
      buf_int++;
    }
    size += num_nodes * sizeof(int);
  }

  /* copy the adjacency info */
  /* send globalID for all adjacencies */
  for (i =  0; i < current_elem->adj_len; i++) {
    if (current_elem->adj[i] != -1 && current_elem->adj_proc[i] == proc) 
      *buf_int = New_Elem_Index[current_elem->adj[i]];
    else
      *buf_int = current_elem->adj[i];
    buf_int++;
    *buf_int = current_elem->adj_proc[i];
    buf_int++;
  }
  size += current_elem->adj_len * 2 * sizeof(int);

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
   * do so by setting the globalID to -1
   */
  current_elem->globalID = -1;
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
  int *buf_int;
  float *buf_float;
  int size, num_nodes;
  int i, j, idx;
  int proc;
  int gid = num_gid_entries-1;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  mesh = (MESH_INFO_PTR) data;
  elem = mesh->elements;
  elem_mig = (ELEM_INFO *) buf;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  if ((idx = in_list(elem_gid[gid], New_Elem_Index_Size, New_Elem_Index)) == -1) {
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
  buf_int = (int *) (buf + size);

  /* copy the connect table */
  if (mesh->num_dims > 0) {
    current_elem->connect = (int *) malloc(num_nodes * sizeof(int));
    if (current_elem->connect == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      *ierr = ZOLTAN_MEMERR;
      return;
    }
    for (i = 0; i < num_nodes; i++) {
      current_elem->connect[i] = *buf_int;
      buf_int++;
    }
    size += num_nodes * sizeof(int);
  }

  /* copy the adjacency info */
  /* globalIDs are received; convert to local IDs when adj elem is local */
  if (current_elem->adj_len > 0) {
    current_elem->adj      = (int *)malloc(current_elem->adj_len * sizeof(int));
    current_elem->adj_proc = (int *)malloc(current_elem->adj_len * sizeof(int));
    if (current_elem->adj == NULL || current_elem->adj_proc == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      *ierr = ZOLTAN_MEMERR;
      return;
    }
    for (i =  0; i < current_elem->adj_len; i++) {
      current_elem->adj[i] = *buf_int;
      buf_int++;
      current_elem->adj_proc[i] = *buf_int;
      buf_int++;
      if (current_elem->adj[i] != -1 && current_elem->adj_proc[i] == proc) 
        current_elem->adj[i] = in_list(current_elem->adj[i], 
                                       New_Elem_Index_Size, New_Elem_Index);
    }
    size += current_elem->adj_len * 2 * sizeof(int);

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
