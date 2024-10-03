// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "dr_const.h"
#include "dr_externs.h"
#include "dr_maps_const.h"
#include "dr_util_const.h"
#include "dr_err_const.h"
#include "dr_par_util_const.h"

#include <limits.h>

// This is a partial port to C++, only changing the code to use
// the C++ bindings for Zoltan.

#include "zoltan_dd_cpp.h"
#include "zoltan_comm_cpp.h"

#define MAP_ALLOC 10

/*
 *  Routines to build elemental communication maps, given a distributed mesh.
 */

struct map_list_head {
  int map_alloc_size;
  ZOLTAN_ID_TYPE *glob_id;
  int *elem_id;
  int *side_id;
  ZOLTAN_ID_TYPE *neigh_id;
};

static void compare_maps_with_ddirectory_results(int, MESH_INFO_PTR);
static void sort_and_compare_maps(int, int, MESH_INFO_PTR, 
  struct map_list_head *, int, int *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int build_elem_comm_maps(int proc, MESH_INFO_PTR mesh)
{
/*
 * Build element communication maps, given a distributed mesh.
 * This routine builds initial communication maps for Chaco input
 * (for Nemesis, initial communication maps are read from the Nemesis file)
 * and rebuilds communication maps after data migration.
 *
 * One communication map per neighboring processor is built.
 * The corresponding maps on neighboring processors
 * must be sorted in the same order, so that neighboring processors do not
 * have to use ghost elements.   For each communication map's pair of
 * processors, the lower-numbered processor determines the order of the
 * elements in the communication map.  The sort key is the elements' global
 * IDs on the lower-number processor; the secondary key is the neighboring
 * elements global IDs.  The secondary key is used when a single element
 * must communicate with more than one neighbor.
 */

const char *yo = "build_elem_comm_maps";
int i, j;
ELEM_INFO *elem;
ZOLTAN_ID_TYPE iadj_elem;
int iadj_proc;
int indx;
int num_alloc_maps;
int max_adj = 0;
int max_adj_per_map;
int cnt, offset;
int *sindex = NULL;
int tmp;
struct map_list_head *tmp_maps = NULL, *map = NULL;

  DEBUG_TRACE_START(proc, yo);

  /*
   *  Free the old maps, if they exist.
   */

  if (mesh->ecmap_id != NULL) {
    safe_free((void **) &(mesh->ecmap_id));
    safe_free((void **) &(mesh->ecmap_cnt));
    safe_free((void **) &(mesh->ecmap_elemids));
    safe_free((void **) &(mesh->ecmap_sideids));
    safe_free((void **) &(mesh->ecmap_neighids));
    mesh->necmap = 0;
  }

  /*
   *  Look for off-processor adjacencies.
   *  Loop over all elements 
   */

  num_alloc_maps = MAP_ALLOC;
  mesh->ecmap_id = (int *) malloc(num_alloc_maps * sizeof(int));
  mesh->ecmap_cnt = (int *) malloc(num_alloc_maps * sizeof(int));
  tmp_maps = (struct map_list_head*) malloc(num_alloc_maps 
                                          * sizeof(struct map_list_head));

  if (mesh->ecmap_id == NULL || mesh->ecmap_cnt == NULL || tmp_maps == NULL) {
    Gen_Error(0, "Fatal:  insufficient memory");
    DEBUG_TRACE_END(proc, yo);
    return 0;
  }

  for (i = 0; i < mesh->num_elems; i++) {
    elem = &(mesh->elements[i]);
    for (j = 0; j < elem->adj_len; j++) {

      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (elem->adj[j] == ZOLTAN_ID_INVALID) continue;

      iadj_elem = elem->adj[j];
      iadj_proc = elem->adj_proc[j];

      if (iadj_proc != proc) {
        /* 
         * Adjacent element is off-processor.
         * Add this element to the temporary data structure for 
         * the appropriate neighboring processor.
         */
        if ((indx = in_list2(iadj_proc, mesh->necmap, mesh->ecmap_id)) == -1) {
          /*
           * Start a new communication map.
           */

          if (mesh->necmap >= num_alloc_maps) {
            num_alloc_maps += MAP_ALLOC;
            mesh->ecmap_id = (int *) realloc(mesh->ecmap_id,
                                            num_alloc_maps * sizeof(int));
            mesh->ecmap_cnt = (int *) realloc(mesh->ecmap_cnt,
                                             num_alloc_maps * sizeof(int));
            tmp_maps = (struct map_list_head *) realloc(tmp_maps,
                               num_alloc_maps * sizeof(struct map_list_head));
            if (mesh->ecmap_id == NULL || mesh->ecmap_cnt == NULL || 
                tmp_maps == NULL) {
              Gen_Error(0, "Fatal:  insufficient memory");
              DEBUG_TRACE_END(proc, yo);
              return 0;
            }
          }
          mesh->ecmap_id[mesh->necmap] = iadj_proc;
          mesh->ecmap_cnt[mesh->necmap] = 0;
          map = &(tmp_maps[mesh->necmap]);
          map->glob_id  = (ZOLTAN_ID_TYPE *) malloc(MAP_ALLOC * sizeof(ZOLTAN_ID_TYPE));
          map->elem_id  = (int *) malloc(MAP_ALLOC * sizeof(int));
          map->side_id  = (int *) malloc(MAP_ALLOC * sizeof(int));
          map->neigh_id = (ZOLTAN_ID_TYPE *) malloc(MAP_ALLOC * sizeof(ZOLTAN_ID_TYPE));
          if (map->glob_id == NULL || map->elem_id == NULL || 
              map->side_id == NULL || map->neigh_id == NULL) {
            Gen_Error(0, "Fatal:  insufficient memory");
            DEBUG_TRACE_END(proc, yo);
            return 0;
          }
          map->map_alloc_size = MAP_ALLOC;
          indx = mesh->necmap;
          mesh->necmap++;
        }
        /* Add to map for indx. */
        map = &(tmp_maps[indx]);
        if (mesh->ecmap_cnt[indx] >= map->map_alloc_size) {
          map->map_alloc_size += MAP_ALLOC;
          map->glob_id  = (ZOLTAN_ID_TYPE *) realloc(map->glob_id, map->map_alloc_size * sizeof(ZOLTAN_ID_TYPE));
          map->elem_id  = (int *) realloc(map->elem_id, 
                                          map->map_alloc_size * sizeof(int));
          map->side_id  = (int *) realloc(map->side_id, 
                                          map->map_alloc_size * sizeof(int));
          map->neigh_id = (ZOLTAN_ID_TYPE *) realloc(map->neigh_id, map->map_alloc_size * sizeof(ZOLTAN_ID_TYPE));
          if (map->glob_id == NULL || map->elem_id == NULL || 
              map->side_id == NULL || map->neigh_id == NULL) {
            Gen_Error(0, "Fatal:  insufficient memory");
            DEBUG_TRACE_END(proc, yo);
            return 0;
          }
        }
        tmp = mesh->ecmap_cnt[indx];
        map->glob_id[tmp] = elem->globalID;
        map->elem_id[tmp] = i;
        map->side_id[tmp] = j+1;  /* side is determined by position in
                                          adj array (+1 since not 0-based). */
        map->neigh_id[tmp] = iadj_elem;
        mesh->ecmap_cnt[indx]++;
        max_adj++;
      }
    }
  }

  /* 
   * If no communication maps, don't need to do anything else. 
   */

  if (mesh->necmap > 0) {

    /*
     * Allocate data structure for element communication map arrays.
     */

    mesh->ecmap_elemids  = (int *) malloc(max_adj * sizeof(int));
    mesh->ecmap_sideids  = (int *) malloc(max_adj * sizeof(int));
    mesh->ecmap_neighids = (ZOLTAN_ID_TYPE *) malloc(max_adj * sizeof(ZOLTAN_ID_TYPE));


    /*
     * Allocate temporary memory for sort index.
     */
    max_adj_per_map = 0;
    for (i = 0; i < mesh->necmap; i++)
      if (mesh->ecmap_cnt[i] > max_adj_per_map)
        max_adj_per_map = mesh->ecmap_cnt[i];
    sindex = (int *) malloc(max_adj_per_map * sizeof(int));

    cnt = 0;
    for (i = 0; i < mesh->necmap; i++) {

      map = &(tmp_maps[i]);
      for (j = 0; j < mesh->ecmap_cnt[i]; j++)
        sindex[j] = j;

      /*
       * Sort the map so that adjacent processors have the same ordering
       * for the communication.  
       * Assume the ordering of the lower-numbered processor in the pair
       * of communicating processors.
       */

      if (proc < mesh->ecmap_id[i]) 
        quicksort_pointer_inc_id_id(sindex, map->glob_id, map->neigh_id,
                                    0, mesh->ecmap_cnt[i]-1);
      else
        quicksort_pointer_inc_id_id(sindex, map->neigh_id, map->glob_id,
                                    0, mesh->ecmap_cnt[i]-1);

      /*
       * Copy sorted data into elem map arrays. 
       */

      offset = cnt;
      for (j = 0; j < mesh->ecmap_cnt[i]; j++) {
        mesh->ecmap_elemids[offset]  = map->elem_id[sindex[j]];
        mesh->ecmap_sideids[offset]  = map->side_id[sindex[j]];
        mesh->ecmap_neighids[offset] = map->neigh_id[sindex[j]];
        offset++;
      }

      cnt += mesh->ecmap_cnt[i];
    }
  }

  /* Free temporary data structure. */
  for (i = 0; i < mesh->necmap; i++) {
    safe_free((void **) &(tmp_maps[i].glob_id));
    safe_free((void **) &(tmp_maps[i].elem_id));
    safe_free((void **) &(tmp_maps[i].side_id));
    safe_free((void **) &(tmp_maps[i].neigh_id));
  }
  safe_free((void **) &tmp_maps);
  safe_free((void **) &sindex);

  if (Test.DDirectory) 
    compare_maps_with_ddirectory_results(proc, mesh);

  DEBUG_TRACE_END(proc, yo);
  return 1;
}

/******************************************************************************/
/******************************************************************************/

static void compare_maps_with_ddirectory_results(
  int proc, 
  MESH_INFO_PTR mesh
)
{
/*
 * Routine to demonstrate the use of the Zoltan Distributed Directory 
 * to build communication maps.  This functionality essentially duplicates 
 * that in build_elem_comm_maps.  It provides a test of the DDirectory,
 * as the maps generated by the directory should match those generated
 * by build_elem_comm_maps.
 * This test is probably more complicated than necessary, but perhaps it
 * will be adapted to become a communication maps tool within Zoltan.
 */
static const int want_size = 4;
int num_elems = mesh->num_elems;
Zoltan_DD *dd = NULL;
ZOLTAN_ID_PTR gids = NULL;
ZOLTAN_ID_PTR lids = NULL;
ZOLTAN_ID_PTR my_gids = NULL;
ZOLTAN_ID_PTR nbor_gids = NULL;
ZOLTAN_ID_PTR nbor_lids = NULL;
ZOLTAN_ID_PTR i_want = NULL;
ZOLTAN_ID_PTR others_want = NULL;
int *ownerlist = NULL;
int *sindex = NULL;
int num_nbor = 0;    /* Number of neighboring elements not on this processor. */
                     /* This counter counts duplicate entries when an element */
                     /* is a neighbor of > 1 element on this processor.       */
int cnt;
int num_others = 0;
int num_maps = 0;
int map_size = 0;
int max_map_size = 0;
int nbor_proc;
int i, j, k, ierr, error = 0, gerror = 0;
ELEM_INFO_PTR current;
struct map_list_head map;
Zoltan_Comm *comm;


  /* Load array of element globalIDs for elements on this processor. */
  /* Count number of neighboring elements not on this processor.     */

  gids = (ZOLTAN_ID_PTR) malloc(sizeof(ZOLTAN_ID_TYPE) * 2 * num_elems);
  if (num_elems > 0 && gids == NULL) {
    Gen_Error(0, "Fatal:  insufficient memory");
    error = 1;
  }
  lids = gids + num_elems;

  for (i = 0; i < num_elems; i++) {
    current = &(mesh->elements[i]);
    gids[i] = current->globalID;
    lids[i] = i;
    for (j = 0; j < current->adj_len; j++) {
      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (current->adj[j] == ZOLTAN_ID_INVALID) continue;
      if (current->adj_proc[j] != proc)
        num_nbor++;
    }
  }

  /* 
   * Create DDirectory and register all owned elements. 
   */

  dd = new Zoltan_DD(zoltan_get_global_comm(), 1, 1, 0, 0, 0);

  ierr = dd->Update(gids, lids, NULL, NULL, num_elems);

  if (ierr) {
    Gen_Error(0, "Fatal:  Error returned by Zoltan_DD::Update");
    error = 1;
  }

  free(gids);

  /*
   * Test copy operator and constructor
   */

  Zoltan_DD ddCopy(*dd);  // uses copy constructor
  Zoltan_DD ddNew;        // uses ordinary constructor

  delete dd;

  ddNew = ddCopy;         // uses copy operator

  dd = &ddNew;

  /* 
   * Use the DDirectory to find owners of off-processor neighboring elements. 
   * Of course, we have this info in ELEM_INFO, but we do the find to test
   * the DDirectory.
   */

  nbor_gids = (ZOLTAN_ID_PTR) malloc(sizeof(ZOLTAN_ID_TYPE) * 3 * num_nbor);
  if (num_nbor > 0 && nbor_gids == NULL) {
    Gen_Error(0, "Fatal:  insufficient memory");
    error = 1;
  }
  nbor_lids = nbor_gids + num_nbor;
  my_gids = nbor_lids + num_nbor;
  ownerlist = (int *) malloc(sizeof(int) * num_nbor);

  /* 
   * Get list of elements whose info is needed. 
   */
  cnt = 0;
  for (i = 0; i < num_elems; i++) {
    current = &(mesh->elements[i]);
    for (j = 0; j < current->adj_len; j++) {
      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (current->adj[j] == ZOLTAN_ID_INVALID) continue;
      if (current->adj_proc[j] != proc) {
        nbor_gids[cnt] = current->adj[j];
        my_gids[cnt] = current->globalID;
        cnt++;
      }
    }
  }

  /* Sanity check */
  if (cnt != num_nbor) {
    Gen_Error(0, "Fatal:  cnt != num_nbor");
    error = 1;
  }

  ierr = dd->Find( nbor_gids, nbor_lids, NULL, NULL, num_nbor, ownerlist);

  if (ierr) {
    Gen_Error(0, "Fatal:  Error returned by Zoltan_DD::Find");
    error = 1;
  }

  /*
   * Check for errors 
   */

  MPI_Allreduce(&error, &gerror, 1, MPI_INT, MPI_SUM, zoltan_get_global_comm());

  if (gerror) {
    Gen_Error(0, "Fatal:  Error returned by DDirectory Test");
    error_report(proc);
    free(nbor_gids);
    free(ownerlist);
    return;
  }

  /* 
   * Use the Communication library to invert this information and build 
   * communication maps.
   * We know what to receive and from where; compute what to send and
   * to whom.
   * Again, this info was computed by build_elem_com_maps, but we are
   * testing DDirectory here.
   */

  /* 
   * Build list of off-proc elements (and their owners)
   * that this proc wants.
   * This list includes duplicate entries when an element 
   * is a neighbor of > 1 element on this processor.
   */

  i_want = (ZOLTAN_ID_PTR) malloc(sizeof(ZOLTAN_ID_TYPE) * want_size * num_nbor);
  if (num_nbor > 0 && i_want == NULL) {
    Gen_Error(0, "Fatal:  insufficient memory");
    return;
  }

  j = 0;
  for (i = 0; i < num_nbor; i++) {
    i_want[j++] = proc;
    i_want[j++] = nbor_gids[i];
    i_want[j++] = nbor_lids[i];
    i_want[j++] = my_gids[i];
  }

  comm = new Zoltan_Comm(num_nbor, ownerlist, zoltan_get_global_comm(), 747, 
                        &num_others);

  /*
   * Test copy operator and constructor
   */

printf("Test comm copy functions\n");
  Zoltan_Comm commCopy(*comm);    // uses copy constructor
  Zoltan_Comm commNew;            // uses ordinary constructor

  delete comm;

  commNew = commCopy;             // uses copy operator

  comm = &commNew;

  /* 
   * Do communication to determine which of this proc's data is wanted by 
   * other procs.
   * This info will determine what is put in this proc's communication maps.
   */

  others_want = (ZOLTAN_ID_PTR) malloc(sizeof(ZOLTAN_ID_TYPE)*want_size*(num_others+1));
  if (others_want == NULL) {
    Gen_Error(0, "Fatal:  insufficient memory");
    return;
  }

  ierr = comm->Do(757, (char *) i_want, want_size * sizeof(ZOLTAN_ID_TYPE), 
                    (char *) others_want);
  if (ierr) {
    Gen_Error(0, "Fatal:  Error returned from Zoltan_Comm::Do");
    return;
  }

  free(i_want);

  /*  
   * Find number of maps and the size of the largest map. 
   * The maps should be grouped by neighboring processor
   * in others_want.
   */

  num_maps = 0;
  map_size = max_map_size = 0;
  for (i = 0, j = 0; i < num_others; i++, j += want_size) {
    nbor_proc = (int)others_want[j];
    map_size++;
    if (i == (num_others - 1) || nbor_proc != (int) others_want[j+want_size]) {
      /* End of map reached */
      num_maps++;
      if (map_size > max_map_size) max_map_size = map_size;
      map_size = 0;
    }
  }

  if (num_maps != mesh->necmap) {
    printf("%d DDirectory Test:  Different number of maps: "
           "%d != %d\n", proc, num_maps, mesh->necmap);
  }

  /* 
   * For each map, 
   *   build a map_list_head for the map;
   *   sort the map_list_head appropriately (as in build_elem_comm_maps);
   *   compare sorted lists with actually communication maps.
   */

  sindex = (int *) malloc(max_map_size * sizeof(int));
  if (max_map_size > 0 && sindex == NULL) {
    Gen_Error(0, "Fatal:  insufficient memory");
    return;
  }

  map.map_alloc_size = max_map_size;
  map.glob_id  = (ZOLTAN_ID_TYPE *) malloc(max_map_size * sizeof(ZOLTAN_ID_TYPE));
  map.elem_id  = (int *) malloc(max_map_size * sizeof(int));
  map.side_id  = (int *) malloc(max_map_size * sizeof(int));
  map.neigh_id  = (ZOLTAN_ID_TYPE *) malloc(max_map_size * sizeof(ZOLTAN_ID_TYPE));
  
  if (max_map_size > 0 && map.neigh_id == NULL) {
    Gen_Error(0, "Fatal:  insufficient memory");
    return;
  }
  
  if (Debug_Driver > 3) {
    /* For high debug levels, serialize the following section so that
     * output of generated map is serialized (and not junked up).
     */
    print_sync_start(proc, 1);
  }

  map_size = 0;
  for (i = 0, j = 0; i < num_others; i++, j += want_size) {
    nbor_proc = (int)others_want[j];
    map.glob_id[map_size] = others_want[j+1];
    map.elem_id[map_size] = others_want[j+2];
    current = &(mesh->elements[map.elem_id[map_size]]);
    for (k = 0; k < current->adj_len; k++) {
      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (current->adj[k] == ZOLTAN_ID_INVALID) continue;
      if (current->adj_proc[k] == nbor_proc && 
          current->adj[k] == others_want[j+3]) {
        map.side_id[map_size] = k + 1;
        map.neigh_id[map_size] = current->adj[k];
        break;
      }
    }
    map_size++;
    if (i == (num_others-1) || nbor_proc != (int) others_want[j+want_size]) {
      /*
       * End of map has been reached.
       * Sort and compare the current map.
       */
      sort_and_compare_maps(proc, nbor_proc, mesh, &map, map_size, sindex);
      
      /*
       * Reinitialize data structures for new map.
       */
      map_size = 0;
    }
  }

  if (Debug_Driver > 3) {
    /* For high debug levels, serialize the previous section so that
     * output of generated map is serialized (and not junked up).
     */
    int nprocs = 0;
    MPI_Comm_size(zoltan_get_global_comm(), &nprocs);
    print_sync_end(proc, nprocs, 1);
  }

  free(map.glob_id);
  free(map.elem_id);
  free(map.neigh_id);
  free(map.side_id);
  free(sindex);
  free(nbor_gids);
  free(ownerlist);
  free(others_want);
}


static void sort_and_compare_maps(
  int proc,
  int nbor_proc,
  MESH_INFO_PTR mesh, 
  struct map_list_head *map, 
  int map_size, 
  int *sindex
)
{
/*
 *  Routine to sort a given communication map for a single neighbor processor
 *  and compare it to the actual communication map for that processor
 *  (generated by build_elem_comm_maps).
 *  If the DDirectory were used to build comm maps, this routine could be
 *  modified to do the assignments to mesh->ecmap_*, rather than do comparisons
 *  with it.
 */
int i, j;
int cnt = 0;
int indx;

  /*
   *  Sort the given map according to element ids.
   *  Primary key is determined as in build_elem_comm_maps.
   */

  for (i = 0; i < map_size; i++)
    sindex[i] = i;

  if (proc < nbor_proc)
    quicksort_pointer_inc_id_id(sindex, map->glob_id, map->neigh_id,
                                0, map_size-1);
  else 
    quicksort_pointer_inc_id_id(sindex, map->neigh_id, map->glob_id,
                                0, map_size-1);
  
  /*
   * Compute offset into mesh communication maps for the given nbor proc.
   */
  if ((indx = in_list2(nbor_proc, mesh->necmap, mesh->ecmap_id)) == -1) {
    printf("%d DDirectory Test:  Comm map for nbor proc %d does not exist\n",
           proc, nbor_proc);
    return;
  }

  cnt = 0;
  for (i = 0; i < indx; i++)
    cnt += mesh->ecmap_cnt[i];

  /*
   * Compare given map to mesh communication map for this nbor proc.
   * If the DDirectory code were used to construct maps, assignments
   * would be done here (rather than comparisons).
   */

  if (map_size != mesh->ecmap_cnt[indx]) {
    printf("%d DDirectory Test:  Different map size for nbor_proc %d: "
           "%d != %d\n", proc, nbor_proc, map_size, mesh->ecmap_cnt[indx]);
    return;
  }

  for (i = 0; i < map_size; i++) {
    j = sindex[i];
    if (map->elem_id[j] != mesh->ecmap_elemids[i+cnt]) {
      printf("%d DDirectory Test: Different element IDs for nbor_proc %d: "
             "%d != %d\n", proc, nbor_proc, map->elem_id[j], 
             mesh->ecmap_elemids[i+cnt]);
    }
  }

  for (i = 0; i < map_size; i++) {
    j = sindex[i];
    if (map->side_id[j] != mesh->ecmap_sideids[i+cnt]) {
      printf("%d DDirectory Test: Different side IDs for nbor_proc %d: "
             "%d != %d\n", proc, nbor_proc, map->side_id[j], 
             mesh->ecmap_sideids[i+cnt]);
    }
  }

  for (i = 0; i < map_size; i++) {
    j = sindex[i];
    if (map->neigh_id[j] != mesh->ecmap_neighids[i+cnt]) {
      printf("%d DDirectory Test: Different neigh IDs for nbor_proc %d: "
             ZOLTAN_ID_SPEC " != " ZOLTAN_ID_SPEC "\n", proc, nbor_proc, map->neigh_id[j], 
             mesh->ecmap_neighids[i+cnt]);
    }
  }

  if (Debug_Driver > 3) {
    printf("%d  *************   DDirectory Map for %d    ***************\n",
            proc, nbor_proc);
    printf("Local ID\tSide ID\tGlobal ID\tNeigh ID\n");
    for (i = 0; i < map_size; i++) {
      j = sindex[i];
      printf("\t%d\t%d\t" ZOLTAN_ID_SPEC "\t" ZOLTAN_ID_SPEC "\n", 
             map->elem_id[j], map->side_id[j], 
             map->glob_id[j], map->neigh_id[j]);
    }
  }
}

