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
 * Questions? Contact Karen Devine      kddevin@sandia.gov
 *                    Erik Boman        egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "zz_rand.h"
#include "zz_util_const.h"
#include "zz_id_const.h"
#include "params_const.h"
#include "all_allo_const.h"
#include "hier.h"
#include "zoltan_comm.h"
#include "key_params.h"
#include "parmetis_interface_params.h"
#include "scotch_interface_params.h"
#include "third_library_params.h"
#include "graph_params.h"
#include "order_params.h"
#include "phg_params.h"
#include "rcb_params.h"
#include "rib_params.h"
#include "hsfc_params.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#define ZOLTAN_MAX_SIBLINGS  32768      /* A sanity check, valid in 2010 AD */

/***********************************************************************/
/*
 * The topology can be provided by giving a platform name in the
 * PLATFORM_NAME parameter.
 *
 * Alternatively the topology can be indicated in a string with the
 * TOPOLOGY parameter.
 *
 * Quad-socket six-core:  TOPOLOGY="4,6".
 * Quad-socket six-core , every core pair shares a cache: TOPOLOGY="4,3,2"
 * We assume that every node has the same structure.
 *
 * If you add a node topology to the zoltan_hier_platform_specs, include
 * a "1" before the
 * node topology numbers, and include it in the level count.  Process are
 * partitioned
 * across nodes before they are partitioned within nodes.
 *
 * The machine name should be lower case.
 */

/*
 *  Platforms whose hierarchical structure is built-in to Zoltan.
 *  zoltan_platform_names should be in the same order as
 *  zoltan_platform_specification
 */

enum zoltan_platform_names {
      GLORY,
      REDSKY,
      CTX,
      ODIN,
      OCTOPI,
      S861036,
      ZOLTAN_HIER_LAST_PLATFORM
};

static zoltan_platform_specification
       zoltan_hier_platform_specs[ZOLTAN_HIER_LAST_PLATFORM]={

{"glory",         /* machine named Glory */
  3,              /* 3-level hierarchy */
  {1, 4, 4},      /* 1 node, 4 sockets, 4 cpus */
  {-1}},          /* dummy initializer to avoid compiler warnings */

{"redsky",        /* machine named RedSky */
  3,              /* 3-level hierarchy */
  {1, 2, 4},      /* 1 node, 2 sockets, 4 cpus */
  {-1}},          /* dummy initializer to avoid compiler warnings */

{"ctx",           /* machine named CTX */
  3,              /* 3-level hierarchy */
  {1, 2, 6},      /* 1 node, 2 sockets, 6 cpus */
  {-1}},          /* dummy initializer to avoid compiler warnings */

{"odin",          /* machine named odin */
  3,              /* 3-level hierarchy */
  {1, 2, 4},      /* 1 node, 2 sockets, 4 cpus */
  {-1}},          /* dummy initializer to avoid compiler warnings */

{"octopi",        /* eight-core machine named octopi */
  2,              /* 2-level hierarchy */
  {2, 4},         /* 2 sockets, 4 cpus */
  {-1}},          /* dummy initializer to avoid compiler warnings */

{"s861036",      /* dual-core machine named s861036 */
  1,             /* 1-level hierarchy */
  {2},           /* 2 cpus */
  {-1}}          /* dummy initializer to avoid compiler warnings */
};

/***********************************************************************/
static int Zoltan_Hier_Assist_Num_Levels(void *data, int *ierr)
{
  zoltan_platform_specification *spec = (zoltan_platform_specification *)data;
  *ierr = ZOLTAN_OK;

  if (spec == NULL){
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  else{
    return spec->numLevels;
  }
}

/***********************************************************************/
static int Zoltan_Hier_Assist_Part_Number(void *data, int level, int *ierr)
{
  zoltan_platform_specification *spec = (zoltan_platform_specification *)data;
  *ierr = ZOLTAN_OK;

  return spec->my_part[level];
}

/***********************************************************************/
static void Zoltan_Hier_Assist_Method(
  void *data,
  int level,
  struct Zoltan_Struct *to,
  int *ierr)
{
  struct Zoltan_Struct *from;
  *ierr = ZOLTAN_OK;

  from = (struct Zoltan_Struct *)data;

  if (from->LB.Imb_Tol_Len > 0){
    memcpy(to->LB.Imbalance_Tol, from->LB.Imbalance_Tol,
           sizeof(float) * from->LB.Imb_Tol_Len);
  }

  to->Debug_Proc = 0;

  strcpy(to->LB.Approach, from->LB.Approach);

  if (from->Seed > 0){
    to->Seed = from->Seed;
    Zoltan_Srand(from->Seed, NULL);
  }

  if ((from->Get_Num_Edges != NULL || from->Get_Num_Edges_Multi != NULL) &&
      (from->Get_Edge_List != NULL || from->Get_Edge_List_Multi != NULL)) {

    Zoltan_Filter_Params(to, from, Graph_Package_params,
                         from->Debug_Level, to->Proc, 0);
    Zoltan_Filter_Params(to, from, ZG_params, from->Debug_Level,
                         to->Proc, 0);
    Zoltan_Filter_Params(to, from, Order_params, from->Debug_Level,
                         to->Proc, 0);
    Zoltan_Filter_Params(to, from, Parmetis_params, from->Debug_Level,
                         to->Proc, 0);
    Zoltan_Filter_Params(to, from, Scotch_params, from->Debug_Level,
                         to->Proc, 0);
    Zoltan_Filter_Params(to, from, PHG_params, from->Debug_Level, to->Proc, 0);
    Zoltan_Set_Param(to, "LB_METHOD", "GRAPH");
  }
  else if (from->Get_Num_Geom != NULL &&
          (from->Get_Geom != NULL || from->Get_Geom_Multi != NULL)) {

    Zoltan_Filter_Params(to, from, RCB_params, from->Debug_Level, to->Proc, 0);
    Zoltan_Filter_Params(to, from, RIB_params, from->Debug_Level, to->Proc, 0);
    Zoltan_Filter_Params(to, from, HSFC_params, from->Debug_Level, to->Proc, 0);
    Zoltan_Set_Param(to, "LB_METHOD", "RCB");
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }

  return;
}

/*****************************************************************************/
/* parameters for the hierarchical partitioning.  */

static PARAM_VARS Hier_params[] = {
  {  "HIER_DEBUG_LEVEL", NULL, "INT", 0},
  {  "HIER_GENERATE_FILES", NULL, "INT", 0},
  {  "HIER_CHECKS", NULL, "INT" , 0},

  {  "HIER_ASSIST", NULL, "INT", 0}, /* If "1", Zoltan determines hierarchy */
     /* If HIER_ASSIST is "1", define either PLATFORM_NAME or TOPOLOGY */

  {  "PLATFORM_NAME", NULL, "STRING", 0},
     /* a name from zoltan_hier_platform_specs above */

  {  "TOPOLOGY", NULL, "STRING", 0},
     /* topology: for example
          double socket, quad core: "2, 4"
          dual processor work station: "2"
          quad socket, each with 2 L3 caches, 3 cores per cache: "4,2,3"
      */

  { "USE_TIMERS", NULL, "INT", 0 },
  {  NULL,              NULL,  NULL, 0 }};


int Zoltan_Hier_Set_Param(
  char *name,                 /* name of variable */
  char *val                   /* value of variable */
)
{
  int status;
  PARAM_UTYPE result;        /* value returned from Zoltan_Check_Param */
  int index;                 /* index returned from Zoltan_Check_Param */

  status = Zoltan_Check_Param(name, val, Hier_params, &result, &index);
  return(status);
}

/*****************************************************************************/
/* static helper function for Zoltan_Hier */
/* compute part sizes for the current level and set them in hierzz */
/* part_sizes: input array of size
   hpp->origzz->Num_Global_Parts * hpp->origzz->Obj_Weight_Dim
   containing the percentage of work to be
   assigned to each final global part.               */
/* returns error condition */

static int set_hier_part_sizes(HierPartParams *hpp, float *part_sizes) {
  int ierr = ZOLTAN_OK;
  float *my_level_part_sizes=NULL, *level_part_sizes=NULL;
  int *part_ids=NULL, *wgt_idx=NULL;
  int i;
  char msg[256];
  int part_weight_dim = hpp->origzz->Obj_Weight_Dim;

  /* when this is called, hpp->num_parts contains the number of
     parts to be computed at this level, and hpp->part_to_compute
     contains the part id to be computed by this process at this
     level, hpp->hier_comm is a communicator for all procs participating
     at this level. */

  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] set_hier_part_sizes at level %d, computing %d parts\n",
           hpp->origzz->Proc, hpp->level, hpp->num_parts);
  }

  /* careful of part_weight_dim of 0 for variable part sizes */
  if (part_weight_dim == 0) part_weight_dim = 1;

  /* allocate an array for input to reduction to compute
     part sizes for this level */
  my_level_part_sizes = (float *)ZOLTAN_MALLOC(hpp->num_parts *
                                               part_weight_dim *
                                               sizeof(float));
  if (!my_level_part_sizes) {
    sprintf(msg, "Out of memory, tried to alloc %u bytes",
        (unsigned int)(hpp->num_parts * part_weight_dim * sizeof(float)));
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, "set_hier_part_sizes", msg);
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  for (i=0; i<hpp->num_parts * part_weight_dim; i++) {
    my_level_part_sizes[i] = 0;
  }

  /* put in my part_sizes for the part I'll be computing */
  for (i=0; i<part_weight_dim; i++) {
    my_level_part_sizes[hpp->part_to_compute*part_weight_dim+i] =
      part_sizes[hpp->origzz->Proc*part_weight_dim+i];
  }

  /* allocate an array for result of reduction of
     part sizes for this level */
  level_part_sizes = (float *)ZOLTAN_MALLOC(hpp->num_parts *
                                            part_weight_dim *
                                            sizeof(float));
  if (!level_part_sizes) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, "set_hier_part_sizes",
                       "Out of memory");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* do the reduction to get global array on each proc */
  MPI_Allreduce(my_level_part_sizes, level_part_sizes,
                hpp->num_parts * part_weight_dim,
                MPI_FLOAT, MPI_SUM, hpp->hier_comm);

  /* allocate and populate extra args to set_part_sizes) */
  part_ids = (int *)ZOLTAN_MALLOC(hpp->num_parts *
                                  part_weight_dim *
                                  sizeof(int));
  wgt_idx = (int *)ZOLTAN_MALLOC(hpp->num_parts *
                                 part_weight_dim *
                                 sizeof(int));
  if (!part_ids || !wgt_idx) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, "set_hier_part_sizes",
                       "Out of memory");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  for (i=0; i<hpp->num_parts * part_weight_dim; i++) {
    part_ids[i] = i/part_weight_dim;
    wgt_idx[i] = i%part_weight_dim;
  }

  if (hpp->output_level >= HIER_DEBUG_ALL) {
    for (i=0; i<hpp->num_parts * part_weight_dim; i++) {
      printf("[%d] setting part_size[%d] to %.3f\n", hpp->origzz->Proc,
             i, level_part_sizes[i]);
    }
  }

  /* set the part sizes in hpp->hierzz */
  Zoltan_LB_Set_Part_Sizes(hpp->hierzz, 1,
                           hpp->num_parts * part_weight_dim,
                           part_ids, wgt_idx, level_part_sizes);

End:
  ZOLTAN_FREE(&my_level_part_sizes);
  ZOLTAN_FREE(&level_part_sizes);
  ZOLTAN_FREE(&part_ids);
  ZOLTAN_FREE(&wgt_idx);

  return ierr;
}

/*****************************************************************************/
/*****************************************************************************/
/** After final partitioning, migrate only the objects to the new owners,    */
/** so they can report them as the objects to be imported to them.           */
/*****************************************************************************/
/*****************************************************************************/

static int final_migrate(
  HierPartParams *hpp,
  int num_export,
  ZOLTAN_ID_TYPE *hier_export_gids,
  ZOLTAN_ID_TYPE *hier_export_lids,
  int *hier_export_procs)
{
  ZOLTAN_COMM_OBJ *plan=NULL;
  ZOLTAN_ID_TYPE *importList=NULL;
  int i, nImports, numGno, next, tag, ierr;
  MPI_Comm comm = hpp->hier_comm;

  tag = 11111;

  ierr = Zoltan_Comm_Create(&plan, num_export, hier_export_procs,
                            comm, tag, &nImports);

  if (ierr != ZOLTAN_OK)
    goto End;

  if (nImports > 0){
    importList = (ZOLTAN_ID_TYPE *)
                  ZOLTAN_MALLOC(nImports*sizeof(ZOLTAN_ID_TYPE));
    if (!importList){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  tag++;
  ierr = Zoltan_Comm_Do(plan, tag, (char *)hier_export_gids,
                        sizeof(ZOLTAN_ID_TYPE), (char *)importList);

  if (ierr != ZOLTAN_OK)
    goto End;

  Zoltan_Comm_Destroy(&plan);

  /* Modify gno list to contain this process' gnos.
   * We don't update the related fields
   *  (vwgt, adjncy, ewgt, adjproc) this time because we are done partitioning.
   */

  numGno = hpp->num_obj + nImports - num_export;
  next=0;

  if (numGno){
    ZOLTAN_GNO_TYPE *gnoList = (ZOLTAN_GNO_TYPE *)
                                ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * numGno);
    if (!gnoList){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i=0; i < num_export; i++){
      hpp->gno[hier_export_lids[i]] = hpp->invalid_gno;
    }

    for (i=0; i < hpp->num_obj; i++){
      if (hpp->gno[i] != hpp->invalid_gno)
        gnoList[next++] = hpp->gno[i];
    }

    for (i=0; i < nImports; i++){
      gnoList[next++] = importList[i];
    }

    ZOLTAN_FREE(&hpp->gno);
    hpp->gno = gnoList;
  }
  else{
    ZOLTAN_FREE(&hpp->gno);
  }

  hpp->num_obj = numGno;

End:
  if (plan) Zoltan_Comm_Destroy(&plan);
  ZOLTAN_FREE(&importList);

  return ierr;
}

/*****************************************************************************/
static int migrate_to_next_subgroups(HierPartParams *hpp, int num_export,
  ZOLTAN_ID_TYPE *hier_export_lids, int *hier_export_procs, MPI_Comm next_comm)
{
  MPI_Comm comm;
  int ids[2];
  int i, j, k, w, ierr;
  int rank, size;
  int nextRank, nextSize, nextGroup;
  int gnos_per_gid, keySize, hashTableSize;
  int nVtx, nEdge;
  int numUniqueNbors;
  int start, edim, vdim, gdim;
  int tag, nNewVtx, nNewAdj;
  int ownerNewGroup, adjNewGroup, adjNewRank;
  int *id_map = NULL, *to_proc=NULL;
  int *adjProcNext=NULL, *owner=NULL, *adjOwner=NULL;
  int *newXadj=NULL, *edgeSizes=NULL, *newAdjProc=NULL;
  int *procList=NULL;
  Zoltan_DD_Directory *dd=NULL;
  ZOLTAN_MAP *nborMap=NULL;
  ZOLTAN_COMM_OBJ *plan=NULL;
  ZZ *zz=NULL;
  char *keyptr=NULL;
  ZOLTAN_GNO_TYPE nbor;
  ZOLTAN_GNO_TYPE *adjNext=NULL, *nborList=NULL;
  ZOLTAN_GNO_TYPE *newVtx=NULL, *newAdjncy=NULL;
  float *ewgtNext=NULL, *ewgt=NULL;
  float *newEwgts=NULL, *newVwgts=NULL;
  double *newGeom=NULL;
  intptr_t value_in, value_out;

  ierr = ZOLTAN_OK;

  zz = hpp->hierzz;
  comm = zz->Communicator;
  vdim = hpp->obj_wgt_dim;
  gdim = hpp->ndims;
  nVtx = hpp->num_obj;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  nextGroup = hpp->part_to_compute;
  MPI_Comm_rank(next_comm, &nextRank);
  MPI_Comm_size(next_comm, &nextSize);

  if ( !(id_map = (int *)ZOLTAN_MALLOC(sizeof(int) * 2 * size))){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  ids[0] = nextGroup;
  ids[1] = nextRank;

  MPI_Allgather(ids, 2, MPI_INT, id_map, 2, MPI_INT, comm);

  /* Global mapping of gno's to their new owner in the current group. */

  gnos_per_gid = sizeof(ZOLTAN_GNO_TYPE) / sizeof(ZOLTAN_ID_TYPE);
  ierr = Zoltan_DD_Create(&dd, comm, gnos_per_gid, 0, 0, nVtx, 0);

  if (ierr != ZOLTAN_OK)
    goto End;

  to_proc = (int *)ZOLTAN_MALLOC(sizeof(int) * nVtx);

  if (nVtx && !to_proc){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  for (i=0; i < nVtx; i++){
    to_proc[i] = rank;
  }

  for (i=0; i < num_export; i++){
    to_proc[hier_export_lids[i]] = hier_export_procs[i];
  }

  ierr = Zoltan_DD_Update(dd, (ZOLTAN_ID_PTR)hpp->gno, NULL, NULL, to_proc,
                          nVtx);

  if (ierr != ZOLTAN_OK)
    goto End;

  /* Export objects and weights */

  tag = 11111;

  ierr = Zoltan_Comm_Create(&plan, nVtx, to_proc, comm, tag, &nNewVtx);

  if (ierr != ZOLTAN_OK)
    goto End;

  if (nNewVtx > 0){
    ierr = ZOLTAN_MEMERR;

    newVtx = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(nNewVtx*sizeof(ZOLTAN_GNO_TYPE));
    if (!newVtx)
      goto End;

    if (vdim > 0){
      newVwgts = (float *)ZOLTAN_MALLOC(nNewVtx * sizeof(float) * vdim);
      if (!newVwgts)
        goto End;
    }

    if (gdim){
      newGeom= (double *)ZOLTAN_MALLOC(nNewVtx * sizeof(double) * gdim);
      if (!newGeom)
        goto End;
    }

    ierr = ZOLTAN_OK;
  }

  tag++;
  ierr = Zoltan_Comm_Do(plan, tag, (char *)hpp->gno, sizeof(ZOLTAN_GNO_TYPE),
                        (char *)newVtx);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&hpp->gno);
  hpp->gno = newVtx;

  if (vdim > 0){
    tag++;
    ierr = Zoltan_Comm_Do(plan, tag, (char *)hpp->vwgt, sizeof(float) * vdim,
                          (char *)newVwgts);

    if (ierr != ZOLTAN_OK)
      goto End;

    ZOLTAN_FREE(&hpp->vwgt);
    hpp->vwgt= newVwgts;
  }

  if (gdim){
    tag++;
    ierr = Zoltan_Comm_Do(plan, tag, (char *)hpp->geom_vec, sizeof(double)*gdim,
                         (char *)newGeom);

    if (ierr != ZOLTAN_OK)
      goto End;

    ZOLTAN_FREE(&hpp->geom_vec);
    hpp->geom_vec= newGeom;
  }

  hpp->num_obj = nNewVtx;

  if (!hpp->use_graph){
    Zoltan_Comm_Destroy(&plan);
    Zoltan_DD_Destroy(&dd);
    ZOLTAN_FREE(&id_map);
    ZOLTAN_FREE(&to_proc);
    return ierr;
  }

  /* Export graph. Include only adjacencies that will be owned by a process
   * in the new process' next sub group.  For adjproc field, use the rank
   * of the owner in the new sub group.
   */

  edim = hpp->edge_wgt_dim;
  nEdge = hpp->xadj[nVtx];

  /* A map of neighbor vertices to their new owners in current group. */

  keySize = sizeof(ZOLTAN_GNO_TYPE);
  hashTableSize = Zoltan_Recommended_Hash_Size(nVtx * 1.1);

  nborMap = Zoltan_Map_Create(zz,
            hashTableSize,
            keySize,
            1,  /* save value of the key, not a pointer to it */
            0); /* don't know number of keys, so allocate dynamically */

  if (!nborMap){
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  value_in = 1;
  for (i=0; i < nEdge; i++){  /* "Find_Add" means if not found, add it */
    ierr = Zoltan_Map_Find_Add(zz, nborMap,
                              (char *)(hpp->adjncy + i),  /* pointer to key */
                              value_in,                   /* value */
                              &value_out);

    if (ierr != ZOLTAN_OK)
      goto End;

    if (value_out == value_in){   /* That id was not found in the map */
      value_in++;
    }
  }

  numUniqueNbors = value_in - 1;

  if (numUniqueNbors > 0){
    nborList = (ZOLTAN_GNO_TYPE *)
                ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * numUniqueNbors);
    procList = (int *)ZOLTAN_MALLOC(sizeof(int) * numUniqueNbors);

    if (!nborList || !procList){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    ierr = Zoltan_Map_First(zz, nborMap, &keyptr, &value_out);

    if (ierr != ZOLTAN_OK)
      goto End;

    nborList[value_out-1] = *((ZOLTAN_GNO_TYPE *)keyptr);

    for (i=1; i < numUniqueNbors; i++){
      ierr = Zoltan_Map_Next(zz, nborMap, &keyptr, &value_out);

      if (ierr != ZOLTAN_OK)
        goto End;

      nborList[value_out-1] = *((ZOLTAN_GNO_TYPE *)keyptr);
    }
  }

  ierr = Zoltan_DD_Find(dd, (ZOLTAN_ID_PTR)nborList, NULL, NULL, procList,
                        numUniqueNbors, NULL);

  if (ierr != ZOLTAN_OK)
    goto End;

  Zoltan_DD_Destroy(&dd);
  ZOLTAN_FREE(&nborList);

  /* Overwrite adjacency lists with new values. */

  start = 0;
  adjNext = hpp->adjncy;
  ewgtNext = hpp->ewgts;
  adjProcNext = hpp->adjproc;
  ewgt = hpp->ewgts;

  edgeSizes = (int *)ZOLTAN_MALLOC(nVtx * sizeof(int));

  if (nVtx && !edgeSizes){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  for (i=0; i < nVtx; i++){
    owner = id_map + to_proc[i]*2;
    ownerNewGroup = *owner;
    edgeSizes[i] = 0;

    for (j=start; j < hpp->xadj[i+1]; j++){
      nbor = hpp->adjncy[j];
      ierr = Zoltan_Map_Find(zz, nborMap, (char *)(&nbor), &value_out);

      if (value_out == ZOLTAN_NOT_FOUND){
        ierr = ZOLTAN_FATAL;
        goto End;
      }

      adjOwner = id_map + procList[value_out-1] * 2;
      adjNewGroup = *adjOwner++;
      adjNewRank= *adjOwner;

      if (adjNewGroup == ownerNewGroup){
        *adjNext++ = nbor;
        *adjProcNext++ = adjNewRank;
        for (w=0; w < edim; w++){
          *ewgtNext++ = ewgt[w];
        }
        edgeSizes[i]++;
      }
      ewgt += edim;
    }

    start = hpp->xadj[i+1];
    hpp->xadj[i+1] = hpp->xadj[i] + edgeSizes[i];
  }

  ZOLTAN_FREE(&procList);
  ZOLTAN_FREE(&hpp->xadj);
  ZOLTAN_FREE(&id_map);
  Zoltan_Map_Destroy(zz, &nborMap);


  newXadj = (int *)ZOLTAN_MALLOC((nNewVtx+1) * sizeof(int));
  if (!newXadj){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  tag++;
  ierr = Zoltan_Comm_Do(plan, tag, (char *)edgeSizes, sizeof(int),
                        (char *)(newXadj + 1));

  if (ierr != ZOLTAN_OK)
    goto End;

  hpp->xadj = newXadj;

  newXadj[0] = 0;
  for (i=1; i <= nNewVtx; i++){
    newXadj[i] += newXadj[i-1];
  }

  nNewAdj = newXadj[nNewVtx];

  tag++;
  ierr = Zoltan_Comm_Resize(plan, edgeSizes, tag, &k);

  if (ierr != ZOLTAN_OK)
    goto End;

  if (k != nNewAdj){
    /* error */
  }

  ZOLTAN_FREE(&edgeSizes);

  if (nNewAdj > 0){
    ierr = ZOLTAN_MEMERR;

    newAdjncy= (ZOLTAN_GNO_TYPE *)
                ZOLTAN_MALLOC(nNewAdj * sizeof(ZOLTAN_GNO_TYPE));
    if (!newAdjncy)
      goto End;

    newAdjProc= (int *)ZOLTAN_MALLOC(nNewAdj* sizeof(int));
    if (!newAdjProc)
      goto End;

    if (edim){
      newEwgts = (float *)ZOLTAN_MALLOC(sizeof(float) * edim * nNewAdj);
      if (!newEwgts)
        goto End;
    }
    ierr = ZOLTAN_OK;
  }

  tag++;
  ierr = Zoltan_Comm_Do(plan, tag, (char *)hpp->adjncy, sizeof(ZOLTAN_GNO_TYPE),
                        (char *)newAdjncy);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&hpp->adjncy);
  hpp->adjncy = newAdjncy;

  tag++;
  ierr = Zoltan_Comm_Do(plan, tag, (char *)hpp->adjproc, sizeof(int),
                        (char *)newAdjProc);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&hpp->adjproc);
  hpp->adjproc= newAdjProc;

  if (edim){
    tag++;
    ierr = Zoltan_Comm_Do(plan, tag, (char *)hpp->ewgts, sizeof(float) * edim ,
                          (char *)newEwgts);

    if (ierr != ZOLTAN_OK)
      goto End;

    Zoltan_Comm_Destroy(&plan);

    ZOLTAN_FREE(&hpp->ewgts);
    hpp->ewgts= newEwgts;
  }


End:

  if (ierr != ZOLTAN_OK){
    ZOLTAN_FREE(&newVtx);
    ZOLTAN_FREE(&newVwgts);
    ZOLTAN_FREE(&newGeom);
    ZOLTAN_FREE(&newXadj);
    ZOLTAN_FREE(&newAdjncy);
    ZOLTAN_FREE(&newAdjProc);
    ZOLTAN_FREE(&newEwgts);
  }

  ZOLTAN_FREE(&nborList);
  ZOLTAN_FREE(&procList);
  ZOLTAN_FREE(&edgeSizes);
  ZOLTAN_FREE(&id_map);
  ZOLTAN_FREE(&to_proc);

  if (plan) Zoltan_Comm_Destroy(&plan);
  if (dd) Zoltan_DD_Destroy(&dd);

  return ierr;
}

/***********************************************************************/
static void view_hierarchy_specification(
  zoltan_platform_specification *spec,
  int rank,
  int verbose
)
{
int i;

  if (verbose){
    if (spec->platform_name){
      printf("%s\n",spec->platform_name);
    }
    printf("Number of siblings at each level: ");
    for (i=0; i < spec->numLevels; i++){
      printf("%d ",spec->num_siblings[i]);
    }
    printf("\n");
  }

  printf("Part for MPI rank %d at each level: ", rank);
  for (i=0; i < spec->numLevels; i++){
    printf("%d ",spec->my_part[i]);
  }
  printf("\n");

  fflush(stdout);
}

/***********************************************************************/
static char *make_platform_name_string()
{
int i;
int len;
char *msg;
char *yo = "make_platform_name_string";


  for (i=0, len=0; i < ZOLTAN_HIER_LAST_PLATFORM; i++){
    len += strlen(zoltan_hier_platform_specs[i].platform_name);
  }

  len += ((ZOLTAN_HIER_LAST_PLATFORM * 3) + 64);

  msg = (char *)ZOLTAN_MALLOC(len);
  if (!msg){
    ZOLTAN_PRINT_ERROR(-1, yo, "Out of memory");
    return NULL;
  }
  msg[0] = '\0';

  for (i=0; i < ZOLTAN_HIER_LAST_PLATFORM; i++){
    strcat(msg, zoltan_hier_platform_specs[i].platform_name);
    strcat(msg, " ");
    if (i % 5  == 0) strcat(msg, "\n");
  }

  if (ZOLTAN_HIER_LAST_PLATFORM % 5)
    strcat(msg, "\n");

  return msg;
}

/*****************************************************************************/
/* Initialize the parameter structure for hierarchical */
static int Zoltan_Hier_Initialize_Params(ZZ *zz, HierPartParams *hpp) {

  char *yo = "Zoltan_Hier_Initialize_Params";
  int assist, i=0, j, len;
  int num_cpus, num_siblings;
  char platform[MAX_PARAM_STRING_LEN];
  char topology[MAX_PARAM_STRING_LEN];
  char *c=NULL;
  div_t result;

  Zoltan_Bind_Param(Hier_params, "HIER_DEBUG_LEVEL",
                    (void *) &hpp->output_level);
  Zoltan_Bind_Param(Hier_params, "HIER_GENERATE_FILES",
                    (void *) &hpp->gen_files);
  Zoltan_Bind_Param(Hier_params, "HIER_CHECKS", (void *) &hpp->checks);
  Zoltan_Bind_Param(Hier_params, "HIER_ASSIST", (void *) &assist);
  Zoltan_Bind_Param(Hier_params, "PLATFORM_NAME", (void *) platform);
  Zoltan_Bind_Param(Hier_params, "TOPOLOGY", (void *) topology);
  Zoltan_Bind_Param(Hier_params, "USE_TIMERS", (void *) &hpp->use_timers);

  /* set default values */
  hpp->output_level = HIER_DEBUG_NONE;
  hpp->checks = 0;
  assist = 0;
  platform[0] = topology[0] = 0;

  /* Get application values of parameters. */
  Zoltan_Assign_Param_Vals(zz->Params, Hier_params, zz->Debug_Level, zz->Proc,
                           zz->Debug_Proc);

  if (!assist)
    return ZOLTAN_OK;

  if (platform[0]){
    len = strlen(platform);
    for (i=0; i < len; i++){
      if (isupper((int)platform[i]))
        platform[i] = (char)tolower((int)platform[i]);
    }

    for (i=0; i < ZOLTAN_HIER_LAST_PLATFORM; i++){
      if (strcmp(platform, zoltan_hier_platform_specs[i].platform_name))
        continue;
      hpp->spec = zoltan_hier_platform_specs + i;
      break;
    }
  }

  if (!hpp->spec && topology[0]){
    hpp->spec = (zoltan_platform_specification *)
                 ZOLTAN_CALLOC(sizeof(zoltan_platform_specification), 1);

    if (!hpp->spec){
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
      return ZOLTAN_MEMERR;
    }

    hpp->spec->platform_name = NULL;

    if (topology[0]){
      hpp->spec->num_siblings[0] = 1; /* the node or machine itself
                                         is the first level */
      i = 1;
      c = topology;
    }

    j = 0;

    while (*c){
      while (*c && !isdigit(*c)) c++;

      if (*c){
        if (i == ZOLTAN_PLATFORM_MAX_LEVELS){
          ZOLTAN_FREE(&(hpp->spec));
          break;
        }

        sscanf(c, "%d", hpp->spec->num_siblings +  i);

        if ((hpp->spec->num_siblings[i] < 1) ||
            (hpp->spec->num_siblings[i] > ZOLTAN_MAX_SIBLINGS)){
          ZOLTAN_FREE(&(hpp->spec));
          break;
        }
        i++;
        j++;
      }

      while (*c && isdigit(*c)) c++;

    }

    hpp->spec->numLevels = i;

    if (j == 0){
      ZOLTAN_FREE(&(hpp->spec));
    }
  }

  if (!hpp->spec){
    if (zz->Proc == 0){
      char *pnames = make_platform_name_string();
      char *msg = NULL;
      if (pnames == NULL){
        ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
        return ZOLTAN_MEMERR;
      }
      i = strlen(pnames) + 2048;
      msg = (char *)ZOLTAN_MALLOC(i);
      if (!msg){
        ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
        return ZOLTAN_MEMERR;
      }
      strcpy(msg,"Error:\n");
      strcat(msg, "HIER_ASSIST requested but insufficient topology "
                  "information provided.\n\n"
                  "Specify PLATFORM_NAME or TOPOLOGY.\n\n");

      strcat(msg,"TOPOLOGY is the number of hardware siblings "
                 "at each level in a topology.\n"
                 "  Ex. TOPOLOGY=\"2, 4\" describes a dual-socket "
                   "quad-core computing cluster.\n"
                 "  Ex. TOPOLOGY=\"4\" describes a quad-core "
                   "desktop computer.\n\n");

      strcat(msg,"Zoltan assumes the run-time system has "
                 "pinned each process to a CPU.\n");
      strcat(msg,"It assumes MPI process ranks map to the topology.  "
                 "(In the 2,4 example,\n");
      strcat(msg,"this means ranks 0-7 are on the same node, "
                 "and 0-3 on the same socket.)\n\n");

      strcat(msg, "PLATFORM_NAME can be one of the following:\n");
      strcat(msg, pnames);

      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);

      ZOLTAN_FREE(&pnames);
      ZOLTAN_FREE(&msg);
    }

    return ZOLTAN_FATAL;
  }

  /*
   * Compute which part my process has at each level.  We are assuming that
   * MPI laid out the process ranks with respect to the topology.  This
   * may not be true and eventually we want a way to determine the
   * topological rank of each process.
   */

  num_cpus = 1;
  for (i = 0; i < hpp->spec->numLevels; i++)
    num_cpus *= hpp->spec->num_siblings[i];

  result = div(zz->Num_Proc, num_cpus);

  hpp->spec->num_siblings[0] = result.quot;

  if (result.rem > 0)
    hpp->spec->num_siblings[0]++;  /* number of nodes */

  for (i=0; i < hpp->spec->numLevels; i++){

    /* total number of objects at this level */
    num_siblings = hpp->spec->num_siblings[i];

    /* total number of cpus within an object at this level */
    num_cpus = 1;

    for (j = hpp->spec->numLevels-1; j > i; j--)
      num_cpus *= hpp->spec->num_siblings[j];

    result = div(zz->Proc, num_cpus);

    result = div(result.quot, num_siblings);

    hpp->spec->my_part[i] = result.rem;
  }

  if (hpp->output_level >= HIER_DEBUG_LIST){
    MPI_Barrier(zoltan_get_global_comm());
    for (i=0; i < zz->Num_Proc; i++){
      if (i == zz->Proc){
        view_hierarchy_specification(hpp->spec, i, (i==0));
      }
      MPI_Barrier(zoltan_get_global_comm());
    }
    MPI_Barrier(zoltan_get_global_comm());
  }

  return ZOLTAN_OK;
}


/***********************************************************************/
/* callbacks registered by hier method, accessing hpp data structure   */
/***********************************************************************/

static int Zoltan_Hier_Num_Obj_Fn(void *data, int *ierr) {
  HierPartParams *hpp = (HierPartParams *)data;

  *ierr = ZOLTAN_OK;
  return hpp->num_obj;
}

/***********************************************************************/
static void Zoltan_Hier_Obj_List_Fn(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  ZOLTAN_ID_TYPE * global_ids,
  ZOLTAN_ID_TYPE  *local_ids,
  int wgt_dim,
  float *obj_wgts,
  int *ierr)
{
  HierPartParams *hpp = (HierPartParams *)data;
  int j;

  *ierr = ZOLTAN_OK;

  for (j=0; j < hpp->num_obj; j++){
    global_ids[j] = (ZOLTAN_ID_TYPE)hpp->gno[j];
    local_ids[j] = j;
  }

  if (wgt_dim > 0){
    memcpy(obj_wgts, hpp->vwgt, sizeof(float) * hpp->num_obj * wgt_dim);
  }

}

/***********************************************************************/
static int Zoltan_Hier_Num_Geom_Fn(void *data, int *ierr) {
  HierPartParams *hpp = (HierPartParams *)data;

  *ierr = ZOLTAN_OK;
  return hpp->ndims;
}

/***********************************************************************/
static void Zoltan_Hier_Geom_Multi_Fn(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_obj,
  ZOLTAN_ID_TYPE *global_id,
  ZOLTAN_ID_TYPE *local_id,
  int num_dim,
  double *coord,
  int *ierr)
{
  HierPartParams *hpp = (HierPartParams *)data;
  double *coord_ptr;
  int i, j, idx;

  *ierr = ZOLTAN_OK;

  for (i=0; i < num_obj; i++){
    idx = local_id[i];
    coord_ptr = hpp->geom_vec + (idx * num_dim);
    for (j=0; j < num_dim; j++){
      *coord++ = coord_ptr[j];
    }
  }
}

/***********************************************************************/
static void Zoltan_Hier_Num_Edges_Multi_Fn(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_obj,
  ZOLTAN_ID_TYPE *global_id,
  ZOLTAN_ID_TYPE *local_id,
  int *num_edges,
  int *ierr
)
{
  HierPartParams *hpp = (HierPartParams *)data;
  int i, idx;

  *ierr = ZOLTAN_OK;

  for (i=0; i < num_obj; i++){
    idx = local_id[i];
    num_edges[i] = hpp->xadj[idx+1] - hpp->xadj[idx];
  }
}

/***********************************************************************/
static void Zoltan_Hier_Edge_List_Multi_Fn(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_obj,
  ZOLTAN_ID_TYPE *global_id,
  ZOLTAN_ID_TYPE *local_id,
  int *num_edges,
  ZOLTAN_ID_TYPE *nbor_global_id,
  int *nbor_procs,
  int wgt_dim,
  float *ewgts,
  int *ierr)
{
  HierPartParams *hpp = (HierPartParams *)data;
  int i, j, k, idx;
  int *out_proc;
  ZOLTAN_ID_TYPE *out_gid;
  float *out_weight, *wgts;

  *ierr = ZOLTAN_OK;

  out_proc = nbor_procs;
  out_gid = nbor_global_id;
  out_weight = ewgts;

  for (i=0; i < num_obj; i++){
    idx = local_id[i];

    for (j= hpp->xadj[idx]; j < hpp->xadj[idx+1]; j++){
      *out_proc++ = hpp->adjproc[j];
      *out_gid++ = hpp->adjncy[j];
      wgts = hpp->ewgts + j*wgt_dim;
      for (k=0; k < wgt_dim; k++){
        *out_weight++ = *wgts++;
      }
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/** Zoltan_Hier: main routine for hierarchical balancing *********************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Hier(
  ZZ *zz,               /* Zoltan structure */
  float *part_sizes,    /* Input:  Array of size
                           zz->LB.Num_Global_Parts * zz->Obj_Weight_Dim
                           containing the percentage of work to be
                           assigned to each part.               */
  int *num_imp,         /* number of objects to be imported */
  ZOLTAN_ID_TYPE **imp_gids,  /* global ids of objects to be imported */
  ZOLTAN_ID_TYPE **imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,      /* list of processors to import from */
  int **imp_to_part,    /* list of parts to which imported objects are
                           assigned.  */
  int *num_exp,         /* number of objects to be exported */
  ZOLTAN_ID_TYPE **exp_gids,  /* global ids of objects to be exported */
  ZOLTAN_ID_TYPE **exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs,      /* list of processors to export to */
  int **exp_to_part     /* list of parts to which exported objects are
                           assigned. */
)
{
  int ierr = ZOLTAN_OK;   /* error flag for initialization checks */
  HierPartParams hpp;     /* hierarchical partitioning parameters,
                             mainly things that will be needed in the
                             callbacks (a pointer to this is passed as
                             the user data) */
  char msg[256];
  char *yo = "Zoltan_Hier";
  int i;
  int num_obj, gno_size_for_dd, graph_type = 0;
  int userDataLen, hier_changes, hier_num_gid_entries, hier_num_lid_entries;
  int hier_num_import_objs, hier_num_export_objs;

  MPI_Datatype mpi_gno_datatype;
  Zoltan_DD_Directory *dd=NULL;
  int *input_parts=NULL;
  int *hier_import_procs=NULL, *hier_import_to_part=NULL;
  int *hier_export_procs=NULL, *hier_export_to_part=NULL;
  int *fromProc = NULL, *toPart = NULL, *fProc = NULL;
  int nImport;
  ZOLTAN_ID_TYPE *global_ids=NULL, *local_ids=NULL;
  ZOLTAN_ID_TYPE *inGids=NULL, *inLids=NULL, *appids=NULL;
  ZOLTAN_ID_TYPE *hier_import_gids=NULL, *hier_import_lids=NULL;
  ZOLTAN_ID_TYPE *hier_export_gids=NULL, *hier_export_lids=NULL;
  ZOLTAN_ID_TYPE *id_list, *gid, *lid;
  ZOLTAN_GNO_TYPE *gnoList = NULL, *vtxdist=NULL;
  ZOLTAN_GNO_TYPE localsize, globalsize, gno1;

  /* Timers; should probably switch to Zoltan_Timer later */
  double timeInitSetup = 0.;
  double timeLevelSetup = 0.;
  double timePartSetup = 0.;
  double timePartition = 0.;
  double timeMigDown = 0.;
  double timeMigFinal = 0.;
  double timeImpLists = 0.;
  double timeStart = 0;
  double timeEnd = 0;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialize return arguments */
  *num_imp   = *num_exp   = -1;
  *imp_gids  = *exp_gids  = NULL;
  *imp_lids  = *exp_lids  = NULL;
  *imp_procs = *exp_procs = NULL;
  *imp_to_part = *exp_to_part = NULL;

  /* Initialize hpp structure */
  hpp.output_level=0;
  hpp.checks=0;
  hpp.gen_files=0;
  hpp.num_levels=0;
  hpp.level=0;
  hpp.hier_comm=0;
  hpp.origzz=zz;
  hpp.hierzz=NULL;
  hpp.part_to_compute=0;
  hpp.num_parts=0;
  hpp.use_geom=0;
  hpp.use_graph=0;
  hpp.num_obj=0;
  hpp.obj_wgt_dim=zz->Obj_Weight_Dim;
  hpp.edge_wgt_dim=zz->Edge_Weight_Dim;
  hpp.gno=NULL;
  hpp.vwgt=NULL;
  hpp.xadj=NULL;
  hpp.adjncy=NULL;
  hpp.ewgts=NULL;
  hpp.adjproc=NULL;
  hpp.ndims=0;
  hpp.geom_vec=NULL;
  hpp.spec=NULL;
  hpp.use_timers=0;

  /* Cannot currently do hierarchical balancing for num_parts != num_procs */
  if ((zz->Num_Proc != zz->LB.Num_Global_Parts) ||
      (!zz->LB.Single_Proc_Per_Part)) {
    ZOLTAN_HIER_ERROR(ZOLTAN_FATAL,
     "number_parts != number_processes not yet supported by LB_METHOD HIER");
  }

  /* Initialize hierarchical partitioning parameters. */
  ierr = Zoltan_Hier_Initialize_Params(zz, &hpp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_HIER_ERROR(ierr, "Zoltan_Hier_Initialize_Params returned error");
  }

  if (hpp.use_timers) {
    MPI_Barrier(zz->Communicator);
    timeStart = Zoltan_Time(zz->Timer);
  }

  if (!hpp.spec){

    /* Make sure we have the callbacks we need */

    if (zz->Get_Hier_Num_Levels == NULL) {
      ZOLTAN_HIER_ERROR(ZOLTAN_FATAL,"Must register ZOLTAN_HIER_NUM_LEVELS_FN");
    }
    if (zz->Get_Hier_Part == NULL) {
      ZOLTAN_HIER_ERROR(ZOLTAN_FATAL, "Must register ZOLTAN_HIER_PART_FN");
    }
    if (zz->Get_Hier_Method == NULL) {
      ZOLTAN_HIER_ERROR(ZOLTAN_FATAL, "Must register ZOLTAN_HIER_METHOD_FN");
    }
  }
  else{

    /* Zoltan defines the callbacks based on the network topology */

    Zoltan_Set_Hier_Num_Levels_Fn(zz, Zoltan_Hier_Assist_Num_Levels,
                                 (void *)hpp.spec);
    Zoltan_Set_Hier_Part_Fn(zz, Zoltan_Hier_Assist_Part_Number,
                           (void *)hpp.spec);
    Zoltan_Set_Hier_Method_Fn(zz, Zoltan_Hier_Assist_Method, (void *)zz);
  }

  /* do we have callbacks to get geometric and/or graph information? */
  /* TODO KDDKDD Set use_geom and use_graph based on the methods requested */
  hpp.use_geom = ((zz->Get_Geom != NULL) ||
                  (zz->Get_Geom_Multi != NULL));
  hpp.use_graph = ((zz->Get_Num_Edges != NULL) ||
                   (zz->Get_Num_Edges_Multi != NULL));

  /* build our initial intermediate representation */

  ierr = Zoltan_Get_Obj_List(zz, &num_obj,
                             &global_ids, &local_ids,
                             hpp.obj_wgt_dim, &hpp.vwgt, &input_parts);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      ZOLTAN_HIER_ERROR(ierr, "Zoltan_Get_Obj_List returned error.");
  }

  hpp.num_obj = num_obj;
  ZOLTAN_FREE(&input_parts);

  /* build a graph */
  if (hpp.use_graph){
    SET_GLOBAL_GRAPH(&graph_type);
    ierr = Zoltan_Build_Graph(zz, &graph_type,
            hpp.checks, hpp.num_obj,
            global_ids, local_ids,       /* input ZOLTAN_ID_TYPEs */
            hpp.obj_wgt_dim, &hpp.edge_wgt_dim,
            &vtxdist, &hpp.xadj, &hpp.adjncy, /* internal ZOLTAN_GNO_TYPEs */
            &hpp.ewgts, &hpp.adjproc);

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      ZOLTAN_HIER_ERROR(ierr, "Zoltan_Build_Graph returned error.");
    }
    hpp.invalid_gno = vtxdist[zz->Num_Proc];
    gno1 = vtxdist[zz->Proc];
    ZOLTAN_FREE(&vtxdist);
  }
  else{
    MPI_Datatype MPI_GNOTYPE = Zoltan_mpi_gno_type();
    ZOLTAN_GNO_TYPE localCount = hpp.num_obj, scanCount=0;
    ZOLTAN_GNO_TYPE totalCount;
    MPI_Scan(&localCount, &scanCount, 1, MPI_GNOTYPE, MPI_SUM, zz->Communicator);
    totalCount = scanCount;
    MPI_Bcast(&totalCount, 1, MPI_GNOTYPE, zz->Num_Proc-1, zz->Communicator);

    hpp.invalid_gno = totalCount;
    gno1 = scanCount - localCount;
  }

  /* Check that the space of global numbers fits in
   * a ZOLTAN_ID_TYPE.  Caller is
   * using tuples of ZOLTAN_ID_TYPE, which we mapped
   * to singleton ZOLTAN_GNO_TYPE
   * global numbers.  But the Zoltan query functions will store them in
   * singleton ZOLTAN_ID_TYPEs.
   */

  if (sizeof(ZOLTAN_ID_TYPE) < sizeof(ZOLTAN_GNO_TYPE)){
    mpi_gno_datatype = Zoltan_mpi_gno_type();
    localsize = hpp.num_obj;

    MPI_Allreduce(&localsize, &globalsize, 1, mpi_gno_datatype, MPI_SUM,
                  zz->Communicator);

    if (globalsize >= ZOLTAN_ID_INVALID){
      if (zz->Proc == 0){
        fprintf(stderr,"data type for ZOLTAN_ID_TYPE is too small\n");
      }
      ierr = ZOLTAN_FATAL;
      goto End;
    }
  }

  if (hpp.num_obj){
    hpp.gno = (ZOLTAN_GNO_TYPE *)
               ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * hpp.num_obj);
    if (!hpp.gno)
      ZOLTAN_HIER_ERROR(ZOLTAN_MEMERR, "Out of memory");
  }

  for (i=0; i < hpp.num_obj; i++) {
    hpp.gno[i] = gno1 + i;
  }

  /* if we're going to need coordinates */
  if (hpp.use_geom) {
    /* Get coordinate information */
    ierr = Zoltan_Get_Coordinates(zz, hpp.num_obj, global_ids,
                                  local_ids, &hpp.ndims, &hpp.geom_vec);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_HIER_ERROR(ierr, "Error returned from Zoltan_Get_Coordinates");
    }
  }

  /* find out how many levels of hierarchy this proc will participate in */
  hpp.num_levels = zz->Get_Hier_Num_Levels(zz->Get_Hier_Num_Levels_Data,
                                           &ierr);

  if (hpp.output_level >= HIER_DEBUG_ALL) {
    printf("HIER: Proc %d to compute %d levels\n", zz->Proc, hpp.num_levels);
  }

  /* initialize our communicator to the "world" as seen by Zoltan */
  MPI_Comm_dup(zz->Communicator, &hpp.hier_comm);

  if (hpp.use_timers) {
    MPI_Barrier(zz->Communicator);
    timeEnd = Zoltan_Time(zz->Timer);
    timeInitSetup += timeEnd - timeStart;
  }

  /* loop over levels of hierarchical balancing to be done */
  for (hpp.level = 0; hpp.level < hpp.num_levels; hpp.level++) {

    if (hpp.use_timers) {
      MPI_Barrier(hpp.hier_comm);
      timeStart = Zoltan_Time(zz->Timer);
    }

    /* determine parts to compute at this level */
    hpp.part_to_compute =
      zz->Get_Hier_Part(zz->Get_Hier_Part_Data, hpp.level, &ierr);
    /* number of parts is one more than the highest part id
       specified on procs in the current hier_comm */
    MPI_Allreduce(&hpp.part_to_compute, &hpp.num_parts, 1, MPI_INT,
                  MPI_MAX, hpp.hier_comm);
    hpp.num_parts++;

    if (hpp.num_parts == 1){
      /*
       * If there is only one part, and we've not done any partitioning
       * yet, we can skip this step.
       */
      MPI_Comm_size(hpp.hier_comm, &i);
      if (i == zz->Num_Proc)
        continue;
    }

    if (hpp.output_level >= HIER_DEBUG_ALL ||
        zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
      printf("HIER: Proc %d computing part %d of %d at level %d\n",
             zz->Proc, hpp.part_to_compute, hpp.num_parts, hpp.level);
    }

    /* should make sure we have reasonable parts to compute */

    hpp.hierzz = NULL;

    if (hpp.use_timers) {
      MPI_Barrier(hpp.hier_comm);
      timeEnd = Zoltan_Time(zz->Timer);
      timeLevelSetup += timeEnd - timeStart;
    }

    if (hpp.num_parts > 1){

      if (hpp.use_timers) {
        MPI_Barrier(hpp.hier_comm);
        timeStart = Zoltan_Time(zz->Timer);
      }

      /* construct appropriate ZZ and input arrays */
      /* create a brand new one */
      hpp.hierzz = Zoltan_Create(hpp.hier_comm);

      /* and copy in some specified params from zz where appropriate */

      /* just copy debug level to child Zoltan_Struct, use can override
         by setting params of the hierzz in the Get_Hier_Method callback */

      hpp.hierzz->Debug_Level = zz->Debug_Level;
      hpp.hierzz->Timer = zz->Timer;
      hpp.hierzz->Deterministic = zz->Deterministic;
      hpp.hierzz->Obj_Weight_Dim = zz->Obj_Weight_Dim;
      hpp.hierzz->Edge_Weight_Dim = zz->Edge_Weight_Dim;

      /* remapping does not make sense for internal steps, only at the end */
      hpp.hierzz->LB.Remap_Flag = 0;

      /* let the application specify any balancing params for this level */
      zz->Get_Hier_Method(zz->Get_Hier_Method_Data, hpp.level,
                          hpp.hierzz, &ierr);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Get_Hier_Method callback returned error.");
      }

      /* set the numbers of parts */
      sprintf(msg, "%d", hpp.num_parts);
      Zoltan_Set_Param(hpp.hierzz, "NUM_GLOBAL_PARTS", msg);

      /* specify the callbacks */

      ierr = Zoltan_Set_Num_Obj_Fn(hpp.hierzz, Zoltan_Hier_Num_Obj_Fn,
                                   (void *) &hpp);

      ierr = Zoltan_Set_Obj_List_Fn(hpp.hierzz, Zoltan_Hier_Obj_List_Fn,
                                    (void *) &hpp);

      if (hpp.use_geom) {

        ierr = Zoltan_Set_Num_Geom_Fn(hpp.hierzz, Zoltan_Hier_Num_Geom_Fn,
                                      (void *) &hpp);

        ierr = Zoltan_Set_Geom_Multi_Fn(hpp.hierzz, Zoltan_Hier_Geom_Multi_Fn,
                                        (void *) &hpp);
      }

      if (hpp.use_graph) {

        ierr = Zoltan_Set_Num_Edges_Multi_Fn(hpp.hierzz,
                                             Zoltan_Hier_Num_Edges_Multi_Fn,
                                             (void *) &hpp);

        ierr = Zoltan_Set_Edge_List_Multi_Fn(hpp.hierzz,
                                             Zoltan_Hier_Edge_List_Multi_Fn,
                                             (void *) &hpp);
      }

      /* specify the GIDs (just the global numbering) */
      Zoltan_Set_Param(hpp.hierzz, "NUM_GID_ENTRIES", "1");
      Zoltan_Set_Param(hpp.hierzz, "NUM_LID_ENTRIES", "1");

      Zoltan_Set_Param(hpp.hierzz, "RETURN_LISTS", "EXPORT");

      /* deal with part sizes, etc */
      /* we have the assumption here that the final result is one
         part per process */

      ierr = set_hier_part_sizes(&hpp, part_sizes);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_HIER_ERROR(ierr, "set_hier_part_sizes returned error");
      }

      if (hpp.gen_files){
        sprintf(msg,"level_%d",hpp.level);
        ierr = Zoltan_Generate_Files(hpp.hierzz, msg, zz->Proc, 0, 1, 0);
        if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
          ZOLTAN_HIER_ERROR(ierr, "Zoltan_Generate_Files returned error.");
        }
      }

      if (hpp.use_timers) {
        MPI_Barrier(hpp.hier_comm);
        timeEnd = Zoltan_Time(zz->Timer);
        timePartSetup += timeEnd - timeStart;
        MPI_Barrier(hpp.hier_comm);
        timeStart = Zoltan_Time(zz->Timer);
      }

      /* call partitioning method to compute the part at this level */

      ierr = Zoltan_LB_Partition(hpp.hierzz, &hier_changes,
                                 &hier_num_gid_entries, &hier_num_lid_entries,
                                 &hier_num_import_objs,
                                 &hier_import_gids, &hier_import_lids,
                                 &hier_import_procs, &hier_import_to_part,
                                 &hier_num_export_objs,
                                 &hier_export_gids, &hier_export_lids,
                                 &hier_export_procs, &hier_export_to_part);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_LB_Partition returned error.");
      }
      if (hpp.use_timers) {
        MPI_Barrier(hpp.hier_comm);
        timeEnd = Zoltan_Time(zz->Timer);
        timePartition += timeEnd - timeStart;
      }
    }
    else{
      hier_changes=0;
      hier_num_gid_entries = 1; hier_num_lid_entries = 0;
      hier_num_import_objs = hier_num_export_objs = 0;
      hier_import_gids = hier_import_lids=NULL;
      hier_import_procs = hier_import_to_part=NULL;
      hier_export_gids = hier_export_lids=NULL;
      hier_export_procs = hier_export_to_part=NULL;
    }

    if (hpp.level < hpp.num_levels - 1){
      /*
       * Compute the next level of sub communicators
       * and migrate objects downward.
       */
      MPI_Comm next_comm;

      if (hpp.use_timers) {
        MPI_Barrier(hpp.hier_comm);
        timeStart = Zoltan_Time(zz->Timer);
      }

      MPI_Comm_split(hpp.hier_comm, hpp.part_to_compute, 0, &next_comm);

      ierr = migrate_to_next_subgroups(&hpp,
                                       hier_num_export_objs, hier_export_lids,
                                       hier_export_procs, next_comm);
      if (ierr != ZOLTAN_OK)
        goto End;

      MPI_Comm_free(&hpp.hier_comm);
      hpp.hier_comm = next_comm;

      if (hpp.use_timers) {
        MPI_Barrier(hpp.hier_comm);
        timeEnd = Zoltan_Time(zz->Timer);
        timeMigDown += timeEnd - timeStart;
      }
    }
    else{
      if (hpp.use_timers) {
        MPI_Barrier(hpp.hier_comm);
        timeStart = Zoltan_Time(zz->Timer);
      }

      /*
       * Migrate objects (without weights or adjacencies) their new owners.  Now
       * each process has a list representing the objects to be imported to it.
       */
      ierr = final_migrate(&hpp,
                           hier_num_export_objs, hier_export_gids,
                            hier_export_lids, hier_export_procs);

      if (ierr != ZOLTAN_OK)
        goto End;

      if (hpp.use_timers) {
        MPI_Barrier(hpp.hier_comm);
        timeEnd = Zoltan_Time(zz->Timer);
        timeMigFinal += timeEnd - timeStart;
      }
      /* Still need hpp.gno to create our import list */
      ZOLTAN_FREE(&hpp.vwgt);
      ZOLTAN_FREE(&hpp.xadj);
      ZOLTAN_FREE(&hpp.adjncy);
      ZOLTAN_FREE(&hpp.ewgts);
      ZOLTAN_FREE(&hpp.adjproc);
      ZOLTAN_FREE(&hpp.geom_vec);
      MPI_Comm_free(&hpp.hier_comm);
      hpp.hier_comm = MPI_COMM_NULL;
    }

    ierr = Zoltan_LB_Free_Part(&hier_import_gids, &hier_import_lids,
                               &hier_import_procs, &hier_import_to_part);
    if (ierr != ZOLTAN_OK){
      ZOLTAN_HIER_ERROR(ierr, "Zoltan_LB_Free_Part returned error.");
    }

    ierr = Zoltan_LB_Free_Part(&hier_export_gids, &hier_export_lids,
                               &hier_export_procs, &hier_export_to_part);
    if (ierr != ZOLTAN_OK){
      ZOLTAN_HIER_ERROR(ierr, "Zoltan_LB_Free_Part returned error.");
    }

    /* clean up hierzz */
    Zoltan_Destroy(&hpp.hierzz);
  }

  /*
   * Fill out import lists.
   *    hpp.num_obj - the number of objects I'm left with after partitioning
   *    num_obj     - the number of objects I had before partitioning
   */
  if (hpp.use_timers) {
    MPI_Barrier(zz->Communicator);
    timeStart = Zoltan_Time(zz->Timer);
  }

  userDataLen = sizeof(ZOLTAN_ID_TYPE) * (zz->Num_GID + zz->Num_LID);
  gno_size_for_dd = sizeof(ZOLTAN_GNO_TYPE) / sizeof(ZOLTAN_ID_TYPE);

  ierr = Zoltan_DD_Create(&dd, zz->Communicator, gno_size_for_dd, 0,
                          userDataLen, num_obj, 0);

  if (ierr != ZOLTAN_OK)
    goto End;


  if (num_obj){
    gnoList = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE)*num_obj);
    if (!gnoList){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  i = (hpp.num_obj > num_obj ? hpp.num_obj : num_obj);

  if (i > 0){
    appids = (ZOLTAN_ID_TYPE *)ZOLTAN_MALLOC(userDataLen * i);
    if (!appids){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  id_list = appids;
  gid = global_ids;
  lid = local_ids;

  for (i=0; i < num_obj; i++){
    gnoList[i] = gno1 + i;
    ZOLTAN_SET_GID(zz, id_list, gid);
    id_list+= zz->Num_GID;
    gid+= zz->Num_GID;

    if (zz->Num_LID){
      ZOLTAN_SET_LID(zz, id_list, lid);
      id_list+= zz->Num_LID;
      lid+= zz->Num_LID;
    }
  }

  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);

  ierr = Zoltan_DD_Update(dd, (ZOLTAN_ID_TYPE *)gnoList, NULL,
                          (char *)appids, NULL, num_obj);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&gnoList);

  if (hpp.num_obj > 0){
    fromProc = (int *)ZOLTAN_MALLOC(sizeof(int) * hpp.num_obj);
    if (!fromProc){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  ierr = Zoltan_DD_Find(dd, (ZOLTAN_ID_TYPE *)hpp.gno, NULL,
                        (char *)appids, NULL, hpp.num_obj, fromProc);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&hpp.gno);
  Zoltan_DD_Destroy(&dd);

  nImport = 0;
  for (i = 0; i < hpp.num_obj; i++)
    if (fromProc[i] != zz->Proc) /* Importing object; assume #parts == #procs */
      nImport++;

  *num_imp = nImport;

  if (nImport > 0){
    if (!Zoltan_Special_Malloc(zz, (void **)imp_gids, nImport,
                               ZOLTAN_SPECIAL_MALLOC_GID)
     || !Zoltan_Special_Malloc(zz, (void **)imp_lids, nImport,
                               ZOLTAN_SPECIAL_MALLOC_LID)
     || !Zoltan_Special_Malloc(zz, (void **)imp_procs, nImport,
                               ZOLTAN_SPECIAL_MALLOC_INT)
     || !Zoltan_Special_Malloc(zz, (void **)imp_to_part, nImport,
                               ZOLTAN_SPECIAL_MALLOC_INT)) {
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    toPart = *imp_to_part;
    inGids = *imp_gids;
    inLids = *imp_lids;
    fProc = *imp_procs;
  }

  gid = inGids;
  lid = inLids;
  id_list = appids;

  nImport = 0;
  for (i=0; i < hpp.num_obj; i++){
    if (fromProc[i] != zz->Proc) { /* Importing obj; assume #parts==#procs */
      toPart[nImport] = zz->Proc;

      ZOLTAN_SET_GID(zz, gid, id_list);
      id_list += zz->Num_GID;
      gid += zz->Num_GID;

      ZOLTAN_SET_LID(zz, lid, id_list);
      id_list += zz->Num_LID;
      lid += zz->Num_LID;

      fProc[nImport] = fromProc[i];
      nImport++;
    }
    else {
      id_list += zz->Num_GID + zz->Num_LID;
    }
  }

  ZOLTAN_FREE(&appids);

  if (hpp.use_timers) {
    MPI_Barrier(zz->Communicator);
    timeEnd = Zoltan_Time(zz->Timer);
    timeImpLists += timeEnd - timeStart;
  }

End:

  if (hpp.use_timers && zz->Proc == 0) {
    printf("HIER TIME InitSetup  %e\n", timeInitSetup);
    printf("HIER TIME LevelSetup %e\n", timeLevelSetup);
    printf("HIER TIME PartSetup  %e\n", timePartSetup);
    printf("HIER TIME Partition  %e\n", timePartition);
    printf("HIER TIME MigDown    %e\n", timeMigDown);
    printf("HIER TIME MigFinal   %e\n", timeMigFinal);
    printf("HIER TIME ImpLists   %e\n", timeImpLists);
  }
  ZOLTAN_FREE(&vtxdist);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&input_parts);
  ZOLTAN_FREE(&gnoList);
  ZOLTAN_FREE(&fromProc);
  ZOLTAN_FREE(&appids);
  ZOLTAN_FREE(&hpp.gno);
  ZOLTAN_FREE(&hpp.vwgt);
  ZOLTAN_FREE(&hpp.xadj);
  ZOLTAN_FREE(&hpp.adjncy);
  ZOLTAN_FREE(&hpp.ewgts);
  ZOLTAN_FREE(&hpp.adjproc);
  ZOLTAN_FREE(&hpp.geom_vec);
  if (dd) Zoltan_DD_Destroy(&dd);
  Zoltan_Destroy(&hpp.hierzz);
  if (hpp.hier_comm != MPI_COMM_NULL) MPI_Comm_free(&hpp.hier_comm);
  if (hpp.spec != NULL) ZOLTAN_FREE(&hpp.spec);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
