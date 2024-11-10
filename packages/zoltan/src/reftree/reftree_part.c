// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include "zz_const.h"
#include "reftree.h"
#include "all_allo_const.h"
#include "params_const.h"

/* TEMP SINGLE WEIGHT
        Currently, it has not been decided how to handle multi-component
        weights.  This code just uses the first component.  To find all
        locations where the code would change if a different reduction
        is defined, search for "TEMP SINGLE WEIGHT" */

/* TEMP HA Support for heterogeneous architectures.
        To find all places where changes may be needed to support
        heterogeneous architectures, search for "TEMP HA".
        Support is provided for different processing power by producing
        unequal sized partitions.
        Support is NOT provided for unequal communication.  I don't know
        how to do that.  Presumably it would involve mapping the partitions
        to the processors in a way that the smaller boundary interfaces go
        to the slower communication channels.  Currently this algorithm
        doesn't even consider the mapping of partitions to processors.  */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Prototypes for functions internal to this file */
static int Zoltan_Reftree_Sum_Weights(ZZ *zz);

static void Zoltan_Reftree_Sum_My_Weights(ZZ *zz, ZOLTAN_REFTREE *subroot, 
       int *count, int wdim);
static void Zoltan_Reftree_Sum_All_Weights(ZZ *zz, ZOLTAN_REFTREE *subroot, int wdim);
static void Zoltan_Reftree_List_Other_Leaves(ZZ *zz, ZOLTAN_REFTREE *subroot, 
       ZOLTAN_ID_PTR list, int *count);
static int Zoltan_Reftree_Partition(ZZ *zz, float *part_sizes, int *num_export, 
       ZOLTAN_ID_PTR *export_global_ids, ZOLTAN_ID_PTR *export_local_ids, 
       int **export_to_partition, int **export_procs);
static int Zoltan_Reftree_Part_Recursive(ZZ *zz, ZOLTAN_REFTREE *subroot, int *part,
       float *current_size, int *num_exp, float *cutoff,
       int num_part);
static int Zoltan_Reftree_Mark_and_Count(ZOLTAN_REFTREE *subroot, int part, 
       int *num_exp, ZZ *zz);
static int Zoltan_Reftree_Export_Lists(ZZ *zz, ZOLTAN_REFTREE *subroot, 
       int *num_export, ZOLTAN_ID_PTR *export_global_ids,
       ZOLTAN_ID_PTR *export_local_ids, int **export_to_partition,
       int **export_procs);
static int export_it(ZOLTAN_REFTREE *subroot, ZZ *zz, int *ierr);
static int get_current_part(ZOLTAN_REFTREE *subroot, ZZ *zz, int *ierr);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Reftree_Part(

  ZZ *zz,                       /* The Zoltan structure */
  float *part_sizes,            /* Input:  Array of size zz->Num_Global_Parts
                                   containing the percentage of work to be
                                   assigned to each partition.               */
  int *num_import,              /* Not computed, set to -1 */
  ZOLTAN_ID_PTR *import_global_ids, /* Not computed */
  ZOLTAN_ID_PTR *import_local_ids,  /* Not computed */
  int **import_procs,           /* Not computed */
  int **import_to_part,         /* Not computed */
  int *num_export,              /* Number of objects to be exported */
  ZOLTAN_ID_PTR *export_global_ids, /* global ids of objects to be exported */
  ZOLTAN_ID_PTR *export_local_ids,  /* local  ids of objects to be exported */
  int **export_procs,           /* list of processors to export to */
  int **export_to_partition     /* list of partitions to export to */
)
{
char *yo = "Zoltan_Reftree_Part";
int ierr;       /* error code returned by called routines */
int final_ierr; /* error code returned by this routine */
double time0 = 0, time1= 0, time2 = 0, time3 =0, time4 =0;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initializations in case of early exit. */
  *num_export = -1;
  *num_import = -1;
  final_ierr = ZOLTAN_OK;

  /*
   * initialize the tree (first call only)
   */

  if (zz->LB.Data_Structure == NULL) {
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) time0 = Zoltan_Time(zz->Timer);
    ierr = Zoltan_Reftree_Init(zz);
    if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                     "Error returned by Zoltan_Reftree_Init.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return(ierr);
    }
    if (ierr==ZOLTAN_WARN) final_ierr = ZOLTAN_WARN;
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) time1 = Zoltan_Time(zz->Timer);
  } else {
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
      time1 = Zoltan_Time(zz->Timer);
      time0 = time1 + 1.0;
    }
  }

  /*
   * build the refinement tree
   */

  ierr = Zoltan_Reftree_Build(zz);
  if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Error returned by Zoltan_Reftree_Build.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }
  if (ierr==ZOLTAN_WARN) final_ierr = ZOLTAN_WARN;
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) time2 = Zoltan_Time(zz->Timer);

  /*
   * sum the weights in the tree
   */

  ierr = Zoltan_Reftree_Sum_Weights(zz);
  if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Error returned by Zoltan_Reftree_Sum_Weights.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }
  if (ierr==ZOLTAN_WARN) final_ierr = ZOLTAN_WARN;
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) time3 = Zoltan_Time(zz->Timer);

  /*
   * determine the new partition
   */

  ierr = Zoltan_Reftree_Partition(zz, part_sizes, num_export, export_global_ids,
                          export_local_ids, export_to_partition, export_procs);
  if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Error returned by Zoltan_Reftree_Partition.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }
  if (ierr==ZOLTAN_WARN) final_ierr = ZOLTAN_WARN;
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) time4 = Zoltan_Time(zz->Timer);

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
    if (time0 <= time1) {
      Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, time1-time0,
                     "REFTREE Time to initialize :");
    }
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, time2-time1, 
                   "REFTREE Time to build tree :");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, time3-time2,
                   "REFTREE Time to sum weights:");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, time4-time3,
                   "REFTREE Time to partition  :");
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return(final_ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_Reftree_Sum_Weights(ZZ *zz)

{
/*
 * Function to sum the weights in the refinement tree.  On input the
 * refinement tree should be valid and have weight set.  On output the
 * values in summed_weight at each node is the sum of the weights in the
 * subtree with that node as root.
 * This function also sets assigned_to_me for interior nodes to be
 * 1 if the entire subtree is assigned to this processor
 * 0 if none of the subtree is assigned to this processor
 * -1 if some of the subtree is assigned to this processor
 */
char *yo = "Zoltan_Reftree_Sum_Weights";
ZOLTAN_REFTREE *root;         /* Root of the refinement tree */
int wdim;                 /* Dimension of the weight array */
int i,j;                  /* loop counters */
int count;                /* counter */
ZOLTAN_ID_PTR leaf_list = NULL;      
                          /* leaves for which some proc requests weight */
ZOLTAN_ID_PTR all_leaflist = NULL;   
                          /* leaf_list from all processors */
int reqsize;              /* length of leaf_list */
int *reqsize_all;         /* reqsize from all processors */
int sum_reqsize;          /* sum of all reqsize */
int *displs;              /* running sum of all reqsize */
int my_start;             /* position in leaf_list of this proc's list */
int nproc;                /* number of processors */
ZOLTAN_REFTREE *node;         /* a node in the refinement tree */
struct Zoltan_Reftree_hash_node **hashtab; /* hash table */
int hashsize;             /* dimension of hash table */
float *send_float;        /* sending message of floats */
float *req_weights;       /* the requested weights */
int num_gid_entries = zz->Num_GID; /* Number of array entries in a global ID */

   ZOLTAN_TRACE_ENTER(zz, yo);

  /*
   * set the root and hash table
   */

  root = ((struct Zoltan_Reftree_data_struct *)zz->LB.Data_Structure)->reftree_root;
  if (root == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Refinement tree not defined.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_FATAL);
  }
  hashtab  = ((struct Zoltan_Reftree_data_struct *)zz->LB.Data_Structure)->hash_table;
  hashsize = ((struct Zoltan_Reftree_data_struct *)zz->LB.Data_Structure)->hash_table_size;

  /*
   * Determine the dimension of the weight array
   */

  if (zz->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = zz->Obj_Weight_Dim;
  }

  /*
   * In the first pass, sum the weights of the nodes that are assigned to
   * this processor, and count the leaves that are not.
   */

  count = 0;
  for (i=0; i<root->num_child; i++) {
    Zoltan_Reftree_Sum_My_Weights(zz,&(root->children[i]),&count,wdim);
  }
  root->assigned_to_me = -1;

  /*
   * Make a list of the leaves that are not assigned to this processor
   */

  if (count == 0)
    leaf_list = ZOLTAN_MALLOC_GID(zz);
  else
    leaf_list = ZOLTAN_MALLOC_GID_ARRAY(zz, count);
  if (leaf_list == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_MEMERR);
  }

  count = 0;
  Zoltan_Reftree_List_Other_Leaves(zz, root,leaf_list,&count);

  /*
   * Get the unknown leaf weights from other processors.
   */

  nproc = zz->Num_Proc;
  reqsize = count;

  /*
   * Build a list of all processor's request list by concatinating them in
   * the order of the processor ranks
   */

  /*
   * Determine the request size of all processors
   */

  reqsize_all = (int *)ZOLTAN_MALLOC(nproc*sizeof(int));
  displs = (int *)ZOLTAN_MALLOC(nproc*sizeof(int));
  if (reqsize_all == NULL || displs == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    Zoltan_Multifree(__FILE__, __LINE__, 3, &displs,
                                            &reqsize_all,
                                            &leaf_list);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_MEMERR);
  }

  MPI_Allgather((void *)&reqsize,1,MPI_INT,(void *)reqsize_all,1,MPI_INT,
                zz->Communicator);
  displs[0] = 0;
  for (i=1; i<nproc; i++) displs[i] = displs[i-1]+reqsize_all[i-1];
  sum_reqsize = displs[nproc-1] + reqsize_all[nproc-1];
  my_start = displs[zz->Proc];

  /*
   * If sum_reqsize is 0, nothing needs to be communciated
   */

  if (sum_reqsize == 0) {
    Zoltan_Multifree(__FILE__, __LINE__, 3, &displs,
                                            &reqsize_all,
                                            &leaf_list);
  }
  else {

  /*
   * Gather the request list from all processors
   */

    all_leaflist = ZOLTAN_MALLOC_GID_ARRAY(zz, sum_reqsize);
    if (all_leaflist == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      Zoltan_Multifree(__FILE__, __LINE__, 4, &all_leaflist,
                                              &displs,
                                              &reqsize_all,
                                              &leaf_list);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return(ZOLTAN_MEMERR);
    }

    /* KDDKDD Changed MPI_BYTE to ZOLTAN_ID_MPI_TYPE  */

    /* Account for number of array entries in an ID. */
    for (i=0; i<nproc; i++) {
      reqsize_all[i] = reqsize_all[i]*num_gid_entries;
      displs[i] = displs[i]*num_gid_entries;
    }

    MPI_Allgatherv((void *)leaf_list,reqsize*num_gid_entries,ZOLTAN_ID_MPI_TYPE,
                   (void *)all_leaflist,reqsize_all,displs,ZOLTAN_ID_MPI_TYPE,
                   zz->Communicator);

    ZOLTAN_FREE(&displs);
    ZOLTAN_FREE(&leaf_list);

    for (i=0; i<nproc; i++) reqsize_all[i] = reqsize_all[i]/num_gid_entries;

  /* 
   * Create a list with the partial sums this processor has
   */

    send_float = (float *) ZOLTAN_MALLOC(sizeof(float)*wdim*sum_reqsize);
    if (send_float == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      Zoltan_Multifree(__FILE__, __LINE__, 3, &send_float,
                                              &all_leaflist,
                                              &reqsize_all);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return(ZOLTAN_MEMERR);
    }

    for (i=0; i<sum_reqsize; i++) {
      node = Zoltan_Reftree_hash_lookup(zz, hashtab,
                                    &(all_leaflist[i*num_gid_entries]),
                                    hashsize);
      if (node == NULL)
         for (j=0; j<wdim; j++) send_float[i*wdim+j] = 0.0;
      else
         for (j=0; j<wdim; j++) send_float[i*wdim+j] = node->my_sum_weight[j];
    }

  /*
   * Sum the weights over all the processors
   */

    if (reqsize == 0)
      req_weights = (float *) ZOLTAN_MALLOC(sizeof(float)*wdim);
    else
      req_weights = (float *) ZOLTAN_MALLOC(sizeof(float)*wdim*reqsize);
    if (req_weights == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      Zoltan_Multifree(__FILE__, __LINE__, 4, &req_weights,
                                              &send_float,
                                              &all_leaflist,
                                              &reqsize_all);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return(ZOLTAN_MEMERR);
    }

    MPI_Reduce_scatter((void *)send_float, (void *)req_weights, reqsize_all,
                       MPI_FLOAT, MPI_SUM, zz->Communicator);

    ZOLTAN_FREE(&send_float);
    ZOLTAN_FREE(&reqsize_all);

  /*
   * Set the weights this processor requested
   */

    for (i=0; i<count; i++) {
      node = Zoltan_Reftree_hash_lookup(zz, hashtab,
                                  &(all_leaflist[(i+my_start)*num_gid_entries]),
                                  hashsize);
      for (j=0; j<wdim; j++) node->summed_weight[j] = req_weights[i*wdim+j];
    }

    ZOLTAN_FREE(&req_weights);
    ZOLTAN_FREE(&all_leaflist);
  }

  /*
   * All the leaves now have summed_weight set.
   * Sum the weights throughout the tree.
   */

  Zoltan_Reftree_Sum_All_Weights(zz,root,wdim);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Reftree_Sum_My_Weights(ZZ *zz, ZOLTAN_REFTREE *subroot, 
       int *count, int wdim)

{
/*
 * Function to recursively sum the weights of the nodes assigned
 * to this processor, and set assigned_to_me for interior nodes
 */
int i, j;          /* loop counter */
int none_assigned; /* flag for no children assigned to this proc */
int all_assigned;  /* flag for all children assigned to this proc */

  if (subroot->num_child == 0) {

  /*
   * If there are no children, then the sum is the weight of this node if
   * it is assigned to this processor, or 0 if it is not.
   */

    if (subroot->assigned_to_me) {
      for (i=0; i<wdim; i++) subroot->my_sum_weight[i] = subroot->weight[i];
      for (i=0; i<wdim; i++) subroot->summed_weight[i] = subroot->weight[i];
    }
    else {
      for (i=0; i<wdim; i++) subroot->my_sum_weight[i] = 0.0;
      for (i=0; i<wdim; i++) subroot->summed_weight[i] = 0.0;
      *count += 1;
    }

  }
  else {

  /*
   * If there are children, sum the weights of the children with the
   * node's weight and set assigned to me.
   */

    if (subroot->assigned_to_me) {
      for (i=0; i<wdim; i++) subroot->my_sum_weight[i] = subroot->weight[i];
    }
    else {
      for (i=0; i<wdim; i++) subroot->my_sum_weight[i] = 0.0;
    }
    none_assigned = 1;
    all_assigned = 1;

    for (j=0; j<subroot->num_child; j++) {
      Zoltan_Reftree_Sum_My_Weights(zz,&(subroot->children[j]),count,wdim);
      for (i=0; i<wdim; i++)
        subroot->my_sum_weight[i] += (subroot->children[j]).my_sum_weight[i];
      if ((subroot->children[j]).assigned_to_me == 1) none_assigned = 0;
      if ((subroot->children[j]).assigned_to_me == 0) all_assigned = 0;
      if ((subroot->children[j]).assigned_to_me == -1) {
        none_assigned = 0;
        all_assigned = 0;
      }
    }
    if (none_assigned)
      subroot->assigned_to_me = 0;
    else if (all_assigned)
      subroot->assigned_to_me = 1;
    else
      subroot->assigned_to_me = -1;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Reftree_Sum_All_Weights(ZZ *zz, ZOLTAN_REFTREE *subroot, int wdim)

{
/*
 * Function to recursively sum the weights of all the nodes in the tree
 * assuming that summed_weight contains partial sums in the leaves
 */
int i, j;   /* loop counter */

  /*
   * If there are no children, do nothing
   */

  if (subroot->num_child != 0) {

  /*
   * If there are children, sum the weights of the children with the
   * weight of this node.
   */

    for (i=0; i<wdim; i++) subroot->summed_weight[i] = subroot->weight[i];

    for (j=0; j<subroot->num_child; j++) {
      Zoltan_Reftree_Sum_All_Weights(zz,&(subroot->children[j]),wdim);
      for (i=0; i<wdim; i++)
        subroot->summed_weight[i] += (subroot->children[j]).summed_weight[i];
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Reftree_List_Other_Leaves(ZZ *zz, ZOLTAN_REFTREE *subroot, 
       ZOLTAN_ID_PTR list, int *count)

{
/*
 * Function to make a list of the leaves not assigned to this processor
 */
int j;   /* loop counter */

  if (subroot->num_child == 0) {

  /*
   * If there are no children, then add it to the list if it is not
   * assigned to this processor
   */

    if (!subroot->assigned_to_me) {
      ZOLTAN_SET_GID(zz, &(list[(*count)*zz->Num_GID]),subroot->global_id);
      *count += 1;
    }

  }
  else {

  /*
   * If there are children, traverse the subtrees
   */

    for (j=0; j<subroot->num_child; j++) {
      Zoltan_Reftree_List_Other_Leaves(zz, &(subroot->children[j]),list,count);
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_Reftree_Partition(ZZ *zz, float *part_sizes, int *num_export, 
       ZOLTAN_ID_PTR *export_global_ids, ZOLTAN_ID_PTR *export_local_ids, 
       int **export_to_partition, int **export_procs)

{
/*
 * Function to partition the grid and fill the export arrays.
 */

char *yo = "Zoltan_Reftree_Partition";
char msg[256];
int num_exp;          /* count the number of export objects */
ZOLTAN_REFTREE *root;     /* root of the tree */
float *cutoff;        /* the relative sizes of the partitions */
int part;             /* partition under construction */
float current_size;   /* amount of weight consumed so far */
int num_part;         /* number of partitions */
int ierr;             /* error flag */
int wdim;             /* Max(zz->Obj_Weight_Dim, 1) */

  ZOLTAN_TRACE_ENTER(zz, yo);

  root = ((struct Zoltan_Reftree_data_struct *)zz->LB.Data_Structure)->reftree_root;

  /*
   * determine the size of the partitions and tolerance interval
   */

  num_part = zz->LB.Num_Global_Parts;

/* Set the cutoff points of summed weights for the end of each partition */

/* TEMP HA
        Determine a partition of unity for the sizes of the partitions
        relative to the processing power of the processors, to support
        heterogeneous architectures.
        cutoff(0) = is the percent of work that should be assigned to
        partition 0
        cutoff(i) - cutoff(i-1) is the percent of work that should be
        assigned to partition i.
        cutoff(num_part-1) = 1.0
        When support for heterogeneous architectures is added to Zoltan,
        there should be a "vector of partition sizes" passed into this
        routine, which would be used to define cutoff.  For now, it is
        just set to be equal sizes.
*/

  cutoff = (float *)ZOLTAN_MALLOC((num_part)*sizeof(float));
  if (cutoff == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_MEMERR);
  }

  cutoff[0] = part_sizes[0];  /* TEMP SINGLE WEIGHT */
  wdim = ((zz->Obj_Weight_Dim == 0) ? 1 : zz->Obj_Weight_Dim);
  for (part = 1; part < num_part; part++) 
     cutoff[part] = cutoff[part-1] + part_sizes[part*wdim]; /* TEMP SINGLE WEIGHT */
  cutoff[num_part-1] = 1.0; /* just to make sure roundoff doesn't bite us */

/* multiply cutoff by the total weight so that cutoff gives the desired
   actual weight of each partition instead of relative sizes */

/* TEMP HA
   It is quite possible that once we know what is passed into this
   routine to determine relative partition sizes, then this can be combined
   with the above partition of unity to simplify the code.  I just did it
   this way to make it clear what is happening.
*/

  for (part=0; part<num_part; part++) {
    cutoff[part] = cutoff[part]*root->summed_weight[0]; /* TEMP SINGLE WEIGHT */
  }

  /*
   * traverse the tree to define the partition and count the number of exports
   */

  num_exp = 0;
  part = 0;
  current_size = 0.0;
  ierr = Zoltan_Reftree_Part_Recursive(zz,root,&part,&current_size,&num_exp,
                                       cutoff,num_part);
  ZOLTAN_FREE(&cutoff);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }

  /*
   * if no exports, we're done
   */

  if (zz->LB.Return_Lists == ZOLTAN_LB_NO_LISTS) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_OK);
  }
  else if (zz->LB.Return_Lists == ZOLTAN_LB_CANDIDATE_LISTS) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Candidate Lists not supported in REFTREE;"
                                     "change RETURN_LISTS parameter.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_FATAL);
  }
  else if (num_exp == 0) {
    *num_export = 0;
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_OK);
  }

  /*
   * allocate space for the export lists
   */

  if (!Zoltan_Special_Malloc(zz,(void **)export_global_ids,num_exp,
                         ZOLTAN_SPECIAL_MALLOC_GID)) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }
  if (!Zoltan_Special_Malloc(zz,(void **)export_local_ids,num_exp,
                         ZOLTAN_SPECIAL_MALLOC_LID)) {
    Zoltan_Special_Free(zz,(void **)export_global_ids,ZOLTAN_SPECIAL_MALLOC_GID);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }
  if (!Zoltan_Special_Malloc(zz,(void **)export_to_partition,num_exp,
                         ZOLTAN_SPECIAL_MALLOC_INT)) {
    Zoltan_Special_Free(zz,(void **)export_global_ids,ZOLTAN_SPECIAL_MALLOC_GID);
    Zoltan_Special_Free(zz,(void **)export_local_ids,ZOLTAN_SPECIAL_MALLOC_LID);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }
  if (!Zoltan_Special_Malloc(zz,(void **)export_procs,num_exp,
                         ZOLTAN_SPECIAL_MALLOC_INT)) {
    Zoltan_Special_Free(zz,(void **)export_to_partition,ZOLTAN_SPECIAL_MALLOC_INT);
    Zoltan_Special_Free(zz,(void **)export_global_ids,ZOLTAN_SPECIAL_MALLOC_GID);
    Zoltan_Special_Free(zz,(void **)export_local_ids,ZOLTAN_SPECIAL_MALLOC_LID);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }

  /*
   * traverse the tree to make the export lists
   */

  *num_export = 0;
  ierr = Zoltan_Reftree_Export_Lists(zz,root,num_export,export_global_ids,
                          export_local_ids,export_to_partition,export_procs);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }

  if (num_exp != *num_export) {
    sprintf(msg, "num_exp = %d not equal to num_export = %d.",
            num_exp,*num_export);
    ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_WARN);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_Reftree_Part_Recursive(ZZ *zz, ZOLTAN_REFTREE *subroot, int *part,
                              float *current_size, int *num_exp, float *cutoff,
                              int num_part)

{
/*
 * function to recursively define the partition and count the number of exports
 */
int i;         /* loop counter */
float newsize; /* size of partition if this subroot gets added to it */
float eps;     /* imbalance tolerance in units of weight */
float imb_tol; /* Zoltan imbalance tolerance */
int ierr;      /* error flag */

  imb_tol = zz->LB.Imbalance_Tol[0];  /* Only one weight currently supported. */
  newsize = *current_size + subroot->summed_weight[0]; /* TEMP SINGLE WEIGHT */
  if (*part != num_part-1)
    eps = (imb_tol - 1.0)*(cutoff[*part+1]-cutoff[*part])/2.0;
  else
    eps = 0.0;

  if (newsize <= cutoff[*part] + eps || *part == num_part-1) {

  /*
   * This subtree fits in the current partition
   */

    subroot->partition = *part;
    *current_size = newsize;

  /*
   * If there are no leaves of this subtree assigned to this processor, there
   * are no exports below this node.
   * Otherwise, traverse the subtree setting partition and counting exports
   */

    if (subroot->assigned_to_me) {
      if (subroot->num_child == 0) {
        if (export_it(subroot,zz,&ierr))
          *num_exp += 1;
        if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) return(ierr);
      }
      else
        for (i=0; i<subroot->num_child; i++) {
          ierr = Zoltan_Reftree_Mark_and_Count(&(subroot->children[i]), *part,
                                               num_exp, zz);
          if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) return(ierr);
        }
    }

  /*
   * See if it is close enough to filling the partition
   */

    if (*part != 0)
      eps = (imb_tol - 1.0)*(cutoff[*part]-cutoff[*part-1])/2.0;
    else
      eps = (imb_tol - 1.0)*(cutoff[*part])/2.0;

    if (*current_size >= cutoff[*part] - eps && *part < num_part-1) {
      *part += 1;
    }
  }
  else {

  /*
   * This subtree is too big for the current partition
   */

  /*
   * If it has children, traverse them.
   */

    if (subroot->num_child != 0) {
      subroot->partition = -1;
      for (i=0; i<subroot->num_child; i++)
        ierr = Zoltan_Reftree_Part_Recursive(zz, &(subroot->children[i]),part,
                                     current_size, num_exp, cutoff, num_part);
    }
    else {

  /*
   * If there are no children, move on to the next partition
   */

      if (*part != num_part-1)
        eps = (imb_tol - 1.0)*(cutoff[*part+1]-cutoff[*part])/2.0;
      else
        eps = 0.0;

      while (newsize > cutoff[*part]+eps && *part < num_part-1) {
        *part += 1;
        if (*part != num_part-1)
         eps = (imb_tol - 1.0)*(cutoff[*part+1]-cutoff[*part])/2.0;
        else
         eps = 0.0;
      }
      subroot->partition = *part;
      *current_size = newsize;
      if (export_it(subroot,zz,&ierr)) *num_exp += 1;
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) return(ierr);
    }
  }
  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_Reftree_Mark_and_Count(ZOLTAN_REFTREE *subroot, int part, 
                                         int *num_exp, ZZ *zz)
{
/*
 * Function to set the partition and count exports
 */
int i, ierr;

  subroot->partition = part;
  if (subroot->num_child == 0) {
    if (export_it(subroot,zz,&ierr))
      *num_exp += 1;
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) return(ierr);
  } else {
    for (i=0; i<subroot->num_child; i++) {
      ierr = Zoltan_Reftree_Mark_and_Count(&(subroot->children[i]), part,
                                           num_exp, zz);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) return(ierr);
    }
  }
  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_Reftree_Export_Lists(ZZ *zz, ZOLTAN_REFTREE *subroot, 
                            int *num_export, ZOLTAN_ID_PTR *export_global_ids,
                            ZOLTAN_ID_PTR *export_local_ids,
                            int **export_to_partition, int **export_procs)
{
/*
 * Function to build the export lists
 */
int i, ierr;

/*
 * if this subtree has no leaves assigned to this processor then there can be
 * no exports below it
 */

  if (!subroot->assigned_to_me) return(ZOLTAN_OK);

  if (subroot->num_child == 0) {

/*
 * if this is a leaf, put it on the export lists if it is to be exported
 */

    if (export_it(subroot,zz,&ierr)) {
      ZOLTAN_SET_GID(zz, &((*export_global_ids)[(*num_export)*zz->Num_GID]),
                     subroot->global_id);
      ZOLTAN_SET_LID(zz, &((*export_local_ids)[(*num_export)*zz->Num_LID]),
                     subroot->local_id);
      (*export_to_partition)[*num_export] = subroot->partition;
      (*export_procs)[*num_export] = Zoltan_LB_Part_To_Proc(zz,subroot->partition,subroot->global_id);
      *num_export += 1;
    }
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) return(ierr);
  }
  else {

/*
 * if it is not a leaf, traverse the subtree
 */

    for (i=0; i<subroot->num_child; i++) {
      ierr = Zoltan_Reftree_Export_Lists(zz, &(subroot->children[i]),num_export,
                                         export_global_ids,export_local_ids,
                                         export_to_partition,export_procs);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) return(ierr);
    }
  }
  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int export_it(ZOLTAN_REFTREE *subroot, ZZ *zz, int *ierr)
{
/*
 * Function to determine if an object belongs on the export list.
 */

int current_part;

/* return TRUE if it is currently assigned to this processor and
   either the new partition is not the old partition or the new
   partition is not assigned to this processor */

  current_part = get_current_part(subroot,zz,ierr);
  if (*ierr != ZOLTAN_OK && *ierr != ZOLTAN_WARN) return(FALSE);

  if ((current_part != subroot->partition ||
       Zoltan_LB_Part_To_Proc(zz,subroot->partition,subroot->global_id) != zz->Proc)
       && subroot->assigned_to_me) return(TRUE);

  return(FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int get_current_part(ZOLTAN_REFTREE *subroot, ZZ *zz, int *ierr)
{
/*
 * Function to return the current partition of an object.
 * If there is no user defined get_partition function, then the returned
 * value only indicates whether or not the partition number is this
 * processor's number.
 */

char *yo = "get_current_part";
int result;

  *ierr = ZOLTAN_OK;

/* if the user registered a partition function, then use it */
  if (zz->Get_Part != NULL) {
    result = zz->Get_Part(zz->Get_Part_Data,zz->Num_GID,zz->Num_LID,
                               subroot->global_id,subroot->local_id, ierr);
    if (*ierr != ZOLTAN_OK && *ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from ZOLTAN_PART_FN");
    }
  }
  else if (zz->Get_Part_Multi != NULL) {
/* Not the best use of Multi function, but best I can do. KDD */
    zz->Get_Part_Multi(zz->Get_Part_Multi_Data,
                            zz->Num_GID,zz->Num_LID,1,
                            subroot->global_id,subroot->local_id,
                            &result,ierr);
    if (*ierr != ZOLTAN_OK && *ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from ZOLTAN_PART_MULTI_FN");
    }
  }
  else {

/* otherwise, return my processor number if the object is assigned to
   this processor, or any other value if it is not */

    if (subroot->assigned_to_me)
      result = zz->Proc;
    else
      result = zz->Proc-1;
  }

  return(result);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
