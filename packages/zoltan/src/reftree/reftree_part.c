#include <stdio.h>
#include "lb_const.h"
#include "reftree_const.h"
#include "all_allo_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Prototypes for functions internal to this file */

void LB_Reftree_Sum_My_Weights(LB *lb, LB_REFTREE *subroot, int *count, int wdim);
void LB_Reftree_Sum_All_Weights(LB *lb, LB_REFTREE *subroot, int wdim);
void LB_Reftree_List_Other_Leaves(LB_REFTREE *subroot, LB_GID *list, int *count);
int LB_Reftree_Partition(LB *lb, int *num_export, LB_GID **export_global_ids,
                         LB_LID **export_local_ids, int **export_procs);
void LB_Reftree_Part_Recursive(LB *lb, LB_REFTREE *subroot, int *part,
                               float *current_size, int *num_exp, float *cutoff,
                               float partition_size);
void LB_Reftree_Mark_and_Count(LB_REFTREE *subroot, int part, int *num_exp);
int LB_Reftree_Export_Lists(LB *lb, LB_REFTREE *subroot, int *num_export,
                            LB_GID **export_global_ids,
                            LB_LID **export_local_ids, int **export_procs);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Reftree_Part(

  LB *lb,                     /* The load-balancing structure */
  int *num_import,            /* Not computed, set to -1 */
  LB_GID **import_global_ids, /* Not computed */
  LB_LID **import_local_ids,  /* Not computed */
  int **import_procs,         /* Not computed */
  int *num_export,            /* Number of objects to be exported */
  LB_GID **export_global_ids, /* global ids of objects to be exported */
  LB_LID **export_local_ids,  /* local  ids of objects to be exported */
  int **export_procs          /* list of processors to export to */
)
{
char *yo = "LB_Reftree_Part";
int ierr;       /* error code returned by called routines */
int final_ierr; /* error code returned by this routine */

  *num_import = -1;
  final_ierr = LB_OK;

  /*
   * initialize the tree (first call only)
   */

  if (lb->Data_Structure == NULL) {
    ierr = LB_Reftree_Init(lb);
    if (ierr==LB_FATAL || ierr==LB_MEMERR) {
      fprintf(stderr, "[%d] Error from %s: Error returned by "
              "LB_Reftree_Init\n",lb->Proc,yo);
      return(ierr);
    }
    if (ierr==LB_WARN) final_ierr = LB_WARN;
  }

  /*
   * build the refinement tree
   */

  ierr = LB_Reftree_Build(lb);
  if (ierr==LB_FATAL || ierr==LB_MEMERR) {
    fprintf(stderr, "[%d] Error from %s: Error returned by "
            "LB_Reftree_Build\n",lb->Proc,yo);
    return(ierr);
  }
  if (ierr==LB_WARN) final_ierr = LB_WARN;

  /*
   * sum the weights in the tree
   */

  ierr = LB_Reftree_Sum_Weights(lb);
  if (ierr==LB_FATAL || ierr==LB_MEMERR) {
    fprintf(stderr, "[%d] Error from %s: Error returned by "
            "LB_Reftree_Sum_Weights\n",lb->Proc,yo);
    return(ierr);
  }
  if (ierr==LB_WARN) final_ierr = LB_WARN;

  /*
   * determine the new partition
   */

  ierr = LB_Reftree_Partition(lb, num_export, export_global_ids,
                              export_local_ids, export_procs);
  if (ierr==LB_FATAL || ierr==LB_MEMERR) {
    fprintf(stderr, "[%d] Error from %s: Error returned by "
            "LB_Reftree_Partition\n",lb->Proc,yo);
    return(ierr);
  }
  if (ierr==LB_WARN) final_ierr = LB_WARN;

  /*
   * delete the tree, except for the first level (initial coarse grid)
   */

  LB_Reftree_Reinitialize(lb);

  return(final_ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Reftree_Sum_Weights(LB *lb)

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
 *
 * TEMP This version uses all-to-all communciation.  Need to find a less
 *      communication intensive approach.
 * Note to myself: The log(p) pairwise approach didn't work because
 * intermediate processors didn't pass on partial sums for the nodes
 * they didn't have.  I believe this can be fixed by added a second
 * log(p) pairwise communication phase to build a list of all processors'
 * requests, so that intermediate processors know to keep the partial
 * sums that they don't need.
 */
char *yo = "LB_Reftree_Sum_Weights";
LB_REFTREE *root;         /* Root of the refinement tree */
int wdim;                 /* Dimension of the weight array */
int i,j,comm_loop,rproc;  /* loop counters */
int count;                /* counter */
LB_GID *leaf_list;        /* list of the leaves not assigned to this proc */
int nproc;                /* number of processors */
int myproc;               /* this processor's processor number */
MPI_Status mess_status;   /* status of an MPI message */
int reqsize;              /* number of leaves for which weights are requested */
LB_GID *recv_gid;         /* received message of LB_GIDs */
LB_REFTREE *node;         /* a node in the refinement tree */
struct LB_reftree_hash_node **hashtab; /* hash tree */
int hashsize;             /* dimension of hash tree */
float *send_float;        /* sending message of floats */
float *recv_float;        /* received message of floats */

  /*
   * set the root and hash table
   */

  root = ((struct LB_reftree_data_struct *)lb->Data_Structure)->reftree_root;
  if (root == NULL) {
    fprintf(stderr, "[%d] Error from %s: Refinement tree not defined.\n",
            lb->Proc, yo);
    return(LB_FATAL);
  }
  hashtab  = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table;
  hashsize = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table_size;

  /*
   * Determine the dimension of the weight array
   */

  if (lb->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = lb->Obj_Weight_Dim;
  }

  /*
   * In the first pass, sum the weights of the nodes that are assigned to
   * this processor, and count the leaves that are not.
   */

  count = 0;
  for (i=0; i<root->num_child; i++) {
    LB_Reftree_Sum_My_Weights(lb,&(root->children[i]),&count,wdim);
  }
  root->assigned_to_me = -1;

  /*
   * Make a list of the leaves that are not assigned to this processor
   */

  if (count == 0)
    leaf_list = (LB_GID *) LB_MALLOC(sizeof(LB_GID));
  else
    leaf_list = (LB_GID *) LB_MALLOC(count*sizeof(LB_GID));
  if (leaf_list == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
            lb->Proc,yo);
    return(LB_MEMERR);
  }

  count = 0;
  LB_Reftree_List_Other_Leaves(root,leaf_list,&count);

  /*
   * Get the unknown leaf weights from other processors.
   */

  nproc = lb->Num_Proc;
  myproc = lb->Proc;

  /*
   * Communication loop
   */

  for (comm_loop=0; comm_loop<nproc; comm_loop++) {

  /*
   * Broadcast the size of the request
   */

    reqsize = count;
    MPI_Bcast((void *)&reqsize, 1, MPI_INT, comm_loop, lb->Communicator);

  /*
   * If the request size is 0, skip this communication
   */

    if (reqsize != 0) {

  /*
   * Allocate space to receive the request
   */

      if (myproc != comm_loop) {
        recv_gid = (LB_GID *)LB_MALLOC(reqsize*sizeof(LB_GID));
        if (recv_gid == NULL) {
          fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
                  lb->Proc,yo);
          LB_FREE(&leaf_list);
          return(LB_MEMERR);
        }
      }

  /*
   * broadcast the list of leaves processor comm_loop needs to know
   */

/* TEMP sending this as bytes won't work if the parallel machine consists
        of machines with different representations (e.g. different byte
        ordering) of whatever LB_GID is, but what else can we do? */

      if (myproc == comm_loop)
        MPI_Bcast((void *)leaf_list,reqsize*sizeof(LB_GID),MPI_BYTE,comm_loop,
                  lb->Communicator);
      else
        MPI_Bcast((void *)recv_gid,reqsize*sizeof(LB_GID),MPI_BYTE,comm_loop,
                  lb->Communicator);

  /*
   *  Reply with any weights I have, and 0. where I don't have it
   */

      if (myproc != comm_loop) {

        send_float = (float *) LB_MALLOC(sizeof(float)*wdim*reqsize);
        if (send_float == NULL) {
          fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
                  lb->Proc,yo);
          LB_FREE(&recv_gid);
          LB_FREE(&leaf_list);
          return(LB_MEMERR);
        }

        for (i=0; i<reqsize; i++) {
          node = LB_Reftree_hash_lookup(hashtab,recv_gid[i],hashsize);
          if (node == NULL)
             for (j=0; j<wdim; j++) send_float[i*wdim+j] = 0.0;
          else
             for (j=0; j<wdim; j++) send_float[i*wdim+j] = node->my_sum_weight[j];
        }

        MPI_Send((void *)send_float,wdim*reqsize,MPI_FLOAT,comm_loop,
                 100+comm_loop,lb->Communicator);
        LB_FREE(&send_float);
        LB_FREE(&recv_gid);
      }

  /*
   * Receive the weights and add them into the tree
   */

      if (myproc == comm_loop) {
        recv_float = (float *)LB_MALLOC(sizeof(float)*wdim*reqsize);
        if (recv_float == NULL) {
          fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
                  lb->Proc,yo);
          LB_FREE(&leaf_list);
          return(LB_MEMERR);
        }

        for (rproc=1; rproc<nproc; rproc++) {
          MPI_Recv((void *)recv_float,reqsize*wdim,MPI_FLOAT,MPI_ANY_SOURCE,
                   100+comm_loop, lb->Communicator,&mess_status);

          for (i=0; i<reqsize; i++) {
            node = LB_Reftree_hash_lookup(hashtab,leaf_list[i],hashsize);
            for (j=0; j<wdim; j++) node->summed_weight[j] += recv_float[i];
          }
        }
        LB_FREE(&recv_float);
      }
    }

  } /* end comm_loop */

  LB_FREE(&leaf_list);

  /*
   * All the leaves now have summed_weight set.
   * Sum the weights throughout the tree.
   */

  LB_Reftree_Sum_All_Weights(lb,root,wdim);

  return(LB_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Reftree_Sum_My_Weights(LB *lb, LB_REFTREE *subroot, int *count, int wdim)

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
   * If there are children, sum the weights of the children and set
   * assigned to me.
   * TEMP Not using the weights of interior nodes because I haven't
   *      determined how to handle interior nodes that are assigned to
   *      other processors.
   */

    for (i=0; i<wdim; i++) subroot->my_sum_weight[i] = 0.0;
    none_assigned = 1;
    all_assigned = 1;

    for (j=0; j<subroot->num_child; j++) {
      LB_Reftree_Sum_My_Weights(lb,&(subroot->children[j]),count,wdim);
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

void LB_Reftree_Sum_All_Weights(LB *lb, LB_REFTREE *subroot, int wdim)

{
/*
 * Function to recursively sum the weights of all the nodes in the tree
 * assuming that summed_weight contains partial sums
 */
int i, j;   /* loop counter */

  /*
   * If there are no children, do nothing
   */

  if (subroot->num_child != 0) {

  /*
   * If there are children, sum the weights of the children.
   * TEMP Not using the weights of interior nodes because I haven't
   *      determined how to handle interior nodes that are assigned to
   *      other processors.
   */

    for (i=0; i<wdim; i++) subroot->summed_weight[i] = 0.0;

    for (j=0; j<subroot->num_child; j++) {
      LB_Reftree_Sum_All_Weights(lb,&(subroot->children[j]),wdim);
      for (i=0; i<wdim; i++)
        subroot->summed_weight[i] += (subroot->children[j]).summed_weight[i];
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Reftree_List_Other_Leaves(LB_REFTREE *subroot, LB_GID *list, int *count)

{
/*
 * Function to make a list of the leaves not assigned to this processor
 */
int i, j;   /* loop counter */

  if (subroot->num_child == 0) {

  /*
   * If there are no children, then add it to the list if it is not
   * assigned to this processor
   */

    if (!subroot->assigned_to_me) {
      LB_SET_GID(list[*count],subroot->global_id);
      *count += 1;
    }

  }
  else {

  /*
   * If there are children, traverse the subtrees
   */

    for (j=0; j<subroot->num_child; j++) {
      LB_Reftree_List_Other_Leaves(&(subroot->children[j]),list,count);
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Reftree_Partition(LB *lb, int *num_export, LB_GID **export_global_ids,
                         LB_LID **export_local_ids, int **export_procs)

{
/*
 * Function to partition the grid and fill the export arrays.
 */

char *yo = "LB_Reftree_Partition";
int num_exp;          /* count the number of export objects */
LB_REFTREE *root;     /* root of the tree */
float partition_size; /* amount of weight for each partition */
float cutoff;         /* the value at which the current partition is full */
int part;             /* partition under construction */
float current_size;   /* amount of weight consumed so far */

/* TEMP don't know what to do with multicomponent weights.  Just using
        the first component */

/* TEMP don't know how to allow a tolerance in the partition sizes */

  root = ((struct LB_reftree_data_struct *)lb->Data_Structure)->reftree_root;

  /*
   * determine the size of the partitions
   */

  partition_size = root->summed_weight[0]/lb->Num_Proc;

  /*
   * traverse the tree to define the partition and count the number of exports
   */

  num_exp = 0;
  part = 0;
  current_size = 0.0;
  cutoff = partition_size;
  LB_Reftree_Part_Recursive(lb,root,&part,&current_size,&num_exp,&cutoff,
                            partition_size);

  /*
   * if no exports, we're done
   */

  if (num_exp == 0) {
    *num_export = num_exp;
    return(LB_OK);
  }

  /*
   * allocate space for the export lists
   */

  if (!LB_Special_Malloc(lb,(void **)export_global_ids,num_exp,
                         LB_SPECIAL_MALLOC_GID))
    return LB_MEMERR;
  if (!LB_Special_Malloc(lb,(void **)export_local_ids,num_exp,
                         LB_SPECIAL_MALLOC_LID)) {
    LB_Special_Free(lb,(void **)export_global_ids,LB_SPECIAL_MALLOC_GID);
    return LB_MEMERR;
  }
  if (!LB_Special_Malloc(lb,(void **)export_procs,num_exp,
                         LB_SPECIAL_MALLOC_INT)) {
    LB_Special_Free(lb,(void **)export_global_ids,LB_SPECIAL_MALLOC_GID);
    LB_Special_Free(lb,(void **)export_local_ids,LB_SPECIAL_MALLOC_LID);
    return LB_MEMERR;
  }

  /*
   * traverse the tree to make the export lists
   */

  *num_export = 0;
  LB_Reftree_Export_Lists(lb,root,num_export,export_global_ids,
                          export_local_ids,export_procs);

  if (num_exp != *num_export) {
    fprintf(stderr, "[%d] Warning from %s: num_exp = %d not equal to "
            "num_export = %d\n",lb->Proc,yo,num_exp,*num_export);
    return(LB_WARN);
  }

  return(LB_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Reftree_Part_Recursive(LB *lb, LB_REFTREE *subroot, int *part,
                              float *current_size, int *num_exp, float *cutoff,
                              float partition_size)

{
/*
 * function to recursively define the partition and count the number of exports
 */
int i;        /* loop counter */

  if (*current_size + subroot->summed_weight[0] <= *cutoff) {

  /*
   * This subtree fits in the current partition
   */

    subroot->partition = *part;
    *current_size += subroot->summed_weight[0];

  /*
   * If this is this processor's partition, there are no exports below this
   * node, so we don't have to traverse the subtree.
   * If there are no leaves of this subtree assigned to this processor, there
   * are no exports below this node.
   * Otherwise, traverse the subtree setting partition and counting exports
   */

    if (*part != lb->Proc && subroot->assigned_to_me) {
      if (subroot->num_child == 0)
        *num_exp += 1;
      else
        for (i=0; i<subroot->num_child; i++)
          LB_Reftree_Mark_and_Count(&(subroot->children[i]), *part, num_exp);
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
        LB_Reftree_Part_Recursive(lb, &(subroot->children[i]), part,
                                 current_size, num_exp, cutoff, partition_size);
    }
    else {

  /*
   * If there are no children, move on to the next partition
   */

      while (*current_size + subroot->summed_weight[0] > *cutoff) {
        *part += 1;
        *cutoff = (*part + 1)*partition_size;
      }
      subroot->partition = *part;
      *current_size += subroot->summed_weight[0];
      if (*part != lb->Proc && subroot->assigned_to_me) *num_exp += 1;
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Reftree_Mark_and_Count(LB_REFTREE *subroot, int part, int *num_exp)

{
/*
 * Function to set the partition and count exports
 */
int i;

  subroot->partition = part;
  if (subroot->num_child == 0) {
    if (subroot->assigned_to_me) *num_exp += 1;
  } else {
    for (i=0; i<subroot->num_child; i++)
      LB_Reftree_Mark_and_Count(&(subroot->children[i]), part, num_exp);
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Reftree_Export_Lists(LB *lb, LB_REFTREE *subroot, int *num_export,
                            LB_GID **export_global_ids,
                            LB_LID **export_local_ids, int **export_procs)
{
/*
 * Function to build the export lists
 */
int i;

/*
 * if this subtree has no leaves assigned to this processor, or if the
 * partition of this subtree is this processor, then there can be no
 * exports below it
 */

  if (!subroot->assigned_to_me || subroot->partition == lb->Proc) return;

  if (subroot->num_child == 0) {

/*
 * if this is a leaf, put it on the export lists
 */

    LB_SET_GID((*export_global_ids)[*num_export],subroot->global_id);
    LB_SET_LID((*export_local_ids)[*num_export],subroot->local_id);
    (*export_procs)[*num_export] = subroot->partition;
    *num_export += 1;

  }
  else {

/*
 * if it is not a leaf, traverse the subtree
 */

    for (i=0; i<subroot->num_child; i++)
      LB_Reftree_Export_Lists(lb, &(subroot->children[i]), num_export,
                              export_global_ids,export_local_ids,export_procs);
  }
}
