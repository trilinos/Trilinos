/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

#include <stdio.h>
#include "lb_const.h"
#include "reftree_const.h"
#include "all_allo_const.h"
#include "params_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Prototypes for functions internal to this file */

void LB_Reftree_Free_Subtree(LB_REFTREE *subroot);
int any_vert_equals(int vert_match, int *vertices, int start, int finish);
int order_tri_bisect(LB *lb, int *vert1, int *order, int *vertices,
                     int *in_vertex, int *out_vertex, LB_REFTREE *subroot);
int order_other_ref(LB *lb, LB_REFTREE *parent, int num_child, int *num_vert,
                    int *vert1, int *vertices, int *order, int *in_vertex,
                    int *out_vertex);
int order_other_ref_recur(int new_entry, int level, int *order, int *on_path,
                          int num_child, int *has_out, int **share_vert,
                          int *solved);
int find_inout(int level, int num_child, int *num_vert, int *vert1,
               int *vertices, int *in_vertex, int *out_vertex, int *order);
int LB_Reftree_Reinit_Coarse(LB *lb);
static int LB_Reftree_Build_Recursive(LB *lb,LB_REFTREE *subroot);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*  Parameters structure for reftree methods */
static PARAM_VARS REFTREE_params[] = {
        { "REFTREE_HASH_SIZE", NULL, "INT" },
        { NULL, NULL, NULL } };

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Set_Reftree_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */

    status = LB_Check_Param(name, val, REFTREE_params, &result, &index);

    return(status);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Reftree_Init(LB *lb)

{
/*
 *  Function to initialize a refinement tree.  This creates the root and
 *  the first level of the tree, which corresponds to the initial coarse grid
 */
char *yo = "LB_Reftree_Init";
struct LB_reftree_data_struct *reftree_data; /* data pointed to by lb */
LB_REFTREE *root;          /* Root of the refinement tree */
struct LB_reftree_hash_node **hashtab; /* hash table */
int nproc;                 /* number of processors */
LB_GID *local_gids;        /* coarse element Global IDs from user */
LB_LID *local_lids;        /* coarse element Local IDs from user */
LB_GID *all_gids;          /* coarse element Global IDs from all procs */
int *assigned;             /* 1 if the element is assigned to this proc */
int *num_vert;             /* number of vertices for each coarse element */
int *vertices;             /* vertices for the coarse elements */
int *in_vertex;            /* "in" vertex for each coarse element */
int *out_vertex;           /* "out" vertex for each coarse element */
int in_order;              /* 1 if user is supplying order of the elements */
int num_obj;               /* number of coarse objects known to this proc */
int *num_obj_all;          /* num_obj from each processor */
int *displs;               /* running sum of num_obj_all */
int sum_num_obj;           /* full sum of num_obj_all */
int total_num_obj;         /* number of objects in the whole coarse grid */
int ierr;                  /* error flag from calls */
int final_ierr;            /* error flag returned by this routine */
int wdim;                  /* dimension of object weights */
int count;                 /* counter for number of objects */
int sum_vert;              /* summation of number of vertices of objects */
int *order;                /* permutation array for ordering coarse elements */
int found;                 /* flag for terminating first/next query loop */
int hashsize;              /* size of the hash table */
int i, j;                  /* loop counters */
unsigned char *p;          /* for setting IDs to NULL */

  final_ierr = LB_OK;

  if (lb->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = lb->Obj_Weight_Dim;
  }

  nproc = lb->Num_Proc;

  /*
   * Allocate the root of the refinement tree for this load balancing structure.
   * If a tree already exists, destroy it first.
   */

  if (lb->Data_Structure != NULL) LB_Reftree_Free_Structure(lb);

  root = (LB_REFTREE *) LB_MALLOC(sizeof(LB_REFTREE));
  if (root == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
            lb->Proc, yo);
    return(LB_MEMERR);
  }

  /*
   * Initialize the root
   */

  p = (unsigned char *)&(root->global_id);
  for (i=0; i<sizeof(LB_GID); i++,p++) *p = (unsigned char)NULL;
  p = (unsigned char *)&(root->local_id);
  for (i=0; i<sizeof(LB_LID); i++,p++) *p = (unsigned char)NULL;
  root->children       = (LB_REFTREE *) NULL;
  root->num_child      = 0;
  root->num_vertex     = 0;
  root->vertices       = (int *) NULL;
  root->in_vertex      = (int) NULL;
  root->out_vertex     = (int) NULL;
  root->assigned_to_me = 0;
  root->partition      = 0;

  root->weight = (float *) LB_MALLOC(wdim*sizeof(float));
  root->summed_weight = (float *) LB_MALLOC(wdim*sizeof(float));
  root->my_sum_weight = (float *) LB_MALLOC(wdim*sizeof(float));
  if (root->weight == NULL || root->summed_weight == NULL ||
      root->my_sum_weight == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
            lb->Proc, yo);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }
  for (i=0; i<wdim; i++) {
    root->weight[i] = 0.0;
    root->summed_weight[i] = 0.0;
    root->my_sum_weight[i] = 0.0;
  }

  /*
   * Allocate and initialize the hash table.
   */

  REFTREE_params[0].ptr = (void *) &hashsize;
  hashsize = DEFAULT_HASH_TABLE_SIZE;
  LB_Assign_Param_Vals(lb->Params, REFTREE_params, lb->Debug_Level, lb->Proc);

  hashtab = (struct LB_reftree_hash_node **)
            LB_MALLOC(sizeof(struct LB_reftree_hash_node *)*hashsize);
  if (hashtab == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
            lb->Proc, yo);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }
  for (i=0; i<hashsize; i++)
    hashtab[i] = (struct LB_reftree_hash_node *)NULL;

  /*
   * set the lb pointer for later access to the refinement tree and hash table
   */

  reftree_data = (struct LB_reftree_data_struct *)
                 LB_MALLOC(sizeof(struct LB_reftree_data_struct));
  reftree_data->reftree_root = root;
  reftree_data->hash_table = hashtab;
  reftree_data->hash_table_size = hashsize;
  lb->Data_Structure = (void *) reftree_data;

  /*
   * Get the list of initial elements known to this processor
   */


  /*
   * Get the number of objects
   */

  if (lb->Get_Num_Coarse_Obj == NULL) {
    fprintf(stderr, "[%d] Error in %s: Must register LB_NUM_COARSE_OBJ_FN.\n",
            lb->Proc, yo);
    LB_Reftree_Free_Structure(lb);
    return(LB_FATAL);
  }

  num_obj = lb->Get_Num_Coarse_Obj(lb->Get_Num_Coarse_Obj_Data, &ierr);
  if (ierr) {
    fprintf(stderr, "[%d] Error in %s:  Error returned from user function"
                    "Get_Num_Coarse_Obj.\n", lb->Proc, yo);
    LB_Reftree_Free_Structure(lb);
    return(ierr);
  }

  /*
   * Get the objects, if the number is not 0
   */

  if (num_obj > 0) {

    num_obj += 1; /* allocate one extra spot for the last call to NEXT_OBJ */
    local_gids = (LB_GID *) LB_MALLOC(num_obj*sizeof(LB_GID));
    local_lids = (LB_LID *) LB_MALLOC(num_obj*sizeof(LB_LID));
    assigned   = (int *) LB_MALLOC(num_obj*sizeof(int));
    num_vert   = (int *) LB_MALLOC(num_obj*sizeof(int));
    vertices   = (int *) LB_MALLOC(MAXVERT*num_obj*sizeof(int));
    in_vertex  = (int *) LB_MALLOC(num_obj*sizeof(int));
    out_vertex = (int *) LB_MALLOC(num_obj*sizeof(int));
    num_obj -= 1;

    if (local_gids == NULL || local_lids == NULL || assigned  == NULL ||
        num_vert   == NULL || vertices   == NULL || in_vertex == NULL ||
        out_vertex == NULL) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&local_gids);
      LB_FREE(&local_lids);
      LB_FREE(&assigned);
      LB_FREE(&num_vert);
      LB_FREE(&vertices);
      LB_FREE(&in_vertex);
      LB_FREE(&out_vertex);
      LB_Reftree_Free_Structure(lb);
      return(LB_MEMERR);
    }

    if (lb->Get_Coarse_Obj_List != NULL) {

  /*
   * Get objects via list
   */

      lb->Get_Coarse_Obj_List(lb->Get_Coarse_Obj_List_Data, local_gids,
                              local_lids, assigned, num_vert, vertices,
                              &in_order, in_vertex, out_vertex, &ierr);
      if (ierr) {
        fprintf(stderr, "[%d] Error in %s:  Error returned from user function"
                        "Get_Coarse_Obj_List.\n", lb->Proc, yo);
        LB_FREE(&local_gids);
        LB_FREE(&local_lids);
        LB_FREE(&assigned);
        LB_FREE(&num_vert);
        LB_FREE(&vertices);
        LB_FREE(&in_vertex);
        LB_FREE(&out_vertex);
        LB_Reftree_Free_Structure(lb);
        return(ierr);
      }

    }

    else if (lb->Get_First_Coarse_Obj != NULL &&
             lb->Get_Next_Coarse_Obj  != NULL) {

  /*
   * Get objects via first/next
   */

      sum_vert = 0;
      count = 0;
      found = lb->Get_First_Coarse_Obj(lb->Get_First_Coarse_Obj_Data,
                                       &local_gids[count], &local_lids[count],
                                       &assigned[count],
                                       &num_vert[count], &vertices[sum_vert],
                                       &in_order,
                                       &in_vertex[count], &out_vertex[count],
                                       &ierr);
      if (ierr) {
        fprintf(stderr, "[%d] Error in %s:  Error returned from user function"
                        "Get_First_Coarse_Obj.\n", lb->Proc, yo);
        LB_FREE(&local_gids);
        LB_FREE(&local_lids);
        LB_FREE(&assigned);
        LB_FREE(&num_vert);
        LB_FREE(&vertices);
        LB_FREE(&in_vertex);
        LB_FREE(&out_vertex);
        LB_Reftree_Free_Structure(lb);
        return(ierr);
      }

      while (found && count <= num_obj) {
        sum_vert += num_vert[count];
        count += 1;
        found = lb->Get_Next_Coarse_Obj(lb->Get_Next_Coarse_Obj_Data,
                                        &local_gids[count], &local_lids[count],
                                        &assigned[count],
                                        &num_vert[count], &vertices[sum_vert],
                                        &in_vertex[count], &out_vertex[count],
                                        &ierr);
        if (ierr) {
          fprintf(stderr, "[%d] Error in %s:  Error returned from user function"
                          "Get_Next_Coarse_Obj.\n", lb->Proc, yo);
          LB_FREE(&local_gids);
          LB_FREE(&local_lids);
          LB_FREE(&assigned);
          LB_FREE(&num_vert);
          LB_FREE(&vertices);
          LB_FREE(&in_vertex);
          LB_FREE(&out_vertex);
          LB_Reftree_Free_Structure(lb);
          return(ierr);
        }
      }
      if (count != num_obj) {
        fprintf(stderr, "[%d] Warning in %s: Number of objects returned by "
                        "First/Next_Coarse_Obj = %d is not equal to the "
                        "number returned by Num_Coarse_Obj = %d\n",
                        lb->Proc, yo, count, num_obj);
        final_ierr = LB_WARN;
      }
    }

    else {
      fprintf(stderr, "Error in %s:  Must define and register either "
                      "LB_COARSE_OBJ_LIST_FN or "
                      "LB_FIRST_COARSE_OBJ_FN/LB_NEXT_COARSE_OBJ_FN pair\n",
                       yo);
      LB_FREE(&local_gids);
      LB_FREE(&local_lids);
      LB_FREE(&assigned);
      LB_FREE(&num_vert);
      LB_FREE(&vertices);
      LB_FREE(&in_vertex);
      LB_FREE(&out_vertex);
      LB_Reftree_Free_Structure(lb);
      return(LB_FATAL);
    }
  } /* endif (num_obj > 0) */

  /*
   * Communicate to get coarse grid objects unknown to this processor.
   */

  /*
   * First determine how many coarse objects are on each processor
   */

  num_obj_all = (int *)LB_MALLOC(nproc*sizeof(int));
  displs = (int *)LB_MALLOC(nproc*sizeof(int));
  if (num_obj_all == NULL || displs == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n", lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_FREE(&num_obj_all);
    LB_FREE(&displs);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }

  MPI_Allgather((void *)&num_obj,1,MPI_INT,(void *)num_obj_all,1,MPI_INT,
                lb->Communicator);
  displs[0] = 0;
  for (i=1; i<nproc; i++) displs[i] = displs[i-1]+num_obj_all[i-1];
  sum_num_obj = displs[nproc-1] + num_obj_all[nproc-1];

  /*
   * Then get the coarse objects from all processors
   */

  all_gids = (LB_GID *) LB_MALLOC(sum_num_obj*sizeof(LB_GID));
  if (all_gids == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n", lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_FREE(&num_obj_all);
    LB_FREE(&displs);
    LB_FREE(&all_gids);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }

/* TEMP sending this as bytes won't work if the parallel machine consists
        of machines with different representations (e.g. different byte
        ordering) of whatever LB_GID is, but what else can we do? */

  for (i=0; i<nproc; i++) {
    num_obj_all[i] = num_obj_all[i]*sizeof(LB_GID);
    displs[i] = displs[i]*sizeof(LB_GID);
  }

  MPI_Allgatherv((void *)local_gids,num_obj*sizeof(LB_GID),MPI_BYTE,
                 (void *)all_gids,num_obj_all,displs,MPI_BYTE,
                 lb->Communicator);

  LB_FREE(&displs);
  LB_FREE(&num_obj_all);

  /*
   * Finally, build a list with each coarse grid element, beginning with
   * the ones this processor knows.  Also set the default order of the
   * elements as given by the user, with processor rank resolving duplicates
   */

  local_gids = (LB_GID *) LB_REALLOC(local_gids,sum_num_obj*sizeof(LB_GID));
  order = (int *) LB_MALLOC(sum_num_obj*sizeof(int));
  if (local_gids == NULL || order == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n", lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_FREE(&all_gids);
    LB_FREE(&order);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }

/* TEMP this is terribly inefficient.  Probably better to sort all_gids to
        identify the duplicates.  Of course, it's not bad if the
        initial grid is really coarse */

  total_num_obj = num_obj;
  count = 0;
  for (i=0; i<sum_num_obj; i++) order[i] = -1;

  for (i=0; i<sum_num_obj; i++) {
    found = 0;
    for (j=0; j<total_num_obj && !found; j++) {
      if (LB_EQ_GID(all_gids[i],local_gids[j])) found = 1;
    }
    if (found) {
      if (order[j-1] == -1) {
        order[j-1] = count;
        count += 1;
      }
    }
    else {
      LB_SET_GID(local_gids[total_num_obj],all_gids[i]);
      order[total_num_obj] = count;
      count += 1;
      total_num_obj += 1;
    }
  }

  if (count != total_num_obj) {
    fprintf(stderr, "[%d] Warning in %s: Number of objects counted while "
                    "setting default order = %d is not equal to the "
                    "number counted while getting objects from other procs "
                    "= %d\n",lb->Proc, yo, count, total_num_obj);
    final_ierr = LB_WARN;
  }

  LB_FREE(&all_gids);

  num_vert = (int *) LB_REALLOC(num_vert,total_num_obj*sizeof(int));
  if (num_vert == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n", lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }

  for (i=num_obj; i<total_num_obj; i++) num_vert[i] = -1;

  /*
   * Determine the order of the coarse grid elements.
   * If the user supplies the order, it was set above.
   */

  if (!in_order) {

  /*
   * TEMP For now, require that the user provide the order.
   */

    fprintf(stderr, "[%d] Warning in %s: Currently not supporting automatic "
                    "determination of the order of the coarse grid objects.  "
                    "Using the order in which they were provided.\n",
                    lb->Proc, yo);
    final_ierr = LB_WARN;

  }

  /*
   * Copy the elements into the child list of the root
   */

  /*
   * Allocate the children of the root
   */

  root->children = (LB_REFTREE *) LB_MALLOC(total_num_obj*sizeof(LB_REFTREE));
  if (root->children == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
            lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_FREE(&order);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }
  root->num_child = total_num_obj;

  /*
   * Make sure the weights have been provided, if needed
   */

  if (lb->Obj_Weight_Dim != 0 && lb->Get_Child_Weight == NULL) {
    fprintf(stderr, "[%d] Error in %s: Must register LB_CHILD_WEIGHT_FN.\n",
            lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_FREE(&order);
    LB_Reftree_Free_Structure(lb);
    return(LB_FATAL);
  }

  /*
   * For each coarse grid object ...
   */

  sum_vert = 0;
  for (i=0; i<total_num_obj; i++) {

  /*
   * allocate memory within the tree node
   */

    root->children[order[i]].weight = (float *) LB_MALLOC(wdim*sizeof(float));
    root->children[order[i]].summed_weight = (float *) LB_MALLOC(wdim*sizeof(float));
    root->children[order[i]].my_sum_weight = (float *) LB_MALLOC(wdim*sizeof(float));
    if (num_vert[i] <= 0)
      root->children[order[i]].vertices = (int *) LB_MALLOC(sizeof(int));
    else
      root->children[order[i]].vertices = (int *) LB_MALLOC(num_vert[i]*sizeof(int));
    if (root->children[order[i]].weight        == NULL ||
        root->children[order[i]].summed_weight == NULL ||
        root->children[order[i]].my_sum_weight == NULL ||
        root->children[order[i]].vertices      == NULL) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&local_gids);
      LB_FREE(&local_lids);
      LB_FREE(&assigned);
      LB_FREE(&num_vert);
      LB_FREE(&vertices);
      LB_FREE(&in_vertex);
      LB_FREE(&out_vertex);
      LB_FREE(&order);
      LB_Reftree_Free_Structure(lb);
      return(LB_MEMERR);
    }

  /*
   * Get the weight
   */

    if (lb->Obj_Weight_Dim == 0) {
  /* if an initial element is a leaf, the weight gets set to 1 later */
       *(root->children[order[i]].weight) = 0.0;
    }
    else if (num_vert[i] == -1) {
  /* if the element is not known to this processor, the weight is 0 */
       *(root->children[order[i]].weight) = 0.0;
    }
    else {
      lb->Get_Child_Weight(lb->Get_Child_Weight_Data, local_gids[i],
                           local_lids[i], lb->Obj_Weight_Dim, 
                           root->children[order[i]].weight, &ierr);
    }
    for (j=0; j<wdim; j++) {
      root->children[order[i]].summed_weight[j] = 0.0;
      root->children[order[i]].my_sum_weight[j] = 0.0;
    }

  /*
   * Copy the vertices
   */

    for (j=0; j<num_vert[i]; j++) 
      root->children[order[i]].vertices[j] = vertices[sum_vert+j];
    if (num_vert[i] > 0) sum_vert += num_vert[i];

  /*
   * Copy from temporary arrays and set empty defaults
   */

    if (num_vert[i] == -1) {
  /* elements not known to this processor have more empty entries */
      LB_SET_GID(root->children[order[i]].global_id,local_gids[i]);
      p = (unsigned char *)&(root->children[order[i]].local_id);
      for (i=0; i<sizeof(LB_LID); i++,p++) *p = (unsigned char)NULL;
      root->children[order[i]].children       = (LB_REFTREE *) NULL;
      root->children[order[i]].num_child      = 0;
      root->children[order[i]].num_vertex     = num_vert[i];
      root->children[order[i]].in_vertex      = 0;
      root->children[order[i]].out_vertex     = 0;
      root->children[order[i]].assigned_to_me = 0;
      root->children[order[i]].partition      = 0;
    }
    else {
      LB_SET_GID(root->children[order[i]].global_id,local_gids[i]);
      LB_SET_LID(root->children[order[i]].local_id,local_lids[i]);
      root->children[order[i]].children       = (LB_REFTREE *) NULL;
      root->children[order[i]].num_child      = 0;
      root->children[order[i]].num_vertex     = num_vert[i];
      root->children[order[i]].in_vertex      = in_vertex[i];
      root->children[order[i]].out_vertex     = out_vertex[i];
      root->children[order[i]].assigned_to_me = assigned[i];
      root->children[order[i]].partition      = 0;
    }

  /*
   * Add it to the hash table
   */

    LB_Reftree_Hash_Insert(&(root->children[order[i]]),hashtab,hashsize);

  }

  /*
   * clean up and return error code
   */

  LB_FREE(&local_gids);
  LB_FREE(&local_lids);
  LB_FREE(&assigned);
  LB_FREE(&num_vert);
  LB_FREE(&vertices);
  LB_FREE(&in_vertex);
  LB_FREE(&out_vertex);
  LB_FREE(&order);
  return(final_ierr);
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Reftree_Build(LB *lb)

{
/*
 * Function to build a refinement tree
 */
char *yo = "LB_Reftree_Build";
LB_REFTREE *root;          /* Root of the refinement tree */
int ierr;                  /* Error code returned by called functions */
int i;                     /* loop counter */

  /*
   * Initialize the tree, if not already there, and set the root.  If already
   * there, reinitialize coarse grid.
   */

  if (lb->Data_Structure == NULL) {
    ierr = LB_Reftree_Init(lb);
    if (ierr==LB_FATAL || ierr==LB_MEMERR) {
      fprintf(stderr,"[%d] Error in %s:  Error returned from LB_Reftree_Init\n",
                      lb->Proc, yo);
      return(ierr);
    }
  }
  else {
    LB_Reftree_Reinit_Coarse(lb);
  }
  root = ((struct LB_reftree_data_struct *)lb->Data_Structure)->reftree_root;

  /*
   * Verify the required child query functions are registered
   */

  if (lb->Get_Num_Child == NULL || lb->Get_Child_List == NULL) {
    fprintf(stderr, "[%d] Error in %s: Must register LB_NUM_CHILD_FN"
            " and LB_CHILD_LIST_FN\n",lb->Proc, yo);
    LB_Reftree_Free_Structure(lb);
    return(LB_FATAL);
  }

  /*
   * For each child of the root, build the subtree rooted at that child
   * and, if the weights are not provided, set its weight if it is a leaf.
   * Skip elements not known to this processor.
   */

  for (i=0; i<root->num_child; i++) {
    if ( (root->children[i]).num_vertex != -1 ) {
      ierr = LB_Reftree_Build_Recursive(lb,&(root->children[i]));
      if (ierr==LB_FATAL || ierr==LB_MEMERR) {
        fprintf(stderr, "[%d] Error in %s:  Error returned from LB_Reftree"
                        "_Build_Recursive\n",lb->Proc, yo);
        return(ierr);
      }
    }
  }

  return(LB_OK);
}

static int LB_Reftree_Build_Recursive(LB *lb,LB_REFTREE *subroot)

{
/*
 * Recursive function to traverse a tree while building it
 */
static int TEMP_first_warning = 1; /* TEMP until ref_type is fully supported */
char *yo = "LB_Reftree_Build_Recursive";
int ierr;                  /* error code called routines */
int final_ierr;            /* error code returned by this routine */
int num_obj;               /* number of children returned by user */
LB_GID *local_gids;        /* coarse element Global IDs from user */
LB_LID *local_lids;        /* coarse element Local IDs from user */
int *assigned;             /* 1 if the element is assigned to this proc */
int *num_vert;             /* number of vertices for each coarse element */
int *vertices;             /* vertices for the coarse elements */
int *in_vertex;            /* "in" vertex for each coarse element */
int *out_vertex;           /* "out" vertex for each coarse element */
LB_REF_TYPE ref_type;      /* type of refinement that creates children */
int *order;                /* order of the children */
int wdim;                  /* dimension for weights */
int i, j;                  /* loop counters */
int sum_vert;              /* running sum of the number of vertices */
int *vert1;                /* array containing the first vertex for each child*/
struct LB_reftree_hash_node **hashtab; /* hash tree */
int hashsize;              /* size of the hash table */

  final_ierr = LB_OK;
  if (lb->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = lb->Obj_Weight_Dim;
  }

  /*
   * Print a warning if a nonexistent subroot is passed in
   */

  if (subroot == NULL) {
    fprintf(stderr, "[%d] Warning in %s:  Called with nonexistent subroot",
                    lb->Proc, yo);
    return(LB_WARN);
  }

  /*
   * Get the number of children of this node
   */

  num_obj = lb->Get_Num_Child(lb->Get_Num_Child_Data, subroot->global_id,
                              subroot->local_id, &ierr);
  if (ierr) {
    fprintf(stderr, "[%d] Error in %s:  Error returned from user function"
                    "Get_Num_Child.\n", lb->Proc, yo);
    LB_Reftree_Free_Structure(lb);
    return(ierr);
  }

  /*
   * If there are no children, set the weight if it is not user provided,
   * and return.  The default is to use 1.0 for leaves and 0.0 for others.
   */

  if (num_obj == 0) {
    if (lb->Obj_Weight_Dim == 0) *(subroot->weight) = 1.0;
    return(LB_OK);
  }

  /*
   * Get the children
   */

  local_gids = (LB_GID *) LB_MALLOC(num_obj*sizeof(LB_GID));
  local_lids = (LB_LID *) LB_MALLOC(num_obj*sizeof(LB_LID));
  assigned   = (int *) LB_MALLOC(num_obj*sizeof(int));
  num_vert   = (int *) LB_MALLOC(num_obj*sizeof(int));
  vertices   = (int *) LB_MALLOC(MAXVERT*num_obj*sizeof(int));
  in_vertex  = (int *) LB_MALLOC(num_obj*sizeof(int));
  out_vertex = (int *) LB_MALLOC(num_obj*sizeof(int));
  vert1      = (int *) LB_MALLOC((num_obj+1)*sizeof(int));
  if (local_gids == NULL || local_lids == NULL || assigned  == NULL ||
      num_vert   == NULL || vertices   == NULL || in_vertex == NULL ||
      out_vertex == NULL || vert1      == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
            lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_FREE(&vert1);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }
  lb->Get_Child_List(lb->Get_Child_List_Data, subroot->global_id,
                     subroot->local_id, local_gids, local_lids, assigned,
                     num_vert, vertices, &ref_type, in_vertex, out_vertex,
                     &ierr);
  if (ierr) {
    fprintf(stderr, "[%d] Error in %s:  Error returned from user function"
                    "Get_Child_List.\n", lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_FREE(&vert1);
    LB_Reftree_Free_Structure(lb);
    return(ierr);
  }

  /*
   * Set the start of the list of vertices for each child
   */

  vert1[0] = 0;
  for (i=0; i<num_obj; i++) vert1[i+1] = vert1[i] + num_vert[i];

  /*
   * Determine the order of the children
   */

  order = (int *) LB_MALLOC(num_obj*sizeof(int));
  if (order == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
            lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_FREE(&vert1);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }

  /*
   * TEMP until code is supplied to support these refinement types
   */

  switch (ref_type) {
  case LB_QUAD_QUAD:
    if (TEMP_first_warning) {
      fprintf(stderr, "[%d] Warning in %s:  Currently not supporting "
                      "automatic ordering of elements for refinement type "
                      "LB_QUAD_QUAD.  Using LB_OTHER_REF.\n",lb->Proc,yo);
      TEMP_first_warning = 0;
      final_ierr = LB_WARN;
    }
    ref_type = LB_OTHER_REF;
    break;
  case LB_HEX3D_OCT:
    if (TEMP_first_warning) {
      fprintf(stderr, "[%d] Warning in %s:  Currently not supporting "
                      "automatic ordering of elements for refinement type "
                      "LB_HEX3D_OCT.  Using LB_OTHER_REF.\n",lb->Proc,yo);
      TEMP_first_warning = 0;
      final_ierr = LB_WARN;
    }
    ref_type = LB_OTHER_REF;
    break;
  }

  /* end TEMP */

  /*
   *  Refinement type dependent determination of the order of the children
   *  and the in/out vertices
   */

  switch (ref_type) {

  case LB_IN_ORDER:
    for (i=0; i<num_obj; i++) order[i] = i;
    break;
  case LB_TRI_BISECT:
    ierr = order_tri_bisect(lb,vert1,order,vertices,in_vertex,out_vertex,
                            subroot);
    break;
  case LB_QUAD_QUAD:
    /* TEMP */
    printf("Oops, still got into case for QUAD_QUAD\n");
    for (i=0; i<num_obj; i++) order[i] = i;
    break;
  case LB_HEX3D_OCT:
    /* TEMP */
    printf("Oops, still got into case for HEX3D_OCT\n");
    for (i=0; i<num_obj; i++) order[i] = i;
    break;
  case LB_OTHER_REF:
    ierr = order_other_ref(lb, subroot, num_obj, num_vert, vert1, vertices,
                           order, in_vertex, out_vertex);
    break;

  /*
   * Default case if a bad value gets returned; use them in order.
   */
  default:
    fprintf(stderr, "[%d] Warning in %s:  Unknown value returned for ref_type"
            " = %d.  Using children in order provided.\n",lb->Proc,yo,ref_type);
    for (i=0; i<num_obj; i++) order[i] = i;
    final_ierr = LB_WARN;
  }

  /*
   * Copy the children into the child list of the subroot
   */

  /*
   * Allocate the children
   */

  if (subroot->children != NULL) {
    fprintf(stderr, "[%d] Warning from %s: children already existed; memory"
                    " leak potential.\n",lb->Proc, yo);
    final_ierr = LB_WARN;
  }

  subroot->children = (LB_REFTREE *) LB_MALLOC(num_obj*sizeof(LB_REFTREE));
  if (subroot->children == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
            lb->Proc, yo);
    LB_FREE(&local_gids);
    LB_FREE(&local_lids);
    LB_FREE(&assigned);
    LB_FREE(&num_vert);
    LB_FREE(&vertices);
    LB_FREE(&in_vertex);
    LB_FREE(&out_vertex);
    LB_FREE(&order);
    LB_FREE(&vert1);
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }
  subroot->num_child = num_obj;

  hashtab  = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table;
  hashsize = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table_size;

  /*
   * For each child ...
   */

  sum_vert = 0;
  for (i=0; i<num_obj; i++) {

  /*
   * allocate memory within the child
   */

    subroot->children[order[i]].weight = (float *) LB_MALLOC(wdim*sizeof(float));
    subroot->children[order[i]].summed_weight = (float *) LB_MALLOC(wdim*sizeof(float));
    subroot->children[order[i]].my_sum_weight = (float *) LB_MALLOC(wdim*sizeof(float));
    if (num_vert[i] <= 0)
      subroot->children[order[i]].vertices = (int *) LB_MALLOC(sizeof(int));
    else
      subroot->children[order[i]].vertices = (int *) LB_MALLOC(num_vert[i]*sizeof(int));
    if (subroot->children[order[i]].weight        == NULL ||
        subroot->children[order[i]].summed_weight == NULL ||
        subroot->children[order[i]].my_sum_weight == NULL ||
        subroot->children[order[i]].vertices      == NULL) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&local_gids);
      LB_FREE(&local_lids);
      LB_FREE(&assigned);
      LB_FREE(&num_vert);
      LB_FREE(&vertices);
      LB_FREE(&in_vertex);
      LB_FREE(&out_vertex);
      LB_FREE(&vert1);
      LB_FREE(&order);
      LB_Reftree_Free_Structure(lb);
      return(LB_MEMERR);
    }

  /*
   * Get the weight
   */

    if (lb->Obj_Weight_Dim == 0) {
       *(subroot->children[order[i]].weight) = 0.0;
    }
    else {
      lb->Get_Child_Weight(lb->Get_Child_Weight_Data, local_gids[i],
                           local_lids[i], lb->Obj_Weight_Dim,
                           subroot->children[order[i]].weight, &ierr);
    }
    for (j=0; j<wdim; j++) {
      subroot->children[order[i]].summed_weight[j] = 0.0;
      subroot->children[order[i]].my_sum_weight[j] = 0.0;
    }

  /*
   * Copy the vertices
   */

    for (j=0; j<num_vert[i]; j++)
      subroot->children[order[i]].vertices[j] = vertices[sum_vert+j];
    if (num_vert[i] > 0) sum_vert += num_vert[i];

  /*
   * Copy from temporary arrays and set empty defaults
   */

    LB_SET_GID(subroot->children[order[i]].global_id,local_gids[i]);
    LB_SET_LID(subroot->children[order[i]].local_id,local_lids[i]);
    subroot->children[order[i]].children       = (LB_REFTREE *) NULL;
    subroot->children[order[i]].num_child      = 0;
    subroot->children[order[i]].num_vertex     = num_vert[i];
    subroot->children[order[i]].in_vertex      = in_vertex[i];
    subroot->children[order[i]].out_vertex     = out_vertex[i];
    subroot->children[order[i]].assigned_to_me = assigned[i];
    subroot->children[order[i]].partition      = 0;

  /*
   * Add it to the hash table
   */

    LB_Reftree_Hash_Insert(&(subroot->children[order[i]]),hashtab,hashsize);

  }

  /*
   * clean up
   */

  LB_FREE(&local_gids);
  LB_FREE(&local_lids);
  LB_FREE(&assigned);
  LB_FREE(&num_vert);
  LB_FREE(&vertices);
  LB_FREE(&in_vertex);
  LB_FREE(&out_vertex);
  LB_FREE(&vert1);
  LB_FREE(&order);

  /*
   * recursively do the children
   */

  for (i=0; i<subroot->num_child; i++) {
    ierr = LB_Reftree_Build_Recursive(lb,&(subroot->children[i]));
    if (ierr) final_ierr = ierr;
  }

  return(final_ierr);

}

/*****************************************************************************/

int any_vert_equals(int vert_match, int *vertices, int start, int finish)
{
/*
 * function to see if one of the vertices between start and finish is vert_match
 */
int i, result;
  result = 0;
  for (i=start; i<=finish; i++) {
    if (vertices[i] == vert_match) result = 1;
  }
  return(result);
}

/*****************************************************************************/

int order_tri_bisect(LB *lb, int *vert1, int *order, int *vertices,
                     int *in_vertex, int *out_vertex, LB_REFTREE *subroot)
{
/*
 * Function to determine the order of the children and in/out vertices
 * when refinement is done by bisecting triangles.  Determine which of
 * the first and second child has the in_vertex and out_vertex, and find the
 * common vertex to go between them.
 */
char *yo = "order_tri_bisect";
int i, j;                  /* loop indices */
int parents_vert[6];       /* cross index between children and parent */
int parent_in;             /* index of the parent in vertex */
int parent_out;            /* index of the parent out vertex */
int parent_third;          /* index of the parent triangle non in/out vert */
int not_parent[2];         /* index of a vertex the parent does not have */
int has_in[2];             /* flag for an element having the parent in vert */
int has_out[2];            /* flag for an element having the parent out vert */
int has_third[2];          /* flag for a triangle having the parent non in/out*/
int bad_case;              /* flag for failing to identify order */

  /* verify that 3 vertices were given for each triangle; if not, punt */
  if (vert1[1] != 3 || vert1[2] != 6) {
    fprintf(stderr, "[%d] Warning in %s:  Incorrect number of vertices "
                    "given for bisected triangles.\n",lb->Proc,yo);
    order[0] = 0;
    order[1] = 1;
    in_vertex[0] = vertices[vert1[0]];
    out_vertex[0] = vertices[vert1[0]+1];
    in_vertex[1] = vertices[vert1[1]];
    out_vertex[1] = vertices[vert1[1]+1];
    return(LB_WARN);
  }

  /* determine the relationship between the parent's vertices and
     the children's vertices */
  for (i=0; i<6; i++) {
    parents_vert[i] = -1;
    for (j=0; j<3; j++) {
      if (vertices[i] == subroot->vertices[j]) parents_vert[i] = j;
    }
  }

  /* determine the location of the parents in and out vertices */
  parent_in = -1; parent_out = -1; parent_third = -1;
  for (i=0; i<3; i++) {
    if (subroot->vertices[i] == subroot->in_vertex) {
      parent_in = i;
    }
    else if (subroot->vertices[i] == subroot->out_vertex) {
      parent_out = i;
    }
    else {
    parent_third = i;
    }
  }
  if (parent_in == -1 || parent_out == -1 || parent_third == -1) {
    /* failed to locate one of them */
    fprintf(stderr, "[%d] Warning in %s:  Could not locate in and out "
                    "vertices in the parent.\n",lb->Proc,yo);
    order[0] = 0;
    order[1] = 1;
    in_vertex[0] = vertices[vert1[0]];
    out_vertex[0] = vertices[vert1[0]+1];
    in_vertex[1] = vertices[vert1[1]];
    out_vertex[1] = vertices[vert1[1]+1];
    return(LB_WARN);
  }

  /* find the vertex that the parent doesn't have */
  if (parents_vert[0] == -1) not_parent[0] = 0;
  if (parents_vert[1] == -1) not_parent[0] = 1;
  if (parents_vert[2] == -1) not_parent[0] = 2;
  if (parents_vert[3] == -1) not_parent[1] = 3;
  if (parents_vert[4] == -1) not_parent[1] = 4;
  if (parents_vert[5] == -1) not_parent[1] = 5;

  /* see which children have which special vertices */
  if (parents_vert[0] == parent_in || parents_vert[1] == parent_in ||
      parents_vert[2] == parent_in) has_in[0] = 1;
  else has_in[0] = 0;
  if (parents_vert[0] == parent_out || parents_vert[1] == parent_out ||
      parents_vert[2] == parent_out) has_out[0] = 1;
  else has_out[0] = 0;
  if (parents_vert[0] == parent_third || parents_vert[1] == parent_third ||
      parents_vert[2] == parent_third) has_third[0] = 1;
  else has_third[0] = 0;
  if (parents_vert[3] == parent_in || parents_vert[4] == parent_in ||
      parents_vert[5] == parent_in) has_in[1] = 1;
  else has_in[1] = 0;
  if (parents_vert[3] == parent_out || parents_vert[4] == parent_out ||
      parents_vert[5] == parent_out) has_out[1] = 1;
  else has_out[1] = 0;
  if (parents_vert[3] == parent_third || parents_vert[4] == parent_third ||
      parents_vert[5] == parent_third) has_third[1] = 1;
  else has_third[1] = 0;

  /* look for the case for this refinement */
  bad_case = 0;
  if (has_in[0]) {
    if (has_out[1]) {
      order[0] = 0; order[1] = 1;
      in_vertex[0] = subroot->vertices[parent_in];
      out_vertex[1] = subroot->vertices[parent_out];
      if (has_third[0] && has_third[1]) {
        out_vertex[0] = subroot->vertices[parent_third];
        in_vertex[1] = subroot->vertices[parent_third];
      }else{
        out_vertex[0] = vertices[not_parent[0]];
        in_vertex[1] = vertices[not_parent[1]];
      }
    }
    else if (has_in[1]) {
      if (has_out[0]) {
        order[0] = 1; order[1] = 0;
        in_vertex[1] = subroot->vertices[parent_in];
        out_vertex[0] = subroot->vertices[parent_out];
        if (has_third[0] && has_third[1]) {
          out_vertex[1] = subroot->vertices[parent_third];
          in_vertex[0] = subroot->vertices[parent_third];
        }else{
          out_vertex[1] = vertices[not_parent[1]];
          in_vertex[0] = vertices[not_parent[0]];
        }
      }else{ /* impossible case, no one has the out vertex */
        bad_case = 1;
        order[0] = 0; order[1] = 1;
        in_vertex[0] = subroot->vertices[parent_in];
        out_vertex[0] = subroot->vertices[parent_third];
        in_vertex[1] = subroot->vertices[parent_third];
        out_vertex[1] = subroot->vertices[parent_in];
      }
    }else{ /* impossible case, second child has neither in nor out */
      bad_case = 1;
      order[0] = 0; order[1] = 1;
      in_vertex[0] = subroot->vertices[parent_in];
      out_vertex[0] = subroot->vertices[parent_third];
      in_vertex[1] = vertices[3];
      out_vertex[1] = vertices[4];
    }
  }
  else if (has_out[0]) {
    if (has_in[1]) {
      order[0] = 1; order[1] = 0;
      in_vertex[1] = subroot->vertices[parent_in];
      out_vertex[0] = subroot->vertices[parent_out];
      if (has_third[0] && has_third[1]) {
        out_vertex[1] = subroot->vertices[parent_third];
        in_vertex[0] = subroot->vertices[parent_third];
      }else{
        out_vertex[1] = vertices[not_parent[1]];
        in_vertex[0] = vertices[not_parent[0]];
      }
    }else{ /* impossible case, no one has the in vertex */
      bad_case = 1;
      order[0] = 0; order[1] = 1;
      in_vertex[0] = subroot->vertices[parent_out];
      out_vertex[0] = subroot->vertices[parent_third];
      in_vertex[1] = subroot->vertices[parent_third];
      out_vertex[1] = subroot->vertices[parent_out];
    }
  }else{ /* impossible case, first child has neither in nor out */
    bad_case = 1;
    order[0] = 0; order[1] = 1;
    in_vertex[0] = vertices[0];
    out_vertex[0] = vertices[1];
    in_vertex[1] = vertices[3];
    out_vertex[1] = vertices[4];
  }
  if (bad_case) {
    fprintf(stderr, "[%d] Warning in %s:  Vertices of children did not "
                    "match the in and out vertices of parent.\n",lb->Proc,yo);
    return(LB_WARN);
  }
  else {
    return(LB_OK);
  }
}

/*****************************************************************************/

int order_other_ref(LB *lb, LB_REFTREE *parent, int num_child, int *num_vert,
                    int *vert1, int *vertices, int *order, int *in_vertex,
                    int *out_vertex)
{
/*
 * Function to determine the order of the children for an undetermined
 * refinement scheme.  This is expensive as it performs a tree search
 * to solve this NP hard problem, but it should work for any refinement.
 */

char *yo = "order_other_ref";
int i, j, vi, vj;   /* loop counters */
int *has_in;        /* flag for children having in vertex */
int *has_out;       /* flag for children having out vertex */
int **share_vert;   /* flag for two elements sharing a vertex */
int solved;         /* flag for having found the solution */
int final_ierr;     /* error code returned */
int *on_path;       /* flag for already placed element on path */

/* TEMP currently looks for a path in which sequential elements share a
        vertex.  May get a better partition by first searching for a path
        in which they share an edge, and then doing this if that fails */

  final_ierr = LB_OK;

  /*
   * Determine which elements contain the in and out vertices of the parent
   */

  has_in = (int *) LB_MALLOC(num_child*sizeof(int));
  has_out = (int *) LB_MALLOC(num_child*sizeof(int));
  if (has_in == NULL || has_out == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n", lb->Proc, yo);
    LB_FREE(&has_in);
    LB_FREE(&has_out);
    return(LB_MEMERR);
  }

  for (i=0; i<num_child; i++) {
    has_in[i] = 0;
    has_out[i] = 0;
    for (j=0; j<num_vert[i] && !has_in[i]; j++)
      if (vertices[vert1[i]+j] == parent->in_vertex) has_in[i] = 1;
    for (j=0; j<num_vert[i] && !has_out[i]; j++)
      if (vertices[vert1[i]+j] == parent->out_vertex) has_out[i] = 1;
  }

  /*
   * Determine which elements share vertices other than the in/out vertices
   */

  share_vert = (int **) LB_MALLOC(num_child*sizeof(int *));
  if (share_vert == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n", lb->Proc, yo);
    LB_FREE(&share_vert);
    LB_FREE(&has_in);
    LB_FREE(&has_out);
    return(LB_MEMERR);
  }
  for (i=0; i<num_child; i++) {
    share_vert[i] = (int *) LB_MALLOC(num_child*sizeof(int));
    if (share_vert[i] == NULL) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n", lb->Proc, yo);
      for (j=0; j<=i; j++) LB_FREE(&(share_vert[j]));
      LB_FREE(&share_vert);
      LB_FREE(&has_in);
      LB_FREE(&has_out);
      return(LB_MEMERR);
    }
  }

  for (i=0; i<num_child; i++) {
    share_vert[i][i] = 1;
    for (j=i+1; j<num_child; j++) {
      share_vert[i][j] = 0;
      share_vert[j][i] = 0;
      for (vi=0; vi<num_vert[i] && !share_vert[i][j]; vi++) {
        for (vj=0; vj<num_vert[j] && !share_vert[i][j]; vj++) {
          if (vertices[vert1[i]+vi] == vertices[vert1[j]+vj]) {
            if (vertices[vert1[i]+vi] != parent->in_vertex &&
                vertices[vert1[i]+vi] != parent->out_vertex) {
              share_vert[i][j] = 1;
              share_vert[j][i] = 1;
            }
          }
        }
      }
    }
  }

  /*
   * Perform tree search to find solution
   */

  solved = 0;
  on_path = (int *) LB_MALLOC(num_child*sizeof(int));
  if (on_path == NULL) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n", lb->Proc, yo);
    for (j=0; j<=i; j++) LB_FREE(&(share_vert[j]));
    LB_FREE(&on_path);
    LB_FREE(&share_vert);
    LB_FREE(&has_in);
    LB_FREE(&has_out);
    return(LB_MEMERR);
  }
  for (i=0; i<num_child; i++) on_path[i]=0;

  /*
   * Try each element with the in vertex to start the path
   */

  for (i=0; i<num_child && !solved; i++) {
    if (has_in[i]) {
      order_other_ref_recur(i,0,order,on_path,num_child,
                            has_out,share_vert,&solved);
    }
  }

  /*
   * This should have found a solution, but if not then use given order
   */

  if (!solved) {
    fprintf(stderr, "[%d] Warning from %s: Couldn't find path through children."
                    "  Using given order.\n", lb->Proc, yo);
    for (i=0; i<num_child; i++) order[i] = i;
    final_ierr = LB_WARN;
  }

  LB_FREE(&on_path);
  LB_FREE(&share_vert);
  LB_FREE(&has_in);
  LB_FREE(&has_out);

  /*
   * Finally, determine the in and out vertices of each child
   */

  in_vertex[order[0]] = parent->in_vertex;
  out_vertex[order[num_child-1]] = parent->out_vertex;
  solved = find_inout(0, num_child, num_vert, vert1, vertices, in_vertex,
                      out_vertex, order);
  if (!solved) {
    fprintf(stderr, "[%d] Warning from %s: Couldn't find good set of in/out"
                    " vertices.  Using first and second.\n", lb->Proc, yo);
    for (i=0; i<num_child; i++) {
      in_vertex[i]  = vertices[vert1[i]];
      out_vertex[i] = vertices[vert1[i]+1];
    }
    final_ierr = LB_WARN;
  }

  return(final_ierr);
}

/*****************************************************************************/

int order_other_ref_recur(int new_entry, int level, int *order, int *on_path,
                          int num_child, int *has_out, int **share_vert,
                          int *solved)
{
/*
 * Recursive routine to search the solution space tree
 */
int i;

  if (level == num_child-1) {

  /*
   * End of a path, success if this element has the out vertex
   */
    if (has_out[new_entry]) {
      order[level] = new_entry;
      *solved = 1;
    }
    else {
      *solved = 0;
    }
  }

  else {

  /*
   * Add this element to the proposed path
   */

    order[level] = new_entry;
    on_path[new_entry] = 1;

  /*
   * Try each element that is not already on the path and shares a vertex
   * with the current new entry
   */

    for (i=0; i<num_child && !(*solved); i++) {
      if (!on_path[i] && share_vert[new_entry][i]) {
        order_other_ref_recur(i, level+1, order, on_path, num_child, has_out,
                              share_vert, solved);
      }
    }

    on_path[new_entry] = 0;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int find_inout(int level, int num_child, int *num_vert, int *vert1,
               int *vertices, int *in_vertex, int *out_vertex, int *order)
{
/*
 * Function to find in and out vertices.
 * On first call, the first in_vertex and last out_vertex should already be set.
 * level should be 0 in the first call.
 */
int i, j;                       /* loop counters */
int solved;                     /* found a solution */

  if (level == num_child-1) {

  /*
   * Last element.  Success if the in vertex is not the last out
   */

    if (in_vertex[order[level]] == out_vertex[order[level]])
      solved = 0;
    else
      solved = 1;

  }
  else {

  /*
   * Not last element.
   * Try each vertex that is not the in vertex, and if the next element in
   * the path shares that vertex, move on to the next element
   */

    solved = 0;
    for (i=0; i<num_vert[order[level]] && !solved; i++) {
      if (vertices[vert1[order[level]]+i] != in_vertex[order[level]]) {
        for (j=0; j<num_vert[order[level+1]] && !solved; j++) {
          if (vertices[vert1[order[level+1]]+j] == vertices[vert1[order[level]]+i]) {
            out_vertex[order[level]]  = vertices[vert1[order[level]]+i];
            in_vertex[order[level+1]] = vertices[vert1[order[level]]+i];
            solved = find_inout(level+1, num_child, num_vert, vert1, vertices,
                                in_vertex, out_vertex, order);
          }
        }
      }
    }

  }

  return(solved);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Reftree_Free_Structure(LB *lb)

{
/*
 *  Function to free all the memory of a refinement tree
 */
struct LB_reftree_data_struct *reftree_data; /* data structure from lb */
LB_REFTREE *root;                            /* Root of the refinement tree */
struct LB_reftree_hash_node **hashtab;       /* hash table */
int hashsize;                                /* dimension of hash table */
int i;                                       /* loop counter */

  reftree_data = (struct LB_reftree_data_struct *)lb->Data_Structure;

  root = reftree_data->reftree_root;

  if (root != NULL) {

  /*
   * Eliminate all the children, recursively
   */

    if (root->children != NULL) {
      for (i=0; i<root->num_child; i++)
        LB_Reftree_Free_Subtree(&(root->children[i]));

      LB_FREE(&(root->children));
    }

  /*
   * Free the memory in root, and root itself
   */

    LB_FREE(&(root->weight));
    LB_FREE(&(root->summed_weight));
    LB_FREE(&(root->my_sum_weight));
    LB_FREE(&root);
  }

  /*
   * Free the memory in the hash table
   */

  hashtab  = reftree_data->hash_table;
  hashsize = reftree_data->hash_table_size;

  if (hashtab != NULL) {
    LB_Reftree_Clear_Hash_Table(hashtab,hashsize);
    LB_FREE(&hashtab);
  }

  LB_FREE(&(lb->Data_Structure));

}

void LB_Reftree_Free_Subtree(LB_REFTREE *subroot)

{
/*
 *  Function to free the memory of a subtree
 */
int i;   /* loop counter */

  if (subroot != NULL) {

  /*
   * Eliminate all the children, recursively
   */

    if (subroot->children != NULL) {
      for (i=0; i<subroot->num_child; i++) 
        LB_Reftree_Free_Subtree(&(subroot->children[i]));
      
      LB_FREE(&(subroot->children));
    }

  /*
   * Free the memory in subroot
   */

    LB_FREE(&(subroot->weight));
    LB_FREE(&(subroot->summed_weight));
    LB_FREE(&(subroot->my_sum_weight));
    LB_FREE(&(subroot->vertices));
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Reftree_Reinitialize(LB *lb)

{
/*
 *  Function to set a refinement tree back to the initial grid (1 tree level)
 */
struct LB_reftree_data_struct *reftree_data; /* data structure from lb */
LB_REFTREE *root;                            /* Root of the refinement tree */
struct LB_reftree_hash_node **hashtab;       /* hash table */
int hashsize;                                /* dimension of hash table */
LB_REFTREE *child;                           /* A child of the root */
int i, j;                                    /* loop counters */


  reftree_data = (struct LB_reftree_data_struct *)lb->Data_Structure;

  root = reftree_data->reftree_root;
  hashtab  = reftree_data->hash_table;
  hashsize = reftree_data->hash_table_size;

  if (root != NULL) {

  /*
   * Clear out the hash table
   */

    LB_Reftree_Clear_Hash_Table(hashtab,hashsize);

  /*
   * Go through the children of the root (initial elements)
   */

    if (root->children != NULL) {
      for (i=0; i<root->num_child; i++) {
        child = &(root->children[i]);

  /*
   * Eliminate each child of this object, and update this object, and
   * add it back to the hash table
   */
  
        if (child->children != NULL) {
          for (j=0; j<child->num_child; j++)
            LB_Reftree_Free_Subtree(&(child->children[j]));

          LB_FREE(&(child->children));
          child->children = (LB_REFTREE *)NULL;
        }
        child->num_child = 0;
        LB_Reftree_Hash_Insert(child,hashtab,hashsize);
      }
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Reftree_Reinit_Coarse(LB *lb)

{
/*
 *  Function to reestablish which coarse grid elements are known to this proc
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

char *yo = "LB_Reftree_Reinit_Coarse";
LB_REFTREE *root;     /* Root of the refinement tree */
struct LB_reftree_hash_node **hashtab; /* hash table */
int hashsize;         /* dimension of hash table */
int i, j;             /* loop counter */
LB_GID *local_gids;   /* coarse element Global IDs from user */
LB_LID *local_lids;   /* coarse element Local IDs from user */
int *assigned;        /* 1 if the element is assigned to this proc */
int *num_vert;        /* number of vertices for each coarse element */
int *vertices;        /* vertices for the coarse elements */
int *in_vertex;       /* "in" vertex for each coarse element */
int *out_vertex;      /* "out" vertex for each coarse element */
LB_GID slocal_gids;   /* coarse element Global IDs from user */
LB_LID slocal_lids;   /* coarse element Local IDs from user */
int sassigned;        /* 1 if the element is assigned to this proc */
int snum_vert;        /* number of vertices for a coarse element */
int sin_vertex;       /* "in" vertex for a coarse element */
int sout_vertex;      /* "out" vertex for a coarse element */
int in_order;         /* 1 if user is supplying order of the elements */
int num_obj;          /* number of coarse objects known to this proc */
int ierr;             /* error flag */
LB_REFTREE *tree_node;/* pointer to an initial grid element in the tree */
int final_ierr;       /* error code returned */
int sum_vert;         /* running total of number of vertices */
int found;            /* flag for another coarse grid element */

  root = ((struct LB_reftree_data_struct *)lb->Data_Structure)->reftree_root;
  hashtab  = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table;
  hashsize = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table_size;
  final_ierr = LB_OK;

  /*
   * Mark all coarse elements as unknown
   */

  for (i=0; i<root->num_child; i++) {
    ((root->children)[i]).num_vertex = -1;
  }

  /*
   * Get the coarse grid objects and update the vertices, whether the element
   * is assigned to this processor, weight and, if not already set, in/out vert
   */

  if (lb->Get_Coarse_Obj_List != NULL) {

  /*
   * Get objects via list
   */

    num_obj = lb->Get_Num_Coarse_Obj(lb->Get_Num_Coarse_Obj_Data, &ierr);
    if (ierr) {
      fprintf(stderr, "[%d] Error in %s:  Error returned from user function"
                      "Get_Num_Coarse_Obj.\n", lb->Proc, yo);
      return(ierr);
    }

    if (num_obj > 0) {
      local_gids = (LB_GID *) LB_MALLOC(num_obj*sizeof(LB_GID));
      local_lids = (LB_LID *) LB_MALLOC(num_obj*sizeof(LB_LID));
      assigned   = (int *) LB_MALLOC(num_obj*sizeof(int));
      num_vert   = (int *) LB_MALLOC(num_obj*sizeof(int));
      vertices   = (int *) LB_MALLOC(MAXVERT*num_obj*sizeof(int));
      in_vertex  = (int *) LB_MALLOC(num_obj*sizeof(int));
      out_vertex = (int *) LB_MALLOC(num_obj*sizeof(int));

      if (local_gids == NULL || local_lids == NULL || assigned  == NULL ||
          num_vert   == NULL || vertices   == NULL || in_vertex == NULL ||
          out_vertex == NULL) {
        fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
                lb->Proc, yo);
        LB_FREE(&local_gids);
        LB_FREE(&local_lids);
        LB_FREE(&assigned);
        LB_FREE(&num_vert);
        LB_FREE(&vertices);
        LB_FREE(&in_vertex);
        LB_FREE(&out_vertex);
        return(LB_MEMERR);
      }

      lb->Get_Coarse_Obj_List(lb->Get_Coarse_Obj_List_Data, local_gids,
                              local_lids, assigned, num_vert, vertices,
                              &in_order, in_vertex, out_vertex, &ierr);
      if (ierr) {
        fprintf(stderr, "[%d] Error in %s:  Error returned from user function"
                        "Get_Coarse_Obj_List.\n", lb->Proc, yo);
        LB_FREE(&local_gids);
        LB_FREE(&local_lids);
        LB_FREE(&assigned);
        LB_FREE(&num_vert);
        LB_FREE(&vertices);
        LB_FREE(&in_vertex);
        LB_FREE(&out_vertex);
        return(ierr);
      }

      sum_vert = 0;
      for (i=0; i<num_obj; i++) {

        tree_node = LB_Reftree_hash_lookup(hashtab,local_gids[i],hashsize);
        if (tree_node == NULL) {
          fprintf(stderr, "[%d] Warning in %s: coarse grid element not"
                          " previously seen.\n", lb->Proc, yo);
          final_ierr = LB_WARN;
        }
        else {
          tree_node->num_vertex = num_vert[i];
          LB_FREE(&(tree_node->vertices));
          if (num_vert[i] <= 0)
            tree_node->vertices = (int *) LB_MALLOC(sizeof(int));
          else
            tree_node->vertices = (int *) LB_MALLOC(num_vert[i]*sizeof(int));
          if (tree_node->vertices == NULL) {
            fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
                    lb->Proc, yo);
            LB_FREE(&local_gids);
            LB_FREE(&local_lids);
            LB_FREE(&assigned);
            LB_FREE(&num_vert);
            LB_FREE(&vertices);
            LB_FREE(&in_vertex);
            LB_FREE(&out_vertex);
            return(LB_MEMERR);
          }

          for (j=0; j<num_vert[i]; j++)
            tree_node->vertices[j] = vertices[sum_vert+j];
          if (num_vert[i] > 0) sum_vert += num_vert[i];

          tree_node->assigned_to_me = assigned[i];
/* TEMP if not provided in_order, then in/out are not returned and must be
        determined */
          if (tree_node->in_vertex == 0) tree_node->in_vertex = in_vertex[i];
          if (tree_node->out_vertex == 0) tree_node->out_vertex = out_vertex[i];
          lb->Get_Child_Weight(lb->Get_Child_Weight_Data, local_gids[i],
                               local_lids[i], lb->Obj_Weight_Dim,
                               tree_node->weight, &ierr);
        }
      }
      LB_FREE(&local_gids);
      LB_FREE(&local_lids);
      LB_FREE(&assigned);
      LB_FREE(&num_vert);
      LB_FREE(&vertices);
      LB_FREE(&in_vertex);
      LB_FREE(&out_vertex);
    }

  }
  else {

  /*
   * Get objects via first/next
   */

    vertices = (int *) LB_MALLOC(MAXVERT*sizeof(int));
    if (tree_node->vertices == NULL) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      return(LB_MEMERR);
    }

    found = lb->Get_First_Coarse_Obj(lb->Get_First_Coarse_Obj_Data,
                                     &slocal_gids, &slocal_lids, &sassigned,
                                     &snum_vert, vertices, &in_order,
                                     &sin_vertex, &sout_vertex, &ierr);
    if (ierr) {
      fprintf(stderr, "[%d] Error in %s:  Error returned from user function"
                      "Get_First_Coarse_Obj.\n", lb->Proc, yo);
      LB_FREE(&vertices);
      return(ierr);
    }
    while (found) {
      tree_node = LB_Reftree_hash_lookup(hashtab,slocal_gids,hashsize);
      if (tree_node == NULL) {
        fprintf(stderr, "[%d] Warning in %s: coarse grid element not"
                        " previously seen.\n", lb->Proc, yo);
        final_ierr = LB_WARN;
      }
      else {
        tree_node->num_vertex = snum_vert;
        LB_FREE(&(tree_node->vertices));
        if (snum_vert <= 0)
          tree_node->vertices = (int *) LB_MALLOC(sizeof(int));
        else
          tree_node->vertices = (int *) LB_MALLOC(snum_vert*sizeof(int));
        if (tree_node->vertices == NULL) {
          fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
                  lb->Proc, yo);
          return(LB_MEMERR);
        }
  
        for (j=0; j<snum_vert; j++)
          tree_node->vertices[j] = vertices[j];
  
        tree_node->assigned_to_me = sassigned;
/* TEMP if not provided in_order, then in/out are not returned and must be
        determined */
        if (tree_node->in_vertex == 0) tree_node->in_vertex = sin_vertex;
        if (tree_node->out_vertex == 0) tree_node->out_vertex = sout_vertex;
        lb->Get_Child_Weight(lb->Get_Child_Weight_Data, slocal_gids,
                             slocal_lids, lb->Obj_Weight_Dim,
                             tree_node->weight, &ierr);
      }

      found = lb->Get_Next_Coarse_Obj(lb->Get_Next_Coarse_Obj_Data,
                                       &slocal_gids, &slocal_lids, &sassigned,
                                       &snum_vert, vertices,
                                       &sin_vertex, &sout_vertex, &ierr);
    }
  }
  return(final_ierr);
}

void LB_Reftree_Print(LB *lb, LB_REFTREE *subroot, int level)
{
/*
 * Print the refinement tree, for debugging
 */

  int i, me;
  int *p;

  if (subroot == NULL) return;

  me = lb->Proc;
  printf("\n");
  printf("[%d] refinement tree node with local id %d on level %d\n",
         me,subroot->local_id,level);
  printf("[%d]   Global ID ",me);
  p = (int *)(&(subroot->global_id));
  for (i=0; i<sizeof(LB_GID)/4; i++,p++) printf("%d ",*p);
  printf("\n");
  printf("[%d]   first weight %f\n",me,subroot->weight[0]);
  printf("[%d]   first summed weight %f\n",me,subroot->summed_weight[0]);
  printf("[%d]   first my_sum weight %f\n",me,subroot->my_sum_weight[0]);
  printf("[%d]   number of vertices %d\n",me,subroot->num_vertex);
  printf("[%d]   vertices ",me);
  for (i=0; i<subroot->num_vertex; i++) printf("%d ",subroot->vertices[i]);
  printf("\n");
  printf("[%d]   in and out vertices %d %d\n",me,subroot->in_vertex,subroot->out_vertex);
  printf("[%d]   assigned_to_me %d\n",me,subroot->assigned_to_me);
  printf("[%d]   partition %d\n",me,subroot->partition);
  printf("[%d]   number of children %d \n",me,subroot->num_child);
  printf("[%d]   children follow.\n",me);
  for (i=0; i<subroot->num_child; i++)
    LB_Reftree_Print(lb,&(subroot->children[i]),level+1);
}
