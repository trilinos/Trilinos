#include <stdio.h>
#include "lb_const.h"
#include "reftree_const.h"
#include "all_allo_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Prototypes for functions internal to this file */

void LB_Reftree_Free_Subtree(LB_REFTREE *subroot);
int any_vert_equals(int vert_match, int *vertices, int start, int finish);

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
LB_GID *local_gids;        /* coarse element Global IDs from user */
LB_LID *local_lids;        /* coarse element Local IDs from user */
int *assigned;             /* 1 if the element is assigned to this proc */
int *num_vert;             /* number of vertices for each coarse element */
int *vertices;             /* vertices for the coarse elements */
int *in_vertex;            /* "in" vertex for each coarse element */
int *out_vertex;           /* "out" vertex for each coarse element */
int in_order;              /* 1 if user is supplying order of the elements */
int num_obj;               /* number of coarse objects known to this proc */
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

  /*
   * Allocate the root of the refinement tree for this load balancing object.
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

  hashsize = HASH_TABLE_SIZE;
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
        found = lb->Get_Next_Coarse_Obj(lb->Get_First_Coarse_Obj_Data,
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
TEMP For now, assume all coarse grid objects are known to all processors.
*/

  total_num_obj = num_obj;
  /* TEMP in what follows, most local_gids etc will be replaced by whatever
          I call them after doing this */

  /*
   * Determine the order of the coarse grid elements.
   */

  order = (int *) LB_MALLOC(total_num_obj*sizeof(int));
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
    LB_Reftree_Free_Structure(lb);
    return(LB_MEMERR);
  }

  if (in_order) {

  /*
   * User supplied order; permutation matrix is the identity
   */

    for (i=0; i<total_num_obj; i++) order[i] = i;
  }

  else {

  /*
   * TEMP For now, require that the user provide the order.
   */

    for (i=0; i<total_num_obj; i++) order[i] = i;

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
    sum_vert += num_vert[i];

  /*
   * Copy from temporary arrays and set empty defaults
   */

    root->children[order[i]].global_id      = local_gids[i];
    root->children[order[i]].local_id       = local_lids[i];
    root->children[order[i]].children       = (LB_REFTREE *) NULL;
    root->children[order[i]].num_child      = 0;
    root->children[order[i]].num_vertex     = num_vert[i];
    root->children[order[i]].in_vertex      = in_vertex[i];
    root->children[order[i]].out_vertex     = out_vertex[i];
    root->children[order[i]].assigned_to_me = assigned[i];
    root->children[order[i]].partition      = 0;

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
   * Initialize the tree, if not already there, and set the root
   */

  if (lb->Data_Structure == NULL) {
    ierr = LB_Reftree_Init(lb);
    if (ierr==LB_FATAL || ierr==LB_MEMERR) {
      fprintf(stderr,"[%d] Error in %s:  Error returned from LB_Reftree_Init\n",
                      lb->Proc, yo);
      return(ierr);
    }
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
   * and, if the weights are not provided, set its weight if it is a leaf
   */

  for (i=0; i<root->num_child; i++) {
    ierr = LB_Reftree_Build_Recursive(lb,&(root->children[i]));
    if (ierr==LB_FATAL || ierr==LB_MEMERR) {
      fprintf(stderr, "[%d] Error in %s:  Error returned from LB_Reftree"
                      "_Build_Recursive\n",lb->Proc, yo);
      return(ierr);
    }
  }

  return(LB_OK);
}

int LB_Reftree_Build_Recursive(LB *lb,LB_REFTREE *subroot)

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
int bad_case;              /* flag for failing to identify order */
int has_in[2];             /* flag for an element having the parent in vert */
int has_out[2];            /* flag for an element having the parent out vert */
int has_third[2];          /* flag for a triangle having the parent non in/out*/
int not_parent[2];         /* index of a vertex the parent does not have */
int parent_in;             /* index of the parent in vertex */
int parent_out;            /* index of the parent out vertex */
int parent_third;          /* index of the parent triangle non in/out vert */
int parents_vert[6];       /* cross index between children and parent */
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
  case LB_OTHER_REF:
    if (TEMP_first_warning) {
      fprintf(stderr, "[%d] Warning in %s:  Currently not supporting "
                      "automatic ordering of elements for refinement type "
                      "LB_OTHER_REF.  Using LB_IN_ORDER.\n",lb->Proc,yo);
      TEMP_first_warning = 0;
      final_ierr = LB_WARN;
    }
    ref_type = LB_IN_ORDER;
    break;
  case LB_QUAD_QUAD:
    if (TEMP_first_warning) {
      fprintf(stderr, "[%d] Warning in %s:  Currently not supporting "
                      "automatic ordering of elements for refinement type "
                      "LB_QUAD_QUAD.  Using LB_IN_ORDER.\n",lb->Proc,yo);
      TEMP_first_warning = 0;
      final_ierr = LB_WARN;
    }
    ref_type = LB_IN_ORDER;
    break;
  case LB_HEX3D_OCT:
    if (TEMP_first_warning) {
      fprintf(stderr, "[%d] Warning in %s:  Currently not supporting "
                      "automatic ordering of elements for refinement type "
                      "LB_HEX3D_OCT.  Using LB_IN_ORDER.\n",lb->Proc,yo);
      TEMP_first_warning = 0;
      final_ierr = LB_WARN;
    }
    ref_type = LB_IN_ORDER;
    break;
  }

  /* end TEMP */

  /*
   *  Refinement type dependent determination of the order of the children
   *  and the in/out vertices
   */

  switch (ref_type) {

  /*
   * If provided in order, order is the identity permutation and
   * in_vertex and out_vertex were provided by the user or not needed.
   * TEMP not quite -- I said I would pick one if -1 is returned
   */
  case LB_IN_ORDER:
    for (i=0; i<num_obj; i++) order[i] = i;
    break;

  /*
   * For newest node bisection of triangles, determine which of the first
   * and second child has the in_vertex and out_vertex, and find the
   * common vertex to go between them.
   */
  case LB_TRI_BISECT:
    /* verify that 3 vertices were given for each triangle; if not, punt */
    if (vert1[1] != 3 || vert1[2] != 6) {
      fprintf(stderr, "[%d] Warning in %s:  Incorrect number of vertices "
                      "given for bisected triangles.\n",lb->Proc,yo);
      final_ierr = LB_WARN;
      order[0] = 0;
      order[1] = 1;
      in_vertex[0] = vertices[vert1[0]];    /* TEMP I can at least try to */
      out_vertex[0] = vertices[vert1[0]+1]; /* find a common vertex for   */
      in_vertex[1] = vertices[vert1[1]];    /* out[0] and in[1]           */
      out_vertex[1] = vertices[vert1[1]+1];
      break;
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
      final_ierr = LB_WARN;
      order[0] = 0;
      order[1] = 1;
      in_vertex[0] = vertices[vert1[0]];    /* TEMP I can at least try to */
      out_vertex[0] = vertices[vert1[0]+1]; /* find a common vertex for   */
      in_vertex[1] = vertices[vert1[1]];    /* out[0] and in[1]           */
      out_vertex[1] = vertices[vert1[1]+1];
      break;
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
      final_ierr = LB_WARN;
    }
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
    /* TEMP */
    printf("Oops, still got into case for OTHER_REF\n");
    for (i=0; i<num_obj; i++) order[i] = i;
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
    sum_vert += num_vert[i];

  /*
   * Copy from temporary arrays and set empty defaults
   */

    subroot->children[order[i]].global_id      = local_gids[i];
    subroot->children[order[i]].local_id       = local_lids[i];
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
