
#include "zz_const.h"
#include "rcb.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RCB_Build_Structure(ZZ *zz, int *num_obj, int *max_obj, int wgtflag,
                               double overalloc, int use_ids, int gen_tree)
{
/*
 *  Function to build the geometry-based data structures for 
 *  Steve Plimpton's RCB implementation.
 */
char *yo = "Zoltan_RCB_Build_Structure";
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */
struct rcb_tree *treeptr;
int i, ierr = 0;

  /*
   * Allocate an RCB data structure for this Zoltan structure.
   * If the previous data structure is still there, free the Dots and IDs first;
   * the other fields can be reused.
   */

  if (zz->LB.Data_Structure == NULL) {
    rcb = (RCB_STRUCT *) ZOLTAN_MALLOC(sizeof(RCB_STRUCT));
    if (rcb == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    zz->LB.Data_Structure = (void *) rcb;
    rcb->Tree_Ptr = NULL;
    rcb->Box = NULL;
    rcb->Global_IDs = NULL;
    rcb->Local_IDs = NULL;
    memset(&(rcb->Dots), 0, sizeof(struct Dot_Struct));

    Zoltan_Initialize_Transformation(&(rcb->Tran));

    rcb->Box = (struct rcb_box *) ZOLTAN_MALLOC(sizeof(struct rcb_box));
    if (rcb->Box == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      Zoltan_RCB_Free_Structure(zz);
      return(ZOLTAN_MEMERR);
    }
    if (gen_tree) {
      rcb->Tree_Ptr = (struct rcb_tree *)
        ZOLTAN_CALLOC(zz->LB.Num_Global_Parts, sizeof(struct rcb_tree));
      if (rcb->Tree_Ptr == NULL) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        Zoltan_RCB_Free_Structure(zz);
        return(ZOLTAN_MEMERR);
      }
      /* initialize Tree_Ptr */
      for (i = 0; i < zz->LB.Num_Global_Parts; i++) {
         treeptr = &(rcb->Tree_Ptr[i]);
         /* initialize dim to -1 to prevent use of cut */
         treeptr->dim = -1;
         treeptr->cut = 0.0;
         treeptr->parent = treeptr->left_leaf = treeptr->right_leaf = 0;
      }
    }
  }
  else {
    rcb = (RCB_STRUCT *) zz->LB.Data_Structure;
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    Zoltan_Free_And_Reset_Dot_Structure(&rcb->Dots);
  }

  ierr = Zoltan_RB_Build_Structure(zz, &(rcb->Global_IDs), &(rcb->Local_IDs),
                               &(rcb->Dots), num_obj, max_obj,
                               &(rcb->Num_Dim),
                               wgtflag, overalloc, use_ids, 1);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Error returned from Zoltan_RB_Build_Structure.");
    Zoltan_RCB_Free_Structure(zz);
    return(ierr);
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_RCB_Free_Structure(ZZ *zz)
{
/*
 * Deallocate the persistent RCB data structures in zz->Structure.
 */
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */

  rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);

  if (rcb != NULL) {
    ZOLTAN_FREE(&(rcb->Tree_Ptr));
    ZOLTAN_FREE(&(rcb->Box));
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    Zoltan_Free_And_Reset_Dot_Structure(&rcb->Dots);
    ZOLTAN_FREE(&(zz->LB.Data_Structure));

  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#define COPY_BUFFER(buf, type, num) \
  if (from->buf) { \
    to->buf = (type *)ZOLTAN_MALLOC((num) * sizeof(type)); \
    if (!to->buf) { \
      Zoltan_RCB_Free_Structure(toZZ); \
      ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory."); \
      return(ZOLTAN_MEMERR); \
    } \
    memcpy(to->buf, from->buf, (num) * sizeof(type)); \
  } \
  else { \
    to->buf = NULL; \
  }

int Zoltan_RCB_Copy_Structure(ZZ *toZZ, ZZ const *fromZZ)
{
  char *yo = "Zoltan_RCB_Copy_Structure";
  RCB_STRUCT *to;
  RCB_STRUCT const *from;

  from = (RCB_STRUCT const *)fromZZ->LB.Data_Structure;
  Zoltan_RCB_Free_Structure(toZZ);

  if (!from){
    return(ZOLTAN_OK);
  }

  to = (RCB_STRUCT *)ZOLTAN_MALLOC(sizeof(RCB_STRUCT));
  if (to == NULL) {
    ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }

  toZZ->LB.Data_Structure = (void *)to;
  *to = *from;

  COPY_BUFFER(Tree_Ptr, struct rcb_tree, fromZZ->LB.Num_Global_Parts);

  COPY_BUFFER(Box, struct rcb_box, 1);

  return ZOLTAN_OK;
}

/*
** For debugging purposes, print out the RCB structure
**   Indicate how many objects you want printed out, -1 for all.
*/
void Zoltan_RCB_Print_Structure(ZZ *zz, int howMany)
{
  int num_obj, i, len;
  RCB_STRUCT *rcb;
  struct rcb_tree r;
  struct rcb_box *b;
  int printed = 0;

  rcb = (RCB_STRUCT *)zz->LB.Data_Structure;
  num_obj = Zoltan_Print_Obj_List(zz, rcb->Global_IDs, rcb->Local_IDs, 
    0, NULL, NULL, howMany);

#if 0               /* TODO64 fix this */
  struct Dot_Struct dot;
  for (i=0; rcb->Dots && (i<num_obj); i++){
    dot = rcb->Dots[i];
    printf("(Dots %d) (%6.4f %6.4f %6.4f) (%6.4f %6.4f %6.4f %6.4f) proc %d, part %d, new part %dn",
     i, dot.X[0], dot.X[1], dot.X[2], 
     dot.Weight[0], dot.Weight[1], dot.Weight[2], dot.Weight[3],
     dot.Proc, dot.Input_Part, dot.Part);
    printed = 1;
  }
  if (!printed){
    printf("Dots: NULL\n");
  }
#endif
  len = zz->LB.Num_Global_Parts;
  printed = 0;

  for (i=0; rcb->Tree_Ptr && (i<len); i++){
    r = rcb->Tree_Ptr[i];
    printf("(Tree %d) cut: %6.4f, dim %d, up %d, left %d, right %d\n",
      i, r.cut, r.dim, r.parent, r.left_leaf, r.right_leaf);
    printed=1;
  }
  if (!printed){
    printf("Tree: NULL\n");
  }

  b = rcb->Box;
  if (b){
     printf("Box: (%6.4f, %6.4f) (%6.4f, %6.4f) (%6.4f, %6.4f)\n",
       b->lo[0], b->hi[0],
       b->lo[1], b->hi[1],
       b->lo[2], b->hi[2]);
  }
  else{
    printf("Box: NULL\n");
  }

  Zoltan_Print_Transformation(&(rcb->Tran));
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
