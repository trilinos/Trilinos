#include "lb_const.h"
#include "all_allo_const.h"
#include "octree_const.h"


static void    LB_POct_DfsTraversal(LB *lb, OCT_Global_Info *OCT_info, 
                                    pOctant oct);
static void    LB_POct_printRegionInfo(LB *lb, OCT_Global_Info * OCT_info, 
                                       pOctant oct);


/*****************************************************************************/
/*
 * void LB_POct_init(int processor_id)
 *
 * sets up global variables for the octree partitioner
 */
OCT_Global_Info *LB_POct_init(struct LB_Struct *lb, int pid, int dim) {
  char *yo = "LB_POct_init";
  OCT_Global_Info *OCT_info;
  if((OCT_info = (OCT_Global_Info *) LB_MALLOC(sizeof(OCT_Global_Info))) == NULL) {
    LB_TRACE_EXIT(lb, yo);
    return NULL;
  }
  lb->Data_Structure = (void *) OCT_info;
  OCT_info->OCT_localpid=pid;
  OCT_info->OCT_rootlist=RL_initRootList();
  LB_Oct_initCounters();
  if((dim > 3) || (dim < 2)) {
    fprintf(stderr,"WARNING: illegal dimension, using default (3D).\n");
    OCT_info->OCT_dimension = 3;
  }
  else
    OCT_info->OCT_dimension = dim;
  return(OCT_info);
}

/*****************************************************************************/
/*
 * LB_POct_new(OCT_Global_Info *)
 *
 * create a new octant on the local processor and return
 * a pointer to it.  It will have no parents or children.
 */
extern pOctant LB_POct_new(OCT_Global_Info *OCT_info) {
  pOctant newoct = LB_Oct_new();
  if(!newoct)
    return NULL;
  newoct->ppid = OCT_info->OCT_localpid;
  LB_Oct_modify_newpid(newoct, OCT_info->OCT_localpid);
  
  return newoct;
}
/*****************************************************************************/
/* 
 * void LB_POct_free(pOctant octant)
 *
 * frees up space in memory
 * this does not delete the attached regions, must specifically call
 * LB_Oct_clearRegions 
 */

/* KDDKDDFREE Changed oct to *oct to allow NULL from LB_FREE to propagate back 
 * KDDKDDFREE to the calling routine. */
void LB_POct_free(OCT_Global_Info *OCT_info, pOctant *oct) {

  pRList RootList = LB_POct_localroots(OCT_info);

  /* traverse through local root list, if octant a local root */

  if(LB_Oct_Ppid(*oct) != OCT_info->OCT_localpid) {
    /* KDDKDDFREE Now passing pointer to OCT_rootlist so that, if
     * KDDKDDFREE head of list is deleted, this pointer can be updated
     * KDDKDDFREE appropriately (i.e., no longer points to deleted entry). */
    RL_delRootOctant(OCT_info, &(OCT_info->OCT_rootlist), *oct);
  }

  /* free up space in memory */
  LB_Oct_free(OCT_info, oct);

}

/*****************************************************************************/
/*
 * LB_POct_setparent(pOctant octant, pOctant parent, int parent_processor_id)
 *
 * sets the parent of the octant. If the parent is offprocessor, then
 * add octant to the local root list
 */
void LB_POct_setparent(OCT_Global_Info *OCT_info, pOctant oct, pOctant parent, int ppid) {

  pRList RootList = LB_POct_localroots(OCT_info);

  if(LB_Oct_Ppid(oct) == OCT_info->OCT_localpid) {
    if(ppid != OCT_info->OCT_localpid) {
      RL_addRootOctant(RootList, oct);     /* was local -- now nonlocal */
    }
  }
  else {
    if(ppid == OCT_info->OCT_localpid) {
      /* KDDKDDFREE Now passing pointer to OCT_rootlist so that, if 
       * KDDKDDFREE head of list is deleted, this pointer can be updated
       * KDDKDDFREE appropriately (i.e., no longer points to deleted entry). */
      RL_delRootOctant(OCT_info, &(OCT_info->OCT_rootlist), oct);     /* was foreign -- now local  */
    }
  }

  oct->ppid=ppid;
  oct->parent=parent;
  if(parent && (ppid == OCT_info->OCT_localpid) && (oct->mapidx < 0))
    oct->mapidx=parent->mapidx;

}

/*****************************************************************************/
/*
 * pRList LB_POct_localroots()
 *
 * return the list of local roots
 */
pRList LB_POct_localroots(OCT_Global_Info *OCT_info)
{
  return(OCT_info->OCT_rootlist);
}

/*****************************************************************************/
/*
 * void LB_POct_printResults()
 *
 * prints out the intermediate results of the octree structure
 */
void LB_POct_printResults(LB *lb, OCT_Global_Info *OCT_info) {
  pRList ptr;                                  /* pointer to local root list */

  ptr = OCT_info->OCT_rootlist;
  /* go through each entry in local root list and travers down subtree */
  while(ptr != NULL) {
    LB_POct_DfsTraversal(lb, OCT_info, ptr->oct);
    ptr = ptr->next;
  }
}

/*****************************************************************************/
/*
 * void LB_POct_DfsTraversal(pOctant octant)
 *
 * traverse through the octree in DFS order to get a printout
 */
static void LB_POct_DfsTraversal(LB *lb, OCT_Global_Info *OCT_info, pOctant oct)
{
  int i;                                                    /* index counter */

  if(oct == NULL)
    return;
  if(LB_Oct_isTerminal(oct))
    LB_POct_printRegionInfo(lb, OCT_info, oct);
  else {
    for(i=0; i<8; i++)
      LB_POct_DfsTraversal(lb, OCT_info, oct->child[i]);
    LB_POct_printRegionInfo(lb, OCT_info, oct);
  }
}

/*****************************************************************************/
/*
 * void LB_POct_printRegionInfo(pOctant octant)
 *
 * prints out region information
 */
static void LB_POct_printRegionInfo(LB *lb, OCT_Global_Info * OCT_info, pOctant oct) {
  pRegion ptr;            /* pointer to iterate through octant's region list */
  pOctant parent;

#if 0
  if(!LB_Oct_isTerminal(oct)) {
    fprintf(stderr, "WARNING: printing regions in a non-terminal octant.\n");
    return;
  }
#endif

  parent = LB_Oct_parent(oct);
  printf("(Proc %d) ocant %d:\n",
	 OCT_info->OCT_localpid, oct->id);
  printf("\tbounds\tmin=%f, %f, %f\n\t\t max %f, %f, %f\n",
	 oct->min[0], oct->min[1], oct->min[2],
	 oct->max[0], oct->max[1], oct->max[2]);
  if(parent != NULL)
    printf("\tparent octant: %d", parent->id);
  else
    printf("\tparent octant: NULL");
  printf(" \tmigrate: from %d to %d\n", OCT_info->OCT_localpid, oct->npid);

  if(!LB_Oct_isTerminal(oct)) {
    return;
  }
  ptr = oct->list;

  if(ptr==NULL)
    printf("\tOctant is EMPTY\n");

  while(ptr != NULL) {
    printf("\tGlobal_ID:");
    LB_PRINT_GID(lb, ptr->Global_ID);
    printf(" Local_ID:");
    LB_PRINT_LID(lb, ptr->Local_ID);
    printf(" Proc:%d coord:(%f, %f, %f)\n", 
	   ptr->Proc, ptr->Coord[0], ptr->Coord[1], ptr->Coord[2]);
    /*
      printf("%lf %lf %lf,  %d -> %d\n", 
      ptr->Coord[0], ptr->Coord[1], ptr->Coord[2],
      OCT_info->OCT_localpid,oct->npid);
    */
    ptr = ptr->next;
  }
}

/*****************************************************************************/
/*
 * pOctant LB_POct_nextDfs(pOctant octant)
 *
 * returns the next octant in a DFS ordering
 */
pOctant LB_POct_nextDfs(OCT_Global_Info *OCT_info, pOctant octant) {
  pOctant parent,                                     /* parent of an octant */
          child;                                      /* child of an octant */
  int pid;
  int i;                                              /* index counter */

  if (!octant)
    return(NULL);

  for (i=0; i<8; i++) {
    child = LB_Oct_child(octant,i);
    pid = LB_Oct_Cpid(octant,i);
    if ((pid == OCT_info->OCT_localpid) && child)
      return(child);          /* Go down */
  }

  parent = LB_Oct_parent(octant);
  pid = LB_Oct_Ppid(octant);
  while ((pid == OCT_info->OCT_localpid) && parent) {
    for (i=octant->which; i<7; i++) {
      child = LB_Oct_child(parent,i+1);
      pid = LB_Oct_Cpid(parent, i+1);
      if ((pid == OCT_info->OCT_localpid) && child)
	return(child);
    }

    octant=parent;                                       /* Go up */
    parent = LB_Oct_parent(octant);
    pid = LB_Oct_Ppid(octant);
  }
  
  return(NULL);         /* no more octants remain in dfs ordering */
}

/*****************************************************************************/
/*
 * int LB_POct_local(pOctant octant, int child_index)
 *
 * returns true if ith child is local
 */
int LB_POct_local(OCT_Global_Info *OCT_info, pOctant octant, int i) {
  if(octant->cpid[i] == OCT_info->OCT_localpid)
    return 1;
  else {
    /*
     *if(octant->cpid[i] == -1)
     * fprintf(stderr, "WARNING: cpid was not set, returning 0 as default.\n");
     */
    return 0;
  }
}

/*****************************************************************************/
/*
 * int LB_POct_delTree(pOctant root)
 *
 * recursivly traverses down root's subtree deleting all the octants
 */
/* KDDKDDFREE Changed root to *root to allow NULL from LB_FREE to propagate 
 * KDDKDDFREE back to calling routine. */
int LB_POct_delTree(OCT_Global_Info *OCT_info, pOctant *root) {
  int i;                                               /* index counter */
  pOctant child;                                       /* child of an octant */
  
  if(*root == NULL)
    return 1;

  if(LB_Oct_isTerminal(*root)) {
    if(LB_Oct_nRegions(*root))
      LB_Oct_clearRegions(*root);
    LB_POct_free(OCT_info, root);
  }
  else {
    for(i=0; i<8; i++) {
      child = LB_Oct_child(*root, i);
      if(child != NULL && LB_POct_local(OCT_info,*root, i)) {
	LB_POct_delTree(OCT_info,&child);
        /* KDDKDDFREE propagate NULL from LB_POct_delTree to root->child */
        (*root)->child[i] = NULL;  
        (*root)->cpid[i]  = -1;
      }
      /* KDDKDDFREE Added this condition so that tests (in other parts of the
       * KDDKDDFREE code) for NULL children work */
      else if (child != NULL) {
        (*root)->child[i] = NULL;
        (*root)->cpid[i]  = -1;
      }
      /* END KDDKDDFREE */
    }
    LB_POct_free(OCT_info, root);
  }
  return 1;
}
