/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "all_allo_const.h"
#include "octree_const.h"


/*****************************************************************************/
/*
 * void Zoltan_Oct_POct_init(int processor_id)
 *
 * sets up global variables for the octree partitioner
 */
OCT_Global_Info *Zoltan_Oct_POct_init(ZZ *zz, int pid, int dim) {
  char *yo = "Zoltan_Oct_POct_init";
  OCT_Global_Info *OCT_info;
  if((OCT_info = (OCT_Global_Info *) ZOLTAN_MALLOC(sizeof(OCT_Global_Info))) == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return NULL;
  }
  zz->LB.Data_Structure = (void *) OCT_info;
  OCT_info->OCT_localpid=pid;
  OCT_info->OCT_rootlist=RL_initRootList();
  Zoltan_Oct_initCounters();
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
 * Zoltan_Oct_POct_new(OCT_Global_Info *)
 *
 * create a new octant on the local processor and return
 * a pointer to it.  It will have no parents or children.
 */
extern pOctant Zoltan_Oct_POct_new(OCT_Global_Info *OCT_info) {
  pOctant newoct = Zoltan_Oct_new();
  if(!newoct)
    return NULL;
  newoct->ppid = OCT_info->OCT_localpid;
  Zoltan_Oct_modify_newpid(newoct, OCT_info->OCT_localpid);
  
  return newoct;
}
/*****************************************************************************/
/* 
 * void Zoltan_Oct_POct_free(pOctant octant)
 *
 * frees up space in memory
 * this does not delete the attached regions, must specifically call
 * Zoltan_Oct_clearRegions 
 */

/* KDDKDDFREE Changed oct to *oct to allow NULL from ZOLTAN_FREE to propagate back 
 * KDDKDDFREE to the calling routine. */
void Zoltan_Oct_POct_free(OCT_Global_Info *OCT_info, pOctant *oct) {

  /* KDDKDDFREE  This variable no longer needed.
  pRList RootList = Zoltan_Oct_POct_localroots(OCT_info);
   * KDDKDDFREE */

  /* traverse through local root list, if octant a local root */

  if(Zoltan_Oct_Ppid(*oct) != OCT_info->OCT_localpid) {
    /* KDDKDDFREE Now passing pointer to OCT_rootlist so that, if
     * KDDKDDFREE head of list is deleted, this pointer can be updated
     * KDDKDDFREE appropriately (i.e., no longer points to deleted entry). */
    RL_delRootOctant(OCT_info, &(OCT_info->OCT_rootlist), *oct);
  }

  /* free up space in memory */
  Zoltan_Oct_free(OCT_info, oct);

}

/*****************************************************************************/
/*
 * Zoltan_Oct_POct_setparent(pOctant octant, pOctant parent, int parent_processor_id)
 *
 * sets the parent of the octant. If the parent is offprocessor, then
 * add octant to the local root list
 */
void Zoltan_Oct_POct_setparent(OCT_Global_Info *OCT_info, pOctant oct, pOctant parent, int ppid) {

  pRList RootList = Zoltan_Oct_POct_localroots(OCT_info);

  if(Zoltan_Oct_Ppid(oct) == OCT_info->OCT_localpid) {
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
 * pRList Zoltan_Oct_POct_localroots()
 *
 * return the list of local roots
 */
pRList Zoltan_Oct_POct_localroots(OCT_Global_Info *OCT_info)
{
  return(OCT_info->OCT_rootlist);
}


/*****************************************************************************/
/*
 * pOctant Zoltan_Oct_POct_nextDfs(pOctant octant)
 *
 * returns the next octant in a DFS ordering
 */
pOctant Zoltan_Oct_POct_nextDfs(OCT_Global_Info *OCT_info, pOctant octant) {
  pOctant parent,                                     /* parent of an octant */
          child;                                      /* child of an octant */
  int pid;
  int i;                                              /* index counter */

  if (!octant)
    return(NULL);

  for (i=0; i<8; i++) {
    child = Zoltan_Oct_child(octant,i);
    pid = Zoltan_Oct_Cpid(octant,i);
    if ((pid == OCT_info->OCT_localpid) && child)
      return(child);          /* Go down */
  }

  parent = Zoltan_Oct_parent(octant);
  pid = Zoltan_Oct_Ppid(octant);
  while ((pid == OCT_info->OCT_localpid) && parent) {
    for (i=octant->which; i<7; i++) {
      child = Zoltan_Oct_child(parent,i+1);
      pid = Zoltan_Oct_Cpid(parent, i+1);
      if ((pid == OCT_info->OCT_localpid) && child)
	return(child);
    }

    octant=parent;                                       /* Go up */
    parent = Zoltan_Oct_parent(octant);
    pid = Zoltan_Oct_Ppid(octant);
  }
  
  return(NULL);         /* no more octants remain in dfs ordering */
}

/*****************************************************************************/
/*
 * int Zoltan_Oct_POct_local(pOctant octant, int child_index)
 *
 * returns true if ith child is local
 */
int Zoltan_Oct_POct_local(OCT_Global_Info *OCT_info, pOctant octant, int i) {
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
 * int Zoltan_Oct_POct_delTree(pOctant root)
 *
 * recursivly traverses down root's subtree deleting all the octants
 */
/* KDDKDDFREE Changed root to *root to allow NULL from ZOLTAN_FREE to propagate 
 * KDDKDDFREE back to calling routine. */
int Zoltan_Oct_POct_delTree(OCT_Global_Info *OCT_info, pOctant *root) {
  int i;                                               /* index counter */
  pOctant child;                                       /* child of an octant */
  
  if(*root == NULL)
    return 1;

  if(Zoltan_Oct_isTerminal(*root)) {
    if(Zoltan_Oct_nRegions(*root))
      Zoltan_Oct_clearRegions(*root);
    Zoltan_Oct_POct_free(OCT_info, root);
  }
  else {
    for(i=0; i<8; i++) {
      child = Zoltan_Oct_child(*root, i);
      if(child != NULL && Zoltan_Oct_POct_local(OCT_info,*root, i)) {
	Zoltan_Oct_POct_delTree(OCT_info,&child);
        /* KDDKDDFREE propagate NULL from Zoltan_Oct_POct_delTree to root->child */
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
    Zoltan_Oct_POct_free(OCT_info, root);
  }
  return 1;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
