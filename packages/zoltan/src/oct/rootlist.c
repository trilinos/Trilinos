/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "zz_const.h"
#include "all_allo_const.h"
#include "rootlist_const.h"
#include "octree_const.h"

static int RL_OctGreaterThan(pOctant oct1, pOctant oct2) {
  if(oct2->type == LOCALOCT)
    return (oct1->id > oct2->id);
  return (oct1->area > oct2->area);
}

static int RL_OctEqual(pOctant oct1, pOctant oct2) {
  if(oct2->type == LOCALOCT)
    return (oct1->id == oct2->id);
  return (oct1->area == oct2->area);
}

pRList  RL_initRootList() {
  pRList sentinal1 = (pRList) ZOLTAN_MALLOC(sizeof(RList));
  pRList sentinal2 = (pRList) ZOLTAN_MALLOC(sizeof(RList));
  sentinal1->oct = NULL;
  sentinal1->next = sentinal2;
  sentinal2->oct = NULL;
  sentinal2->next = sentinal1;
  return sentinal1;
}

pOctant RL_nextRootOctant(pRList *rlist) {
  if(*rlist == NULL)
    return NULL;
  *rlist = (*rlist)->next;
  return (*rlist)->oct;
}

int RL_addRootOctant(pRList rlist, pOctant oct) {
  pRList  prev;
  pRList  node;
  pOctant rootoct;

  if((oct == NULL) || (rlist == NULL) || (rlist->oct != NULL))
    return -1;

  node = (pRList) ZOLTAN_MALLOC(sizeof(RList));

  prev = rlist;
  while((rootoct = RL_nextRootOctant(&rlist))) {
    if(!rootoct || RL_OctGreaterThan(rootoct,oct)) 
      break;
    prev = rlist;
  }
  prev->next = node;
  node->next = rlist;
  node->oct = oct;
  return 0;
}

/* KDDKDDFREE changed rlist to *rootlist so that if head of list is deleted,
 * KDDKDDFREE a new head pointer can be propagated back to the calling routine.
 */
int RL_delRootOctant(OCT_Global_Info *OCT_info, pRList *rootlist, pOctant oct) {
  pRList  rlist = *rootlist;
  pRList  temp;
  pRList  prev;
  pOctant rootoct;
  int     result = -1;

  if((oct == NULL) || (rlist == NULL) || (rlist->oct != NULL))
    return -1;
  
  prev = rlist;
  while((rootoct = RL_nextRootOctant(&rlist))) {
    while(rootoct && RL_OctEqual(rootoct,oct)) {
      temp = rlist;
      rootoct = RL_nextRootOctant(&rlist);
      prev->next = rlist;
      /* KDDKDDFREE  Update *rootlist (head of list) if the lead entry is
       * KDDKDDFREE  being deleted.  */
      if (temp == *rootlist) 
        *rootlist = rlist;
      /* KDDKDDFREE  If there is only one item in list and it will be deleted,
       * KDDKDDFREE  set head of list to NULL. */
      if (temp == prev && temp == rlist)
        *rootlist = NULL;
      /* END KDDKDDFREE */
      ZOLTAN_FREE(&temp);
      result = 0;
    }
    if(result == 0)
      break;
    prev = rlist;
  }

  return result;
}

/* KDDKDDFREE changed rlist to *rlist to allow NULL from ZOLTAN_FREE to propagate 
 * KDDKDDFREE back to the calling routine.  */
int RL_clearRootOctants(pRList *rlist) {
  pRList  head;
  if((*rlist == NULL) || ((*rlist)->oct != NULL))
    return -1;
  head = *rlist;
  *rlist = (*rlist)->next;
  while((*rlist)->next != head) {  
    head->next = (*rlist)->next;
    ZOLTAN_FREE(rlist);
    *rlist = head->next;
  }
  return ZOLTAN_OK;
}

/* KDDKDDFREE changed rlist to *rlist to allow NULL from ZOLTAN_FREE to propagate 
 * KDDKDDFREE back to the calling routine.  */
int RL_freeList(pRList *rlist) {
  pRList  head;
  if((*rlist == NULL) || ((*rlist)->oct != NULL))
    return -1;
  head = *rlist;
  *rlist = (*rlist)->next;
  while(*rlist != head) {  
    head->next = (*rlist)->next;
    ZOLTAN_FREE(rlist);
    *rlist = head->next;
  }
  ZOLTAN_FREE(rlist);
  return ZOLTAN_OK;
}

int RL_numRootOctants(pRList rlist) {
  int nroots = 0;

  if((rlist == NULL) || (rlist->oct != NULL))
    return -1;

  while(RL_nextRootOctant(&rlist)) 
    nroots++;

  return nroots;
}

int RL_printRootOctants(pRList rlist) {
  pOctant rootoct;
  COORD rmin, rmax;
  if((rlist == NULL) || (rlist->oct != NULL))
    return -1;

  while((rootoct = RL_nextRootOctant(&rlist))) {
    Zoltan_Oct_bounds(rootoct,rmin,rmax);
    fprintf(stderr, "RL_printRootOctants ppid %d id %d area %f  npid %d ", rootoct->ppid, rootoct->id, rootoct->area, rootoct->npid);   
    fprintf(stderr,"min box %f %f %f ",rmin[0], rmin[1], rmin[2]);
    fprintf(stderr,"max box %f %f %f\n",rmax[0], rmax[1], rmax[2]);
  }
  return 0;
}
