#include "lb_const.h"
#include "all_allo_const.h"
#include "rootlist_const.h"

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
  pRList sentinal1 = (pRList) LB_MALLOC(sizeof(RList));
  pRList sentinal2 = (pRList) LB_MALLOC(sizeof(RList));
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

  node = (pRList) LB_MALLOC(sizeof(RList));

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

int RL_delRootOctant(pRList rlist, pOctant oct) {
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
      LB_FREE(&temp);
      result = 0;
    }
    if(result == 0)
      break;
    prev = rlist;
  }

  return result;
}

int RL_clearRootOctants(pRList rlist) {
  pRList  head;
  if((rlist == NULL) || (rlist->oct != NULL))
    return -1;
  head = rlist;
  rlist = rlist->next;
  while(rlist->next != head) {  
    head->next = rlist->next;
    LB_FREE(&rlist);
    rlist = head->next;
  }
  return LB_OK;
}

int RL_freeList(pRList rlist) {
  pRList  head;
  if((rlist == NULL) || (rlist->oct != NULL))
    return -1;
  head = rlist;
  rlist = rlist->next;
  while(rlist != head) {  
    head->next = rlist->next;
    LB_FREE(&rlist);
    rlist = head->next;
  }
  LB_FREE(&rlist);
  return LB_OK;
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
    LB_Oct_bounds(rootoct,rmin,rmax);
    fprintf(stderr, "RL_printRootOctants ppid %d id %d area %f  npid %d ", rootoct->ppid, rootoct->id, rootoct->area, rootoct->npid);   
    fprintf(stderr,"min box %f %f %f ",rmin[0], rmin[1], rmin[2]);
    fprintf(stderr,"max box %f %f %f\n",rmax[0], rmax[1], rmax[2]);
  }
  return 0;
}
