/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "octree_const.h"
#include "costs.h"
#include "octupdate_const.h"

/*
 * costs.c
 *
 * Routines to do a DFS on an octree, calculating the costs for
 * all subtrees.
 */

static void  LB_costs_init(OCT_Global_Info *OCT_info,pOctant octant);
static float LB_costs_subtree_compute(OCT_Global_Info *OCT_info,pOctant octant, int *seq);
static float LB_costs_weight(pOctant octant);

/*
 * void LB_costs_init(pOctant octant)
 *
 * initialize costs for the subtree rooted at octant
 * ATTN: This function may not be necessary anymore
 */
static void LB_costs_init(OCT_Global_Info *OCT_info,pOctant octant) {

  pOctant children[8];                 /* children of the octant */
  int i;                               /* index counter */

  LB_Oct_modify_newpid(octant, OCT_info->OCT_localpid);
  LB_Oct_modify_cost(octant, 0);
  
  if (!LB_Oct_isTerminal(octant)) {
    LB_Oct_children(octant,children);
    for (i=0; i<8; i++) {
      if(LB_POct_local(OCT_info, octant, i))
	LB_costs_init(OCT_info,children[i]);
    }
  }
}

/*
 * void LB_costs_free(pOctant octant)
 *
 * deletes the cost associated with the octant
 */
void LB_costs_free(OCT_Global_Info *OCT_info,pOctant octant) {
  pOctant children[8];                       /* children of the octant */
  int i;                                     /* index counter */

  LB_Oct_modify_cost(octant, 0);

  if (!LB_Oct_isTerminal(octant)) {
    LB_Oct_children(octant,children);
    for (i=0; i<8; i++)
      if(LB_POct_local(OCT_info, octant, i))
	LB_costs_free(OCT_info,children[i]);
  }
}

/*
 * float LB_costs_subtree_compute(pOctant octant, int sequence_number)
 *
 * Do a DFS on the octree, calculating the costs for any given subtree.
 * (Subtree is defined as any octant and all its descendants.)
 *
 * Tag every octant visited with a sequence number, automatically
 * incrementing the sequence number.
 *
 * NOTE: must call LB_costs_init() first
 */
static float LB_costs_subtree_compute(OCT_Global_Info *OCT_info,pOctant octant, int *seq) {
  pOctant children[8];                       /* the children of the octant */
  float c = 0;                               /* cost of each subtree */
  int i = 0;                                 /* index counter */

/* #ifdef LGG_MIGOCT */
  LB_Oct_setID(octant,(*seq)++);               /* set new ID for local ordering */
/* #endif  */ /* LGG_MIGOCT */

  if (!LB_Oct_isTerminal(octant)) {
    /* get the children of each octant */
    LB_Oct_children(octant,children);   
    /* sum the cost for each child to get octant's cost */
    for (i=0; i<8; i++) 
      if(children[i] && LB_POct_local(OCT_info, octant, i))
	c += LB_costs_subtree_compute(OCT_info,children[i], seq); 
  }
  else                                                           /* terminal */
    c=LB_costs_weight(octant);

  /* set the cost data to the octant */
  LB_Oct_modify_cost(octant, c);
  return(c);
}  

/* 
 * float LB_costs_value(pOctant octant)
 *
 * return cost for an octant
 */
float LB_costs_value(pOctant oct)
{ return(oct->cost); }

/*
 * floatg LB_costs_global_compute()                 (was: all_subtree_costs())
 *
 * return costs of all subtrees.
 *
 * Note: MUST set sequence numbers first, otherwise
 * it will not know what order to traverse the different
 * subtrees in. 
 */
float LB_costs_global_compute(OCT_Global_Info *OCT_info) {
  pRList  RootList;                           /* list of all local roots */
  pOctant RootOct;
  int     seq = 0;                            /* sequencing number */
  float   totcost= 0;                         /* total cost local octree */

  /* initialize octants for COST and NPID data tags */
  /* and calculate cost of all the subtree */
/*   LB_POct_printResults(OCT_info); */
  RootList = LB_POct_localroots(OCT_info);
  while((RootOct = RL_nextRootOctant(&RootList))) {
    LB_costs_init(OCT_info, RootOct);
    totcost+=LB_costs_subtree_compute(OCT_info, RootOct, &seq);
  }
  return(totcost);
}

/*
 * float LB_costs_weight(pOctant octant)
 *
 * calculates the cost of an octant. returns the cost..
 */
static float LB_costs_weight(pOctant octant) {
  pRegion region;                                  /* a region from the list */
  float cost;                                      /* cost of the octant */

  region = LB_Oct_regionlist(octant);
  cost=0;

  /* iterate through the region list, summing each of their weights */
  while (region != NULL) {
    cost += region->Weight;
    region = region->next;
  }
  return cost;
}
