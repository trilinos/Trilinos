#include "octant_const.h"
#include "costs.h"
#include "octupdate_const.h"

/*
 * costs.c
 *
 * Routines to do a DFS on an octree, calculating the costs for
 * all subtrees.
 */


/*
 * costs_init(octant)
 *
 * initialize costs for the subtree rooted at octant
 * ATTN: This function may not be necessary anymore
 */
void costs_init(pOctant octant) {
  pOctant children[8];
  float *data;
  int i;

  /*  OEN_attachDataI(octant, "NPID", -1);
      SAFE_MALLOC(data,float*,sizeof(float));
      OEN_attachDataP(octant, "COST", data);
   */
  POC_modify_newpid(octant, OCT_localpid);
  POC_modify_cost(octant, 0);
  
  if (!POC_isTerminal(octant)) {
    POC_children(octant,children);
    for (i=0; i<8; i++) {
      /* if (POC_local(children[i])) */
      if(POC_local(octant, i))
	costs_init(children[i]);
    }
  }
}

/*
 * costs_free(octant)
 *
 */
void costs_free(pOctant octant) {
  float *data;

  pOctant children[8];
  int i;

  /* should free seq number? */

  /* POC_modify_newpid(octant, OCT_localpid); */
  POC_modify_cost(octant, 0);

  if (!POC_isTerminal(octant)) {
    POC_children(octant,children);
    for (i=0; i<8; i++)
      /* if (POC_local(children[i])) */
      if(POC_local(octant, i))
	costs_free(children[i]);
  }
}

/*
 * costs_subtree_compute(pmeshpb,octant,seq)
 *
 * Do a DFS on the octree, calculating the costs for any given subtree.
 * (Subtree is defined as any octant and all its descendants.)
 *
 * Tag every octant visited with a sequence number, automatically
 * incrementing the sequence number.
 *
 * NOTE: must call costs_init() first
 */
/* REMOVED: parameter "type" */
float costs_subtree_compute(pOctant octant, int *seq) {
  pOctant children[8];                       /* the children of the octant */
  float c;                                   /* cost of each subtree */
  float *data;                               /* COST data attached to octant */
  int i;                                     /* index counter */

  POC_setID(octant,(*seq)++);               /* set new ID for local ordering */
  c=0;                                           /* initialize cost variable */

  if (!POC_isTerminal(octant)) {
    /* get the children of each octant */
    POC_children(octant,children);
    
    /* sum the cost for each child to get octant's cost */
    for (i=0; i<8; i++) {
      /* if (children[i] && (POC_local(children[i]))) */
      if(children[i] && POC_local(octant, i))
	c += costs_subtree_compute(children[i], seq);
    }
  }
  else                                           /* terminal */
    c=costs_weight(octant);

  /* attach the COST data to the octant */
  /* *((float *)OEN_dataP(octant,"COST"))=c; */
  POC_modify_cost(octant, c);
  return(c);
}  

/* 
 * costs_value(octant)
 *
 * return cost for an octant
 */
float costs_value(pOctant oct)
{ return(oct->cost); }

/*
 * costs_global_compute()    (was: all_subtree_costs())
 *
 * return costs of all subtrees.
 *
 * Note: MUST set sequence numbers first, otherwise
 * it will not know what order to traverse the different
 * subtrees in. 
 */
float costs_global_compute() {
  int seq;                                    /* sequencing number */
  float totcost;                              /* total cost local octree */
  int i;                                      /* index counter */
  int nroot;                                  /* number of local roots */
  pOctant *root;                              /* root of a subtree */
  void *temp;                                 /* temp var used for iterating */
  pRList localroots;                          /* list of all local roots */
  pOctant lr,                                 /* a local root */
          oct;                                /* an octant of a subtree */
  /* REMOVED : int type=bc_getparam("COSTTYPE"); */

  /* initialize variables */
  seq=0;
  totcost=0;
  
  /* get the roots in order */
  oct_roots_in_order(&root,&nroot);
  {
    for (i=0; i<nroot; i++) {
      /* initialize octants for COST and NPID data tags */
      costs_init(root[i]);
      /* calculate cost of all the subtree */
      totcost+=costs_subtree_compute(root[i], &seq);
      fprintf(stderr, "Computing costs on local root %d.%d seq=%d tot=%f\n",
	     OCT_localpid, POC_id(root[i]), seq, totcost);
    }
  }

  free(root);

  return(totcost);
}

float costs_weight(pOctant octant) {
  pRegion region;                      /* a region from the list */
  void *temp;                          /* temp var used for iterations */
  float cost;                          /* cost of the octant */

  region = POC_regionlist(octant);
  cost=0;

  temp=NULL;
  while (region != NULL) {
    cost += region->Weight;
    region = region->next;
  }
  
  return cost;
}
