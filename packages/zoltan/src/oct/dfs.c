#include <stdio.h>
#include <stdlib.h>
#include "octant_const.h"
#include "dfs.h"
#include "costs_const.h"
#include "msg_const.h"
#include "migoct_const.h"
#include "util_const.h"

/* type definition for coordinants */
typedef double Coord[3];

/* determines criterion for visiting octants */
static int CLOSE = 0;
int count;
/* static int dimensions = 8; */

void dfs_set_visit_criterion(int visit) {
  if((visit < 0) || (visit > 1)) {
    fprintf(stderr, "%s\n%s\n", 
	 "Warning: invalid visit criterion passed to dfs_set_visit_criterion.",

	 "Using default method of visiting octants. (CLOSE = 0)");
  }
  else {
    CLOSE = visit;
  }
}

/*
 * void dfs_partition()
 * 
 * This function calls the different subfunctions to partition the octree 
 *
 */
void dfs_partition(int *counter) {
  float mycost,                     /* cost of the octant */
        globalcost,                 /* costs of all the octants */
        prefcost;                   /* sum of costs from previous processors */
  int nprevoct;                     /* the number of previous octants */
  int id;                           /* octant id */
  void *temp;                       /* temp variable for iterations */
  pRList localroots;                /* list of the local roots */
  pOctant lr,                       /* pointer to local root octant */
          oct;                      /* pointer to an octant */

  count = 0;
  /* Note: set root oct id's first so costs_global_compute() can know order */
  mycost=costs_global_compute();                 /* sets ids in order from 0 */
  /* gets the number of octants from the previous processors */
  nprevoct=msg_int_scan(POC_nOctants());

  /* add nprevoct to the id counter */
  localroots = POC_localroots();                  /* get the local root list */
  /* iterate through, and advance each id by nprevocts */
  /* this is trying to make all the octant id's to be unquie globally */
  while(localroots != NULL) {
    dfs_SetIds(localroots->oct, nprevoct);
    fprintf(stderr, "(%d) set local root id %d\n", msg_mypid, 
	    POC_id(localroots->oct));
    localroots = localroots->next;
  }

  /* Sum a value from each processor, and return sum to all processors */
  globalcost=msg_double_sum(mycost);
  prefcost=msg_float_scan(mycost);
  
  /* Initialize static vars */
  optcost=globalcost/msg_nprocs;                   /* Optimal partition size */
  partition=prefcost/optcost;                /* Start work on this partition */
  if (partition==msg_nprocs)
    partition=msg_nprocs-1;

  total=partition*optcost;               /* Total cost of all previous parts */
  pcost=prefcost-partition*optcost;                /* Current partition cost */

  pmass=0.0;                            /* initialize octant volume variable */
  vector_set_comp(pcoord,0,0,0);

#ifdef DEBUG 
  PRINT_IN_ORDER()
    printf("Global %6.1f  Pref %6.1f  Opt %6.1f  Pnum %2d  Pcost %8.2f\n",
	   globalcost,prefcost,optcost,partition,pcost);
#endif

  visit_all_subtrees();

  /* PRINT_IN_ORDER()
    printf("Partitions = %d\n",partition); */
  (*counter) = count;
}

int dfs_SetIds(pOctant oct, int nprevoct) {
  int id,                                      /* the id number of an octant */
      i;                                       /* index counter */
  pOctant child;                               /* ith child of an octant */

  if(POC_isTerminal(oct)) {
    id = POC_id(oct);
    POC_setID(oct, id+nprevoct);                     /* now have global id's */
  }
  else {
    for(i=0; i<8; i++) {
      child = POC_child(oct, i);
      dfs_SetIds(child, nprevoct);
    }
    id = POC_id(oct);
    POC_setID(oct, id+nprevoct);                     /* now have global id's */
  }
  return 0;
}

void visit_all_subtrees() {
  pRList localroots;                    /* list of all local roots */

  /* get the list of all the local roots */
  localroots = POC_localroots();

  /* iterate through each root in localroot list */
  while (localroots != NULL) {
    fprintf(stderr, "(%d) local root %d\n", msg_mypid, 
	    POC_id(localroots->oct));
    visit(localroots->oct);   /* mark each subtree as being visited */
#if 0
    costs_free(localroots->oct); /* clear attached cost data (is it needed?) */
#endif
    localroots = localroots->next;
  }
}

/*
 * visit(oct)
 * 
 * This routine references the following (static) global variables:
 *
 *   partition - (RW) number of the partition we are currently working on
 *   total     - (RW) total cost of all *previous* partitions
 *   pcost     - (RW) partition cost for current partition
 *   optcost   - (RO) optimal partition cost
 */
void visit(pOctant octant) {
  float cost;                 /* Cost of this octant */
  float togo;                 /* Remaining room in current partition */
  float behind;               /* How many to make up for from all prev parts */
  pOctant children[8];        /* children of the octant */
  int i;                      /* index counter */
  double origin[3];           /* center of the octant */
  double size[3];             /* size of the octant */
  double volume;              /* volume of the octant */
  double prod[3];             /* product of octant origin and its volume */

  count++;
  cost = costs_value(octant);                  /* get the cost of the octant */
  behind = partition * optcost - total;          /* calcuate how much behind */

  /* If octant does not overflow the current partition, then use it. */
  if( cost==0 || (pcost+cost) <= (optcost+behind)) {
    tag_subtree(octant,partition);
    pcost+=cost;
    
    POC_origin_volume(octant, origin, &volume);
    
    vector_cmult(prod,volume,origin);
    pmass+=volume;
    vector_add(pcoord,pcoord,prod);
    
    return;
  }

  /* 
   * Can't use entire octant because it is too big. If it has suboctants, 
   * visit them.
   */
  
  if (!POC_isTerminal(octant)) {
    POC_modify_newpid(octant, partition);                         /* Nonterm */
    POC_children(octant,children);

    /* currently CLOSE is defined to be 0, a functionality not used */
    if (CLOSE) {
      i=0;
      visit_by_dist(octant,children);
    }
    else
      for (i=0; i<8; i++)                    /* Simple - just visit in order */
	/* if (children[i] && POC_local(children[i])) */
	if(children[i] && POC_local(octant,i))
	  visit(children[i]);
    return;
  }
  
  /* 
   * No suboctants!
   * We've hit bottom - have to decide whether to add to
   * the current partition or start a new one.
   */
  togo = behind + optcost - pcost;

  if ((cost-togo) >= togo ) {
    /*
     * End current partition and start new one. We are more "over" than "under"
     */
    partition++;                                /* Move on to next partition */
    total += pcost;
    pcost = 0;
    pmass = 0;
    vector_set_comp(pcoord,0,0,0);
  }

  /*** Add terminal octant to current partition */
  POC_modify_newpid(octant, partition);
  pcost += cost;

  /* POC_size_and_origin(octant, size, origin); */
  POC_origin_volume(octant, origin, &volume);

  vector_cmult(prod,volume,origin);
  pmass += volume;
  vector_add(pcoord,pcoord,prod);
}

void tag_subtree(pOctant octant, int partition) {
  pOctant children[8];                                 /* children of octant */
  int i;                                               /* index counter */

  /* modify NPID so octant know where to migrate to */
  POC_modify_newpid(octant, partition);

  if (POC_isTerminal(octant))
    return;

  /* if octant has children, have to tag them too */
  POC_children(octant,children);
  
  for (i=0; i<8; i++)                        /* Simple - just visit in order */
    /* if (children[i] && POC_local(children[i])) */
    if(children[i] && POC_local(octant,i))
      tag_subtree(children[i],partition);
}

void dfs_migrate(LB_TAG **export_tags, int *nsentags,
		 LB_TAG **import_tags, int *nrectags) 
{
  pRList lroots;                              /* list of all local roots */
  pOctant root;                               /* root of a subtree */
  pOctant oct;                                /* octree octant */
  pOctant *docts = NULL;                      /* array of octants being sent */
  int *dpids = NULL;                          /* array of octant pids */
  int dcount;                                 /* count of octants being sent */
  int pid;                                    /* processor id */
  int i;                                      /* index counter */

  i = POC_nOctants();
  if(POC_nOctants()) {
    docts = (pOctant *)malloc(sizeof(pOctant)*POC_nOctants());
    dpids = (int *)malloc(sizeof(int)*POC_nOctants());
    if((docts == NULL) || (dpids == NULL)) {
      fprintf(stderr, "Dfs_Migrate: cannot allocate arrays.\n");
      abort();
    }
  }

  dcount=0;

  lroots = POC_localroots();
  while(lroots != NULL) {
    for (oct=lroots->oct; oct; oct=POC_nextDfs(oct)) {
      pid = POC_data_newpid(oct);
      if (pid<0 || pid>=msg_nprocs) {
	fprintf(stderr,"%d dfs_migrate: bad dest pid %d\n",
		msg_mypid,pid);
	abort();
      }
      docts[dcount]=oct;
      dpids[dcount++]=pid;
    }
    lroots = lroots->next;
  }

  if (dcount!=POC_nOctants()) {
    fprintf(stderr, "ERROR: in dfs_migrate, octant count mismatch\n");
    abort();
  }

  Migrate_Objects(docts, dpids, dcount, export_tags, nsentags, 
		  import_tags, nrectags);

  free(docts);
  free(dpids);
  
#if 0
  PRINT_IN_ORDER()
    printf("dfs_migrate: sending %d regions\n",nregions);
  {
    pMesh mesh;
    int local_meshregions, global_meshregions;
    int global_migrate;

    mesh=pmdb_get_pmesh(pmeshpb);
    local_meshregions=M_nRegion(mesh);
    global_meshregions=msg_int_sum(local_meshregions);

    global_migrate=msg_int_sum(nregions);
    
    if(msg_mypid == 0) { 
      printf("OCTPART volume: %8d of %8d = %.3f\n",
	     global_migrate,global_meshregions,
	     (double)global_migrate/global_meshregions);
    }
  }
  
  migreg_migrate_regions(pmeshpb,migregions,nregions,TRUE);
  free(migregions);
#endif
}

/*
 * void visit_by_dist()
 *
 * tries to find the closest child to add to the partition 
 */
void visit_by_dist(pOctant octant, pOctant children[8])
{
  double min[3],                /* min bounds of the octant */
         max[3];                /* max bounds of the octant */
  Coord cmin[8],                /* array of min bounds for octant's children */
        cmax[8];                /* array of max bounds for octant's children */
  double origin[3];             /* the origin of the octant */
  Coord corigin[8];             /* array of origin pnts of octant's children */
  double pcentroid[3];          /* centroid of the octant */
  int i;                        /* index counter */
  int minchild;                 /* lowest numbered child */
  double dist,                  /* distance */
         mindist;               /* lowest distance */
  int visited[8];               /* flag showing which child has been visited */

  /* initializing data */
  mindist=0;
  pcentroid[0] = pcentroid[1] = pcentroid[2] = 0;

  /* get the bounds of the octant */
  POC_bounds(octant,min,max);

  /* use bounds to find octant's origin */
  bounds_to_origin(min,max,origin);

  /* get bounds, and origin for each of the children */
  for (i=0; i<8; i++) {
    visited[i]=0;
    child_bounds(min,max,origin,i,cmin[i],cmax[i]);
    bounds_to_origin(cmin[i],cmax[i],corigin[i]);
  }
  
 
  for(minchild=0; minchild>=0; ) {        /* Visit child closest to centroid */
    minchild= -1;
    
    if (pmass>0)
      vector_divc(pcentroid,pcoord,pmass);
    
    /* for each of the child, find the one with the closest distance */
    for (i=0; i<8; i++)
      /* if ((POC_local(children[i])) && (!visited[i])) { */
      if(POC_local(octant, i) && !visited[i]) {
	dist=vector_dist(pcentroid,corigin[i]);
	if (pmass==0 || minchild<0 || dist<mindist) {
	  mindist=dist;
	  minchild=i;
	}
      }

    if (minchild>=0) {
      /* visit that child, so that it can be pu into the partition */
      visit(children[minchild]);
      visited[minchild]=1;          /* mark the child as having been visited */
    }
  }
}

