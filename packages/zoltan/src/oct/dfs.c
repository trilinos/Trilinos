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
#include "octree_const.h"
#include "dfs.h"
#include "costs_const.h"
#include "msg_const.h"
#include "migoct_const.h"
#include "migtags_const.h"
#include "oct_util_const.h"

static int CLOSE = 0;           /* determines criterion for visiting octants */
static int DFS_Part_Count;      /* count of number of times Zoltan_Oct_dfs_partition
                                   was called */
static int partition;              /* Partition number we are working on */
static float total;                /* Cost of all complete partitions so far */
static float pcost;                /* Current partition cost */
static float optcost;              /* Optimal partition cost */
static float pmass;                /* octant volume for partition */
static double pcoord[3];           /* Sum of octant position-volume products */

static void Zoltan_Oct_visit(ZZ *zz,pOctant octant);
static void Zoltan_Oct_visit_all_subtrees(ZZ *zz);
static void Zoltan_Oct_tag_subtree(OCT_Global_Info *OCT_info,pOctant octant, int part);
static void Zoltan_Oct_visit_by_dist(ZZ *zz,pOctant octant, pOctant children[8]);
static int Zoltan_Oct_dfs_SetIds(OCT_Global_Info *OCT_info, pOctant oct, int nprevoct);



/*****************************************************************************/
/*
 * void Zoltan_Oct_dfs_partition()
 * 
 * This function calls the different subfunctions to partition the octree 
 */
void Zoltan_Oct_dfs_partition(ZZ *zz, int *counter, float *c1) {
  float mycost;                     /* cost of the octant */
  float globalcost;                 /* costs of all the octants */
  float prefcost;                   /* sum of costs from previous processors */
/* #ifdef LGG_MIGOCT */
  int nprevoct;                     /* the number of previous octants */
  pRList RootList;                  /* list of the local roots */
  pOctant RootOct;
/* #endif */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);

  DFS_Part_Count = 0;
  *c1 = mycost = Zoltan_Oct_costs_global_compute(OCT_info);
 
/* #ifdef LGG_MIGOCT */
  /* gets the number of octants from the previous processors */
  nprevoct=Zoltan_Oct_msg_int_scan(zz->Communicator, zz->Proc, Zoltan_Oct_nOctants());

  /* iterate through, and advance each id by nprevocts */
  /* this is trying to make all the octant id's to be unquie globally */
  RootList = Zoltan_Oct_POct_localroots(OCT_info);   
  while((RootOct = RL_nextRootOctant(&RootList))) 
    Zoltan_Oct_dfs_SetIds(OCT_info, RootOct, nprevoct);
  
/* #endif */ /* LGG_MIGOCT */

  /* Sum a value from each processor, and return sum to all processors */
  MPI_Allreduce(&mycost,&globalcost,1,MPI_FLOAT,MPI_SUM,zz->Communicator);
  prefcost=Zoltan_Oct_msg_float_scan(zz->Communicator, zz->Proc, mycost);
  
  /* Initialize static vars */
  optcost=globalcost/zz->Num_Proc;               /* Optimal partition size */
  if(optcost > 0)
    partition=(int)(prefcost/optcost);        /* Start work on this partition */
  else
    partition=0;
  if (partition==zz->Num_Proc)
    partition=zz->Num_Proc-1;

  total=partition*optcost;               /* Total cost of all previous parts */
  pcost=prefcost-partition*optcost;                /* Current partition cost */

  pmass=0.0;                            /* initialize octant volume variable */
  vector_set_comp(pcoord,0,0,0);

  Zoltan_Oct_visit_all_subtrees(zz);

  (*counter) = DFS_Part_Count;
}


/* #ifdef LGG_MIGOCT */
/*
 * int Zoltan_Oct_dfs_SetIds(pOctant octant, int number_of_previous_octants)
 *
 * sets the ids of all the octants so that there is a global numbering
 */
static int Zoltan_Oct_dfs_SetIds(OCT_Global_Info *OCT_info, pOctant oct, int nprevoct) {
  int id,                                      /* the id number of an octant */
      i;                                       /* index counter */
  int pid;
  pOctant child;                               /* ith child of an octant */

  if(Zoltan_Oct_isTerminal(oct)) {
    id = Zoltan_Oct_id(oct);
    Zoltan_Oct_setID(oct, id+nprevoct);                     /* now have global id's */
  }
  else {
    for(i=0; i<8; i++) {
      child = Zoltan_Oct_child(oct, i);
      pid = Zoltan_Oct_Cpid(oct,i);
      if ((pid == OCT_info->OCT_localpid) && child != NULL)
	Zoltan_Oct_dfs_SetIds(OCT_info,child, nprevoct);
    }
    id = Zoltan_Oct_id(oct);
    Zoltan_Oct_setID(oct, id+nprevoct);                     /* now have global id's */
  }
  return 0;
}
/* #endif */ /* LGG_MIGOCT */

/*****************************************************************************/
/*
 * void Zoltan_Oct_visit_all_subtrees()
 *
 * visits each of the subtrees that are on the local processor
 */
static void Zoltan_Oct_visit_all_subtrees(ZZ *zz) {
  pRList  RootList;                           /* list of all local roots */
  pOctant RootOct;
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);

  /* get the list of all the local roots */
  /* iterate through each root in localroot list */ 
  /* mark each subtree as being Zoltan_Oct_visited */
  /* and free the costs */
  RootList = Zoltan_Oct_POct_localroots(OCT_info);
  while ((RootOct = RL_nextRootOctant(&RootList))) {
    Zoltan_Oct_visit(zz, RootOct);  
    Zoltan_Oct_costs_free(OCT_info, RootOct);
  }
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_visit(pOctant octant)
 * 
 * This routine references the following (static) global variables:
 *
 *   partition - (RW) number of the partition we are currently working on
 *   total     - (RW) total cost of all *previous* partitions
 *   pcost     - (RW) partition cost for current partition
 *   optcost   - (RO) optimal partition cost
 */
static void Zoltan_Oct_visit(ZZ *zz, pOctant octant) {
  float cost;                 /* Cost of this octant */
  float togo;                 /* Remaining room in current partition */
  float behind;               /* How many to make up for from all prev parts */
  pOctant children[8];        /* children of the octant */
  int i;                      /* index counter */
  COORD origin;               /* center of the octant */
  double volume;              /* volume of the octant */
  double prod[3];             /* product of octant origin and its volume */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);

  DFS_Part_Count++;
  cost = Zoltan_Oct_costs_value(octant);               /* get the cost of the octant */
  behind = partition * optcost - total;          /* calcuate how much behind */

  /* If octant does not overflow the current partition, then use it. */
  if( cost==0 || (pcost+cost) <= (optcost+behind)) {
    Zoltan_Oct_tag_subtree(OCT_info,octant,partition);
    pcost+=cost;
    
    Zoltan_Oct_origin_volume(octant, origin, &volume);
    
    vector_cmult(prod,volume,origin);
    pmass+=volume;
    vector_add(pcoord,pcoord,prod);
    
    return;
  }

  /* 
   * Can't use entire octant because it is too big. If it has suboctants, 
   * visit them.
   */
  
  if (!Zoltan_Oct_isTerminal(octant)) {
    Zoltan_Oct_modify_newpid(octant, partition);                         /* Nonterm */
    Zoltan_Oct_children(octant,children);

    /* currently CLOSE is defined to be 0, a functionality not used */
    if (CLOSE) {
      i=0;
      Zoltan_Oct_visit_by_dist(zz, octant, children);
    }
    else
      for (i=0; i<8; i++)                    /* Simple - just visit in order */
	if(children[i] && Zoltan_Oct_POct_local(OCT_info, octant,i))
	  Zoltan_Oct_visit(zz,children[i]);
    return;
  }
  
  /* 
   * No suboctants!
   * We've hit bottom - have to decide whether to add to
   * the current partition or start a new one.
   */
  togo = behind + optcost - pcost;

  if ((cost-togo) >= togo) {
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
  Zoltan_Oct_modify_newpid(octant, partition);
  pcost += cost;

  Zoltan_Oct_origin_volume(octant, origin, &volume);

  vector_cmult(prod,volume,origin);
  pmass += volume;
  vector_add(pcoord,pcoord,prod);
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_tag_subtree(pOctant octant, int partition_number)
 *
 * marks all the octants within the subtree to be in the current partition
 */
static void Zoltan_Oct_tag_subtree(OCT_Global_Info *OCT_info,pOctant octant, int part) {
  pOctant children[8];                                 /* children of octant */
  int i;                                               /* index counter */

  /* modify NPID so octant know where to migrate to */
  Zoltan_Oct_modify_newpid(octant, part);

  if (Zoltan_Oct_isTerminal(octant))
    return;

  /* if octant has children, have to tag them too */
  Zoltan_Oct_children(octant,children);
  
  for (i=0; i<8; i++)                        /* Simple - just visit in order */
    /* if (children[i] && Zoltan_Oct_local(OCT_info, children[i])) */
    if(children[i] && Zoltan_Oct_POct_local(OCT_info, octant,i))
      Zoltan_Oct_tag_subtree(OCT_info,children[i],part);
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_dfs_migrate()
 *
 * sets up information so the migrate octant routines can create the
 * proper export_tags and import_tags arrays
 */
void Zoltan_Oct_dfs_migrate(ZZ *zz, int *nsentags,
		    pRegion *import_regs, int *nrectags, 
		    float *c2, float *c3, int *counter3, int *counter4) 
{
  pRList RootList;                  /* list of the local roots */
  pOctant oct;                                /* octree octant */
  pOctant *docts = NULL;                      /* array of octants being sent */
  int *dpids = NULL;                          /* array of octant pids */
  int dcount;                                 /* count of octants being sent */
  int pid;                                    /* processor id */
  int nrecocts;
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);
  char *yo = "Zoltan_Oct_dfs_migrate";

  if(Zoltan_Oct_nOctants()) {        /* allocate space for octants being migrated */
    docts = (pOctant *) ZOLTAN_MALLOC(Zoltan_Oct_nOctants() * sizeof(pOctant));
    if(!docts) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "cannot allocate arrays.");
      abort();
    }
    dpids = (int *) ZOLTAN_MALLOC(Zoltan_Oct_nOctants() * sizeof(int));
    if(!dpids) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "cannot allocate arrays.");
      ZOLTAN_FREE(&docts);
      abort();
    }
  }

  dcount=0;

  /* go through the local octants and make sure each has a valid npid */
  RootList = Zoltan_Oct_POct_localroots(OCT_info);  

  while((oct = RL_nextRootOctant(&RootList))) 
    while(oct) {
      pid = Zoltan_Oct_data_newpid(oct);
      if (pid<0 || pid>=zz->Num_Proc) {
	fprintf(stderr,"%d Zoltan_Oct_dfs_migrate: bad dest pid %d\n", zz->Proc, pid);
	abort();
      }
      if (dcount<Zoltan_Oct_nOctants()) {   
	docts[dcount]=oct;
	dpids[dcount]=pid;
      }
      oct=Zoltan_Oct_POct_nextDfs(OCT_info, oct);
      dcount++;
    }

  if (dcount!=Zoltan_Oct_nOctants()) {
    fprintf(stderr, "ERROR: in Zoltan_Oct_dfs_migrate, octant count mismatch (I counted %d but there should be %d)\n",dcount, Zoltan_Oct_nOctants());
 /*    dcount=Zoltan_Oct_nOctants();  */
 /*    abort(); */
  }

  /* setup the import_regs */
  Zoltan_Oct_migrate_objects(zz, docts, dpids, dcount, nsentags,
                     import_regs, nrectags, c2, c3, counter3, counter4);
  Zoltan_Oct_migrate_octants(zz, dpids, docts, dcount, &nrecocts);

  ZOLTAN_FREE(&docts);
  ZOLTAN_FREE(&dpids);
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_visit_by_dist()
 *
 * tries to find the closest child to add to the partition 
 */
static void Zoltan_Oct_visit_by_dist(ZZ *zz,pOctant octant, pOctant children[8])
{
  COORD min,                    /* min bounds of the octant */
        max;                    /* max bounds of the octant */
  COORD cmin[8],                /* array of min bounds for octant's children */
        cmax[8];                /* array of max bounds for octant's children */
  COORD origin;                 /* the origin of the octant */
  COORD corigin[8];             /* array of origin pnts of octant's children */
  COORD pcentroid;              /* centroid of the octant */
  int i;                        /* index counter */
  int minchild;                 /* lowest numbered child */
  double dist;                  /* distance */
  double mindist;               /* lowest distance */
  int visited[8];               /* flag showing which child has been visited */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);

  /* initializing data */
  mindist=0;
  pcentroid[0] = pcentroid[1] = pcentroid[2] = 0;

  /* get the bounds of the octant */
  Zoltan_Oct_bounds(octant,min,max);

  /* use bounds to find octant's origin */
  Zoltan_Oct_bounds_to_origin(min,max,origin);

  /* get bounds, and origin for each of the children */
  for (i=0; i<8; i++) {
    visited[i]=0;
    Zoltan_Oct_child_bounds(min,max,origin,i,cmin[i],cmax[i]);
    Zoltan_Oct_bounds_to_origin(cmin[i],cmax[i],corigin[i]);
  }
  
  /* Visit child closest to centroid */
  for(minchild=0; minchild>=0; ) {      
    minchild= -1;
    
    if (pmass>0)
      vector_divc(pcentroid,pcoord,pmass);
    
    /* for each of the child, find the one with the closest distance */
    for (i=0; i<8; i++)
      /* if ((Zoltan_Oct_POct_local(children[i])) && (!visited[i])) { */
      if(Zoltan_Oct_POct_local(OCT_info, octant, i) && !visited[i]) {
	dist=vector_dist(pcentroid,corigin[i]);
	if (pmass==0 || minchild<0 || dist<mindist) {
	  mindist=dist;
	  minchild=i;
	}
      }

    if (minchild>=0) {
      /* visit that child, so that it can be pu into the partition */
      Zoltan_Oct_visit(zz, children[minchild]);
      /* mark the child as having been visited */
      visited[minchild]=1;        
    }
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
