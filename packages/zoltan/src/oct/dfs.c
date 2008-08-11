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

static int DFS_Part_Count;      /* count of number of times Zoltan_Oct_dfs_partition
                                   was called */
static int partition;              /* Partition number we are working on */
static float total;                /* Cost of all complete partitions so far */
static float pcost;                /* Current partition cost */
static float optcost;              /* Optimal partition cost */
static float pmass;                /* octant volume for partition */
static float globalcost;
static double pcoord[3];           /* Sum of octant position-volume products */
static float tmpcost;
static float optsize;

static void Zoltan_Oct_visit(ZZ *zz,pOctant octant, float *part_sizes);
static void Zoltan_Oct_visit_all_subtrees(ZZ *zz, float *part_sizes);
static void Zoltan_Oct_tag_subtree(OCT_Global_Info *OCT_info,pOctant octant, int part);
static int Zoltan_Oct_dfs_SetIds(OCT_Global_Info *OCT_info, pOctant oct, int nprevoct);



/****************************************************************************/
/*
 * void Zoltan_Oct_dfs_partition()
 * 
 * This function calls the different subfunctions to partition the octree 
 */
void Zoltan_Oct_dfs_partition(ZZ *zz, int *counter, float *c1,
			      float *part_sizes) {
  float mycost;                    /* cost of the octant */
  /*float globalcost;*/            /* costs of all the octants */
  float prefcost;                  /* sum of costs from previous processors */
  int nprevoct;                    /* the number of previous octants */
  pRList RootList;                 /* list of the local roots */
  pOctant RootOct;
  float prevwork;
  float pastwork;

  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);

  DFS_Part_Count = 0;
  *c1 = mycost = Zoltan_Oct_costs_global_compute(OCT_info);
 
  /* gets the number of octants from the previous processors */
  nprevoct=Zoltan_Oct_msg_int_scan(zz->Communicator, zz->Proc, 
				   Zoltan_Oct_nOctants());

  /* iterate through, and advance each id by nprevocts */
  /* this is trying to make all the octant id's to be unquie globally */
  RootList = Zoltan_Oct_POct_localroots(OCT_info);   
  while((RootOct = RL_nextRootOctant(&RootList))) 
    Zoltan_Oct_dfs_SetIds(OCT_info, RootOct, nprevoct);
  
  /* Sum a value from each processor, and return sum to all processors */
  MPI_Allreduce(&mycost,&globalcost,1,MPI_FLOAT,MPI_SUM,zz->Communicator);
  prefcost=Zoltan_Oct_msg_float_scan(zz->Communicator, zz->Proc, mycost);
  
  /* Initialize static vars */
  optcost=globalcost/zz->Num_Proc;                /* Optimal partition size */
#if 0
  if(optcost > 0)
    partition=(int)(prefcost/optcost);      /* Start work on this partition */
  else
    partition=0;
  if (partition==zz->Num_Proc)
    partition=zz->Num_Proc-1;
#else
  tmpcost=0;
  prevwork = prefcost;
  if(optcost > 0) {
    partition=(int)(prefcost/optcost);      /* Start work on this partition */
    partition = 0;
    optsize = part_sizes[partition]*globalcost;
    pastwork = part_sizes[partition]*globalcost;
    
    while(pastwork < prevwork) {
      tmpcost += part_sizes[partition];
      partition++;
      optsize = part_sizes[partition]*globalcost;
      pastwork += (part_sizes[partition]*globalcost);
    }
    if((pastwork == prevwork) && (partition < (zz->Num_Proc - 1))){
      tmpcost += part_sizes[partition];
      partition++;
      optsize = part_sizes[partition]*globalcost;
      pastwork += (part_sizes[partition]*globalcost);
    }
  }
  else
    partition=0;
  while((part_sizes[partition] == 0) && (partition < (zz->Num_Proc - 1))) {
    partition++;
    optsize = part_sizes[partition]*globalcost;
  }
#endif
  
  /*optcost = part_sizes[partition]*globalcost;*/
  /*total=partition*optcost;*/          /* Total cost of all previous parts */
  total=tmpcost*globalcost;             /* Total cost of all previous parts */
  /*pcost=prefcost-(partition*optcost);*/         /* Current partition cost */
  pcost=prefcost-(tmpcost*globalcost);            /* Current partition cost */

  pmass=0.0;                           /* initialize octant volume variable */
  vector_set_comp(pcoord,0,0,0);

  Zoltan_Oct_visit_all_subtrees(zz, part_sizes);

  (*counter) = DFS_Part_Count;
}

/*
 * int Zoltan_Oct_dfs_SetIds(pOctant octant, int number_of_previous_octants)
 *
 * sets the ids of all the octants so that there is a global numbering
 */
static int Zoltan_Oct_dfs_SetIds(OCT_Global_Info *OCT_info, pOctant oct,
				 int nprevoct) {
  int id,                                     /* the id number of an octant */
      i;                                      /* index counter */
  int pid;
  pOctant child;                              /* ith child of an octant */

  if(Zoltan_Oct_isTerminal(oct)) {
    id = Zoltan_Oct_id(oct);
    Zoltan_Oct_setID(oct, id+nprevoct);             /* now have global id's */
  }
  else {
    for(i=0; i<8; i++) {
      child = Zoltan_Oct_child(oct, i);
      pid = Zoltan_Oct_Cpid(oct,i);
      if ((pid == OCT_info->OCT_localpid) && child != NULL)
	Zoltan_Oct_dfs_SetIds(OCT_info,child, nprevoct);
    }
    id = Zoltan_Oct_id(oct);
    Zoltan_Oct_setID(oct, id+nprevoct);             /* now have global id's */
  }
  return 0;
}

/****************************************************************************/
/*
 * void Zoltan_Oct_visit_all_subtrees()
 *
 * visits each of the subtrees that are on the local processor
 */
static void Zoltan_Oct_visit_all_subtrees(ZZ *zz, float *part_sizes) {
  pRList  RootList;                           /* list of all local roots */
  pOctant RootOct;
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);

  /* get the list of all the local roots */
  /* iterate through each root in localroot list */ 
  /* mark each subtree as being Zoltan_Oct_visited */
  /* and free the costs */
  RootList = Zoltan_Oct_POct_localroots(OCT_info);
  while ((RootOct = RL_nextRootOctant(&RootList))) {
    Zoltan_Oct_visit(zz, RootOct, part_sizes);  
    Zoltan_Oct_costs_free(OCT_info, RootOct);
  }
}

/****************************************************************************/
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
static void Zoltan_Oct_visit(ZZ *zz, pOctant octant, float *part_sizes) {
  float cost;                /* Cost of this octant */
  float togo;                /* Remaining room in current partition */
  float behind;              /* How many to make up for from all prev parts */
  pOctant children[8];       /* children of the octant */
  int i;                     /* index counter */
  COORD origin;              /* center of the octant */
  double volume;             /* volume of the octant */
  double prod[3];            /* product of octant origin and its volume */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);

  DFS_Part_Count++;
  cost = Zoltan_Oct_costs_value(octant);      /* get the cost of the octant */
  /*behind = partition * optcost - total;*/     /* calcuate how much behind */
  behind = (tmpcost*globalcost) - total;        /* calcuate how much behind */
  
  if(0)
    fprintf(stderr,"LGG[%d] pc=%f, c=%f, ps=%f, b=%f\n", partition, pcost,
	    cost, optsize, behind);
  /* If octant does not overflow the current partition, then use it. */
  /*if( cost==0 || (pcost+cost) <= (optcost+behind)) {*/
  if(cost==0 || ((pcost+cost) <= (optsize+behind))) {
    Zoltan_Oct_tag_subtree(OCT_info,octant,partition);
    /*fprintf(stderr,"LGG[%d] pc=%f, c=%f, ps=%f, b=%f\n", partition, pcost,
	    cost, optsize, behind);*/
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
    Zoltan_Oct_modify_newpid(octant, partition);                 /* Nonterm */
    Zoltan_Oct_children(octant,children);

    for (i=0; i<8; i++)                    /* Simple - just visit in order */
      if(children[i] && Zoltan_Oct_POct_local(OCT_info, octant,i))
        Zoltan_Oct_visit(zz,children[i],part_sizes);
    return;
  }
  
  /* 
   * No suboctants!
   * We've hit bottom - have to decide whether to add to
   * the current partition or start a new one.
   */
  togo = behind + optsize - pcost;

  /*printf("proc=%d, part=%d, b=%f, pcost=%f, cost=%f, os=%f\n",
	 zz->Proc, partition, behind, pcost, cost, optsize);*/

  if ((cost-togo) >= togo) {
    /*printf("proc=%d, part=%d, togo=%f, pcost=%f, cost=%f, g=%f\n",
	   zz->Proc, partition, togo, pcost, cost, globalcost);*/
    /*
     * End current part and start new one. We are more "over" than "under"
     */
    tmpcost += part_sizes[partition];
    partition++;                               /* Move on to next partition */
    while((part_sizes[partition] == 0) && (partition < (zz->Num_Proc - 1)))
      partition++;
    optsize = part_sizes[partition]*globalcost;
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

/****************************************************************************/
/*
 * void Zoltan_Oct_tag_subtree(pOctant octant, int partition_number)
 *
 * marks all the octants within the subtree to be in the current partition
 */
static void Zoltan_Oct_tag_subtree(OCT_Global_Info *OCT_info,pOctant octant, 
				   int part) {
  pOctant children[8];                                /* children of octant */
  int i;                                              /* index counter */

  /* modify NPID so octant know where to migrate to */
  Zoltan_Oct_modify_newpid(octant, part);

  if (Zoltan_Oct_isTerminal(octant))
    return;

  /* if octant has children, have to tag them too */
  Zoltan_Oct_children(octant,children);
  
  for (i=0; i<8; i++)                       /* Simple - just visit in order */
    /* if (children[i] && Zoltan_Oct_local(OCT_info, children[i])) */
    if(children[i] && Zoltan_Oct_POct_local(OCT_info, octant,i))
      Zoltan_Oct_tag_subtree(OCT_info,children[i],part);
}

/****************************************************************************/
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
  pRList RootList;                           /* list of the local roots */
  pOctant oct;                               /* octree octant */
  pOctant *docts = NULL;                     /* array of octants being sent */
  int *dpids = NULL;                         /* array of octant pids */
  int dcount;                                /* count of octants being sent */
  int pid;                                   /* processor id */
  int nrecocts;
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);
  char *yo = "Zoltan_Oct_dfs_migrate";

  if(Zoltan_Oct_nOctants()) {  /* allocate space for octants being migrated */
    docts = (pOctant *)ZOLTAN_MALLOC(Zoltan_Oct_nOctants() * sizeof(pOctant));
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

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
