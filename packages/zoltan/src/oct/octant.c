#include<stdlib.h>
#include<stdio.h>
#include "octant_const.h"
#include "util_const.h"

/* static int dimensions = 8; */

/* WARNING: GLOBAL VARIABLES... BE CAREFUL WHEN USING */
pRList OCT_rootlist;          /* list of all the local roots          */
int OCT_count;                /* count of all local octants           */
int OCT_localpid;             /* the processor id                     */
int idcount;                  /* count for id's, help with uniqueness */
double gmin[3];               /* global root's min bounds */
double gmax[3];               /* global root's max bounds */
int dim = 8;
/* Map *array; */
/* static int msg_nprocs; */         /* total number of processors available */

void POC_init(int pid) {
  OCT_count=0;
  OCT_localpid=pid;
  OCT_rootlist=NULL;
  idcount=0;
}

pOctant POC_malloc() {
  pOctant ptr;
  
  ptr=(pOctant)malloc(sizeof(Octant));
  if (!ptr) {
    perror("POC_malloc: malloc failure");
    abort();
  }
  OCT_count++;
  return(ptr);
}

/*
 * POC_new()
 *
 * allocate a new Paroct on the local processor and return
 * a pointer to it.  It will have no parents or children.
 */
pOctant POC_new() {
  pOctant new;
  int i;
  
  new=POC_malloc();
  for (i=0; i<8; i++) {
    new->child[i] = NULL;
    new->cpid[i] = -1;
  }
  new->parent= NULL;
  new->ppid = OCT_localpid;
  new->id=(idcount++);
  new->which = -1;
  new->numChild = 0;
  new->list=NULL;
  new->cost = 0;
  new->npid = OCT_localpid;
  return(new);
}

/* 
 * this does not delete the attached regions, must specifically call
 * POC_clearRegions 
 */
void POC_free(pOctant oct) {
  pRList ptr, prev;
  pRegion c;
  
  ptr = OCT_rootlist;
  prev = NULL;
  if (oct->ppid != OCT_localpid) {
    while(ptr != NULL) {
      if(ptr->oct->id == oct->id) {
	if(ptr == OCT_rootlist)
	  OCT_rootlist = ptr->next;
	else
	  prev->next = ptr->next;
	ptr->oct = NULL;
	free(ptr);
	ptr = NULL;
      }
      else {
	prev = ptr;
	ptr = ptr->next;
      }
    }
  }
/********  
  c = oct->list;
  while(c != NULL) {
    oct->list = c->next;
    free(c);
    c = oct->list;
  }
********/
  oct->list = NULL;
  OCT_count--;
  free(oct);
}

void POC_setID(pOctant oct, int id) 
{ oct->id=id; }

int POC_id(pOctant oct) 
{ return(oct->id); }

void POC_setparent(pOctant oct, pOctant parent, int ppid) {
  pRList tmp,                          /* temp variable used for iterating */
         prev;                         /* pointer to previous root looked at */

  if (oct->ppid == OCT_localpid) {                              /* was local */
    if (ppid != OCT_localpid) {                              /* now nonlocal */
      tmp = (pRList)malloc(sizeof(RList));
      if(tmp == NULL) {
	fprintf(stderr, "%s: %s\n", "POC_setparent",
		"Could not allocalte new rootlist entry.");
	abort();
      }
      tmp->oct = oct;
      prev = OCT_rootlist;
      if(prev != NULL) {
	while(prev->next != NULL)
	  prev = prev->next;
	tmp->next = NULL;
	prev->next = tmp;
      }
      else {
	tmp->next = OCT_rootlist;
	OCT_rootlist = tmp;
      }
    }
  }
  else                                                        /* was foreign */
    if (ppid == OCT_localpid) {                               /* now local   */
      tmp = OCT_rootlist; 
      while(tmp != NULL) {
	if(tmp->oct->id == oct->id) {
	  if(tmp == OCT_rootlist)
	    OCT_rootlist = tmp->next;
	  else
	    prev->next = tmp->next;
	  tmp->oct = NULL;
	  free(tmp);
	  tmp = NULL;
	}
	else {
	  prev = tmp;
	  tmp = tmp->next;
	}
      }
      /* PList_remItem(OCT_rootlist,oct); */
    }
  oct->ppid=ppid;
  oct->parent=parent;
}

void POC_setchildnum(pOctant oct, int childnum)
{ oct->which=childnum; }

int POC_childnum(pOctant oct)
{ return(oct->which); }

void POC_setchild(pOctant oct, int i, pOctant child) {
  oct->child[i]=child;         /* need to make sure child's info is correct */
  oct->numChild++;
}
void POC_setchildren(pOctant oct, pOctant children[8], int cpids[8]) {
  int i;                                                    /* index counter */
  
  for (i=0; i<8; i++)
    POC_setchild(oct,i,children[i]);
}

void POC_setbounds(pOctant oct, double min[3], double max[3])
{ vector_set(oct->min,min); vector_set(oct->max,max); }

void POC_bounds(pOctant oct, double min[3], double max[3])
{ vector_set(min,oct->min); vector_set(max,oct->max); }

pOctant POC_parent(pOctant oct) 
{ return(oct->parent); }

/*
 * POC_child(oct,i,child,cpid)
 * return the ith child of oct, and it's processor id cpid.
 */
pOctant POC_child(pOctant oct, int i) 
{ return(oct->child[i]); }

/*
 * POC_children(oct,children,cpids)
 * fill in values of all an octants children
 */
int POC_children(pOctant oct, pOctant children[8]) {
  int i;                                                    /* index counter */
  
  for (i=0; i<8; i++)
    children[i] = POC_child(oct,i);
  return(oct->numChild);
}
  
/*
 * POC_isTerminal(oct)
 * returns TRUE if the octant is terminal (has NO children)
 */
int POC_isTerminal(pOctant oct) {
  int i,                                         /* index counter */
      child;                                     /* flag identifying a child */

  child=0;
  for (i=0; (!child) && (i<8) ; i++)
    child=(int)oct->child[i];
  if(child == 0)
    return(1);
  else
    return(0);
}

/*
 * POC_regionlist(oct)
 * get a copy of the octant's region list
 * NOTE: caller must free with PList_delete() when done
 */
pRegion POC_regionlist(pOctant oct) {
  if (!POC_isTerminal(oct))
    abort();
  return(oct->list);
}

/*
 * POC_addRegion(oct,region)
 * add a region to oct's list
 */
void POC_addRegion(pOctant oct, pRegion region) { 
  pRegion entry;

  entry = (pRegion)malloc(sizeof(Region));
  vector_set(entry->Coord, region->Coord);
  entry->Weight = region->Weight;
  entry->Tag.Global_ID = region->Tag.Global_ID;
  entry->Tag.Local_ID = region->Tag.Local_ID;
  entry->Tag.Proc = region->Tag.Proc;

  entry->next = oct->list; 
  oct->list = entry;
}

/*
 * POC_remRegion(oct,region)
 * remove a region from an oct's list
 */
void POC_remRegion(pOctant oct, pRegion region) { 
  pRegion tmp,                       /* temp var used for iterating */
          prev;                      /* pointer to previous region looked at */
  
  tmp = oct->list;
  while(tmp != NULL) {
    if(tmp->Tag.Global_ID == region->Tag.Global_ID) {
      if(tmp == oct->list)
	oct->list = tmp->next;
      else
	prev->next = tmp->next;
      free(tmp);
      tmp = NULL;
    }
    else {
      prev = tmp;
      tmp = tmp->next;
    }
  }
}

/*
 * POC_clearRegions(oct)
 * erase all of a oct's regions
 */
void POC_clearRegions(pOctant oct) { 
  pRegion ptr;                                         /* pointer to regions */
  
  ptr = oct->list;
  while(ptr != NULL) {
    oct->list = ptr->next;
    free(ptr);
    ptr = oct->list;
  }
  oct->list=NULL;
}

/*
 * POC_nRegions(oct)
 * return the number of regions in the octant's list
 */
int POC_nRegions(pOctant oct) {
  pRegion ptr;
  int count;

  if (!POC_isTerminal(oct))
    abort();
  count = 0;
  ptr = oct->list;
  while(ptr != NULL) {
    count ++;
    ptr = ptr->next;
  }
  return(count);
}

/*
 * POC_localroots()
 * return the list of local roots
 * caller must free space with PList_delete() when done
 */
pRList POC_localroots()
{ return(OCT_rootlist); }

/*
 * attached data routines 
 */
void POC_modify_cost(pOctant oct, float cost)
{ oct->cost = cost; }

void POC_modify_newpid(pOctant oct, int newpid) { 
  oct->npid = newpid; 
}

float POC_data_cost(pOctant oct)
{ return(oct->cost); }

int POC_data_newpid(pOctant oct)
{ return(oct->npid); }

/*
 * POC_nlocal(oct)
 * return the number of local leaves in the subtree
 */
int POC_nlocal(pOctant oct) {
  int i;                                    /* index counter */
  pOctant child;                            /* child of an octant */
  int total;                                /* total number of local octants */

  if (POC_isTerminal(oct))
    return(1);
  total=0;
  for (i=0; i<8; i++) {
    child = POC_child(oct,i);
    if (child)
      total+=POC_nlocal(child);
  }
  return(total);
}

int POC_nOctants()
{ return (OCT_count); }

/*
 * void POC_origin_volume(pOctant oct, Coord origin, double *volume)
 */
void POC_origin_volume(pOctant oct, double origin[3], double *volume) {
  double min[3],                                     /* octant minimum bound */
         max[3];                                     /* octant maximum bound */
  double size[3];                                    /* size of the octant */

  POC_bounds(oct,min,max);
  bounds_to_origin_size(min,max,origin,size);  
  *volume=size[0]*size[1]*size[2];
}

void POC_printResults() {
  pRList ptr;

  ptr = OCT_rootlist;
  while(ptr != NULL) {
    POC_DfsTraversal(ptr->oct);
    ptr = ptr->next;
  }
}

void POC_DfsTraversal(pOctant oct) {
  int i;

  if(POC_isTerminal(oct))
    POC_printRegionInfo(oct);
  else
    for(i=0; i<8; i++)
      POC_DfsTraversal(oct->child[i]);
}

void POC_printRegionInfo(pOctant oct) {
  pRegion ptr;
  
  if(!POC_isTerminal(oct)) {
    fprintf(stderr, "WARNING: printing regions in a non-terminal octant.\n");
    return;
  }
  printf("(%d) ocant %d with bounds min=%lf, %lf, %lf, & max %lf, %lf, %lf\n",
	 OCT_localpid, oct->id,
	 oct->min[0], oct->min[1], oct->min[2],
	 oct->max[0], oct->max[1], oct->max[2]);
  printf("from %d to %d\n", OCT_localpid, oct->npid);
  ptr = oct->list;
  if(ptr==NULL)
    printf("\tEMPTY\n");
  while(ptr != NULL) {
    printf("%d\t%d %d %d %d coord:(%lf, %lf, %lf)\n", ptr->attached,
	   ptr->Tag.Global_ID, ptr->Tag.Local_ID, ptr->Tag.Proc, oct->npid,
	   ptr->Coord[0], ptr->Coord[1], ptr->Coord[2]);
    ptr = ptr->next;
  }
  printf("\n\n");
}

pOctant POC_nextDfs(pOctant octant) {
  pOctant parent,
          child;
  int ppid,
      childnum,
      cpid;
  int i;

  if (!octant)
    return(NULL);

  for (i=0; i<8; i++) {
    child = POC_child(octant,i);
    if (child)
      return(child);          /* Go down */
  }

  parent = POC_parent(octant);
  while (parent) {
    for (i=octant->which; i<7; i++) {
      child = POC_child(parent,i+1);
      if (child)
	return(child);
    }

    octant=parent;                                       /* Go up */
    parent = POC_parent(octant);
  }
  
  return(NULL);         /* no more octants remain in dfs ordering */
}

int POC_local(pOctant octant, int i) {
  if(octant->cpid[i] == OCT_localpid)
    return 1;
  else {
    if(octant->cpid[i] == -1)
      fprintf(stderr, "WARNING: cpid was not set, returning 0 as default.\n");
    return 0;
  }
}

void POC_setCpid(pOctant octant, int i, int cpid) {
  octant->cpid[i] = cpid;
}

int POC_delTree(pOctant root) {
  int i;
  pOctant child;

  if(POC_isTerminal(root))
    POC_free(root);
  else {
    for(i=0; i<8; i++) {
      child = POC_child(root, i);
      if(child)
	POC_delTree(child);
    }
    POC_free(root);
  }
  return 1;
}
