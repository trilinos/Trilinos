/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "octant_const.h"
#include "oct_util_const.h"

/*****************************************************************************/
static int OCT_count;                /* count of all local octants           */
static int OCT_idcount;              /* count for id's, help with uniqueness */
static pOctant LB_Oct_mallocremote();
static pOctant LB_Oct_malloc();
static int     LB_Oct_nlocal(pOctant oct);
/*****************************************************************************/
/*
 * pOctant LB_Oct_mallocremote()
 *
 * creates space in memory for a new octant
 */
static pOctant LB_Oct_mallocremote() {
  pOctant ptr;                                  /* pointer to the new octant */
  
  ptr=(pOctant) LB_MALLOC(sizeof(Octant));                /* allocate space */
  if (!ptr) 
    return NULL;

  return(ptr);
}
/*****************************************************************************/
/*
 * pOctant LB_Oct_malloc()
 *
 * creates space in memory for a new octant
 */
static pOctant LB_Oct_malloc() {
  pOctant ptr;                                  /* pointer to the new octant */
  
  ptr=(pOctant) LB_MALLOC(sizeof(Octant));                /* allocate space */
  if (!ptr) 
    return NULL;

  OCT_count++;                                      /* increase octant count */
  return(ptr);
}

void LB_Oct_initCounters() {
  OCT_count = 0;
  OCT_idcount = 0;
}

/*****************************************************************************/
/*
 * LB_Oct_newremote()
 *
 * create a new octant on the local processor and return
 * a pointer to it.  It will have no parents or children.
 */
pOctant LB_Oct_newremote() {
  pOctant new_ptr;                              /* pointer to the new octant */
  int i;                                        /* index counter */
  
  new_ptr=LB_Oct_mallocremote();                /* create space for the new octant */
  if(!new_ptr)
    return NULL;
  new_ptr->type = REMOTEOCT;
  /* Null out child pointers, and invalidate cpids */
  for (i=0; i<8; i++) {
    new_ptr->child[i] = NULL;
    new_ptr->cpid[i] = -1;
  }
  /* setup default information about octant */
  new_ptr->parent= NULL;
  new_ptr->remoteptr = NULL;
  new_ptr->ppid = -1;
  new_ptr->id = -1;
  new_ptr->which = -1;
  new_ptr->numChild = 0;
  new_ptr->list=NULL;
  new_ptr->cost = 0;
  new_ptr->dir = 0;
  new_ptr->mapidx = -1;
  new_ptr->area = 0;
  new_ptr->npid = -1;
  /* new_ptr->orientation = -1; */
  return(new_ptr);
}

/*****************************************************************************/
/*
 * LB_Oct_new()
 *
 * create a new octant on the local processor and return
 * a pointer to it.  It will have no parents or children.
 */
pOctant LB_Oct_new() {
  pOctant new_ptr;                              /* pointer to the new octant */
  int i;                                        /* index counter */
  
  new_ptr=LB_Oct_malloc();                   /* create space for the new octant */
  if(!new_ptr)
    return NULL;

  new_ptr->type = LOCALOCT;
  /* Null out child pointers, and invalidate cpids */
  for (i=0; i<8; i++) {
    new_ptr->child[i] = NULL;
    new_ptr->cpid[i] = -1;
  }
  /* setup default information about octant */
  new_ptr->parent= NULL;
  new_ptr->remoteptr = NULL;
  new_ptr->ppid = -1;
  new_ptr->id = (OCT_idcount++);
  new_ptr->which = -1;
  new_ptr->numChild = 0;
  new_ptr->list=NULL;
  new_ptr->cost = 0;
  new_ptr->dir = 0;
  new_ptr->mapidx = -1;
  new_ptr->area = 0;
  new_ptr->npid = -1;
  /* new_ptr->orientation = -1; */
  return(new_ptr);
}

void LB_Oct_free(pOctant oct) {
  oct->list = NULL;
  if(oct->type == LOCALOCT)
    OCT_count--;
  LB_FREE(&oct);
}

/*****************************************************************************/
/*#ifdef LGG_MIGOCT*/
/*
 * ATTN!! :
 * Manually setting the octant's id number is necessary if octants are 
 * being migrated this is not the case with the current implementation 
 * of the code references are left here in case it does become necessary
 * to do so again.
 */
/*
 * void LB_Oct_setID(pOctant octant, int id)
 *
 * sets the id of the octant
 */
void LB_Oct_setID(pOctant oct, int id) {
  oct->id=id;
}
/* #endif */ /* LGG_MIGOCT */ 

/*
 * void LB_Oct_setID(pOctant octant, int id)
 *
 * sets the id of the octant
 */
void LB_Oct_setMapIdx(pOctant oct, int idx) {
  oct->mapidx=idx;
}

void LB_Oct_setDir(pOctant oct, int dir) {
  oct->dir = dir;
}

int LB_Oct_mapidx(pOctant oct) {
  return oct->mapidx;
}

int LB_Oct_dir(pOctant oct) {
  return oct->dir;
}

/*****************************************************************************/
/*
 * void LB_Oct_id(pOctant octant)
 *
 * gets the id of the octant 
 */
int LB_Oct_id(pOctant oct) 
{ return(oct->id); }

/*****************************************************************************/
/*
 * void LB_Oct_setchildnum(pOctant octant, int childnumber)
 *
 * sets the child number of the octant
 */
void LB_Oct_setchildnum(pOctant oct, int childnum)
{ oct->which=childnum; }

/*****************************************************************************/
/*
 * void LB_Oct_setcildren(pOctant octant, pOctant children, int *cpids)
 *
 * sets the child pointers to point to the children
 */
void LB_Oct_setchildren(pOctant oct, pOctant *children, int *cpids) {
  int i;
  for(i = 0; i < 8; i++) {
    LB_Oct_setchild(oct, i, children[i]);
    LB_Oct_setCpid(oct, i, cpids[i]);
  }
}

/*****************************************************************************/
/*
 * void LB_Oct_setcild(pOctant octant, int childnumber, pOctant child)
 *
 * sets the ith child pointer to point to the child
 */
void LB_Oct_setchild(pOctant oct, int i, pOctant child) {
  oct->child[i]=child;         /* need to make sure child's info is correct */
  if( oct->numChild < 8 )
    oct->numChild++;
}

/*****************************************************************************/
/*
 * void LB_Oct_setbounds(pOctant octant, COORD minimum_bounds,
 *                    COORD maximum_bounds)
 *
 * sets the min and max bounds of an octant 
 */
void LB_Oct_setbounds(pOctant oct, COORD min, COORD max)
{ vector_set(oct->min,min); vector_set(oct->max,max); }

/*****************************************************************************/
/*
 * void LB_Oct_bounds(pOctant octant, COORD minimum, COORD maximum)
 *
 * gets the min and max bounds of an octant
 */
void LB_Oct_bounds(pOctant oct, COORD min, COORD max)
{ vector_set(min,oct->min); vector_set(max,oct->max); }

/*****************************************************************************/
/*
 * pOctant LB_Oct_parent(pOctant octant)
 *
 * returns a pointer to the parent of the octant
 */
pOctant LB_Oct_parent(pOctant oct) 
{ return(oct->parent); }

/*****************************************************************************/
/*
 * pOctant LB_Oct_child(pOctant octant, int child_index)
 *
 * return the ith child of oct
 */
pOctant LB_Oct_child(pOctant oct, int i) 
{ return(oct->child[i]); }

/*****************************************************************************/
/*
 * LB_Oct_cpids(pOctant octant, pOctant children[8])
 *
 * fill in values of all an octants children pids
 */
void LB_Oct_cpids(pOctant oct, int cpids[8]) {
  int i;                                                    /* index counter */
  
  for (i=0; i<8; i++)
    cpids[i] = oct->cpid[i];
}
  
/*****************************************************************************/
/*
 * LB_Oct_children(pOctant octant, pOctant children[8])
 *
 * fill in values of all an octants children
 */
int LB_Oct_children(pOctant oct, pOctant children[8]) {
  int i;                                                    /* index counter */
  
  for (i=0; i<8; i++)
    children[i] = LB_Oct_child(oct,i);
  return(oct->numChild);
}
  
/*****************************************************************************/
/*
 * int LB_Oct_isTerminal(pOctant octant)
 *
 * returns TRUE if the octant is terminal (has NO children)
 */
int LB_Oct_isTerminal(pOctant oct) {
  int i,                                         /* index counter */
      child;                                     /* flag identifying a child */

  child=0;
  for (i=0; (!child) && (i<8) ; i++)
    if (oct->child[i] != NULL)
      child=1;
  if(child == 0)
    return(1);
  else
    return(0);
}

/*****************************************************************************/
/*
 * pRegion LB_Oct_regionlist(pOctant octant)
 * get a copy of the octant's region list
 */
pRegion LB_Oct_regionlist(pOctant oct) {
  if (!LB_Oct_isTerminal(oct))
    abort();
  return(oct->list);
}

/*****************************************************************************/
/*
 * void LB_Oct_addRegion(pOctant octant, pRegion region)
 * add a region to oct's list
 */
int LB_Oct_addRegion(LB *lb, pOctant oct, pRegion region) { 
  char *yo = "LB_Oct_addRegion";
  pRegion entry;                      /* pointer to new entry in region list */

  if(oct == NULL) 
    return LB_WARN;

  entry = (pRegion) LB_MALLOC(sizeof(Region));   /* malloc space for region */
  if(entry == NULL) {
    LB_PRINT_ERROR(lb->Proc, yo, "Cannot allocated memory for region.");
    return LB_MEMERR;
  }

  entry->Global_ID = LB_MALLOC_GID(lb);
  entry->Local_ID = LB_MALLOC_LID(lb);
  /* copy region information into the entry */
  vector_set(entry->Coord, region->Coord);
  entry->Weight = region->Weight;
  LB_SET_GID(lb, entry->Global_ID, region->Global_ID);
  LB_SET_LID(lb, entry->Local_ID, region->Local_ID);
  entry->Proc = region->Proc;

  /* attach region to region list */
  entry->next = oct->list; 
  oct->list = entry;
  return LB_OK;
}

/*****************************************************************************/
/*
 * void LB_Oct_clearRegions(pOctant octant)
 * erase all of a oct's regions
 */
void LB_Oct_clearRegions(pOctant oct) { 
  pRegion ptr;                                         /* pointer to regions */
  
  ptr = oct->list;
  while(ptr != NULL) {
    oct->list = ptr->next;
    LB_FREE(&(ptr->Global_ID));
    LB_FREE(&(ptr->Local_ID));
    LB_FREE(&ptr);
    ptr = oct->list;
  }
  oct->list=NULL;
}

/*****************************************************************************/
/*
 * int LB_Oct_nRegions(pOctant octant)
 * return the number of regions in the octant's list
 */
int LB_Oct_nRegions(pOctant oct) {
  pRegion ptr;                     /* pointer to iterate through region list */
  int count;                       /* count of number of regions in list */

  if (oct == NULL) 
    return 0;

  if(!LB_Oct_isTerminal(oct)) 
    return 0;

  count = 0;
  ptr = oct->list;
  while(ptr != NULL) {
    count ++;
    ptr = ptr->next;
  }
  return(count);
}
/*****************************************************************************/
/***                        attached data routines                         ***/
/*
 * void LB_Oct_modify_cost(pOctant octant, float cost)
 *
 * modifies the cost field of the octant
 */
void LB_Oct_modify_cost(pOctant oct, float cost)
{
  oct->cost = cost;
}

/*****************************************************************************/
/*
 * void LB_Oct_modify_newpid(pOctant octant, int new_processor_id)
 *
 * modifies the npid field of the octant, telling where the octant should
 * migrate to
 */
void LB_Oct_modify_newpid(pOctant oct, int newpid)
{
  oct->npid = newpid;
}

/*****************************************************************************/
/*
 * int LB_Oct_data_newpid(pOctant octant)
 *
 * returns the new processor id of the octant 
 */
int LB_Oct_data_newpid(pOctant oct)
{ 
  return(oct->npid);
}

/*****************************************************************************/
/*
 * int LB_Oct_nlocal(pOctant octant)
 *
 * return the number of local leaves in the subtree
 */
static int LB_Oct_nlocal(pOctant oct) {
  int i;                                    /* index counter */
  pOctant child;                            /* child of an octant */
  int total;                                /* total number of local octants */

  if (LB_Oct_isTerminal(oct))
    return(1);
  total=0;
  for (i=0; i<8; i++) {
    child = LB_Oct_child(oct,i);
    if (child)
      total+=LB_Oct_nlocal(child);
  }
  return(total);
}

/*****************************************************************************/
/*
 * int LB_Oct_nOctant()
 *
 * returns the number of octants on local processor
 */
int LB_Oct_nOctants()
{ return (OCT_count); }

/*****************************************************************************/
/*
 * void LB_Oct_origin_volume(pOctant oct, COORD origin, double *volume)
 *
 * gets the origin and volume of the octant
 */
void LB_Oct_origin_volume(pOctant oct, COORD origin, double *volume) {
  COORD min,                                         /* octant minimum bound */
        max;                                         /* octant maximum bound */
  double size[3];                                    /* size of the octant */

  LB_Oct_bounds(oct,min,max);
  LB_bounds_to_origin_size(min,max,origin,size);  
  *volume=size[0]*size[1]*size[2];
}

int LB_Oct_Cpid(pOctant octant, int i)
{
  return octant->cpid[i];
}

int LB_Oct_Ppid(pOctant octant)
{
  return octant->ppid;
}

int LB_Oct_childnum(pOctant octant)
{
  return octant->which;
}

/*****************************************************************************/
/*
 * void LB_Oct_setCpid(pOctant octant, int child_index, int child_processor_id)
 *
 * sets the child processor id field of an octant
 */
void LB_Oct_setCpid(pOctant octant, int i, int cpid) {
  octant->cpid[i] = cpid;
}
