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
static pOctant Zoltan_Oct_mallocremote();
static pOctant Zoltan_Oct_malloc();
static int     Zoltan_Oct_nlocal(pOctant oct);
/*****************************************************************************/
/*
 * pOctant Zoltan_Oct_mallocremote()
 *
 * creates space in memory for a new octant
 */
static pOctant Zoltan_Oct_mallocremote() {
  pOctant ptr;                                  /* pointer to the new octant */
  
  ptr=(pOctant) ZOLTAN_MALLOC(sizeof(Octant));                /* allocate space */
  if (!ptr) 
    return NULL;

  return(ptr);
}
/*****************************************************************************/
/*
 * pOctant Zoltan_Oct_malloc()
 *
 * creates space in memory for a new octant
 */
static pOctant Zoltan_Oct_malloc() {
  pOctant ptr;                                  /* pointer to the new octant */
  
  ptr=(pOctant) ZOLTAN_MALLOC(sizeof(Octant));                /* allocate space */
  if (!ptr) 
    return NULL;

  OCT_count++;                                      /* increase octant count */
  return(ptr);
}

void Zoltan_Oct_initCounters() {
  OCT_count = 0;
  OCT_idcount = 0;
}

/*****************************************************************************/
/*
 * Zoltan_Oct_newremote()
 *
 * create a new octant on the local processor and return
 * a pointer to it.  It will have no parents or children.
 */
pOctant Zoltan_Oct_newremote() {
  pOctant new_ptr;                              /* pointer to the new octant */
  int i;                                        /* index counter */
  
  new_ptr=Zoltan_Oct_mallocremote();                /* create space for the new octant */
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
 * Zoltan_Oct_new()
 *
 * create a new octant on the local processor and return
 * a pointer to it.  It will have no parents or children.
 */
pOctant Zoltan_Oct_new() {
  pOctant new_ptr;                              /* pointer to the new octant */
  int i;                                        /* index counter */
  
  new_ptr=Zoltan_Oct_malloc();                   /* create space for the new octant */
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

/* KDDKDDFREE changed oct to *oct to allow NULL ptr from ZOLTAN_FREE to 
 * KDDKDDFREE propagate back to the calling routine. */
void Zoltan_Oct_free(OCT_Global_Info *OCT_info, pOctant *oct) {
/* KDDKDDFREE Bookkeeping -- set parent's child pointers and children's
 * KDDKDDFREE parent pointers to NULL when delete an octant.  Following
 * KDDKDDFREE pointers that are no longer valid but were not set to NULL
 * KDDKDDFREE caused many problems on the DEC cluster.  1/2001. */
int i;
pOctant parent;


  /* Set parent's child pointer to NULL */
  parent = (*oct)->parent;
  if (parent != NULL && (*oct)->ppid == OCT_info->OCT_localpid) {
    i = Zoltan_Oct_childnum(*oct);
    parent->child[i] = NULL;
    parent->cpid[i] = -1;
  }

  /* Set children't parent pointer to NULL. */
  for (i = 0; i < 8; i++) 
    if ((*oct)->child[i] != NULL && (*oct)->cpid[i] == OCT_info->OCT_localpid) {
      (*oct)->child[i]->parent = NULL;
      (*oct)->child[i]->ppid = -1;
    }
/* END KDDKDDFREE */

  (*oct)->list = NULL;
  if((*oct)->type == LOCALOCT)
    OCT_count--;
  ZOLTAN_FREE(oct);
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
 * void Zoltan_Oct_setID(pOctant octant, int id)
 *
 * sets the id of the octant
 */
void Zoltan_Oct_setID(pOctant oct, int id) {
  oct->id=id;
}
/* #endif */ /* LGG_MIGOCT */ 

/*
 * void Zoltan_Oct_setID(pOctant octant, int id)
 *
 * sets the id of the octant
 */
void Zoltan_Oct_setMapIdx(pOctant oct, int idx) {
  oct->mapidx=idx;
}

void Zoltan_Oct_setDir(pOctant oct, int dir) {
  oct->dir = dir;
}

int Zoltan_Oct_mapidx(pOctant oct) {
  return oct->mapidx;
}

int Zoltan_Oct_dir(pOctant oct) {
  return oct->dir;
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_id(pOctant octant)
 *
 * gets the id of the octant 
 */
int Zoltan_Oct_id(pOctant oct) 
{ return(oct->id); }

/*****************************************************************************/
/*
 * void Zoltan_Oct_setchildnum(pOctant octant, int childnumber)
 *
 * sets the child number of the octant
 */
void Zoltan_Oct_setchildnum(pOctant oct, int childnum)
{ oct->which=childnum; }

/*****************************************************************************/
/*
 * void Zoltan_Oct_setcildren(pOctant octant, pOctant children, int *cpids)
 *
 * sets the child pointers to point to the children
 */
void Zoltan_Oct_setchildren(pOctant oct, pOctant *children, int *cpids) {
  int i;
  for(i = 0; i < 8; i++) {
    Zoltan_Oct_setchild(oct, i, children[i]);
    Zoltan_Oct_setCpid(oct, i, cpids[i]);
  }
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_setcild(pOctant octant, int childnumber, pOctant child)
 *
 * sets the ith child pointer to point to the child
 */
void Zoltan_Oct_setchild(pOctant oct, int i, pOctant child) {
  oct->child[i]=child;         /* need to make sure child's info is correct */
  if( oct->numChild < 8 )
    oct->numChild++;
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_setbounds(pOctant octant, COORD minimum_bounds,
 *                    COORD maximum_bounds)
 *
 * sets the min and max bounds of an octant 
 */
void Zoltan_Oct_setbounds(pOctant oct, COORD min, COORD max)
{ vector_set(oct->min,min); vector_set(oct->max,max); }

/*****************************************************************************/
/*
 * void Zoltan_Oct_bounds(pOctant octant, COORD minimum, COORD maximum)
 *
 * gets the min and max bounds of an octant
 */
void Zoltan_Oct_bounds(pOctant oct, COORD min, COORD max)
{ vector_set(min,oct->min); vector_set(max,oct->max); }

/*****************************************************************************/
/*
 * pOctant Zoltan_Oct_parent(pOctant octant)
 *
 * returns a pointer to the parent of the octant
 */
pOctant Zoltan_Oct_parent(pOctant oct) 
{ return(oct->parent); }

/*****************************************************************************/
/*
 * pOctant Zoltan_Oct_child(pOctant octant, int child_index)
 *
 * return the ith child of oct
 */
pOctant Zoltan_Oct_child(pOctant oct, int i) 
{ return(oct->child[i]); }

/*****************************************************************************/
/*
 * Zoltan_Oct_cpids(pOctant octant, pOctant children[8])
 *
 * fill in values of all an octants children pids
 */
void Zoltan_Oct_cpids(pOctant oct, int cpids[8]) {
  int i;                                                    /* index counter */
  
  for (i=0; i<8; i++)
    cpids[i] = oct->cpid[i];
}
  
/*****************************************************************************/
/*
 * Zoltan_Oct_children(pOctant octant, pOctant children[8])
 *
 * fill in values of all an octants children
 */
int Zoltan_Oct_children(pOctant oct, pOctant children[8]) {
  int i;                                                    /* index counter */
  
  for (i=0; i<8; i++)
    children[i] = Zoltan_Oct_child(oct,i);
  return(oct->numChild);
}
  
/*****************************************************************************/
/*
 * int Zoltan_Oct_isTerminal(pOctant octant)
 *
 * returns TRUE if the octant is terminal (has NO children)
 */
int Zoltan_Oct_isTerminal(pOctant oct) {
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
 * pRegion Zoltan_Oct_regionlist(pOctant octant)
 * get a copy of the octant's region list
 */
pRegion Zoltan_Oct_regionlist(pOctant oct) {
  if (!Zoltan_Oct_isTerminal(oct))
    abort();
  return(oct->list);
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_addRegion(pOctant octant, pRegion region)
 * add a region to oct's list
 */
int Zoltan_Oct_addRegion(ZZ *zz, pOctant oct, pRegion region) { 
  char *yo = "Zoltan_Oct_addRegion";
  pRegion entry;                      /* pointer to new entry in region list */

  if(oct == NULL) 
    return ZOLTAN_WARN;

  entry = (pRegion) ZOLTAN_MALLOC(sizeof(Region));   /* malloc space for region */
  if(entry == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Cannot allocated memory for region.");
    return ZOLTAN_MEMERR;
  }

  entry->Global_ID = ZOLTAN_MALLOC_GID(zz);
  entry->Local_ID = ZOLTAN_MALLOC_LID(zz);
  /* copy region information into the entry */
  vector_set(entry->Coord, region->Coord);
  entry->Weight = region->Weight;
  ZOLTAN_SET_GID(zz, entry->Global_ID, region->Global_ID);
  ZOLTAN_SET_LID(zz, entry->Local_ID, region->Local_ID);
  entry->Proc = region->Proc;

  /* attach region to region list */
  entry->next = oct->list; 
  oct->list = entry;
  return ZOLTAN_OK;
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_clearRegions(pOctant octant)
 * erase all of a oct's regions
 */
void Zoltan_Oct_clearRegions(pOctant oct) { 
  pRegion ptr;                                         /* pointer to regions */
  
  ptr = oct->list;
  while(ptr != NULL) {
    oct->list = ptr->next;
    ZOLTAN_FREE(&(ptr->Global_ID));
    ZOLTAN_FREE(&(ptr->Local_ID));
    ZOLTAN_FREE(&ptr);
    ptr = oct->list;
  }
  oct->list=NULL;
}

/*****************************************************************************/
/*
 * int Zoltan_Oct_nRegions(pOctant octant)
 * return the number of regions in the octant's list
 */
int Zoltan_Oct_nRegions(pOctant oct) {
  pRegion ptr;                     /* pointer to iterate through region list */
  int count;                       /* count of number of regions in list */

  if (oct == NULL) 
    return 0;

  if(!Zoltan_Oct_isTerminal(oct)) 
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
 * void Zoltan_Oct_modify_cost(pOctant octant, float cost)
 *
 * modifies the cost field of the octant
 */
void Zoltan_Oct_modify_cost(pOctant oct, float cost)
{
  oct->cost = cost;
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_modify_newpid(pOctant octant, int new_processor_id)
 *
 * modifies the npid field of the octant, telling where the octant should
 * migrate to
 */
void Zoltan_Oct_modify_newpid(pOctant oct, int newpid)
{
  oct->npid = newpid;
}

/*****************************************************************************/
/*
 * int Zoltan_Oct_data_newpid(pOctant octant)
 *
 * returns the new processor id of the octant 
 */
int Zoltan_Oct_data_newpid(pOctant oct)
{ 
  return(oct->npid);
}

/*****************************************************************************/
/*
 * int Zoltan_Oct_nlocal(pOctant octant)
 *
 * return the number of local leaves in the subtree
 */
static int Zoltan_Oct_nlocal(pOctant oct) {
  int i;                                    /* index counter */
  pOctant child;                            /* child of an octant */
  int total;                                /* total number of local octants */

  if (Zoltan_Oct_isTerminal(oct))
    return(1);
  total=0;
  for (i=0; i<8; i++) {
    child = Zoltan_Oct_child(oct,i);
    if (child)
      total+=Zoltan_Oct_nlocal(child);
  }
  return(total);
}

/*****************************************************************************/
/*
 * int Zoltan_Oct_nOctant()
 *
 * returns the number of octants on local processor
 */
int Zoltan_Oct_nOctants()
{ return (OCT_count); }

/*****************************************************************************/
/*
 * void Zoltan_Oct_origin_volume(pOctant oct, COORD origin, double *volume)
 *
 * gets the origin and volume of the octant
 */
void Zoltan_Oct_origin_volume(pOctant oct, COORD origin, double *volume) {
  COORD min,                                         /* octant minimum bound */
        max;                                         /* octant maximum bound */
  double size[3];                                    /* size of the octant */

  Zoltan_Oct_bounds(oct,min,max);
  Zoltan_Oct_bounds_to_origin_size(min,max,origin,size);  
  *volume=size[0]*size[1]*size[2];
}

int Zoltan_Oct_Cpid(pOctant octant, int i)
{
  return octant->cpid[i];
}

int Zoltan_Oct_Ppid(pOctant octant)
{
  return octant->ppid;
}

int Zoltan_Oct_childnum(pOctant octant)
{
  return octant->which;
}

/*****************************************************************************/
/*
 * void Zoltan_Oct_setCpid(pOctant octant, int child_index, int child_processor_id)
 *
 * sets the child processor id field of an octant
 */
void Zoltan_Oct_setCpid(pOctant octant, int i, int cpid) {
  octant->cpid[i] = cpid;
}
