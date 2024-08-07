// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "zz_util_const.h"
#include "phg_tree.h"
#include "phg.h"
#include <limits.h>


#define SET_MIN_NODE(ptr, offset, val) (ptr)[2*(offset)]=-(val)
#define SET_MAX_NODE(ptr, offset, val) (ptr)[2*(offset)+1]=(val)
#define IS_EMPTY(interval) ((interval)[1] == -1)
#define IS_INCLUDED(interval1,interval2) ((((interval2)[1] <= (interval1)[1]))&&((interval2)[1] <= (interval1)[1]))


/* Internal use only */
void* Zoltan_PHG_LB_Data_alloc();
int Zoltan_PHG_Timers_copy(ZZ* zz, struct phg_timer_indices* ftimers);
int Zoltan_PHG_Tree_copy(ZZ* zz, Zoltan_PHG_Tree* ftimers);
int Zoltan_PHG_Tree_init(Zoltan_PHG_Tree* ftimers);


void*
Zoltan_PHG_LB_Data_alloc()
{
  Zoltan_PHG_LB_Data * ptr;

  ptr = (Zoltan_PHG_LB_Data*) ZOLTAN_MALLOC(sizeof(Zoltan_PHG_LB_Data));
  if (ptr == NULL)
    return NULL;

  ptr->timers = NULL;
  ptr->tree = NULL;
  return (ptr);
}

struct phg_timer_indices *
Zoltan_PHG_LB_Data_timers(ZZ const * zz)
{
  Zoltan_PHG_LB_Data* ptr = (Zoltan_PHG_LB_Data*)zz->LB.Data_Structure;

  if (!ptr) return NULL;
  return (ptr->timers);
}

Zoltan_PHG_Tree *
Zoltan_PHG_LB_Data_tree(ZZ const * zz)
{
  Zoltan_PHG_LB_Data* ptr = (Zoltan_PHG_LB_Data*)zz->LB.Data_Structure;
  if (!ptr) return NULL;
  return (ptr->tree);
}

void
Zoltan_PHG_LB_Data_free_timers(ZZ* zz)
{
  Zoltan_PHG_LB_Data * ptr = zz->LB.Data_Structure;
  if (ptr == NULL)
    return;
  ZOLTAN_FREE(&(ptr->timers));
}

void
Zoltan_PHG_LB_Data_free_tree(ZZ* zz)
{
  Zoltan_PHG_LB_Data * ptr = zz->LB.Data_Structure;
  if ((ptr == NULL) || (ptr->tree == NULL))
    return;
  ptr->tree->array += 2;
  ZOLTAN_FREE(&(ptr->tree->array));
  ZOLTAN_FREE(&(ptr->tree));
}

/*****************************************************************************/

void Zoltan_PHG_Free_Structure(ZZ *zz)
{
  /* frees all data associated with LB.Data_Structure for hypergraphs */

  Zoltan_PHG_LB_Data_free_timers(zz);
  Zoltan_PHG_LB_Data_free_tree(zz);
  ZOLTAN_FREE(&(zz->LB.Data_Structure));
}

/*****************************************************************************/

int Zoltan_PHG_Copy_Structure(ZZ *to, ZZ const *from)
{
  /* Copies all data associated with LB.Data_Structure for hypergraphs */
  /* For now, that is just timing indices */
  if (from->LB.Data_Structure) {
    int ierr = ZOLTAN_OK;

    Zoltan_PHG_Free_Structure(to);
    ierr = Zoltan_PHG_Timers_copy(to, Zoltan_PHG_LB_Data_timers(from));
    if (ierr != ZOLTAN_OK)
      return (ierr);
    ierr = Zoltan_PHG_Tree_copy(to, Zoltan_PHG_LB_Data_tree(from));

  }
  return ZOLTAN_OK;
}



  /* The tree is a list of couple (-min, max) but only declare as an array of int */

/* Create tree structure */
int
Zoltan_PHG_Tree_create(int part_number, ZZ* zz)
{
  int part2; /* Power of 2 parts */
  Zoltan_PHG_Tree *tree;
  Zoltan_PHG_LB_Data * ptr = zz->LB.Data_Structure;

  if (ptr == NULL) {
    ptr = Zoltan_PHG_LB_Data_alloc();
    zz->LB.Data_Structure = (void*)ptr;  /* TODO simpleGRAPH test this does not get freed */
  }
  if (ptr == NULL)
    return (ZOLTAN_MEMERR);

  if (ptr->tree != NULL) {
    Zoltan_PHG_LB_Data_free_tree(zz);
  }

  tree = ptr->tree = (Zoltan_PHG_Tree*) ZOLTAN_MALLOC(sizeof(Zoltan_PHG_Tree));
  if (ptr->tree == NULL)
    return (ZOLTAN_MEMERR);
  if (part_number == 0)
    return ZOLTAN_OK;

  /* Round up to the next highest power of 2 */
  /* From http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2 */
  part2 = part_number;
  part2--;
  part2 |= part2 >> 1;
  part2 |= part2 >> 2;
  part2 |= part2 >> 4;
  part2 |= part2 >> 8;
  part2 |= part2 >> 16;

  part2++;

  tree->size = 2*part2-1;
  return (Zoltan_PHG_Tree_init(tree));
}


Zoltan_PHG_Tree *
get_tree(ZZ* zz)
{
  Zoltan_PHG_LB_Data * ptr = zz->LB.Data_Structure;

  return (ptr->tree);
}


int
Zoltan_PHG_Tree_init(Zoltan_PHG_Tree* tree)
{
  int i;

  tree->array = (int*) ZOLTAN_MALLOC(sizeof(int)*2*(tree->size));
  if (tree->array == NULL)
    return ZOLTAN_MEMERR;
  /* TRICK: we store -low, thus we can use MPI_MAX for low and high */
  for (i = 0 ; i < tree->size ; ++i) {
    SET_MIN_NODE(tree->array, i, tree->size + 1);
    SET_MAX_NODE(tree->array, i, -1);
  }

  tree->array -= 2; /* Begin at offset 1 */

  return ZOLTAN_OK;
}

int
Zoltan_PHG_Tree_copy(ZZ* zz, Zoltan_PHG_Tree* ftree)
{
  int ierr;
  Zoltan_PHG_LB_Data * ptr = zz->LB.Data_Structure;
  Zoltan_PHG_Tree *ttree;

  if (ptr == NULL)
    return (ZOLTAN_OK);

  Zoltan_PHG_LB_Data_free_tree(zz);

  if (ftree == NULL)
    return (ZOLTAN_OK);

  ttree = ptr->tree = (Zoltan_PHG_Tree*) ZOLTAN_MALLOC(sizeof(Zoltan_PHG_Tree));
  if (ptr->tree == NULL)
    return (ZOLTAN_MEMERR);

  ttree->size = ftree->size;
  ierr = Zoltan_PHG_Tree_init(ttree);

  if (ierr != ZOLTAN_OK)
    return (ierr);
  memcpy (ttree->array+2, ftree->array+2, 2*ttree->size*sizeof(int));
  return (ZOLTAN_OK);
}


/* Build a centralized tree */
int
Zoltan_PHG_Tree_centralize(ZZ *zz)
{
  Zoltan_PHG_Tree* tree;

  tree = Zoltan_PHG_LB_Data_tree(zz);
  if (tree == NULL)
    return ZOLTAN_OK;

  Zoltan_AllReduceInPlace(tree->array + 2, 2*tree->size , MPI_INT, MPI_MAX, zz->Communicator);

  /* TRICK: we store -low, thus we can use MPI_MAX for low and high */
  return ZOLTAN_OK;
}

void
Zoltan_PHG_Tree_Set(ZZ* zz, int father, int lo, int hi)
{
  Zoltan_PHG_Tree* tree;

  tree = Zoltan_PHG_LB_Data_tree(zz);
  if (tree == NULL)
    return;
  SET_MIN_NODE(tree->array, father, lo);
  SET_MAX_NODE(tree->array, father, hi);
}


int
find_interval_in_tree(Zoltan_PHG_Tree *tree, int *interval)
{
  int node = 1;
  int tree_size;

  if (IS_EMPTY(interval))
    return (-1);

  tree_size = get_tree_size(tree);
  if (-interval[0] == interval[1]) /* a leaf */
    return (interval[1]+(tree_size+1)/2);

  while (2*node <= tree_size) {
    if (!(IS_INCLUDED(&(tree->array[2*node]), interval))) {
      return (node/2);
    }
    node *= 2;
    if (-interval[0] > tree->array[2*node+1]) /* We need to search in the highest interval */
      node ++;
  }

  return (node/2);
}

static int
numerote(Zoltan_PHG_Tree *tree, int node, int part, int *partnumber)
{
  int partnum = part;
  if (2*node <= get_tree_size(tree)) { /* Not a leaf */
    partnum = numerote(tree, 2*node, part, partnumber);
    partnum = numerote(tree, 2*node+1, partnum, partnumber);
  }
  partnumber[node] = partnum;

  if (!IS_EMPTY(&(tree->array[2*node])))
    partnum++;
  return (partnum);
}

int *
compute_part_number(Zoltan_PHG_Tree *tree)
{
  int *partnumber;
  int tree_size;

  tree_size = get_tree_size(tree);
  partnumber = (int*)ZOLTAN_CALLOC(tree_size, sizeof(int));
  partnumber -= 1;
  numerote(tree, 1, 0, partnumber);

  return (partnumber);
}



int
Zoltan_PHG_Timers_copy(ZZ* zz, struct phg_timer_indices* ftimers)
{
  int ierr;
  struct phg_timer_indices *ttimer;

  Zoltan_PHG_LB_Data_free_timers(zz);
  if (ftimers == NULL)
    return (ZOLTAN_OK);
  ierr = Zoltan_PHG_Timers_init(zz);
  if (ierr != ZOLTAN_OK)
    return (ierr);
  ttimer = Zoltan_PHG_LB_Data_timers(zz);
  memcpy(ttimer, ftimers, sizeof(struct phg_timer_indices));
  return (ZOLTAN_OK);
}


/*****************************************************************************/
int
Zoltan_PHG_Timers_init(ZZ* zz)
{
  int i;
  int lenght;
  struct phg_timer_indices *timer;

  Zoltan_PHG_LB_Data* ptr = (Zoltan_PHG_LB_Data*)zz->LB.Data_Structure;

  if (!ptr) {
    ptr = Zoltan_PHG_LB_Data_alloc();
    zz->LB.Data_Structure = (void*)ptr;
  }
  if (!ptr)
    return (ZOLTAN_MEMERR);

  if (ptr->timers == NULL)
    ptr->timers = (struct phg_timer_indices*) ZOLTAN_MALLOC(sizeof(struct phg_timer_indices));

  timer= ptr->timers;
  if (!timer)
    return (ZOLTAN_MEMERR);

  lenght = sizeof(struct phg_timer_indices)/sizeof(int);

  for (i=0 ; i < lenght ; ++i)
    ((int*)timer)[i] = -1;

  return (ZOLTAN_OK);
/*   timer->all = */
/*   timer->build = */
/*   timer->setupvmap = */
/*   timer->parkway = */
/*   timer->patoh = */
/*   timer->retlist = */
/*   timer->finaloutput = */
/*   timer->match = */
/*   timer->coarse = */
/*   timer->refine = */
/*   timer->coarsepart = */
/*   timer->project = */
/*   timer->procred = */
/*   timer->vcycle = */
/*   timer->comerge = */
/*   timer->coshuffle = */
/*   timer->coremove = */
/*   timer->cotheend = */
/*   timer->rdrdivide = */
/*   timer->rdbefore = */
/*   timer->rdafter = */
/*   timer->rdsplit = */
/*   timer->rdredist = */
/*   timer->rdsend = */
/*   timer->rdwait = */
/*   timer->rfrefine = */
/*   timer->rfpins = */
/*   timer->rfiso = */
/*   timer->rfgain = */
/*   timer->rfheap = */
/*   timer->rfpass = */
/*   timer->rfroll = */
/*   timer->rfnonroot = */
/*   timer->cpart = */
/*   timer->cpgather = */
/*   timer->cprefine = -1; */
/*   for (i = 0; i < 7; i++) timer->matchstage[i] = -1; */
}

/******************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
