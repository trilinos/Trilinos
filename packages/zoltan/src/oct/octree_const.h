/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __OCT_OCTREE_CONST_H
#define __OCT_OCTREE_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "octant_const.h"
#include "rootlist_const.h"


typedef struct {
  int npid;
  COORD min;
  COORD max;
  pRList list;
} Map;

typedef struct {
  pRList OCT_rootlist;          /* list of all the local roots          */
  Map    *map;
  int mapsize;
  int OCT_localpid;             /* the processor id                     */
  COORD OCT_gmin;
  COORD OCT_gmax;
  int OCT_dimension;
  int GRAY;
  int HILBERT;
} OCT_Global_Info;

/* extern Map *array; */

#define vector_set(r,a)     \
        ( (r)[0]=(a)[0], (r)[1]=(a)[1], (r)[2]=(a)[2] )
#define vector_eq(a,b)      \
        ( (a)[0]==(b)[0] && (a)[1]==(b)[1] && (a)[2]==(b)[2] )
#define vector_lt(a,b)      \
        ( (a)[0]<(b)[0] && (a)[1]<(b)[1] && (a)[2]<(b)[2] )
#define vector_gt(a,b)      \
        ( (a)[0]>(b)[0] && (a)[1]>(b)[1] && (a)[2]>(b)[2] )
#define vector_add(r,a,b)     \
        ( (r)[0]=(a)[0]+(b)[0],     \
          (r)[1]=(a)[1]+(b)[1],     \
          (r)[2]=(a)[2]+(b)[2] )
#define vector_set_comp(r,x,y,z) \
        ( (r)[0]=(x), (r)[1]=(y), (r)[2]=(z) )
#define vector_cmult(r,c,a)       \
        ( (r)[0]=c*(a)[0],        \
          (r)[1]=c*(a)[1],        \
          (r)[2]=c*(a)[2] )
#define vector_divc(r,a,c)        \
        ( (r)[0]=(a)[0]/c,        \
          (r)[1]=(a)[1]/c,        \
          (r)[2]=(a)[2]/c )
#define vector_dist(a,b) \
        ( sqrt( ((a)[0]-(b)[0])*((a)[0]-(b)[0]) + \
                ((a)[1]-(b)[1])*((a)[1]-(b)[1]) + \
                ((a)[2]-(b)[2])*((a)[2]-(b)[2]) ) )

/* Octree functions */

extern OCT_Global_Info *Zoltan_Oct_POct_init(ZZ *, int pid, int dim);
extern pOctant Zoltan_Oct_POct_new(OCT_Global_Info *);
extern void    Zoltan_Oct_POct_free(OCT_Global_Info *OCT_info, pOctant *oct);
extern void    Zoltan_Oct_POct_setparent(OCT_Global_Info *OCT_info,pOctant oct, pOctant parent, int ppid);
extern pRList  Zoltan_Oct_POct_localroots(OCT_Global_Info *);
extern pOctant Zoltan_Oct_POct_nextDfs(OCT_Global_Info *OCT_info, pOctant octant);
extern int     Zoltan_Oct_POct_local(OCT_Global_Info *OCT_info, pOctant octant, int i);
extern int     Zoltan_Oct_POct_delTree(OCT_Global_Info *OCT_info,pOctant *root);

/* KDDKDDFREE moved to octree_const.h to allow OCT_Global_Info arg. */
extern void    Zoltan_Oct_free(OCT_Global_Info *OCT_info, pOctant *oct);
extern int     RL_delRootOctant(OCT_Global_Info *OCT_info, pRList *rootlist, pOctant oct);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /*__OCT_OCTREE_CONST_H*/
