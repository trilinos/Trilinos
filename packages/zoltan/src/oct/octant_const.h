/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __OCTANT_CONST_H
#define __OCTANT_CONST_H

#include "lb_const.h"

typedef double COORD[3];
typedef enum {LOCALOCT, REMOTEOCT} OctType;
typedef struct Region_Node* pRegion;   /* typedef for a pointer to a region  */
typedef struct Region_Node {           /* region = area in 3-space           */
  struct Region_Node *next;            /* pointer to next region in list     */
  COORD  Coord;                        /* centroid location of region        */
  double Weight;                       /* weight of Region - default is 1    */
  LB_ID_PTR Global_ID;                 /* Global ID for the object.          */
  LB_ID_PTR Local_ID;                  /* Local ID for the object.           */
  int Proc;                            /* Processor originally owning the obj*/
  int    attached;                     /* flag to see if region was attached */
} Region;

typedef struct _Octant* pOctant;     /* typedef for a pointer to an octant   */
typedef struct _Octant {             /* octant tree node that has 8 children */
  OctType type;
  pRegion list;                      /* list of regions associated to octant */
  struct _Octant *child[8];          /* array of octant's children           */
  struct _Octant *parent;            /* parent of the octant                 */
  struct _Octant *remoteptr;         /* parent of the octant                 */
  COORD  min;                        /* minimum bounds of an octant          */
  COORD  max;                        /* max bounds of an octant              */
  int ppid;                          /* parent pid, -1 mean a local root     */
  int id;                            /* octant's id number                   */
  int which;                         /* which child of parent                */
  int numChild;                      /* number of children, 0 == terminal    */
  float cost;                        /* cost of the octant                   */
  int npid;                          /* where to migrate to                  */
  int cpid[8];                       /* the pid of the children              */
  int dir;
  int mapidx;
  double area;
  /* int orientation;*/              /* octant traversal orientation         */
} Octant;

/* #ifdef LGG_MIGOCT */
extern void    LB_Oct_setID(pOctant oct, int id);
/* #endif  *//* LGG_MIGOCT */
extern void    LB_Oct_initCounters();
extern pOctant LB_Oct_newremote();
extern pOctant LB_Oct_new();
extern void    LB_Oct_free(pOctant oct);
extern int     LB_Oct_id(pOctant oct);
extern int     LB_Oct_dir(pOctant oct);
extern int     LB_Oct_mapidx(pOctant oct);
extern void    LB_Oct_setDir(pOctant oct, int dir);
extern void    LB_Oct_setchildnum(pOctant oct, int childnum);
extern void    LB_Oct_setchildren(pOctant oct, pOctant *children, int *cpids);
extern void    LB_Oct_setchild(pOctant oct, int i, pOctant child);
extern void    LB_Oct_setbounds(pOctant oct, COORD min, COORD max);
extern void    LB_Oct_setMapIdx(pOctant oct, int idx);
extern void    LB_Oct_bounds(pOctant oct, COORD min, COORD max);
extern pOctant LB_Oct_parent(pOctant oct);
extern pOctant LB_Oct_child(pOctant oct, int i);
extern void    LB_Oct_cpids(pOctant oct, int cpids[8]);
extern int     LB_Oct_children(pOctant oct, pOctant children[8]);
extern int     LB_Oct_isTerminal(pOctant oct);
extern pRegion LB_Oct_regionlist(pOctant oct);
extern int     LB_Oct_addRegion(LB *lb, pOctant oct, pRegion region);
extern void    LB_Oct_clearRegions(pOctant oct);
extern int     LB_Oct_nRegions(pOctant oct);
extern void    LB_Oct_modify_cost(pOctant oct, float cost);
extern void    LB_Oct_modify_newpid(pOctant oct, int newpid);
extern int     LB_Oct_data_newpid(pOctant oct);
extern int     LB_Oct_nOctants(void);
extern void    LB_Oct_origin_volume(pOctant oct,COORD origin,double *volume);
extern int     LB_Oct_Ppid(pOctant octant);
extern int     LB_Oct_Cpid(pOctant octant, int i);
extern int     LB_Oct_childnum(pOctant octant);
extern void    LB_Oct_setCpid(pOctant octant, int i, int cpid); 
/* extern pOctant PLB_Oct_nextDfs(pOcant octant);
 * extern void    PLB_Oct_setOrientation(pOctant octant, int orientation);
 * extern int     PLB_Oct_getOrientation(pOctant octant);
 */

#endif
