/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "SFC.h"
#include "lb_const.h"
#include "octant_const.h"
#include "oct_util_const.h"
#include "dfs_const.h"
#include <values.h>

int LB_get_child_dir(OCT_Global_Info *OCT_info, int dir, int cnum) {
  if(OCT_info->HILBERT) {
    if(OCT_info->OCT_dimension == 3) 
      return Hilbert3dRuleMap[dir][cnum];
    else
      return Hilbert2dRuleMap[dir][cnum];
  }
  else
    if(OCT_info->GRAY) {
      if(OCT_info->OCT_dimension == 3) 
	return Gray3dRuleMap[dir][cnum];
      else
	return Gray2dRuleMap[dir][cnum];
    }
  return 0;
}

static int LB_convert_idx_to_map(OCT_Global_Info *OCT_info, int dir, int cnum) {
  int result = cnum;
  if(OCT_info->HILBERT) {
    if(OCT_info->OCT_dimension == 3) {
      result = Hilbert3dIndexToMap[dir][cnum];
    }
    else
      result = Hilbert2dIndexToMap[dir][cnum];
  }
  else 
    if(OCT_info->GRAY) {
      if(OCT_info->OCT_dimension == 3) 
	result = Gray3dIndexToMap[dir][cnum];
      else
	result = Gray2dIndexToMap[dir][cnum];
    }
  return result;
}

int LB_convert_idx_from_map(OCT_Global_Info *OCT_info, int dir, int cnum) {
  int result = cnum;
  if(OCT_info->HILBERT) {
    if(OCT_info->OCT_dimension == 3) {
      result = Hilbert3dIndexFromMap[dir][cnum];
    }
    else
      result = Hilbert2dIndexFromMap[dir][cnum];
  }
  else 
    if(OCT_info->GRAY) {
      if(OCT_info->OCT_dimension == 3) 
	result = Gray3dIndexFromMap[dir][cnum];
      else
	result = Gray2dIndexFromMap[dir][cnum];
    }
  return result;
}

void LB_set_method(OCT_Global_Info *OCT_info, int method_number) {
  OCT_info->GRAY = OCT_info->HILBERT = 0;
  if(method_number == 1) {
    OCT_info->GRAY = 1;
  }
  else if(method_number == 2) {
    OCT_info->HILBERT = 1;
  }
}

/*
 * LB_in_box
 *
 * Returns true if pt lies withing the rectangular region
 * bounded by coordinate minimums "lower" and coordinate maximums
 * "upper"
 *
 */
int LB_in_box(OCT_Global_Info *OCT_info, COORD pt, COORD lower, COORD upper) {
  int i;
  int j;

  if(OCT_info->OCT_dimension == 2)
    j = 2;                                            /* ignore z coord info */
  else
    j = 3;

  for (i=0; i<j; i++)
    if (pt[i]<=lower[i] || pt[i]>upper[i]) {
      return(0);
    }
  
  return(1);
}


/*
 * LB_in_box_closure
 *
 * Returns true if pt lies withing the rectangular region
 * bounded by coordinate minimums "lower" and coordinate maximums
 * "upper"
 *
 */
int LB_in_box_closure(OCT_Global_Info *OCT_info, COORD pt, COORD lower, COORD upper) {
  int i;
  int j;

  if(OCT_info->OCT_dimension == 2)
    j = 2;                                            /* ignore z coord info */
  else
    j = 3;

  for (i=0; i<j; i++)
    if (pt[i]<lower[i] || pt[i]>upper[i]) {
      return(0);
    }
  
  return(1);
}


/*
 * LB_bounds_to_origin_size()
 * 
 * convert octant bounds to an origin and size
 * NOTE: do not use origin and size for membership
 * in the octant, due to roundoff
 *
 */
void LB_bounds_to_origin_size(COORD min, COORD max, COORD origin, COORD size) {
  int i;

  for (i=0; i<3; i++) {
    size[i]=max[i]-min[i];
    origin[i]=min[i]+size[i]/2.0;
  }
}


/*
 * LB_bounds_to_origin()
 * 
 * convert octant bounds to an origin
 * NOTE: do not use origin and size for membership
 * in the octant, due to roundoff
 *
 */
void LB_bounds_to_origin(COORD min, COORD max, COORD origin)
{
  int i;
  COORD size;

  for (i=0; i<3; i++) {
    size[i]=max[i]-min[i];
    origin[i]=min[i]+size[i]/2.0;
  }
}

void LB_child_bounds_wrapper(OCT_Global_Info *OCT_info,
                             pOctant oct, COORD cmin[], COORD cmax[]) {
  int i, j;
  COORD min,
        max;
  COORD origin;
  /* get the bounds of an octant */
  LB_Oct_bounds(oct,min,max);
  /* calculate the origin from the bounds */
  LB_bounds_to_origin(min,max,origin);

  /* KDDKDD 3/01  Added special cases depending on OCT_info->OCT_dimension.
   * KDDKDD 3/01  When OCT_dimension == 2, LB_convert_idx_to_map was called
   * KDDKDD 3/01  with values >= 4, which are not supported in the GRAY and
   * KDDKDD 3/01  HILBERT maps used in LB_convert_idx_to_map.
   * KDDKDD 3/01  This change calls LB_convert_idx_to_map with only valid 
   * KDDKDD 3/01  values, and sets cmin and cmax so that the loop following 
   * KDDKDD 3/01  the call to LB_child_bounds_wrapper "continues" for children
   * KDDKDD 3/01  4-7, rather than initializing them.
   */

  if (OCT_info->OCT_dimension == 3) {
    for(i=0; i<8; i++) {
      LB_child_bounds(min, max, origin, LB_convert_idx_to_map(OCT_info, oct->dir, i), cmin[i], cmax[i]);  
    }  
  }
  else if (OCT_info->OCT_dimension == 2) {
    for(i=0; i<4; i++) {
      LB_child_bounds(min, max, origin, LB_convert_idx_to_map(OCT_info, oct->dir, i), cmin[i], cmax[i]);  
    }  
    for(i=4; i<8; i++) {
      for(j=0;j<3;j++) {
        cmin[i][j] = MAXDOUBLE;
        cmax[i][j] = -MAXDOUBLE;
      }
    }
  }
  else {
    fprintf(stderr, "LB_child_bounds_wrapper:  Invalid OCT_dimension\n");
    abort();
  }
}

/*
 * LB_child_bounds(pmin,pmax,porigin,childnum,cmin,cmax)
 *
 * given bounds and origin of a parent, compute a child's bounds
 * NOTE: relies on child octant numbering. Assumes in z-curve ordering, so
 *       needs conversion if using Gray code
 */
void LB_child_bounds(COORD pmin, COORD pmax, COORD porigin, int cnum,
		     COORD cmin, COORD cmax) {
  int i;                                  /* index counter */
  int place;                              /* place currently being looked at */

  place=1;
  for (i=0; i<3; i++) {
    if (cnum & place) {
      cmin[i]=porigin[i];
      cmax[i]=pmax[i];
    }
    else {
      cmin[i]=pmin[i];
      cmax[i]=porigin[i];
    }
    place*=2;
  }
}

int LB_child_which_wrapper(OCT_Global_Info *OCT_info,pOctant oct, COORD point) {
  COORD min,max;
  COORD origin;
  LB_Oct_bounds(oct,min,max);                   /* get the bounds of the octant */
  LB_bounds_to_origin(min,max,origin);         /* convert to bound to origin */
  return LB_convert_idx_from_map(OCT_info, oct->dir, LB_child_which(OCT_info,origin, point));
}

/*
 * int LB_child_which(origin,point)
 *
 * Tell which child of octant having bounds will contain
 * the point.  (If point lies outside the octant, will return
 * the child number of the closest child.)
 * Note: Finds number in z-curve ordering. Needs to be converted if using
 *       Gray code or Hilbert code.
 */
int LB_child_which(OCT_Global_Info *OCT_info,COORD porigin, COORD point) {
  int cnum;                                                 /* child number */
  int i;                                                    /* index counter */
  int j;
  
  if(OCT_info->OCT_dimension == 2)
    j = 1;                                       /* ignore z coordinate info */
  else
    j = 2;

  cnum=0;
  for (i=j; i>=0; i--) {
    cnum= cnum<<1;
    
    if (point[i]>porigin[i])
      cnum|=1;
  }

  return(cnum);
}

/*****************************************************************************/
void LB_OCT_Free_Structure(LB *lb)
{
/*
 * Deallocate the persistent OCT data structures in lb->Structure.
 */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *) (lb->Data_Structure);
  pRList  RootList;                           /* list of all local roots */
  pOctant RootOct;
  int i;
  if (OCT_info != NULL) {
    /* free octree */
    /* KDDKDD In distribution from RPI, the following frees were commented out.
     * KDDKDD I do not think they should be commented out, as Debug_Memory
     * KDDKDD reports lots of memory leaks.
     * KDDKDD I have un-commented them.  12/2000
     */
    /* KDDKDDFREE Rearranged the frees because the global octants were trying to
     * KDDKDDFREE access children that were already deleted in LB_POct_delTree.
     * KDDKDDFREE Now delete global octants first, then the tree.  1/2001
     */
    /* free global octants */
    for(i = 0; i < OCT_info->mapsize; i++) {
      RootList = OCT_info->map[i].list;
      while((RootOct = RL_nextRootOctant(&RootList))) {
        LB_Oct_free(OCT_info, &RootOct);
        /* KDDKDDFREE Set RootList's oct ptr to NULL as the octant is deleted.*/
        RootList->oct = NULL;
      }
      RL_freeList(&(OCT_info->map[i].list));
    }
    /* free map array */
    ZOLTAN_FREE(&(OCT_info->map));
    /* Free the octree */
    RootList = LB_POct_localroots(OCT_info); 
    while ((RootOct = RL_nextRootOctant(&RootList))){
      /* KDDKDDFREE  The LB_POct_delTree also frees RootList, since
       * KDDKDDFREE  LB_POct_free frees the root list entry for root octants. 
       * KDDKDDFREE  Need to reset RootList in each iteration. */
      LB_POct_delTree(OCT_info, &RootOct);
      RootList = LB_POct_localroots(OCT_info);
    }
    RL_freeList(&(OCT_info->OCT_rootlist));
    /* KDDKDD End of previously commented-out section. 12/2000 */

    ZOLTAN_FREE(&(lb->Data_Structure));
  }
}
