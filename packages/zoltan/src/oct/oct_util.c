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
#include "hilbert_const.h"
#include "dfs_const.h"
#include "all_allo_const.h"

/* function prototype */
static int compare( unsigned * x, unsigned * y);
static int hilbert_bounds(OCT_Global_Info *OCT_info,COORD min, COORD max, COORD cmin[], COORD cmax[]);
static int hilbert2d_bounds(OCT_Global_Info *OCT_info,COORD min, COORD max, COORD cmin[], COORD cmax[]);

void LB_set_method(OCT_Global_Info *OCT_info, int method_number) {
  
  OCT_info->GRAY = OCT_info->HILBERT = 0;
  if(method_number == 1)
    OCT_info->GRAY = 1;
  else if(method_number == 2)
    OCT_info->HILBERT = 1;
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
      /* fprintf(stderr,"(%d) %lf %lf, %lf %lf, %lf %lf -> %d\n", 
	         OCT_localpid, pt[0], pt[1], lower[0], 
                 lower[1], upper[0], upper[1], i); */
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
  int cnum;                               /* child number */
  int i;
  COORD min,
        max;
  COORD origin;

  /* get the bounds of an octant */
  POC_bounds(oct,min,max);
  /* calculate the origin from the bounds */
  LB_bounds_to_origin(min,max,origin);

  if(OCT_info->HILBERT) {
    /* fprintf(stderr,"Using Hilbert\n"); */ 
    if(OCT_info->OCT_dimension == 3)
      hilbert_bounds(OCT_info,min, max, cmin, cmax);
    else {
      hilbert2d_bounds(OCT_info,min, max, cmin, cmax);
      for(i=4; i<8; i++) {
	LB_child_bounds(min, max, origin, i, cmin[i], cmax[i]);
      }
    }
  }
  else if(OCT_info->GRAY) {
    /* fprintf(stderr,"Using Gray\n"); */
    for(i=0; i<8; i++) {
      cnum = LB_convert_from_gray(i);       /* used for renumbering children */
      LB_child_bounds(min, max, origin, cnum, cmin[i], cmax[i]);
    }
  }
  else {
    /* fprintf(stderr,"Using Default\n"); */
    for(i=0; i<8; i++) {
      LB_child_bounds(min, max, origin, i, cmin[i], cmax[i]);
    }
  }

  /* LB_child_bounds(min, max, origin, cnum, cmin, cmax); */
}

static int compare( unsigned * x, unsigned * y)
{ return ( *x == *y ) ? 0 : ( ( *x < *y ) ? -1 : 1 ); }

static int hilbert2d_bounds(OCT_Global_Info *OCT_info,
                            COORD min, COORD max, COORD cmin[], COORD cmax[])
{
  unsigned ihsfc[4][2];
  unsigned Nword = 1;
  COORD child_min[4],
        child_max[4];
  COORD centroid,
        origin;
  double center[2];
  int i, k, j;
  
  LB_bounds_to_origin(min,max,origin);
  for(i=0; i<4; i++) {
    LB_child_bounds(min, max, origin, i, child_min[i], child_max[i]);
    LB_bounds_to_origin(child_min[i], child_max[i], centroid);
    center[0] = (centroid[0] - OCT_info->OCT_gmin[0])/(OCT_info->OCT_gmax[0] - OCT_info->OCT_gmin[0]);
    center[1] = (centroid[1] - OCT_info->OCT_gmin[1])/(OCT_info->OCT_gmax[1] - OCT_info->OCT_gmin[1]);
    LB_fhsfc2d(center, &Nword, ihsfc[i]);
    ihsfc[i][1] = i;
  }

  qsort(ihsfc, 4, 2*sizeof(unsigned),
       (int (*)(const void *, const void *))compare );

  for(i=0; i<4; i++) {
    k=ihsfc[i][1];
    for(j=0; j<3; j++) {
      cmin[i][j] = child_min[k][j];
      cmax[i][j] = child_max[k][j];
    }
  }
  return(k);
}

static int hilbert_bounds(OCT_Global_Info *OCT_info,
                          COORD min, COORD max, COORD cmin[], COORD cmax[])
{
  unsigned ihsfc[8][2];
  unsigned Nword = 1;
  COORD child_min[8],
        child_max[8];
  COORD centroid,
        origin,
        center;
  int i, k, j;
  
  LB_bounds_to_origin(min,max,origin);
  for(i=0; i<8; i++) {
    LB_child_bounds(min, max, origin, i, child_min[i], child_max[i]);
    LB_bounds_to_origin(child_min[i], child_max[i], centroid);
    center[0] = (centroid[0] - OCT_info->OCT_gmin[0])/(OCT_info->OCT_gmax[0] - OCT_info->OCT_gmin[0]);
    center[1] = (centroid[1] - OCT_info->OCT_gmin[1])/(OCT_info->OCT_gmax[1] - OCT_info->OCT_gmin[1]);
    center[2] = (centroid[2] - OCT_info->OCT_gmin[2])/(OCT_info->OCT_gmax[2] - OCT_info->OCT_gmin[2]);
    LB_fhsfc3d(center, &Nword, ihsfc[i]);
    ihsfc[i][1] = i;
  }

  qsort(ihsfc, 8, 2*sizeof(unsigned), 
       (int (*)(const void *, const void *))compare );

  for(i=0; i<8; i++) {
    k=ihsfc[i][1];
    for(j=0; j<3; j++) {
      cmin[i][j] = child_min[k][j];
      cmax[i][j] = child_max[k][j];
    }
  }
  return(k);
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
  int cnum;                               /* child number */
  int result;
  COORD min,
        max;
  COORD origin;

  POC_bounds(oct,min,max);                   /* get the bounds of the octant */
  LB_bounds_to_origin(min,max,origin);         /* convert to bound to origin */
  
  cnum = LB_child_which(OCT_info,origin, point);       /* find closest child to coord */
  
  if(OCT_info->HILBERT) {
    /* fprintf(stderr,"Using Hilbert\n"); */
    if(OCT_info->OCT_dimension == 3) 
      result = LB_change_to_hilbert(OCT_info,min, max, origin, cnum);
    else 
      result = LB_change_to_hilbert2d(OCT_info,min, max, origin, cnum);

    /* fprintf(stderr,"%d -> %d\n", cnum, result); */
  }
  else if(OCT_info->GRAY) {
    /* fprintf(stderr,"Using Gray\n"); */
    result = LB_convert_to_gray(cnum);
  }
  else {
    /* fprintf(stderr,"Using Default\n"); */
    result = cnum;
  }

  return(result);
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


/*
 * int LB_convert_from_gray(int input)
 *
 * takes a Gray code number and converts it into Binary
 */
int LB_convert_from_gray(int input) {
  int prev,                      /* previous result */
      num,                       /* bit currently being examines */
      res,                       /* result to be returned */
      mask;                      /* used to mask out bits not being coverted */
  int i;                         /* index counter */

  /* initalize variables */
  mask = 0x04;
  res = prev = 0;

  for(i=2; i>=0; i--) {
    res = res << 1;
    num = input & mask;                            /* get one bit from input */
    num = num >> i;
    res |= num ^ prev;                             /* exclusive or with prev */
    mask = mask >> 1;                                          /* shift mask */
    prev = num;
  }

  /* return the result */
  return(res);
}

/*
 * int LB_convert_to_gray(int input)
 *
 * takes in a (binary) number and returns the Gray code of that number
 */
int LB_convert_to_gray(int input) {
  int prev,                      /* previous result */
      num,                       /* bit currently being examines */
      res,                       /* result to be returned */
      mask;                      /* used to mask out bits not being coverted */
  int i;                         /* index counter */

  /* initalize variables */
  mask = 0x04;
  res = prev = 0;

  for(i=2; i>=0; i--) {
    res = res << 1;
    num = input & mask;                            /* get one bit from input */
    num = num >> i;
    res |= num ^ prev;                             /* exclusive or with prev */
    mask = mask >> 1;                                          /* shift mask */
    prev = num ^ prev;
  }

  /* return the result */
  return(res);
}


/*
 * int LB_child_orientation(int parent_orientation, int child_number)
 *
 * returns the orientation the the child requested 
 */
int LB_child_orientation(int o, int cnum) {
  if((cnum == 1) || (cnum == 6) || (cnum < 0))
    return o;
  else if(cnum == 7)
    return ((o+2) % 4);
  else
    return((o ^ 0x01));

#if 0
  if((cnum == 1) || (cnum == 2) || (cnum < 0))
    return o;
  else if(cnum == 3)
    return ((o+2) % 4);
  else
    return((o ^ 0x01));
#endif
}

int LB_change_to_hilbert2d(OCT_Global_Info *OCT_info,
                           COORD min, COORD max, COORD origin, int cnum) {
  unsigned ihsfc[4][2];
  unsigned Nword = 1;
  COORD child_min[4],
        child_max[4];
  COORD centroid;
  double center[2];
  int i;
  
  for(i=0; i<4; i++) {
    LB_child_bounds(min, max, origin, i, child_min[i], child_max[i]);
    LB_bounds_to_origin(child_min[i], child_max[i], centroid);
    center[0] = (centroid[0] - OCT_info->OCT_gmin[0])/(OCT_info->OCT_gmax[0] - OCT_info->OCT_gmin[0]);
    center[1] = (centroid[1] - OCT_info->OCT_gmin[1])/(OCT_info->OCT_gmax[1] - OCT_info->OCT_gmin[1]);
    LB_fhsfc2d(center, &Nword, ihsfc[i]);
    ihsfc[i][1] = i;
  }

  qsort(ihsfc, 4, 2*sizeof(unsigned),
       (int (*)(const void *, const void *))compare );

  for(i=0; i<4; i++) {
    /* 
     * fprintf(stderr,"ihsfc[%d][1] = %d, cnum = %d\n", i, ihsfc[i][1], cnum);
     */
    if(ihsfc[i][1] == cnum)
      break;
  }
  
  if(i>=4) {
    fprintf(stderr,"ERROR in LB_child_which_wrapper, invalid cnum %d\n", i);
    abort();
  }
  
  return(i);
}

int LB_change_to_hilbert(OCT_Global_Info *OCT_info,
                         COORD min, COORD max, COORD origin, int cnum) {
  unsigned ihsfc[8][2];
  unsigned Nword = 1;
  COORD child_min[8],
        child_max[8];
  COORD centroid,
        center;
  int i;
  
  for(i=0; i<8; i++) {
    LB_child_bounds(min, max, origin, i, child_min[i], child_max[i]);
    LB_bounds_to_origin(child_min[i], child_max[i], centroid);
    center[0] = (centroid[0] - OCT_info->OCT_gmin[0])/(OCT_info->OCT_gmax[0] - OCT_info->OCT_gmin[0]);
    center[1] = (centroid[1] - OCT_info->OCT_gmin[1])/(OCT_info->OCT_gmax[1] - OCT_info->OCT_gmin[1]);
    center[2] = (centroid[2] - OCT_info->OCT_gmin[2])/(OCT_info->OCT_gmax[2] - OCT_info->OCT_gmin[2]);
    LB_fhsfc3d(center, &Nword, ihsfc[i]);
    ihsfc[i][1] = i;
  }

  qsort(ihsfc, 8, 2*sizeof(unsigned),
       (int (*)(const void *, const void *))compare );

  for(i=0; i<8; i++) {
    /*
     * fprintf(stderr,"ihsfc[%d][1] = %d, cnum = %d\n", i, ihsfc[i][1], cnum);
     */
    if(ihsfc[i][1] == cnum)
      break;
  }
  
  if(i>=8) {
    fprintf(stderr,"ERROR in LB_child_which_wrapper, invalid cnum %d\n", i);
    abort();
  }
  
  return(i);
}

/*****************************************************************************/
void LB_OCT_Free_Structure(LB *lb)
{
/*
 * Deallocate the persistent RCB data structures in lb->Structure.
 */
OCT_Global_Info *OCT_info = (OCT_Global_Info *) (lb->Data_Structure);

  if (OCT_info != NULL) {
    LB_FREE(&(lb->Data_Structure));
  }
}
