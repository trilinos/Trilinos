/* #include <malloc.h> */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "octant_const.h"
#include "util_const.h"
#include "hilbert_const.h"

typedef double Coord[3];

void set_method(double method_number) {
  
  GRAY = HILBERT = 0;
  if(method_number == 1)
    GRAY = 1;
  else if(method_number == 2)
    HILBERT = 1;
}

/*
 * void *my_malloc(int size)
 *
 * mallocs space in memory. checks for bad input and results
 */
void *my_malloc(int size) {
  char *ret;

  if (!size) {
    fprintf(stderr,"my_malloc: tried to malloc zero sized block\n");
    abort();
  }

  ret=malloc((unsigned)size);
  if (!ret) {
    fprintf(stderr,"my_malloc: failed trying to malloc %d bytes\n",size);
    perror("my_malloc:");
    abort();
  }
  return((void *)ret);
}


/*
 * in_box
 *
 * Returns true if pt lies withing the rectangular region
 * bounded by coordinate minimums "lower" and coordinate maximums
 * "upper"
 *
 */
int in_box(Coord pt, Coord lower, Coord upper) {
  int i;
  int j;

  if(dimension == 2)
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
 * bounds_to_origin_size()
 * 
 * convert octant bounds to an origin and size
 * NOTE: do not use origin and size for membership
 * in the octant, due to roundoff
 *
 */
void bounds_to_origin_size(Coord min, Coord max, Coord origin, Coord size) {
  int i;

  for (i=0; i<3; i++) {
    size[i]=max[i]-min[i];
    origin[i]=min[i]+size[i]/2.0;
  }
}


/*
 * bounds_to_origin()
 * 
 * convert octant bounds to an origin
 * NOTE: do not use origin and size for membership
 * in the octant, due to roundoff
 *
 */
void bounds_to_origin(Coord min, Coord max, Coord origin)
{
  int i;
  Coord size;

  for (i=0; i<3; i++) {
    size[i]=max[i]-min[i];
    origin[i]=min[i]+size[i]/2.0;
  }
}

void child_bounds_wrapper(pOctant oct, int input, Coord cmin, Coord cmax) {
  int cnum;                               /* child number */
  Coord min,
        max;
  Coord origin;

  /* get the bounds of an octant */
  POC_bounds(oct,min,max);
  /* calculate the origin from the bounds */
  bounds_to_origin(min,max,origin);

  if(HILBERT) {
    /* fprintf(stderr,"Using Hilbert\n"); */
    if(dimension == 3)
      cnum = hilbert_bounds(min, max, input);
    else {
      if(input < 4)
	cnum = hilbert2d_bounds(min,max,input);
      else
	cnum = input; 
    }
  }
  else if(GRAY) {
    /* fprintf(stderr,"Using Gray\n"); */
    cnum = convert_from_gray(input);        /* used for renumbering children */
  }
  else {
    /* fprintf(stderr,"Using Default\n"); */
    cnum = input;
  }

  child_bounds(min, max, origin, cnum, cmin, cmax);
}

extern int compare( unsigned * x, unsigned * y)
{ return ( *x == *y ) ? 0 : ( ( *x < *y ) ? -1 : 1 ); }

int hilbert2d_bounds(Coord min, Coord max, int cnum) {
  unsigned ihsfc[4][2];
  unsigned Nword = 1;
  Coord child_min[4],
        child_max[4];
  Coord centroid,
        origin;
  double center[2];
  int i, k;
  
  bounds_to_origin(min,max,origin);
  for(i=0; i<4; i++) {
    child_bounds(min, max, origin, i, child_min[i], child_max[i]);
    bounds_to_origin(child_min[i], child_max[i], centroid);
    center[0] = (centroid[0] - gmin[0])/(gmax[0] - gmin[0]);
    center[1] = (centroid[1] - gmin[1])/(gmax[1] - gmin[1]);
    fhsfc2d(center, &Nword, ihsfc[i]);
    ihsfc[i][1] = i;
  }

  qsort(ihsfc, 4, 2*sizeof(unsigned), (int (*)())compare );

  k=ihsfc[cnum][1];
  return(k);
}

int hilbert_bounds(Coord min, Coord max, int cnum) {
  unsigned ihsfc[8][2];
  unsigned Nword = 1;
  Coord child_min[8],
        child_max[8];
  Coord centroid,
        origin,
        center;
  int i, k;
  
  bounds_to_origin(min,max,origin);
  for(i=0; i<8; i++) {
    child_bounds(min, max, origin, i, child_min[i], child_max[i]);
    bounds_to_origin(child_min[i], child_max[i], centroid);
    center[0] = (centroid[0] - gmin[0])/(gmax[0] - gmin[0]);
    center[1] = (centroid[1] - gmin[1])/(gmax[1] - gmin[1]);
    center[2] = (centroid[2] - gmin[2])/(gmax[2] - gmin[2]);
    fhsfc3d(center, &Nword, ihsfc[i]);
    ihsfc[i][1] = i;
  }

  qsort(ihsfc, 8, 2*sizeof(unsigned), (int (*)())compare );

  k=ihsfc[cnum][1];
  return(k);
}

/*
 * child_bounds(pmin,pmax,porigin,childnum,cmin,cmax)
 *
 * given bounds and origin of a parent, compute a child's bounds
 * NOTE: relies on child octant numbering. Assumes in z-curve ordering, so
 *       needs conversion if using Gray code
 */
void child_bounds(Coord pmin, Coord pmax, Coord porigin, int cnum,
		  Coord cmin, Coord cmax) {
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

int child_which_wrapper(pOctant oct, Coord point) {
  int cnum;                               /* child number */
  int result;
  Coord min,
        max;
  Coord origin;

  POC_bounds(oct,min,max);                   /* get the bounds of the octant */
  bounds_to_origin(min,max,origin);            /* convert to bound to origin */
  
  cnum = child_which(origin, point);          /* find closest child to coord */
  
  if(HILBERT) {
    /* fprintf(stderr,"Using Hilbert\n"); */
    if(dimension == 3) 
      result = change_to_hilbert(min, max, origin, cnum);
    else 
      result = change_to_hilbert2d(min, max, origin, cnum);

    /* fprintf(stderr,"%d -> %d\n", cnum, result); */
  }
  else if(GRAY) {
    /* fprintf(stderr,"Using Gray\n"); */
    result = convert_to_gray(cnum);
  }
  else {
    /* fprintf(stderr,"Using Default\n"); */
    result = cnum;
  }

  return(result);
}

/*
 * int child_which(origin,point)
 *
 * Tell which child of octant having bounds will contain
 * the point.  (If point lies outside the octant, will return
 * the child number of the closest child.)
 * Note: Finds number in z-curve ordering. Needs to be converted if using
 *       Gray code or Hilbert code.
 */
int child_which(Coord porigin, Coord point) {
  int cnum;                                                 /* child number */
  int i;                                                    /* index counter */
  int j;
  
  if(dimension == 2)
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

double dist_point_box(Coord point, Coord min, Coord max) {
  Coord dist;
  double pdist;
  int zero;
  int i;

  if (in_box(point,min,max))
    return(0);

  zero=0;
  pdist=0;
  for (i=0; i<3; i++) {
    if (point[i]<min[i]) {
      pdist=min[i]-point[i];
      dist[i]=pdist;
    }
    else
      if (point[i]>max[i]) {
	pdist=point[i]-max[i];
	dist[i]=pdist;
      }
      else {
	dist[i]=0;
	zero++;
      }
  }

  if (zero<2)
    pdist=sqrt((dist[0]*dist[0]) + (dist[1]*dist[1]) + (dist[2]*dist[2]));

  return(pdist);
}

/*
 * int convert_from_gray(int input)
 *
 * takes a Gray code number and converts it into Binary
 */
int convert_from_gray(int input) {
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
 * int convert_to_gray(int input)
 *
 * takes in a (binary) number and returns the Gray code of that number
 */
int convert_to_gray(int input) {
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

/* lookup table for traversal informaion based on orientation */
/* unsigned int ltable[4] = {0x099, 0x05A, 0x681, 0x642}; */
unsigned int ltable[4] = {0x0B4BD9, 0x06CDDA, 0x6B7B01, 0x66FD02};

/* 000 010 011 001 */

/* 000 010 110 010 101 111 011 001 */
/* 000 111 000 000 000 000 000 000 */

/*
 * int convert_from_hilbert(int number, int orientation)
 *
 * Takes in a number and the traversal orientation and returns the nth
 * number in the orientation's ordering 
 */
int convert_from_hilbert(int n, int o) {
  unsigned int tmp;           /* temporary variable for intermediate results */
  int i,                      /* index counter */
      k;                      /* shifting index */
  int result;                 /* value to be return to the caller */
  unsigned int mask = 0xE00;  /* mask to look at only specific entries */
  unsigned int zeromask = 0xE00000;

  tmp = ltable[o];
  
  for(i=0; i<n; i++)
    zeromask = (zeromask >> 3);
  
  k = 21 - (3*n);                    /* number of times to right shift tmp */
  tmp = tmp & zeromask;
  result = tmp >> k;                                     /* get the result */
  
  /* fprintf(stderr,"n = %d, i = %d, mask = %x, tmp = %x, result = %d\n", 
   *   n, i, zeromask, tmp, result);
   */
  
  return (result);
#if 0
  /* only works in 2D, cannot handle higher octants */
  if (n >= 4)
    return n;

  tmp = ltable[o];                     /* get the ordering from lookup table */
  for(i=0; i<n; i++)                  /* shift mask to get appropriate entry */
    mask = mask >> 3;
  
  k = 9 - (3*n);                       /* number of times to right shift tmp */
  tmp = tmp & mask;
  result = tmp >> k;                                       /* get the result */
  
  return (result);
#endif
}

/*
 * int convert_to_hilbert(int number, int orientation)
 *
 * converts from the default ordering to the hilbert ordering. Only works in
 * 2D.
 */
int convert_to_hilbert(int n, int o) {
  unsigned int mask = 0xE00;        /* mask to look at only specific entries */
  int i,                            /* index counter */
      k;                            /* shifting counter */
  unsigned int test;                /* holds intermediate results */
  unsigned int zeromask = 0xE00000;            /* 111000000000000000000000 */

  /* 000 010 110 010 101 111 011 001 */
  /* 111 000 000 000 000 000 000 000 */

  for(i=0; ; i++) {
    test = ltable[o] & zeromask;
    k = 21 - (3*i);
    test = test >> k;
    zeromask = zeromask >> 3;
    if (test == n)
      break;
    if(i > 8) {
      fprintf(stderr, "ERROR in conversion\n");
      abort();
    }
  }
  
  /* fprintf(stderr,"n = %d -> i= %d\n", n, i); */
  return (i);

#if 0
  /* only works in 2D, cannot handle higher octants */
  if (n >= 4)
    return n;

  for(i=0; ; i++) {
    test = ltable[o] & mask;
    k = 9 - (3*i);
    test = test >> k;
    mask = mask >> 3;
    if (test == n)
      break;
    if(i > 4) {
      fprintf(stderr, "ERROR in conversion\n");
      abort();
    }
  }
  return (i);
#endif
}

/*
 * int child_orientation(int parent_orientation, int child_number)
 *
 * returns the orientation the the child requested 
 */
int child_orientation(int o, int cnum) {
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

int change_to_hilbert2d(Coord min, Coord max, Coord origin, int cnum) {
  unsigned ihsfc[4][2];
  unsigned Nword = 1;
  Coord child_min[4],
        child_max[4];
  Coord centroid;
  double center[2];
  int i;
  
  for(i=0; i<4; i++) {
    child_bounds(min, max, origin, i, child_min[i], child_max[i]);
    bounds_to_origin(child_min[i], child_max[i], centroid);
    center[0] = (centroid[0] - gmin[0])/(gmax[0] - gmin[0]);
    center[1] = (centroid[1] - gmin[1])/(gmax[1] - gmin[1]);
    fhsfc2d(center, &Nword, ihsfc[i]);
    ihsfc[i][1] = i;
  }

  qsort(ihsfc, 4, 2*sizeof(unsigned), (int (*)())compare);

  for(i=0; i<4; i++) {
    /* 
     * fprintf(stderr,"ihsfc[%d][1] = %d, cnum = %d\n", i, ihsfc[i][1], cnum);
     */
    if(ihsfc[i][1] == cnum)
      break;
  }
  
  if(i>=4) {
    fprintf(stderr,"ERROR in child_which_wrapper, invalid cnum %d\n", i);
    abort();
  }
  
  return(i);
}

int change_to_hilbert(Coord min, Coord max, Coord origin, int cnum) {
  unsigned ihsfc[8][2];
  unsigned Nword = 1;
  Coord child_min[8],
        child_max[8];
  Coord centroid,
        center;
  int i;
  
  for(i=0; i<8; i++) {
    child_bounds(min, max, origin, i, child_min[i], child_max[i]);
    bounds_to_origin(child_min[i], child_max[i], centroid);
    center[0] = (centroid[0] - gmin[0])/(gmax[0] - gmin[0]);
    center[1] = (centroid[1] - gmin[1])/(gmax[1] - gmin[1]);
    center[2] = (centroid[2] - gmin[2])/(gmax[2] - gmin[2]);
    fhsfc3d(center, &Nword, ihsfc[i]);
    ihsfc[i][1] = i;
  }

  qsort(ihsfc, 8, 2*sizeof(unsigned), (int (*)())compare);

  for(i=0; i<8; i++) {
    /*
     * fprintf(stderr,"ihsfc[%d][1] = %d, cnum = %d\n", i, ihsfc[i][1], cnum);
     */
    if(ihsfc[i][1] == cnum)
      break;
  }
  
  if(i>=8) {
    fprintf(stderr,"ERROR in child_which_wrapper, invalid cnum %d\n", i);
    abort();
  }
  
  return(i);
}
