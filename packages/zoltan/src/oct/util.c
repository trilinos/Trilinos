/* #include <malloc.h> */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util_const.h"

typedef double Coord[3];

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

  for (i=0; i<3; i++)
    if (pt[i]<=lower[i] || pt[i]>upper[i])
      return(0);

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



/*
 * child_bounds(pmin,pmax,porigin,childnum,cmin,cmax)
 *
 * given bounds and origin of a parent, compute a child's bounds
 * 
 *
 * NOTE: relies on child octant numbering
 */

void child_bounds(Coord pmin, Coord pmax, Coord porigin, int cnum,
		  Coord cmin, Coord cmax) {
  int i;
  int place;

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




/*
 * int child_which(origin,point)
 *
 * Tell which child of octant having bounds will contain
 * the point.  (If point lies outside the octant, will return
 * the child number of the closest child.)
 *
 */
int child_which(Coord porigin, Coord point) {
  int cnum;
  int i;

  cnum=0;
  for (i=2; i>=0; i--) {
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
