/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

#ifndef __UTIL_CONST_H
#define __UTIL_CONST_H

#ifndef lint
static char *cvs_utilconsth_id = "$Id$";
#endif

extern void   set_method(double method_number);
extern void   *my_malloc(int size);
extern int    in_box(COORD pt, COORD lower, COORD upper);
extern void   bounds_to_origin_size(COORD min, COORD max,
				    COORD origin, double size[3]);
extern void   bounds_to_origin(COORD min, COORD max, 
			       COORD origin);
extern void   child_bounds_wrapper(pOctant oct, COORD cmin[], COORD cmax[]);
extern void   child_bounds(COORD pmin, COORD pmax, COORD porigin,
			   int cnum, COORD cmin, COORD cmax);

extern int    compare(unsigned *x, unsigned *y);
extern int    hilbert_bounds(COORD min, COORD max, COORD cmin[], COORD cmax[]);
extern int    hilbert2d_bounds(COORD min, COORD max, 
			       COORD cmin[], COORD cmax[]);

extern int    child_which_wrapper(pOctant oct, COORD point);
extern int    child_which(COORD origin, COORD point);
extern double dist_point_box(COORD point, COORD min, COORD max);
extern int    convert_to_gray(int input);
extern int    convert_from_gray(int input);

extern int    convert_to_hilbert(int n, int o);
extern int    convert_from_hilbert(int n, int o);
extern int    child_orientation(int o, int cnum);
extern int    change_to_hilbert2d(COORD min, COORD max, 
				  COORD origin, int cnum);
extern int    change_to_hilbert(COORD min, COORD max, COORD origin,
				int cnum);

#endif
