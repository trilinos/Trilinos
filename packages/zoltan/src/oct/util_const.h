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

extern void   LB_set_method(int method_number);
extern int    LB_in_box(COORD pt, COORD lower, COORD upper);
extern void   LB_bounds_to_origin_size(COORD min, COORD max,
				       COORD origin, double size[3]);
extern void   LB_bounds_to_origin(COORD min, COORD max, 
			          COORD origin);
extern void   LB_child_bounds_wrapper(pOctant oct, COORD cmin[], COORD cmax[]);
extern void   LB_child_bounds(COORD pmin, COORD pmax, COORD porigin,
			      int cnum, COORD cmin, COORD cmax);

extern int    LB_child_which_wrapper(pOctant oct, COORD point);
extern int    LB_child_which(COORD origin, COORD point);
extern double LB_dist_point_box(COORD point, COORD min, COORD max);
extern int    LB_convert_to_gray(int input);
extern int    LB_convert_from_gray(int input);

extern int    LB_convert_to_hilbert(int n, int o);
extern int    LB_convert_from_hilbert(int n, int o);
extern int    LB_child_orientation(int o, int cnum);
extern int    LB_change_to_hilbert2d(COORD min, COORD max, 
				     COORD origin, int cnum);
extern int    LB_change_to_hilbert(COORD min, COORD max, COORD origin,
				   int cnum);

#endif
