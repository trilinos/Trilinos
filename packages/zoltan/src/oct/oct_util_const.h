/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __OCT_UTIL_CONST_H
#define __OCT_UTIL_CONST_H

#include "octant_const.h"

extern void   LB_set_method(OCT_Global_Info *OCT_info,int method_number);
extern int    LB_in_box(OCT_Global_Info *OCT_info,COORD pt, COORD lower, COORD upper);
extern void   LB_bounds_to_origin_size(COORD min, COORD max,
				       COORD origin, double size[3]);
extern void   LB_bounds_to_origin(COORD min, COORD max, 
			          COORD origin);
extern void   LB_child_bounds_wrapper(OCT_Global_Info *OCT_info,pOctant oct, COORD cmin[], COORD cmax[]);
extern void   LB_child_bounds(COORD pmin, COORD pmax, COORD porigin,
			      int cnum, COORD cmin, COORD cmax);

extern int    LB_child_which_wrapper(OCT_Global_Info *OCT_info,pOctant oct, COORD point);
extern int    LB_child_which(OCT_Global_Info *OCT_info,COORD origin, COORD point);
extern int    LB_convert_to_gray(int input);

extern int    LB_child_orientation(int o, int cnum);
extern int    LB_change_to_hilbert2d(OCT_Global_Info *OCT_info,COORD min, COORD max, 
				     COORD origin, int cnum);
extern int    LB_change_to_hilbert(OCT_Global_Info *OCT_info,COORD min, COORD max, COORD origin,
				   int cnum);
extern void LB_OCT_Free_Structure(LB *lb);

#endif
