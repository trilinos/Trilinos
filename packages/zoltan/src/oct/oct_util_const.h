/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __OCT_UTIL_CONST_H
#define __OCT_UTIL_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "octree_const.h"
extern int    Zoltan_Oct_get_child_dir(OCT_Global_Info *OCT_info, int dir, int cnum);
extern int    Zoltan_Oct_convert_idx_from_map(OCT_Global_Info *OCT_info, int dir, int cnum);

extern void   Zoltan_Oct_set_method(OCT_Global_Info *OCT_info,int method_number);
extern int    Zoltan_Oct_in_box(OCT_Global_Info *OCT_info,COORD pt, COORD lower, COORD upper);
extern int    Zoltan_Oct_in_box_closure(OCT_Global_Info *OCT_info,COORD pt, COORD lower, COORD upper);
extern void   Zoltan_Oct_bounds_to_origin_size(COORD min, COORD max,
				       COORD origin, double size[3]);
extern void   Zoltan_Oct_bounds_to_origin(COORD min, COORD max, 
			          COORD origin);
extern void   Zoltan_Oct_child_bounds_wrapper(OCT_Global_Info *OCT_info,pOctant oct, COORD cmin[], COORD cmax[]);
extern void   Zoltan_Oct_child_bounds(COORD pmin, COORD pmax, COORD porigin,
			      int cnum, COORD cmin, COORD cmax);

extern int    Zoltan_Oct_child_which_wrapper(OCT_Global_Info *OCT_info,pOctant oct, COORD point);
extern int    Zoltan_Oct_child_which(OCT_Global_Info *OCT_info,COORD origin, COORD point);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
