/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef _DR_LOADBAL_CONST_H_
#define _DR_LOADBAL_CONST_H_

extern int run_zoltan(int, PROB_INFO_PTR, MESH_INFO_PTR); 
extern int migrate_elements(int, MESH_INFO_PTR, struct LB_Struct *, 
                            int, int, int, 
                            LB_ID_PTR, LB_ID_PTR, int *, int, LB_ID_PTR,
                            LB_ID_PTR, int *);

#endif /* _DR_LOADBAL_CONST_H_ */
