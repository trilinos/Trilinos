/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * Zoltan is distributed under the GNU Lesser General Public License 2.1.    * 
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __LB_UTIL_CONST_H
#define __LB_UTIL_CONST_H

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

extern void LB_perform_error_checking(LB *);
extern void LB_Get_Obj_List(LB *, LB_GID *, LB_LID *, int, float *, int *);
extern int LB_pad_for_alignment(int);
extern unsigned int LB_hashf(LB_GID, int);
extern int LB_clean_string(char *, char **);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#endif
