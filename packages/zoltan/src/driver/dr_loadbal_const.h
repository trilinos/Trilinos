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
 * $Name$
 *====================================================================*/

#ifndef _DR_LOADBAL_CONST_H_
#define _DR_LOADBAL_CONST_H_

extern int run_zoltan(int, PROB_INFO_PTR, MESH_INFO_PTR); 
extern int migrate_elements(int, MESH_INFO_PTR, struct LB_Struct *, int,
                            LB_GID *, LB_LID *, int *, int, LB_GID *,
                            LB_LID *, int *);

#endif /* _DR_LOADBAL_CONST_H_ */
