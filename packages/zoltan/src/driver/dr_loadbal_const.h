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

extern int run_zoltan(int, PROB_INFO_PTR, ELEM_INFO **);
extern int migrate_elements(int, ELEM_INFO **, struct LB_Struct *, int,
                            LB_GID *, LB_LID *, int *, int, LB_GID *,
                            LB_LID *, int *);

#endif /* _DR_LOADBAL_CONST_H_ */
