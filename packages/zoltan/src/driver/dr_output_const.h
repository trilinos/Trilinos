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
#ifndef lint
static char *cvs_dr_output_const = "$Id$";
#endif

#ifndef _DR_OUTPUT_CONST_H_
#define _DR_OUTPUT_CONST_H_

void print_distributed_mesh(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  ELEM_INFO *elements
);

#endif /* _DR_OUTPUT_CONST_H_ */
