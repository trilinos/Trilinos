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

#ifndef _DR_OUTPUT_CONST_H_
#define _DR_OUTPUT_CONST_H_

#include "dr_input_const.h"

extern void print_distributed_mesh(
  int Proc,
  int Num_Proc,
  ELEM_INFO *elements
);

extern int output_results(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  ELEM_INFO elements[]);

#endif /* _DR_OUTPUT_CONST_H_ */
