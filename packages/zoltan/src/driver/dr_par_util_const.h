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

#ifndef _DR_PAR_UTIL_CONST_H_
#define _DR_PAR_UTIL_CONST_H_

extern
void print_sync_start (
  int proc,
  int do_print_line
);

extern
void print_sync_end (
  int proc,
  int nprocs,
  int do_print_line
);

extern
void boundary_exchange(
  int vec_len,
  int *send_vec,
  int *recv_vec
);
/* Function prototypes */
#endif /* _DR_PAR_UTIL_CONST_H_ */
