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

#ifndef __PAR_CONST_H
#define __PAR_CONST_H

#ifndef lint
static char *cvs_par_const_h = "$Id$";
#endif

#ifdef LB_MPI
#include "mpi.h"
#endif

extern int LB_Proc;
extern int LB_Num_Proc;

extern void LB_print_sync_start(int);
extern void LB_print_sync_end(int);

#endif
