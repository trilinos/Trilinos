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

#ifndef __MSG_CONST_H
#define __MSG_CONST_H

#ifndef lint
static char *cvs_msgconsth_id = "$Id$";
#endif


extern int    msg_int_scan(MPI_Comm, int, int value);
extern float  msg_float_scan(MPI_Comm, int, float value);

#endif
