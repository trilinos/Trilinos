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

#include <mpi.h>

/*
 * Function prototype for utility to find the median
 * of an array of weighted coordinates.
 */
extern int LB_find_median(
  double *dots,         /* array of coordinates                              */
  double *wgts,         /* array of weights associated with dots             */
  int *dotmark,         /* returned list of which side of the median
                           each dot is on:
                                0 - dot is < valuehalf
                                1 - dot is > valuehalf                       */
  int dotnum,           /* number of dots (length of three previous arrays   */
  int proc,             /* this proc number (rank)                           */
  double fractionlo,    /* fraction of weight that should be in bottom half  */
  MPI_Comm local_comm,  /* MPI communicator on which to find median          */
  double *valuehalf,    /* on entry - first guess at median (if first_guess set)                           on exit - the median value                        */
  int first_guess,      /* if set, use value in valuehalf as first guess     */
  int *counter          /* returned for stats, # of median interations       */
);

extern void LB_print_sync_start(LB *, int);
extern void LB_print_sync_end(LB *, int);
extern void LB_Print_Stats (LB *, double, char *);

#endif
