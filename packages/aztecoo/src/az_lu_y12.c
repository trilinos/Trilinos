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

/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include "az_aztec.h"
#ifdef HAVE_AZLU
#include "az_y12m_wrappers.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_fact_lu(double b[], AZ_MATRIX *A_overlapped, double *aflag, 
                double *pivot, int *rnr, int *ha, int *iflag, int *z, 
                int *ifail, int *nn, int *iha, int *n)

/*******************************************************************************

  Routines which call the correct sequence of routines for the sparse direct
  solver: y12m

  The subroutine 'AZ_factorize' first factorizes the matrix.
  The second subroutine 'backsolve' uses precomputed factors to
  perform a back solve.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   Aztec User's Guide).

  aflag:

  pivot:

  b:

  snr, rnr:

  ha:

  iflag:

  z:

  ifail:

  nn:

  n:

  iha:

  nn1:

*******************************************************************************/
{
   double *val;
   int *snr,  *nn1;

   val   = A_overlapped->val;
   snr   = A_overlapped->bindx;
   nn1   = nn;

  *ifail = 0;

  Y12MBF_F77(n, z, val, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail);

  if (*ifail == 0) 
     Y12MCF_F77(n, z, val, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, 
             ifail);


  if (*ifail != 0) {
    (void) fprintf(stderr, "direct: ifail is not zero (%d)\n", *ifail);
    switch(*ifail) {
    case 4:
      (void) fprintf(stderr, "Large growth factor\n");
    break;
    case 3:
      (void) fprintf(stderr, "Matrix may be singular\n");
    break;
    case 5:
      (void) fprintf(stderr, "Allocated space not large enough\n");
    break;
    case 6:
      (void) fprintf(stderr, "Allocated space not large enough\n");
    break;
    case 7:
      (void) fprintf(stderr, "Either the matrix may be singular\n"
                     "or the drop tolerance may be too high\n");
    break;
    case 8:
      (void) fprintf(stderr, "Either the matrix may be singular\n"
                     "or the drop tolerance may be too high\n");
    break;
    case 11:
      (void) fprintf(stderr, "two elements in the same (i,j) position\n");
    break;
    case 12:
      (void) fprintf(stderr, "System has less than two rows\n");
    break;
    case 17:
      (void) fprintf(stderr, "A row without nonzero elements found\n");
    break;
    case 18:
      (void) fprintf(stderr, "A column without nonzero elements found\n");
    break;
    case 24:
      (void) fprintf(stderr, "A column index exceeds matrix dimension \n");
    break;
    case 25:
      (void) fprintf(stderr, "A row index exceeds matrix dimension \n");
    break;
    default:
      break;
    }
    (void) fprintf(stderr, "        Check y12m manual for more information.\n");
      exit(-1);
  }

} /* AZ_factorize*/

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_backsolve(double val[], double pivot[], double b[], int snr[], int ha[],
                  int iflag[], int *ifail, int *nn, int *n, int *iha)

/*******************************************************************************

  Routine which calls the correct sequence of routines for the sparse direct
  solver: y12m

  The subroutine 'backsolve' performs a backsolve using precomputed factors.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   Aztec User's Guide).

  pivot:

  b:

  snr:

  ha:

  iflag:

  ifail:

  nn:

  n:

  iha:

*******************************************************************************/

{

  Y12MDF_F77(n, val, nn, b, pivot, snr, ha, iha, iflag, ifail);

} /* backsolve */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_msr2lu(int oldN, AZ_MATRIX *A_overlapped, int *rnr)

/*******************************************************************************

  This routine converts an MSR matrix into a format suitable for the lu solver.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  oldN:

  val:             Array containing the nonzero entries of the matrix (see file
                   Aztec User's Guide).

  newa:

  snr, rnr:

  z:

  bindx:           MSR index vector.

*******************************************************************************/

{

  /* local variables */

  int i, j, start, N_nz;
  int *bindx;
  double *val;

  /**************************** execution begins ******************************/

  bindx = A_overlapped->bindx;
  val   = A_overlapped->val;

  start = bindx[0];
  N_nz = bindx[oldN];

  for (i = 0 ; i < oldN ; i++ ) {
     for ( j = bindx[i] ; j < bindx[i+1] ; j++ ) rnr[j-1]  = i + 1;
     bindx[i] = i + 1;
     rnr[i]   = i + 1;
  }
  for (i = start ; i < N_nz; i++ ) {
     bindx[i-1]  = bindx[i] + 1;
     val[i-1] = val[i];
  }

} /* AZ_msr2lu */

#endif /* HAVE_AZLU */
