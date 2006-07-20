/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

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
    (void) AZ_printf_err( "direct: ifail is not zero (%d)\n", *ifail);
    switch(*ifail) {
    case 4:
      (void) AZ_printf_err( "Large growth factor\n");
    break;
    case 3:
      (void) AZ_printf_err( "Matrix may be singular\n");
    break;
    case 5:
      (void) AZ_printf_err( "Allocated space not large enough\n");
    break;
    case 6:
      (void) AZ_printf_err( "Allocated space not large enough\n");
    break;
    case 7:
      (void) AZ_printf_err( "Either the matrix may be singular\n"
                     "or the drop tolerance may be too high\n");
    break;
    case 8:
      (void) AZ_printf_err( "Either the matrix may be singular\n"
                     "or the drop tolerance may be too high\n");
    break;
    case 11:
      (void) AZ_printf_err( "two elements in the same (i,j) position\n");
    break;
    case 12:
      (void) AZ_printf_err( "System has less than two rows\n");
    break;
    case 17:
      (void) AZ_printf_err( "A row without nonzero elements found\n");
    break;
    case 18:
      (void) AZ_printf_err( "A column without nonzero elements found\n");
    break;
    case 24:
      (void) AZ_printf_err( "A column index exceeds matrix dimension \n");
    break;
    case 25:
      (void) AZ_printf_err( "A row index exceeds matrix dimension \n");
    break;
    default:
      break;
    }
    (void) AZ_printf_err( "        Check y12m manual for more information.\n");
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
