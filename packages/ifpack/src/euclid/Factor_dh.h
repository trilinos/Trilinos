/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef FACTOR_DH
#define FACTOR_DH

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "euclid_common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  struct _factor_dh
  {
    /* dimensions of local rectangular submatrix; global matrix is n*n */
    int m, n;

    int id;			/* this subdomain's id after reordering */
    int beg_row;		/* global number of 1st locally owned row */
    int first_bdry;		/* local number of first boundary row */
    int bdry_count;		/* m - first_boundary */

    /* if true, factorization was block jacobi, in which case all
       column indices are zero-based; else, they are global.
     */
    bool blockJacobi;

    /* sparse row-oriented storage for locally owned submatrix */
    int *rp;
    int *cval;
    REAL_DH *aval;
    int *fill;
    int *diag;
    int alloc;			/* currently allocated length of cval, aval, and fill arrays */

    /* used for PILU solves (Apply) */
    int num_recvLo, num_recvHi;
    int num_sendLo, num_sendHi;	/* used in destructor */
    double *work_y_lo;		/* recv values from lower nabors; also used as
				   work vector when solving Ly=b for y.
				 */
    double *work_x_hi;		/* recv values from higher nabors; also used as
				   work vector when solving Ux=y for x.
				 */
    double *sendbufLo, *sendbufHi;
    int *sendindLo, *sendindHi;
    int sendlenLo, sendlenHi;
    bool solveIsSetup;
    Numbering_dh numbSolve;

    MPI_Request recv_reqLo[MAX_MPI_TASKS], recv_reqHi[MAX_MPI_TASKS];	/* used for persistent comms */
    MPI_Request send_reqLo[MAX_MPI_TASKS], send_reqHi[MAX_MPI_TASKS];	/* used for persistent comms */
    MPI_Request requests[MAX_MPI_TASKS];
    MPI_Status status[MAX_MPI_TASKS];

    bool debug;
  };

  extern void Factor_dhCreate (Factor_dh * mat);
  extern void Factor_dhDestroy (Factor_dh mat);

  extern void Factor_dhTranspose (Factor_dh matIN, Factor_dh * matOUT);

  extern void Factor_dhInit (void *A, bool fillFlag, bool avalFlag,
			     double rho, int id, int beg_rowP, Factor_dh * F);

  extern void Factor_dhReallocate (Factor_dh F, int used, int additional);
  /* ensures fill, cval, and aval arrays can accomodate
     at least "c" additional entrie
   */

  /* adopted from ParaSails, by Edmond Chow */
  extern void Factor_dhSolveSetup (Factor_dh mat, SubdomainGraph_dh sg);


  extern void Factor_dhSolve (double *rhs, double *lhs, Euclid_dh ctx);
  extern void Factor_dhSolveSeq (double *rhs, double *lhs, Euclid_dh ctx);

  /* functions for monitoring stability */
  extern double Factor_dhCondEst (Factor_dh mat, Euclid_dh ctx);
  extern double Factor_dhMaxValue (Factor_dh mat);
  extern double Factor_dhMaxPivotInverse (Factor_dh mat);

  extern int Factor_dhReadNz (Factor_dh mat);
  extern void Factor_dhPrintTriples (Factor_dh mat, char *filename);

  extern void Factor_dhPrintGraph (Factor_dh mat, char *filename);
  /* seq only */


  extern void Factor_dhPrintDiags (Factor_dh mat, FILE * fp);
  extern void Factor_dhPrintRows (Factor_dh mat, FILE * fp);
  /* prints local matrix to logfile, if open */

#ifdef __cplusplus
}
#endif
#endif
