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

#ifndef MAT_DH_DH
#define MAT_DH_DH

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "euclid_common.h"

  /* this stuff for experimental internal timing */
#define MAT_DH_BINS      10
#define MATVEC_TIME       0	/* time to actually perform matvec */
#define MATVEC_MPI_TIME   1	/* time for comms + vector copying needed */
#define MATVEC_MPI_TIME2  5	/* time for comms, + vector copying needed */
#define MATVEC_TOTAL_TIME 2	/* MATVEC_TIME+MATVEC_MPI_TIME */
#define MATVEC_RATIO      3	/* computation/communication ratio */
#define MATVEC_WORDS      4	/* total words sent to other procs. */

#ifdef __cplusplus
extern "C"
{
#endif

  struct _mat_dh
  {
    int m, n;			/* dimensions of local rectangular submatrix;
				 * the global matrix is n by n.
				 */
    int beg_row;		/* global number of 1st locally owned row */
    int bs;			/* block size */

    /* sparse row-oriented storage for locally owned submatrix */
    int *rp;
    int *len;			/* length of each row; only used for MPI triangular solves */
    int *cval;
    int *fill;
    int *diag;
    double *aval;
    bool owner;			/* for MPI triangular solves */

    /* working space for getRow */
    int len_private;
    int rowCheckedOut;
    int *cval_private;
    double *aval_private;

    /* row permutations to increase positive definiteness */
    int *row_perm;

    /* for timing matvecs in experimental studies */
    double time[MAT_DH_BINS];
    double time_max[MAT_DH_BINS];
    double time_min[MAT_DH_BINS];
    bool matvec_timing;

    /* used for MatVecs */
    int num_recv;
    int num_send;		/* used in destructor */
    MPI_Request *recv_req;
    MPI_Request *send_req;
    double *recvbuf, *sendbuf;
    int *sendind;
    int sendlen;
    int recvlen;
    bool matvecIsSetup;
    Numbering_dh numb;
    MPI_Status *status;

    bool debug;
  };

  extern void Mat_dhCreate (Mat_dh * mat);
  extern void Mat_dhDestroy (Mat_dh mat);

  extern void Mat_dhTranspose (Mat_dh matIN, Mat_dh * matOUT);
  extern void Mat_dhMakeStructurallySymmetric (Mat_dh A);

  /* adopted from ParaSails, by Edmond Chow */
  extern void Mat_dhMatVecSetup (Mat_dh mat);
  extern void Mat_dhMatVecSetdown (Mat_dh mat);

/*========================================================================*/
/* notes: if not compiled with OpenMP, Mat_dhMatVec() and Mat_dhMatVec_omp()
          perform identically; similarly for Mat_dhMatVec_uni()
          and Mat_dhMatVec_uni_omp()
*/

  extern void Mat_dhMatVec (Mat_dh mat, double *lhs, double *rhs);
  /* unthreaded MPI version */

  extern void Mat_dhMatVec_omp (Mat_dh mat, double *lhs, double *rhs);
  /* OpenMP/MPI version */

  extern void Mat_dhMatVec_uni (Mat_dh mat, double *lhs, double *rhs);
  /* unthreaded, single-task version */

  extern void Mat_dhMatVec_uni_omp (Mat_dh mat, double *lhs, double *rhs);
  /* OpenMP/single primary task version */


  extern int Mat_dhReadNz (Mat_dh mat);

  /* for next five, SubdomainGraph_dh() may be NULL; if not null,
     caller must ensure it has been properly initialized;
     if not null, matrix is permuted before printing.

     note: use "-matlab" when calling Mat_dhPrintTriples, to
     insert small value in place of 0.

     Mat_dhPrintCSR only implemented for single cpu, no reordering.
   */
  extern void Mat_dhPrintGraph (Mat_dh mat, SubdomainGraph_dh sg, FILE * fp);
  extern void Mat_dhPrintRows (Mat_dh mat, SubdomainGraph_dh sg, FILE * fp);

  extern void Mat_dhPrintCSR (Mat_dh mat, SubdomainGraph_dh sg,
			      char *filename);
  extern void Mat_dhPrintTriples (Mat_dh mat, SubdomainGraph_dh sg,
				  char *filename);
  extern void Mat_dhPrintBIN (Mat_dh mat, SubdomainGraph_dh sg,
			      char *filename);

  extern void Mat_dhReadCSR (Mat_dh * mat, char *filename);
  extern void Mat_dhReadTriples (Mat_dh * mat, int ignore, char *filename);
  extern void Mat_dhReadBIN (Mat_dh * mat, char *filename);


  extern void Mat_dhPermute (Mat_dh Ain, int *pIN, Mat_dh * Bout);
  /* for single cpu only! */

  extern void Mat_dhFixDiags (Mat_dh A);
  /* inserts diagonal if not explicitly present;
     sets diagonal value in row i to sum of absolute
     values of all elts in row i.
   */

  extern void Mat_dhPrintDiags (Mat_dh A, FILE * fp);

  extern void Mat_dhGetRow (Mat_dh B, int globalRow, int *len, int **ind,
			    double **val);
  extern void Mat_dhRestoreRow (Mat_dh B, int row, int *len, int **ind,
				double **val);

  /* partition matrix into "k" blocks.  User must free storage. */
  extern void Mat_dhPartition (Mat_dh mat, int k, int **beg_rowOUT,
			       int **row_countOUT, int **n2oOUT,
			       int **o2nOUT);




  extern void Mat_dhZeroTiming (Mat_dh mat);
  extern void Mat_dhReduceTiming (Mat_dh mat);


  extern void Mat_dhRowPermute (Mat_dh);

  extern void dldperm (int job, int n, int nnz, int colptr[], int adjncy[],
		       double nzval[], int *perm, double u[], double v[]);


#ifdef __cplusplus
}
#endif
#endif
