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

#ifndef EXTERNAL_ROWS_DH_H
#define EXTERNAL_ROWS_DH_H

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

  extern void ExternalRows_dhCreate (ExternalRows_dh * er);
  extern void ExternalRows_dhDestroy (ExternalRows_dh er);
  extern void ExternalRows_dhInit (ExternalRows_dh er, Euclid_dh ctx);
  extern void ExternalRows_dhRecvRows (ExternalRows_dh extRows);
  extern void ExternalRows_dhSendRows (ExternalRows_dh extRows);
  extern void ExternalRows_dhGetRow (ExternalRows_dh er, int globalRow,
				     int *len, int **cval, int **fill,
				     REAL_DH ** aval);

  struct _extrows_dh
  {
    SubdomainGraph_dh sg;	/* not owned! */
    Factor_dh F;		/* not owned! */

    MPI_Status status[MAX_MPI_TASKS];
    MPI_Request req1[MAX_MPI_TASKS];
    MPI_Request req2[MAX_MPI_TASKS];
    MPI_Request req3[MAX_MPI_TASKS];
    MPI_Request req4[MAX_MPI_TASKS];
    MPI_Request cval_req[MAX_MPI_TASKS];
    MPI_Request fill_req[MAX_MPI_TASKS];
    MPI_Request aval_req[MAX_MPI_TASKS];

    /*------------------------------------------------------------------------
     *  data structures for receiving, storing, and accessing external rows 
     *  from lower-ordered nabors
     *------------------------------------------------------------------------*/
    /* for reception of row counts, row numbers, and row lengths: */
    int rcv_row_counts[MAX_MPI_TASKS];	/* P_i will send rcv_row_counts[i] rows */
    int rcv_nz_counts[MAX_MPI_TASKS];	/* P_i's rows contain rcv_nz_counts[i] nonzeros */
    int *rcv_row_lengths[MAX_MPI_TASKS];	/* rcv_row_lengths[i][] lists the length of each row */
    int *rcv_row_numbers[MAX_MPI_TASKS];	/* rcv_row_lengths[i][] lists the length of each row */

    /* for reception of the actual rows: */
    int *cvalExt;
    int *fillExt;
    REAL_DH *avalExt;

    /* table for accessing the rows */
    Hash_dh rowLookup;

    /*--------------------------------------------------------------------------
     *  data structures for sending boundary rows to higher-ordered nabors
     *--------------------------------------------------------------------------*/
    /* for sending row counts, numbers, and lengths: */
    int *my_row_counts;		/* my_row_counts[i] = nzcount in upper tri portion o */
    int *my_row_numbers;	/* my_row_numbers[i] = global row number of local ro */

    /* for sending the actual rows: */
    int nzSend;			/* total entries in upper tri portions of bdry rows */
    int *cvalSend;
    int *fillSend;
    REAL_DH *avalSend;

    bool debug;
  };
#ifdef __cplusplus
}
#endif

#endif
