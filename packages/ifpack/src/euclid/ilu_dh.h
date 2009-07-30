/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef ILU_MPI_DH
#define ILU_MPI_DH

#include "euclid_common.h"

#ifdef __cplusplus
extern "C" {
#endif

void reallocate_private(int row, int newEntries, int *nzHave,
                int **rp, int **cval, float **aval, double **avalD, int **fill);

extern void ilu_mpi_pilu(Euclid_dh ctx);
  /* driver for comms intermingled with factorization */


extern void iluk_mpi_pilu(Euclid_dh ctx);
  /* the factorization algorithm */

extern void compute_scaling_private(int row, int len, double *AVAL, Euclid_dh ctx);

extern void iluk_mpi_bj(Euclid_dh ctx);

extern void iluk_seq(Euclid_dh ctx);
extern void iluk_seq_block(Euclid_dh ctx);
  /* for sequential or parallel block jacobi.  If used
     for block jacobi, column indices are referenced to 0
     on return; make sure and add beg_row to these values
     before printing the matrix!

     1st version is for single precision, 2nd is for double.
   */

extern void ilut_seq(Euclid_dh ctx);

#ifdef __cplusplus
}
#endif

#endif

