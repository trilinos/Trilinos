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

#include "ifp_c_wrappers.h"
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

struct AZ_IFPACK_STRUCT {                  
/* Data passing structure. This user */
/* defined data structure is used to pass information through   */
/* Aztec and back into the user's matrix-vector product. */
  int nr, nc;
  void *precon, *bmat;
  double *input_vector;
  int user_precon;
  void *user_aux_ptr;
};

typedef struct AZ_IFPACK_STRUCT AZ_IFPACK;

void AZ_ifpack_prec_create(double *x, double *b,
			   int *options, double *params,
			   int *proc_config,
			   AZ_MATRIX *Amat, AZ_PRECOND **Prec);

void AZ_ifpack_iterate(double *x, double *b,
               int *options, double *params, 
               double *status, int *proc_config,
               AZ_MATRIX *Amat);

void AZ_ifpack_precon(double x[], int *, int *,
                     double *, AZ_MATRIX *Amat, AZ_PRECOND *prec);

void AZ_ifpack_prec_destroy(int *options, double *params,
                 int *proc_config, AZ_MATRIX *Amat, AZ_PRECOND *Prec);

void az2ifp_blockmatrix (void **bmat, AZ_MATRIX *Amat);

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif
