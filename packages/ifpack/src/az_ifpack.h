/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
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
