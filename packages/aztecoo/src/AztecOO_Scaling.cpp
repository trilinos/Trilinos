
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

#ifdef AZTEC_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif

#include "AztecOO_Scaling.h"
#include "AztecOO.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"

int AZOO_Scale_Jacobi(int action,
                      Epetra_RowMatrix* A,
                      double b[],
                      double x[]);

int AztecOO_scale_epetra(int action,
                         AZ_MATRIX* Amat,
                         int options[],
                         double b[],
                         double x[],
                         int proc_config[])
{
  AztecOO::MatrixData* Data = (AztecOO::MatrixData*)AZ_get_matvec_data(Amat);
  Epetra_RowMatrix* A = (Epetra_RowMatrix*)(Data->A);

  int returnValue = 0;

  if (options[AZ_scaling] == AZ_Jacobi) {
    AZOO_Scale_Jacobi(action, A, b, x);
  }
  else {
    returnValue = -1;
  }

  return(returnValue);
}

int AZOO_Scale_Jacobi(int action,
                      Epetra_RowMatrix* A,
                      double b[],
                      double x[])
{
  if (action == AZ_INVSCALE_SOL || action == AZ_SCALE_SOL) return(0);

  int numMyRows = A->NumMyRows();

  const Epetra_Map& rowmap = A->RowMatrixRowMap();

  //TEMPORARY: we're getting the diagonal vector every time. remember to
  //implement a way to save the vector for reuse.
  Epetra_Vector vec(rowmap);
  int err = A->ExtractDiagonalCopy(vec);
  if (err != 0) {
    return(err);
  }

  double* vec_vals = NULL;
  vec.ExtractView(&vec_vals);

  //now invert each entry of the diagonal vector...
  int i;
  for(i=0; i<numMyRows; ++i) {
    if (fabs(vec_vals[i]) > 1.e-20) vec_vals[i] = 1.0 / vec_vals[i];
    else                             vec_vals[i] = 1.0;
  }

  if (action == AZ_SCALE_MAT_RHS_SOL) {
    A->LeftScale(vec);
  }

  if (action == AZ_SCALE_MAT_RHS_SOL || action == AZ_SCALE_RHS) {
    for(i=0; i<numMyRows; ++i) {
      b[i] *= vec_vals[i];
    }
  }

  if (action == AZ_INVSCALE_RHS) {
    for(i=0; i<numMyRows; ++i) {
      b[i] /= vec_vals[i];
    }
  }

  return(0);
}

