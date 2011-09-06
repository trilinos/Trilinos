
/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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

//---------------------------------------------------------------
//Prototypes for utility functions that are used internally by
//AztecOO_scale_epetra.

int AZOO_Scale(int action,
               Epetra_RowMatrix* A,
               double b[],
               double x[],
               int options[],
               AZ_SCALING* scaling);

Epetra_Vector* AZOO_create_scaling_vector(Epetra_RowMatrix* A,
                                          int scaling_type);

//---------------------------------------------------------------
//Scaling function which the AztecOO class registers on the
//AZ_SCALING struct as a call-back before calling AZ_iterate.

int AztecOO_scale_epetra(int action,
                         AZ_MATRIX* Amat,
                         int options[],
                         double b[],
                         double x[],
                         int proc_config[],
                         AZ_SCALING* scaling)
{
  (void)proc_config;
  AztecOO::MatrixData* Data = (AztecOO::MatrixData*)AZ_get_matvec_data(Amat);
  Epetra_RowMatrix* A = (Epetra_RowMatrix*)(Data->A);

  int returnValue = 0;

  if (options[AZ_scaling] == AZ_Jacobi ||
      options[AZ_scaling] == AZ_row_sum ||
      options[AZ_scaling] == AZ_sym_diag) {
    returnValue = AZOO_Scale(action, A, b, x, options, scaling);
  }
  else {
    returnValue = -1;
  }

  return(returnValue);
}

int AZOO_Scale(int action,
               Epetra_RowMatrix* A,
               double b[],
               double x[],
               int options[],
               AZ_SCALING* scaling)
{
  (void)x;
  //This function performs Jacobi, row-sum or sym_diag scaling, and
  //basically mirrors the functionality provided by the
  //functions AZ_block_diagonal_scaling, AZ_row_sum_scaling and
  //AZ_sym_diagonal_scaling in the file az_scaling.c.

  if (action == AZ_DESTROY_SCALING_DATA) {
    if (scaling->scaling_data != 0) {
      Epetra_Vector* vec = (Epetra_Vector*)(scaling->scaling_data);
      delete vec;
      scaling->scaling_data = 0;
    }
    return(0);
  }

  if (options[AZ_scaling] != AZ_sym_diag) {
    if (action == AZ_INVSCALE_SOL || action == AZ_SCALE_SOL) {
      return(0);
    }
  }

  int numMyRows = A->NumMyRows();

  //const Epetra_Map& rowmap = A->RowMatrixRowMap();

  Epetra_Vector* vec = NULL;

  //Either get the data stored in scaling->scaling_data, or
  //create a new scaling vector, depending on the value of
  //options[AZ_pre_calc] and action.

  if (options[AZ_pre_calc] == AZ_reuse) {
    if (scaling->scaling_data == NULL) {
      if (options[AZ_output] != AZ_none) {
        cerr << "AZOO_Scale ERROR, AZ_reuse requested, but"
           << " scaling->scaling_data==NULL"<<endl;
      }
      return(-1);
    }

    vec = (Epetra_Vector*)(scaling->scaling_data);
  }
  else {
    if (action == AZ_SCALE_MAT_RHS_SOL) {
      vec = AZOO_create_scaling_vector(A, options[AZ_scaling]);
      if (vec == NULL) {
        if (options[AZ_output] != AZ_none) {
          cerr << "AZOO_create_scaling_vector ERROR"<<endl;
        }
        return(-1);
      }
      if (scaling->scaling_data != NULL) {
        Epetra_Vector* oldvec =
          static_cast<Epetra_Vector*>(scaling->scaling_data);
        delete oldvec;
      }
      scaling->scaling_data = (void*)vec;
    }
    else {
      if (scaling->scaling_data != NULL) {
        vec = (Epetra_Vector*)(scaling->scaling_data);
      }

      if (scaling->scaling_data == NULL || vec == NULL) {
        if (options[AZ_output] != AZ_none) {
          cerr << "AZOO_Scale ERROR, vec == NULL or"
             << " scaling->scaling_data==NULL"<<endl;
        }
        return(-1);
      }
    }
  }

  double* vec_vals = NULL;
  vec->ExtractView(&vec_vals);

  if (action == AZ_SCALE_MAT_RHS_SOL && options[AZ_pre_calc] <= AZ_recalc) {
    A->LeftScale(*vec);

    if (options[AZ_scaling] == AZ_sym_diag) {
      A->RightScale(*vec);
    }
  }

  if (action == AZ_SCALE_MAT_RHS_SOL) {
    if (options[AZ_scaling] == AZ_sym_diag) {
      for(int i=0; i<numMyRows; ++i) {
        b[i] *= vec_vals[i];
        x[i] /= vec_vals[i];
      }
    }
    else {
      for(int i=0; i<numMyRows; ++i) {
        b[i] *= vec_vals[i];
      }
    }
  }

  if (action == AZ_SCALE_SOL) {
    for(int i=0; i<numMyRows; ++i) {
      x[i] /= vec_vals[i];
    }
  }

  if (action == AZ_INVSCALE_SOL) {
    for(int i=0; i<numMyRows; ++i) {
      x[i] *= vec_vals[i];
    }
  }

  if (action == AZ_INVSCALE_RHS) {
    for(int i=0; i<numMyRows; ++i) {
      b[i] /= vec_vals[i];
    }
  }

  if (action == AZ_SCALE_RHS) {
    for(int i=0; i<numMyRows; ++i) {
      b[i] *= vec_vals[i];
    }
  }

  return(0);
}

Epetra_Vector* AZOO_create_scaling_vector(Epetra_RowMatrix* A,
                                          int scaling_type)
{
  //This function creates a new Epetra_Vector, and fills it
  //with the inverse of A's diagonal elements or row-sums,
  //depending on the value of scaling_type.

  Epetra_Vector* vec = new Epetra_Vector(A->RowMatrixRowMap());

  if (scaling_type == AZ_Jacobi || scaling_type == AZ_sym_diag) {
    int err = A->ExtractDiagonalCopy(*vec);
    if (err != 0) {
      delete vec; vec = 0;
      return(vec);
    }

    double* vec_vals = NULL;
    vec->ExtractView(&vec_vals);

    //now invert each entry of the diagonal vector
    for(int i=0; i<A->RowMatrixRowMap().NumMyElements(); ++i) {
      double sval = fabs(vec_vals[i]);
      if (sval > Epetra_MinDouble) {
        vec_vals[i] = scaling_type == AZ_sym_diag ?
                    1.0/sqrt(sval) : 1.0/sval;
      }
      else {
        vec_vals[i] = 1.0;
      }
    }
  }
  else if (scaling_type == AZ_row_sum) {
    int err = A->InvRowSums(*vec);
    if (err != 0) {
      if (err == 1) {
        //err==1 indicates a zero row-sum was found. What should
        //we do in this case? For now, ignore it...
      }
      else {
        delete vec; vec = 0;
        return(vec);
      }
    }
  }
  else {
    delete vec;
    vec = 0;
  }

  return(vec);
}

