
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

//---------------------------------------------------------------
//Prototypes for utility functions that are used internally by
//AztecOO_scale_epetra.

int AZOO_Scale_Jacobi_or_row_sum(int action,
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
  AztecOO::MatrixData* Data = (AztecOO::MatrixData*)AZ_get_matvec_data(Amat);
  Epetra_RowMatrix* A = (Epetra_RowMatrix*)(Data->A);

  int returnValue = 0;

  if (options[AZ_scaling] == AZ_Jacobi ||
      options[AZ_scaling] == AZ_row_sum) {
    returnValue = AZOO_Scale_Jacobi_or_row_sum(action, A, b, x,
                                               options, scaling);
  }
  else {
    returnValue = -1;
  }

  return(returnValue);
}

int AZOO_Scale_Jacobi_or_row_sum(int action,
                                 Epetra_RowMatrix* A,
                                 double b[],
                                 double x[],
                                 int options[],
                                 AZ_SCALING* scaling)
{
  //This function performs either Jacobi or row-sum scaling, and
  //basically mirrors the functionality provided by the
  //functions AZ_block_diagonal_scaling and AZ_row_sum_scaling
  //in the file az_scaling.c.

  if (action == AZ_INVSCALE_SOL || action == AZ_SCALE_SOL) return(0);

  if (action == AZ_DESTROY_SCALING_DATA) {
    if (scaling->scaling_data != 0) {
      Epetra_Vector* vec = (Epetra_Vector*)(scaling->scaling_data);
      delete vec;
      scaling->scaling_data = 0;
    }
  }

  int numMyRows = A->NumMyRows();

  //const Epetra_Map& rowmap = A->RowMatrixRowMap();

  Epetra_Vector* vec = NULL;

  //Either get the data stored in scaling->scaling_data, or
  //create a new scaling vector, depending on the value of
  //options[AZ_pre_calc].

  if (options[AZ_pre_calc] == AZ_reuse) {
    if (scaling->scaling_data == NULL) {
      if (options[AZ_output] != AZ_none) {
        cerr << "AZOO_Scale_Jacobi ERROR, AZ_reuse requested, but"
           << " scaling->scaling_data==NULL"<<endl;
      }
      return(-1);
    }

    vec = (Epetra_Vector*)(scaling->scaling_data);
  }
  else {
    vec = AZOO_create_scaling_vector(A, options[AZ_scaling]);
    if (vec == NULL) {
      if (options[AZ_output] != AZ_none) {
        cerr << "AZOO_create_scaling_vector ERROR"<<endl;
      }
      return(-1);
    }
  }

  double* vec_vals = NULL;
  vec->ExtractView(&vec_vals);

  if (action == AZ_SCALE_MAT_RHS_SOL) {
    A->LeftScale(*vec);
  }

  if (action == AZ_SCALE_MAT_RHS_SOL || action == AZ_SCALE_RHS) {
    for(int i=0; i<numMyRows; ++i) {
      b[i] *= vec_vals[i];
    }
  }

  if (action == AZ_INVSCALE_RHS) {
    for(int i=0; i<numMyRows; ++i) {
      b[i] /= vec_vals[i];
    }
  }

  //if options[AZ_keep_info]==1, store the scaling vector
  //in the scaling->scaling_data pointer for later reuse.
  //
  //(What should we do if AZ_keep_info==1 and the
  //scaling->scaling_data pointer already contains data?
  //For now we'll assume that if that's the case, then the
  //scaling vector we're about to store is the same one that's
  //already there...)

  if (options[AZ_keep_info] == 1) {
    scaling->scaling_data = (void*)vec;
  }
  else {
    delete vec;
    scaling->scaling_data = 0;
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

  if (scaling_type == AZ_Jacobi) {
    int err = A->ExtractDiagonalCopy(*vec);
    if (err != 0) {
      delete vec; vec = 0;
      return(vec);
    }

    double* vec_vals = NULL;
    vec->ExtractView(&vec_vals);

    //now invert each entry of the diagonal vector
    for(int i=0; i<A->RowMatrixRowMap().NumMyElements(); ++i) {
      if (fabs(vec_vals[i]) > Epetra_MinDouble) {
        vec_vals[i] = 1.0/vec_vals[i];
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

