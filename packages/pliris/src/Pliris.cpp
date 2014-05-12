/*
//@HEADER
// ************************************************************************
//
//               Pliris: Parallel Dense Solver Package
//                 Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER
*/

#include "Pliris.h"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "distribute.h"
#include "xlu_solve.h"
#include "permute.h"
#include "x_factor.h"
#include "x_solve.h"


//=============================================================================

Pliris::Pliris(Epetra_Vector * A,
		Epetra_MultiVector * X,
	        Epetra_MultiVector * B)  {



  inConstructor_ = true;

  SetMatrix(A);

  SetLHS(X);

  SetRHS(B);

  inConstructor_ = false;
}

Pliris::Pliris( ) {

}


//=============================================================================
Pliris::~Pliris(void) {


}

//=============================================================================
 int Pliris::GetDistribution( int*  nprocs_row,
                               int* number_of_unknowns,
			       int* nrhs,
                               int* my_rows,
                               int* my_cols,
                               int* my_first_row,
                               int* my_first_col,
                               int* my_rhs,
                               int* my_row,
                               int* my_col ) {



  // This function echoes the multiprocessor distribution of the matrix

  distmat_ ( nprocs_row,
             number_of_unknowns,
             nrhs,
             my_rows,
             my_cols,
             my_first_row,
             my_first_col,
             my_rhs,
             my_row,
             my_col);


  return(0);

}

//=============================================================================
int Pliris::FactorSolve( Epetra_Vector* A,
                         int my_rows,
                         int my_cols,
                         int* matrix_size,
                         int* num_procsr,
                         int* num_rhs,
                         double* secs){



      SetMatrix(A);


      dlusolve_(a_,
                matrix_size,
                num_procsr,
                (a_ + my_rows*my_cols),
                num_rhs,
                secs);

     return(0);

}

//=============================================================================
int Pliris::FactorSolve( Epetra_SerialDenseVector* A,
                         int my_rows,
                         int my_cols,
                         int* matrix_size,
                         int* num_procsr,
                         int* num_rhs,
                         double* secs){

     // Associate the matrix with a_

     SetMatrix(A);

     //  Factor and solve the matrix equation

     dlusolve_  (a_,
                 matrix_size,
                 num_procsr,
                 (a_ + my_rows*my_cols),
                 num_rhs,
                 secs);

     return(0);

}

//  Factor the matrix for later solve

//=============================================================================

int Pliris::Factor( Epetra_Vector* A,
                     int* matrix_size,
                     int* num_procsr,
                     int* permute,
                     double* secs){


     SetMatrix(A);

     // Factor the Matrix

     dfactor_(a_,
              matrix_size,
              num_procsr,
              permute,
              secs);


     // Permute the lower triangular matrix

     dpermute_(a_, permute);

     return(0);
}

// Solve the Matrix after previously factoring

//=============================================================================

int Pliris::Solve(int* permute,
	          int* num_rhs){

  dsolve_(a_,
	  permute,
	  b_,
	  num_rhs);

  return(0);

}


//=============================================================================
int Pliris::SetLHS(Epetra_MultiVector * X) {

  if (X == 0 && inConstructor_ == true) return(0);
  if (X == 0) EPETRA_CHK_ERR(-1);
  X_ = X;
  X_->ExtractView(&x_, &x_LDA_);
  return(0);
}
//=============================================================================
int Pliris::SetRHS(Epetra_MultiVector * B) {

  if (B == 0 && inConstructor_ == true) return(0);
  if (B == 0) EPETRA_CHK_ERR(-1);
  B_ = B;
  B_->ExtractView(&b_, &b_LDA_);

  return(0);
}
//=============================================================================
int Pliris::SetMatrix(Epetra_Vector * A) {

  if (A == 0 && inConstructor_ == true) return(0);
  if (A == 0) EPETRA_CHK_ERR(-1);
  A_ = A;
  A_->ExtractView(&a_);

 return(0);
}

//=============================================================================
int Pliris::SetMatrix(Epetra_SerialDenseVector * AA_) {

  a_ = AA_-> Values();

 return(0);

}
