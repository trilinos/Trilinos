// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#ifndef TFLOP

#include "KundertOO.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsKundertSparse.h"


//=============================================================================
KundertOO::KundertOO(Epetra_RowMatrix * A, 
		 Epetra_MultiVector * X,
		 Epetra_MultiVector * B) {
  //  AllocAzArrays();
  SetKundertDefaults();

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  SetUserMatrix(A);
  
  SetLHS(X);
  SetRHS(B);
  inConstructor_ = false;
}

//=============================================================================
KundertOO::KundertOO() {
  //  AllocAzArrays();
  SetKundertDefaults();
}

//=============================================================================
KundertOO::~KundertOO(void) {
  //  DeleteMemory();
  //  DeleteAzArrays();
}

//=============================================================================
int KundertOO::SetUserMatrix(Epetra_RowMatrix * UserMatrix) {

  if (UserMatrix == 0 && inConstructor_ == true) return(0);
  if (UserMatrix == 0) EPETRA_CHK_ERR(-1);

  UserMatrix_ = UserMatrix;

  return(0);
}

//=============================================================================
int KundertOO::SetLHS(Epetra_MultiVector * X) {

  if (X == 0 && inConstructor_ == true) return(0);
  if (X == 0) EPETRA_CHK_ERR(-1);
  X_ = X;
  X_->ExtractView(&x_, &x_LDA_);
  return(0);
}
//=============================================================================
int KundertOO::SetRHS(Epetra_MultiVector * B) {

  if (B == 0 && inConstructor_ == true) return(0);
  if (B == 0) EPETRA_CHK_ERR(-1);
  B_ = B;
  B_->ExtractView(&b_, &b_LDA_);

  return(0);
}

int KundertOO::SetKundertDefaults() {

 UserOperator_ = 0;
 UserMatrix_ = 0;
 // PrecOperator_ = 0;
 // PrecMatrix_ = 0;
 X_ = 0;
 B_ = 0;
 
 x_LDA_ = 0;
 x_ = 0;
 b_LDA_ = 0;
 b_ = 0;
 Transpose_ = false ; 

 return(0);

}

//=============================================================================

int KundertOO::Solve() { 


  Epetra_RowMatrix   *matAA = (GetUserMatrix()) ; 
  Epetra_CrsMatrix   *matA = dynamic_cast<Epetra_CrsMatrix*>(matAA) ; 
  assert( matA != NULL ) ; 
  Epetra_MultiVector   *vecX = GetLHS() ; 
  Epetra_MultiVector   *vecB = GetRHS() ; 
  assert( vecX->NumVectors() == 1 ) ; 
  assert( vecB->NumVectors() == 1 ) ; 
  int iam =  matA->Comm().MyPID() ;
  assert ( ! GetTrans() ) ;



  Epetra_LinearProblem problem(matA, vecX, vecB);

  Epetra_CrsKundertSparse solver(&problem);

  solver.Solve();

  return(1) ; 
}
#endif
