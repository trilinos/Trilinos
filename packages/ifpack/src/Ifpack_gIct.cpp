/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_gIct.h"
#include "Ifpack_gIct_Utils.h"
#include "Ifpack_Condest.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Util.h"

#ifdef HAVE_IFPACK_TEUCHOS
#include "Teuchos_ParameterList.hpp"
#endif

//==============================================================================
Ifpack_gIct::Ifpack_gIct(Epetra_RowMatrix* A) :
  A_(*A),
  Comm_(A->Comm()),
  IsInitialized_(false),
  IsComputed_(false),
  Condest_(-1.0),
  Athresh_(0.0),
  Rthresh_(1.0),
  Droptol_(0.0),
  Lfil_(0),
  Aict_(0),
  Lict_(0),
  Ldiag_(0),
  U_(0),
  D_(0)
{
#ifdef HAVE_IFPACK_TEUCHOS
  Teuchos::ParameterList List;
  SetParameters(List);
#endif

}
//==============================================================================
Ifpack_gIct::~Ifpack_gIct(){


  delete U_;
  delete D_; // Diagonal is stored separately.  We store the inverse.

  if (Lict_ != 0) {
    Ifpack_AIJMatrix * Lict = (Ifpack_AIJMatrix *) Lict_;
    delete [] Lict->ptr;
    delete [] Lict->col;
    delete [] Lict->val;
    delete Lict;
  }
  if (Aict_ != 0) {
    Ifpack_AIJMatrix * Aict = (Ifpack_AIJMatrix *) Aict_;
    delete Aict;
  }
  if (Ldiag_ != 0) delete [] Ldiag_;

  IsInitialized_ = false;
  IsComputed_ = false;
}

#ifdef HAVE_IFPACK_TEUCHOS
//==========================================================================
int Ifpack_gIct::SetParameters(Teuchos::ParameterList& List)
{

  Lfil_ = List.get("fact: level-of-fill",Lfil_);
  Athresh_ = List.get("fact: absolute threshold", Athresh_);
  Rthresh_ = List.get("fact: relative threshold", Rthresh_);
  Droptol_ = List.get("fact: drop tolerance", Droptol_);

  // set label
  sprintf(Label_, "graph-based ICT (fill=%d, drop=%f)",
	  Lfil_, Droptol_);
  return(0);
}
#endif

//==========================================================================
int Ifpack_gIct::Initialize()
{

  IsInitialized_ = false;
  // FIXME: construction of ILUK graph must be here
  
  IsInitialized_ = true;
  return(0);
}

//==========================================================================
int Ifpack_gIct::ComputeSetup()
{
  // (re)allocate memory for ICT factors
  if (U_) 
    delete U_;
  if (D_)
    delete D_;

  U_ = new Epetra_CrsMatrix(Copy, Matrix().RowMatrixRowMap(), 
			    Matrix().RowMatrixRowMap(), 0);
  D_ = new Epetra_Vector(Matrix().RowMatrixRowMap());

  if (U_ == 0 || D_ == 0)
    IFPACK_CHK_ERR(-1);

  int ierr = 0;
  int i, j;
  int * InI=0, * LI=0, * UI = 0;
  double * InV=0, * LV=0, * UV = 0;
  int NumIn, NumL, NumU;
  bool DiagFound;
  int NumNonzeroDiags = 0;

  // Get Maximun Row length
  int MaxNumEntries = Matrix().MaxNumEntries();

  InI = new int[MaxNumEntries]; // Allocate temp space
  UI  = new int[MaxNumEntries];
  InV = new double[MaxNumEntries];
  UV  = new double[MaxNumEntries];

  double *DV;
  ierr = D_->ExtractView(&DV); // Get view of diagonal

  // First we copy the user's matrix into diagonal vector and U, regardless of fill level

  int NumRows = Matrix().NumMyRows();

  for (i=0; i< NumRows; i++) {

    Matrix().ExtractMyRowCopy(i, MaxNumEntries, NumIn, InV, InI); // Get Values and Indices
    
    // Split into L and U (we don't assume that indices are ordered).
    NumL = 0; 
    NumU = 0; 
    DiagFound = false;
    
    for (j=0; j< NumIn; j++) {
      int k = InI[j];

      if (k==i) {
	DiagFound = true;
	// Store perturbed diagonal in Epetra_Vector D_
	DV[i] += Rthresh_ * InV[j] + EPETRA_SGN(InV[j]) * Athresh_; 
      }

      else if (k < 0) return(-1); // Out of range
      else if (i<k && k<NumRows) {
	UI[NumU] = k;
	UV[NumU] = InV[j];
	NumU++;
      }
    }
    
    // Check in things for this row of L and U

    if (DiagFound) NumNonzeroDiags++;
    if (NumU) U_->InsertMyValues(i, NumU, UV, UI);
    
  }

  delete [] UI;
  delete [] UV;
  delete [] InI;
  delete [] InV;

  U_->TransformToLocal(&Matrix().OperatorDomainMap(), 
		       &Matrix().OperatorRangeMap());

  int ierr1 = 0;
  if (NumNonzeroDiags<U_->NumMyRows()) ierr1 = 1;
  Matrix().Comm().MaxAll(&ierr1, &ierr, 1);
  IFPACK_CHK_ERR(ierr);

  return(0);
}

//==========================================================================
int Ifpack_gIct::Compute() {

  if (!IsInitialized()) 
    IFPACK_CHK_ERR(Initialize());

  IsComputed_ = false;

  // copy matrix into L and U.
  IFPACK_CHK_ERR(ComputeSetup());

  int i;

  int m, n, nz, Nrhs, ldrhs, ldlhs;
  int * ptr=0, * ind;
  double * val, * rhs, * lhs;

  int ierr = Epetra_Util_ExtractHbData(U_, 0, 0, m, n, nz, ptr, ind,
				       val, Nrhs, rhs, ldrhs, lhs, ldlhs);
  if (ierr < 0) 
    IFPACK_CHK_ERR(ierr);

  Ifpack_AIJMatrix * Aict;
  if (Aict_==0) {
    Aict = new Ifpack_AIJMatrix;
    Aict_ = (void *) Aict;
  }
  else Aict = (Ifpack_AIJMatrix *) Aict_;
  Ifpack_AIJMatrix * Lict;
  if (Lict_==0) {
    Lict = new Ifpack_AIJMatrix;
    Lict_ = (void *) Lict;
  }
  else Lict = (Ifpack_AIJMatrix *) Lict_;
  Aict->val = val;
  Aict->col = ind;
  Aict->ptr = ptr;
  double *DV;
  EPETRA_CHK_ERR(D_->ExtractView(&DV)); // Get view of diagonal
    
  crout_ict(m, Aict, DV, Droptol_, Lfil_, Lict, &Ldiag_);

  // Get rid of unnecessary data
  delete [] ptr;
  delete U_;
  delete D_;

  // Create Epetra View of L from crout_ict
  U_ = new Epetra_CrsMatrix(View, A_.RowMatrixRowMap(), A_.RowMatrixRowMap(),0);
  D_ = new Epetra_Vector(View, A_.RowMatrixRowMap(), Ldiag_);

  ptr = Lict->ptr;
  ind = Lict->col;
  val = Lict->val;
    
  for (i=0; i< m; i++) {
    int NumEntries = ptr[i+1]-ptr[i];
    int * Indices = ind+ptr[i];
    double * Values = val+ptr[i];
    U_->InsertMyValues(i, NumEntries, Values, Indices);
  }

  U_->TransformToLocal(&(A_.OperatorDomainMap()), &(A_.OperatorRangeMap()));
  
  D_->Reciprocal(*D_); // Put reciprocal of diagonal in this vector
  // Add up flops
 
  double current_flops = 2 * nz; // Just an estimate
  double total_flops = 0;
    
  A_.Comm().SumAll(&current_flops, &total_flops, 1); // Get total madds across all PEs

  // Now count the rest
  total_flops += (double) U_->NumGlobalNonzeros(); // Accounts for multiplier above
  total_flops += (double) D_->GlobalLength(); // Accounts for reciprocal of diagonal

  IsComputed_ = true;

  return(0);

}

//=============================================================================
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_gIct::ApplyInverse(const Epetra_MultiVector& X, 
			     Epetra_MultiVector& Y) const
{

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK_CHK_ERR(-1); // Return error: X and Y not the same size

  if (!IsComputed())
    IFPACK_CHK_ERR(-1); // compute preconditioner first

  bool Upper = true;
  bool UnitDiagonal = true;

  Epetra_MultiVector Y1(Y);

  U_->Solve(Upper, true, UnitDiagonal, X, Y);
  Y.Multiply(1.0, *D_, Y, 0.0); // y = D*y (D_ has inverse of diagonal)
  U_->Solve(Upper, false, UnitDiagonal, Y, Y); // Solve Uy = y
  return(0);

}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_gIct::Apply(const Epetra_MultiVector& X, 
		      Epetra_MultiVector& Y) const 
{

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK_CHK_ERR(-1); // Return error: X and Y not the same size

  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

  Epetra_MultiVector * X1 = (Epetra_MultiVector *) &X;
  Epetra_MultiVector * Y1 = (Epetra_MultiVector *) &Y;

  U_->Multiply(false, *X1, *Y1);
  Y1->Update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
  Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
  Epetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
  U_->Multiply(true, Y1temp, *Y1);
  Y1->Update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
  return(0);
}
//=============================================================================
double Ifpack_gIct::Condest(const Ifpack_CondestType CT, 
                            const int MaxIters, const double Tol,
                            Epetra_RowMatrix* Matrix)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix);

  return(Condest_);
}
