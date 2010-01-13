//@HEADER
// ***********************************************************************
// 
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
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

#include "Tifpack_OverlapSolveObject.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Flops.hpp"

//==============================================================================
Tifpack_OverlapSolveObject::Tifpack_OverlapSolveObject(char * Label, const Tpetra_Comm & Comm) 
  : Label_(Label),
    L_(0),
    UseLTrans_(false),
    D_(0),
    U_(0),
    UseUTrans_(false),
    UseTranspose_(false),
    Comm_(Comm),
    Condest_(-1.0),
    OverlapMode_(Zero)
{
}

//==============================================================================
Tifpack_OverlapSolveObject::Tifpack_OverlapSolveObject(const Tifpack_OverlapSolveObject & Source) 
  : Label_(Source.Label_),
    L_(Source.L_),
    UseLTrans_(Source.UseLTrans_),
    D_(Source.D_),
    U_(Source.U_),
    UseUTrans_(Source.UseUTrans_),
    UseTranspose_(Source.UseTranspose_),
    Comm_(Source.Comm_),
    Condest_(Source.Condest_),
    OverlapMode_(Source.OverlapMode_)
{
}
//==============================================================================
Tifpack_OverlapSolveObject::~Tifpack_OverlapSolveObject(){
}
//==============================================================================
//=============================================================================
int Tifpack_OverlapSolveObject::Solve(bool Trans, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
				Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const {
//
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//

  // First generate X and Y as needed for this function
  Tpetra_MultiVector * X1 = 0;
  Tpetra_MultiVector * Y1 = 0;
  EPETRA_CHK_ERR(SetupXY(Trans, X, Y, X1, Y1));

  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

  if (!Trans) {

    EPETRA_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, *X1, *Y1));
    EPETRA_CHK_ERR(Y1->Multiply(1.0, *D_, *Y1, 0.0)); // y = D*y (D_ has inverse of diagonal)
    EPETRA_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, *Y1, *Y1)); // Solve Uy = y
    if (L_->Exporter()!=0) {EPETRA_CHK_ERR(Y.Export(*Y1,*L_->Exporter(), OverlapMode_));} // Export computed Y values if needed
  }
  else {
    EPETRA_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, *X1, *Y1)); // Solve Uy = y
    EPETRA_CHK_ERR(Y1->Multiply(1.0, *D_, *Y1, 0.0)); // y = D*y (D_ has inverse of diagonal)
    EPETRA_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, *Y1, *Y1));
    if (U_->Importer()!=0) {EPETRA_CHK_ERR(Y.Export(*Y1,*U_->Importer(), OverlapMode_));} // Export computed Y values if needed
  } 

  return(0);
}
//=============================================================================
int Tifpack_OverlapSolveObject::Multiply(bool Trans, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
				Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const {
//
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//
    
  // First generate X and Y as needed for this function
  Tpetra_MultiVector * X1 = 0;
  Tpetra_MultiVector * Y1 = 0;
  EPETRA_CHK_ERR(SetupXY(Trans, X, Y, X1, Y1));

  if (!Trans) {
    EPETRA_CHK_ERR(U_->Multiply(Trans, *X1, *Y1)); // 
    EPETRA_CHK_ERR(Y1->Update(1.0, *X1, 1.0)); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    EPETRA_CHK_ERR(Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0)); // y = D*y (D_ has inverse of diagonal)
    Tpetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
    EPETRA_CHK_ERR(L_->Multiply(Trans, Y1temp, *Y1));
    EPETRA_CHK_ERR(Y1->Update(1.0, Y1temp, 1.0)); // (account for implicit unit diagonal)
    if (L_->Exporter()!=0) {EPETRA_CHK_ERR(Y.Export(*Y1,*L_->Exporter(), OverlapMode_));} // Export computed Y values if needed
  }
  else {

    EPETRA_CHK_ERR(L_->Multiply(Trans, *X1, *Y1));
    EPETRA_CHK_ERR(Y1->Update(1.0, *X1, 1.0)); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    EPETRA_CHK_ERR(Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0)); // y = D*y (D_ has inverse of diagonal)
    Tpetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
    EPETRA_CHK_ERR(U_->Multiply(Trans, Y1temp, *Y1));
    EPETRA_CHK_ERR(Y1->Update(1.0, Y1temp, 1.0)); // (account for implicit unit diagonal)
    if (L_->Exporter()!=0) {EPETRA_CHK_ERR(Y.Export(*Y1,*L_->Exporter(), OverlapMode_));}
  } 
  return(0);
} 
//=========================================================================
int Tifpack_OverlapSolveObject::SetupXY(bool Trans, 
				       const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xin, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Yin,
				       Tpetra_MultiVector * & Xout, Tpetra_MultiVector * & Yout) const {

  // Generate an X and Y suitable for performing Solve() and Multiply() methods
  
  if (Xin.NumVectors()!=Yin.NumVectors()) EPETRA_CHK_ERR(-1); // Return error: X and Y not the same size
  
  //cout << "Xin = " << Xin << endl;
  Xout = (Tpetra_MultiVector *) &Xin;
  Yout = (Tpetra_MultiVector *) &Yin;
  return(0);
}
//=============================================================================
int Tifpack_OverlapSolveObject::Condest(bool Trans, double & ConditionNumberEstimate) const {

  if (Condest_>=0.0) {
    ConditionNumberEstimate = Condest_;
    return(0);
  }
  // Create a vector with all values equal to one
  Tpetra_Vector Ones(U_->DomainMap());
  Tpetra_Vector OnesResult(L_->RangeMap());
  Ones.PutScalar(1.0);

  EPETRA_CHK_ERR(Solve(Trans, Ones, OnesResult)); // Compute the effect of the solve on the vector of ones
  EPETRA_CHK_ERR(OnesResult.Abs(OnesResult)); // Make all values non-negative
  EPETRA_CHK_ERR(OnesResult.MaxValue(&ConditionNumberEstimate)); // Get the maximum value across all processors
  Condest_ = ConditionNumberEstimate; // Save value for possible later calls
  return(0);
}
