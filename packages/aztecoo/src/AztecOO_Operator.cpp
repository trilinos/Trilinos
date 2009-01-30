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

#include "AztecOO_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"

//==============================================================================
AztecOO_Operator::AztecOO_Operator(AztecOO * solver, int numIters, double tol) 
  : solver_(solver),
    NumIters_(numIters),
    Tol_(tol),
    Label_(0) {

  Label_ = "AztecOO Operator";
}
//==============================================================================
AztecOO_Operator::~AztecOO_Operator() {
}
//==============================================================================
int AztecOO_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  if (!X.Map().SameAs(OperatorDomainMap())) EPETRA_CHK_ERR(-1);
  if (!Y.Map().SameAs(OperatorRangeMap())) EPETRA_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);

  Epetra_MultiVector xtmp(X); // Make copy of X (needed in case X is scaled in solver or if X = Y

  Y.PutScalar(0.0); // Always start with Y = 0

  solver_->SetRHS(&xtmp); // Set RHS to the input X vector copy
  solver_->SetLHS(&Y);

  // Finally do iterations (set tolerance to zero to force all iterations to be done)
  int ierr = solver_->recursiveIterate(NumIters_, Tol_);
  //int ierr = solver_->recursiveIterate(NumIters_, Tol_);

  if (ierr==1) ierr = 0; // We force maxits, don't report as an error

  return(ierr);
}
