/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "AztecOO_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_BlockMap.h"

//==============================================================================
AztecOO_Operator::AztecOO_Operator(AztecOO * solver, int NumIters) 
  : solver_(solver),
    NumIters_(NumIters),
    Label_(0) {

  Label_ = "AztecOO Operator";
}
//==============================================================================
AztecOO_Operator::~AztecOO_Operator() {
}
//==============================================================================
int AztecOO_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  if (!X.Map().SameAs(DomainMap())) EPETRA_CHK_ERR(-1);
  if (!Y.Map().SameAs(RangeMap())) EPETRA_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);

  Epetra_MultiVector xtmp(X); // Make copy of X (needed in case X is scaled in solver or if X = Y

  Y.PutScalar(0.0); // Always start with Y = 0

  solver_->SetRHS(&xtmp); // Set RHS to the input X vector copy
  solver_->SetLHS(&Y);

  // Finally do iterations (set tolerance to zero to force all iterations to be done)
  int ierr = solver_->recursiveIterate(NumIters_, 0.0);

  if (ierr==AZ_maxits) ierr = 0; // We force maxits, don't report as an error

  return(ierr);
}
