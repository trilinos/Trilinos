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

#ifdef ML_WITH_EPETRA

#ifdef ML_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif
#include "ml_include.h"

#include "ml_epetra_operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

//==============================================================================
// constructor -- it's presumed that the user has constructed the ML
// object somewhere else
Epetra_ML_Operator::Epetra_ML_Operator(ML *ml_handle, const Epetra_Comm &myComm,const Epetra_Map &dm, const Epetra_Map &rm)
  : solver_(ml_handle),
    Label_(0),
    Comm_(myComm),
    DomainMap_(dm),
    RangeMap_(rm) {
  Label_ = "Epetra ML_Operator";
}
//==============================================================================
Epetra_ML_Operator::~Epetra_ML_Operator() {
}
//==============================================================================
int Epetra_ML_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  if (!X.Map().SameAs(OperatorDomainMap())) EPETRA_CHK_ERR(-1);
  if (!Y.Map().SameAs(OperatorRangeMap())) EPETRA_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);

  Epetra_MultiVector xtmp(X); // Make copy of X (needed in case X is scaled
                              // in solver or if X = Y
  Y.PutScalar(0.0); // Always start with Y = 0
  // ML_iterate doesn't handle multivectors, so extract and iterate one at
  // a time on them.
  double **xvectors;
  double **yvectors;
  int ierr = xtmp.ExtractView(&xvectors);
  ierr = Y.ExtractView(&yvectors);

  //note: solver_ is the ML handle
  for (int i=0; i < X.NumVectors(); i++)
  {
    switch(solver_->ML_scheme) {
      case(ML_MGFULLV):
        ML_Solve_MGFull(solver_,
                        xvectors[i],  //rhs
                        yvectors[i]); //solution
        break;
      case(ML_SAAMG): //Marian Brezina's solver
        ML_Solve_AMGV(solver_,
                      xvectors[i],  //rhs
                      yvectors[i]); //solution
        break;
      default:
        ML_Solve_MGV(solver_,
                     xvectors[i],  //rhs
                     yvectors[i]); //solution
    }
  }

}

#endif //ifdef ML_WITH_EPETRA
