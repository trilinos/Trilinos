/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_common.h"
#ifdef ML_WITH_EPETRA

#ifdef ML_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif
#include "ml_include.h"

#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace ML_Epetra 
{
  
//==============================================================================
// constructor -- it's presumed that the user has constructed the ML
// object somewhere else
MultiLevelOperator::MultiLevelOperator(ML *ml_handle, const Epetra_Comm &myComm,const Epetra_Map &dm, const Epetra_Map &rm)
  : solver_(ml_handle),
    Label_(0),
    Comm_(myComm),
    DomainMap_(dm),
    RangeMap_(rm) {
  ownership_ = false;
  Label_ = "Epetra ML::MultilLevelOperator";
}
//==============================================================================
MultiLevelOperator::~MultiLevelOperator() {
  if (ownership_ == true) {
    ML_Destroy(&solver_);
  }
}
//==============================================================================
#ifndef WKC
int MultiLevelOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  if (!X.Map().SameAs(OperatorDomainMap())) 
    ML_CHK_ERR(-1);
  if (!Y.Map().SameAs(OperatorRangeMap())) 
    ML_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) 
    ML_CHK_ERR(-3);

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

  return 0;
}

#else
//int MultiLevelOperator::ApplyInverse(const Epetra_MultiVector& X,
//                Epetra_MultiVector& Y , int iBlockSize) const {
int MultiLevelOperator::ApplyInverse(const Epetra_MultiVector& X,
                Epetra_MultiVector& Y ) const {


  if (!X.Map().SameAs(OperatorDomainMap())) 
    ML_CHK_ERR(-1);
  if (!Y.Map().SameAs(OperatorRangeMap())) 
    ML_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) 
    ML_CHK_ERR(-3);

#ifdef ML_MPI
  MPI_Pcontrol ( 8 , "entry" , 2 , 0 , 0 );
#endif

  Epetra_MultiVector xtmp(X); // Make copy of X (needed in case X is scaled
                              // in solver or if X = Y
  Y.PutScalar(0.0); // Always start with Y = 0

  // ML_iterate doesn't handle multivectors, so extract and iterate one at
  // a time on them.

   int iBlockSize = WKC;
   for ( int i = 0 ; i != (X.NumVectors()/iBlockSize) ; i++ )
      {
      int  iOffset = i * iBlockSize;

      Epetra_MultiVector  cur_x ( View , xtmp , iOffset , iBlockSize );
      Epetra_MultiVector  cur_y ( View , Y , iOffset , iBlockSize );
      ML_Solve_MGV(solver_, cur_x , cur_y );
      }

   if ( X.NumVectors() % iBlockSize )
      {
      int  iOffset = (X.NumVectors()/iBlockSize) * iBlockSize;

      Epetra_MultiVector  cur_x ( View , xtmp , iOffset , X.NumVectors() % iBlockSize );
      Epetra_MultiVector  cur_y ( View , Y , iOffset , Y.NumVectors() % iBlockSize );

      ML_Solve_MGV(solver_, xtmp , Y );

      }


  return 0;
}


//==============================================================================
int MultiLevelOperator::ApplyInverse_WKC(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  if (!X.Map().SameAs(OperatorDomainMap())) 
    ML_CHK_ERR(-1);
  if (!Y.Map().SameAs(OperatorRangeMap())) 
    ML_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) 
    ML_CHK_ERR(-3);

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

  return 0;
}
#endif
 
}

#endif //ifdef ML_WITH_EPETRA
