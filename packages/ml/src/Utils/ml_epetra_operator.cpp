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

#include "ml_epetra_operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace Epetra_ML 
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

  return 0;
}

#else
int MultiLevelOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVe
ctor& Y , int iBlockSize) const {


  if (!X.Map().SameAs(OperatorDomainMap())) EPETRA_CHK_ERR(-1);
  if (!Y.Map().SameAs(OperatorRangeMap())) EPETRA_CHK_ERR(-2);
  if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);

  MPI_Pcontrol ( 8 , "entry" , 2 , 0 , 0 );

  Epetra_MultiVector xtmp(X); // Make copy of X (needed in case X is scaled
                              // in solver or if X = Y
  Y.PutScalar(0.0); // Always start with Y = 0

  // ML_iterate doesn't handle multivectors, so extract and iterate one at
  // a time on them.
  double **xvectors;
  double **yvectors;

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

      Epetra_MultiVector  cur_x ( View , xtmp , iOffset , X.NumVectors() % iBloc
kSize );
      Epetra_MultiVector  cur_y ( View , Y , iOffset , Y.NumVectors() % iBlockSi
ze );

      ML_Solve_MGV(solver_, xtmp , Y );

      }


  return 0;
}


//==============================================================================
int MultiLevelOperator::ApplyInverse_WKC(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

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

  return 0;
}
#endif
 
}

#endif //ifdef ML_WITH_EPETRA
