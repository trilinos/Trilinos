/*
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
*/

#include "Thyra_AmesosLinearOpWithSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Epetra_MultiVector.h"

namespace Thyra {

// Constructors/initializers/accessors

AmesosLinearOpWithSolve::AmesosLinearOpWithSolve()
{}

AmesosLinearOpWithSolve::AmesosLinearOpWithSolve(
  const Teuchos::RefCountPtr<const EpetraLinearOpBase>    &epetraFwdOp
  ,const Teuchos::RefCountPtr<Epetra_LinearProblem>       &epetraLP
  ,const Teuchos::RefCountPtr<Amesos_BaseSolver>          &amesosSolver
  )
{
  this->initialize(epetraFwdOp,epetraLP,amesosSolver);
}

void AmesosLinearOpWithSolve::initialize(
  const Teuchos::RefCountPtr<const EpetraLinearOpBase>    &epetraFwdOp
  ,const Teuchos::RefCountPtr<Epetra_LinearProblem>       &epetraLP
  ,const Teuchos::RefCountPtr<Amesos_BaseSolver>          &amesosSolver
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(epetraFwdOp.get()==NULL);
  TEST_FOR_EXCEPT(epetraLP.get()==NULL);
  TEST_FOR_EXCEPT(amesosSolver.get()==NULL);
  TEST_FOR_EXCEPT(epetraLP->GetOperator()!=epetraFwdOp->epetra_op().get());
  TEST_FOR_EXCEPT(epetraLP->GetLHS()!=NULL);
  TEST_FOR_EXCEPT(epetraLP->GetRHS()!=NULL);
#endif
  epetraFwdOp_  = epetraFwdOp;
  epetraLP_     = epetraLP;
  amesosSolver_ = amesosSolver;
}

Teuchos::RefCountPtr<const EpetraLinearOpBase>
AmesosLinearOpWithSolve::extract_epetraFwdOp()
{
  Teuchos::RefCountPtr<const EpetraLinearOpBase> _epetraFwdOp = epetraFwdOp_;
  epetraFwdOp_ = Teuchos::null;
  return _epetraFwdOp;
}

void AmesosLinearOpWithSolve::reset_epetraFwdOp( const Teuchos::RefCountPtr<const EpetraLinearOpBase> &epetraFwdOp )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(get_epetraLP()->GetOperator()!=epetraFwdOp->epetra_op().get());
#endif
  epetraFwdOp_ = epetraFwdOp;
}

void AmesosLinearOpWithSolve::uninitialize(
  Teuchos::RefCountPtr<const EpetraLinearOpBase>    *epetraFwdOp
  ,Teuchos::RefCountPtr<Epetra_LinearProblem>       *epetraLP
  ,Teuchos::RefCountPtr<Amesos_BaseSolver>          *amesosSolver
  )
{

  if(epetraFwdOp)  *epetraFwdOp  = epetraFwdOp_;
  if(epetraLP)     *epetraLP     = epetraLP_;
  if(amesosSolver) *amesosSolver = amesosSolver_;

  epetraFwdOp_  = Teuchos::null;
  epetraLP_     = Teuchos::null;
  amesosSolver_ = Teuchos::null;

}

// Overridden from LinearOpBase

Teuchos::RefCountPtr< const VectorSpaceBase<double> >
AmesosLinearOpWithSolve::range() const
{
  return epetraFwdOp_->range();
}

Teuchos::RefCountPtr< const VectorSpaceBase<double> >
AmesosLinearOpWithSolve::domain() const
{
  return epetraFwdOp_->domain();
}

Teuchos::RefCountPtr<const LinearOpBase<double> >
AmesosLinearOpWithSolve::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

// Overridden from Teuchos::Describable

std::string AmesosLinearOpWithSolve::description() const
{
  std::ostringstream oss;
  oss << "Thyra::AmesosLinearOpWithSolve";
  if(amesosSolver_.get()) {
    oss << "(epetraFwdOp=\'"<<typeid(*epetraFwdOp_->epetra_op()).name()<<"\'"
        << ",amesosSolver=\'"<<typeid(*amesosSolver_).name()<<"\')";
  }
  return oss.str();
}

// protected

// Overridden from SingleScalarLinearOpBase

bool AmesosLinearOpWithSolve::opSupported(ETransp M_trans) const
{
  return true; // ToDo: Determine if the Epetra_Operator supports adjoints or not!
}

void AmesosLinearOpWithSolve::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<double>    &X
  ,MultiVectorBase<double>          *Y
  ,const double                     alpha
  ,const double                     beta
  ) const
{
  Thyra::apply( *epetraFwdOp_, M_trans, X, Y, alpha, beta );
}

// Overridden from SingleScalarLinearOpWithSolveBase

bool AmesosLinearOpWithSolve::solveSupportsTrans(ETransp M_trans) const
{
  return true; // ToDo: Determine if the solver supports adjoints or not!
}

bool AmesosLinearOpWithSolve::solveSupportsSolveTolType(ETransp M_trans, ESolveTolType solveTolType) const
{
  return true; // I am a direct solver so I should be able to do it all!
}

// Overridden from SingleRhsLinearOpWithSolveBase

void AmesosLinearOpWithSolve::solve(
  const ETransp                              M_trans
  ,const MultiVectorBase<double>             &B
  ,MultiVectorBase<double>                   *X
  ,const int                                 numBlocks
  ,const BlockSolveCriteria<double>          blockSolveCriteria[]
  ,SolveStatus<double>                       blockSolveStatus[]
  ) const
{
  typedef SolveCriteria<double>  SC;
  typedef SolveStatus<double>    SS;
#ifdef _DEBUG
  TEST_FOR_EXCEPT(X==NULL);
  TEST_FOR_EXCEPT(blockSolveCriteria==NULL && blockSolveStatus!=NULL);
#endif
  //
  // Get the op(...) range and domain maps
  //
  const Epetra_Map
    &opRangeMap  = ( real_trans(M_trans) == NOTRANS
                     ? epetraFwdOp_->epetra_op()->OperatorRangeMap()
                     : epetraFwdOp_->epetra_op()->OperatorDomainMap() ),
    &opDomainMap = ( real_trans(M_trans) == NOTRANS
                     ? epetraFwdOp_->epetra_op()->OperatorDomainMap()
                     : epetraFwdOp_->epetra_op()->OperatorRangeMap() );
  //
  // Get Epetra_MultiVector views of B and X
  //
  Teuchos::RefCountPtr<const Epetra_MultiVector>
    epetra_B = get_Epetra_MultiVector(opRangeMap,Teuchos::rcp(&B,false));
  Teuchos::RefCountPtr<Epetra_MultiVector>
    epetra_X = get_Epetra_MultiVector(opDomainMap,Teuchos::rcp(X,false));
  //
  // Set B and X in the linear problem
  //
  epetraLP_->SetLHS(&*epetra_X);
  epetraLP_->SetRHS(const_cast<Epetra_MultiVector*>(&*epetra_B)); // Should be okay but cross your fingers!
  //
  // Solve the linear system
  //
  const bool oldUseTranspose = amesosSolver_->UseTranspose();
  amesosSolver_->SetUseTranspose(real_trans(M_trans)==TRANS);
  TEST_FOR_EXCEPTION(
    0!=amesosSolver_->Solve(), CatastrophicSolveFailure
    ,"Error, the Amesos solver of type \'"<<typeid(*amesosSolver_).name()<<"\' could not perform the solve!"
    );
  amesosSolver_->SetUseTranspose(oldUseTranspose);
  //
  // Unset B and X in the linear problem
  //
  epetraLP_->SetLHS(NULL);
  epetraLP_->SetRHS(NULL);
  //
  // Set the solve status if requested
  //
  if(numBlocks && blockSolveStatus) {
    for( int i = 0; i < numBlocks; ++i ) {
      blockSolveStatus[i].solveStatus
        = (blockSolveCriteria[i].solveCriteria.requestedTol!=SC::unspecifiedTolerance()
           ? SOLVE_STATUS_CONVERGED : SOLVE_STATUS_UNKNOWN );
      blockSolveStatus[i].achievedTol = SS::unknownTolerance();
      blockSolveStatus[i].iterations = 1;
    }
  }
}

}	// end namespace Thyra
