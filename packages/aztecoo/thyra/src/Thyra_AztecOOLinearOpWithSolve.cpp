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

#include "Thyra_AztecOOLinearOpWithSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

namespace {

inline
Teuchos::ETransp convert( Thyra::ETransp trans_in )
{
	Teuchos::ETransp  trans_out;
	switch(trans_in) {
		case Thyra::NOTRANS:
			trans_out = Teuchos::NO_TRANS;
			break;
		case Thyra::TRANS:
			trans_out = Teuchos::TRANS;
			break;
		default:
			TEST_FOR_EXCEPT(true); // Should never get here!
	}
	return trans_out;
}

} // namespace

namespace Thyra {

// Constructors/initializers/accessors

AztecOOLinearOpWithSolve::AztecOOLinearOpWithSolve(
  const int       fwdDefaultMaxIterations
  ,const double   fwdDefaultTol
  ,const int      adjDefaultMaxIterations
  ,const double   adjDefaultTol
  )
  :fwdDefaultMaxIterations_(fwdDefaultMaxIterations)
  ,fwdDefaultTol_(fwdDefaultTol)
  ,adjDefaultMaxIterations_(adjDefaultMaxIterations)
  ,adjDefaultTol_(adjDefaultTol)
{}

void AztecOOLinearOpWithSolve::initialize(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >                 &fwdOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<double> >                &precOp
  ,const EPreconditionerInputType                                         precOpType
  ,const Teuchos::RefCountPtr<AztecOO>                                    &aztecFwdSolver
  ,const bool                                                             allowInexactFwdSolve
  ,const Teuchos::RefCountPtr<AztecOO>                                    &aztecAdjSolver
  ,const bool                                                             allowInexactAdjSolve
  ,const double                                                           aztecSolverScalar
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(fwdOp.get()==NULL);
  TEST_FOR_EXCEPT(aztecFwdSolver.get()==NULL);
#endif
  fwdOp_ = fwdOp;
  precOp_ = precOp;
  precOpType_ = precOpType;
  aztecFwdSolver_ = aztecFwdSolver;
  allowInexactFwdSolve_ = allowInexactFwdSolve;
  aztecAdjSolver_ = aztecAdjSolver;
  allowInexactAdjSolve_ = allowInexactAdjSolve;
  aztecSolverScalar_ = aztecSolverScalar;
}

Teuchos::RefCountPtr<const LinearOpBase<double> >
AztecOOLinearOpWithSolve::extract_fwdOp()
{
  Teuchos::RefCountPtr<const LinearOpBase<double> > _fwdOp = fwdOp_;
  fwdOp_ = Teuchos::null;
  return _fwdOp;
}

Teuchos::RefCountPtr<const LinearOpBase<double> >
AztecOOLinearOpWithSolve::extract_precOp()
{
  Teuchos::RefCountPtr<const LinearOpBase<double> > _precOp = precOp_;
  precOp_ = Teuchos::null;
  return _precOp;
}

EPreconditionerInputType
AztecOOLinearOpWithSolve::extract_precOpType()
{
  return precOpType_;
}

void AztecOOLinearOpWithSolve::uninitialize(
  Teuchos::RefCountPtr<const LinearOpBase<double> >                 *fwdOp
  ,Teuchos::RefCountPtr<const LinearOpBase<double> >                *precOp
  ,EPreconditionerInputType                                         *precOpType
  ,Teuchos::RefCountPtr<AztecOO>                                    *aztecFwdSolver
  ,bool                                                             *allowInexactFwdSolve
  ,Teuchos::RefCountPtr<AztecOO>                                    *aztecAdjSolver
  ,bool                                                             *allowInexactAdjSolve
  ,double                                                           *aztecSolverScalar
  )
{
  if(fwdOp) *fwdOp = fwdOp_;
  if(precOp) *precOp = precOp_;
  if(precOpType) *precOpType = precOpType_;
  if(aztecFwdSolver) *aztecFwdSolver = aztecFwdSolver_;
  if(allowInexactFwdSolve) *allowInexactFwdSolve = allowInexactFwdSolve_;
  if(aztecAdjSolver) *aztecAdjSolver = aztecAdjSolver_;
  if(allowInexactAdjSolve) *allowInexactAdjSolve = allowInexactAdjSolve_;
  if(aztecSolverScalar) *aztecSolverScalar = aztecSolverScalar_;

  fwdOp_ = Teuchos::null;
  aztecFwdSolver_ = Teuchos::null;
  allowInexactFwdSolve_ = false;
  aztecAdjSolver_ = Teuchos::null;
  allowInexactAdjSolve_ = false;
  aztecSolverScalar_ = 0.0;
}

// Overridden from LinearOpBase

Teuchos::RefCountPtr< const VectorSpaceBase<double> >
AztecOOLinearOpWithSolve::range() const
{
  return ( fwdOp_.get() ? fwdOp_->range() : Teuchos::null );
}

Teuchos::RefCountPtr< const VectorSpaceBase<double> >
AztecOOLinearOpWithSolve::domain() const
{
  return  ( fwdOp_.get() ? fwdOp_->domain() : Teuchos::null );
}

Teuchos::RefCountPtr<const LinearOpBase<double> >
AztecOOLinearOpWithSolve::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

// Overridden from Teuchos::Describable

std::string AztecOOLinearOpWithSolve::description() const
{
  std::ostringstream oss;
  oss << "Thyra::AztecOOLinearOpWithSolve";
  if(fwdOp_.get()) {
    oss << "(";
    oss << "fwdOp=\'"<<fwdOp_->description()<<"\'";
    oss << ")";
  }
  return oss.str();
}

// ToDo: Add more detailed describe() function override to show all of the good stuff!

// protected

// Overridden from SingleScalarLinearOpBase

bool AztecOOLinearOpWithSolve::opSupported(ETransp M_trans) const
{
  return ::Thyra::opSupported(*fwdOp_,M_trans);
}

// Overridden from SingleRhsLinearOpBase

void AztecOOLinearOpWithSolve::apply(
  const ETransp                M_trans
  ,const VectorBase<double>    &x
  ,VectorBase<double>          *y
  ,const double                alpha
  ,const double                beta
  ) const
{
  Thyra::apply( *fwdOp_, M_trans, x, y, alpha, beta );
}

// Overridden from SingleScalarLinearOpWithSolveBase

bool AztecOOLinearOpWithSolve::solveSupportsTrans(ETransp M_trans) const
{
  if(real_trans(M_trans)==NOTRANS) return true;
  return (aztecAdjSolver_.get()!=NULL);
}

bool AztecOOLinearOpWithSolve::solveSupportsSolveTolType(ETransp M_trans, ESolveTolType solveTolType) const
{
  if(real_trans(M_trans)==NOTRANS) {
    return ( solveTolType==SOLVE_TOL_DEFAULT || (solveTolType==SOLVE_TOL_REL_RESIDUAL_NORM && allowInexactFwdSolve_) );
  }
  // TRANS
  return ( aztecAdjSolver_.get()!=NULL && ( solveTolType==SOLVE_TOL_DEFAULT || (solveTolType==SOLVE_TOL_REL_RESIDUAL_NORM && allowInexactFwdSolve_) ) );
}

int AztecOOLinearOpWithSolve::defaultSolveMaxIterations(ETransp M_trans, ESolveTolType solveTolType) const
{
  return ( real_trans(M_trans) == NOTRANS ? fwdDefaultMaxIterations() : adjDefaultMaxIterations() );
}

// Overridden from SingleRhsLinearOpWithSolveBase

SolveStatus<double>
AztecOOLinearOpWithSolve::solve(
  const ETransp                         M_trans
  ,const VectorBase<double>             &b
  ,VectorBase<double>                   *x
  ,const SolveCriteria<double>          *solveCriteria
  ) const
{
  typedef SolveCriteria<double>  SC;
  typedef SolveStatus<double>    SS;
  //
  // Validate input
  //
//#ifdef _DEBUG
  TEST_FOR_EXCEPT(!this->solveSupportsTrans(M_trans));
  TEST_FOR_EXCEPT( solveCriteria && !this->solveSupportsSolveTolType(M_trans,solveCriteria->solveTolType) );
  TEST_FOR_EXCEPT(x==NULL);
//#endif
  //
  // Get the transpose argument
  //
  const ETransp aztecOpTransp = real_trans(M_trans);
  //
  // Get the solver, operator, and preconditioner that we will use
  //
  Teuchos::RefCountPtr<AztecOO>
    aztecSolver = ( aztecOpTransp == NOTRANS ? aztecFwdSolver_  : aztecAdjSolver_ );
  const Epetra_Operator
    *aztecOp = aztecSolver->GetUserOperator();
  //
  // Get the op(...) range and domain maps
  //
  const Epetra_Map
    &opRangeMap  = aztecOp->OperatorRangeMap(),
    &opDomainMap = aztecOp->OperatorDomainMap();
  //
  // Get Epetra_Vector views of b and x
  //
  Teuchos::RefCountPtr<const Epetra_Vector>
    epetra_b = get_Epetra_Vector(opRangeMap,Teuchos::rcp(&b,false));
  Teuchos::RefCountPtr<Epetra_Vector>
    epetra_x = get_Epetra_Vector(opDomainMap,Teuchos::rcp(x,false));
  //
  // Set the RHS and LHS
  //
  aztecSolver->SetRHS( const_cast<Epetra_Vector*>(&*epetra_b) ); // Should be okay?
  aztecSolver->SetLHS( &*epetra_x );
  //
  // Get the convergence criteria
  //
  double tol            = ( aztecOpTransp==NOTRANS ? fwdDefaultTol()           : adjDefaultTol()           );
  int    maxIterations  = ( aztecOpTransp==NOTRANS ? fwdDefaultMaxIterations() : adjDefaultMaxIterations() );
  if( solveCriteria ) {
    if( solveCriteria->requestedTol != SC::unspecifiedTolerance() )
      tol = solveCriteria->requestedTol;
    if( solveCriteria->maxIterations != SC::unspecifiedMaxIterations() )
      maxIterations = solveCriteria->maxIterations;
  }
  //
  // Solve the linear system
  //
	aztecSolver->Iterate( maxIterations, tol ); // We ignore the returned status (see below)
	//
	// Scale the solution
	//
  if(aztecSolverScalar_ != 1.0)
    epetra_x->Scale(1.0/aztecSolverScalar_);
  //
  // Release the Epetra_Vector views of x and b
  //
  epetra_x = Teuchos::null;
  epetra_b = Teuchos::null;
  //
  // Set the return solve status
  //
	const int     iterations  = aztecSolver->NumIters();
	const double  achievedTol = aztecSolver->ScaledResidual();
	const double *AZ_status   = aztecSolver->GetAztecStatus();
  std::ostringstream oss;
  bool converged = false;
  if(AZ_status[AZ_why]==AZ_normal)           { oss << "Aztec returned AZ_normal"; converged = true; }
  else if(AZ_status[AZ_why]==AZ_param)       oss << "Aztec returned AZ_param";
  else if(AZ_status[AZ_why]==AZ_breakdown)   oss << "Aztec returned AZ_breakdown";
  else if(AZ_status[AZ_why]==AZ_loss)        oss << "Aztec returned AZ_loss";
  else if(AZ_status[AZ_why]==AZ_ill_cond)    oss << "Aztec returned AZ_ill_cond";
  else if(AZ_status[AZ_why]==AZ_maxits)      oss << "Aztec returned AZ_maxits";
  else                                       oss << "Aztec returned an unknown status?";
  SolveStatus<double> solveStatus;
  solveStatus.solveStatus = SOLVE_STATUS_UNKNOWN;
  solveStatus.achievedTol = achievedTol;
  // Note, achieveTol may actually be greater than tol due to ill conditioning and roundoff!
  solveStatus.iterations = iterations;
  solveStatus.message = oss.str();
  if( solveCriteria && solveCriteria->solveTolType == SOLVE_TOL_REL_RESIDUAL_NORM ) {
    // This is for no left preconditioning and no left scaling only!
    solveStatus.solveStatus = ( converged ? SOLVE_STATUS_CONVERGED : SOLVE_STATUS_UNCONVERGED );
  }
  return solveStatus;
}

}	// end namespace Thyra
