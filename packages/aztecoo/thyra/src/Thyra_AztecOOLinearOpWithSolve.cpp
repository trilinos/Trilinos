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

#ifndef __sun

#include "Thyra_AztecOOLinearOpWithSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Time.hpp"

namespace {

Teuchos::RefCountPtr<Teuchos::Time> overallSolveTimer, individualSolveTimer;

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
  ,const bool     outputEveryRhs
  )
  :fwdDefaultMaxIterations_(fwdDefaultMaxIterations)
  ,fwdDefaultTol_(fwdDefaultTol)
  ,adjDefaultMaxIterations_(adjDefaultMaxIterations)
  ,adjDefaultTol_(adjDefaultTol)
  ,outputEveryRhs_(outputEveryRhs)
  ,isExternalPrec_(false)
  ,allowInexactFwdSolve_(false)
  ,allowInexactAdjSolve_(false)
  ,aztecSolverScalar_(0.0)
{
  initializeTimers();
}

void AztecOOLinearOpWithSolve::initialize(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >                 &fwdOp
  ,const Teuchos::RefCountPtr<const PreconditionerBase<double> >          &prec
  ,const bool                                                             isExternalPrec
  ,const Teuchos::RefCountPtr<const LinearOpBase<double> >                &approxFwdOp
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
  isExternalPrec_ = isExternalPrec;
  prec_ = prec;
  approxFwdOp_ = approxFwdOp;
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

Teuchos::RefCountPtr<const PreconditionerBase<double> >
AztecOOLinearOpWithSolve::extract_prec()
{
  Teuchos::RefCountPtr<const PreconditionerBase<double> > _prec = prec_;
  prec_ = Teuchos::null;
  return _prec;
}

bool AztecOOLinearOpWithSolve::isExternalPrec() const
{
  return isExternalPrec_;
}

Teuchos::RefCountPtr<const LinearOpBase<double> >
AztecOOLinearOpWithSolve::extract_approxFwdOp()
{
  Teuchos::RefCountPtr<const LinearOpBase<double> > _approxFwdOp = approxFwdOp_;
  approxFwdOp_ = Teuchos::null;
  return _approxFwdOp;
}

void AztecOOLinearOpWithSolve::uninitialize(
  Teuchos::RefCountPtr<const LinearOpBase<double> >                 *fwdOp
  ,Teuchos::RefCountPtr<const PreconditionerBase<double> >          *prec
  ,bool                                                             *isExternalPrec
  ,Teuchos::RefCountPtr<const LinearOpBase<double> >                *approxFwdOp
  ,Teuchos::RefCountPtr<AztecOO>                                    *aztecFwdSolver
  ,bool                                                             *allowInexactFwdSolve
  ,Teuchos::RefCountPtr<AztecOO>                                    *aztecAdjSolver
  ,bool                                                             *allowInexactAdjSolve
  ,double                                                           *aztecSolverScalar
  )
{
  if(fwdOp) *fwdOp = fwdOp_;
  if(prec) *prec = prec_;
  if(isExternalPrec) *isExternalPrec = isExternalPrec_;
  if(approxFwdOp) *approxFwdOp = approxFwdOp_;
  if(aztecFwdSolver) *aztecFwdSolver = aztecFwdSolver_;
  if(allowInexactFwdSolve) *allowInexactFwdSolve = allowInexactFwdSolve_;
  if(aztecAdjSolver) *aztecAdjSolver = aztecAdjSolver_;
  if(allowInexactAdjSolve) *allowInexactAdjSolve = allowInexactAdjSolve_;
  if(aztecSolverScalar) *aztecSolverScalar = aztecSolverScalar_;

  fwdOp_ = Teuchos::null;
  prec_ = Teuchos::null;
  isExternalPrec_ = false; // Just to make unique
  approxFwdOp_ = Teuchos::null;
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
    oss << "{";
    oss << "fwdOp=\'"<<fwdOp_->description()<<"\'";
    oss << "}";
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

bool AztecOOLinearOpWithSolve::solveSupportsSolveMeasureType(ETransp M_trans, const SolveMeasureType& solveMeasureType) const
{
  if(real_trans(M_trans)==NOTRANS) {
    return (
      solveMeasureType.useDefault()
      ||
      ( solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS)
        &&
        allowInexactFwdSolve_
        )
      );
  }
  // TRANS
  return (
    aztecAdjSolver_.get()!=NULL
    &&
    (
      solveMeasureType.useDefault()
      ||
      ( solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS)
        &&
        allowInexactAdjSolve_
        )
      )
    );
}

// Overridden from SingleRhsLinearOpWithSolveBase

void AztecOOLinearOpWithSolve::solve(
  const ETransp                         M_trans
  ,const MultiVectorBase<double>        &B
  ,MultiVectorBase<double>              *X
  ,const int                            numBlocks
  ,const BlockSolveCriteria<double>     blockSolveCriteria[]
  ,SolveStatus<double>                  blockSolveStatus[]
  ) const
{
  using Teuchos::OSTab;
  typedef SolveCriteria<double>  SC;
  typedef SolveStatus<double>    SS;

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);
  Teuchos::TimeMonitor timeMonitor(*overallSolveTimer);

  Teuchos::RefCountPtr<Teuchos::FancyOStream>  out = this->getOStream();
  Teuchos::EVerbosityLevel                     verbLevel = this->getVerbLevel();
  OSTab tab = this->getOSTab();
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
    *out << "\nSolving block system using AztecOO ...\n\n";

  //
  // Validate input
  //
  TEST_FOR_EXCEPT(numBlocks > 1); // ToDo: Deal with multiple solve criteria later if needed
//#ifdef _DEBUG
  TEST_FOR_EXCEPT(!this->solveSupportsTrans(M_trans));
  TEST_FOR_EXCEPT( numBlocks && blockSolveCriteria && !this->solveSupportsSolveMeasureType(M_trans,blockSolveCriteria[0].solveCriteria.solveMeasureType) );
  TEST_FOR_EXCEPT(X==NULL);
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
  // Get the convergence criteria
  //
  double tol            = ( aztecOpTransp==NOTRANS ? fwdDefaultTol()           : adjDefaultTol()           );
  int    maxIterations  = ( aztecOpTransp==NOTRANS ? fwdDefaultMaxIterations() : adjDefaultMaxIterations() );
  bool   isDefaultSolveCriteria = true;
  if( numBlocks && blockSolveCriteria ) {
    if( blockSolveCriteria[0].solveCriteria.requestedTol != SC::unspecifiedTolerance() ) {
      tol = blockSolveCriteria[0].solveCriteria.requestedTol;
      isDefaultSolveCriteria = true;
    }
  }
  //
  // Get Epetra_MultiVector views of B and X
  //
  Teuchos::RefCountPtr<const Epetra_MultiVector>
    epetra_B = get_Epetra_MultiVector(opRangeMap,Teuchos::rcp(&B,false));
  Teuchos::RefCountPtr<Epetra_MultiVector>
    epetra_X = get_Epetra_MultiVector(opDomainMap,Teuchos::rcp(X,false));
  //
  // Use AztecOO to solve each RHS one at a time (which is all that I can do anyway)
  //
  int totalIterations = 0;
  SolveStatus<double> solveStatus;
  solveStatus.solveStatus = SOLVE_STATUS_CONVERGED;
  solveStatus.achievedTol = -1.0;
  const int m = epetra_B->NumVectors();
  for( int j = 0; j < m; ++j ) {
    Teuchos::TimeMonitor timeMonitor(*individualSolveTimer);
    //
    // Get Epetra_Vector views of B(:,j) and X(:,j)
    //
    const Epetra_Vector   *epetra_b_j = (*epetra_B)(j);
    Epetra_Vector         *epetra_x_j = (*epetra_X)(j);
    TEST_FOR_EXCEPT(!epetra_b_j);
    TEST_FOR_EXCEPT(!epetra_x_j);
    //
    // Set the RHS and LHS
    //
    aztecSolver->SetRHS( const_cast<Epetra_Vector*>(epetra_b_j) ); // Should be okay?
    aztecSolver->SetLHS( epetra_x_j );
    //
    // Solve the linear system
    //
    timer.start(true);
    aztecSolver->Iterate( maxIterations, tol ); // We ignore the returned status but get it below
    timer.stop();
    //
    // Set the return solve status
    //
    const int     iterations  = aztecSolver->NumIters();
    const double  achievedTol = aztecSolver->ScaledResidual();
    const double  *AZ_status  = aztecSolver->GetAztecStatus();
    std::ostringstream oss;
    bool converged = false;
    if(AZ_status[AZ_why]==AZ_normal)           { oss << "Aztec returned AZ_normal."; converged = true; }
    else if(AZ_status[AZ_why]==AZ_param)       oss << "Aztec returned AZ_param.";
    else if(AZ_status[AZ_why]==AZ_breakdown)   oss << "Aztec returned AZ_breakdown.";
    else if(AZ_status[AZ_why]==AZ_loss)        oss << "Aztec returned AZ_loss.";
    else if(AZ_status[AZ_why]==AZ_ill_cond)    oss << "Aztec returned AZ_ill_cond.";
    else if(AZ_status[AZ_why]==AZ_maxits)      oss << "Aztec returned AZ_maxits.";
    else                                       oss << "Aztec returned an unknown status?";
    oss << "  Iterations = " << iterations << ".";
    oss << "  Achieved Tolerance = " << achievedTol << ".";
    oss << "  Total time = " << timer.totalElapsedTime() << " sec.";
    if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE) && outputEveryRhs())
      *Teuchos::OSTab(out).getOStream() << "j="<<j<<": " << oss.str() << "\n";
    //
    totalIterations += iterations;
    solveStatus.achievedTol = TEUCHOS_MAX(solveStatus.achievedTol,achievedTol);
    // Note, achieveTol may actually be greater than tol due to ill conditioning and roundoff!
    solveStatus.message = oss.str();
    if( isDefaultSolveCriteria ) {
      switch(solveStatus.solveStatus) {
        case SOLVE_STATUS_UNKNOWN:
          // Leave overall unknown!
          break;
        case SOLVE_STATUS_CONVERGED:
          solveStatus.solveStatus = ( converged ? SOLVE_STATUS_CONVERGED : SOLVE_STATUS_UNCONVERGED );
          break;
        case SOLVE_STATUS_UNCONVERGED:
          // Leave overall unconverged!
          break;
        default:
          TEST_FOR_EXCEPT(true); // Should never get here!
      }
    }
  }
  //
  // Scale the solution
  //
  if(aztecSolverScalar_ != 1.0)
    epetra_X->Scale(1.0/aztecSolverScalar_);
  //
  // Release the Epetra_MultiVector views of X and B
  //
  epetra_X = Teuchos::null;
  epetra_B = Teuchos::null;
  //
  // Update the overall solve criteria
  //
  totalTimer.stop();
  if( numBlocks && blockSolveStatus ) {
    std::ostringstream oss;
    oss
      << "AztecOO solver "
      << ( solveStatus.solveStatus==SOLVE_STATUS_CONVERGED ? "converged" : "unconverged" )
      << " on m = "<<m<<" RHSs using " << totalIterations << " cumulative iterations"
      << " for an average of " << (totalIterations/m) << " iterations/RHS and"
      << " total CPU time of "<<totalTimer.totalElapsedTime()<<" sec.";
    blockSolveStatus[0].message     = oss.str();
    if(isDefaultSolveCriteria) {
      blockSolveStatus[0].solveStatus = SOLVE_STATUS_UNKNOWN;
      blockSolveStatus[0].achievedTol = SS::unknownTolerance();
    }
    else {
      blockSolveStatus[0].solveStatus = solveStatus.solveStatus;
      blockSolveStatus[0].achievedTol = solveStatus.achievedTol;
    }
  }
  //
  // Report the overall time
  //
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal solve time = "<<totalTimer.totalElapsedTime()<<" sec\n";
}

// private

void AztecOOLinearOpWithSolve::initializeTimers()
{
  if(!overallSolveTimer.get()) {
    overallSolveTimer    = Teuchos::TimeMonitor::getNewTimer("AztecOOLOWS");
    individualSolveTimer = Teuchos::TimeMonitor::getNewTimer("AztecOOLOWS:SingleSolve");
  }
}

}	// end namespace Thyra

#endif // __sun





