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
#include "Thyra_LinearOpWithSolveHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_implicit_cast.hpp"


namespace {


inline
Teuchos::ETransp convert( Thyra::EOpTransp trans_in )
{
  Teuchos::ETransp trans_out;
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


// This class sets some solve instance specific state and then sets it back to
// the default state on destruction.  But using the destructor to unset the
// state we can be sure that the state is rest correctly even if an exception
// is thrown.
class SetAztecSolveState {
public:
  SetAztecSolveState( 
    const Teuchos::RCP<AztecOO> &aztecSolver,
    const Teuchos::RCP<Teuchos::FancyOStream> &fancyOStream,
    const Teuchos::EVerbosityLevel verbLevel,
    const Thyra::SolveMeasureType &solveMeasureType
    );
  ~SetAztecSolveState();
private:
  Teuchos::RCP<AztecOO> aztecSolver_;
  Teuchos::RCP<Teuchos::FancyOStream> fancyOStream_;
  Teuchos::EVerbosityLevel verbLevel_;
  int outputFrequency_;
  int convergenceTest_;
  SetAztecSolveState(); // Not defined and not to be called!
};


SetAztecSolveState::SetAztecSolveState( 
  const Teuchos::RCP<AztecOO> &aztecSolver,
  const Teuchos::RCP<Teuchos::FancyOStream> &fancyOStream,
  const Teuchos::EVerbosityLevel verbLevel,
  const Thyra::SolveMeasureType &solveMeasureType
  )
  :aztecSolver_(aztecSolver.assert_not_null())
{
  
  // Output state
  verbLevel_ = verbLevel;
  if ( Teuchos::VERB_NONE != verbLevel_ ) {
    if (!is_null(fancyOStream)) {
      // AztecOO puts in two tabs before it prints anything.  Therefore,
      // there is not much that we can do to improve the layout of the
      // indentation so just leave it!
      fancyOStream_= Teuchos::tab(
        fancyOStream.create_weak(),
        0, // Don't indent since AztecOO already puts in two tabs (not spaces!)
        Teuchos::implicit_cast<std::string>("AZTECOO")
        );
      aztecSolver_->SetOutputStream(*fancyOStream_);
      aztecSolver_->SetErrorStream(*fancyOStream_);
      // Note, above we can not save the current output and error streams
      // since AztecOO does not define functions to get them.  In the
      // future, AztecOO should define these functions if we are to avoid
      // treading on each others print statements.  However, since the
      // AztecOO object is most likely owned by these Thyra wrappers, this
      // should not be a problem.
    }
  }
  else {
    outputFrequency_ = aztecSolver_->GetAllAztecOptions()[AZ_output];
    aztecSolver_->SetAztecOption(AZ_output,0);
  }

  // Convergence test
  convergenceTest_ = aztecSolver_->GetAztecOption(AZ_conv);
  if (solveMeasureType.useDefault())
  {
    // Just use the default solve measure type already set!
  }
  else if (
    solveMeasureType(
      Thyra::SOLVE_MEASURE_NORM_RESIDUAL,
      Thyra::SOLVE_MEASURE_NORM_RHS
      )
    )
  {
    aztecSolver_->SetAztecOption(AZ_conv,AZ_rhs);
  }
  else if (
    solveMeasureType(
      Thyra::SOLVE_MEASURE_NORM_RESIDUAL,
      Thyra::SOLVE_MEASURE_NORM_INIT_RESIDUAL
      )
    )
  {
    aztecSolver_->SetAztecOption(AZ_conv,AZ_r0);
  }
  else {
    TEST_FOR_EXCEPT("Invalid solve measure type, you should never get here!");
  }

}


SetAztecSolveState::~SetAztecSolveState()
{
      
  // Output state
  if ( Teuchos::VERB_NONE != verbLevel_ ) {
    if (!is_null(fancyOStream_)) {
      aztecSolver_->SetOutputStream(std::cout);
      aztecSolver_->SetErrorStream(std::cerr);
      *fancyOStream_ << "\n";
    }
  }
  else {
    aztecSolver_->SetAztecOption(AZ_output,outputFrequency_);
  }

  // Convergence test
  aztecSolver_->SetAztecOption(AZ_conv,convergenceTest_);
      
}


} // namespace


namespace Thyra {


// Constructors/initializers/accessors


AztecOOLinearOpWithSolve::AztecOOLinearOpWithSolve(
  const int fwdDefaultMaxIterations
  ,const double fwdDefaultTol
  ,const int adjDefaultMaxIterations
  ,const double adjDefaultTol
  ,const bool outputEveryRhs
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
{}


void AztecOOLinearOpWithSolve::initialize(
  const RCP<const LinearOpBase<double> > &fwdOp
  ,const RCP<const LinearOpSourceBase<double> > &fwdOpSrc
  ,const RCP<const PreconditionerBase<double> > &prec
  ,const bool isExternalPrec
  ,const RCP<const LinearOpSourceBase<double> > &approxFwdOpSrc
  ,const RCP<AztecOO> &aztecFwdSolver
  ,const bool allowInexactFwdSolve
  ,const RCP<AztecOO> &aztecAdjSolver
  ,const bool allowInexactAdjSolve
  ,const double aztecSolverScalar
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(fwdOp.get()==NULL);
  TEST_FOR_EXCEPT(fwdOpSrc.get()==NULL);
  TEST_FOR_EXCEPT(aztecFwdSolver.get()==NULL);
#endif
  fwdOp_ = fwdOp;
  fwdOpSrc_ = fwdOpSrc;
  isExternalPrec_ = isExternalPrec;
  prec_ = prec;
  approxFwdOpSrc_ = approxFwdOpSrc;
  aztecFwdSolver_ = aztecFwdSolver;
  allowInexactFwdSolve_ = allowInexactFwdSolve;
  aztecAdjSolver_ = aztecAdjSolver;
  allowInexactAdjSolve_ = allowInexactAdjSolve;
  aztecSolverScalar_ = aztecSolverScalar;
  const std::string fwdOpLabel = fwdOp_->getObjectLabel();
  if (fwdOpLabel.length())
    this->setObjectLabel( "lows("+fwdOpLabel+")" );
}


RCP<const LinearOpSourceBase<double> >
AztecOOLinearOpWithSolve::extract_fwdOpSrc()
{
  RCP<const LinearOpSourceBase<double> >
    _fwdOpSrc = fwdOpSrc_;
  fwdOpSrc_ = Teuchos::null;
  return _fwdOpSrc;
}


RCP<const PreconditionerBase<double> >
AztecOOLinearOpWithSolve::extract_prec()
{
  RCP<const PreconditionerBase<double> >
    _prec = prec_;
  prec_ = Teuchos::null;
  return _prec;
}


bool AztecOOLinearOpWithSolve::isExternalPrec() const
{
  return isExternalPrec_;
}


RCP<const LinearOpSourceBase<double> >
AztecOOLinearOpWithSolve::extract_approxFwdOpSrc()
{
  RCP<const LinearOpSourceBase<double> >
    _approxFwdOpSrc = approxFwdOpSrc_;
  approxFwdOpSrc_ = Teuchos::null;
  return _approxFwdOpSrc;
}


void AztecOOLinearOpWithSolve::uninitialize(
  RCP<const LinearOpBase<double> > *fwdOp,
  RCP<const LinearOpSourceBase<double> > *fwdOpSrc,
  RCP<const PreconditionerBase<double> > *prec,
  bool *isExternalPrec,
  RCP<const LinearOpSourceBase<double> > *approxFwdOpSrc,
  RCP<AztecOO> *aztecFwdSolver,
  bool *allowInexactFwdSolve,
  RCP<AztecOO> *aztecAdjSolver,
  bool *allowInexactAdjSolve,
  double *aztecSolverScalar
  )
{
  if (fwdOp) *fwdOp = fwdOp_;
  if (fwdOpSrc) *fwdOpSrc = fwdOpSrc_;
  if (prec) *prec = prec_;
  if (isExternalPrec) *isExternalPrec = isExternalPrec_;
  if (approxFwdOpSrc) *approxFwdOpSrc = approxFwdOpSrc_;
  if (aztecFwdSolver) *aztecFwdSolver = aztecFwdSolver_;
  if (allowInexactFwdSolve) *allowInexactFwdSolve = allowInexactFwdSolve_;
  if (aztecAdjSolver) *aztecAdjSolver = aztecAdjSolver_;
  if (allowInexactAdjSolve) *allowInexactAdjSolve = allowInexactAdjSolve_;
  if (aztecSolverScalar) *aztecSolverScalar = aztecSolverScalar_;

  fwdOp_ = Teuchos::null;
  fwdOpSrc_ = Teuchos::null;
  prec_ = Teuchos::null;
  isExternalPrec_ = false; // Just to make unique
  approxFwdOpSrc_ = Teuchos::null;
  aztecFwdSolver_ = Teuchos::null;
  allowInexactFwdSolve_ = false;
  aztecAdjSolver_ = Teuchos::null;
  allowInexactAdjSolve_ = false;
  aztecSolverScalar_ = 0.0;
}


// Overridden from LinearOpBase


RCP< const VectorSpaceBase<double> >
AztecOOLinearOpWithSolve::range() const
{
  return ( fwdOp_.get() ? fwdOp_->range() : Teuchos::null );
}


RCP< const VectorSpaceBase<double> >
AztecOOLinearOpWithSolve::domain() const
{
  return ( fwdOp_.get() ? fwdOp_->domain() : Teuchos::null );
}


RCP<const LinearOpBase<double> >
AztecOOLinearOpWithSolve::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}


// Overridden from Teuchos::Describable


std::string AztecOOLinearOpWithSolve::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (fwdOp_.get()) {
    oss << "{";
    oss << "fwdOp="<<fwdOp_->description()<<"";
    oss << "}";
  }
  return oss.str();
}


void AztecOOLinearOpWithSolve::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::OSTab;
  using Teuchos::typeName;
  using Teuchos::describe;
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      out
        << Teuchos::Describable::description() << "{"
        << "rangeDim=" << this->range()->dim()
        << ",domainDim="<< this->domain()->dim() << "}\n";
      OSTab tab(out);
      if (!is_null(fwdOp_)) {
        out << "fwdOp = " << describe(*fwdOp_,verbLevel);
      }
      if (!is_null(prec_)) {
        out << "prec = " << describe(*prec_,verbLevel);
      }
      if (!is_null(aztecFwdSolver_)) {
        if (aztecFwdSolver_->GetUserOperator())
          out
            << "Aztec Fwd Op = "
            << typeName(*aztecFwdSolver_->GetUserOperator()) << "\n";
        if (aztecFwdSolver_->GetUserMatrix())
          out
            << "Aztec Fwd Mat = "
            << typeName(*aztecFwdSolver_->GetUserMatrix()) << "\n";
        if (aztecFwdSolver_->GetPrecOperator())
          out
            << "Aztec Fwd Prec Op = "
            << typeName(*aztecFwdSolver_->GetPrecOperator()) << "\n";
        if (aztecFwdSolver_->GetPrecMatrix())
          out
            << "Aztec Fwd Prec Mat = "
            << typeName(*aztecFwdSolver_->GetPrecMatrix()) << "\n";
      }
      if (!is_null(aztecAdjSolver_)) {
        if (aztecAdjSolver_->GetUserOperator())
          out
            << "Aztec Adj Op = "
            << typeName(*aztecAdjSolver_->GetUserOperator()) << "\n";
        if (aztecAdjSolver_->GetUserMatrix())
          out
            << "Aztec Adj Mat = "
            << typeName(*aztecAdjSolver_->GetUserMatrix()) << "\n";
        if (aztecAdjSolver_->GetPrecOperator())
          out
            << "Aztec Adj Prec Op = "
            << typeName(*aztecAdjSolver_->GetPrecOperator()) << "\n";
        if (aztecAdjSolver_->GetPrecMatrix())
          out
            << "Aztec Adj Prec Mat = "
            << typeName(*aztecAdjSolver_->GetPrecMatrix()) << "\n";
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// protected


// Overridden from LinearOpBase


bool AztecOOLinearOpWithSolve::opSupportedImpl(EOpTransp M_trans) const
{
  return ::Thyra::opSupported(*fwdOp_,M_trans);
}


void AztecOOLinearOpWithSolve::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<double> &X,
  const Ptr<MultiVectorBase<double> > &Y,
  const double alpha,
  const double beta
  ) const
{
  Thyra::apply( *fwdOp_, M_trans, X, Y, alpha, beta );
}


// Overridden from LinearOpWithSolveBase


bool
AztecOOLinearOpWithSolve::solveSupportsImpl(EOpTransp M_trans) const
{
  if (real_trans(M_trans)==NOTRANS) return true;
  return (nonnull(aztecAdjSolver_));
}


bool
AztecOOLinearOpWithSolve::solveSupportsSolveMeasureTypeImpl(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType
  ) const
{
  if (real_trans(M_trans)==NOTRANS) {
    if (solveMeasureType.useDefault())
    {
      return true;
    }
    else if (
      solveMeasureType(
        SOLVE_MEASURE_NORM_RESIDUAL,
        SOLVE_MEASURE_NORM_RHS
        )
      &&
      allowInexactFwdSolve_
      )
    {
      return true;
    }
    else if (
      solveMeasureType(
        SOLVE_MEASURE_NORM_RESIDUAL,
        SOLVE_MEASURE_NORM_INIT_RESIDUAL
        )
      &&
      allowInexactFwdSolve_
      )
    {
      return true;
    }
  }
  else {
    // TRANS
    if (aztecAdjSolver_.get()==NULL)
    {
      return false;
    }
    else if (solveMeasureType.useDefault())
    {
      return true;
    }
    else if (
      solveMeasureType(
        SOLVE_MEASURE_NORM_RESIDUAL,
        SOLVE_MEASURE_NORM_RHS
        )
      &&
      allowInexactFwdSolve_
      )
    {
      return true;
    }
    else if (
      solveMeasureType(
        SOLVE_MEASURE_NORM_RESIDUAL,
        SOLVE_MEASURE_NORM_INIT_RESIDUAL
        )
      &&
      allowInexactFwdSolve_
      )
    {
      return true;
    }
  }
  // If you get here then we don't support the solve measure type!
  return false;
}


// Overridden from LinearOpWithSolveBase


SolveStatus<double>
AztecOOLinearOpWithSolve::solveImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<double> &B,
  const Ptr<MultiVectorBase<double> > &X,
  const Ptr<const SolveCriteria<double> > solveCriteria
  ) const
{

  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::rcpFromPtr;
  using Teuchos::OSTab;
  typedef SolveCriteria<double> SC;
  typedef SolveStatus<double> SS;

#ifdef STRATIMIKOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stratimikos: AztecOOLOWS");
#endif
  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  OSTab tab = this->getOSTab();
  if (out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
    *out << "\nSolving block system using AztecOO ...\n\n";

  //
  // Validate input
  //
  TEUCHOS_ASSERT(this->solveSupportsImpl(M_trans));
  SolveMeasureType solveMeasureType;
  if (nonnull(solveCriteria)) {
    solveMeasureType = solveCriteria->solveMeasureType;
    assertSupportsSolveMeasureType(*this, M_trans, solveMeasureType);
  }

  //
  // Get the transpose argument
  //
  const EOpTransp aztecOpTransp = real_trans(M_trans);

  //
  // Get the solver, operator, and preconditioner that we will use
  //
  RCP<AztecOO>
    aztecSolver = ( aztecOpTransp == NOTRANS ? aztecFwdSolver_  : aztecAdjSolver_ );
  const Epetra_Operator
    *aztecOp = aztecSolver->GetUserOperator();

  //
  // Get the op(...) range and domain maps
  //
  const Epetra_Map
    &opRangeMap = aztecOp->OperatorRangeMap(),
    &opDomainMap = aztecOp->OperatorDomainMap();

  //
  // Get the convergence criteria
  //
  double tol = ( aztecOpTransp==NOTRANS ? fwdDefaultTol() : adjDefaultTol() );
  int maxIterations = ( aztecOpTransp==NOTRANS
    ? fwdDefaultMaxIterations() : adjDefaultMaxIterations() );
  bool isDefaultSolveCriteria = true;
  if (nonnull(solveCriteria)) {
    if ( solveCriteria->requestedTol != SC::unspecifiedTolerance() ) {
      tol = solveCriteria->requestedTol;
      isDefaultSolveCriteria = false;
    }
  }

  //
  // Get Epetra_MultiVector views of B and X
  //

  RCP<const Epetra_MultiVector> epetra_B;
  RCP<Epetra_MultiVector> epetra_X;

  const EpetraOperatorWrapper* opWrapper =
    dynamic_cast<const EpetraOperatorWrapper*>(aztecOp);

  if (opWrapper == 0) {
    epetra_B = get_Epetra_MultiVector(opRangeMap, rcpFromRef(B));
    epetra_X = get_Epetra_MultiVector(opDomainMap, rcpFromPtr(X));
  }

  //
  // Use AztecOO to solve each RHS one at a time (which is all that I can do anyway)
  //

  int totalIterations = 0;
  SolveStatus<double> solveStatus;
  solveStatus.solveStatus = SOLVE_STATUS_CONVERGED;
  solveStatus.achievedTol = -1.0;

  /* Get the number of columns in the multivector. We use Thyra
   * functions rather than Epetra functions to do this, as we
   * might not yet have created an Epetra multivector. - KL */
  //const int m = epetra_B->NumVectors();
  const int m = B.domain()->dim();

  for( int j = 0; j < m; ++j ) {

#ifdef STRATIMIKOS_TEUCHOS_TIME_MONITOR
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Stratimikos: AztecOOLOWS:SingleSolve", SingleSolve);
#endif

    //
    // Get Epetra_Vector views of B(:,j) and X(:,j)
    // How this is done will depend on whether we have a true Epetra operator
    // or we are wrapping a general Thyra operator in an Epetra operator.
    //

    // We need to declare epetra_x_j as non-const because when we have a phony
    // Epetra operator we'll have to copy a thyra vector into it.
    RCP<Epetra_Vector> epetra_b_j;
    RCP<Epetra_Vector> epetra_x_j;

    if (opWrapper == 0) {
      epetra_b_j = rcpFromRef(*const_cast<Epetra_Vector*>((*epetra_B)(j)));
      epetra_x_j = rcpFromRef(*(*epetra_X)(j));
    }
    else {
      if (is_null(epetra_b_j)) {
        epetra_b_j = rcp(new Epetra_Vector(opRangeMap));
        epetra_x_j = rcp(new Epetra_Vector(opDomainMap));
      }
      opWrapper->copyThyraIntoEpetra(*B.col(j), *epetra_b_j);
      opWrapper->copyThyraIntoEpetra(*X->col(j), *epetra_x_j);
    }

    //
    // Set the RHS and LHS
    //

    aztecSolver->SetRHS(&*epetra_b_j);
    aztecSolver->SetLHS(&*epetra_x_j);

    //
    // Solve the linear system
    //
    timer.start(true);
    {
      SetAztecSolveState
        setAztecSolveState(aztecSolver,out,verbLevel,solveMeasureType);
      aztecSolver->Iterate( maxIterations, tol );
      // NOTE: We ignore the returned status but get it below
    }
    timer.stop();

    //
    // Scale the solution 
    // (Originally, this was at the end of the loop after all columns had been
    // processed. It's moved here because we need to do it before copying the
    // solution back into a Thyra vector. - KL
    //
    if (aztecSolverScalar_ != 1.0)
      epetra_x_j->Scale(1.0/aztecSolverScalar_);

    //
    // If necessary, convert the solution back to a non-epetra vector
    //
    if (opWrapper != 0) {
      opWrapper->copyEpetraIntoThyra(*epetra_x_j, X->col(j).ptr());
    }

    //
    // Set the return solve status
    //

    const int iterations = aztecSolver->NumIters();
    const double achievedTol = aztecSolver->ScaledResidual();
    const double *AZ_status = aztecSolver->GetAztecStatus();
    std::ostringstream oss;
    bool converged = false;
    if (AZ_status[AZ_why]==AZ_normal) { oss << "Aztec returned AZ_normal."; converged = true; }
    else if (AZ_status[AZ_why]==AZ_param) oss << "Aztec returned AZ_param.";
    else if (AZ_status[AZ_why]==AZ_breakdown) oss << "Aztec returned AZ_breakdown.";
    else if (AZ_status[AZ_why]==AZ_loss) oss << "Aztec returned AZ_loss.";
    else if (AZ_status[AZ_why]==AZ_ill_cond) oss << "Aztec returned AZ_ill_cond.";
    else if (AZ_status[AZ_why]==AZ_maxits) oss << "Aztec returned AZ_maxits.";
    else oss << "Aztec returned an unknown status?";
    oss << "  Iterations = " << iterations << ".";
    oss << "  Achieved Tolerance = " << achievedTol << ".";
    oss << "  Total time = " << timer.totalElapsedTime() << " sec.";
    if (out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE) && outputEveryRhs())
      Teuchos::OSTab(out).o() << "j="<<j<<": " << oss.str() << "\n";
    
    solveStatus.achievedTol = TEUCHOS_MAX(solveStatus.achievedTol, achievedTol);
    // Note, achieveTol may actually be greater than tol due to ill conditioning and roundoff!

    totalIterations += iterations;

    solveStatus.message = oss.str();
    if ( isDefaultSolveCriteria ) {
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

  aztecSolver->UnsetLHSRHS();
  
  //
  // Release the Epetra_MultiVector views of X and B
  //
  epetra_X = Teuchos::null;
  epetra_B = Teuchos::null;

  //
  // Update the overall solve criteria
  //
  totalTimer.stop();
  SolveStatus<double> overallSolveStatus;
  if (isDefaultSolveCriteria) {
    overallSolveStatus.solveStatus = SOLVE_STATUS_UNKNOWN;
    overallSolveStatus.achievedTol = SS::unknownTolerance();
  }
  else {
    overallSolveStatus.solveStatus = solveStatus.solveStatus;
    overallSolveStatus.achievedTol = solveStatus.achievedTol;
  }
  std::ostringstream oss;
  oss
    << "AztecOO solver "
    << ( overallSolveStatus.solveStatus==SOLVE_STATUS_CONVERGED ? "converged" : "unconverged" )
    << " on m = "<<m<<" RHSs using " << totalIterations << " cumulative iterations"
    << " for an average of " << (totalIterations/m) << " iterations/RHS and"
    << " total CPU time of "<<totalTimer.totalElapsedTime()<<" sec.";
  overallSolveStatus.message = oss.str();

  //
  // Report the overall time
  //
  if (out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal solve time = "<<totalTimer.totalElapsedTime()<<" sec\n";

  return overallSolveStatus;

}


}	// end namespace Thyra
