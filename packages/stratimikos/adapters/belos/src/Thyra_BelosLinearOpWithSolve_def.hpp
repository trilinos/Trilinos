
#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_BelosLinearOpWithSolve_decl.hpp"
#include "Thyra_GeneralSolveCriteriaBelosStatusTest.hpp"
#include "Thyra_LinearOpWithSolveHelpers.hpp"
#include "Teuchos_DebugDefaultAsserts.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace {
  // Set the Belos solver's parameter list to scale its residual norms
  // in the specified way.
  //
  // We break this out in a separate function because the parameters
  // to set depend on which parameters the Belos solver supports.  Not
  // all Belos solvers support both the "Implicit Residual Scaling"
  // and "Explicit Residual Scaling" parameters, so we have to check
  // the solver's list of valid parameters for the existence of these.
  //
  // Scaling options: Belos lets you decide whether the solver will
  // scale residual norms by the (left-)preconditioned initial
  // residual norms (residualScalingType = "Norm of Initial
  // Residual"), or by the unpreconditioned initial residual norms
  // (residualScalingType = "Norm of RHS").  Usually you want to scale
  // by the unpreconditioned initial residual norms.  This is because
  // preconditioning is just an optimization, and you really want to
  // make ||B - AX|| small, rather than ||M B - M (A X)||.  If you're
  // measuring ||B - AX|| and scaling by the initial residual, you
  // should use the unpreconditioned initial residual to match it.
  //
  // Note, however, that the implicit residual test computes
  // left-preconditioned residuals, if a left preconditioner was
  // provided.  That's OK because when Belos solvers (at least the
  // GMRES variants) are given a left preconditioner, they first check
  // the implicit residuals.  If those converge, they then check the
  // explicit residuals.  The explicit residual test does _not_ apply
  // the left preconditioner when computing the residual.  The
  // implicit residual test is just an optimization so that Belos
  // doesn't have to compute explicit residuals B - A*X at every
  // iteration.  This is why we use the same scaling factor for both
  // the implicit and explicit residuals.
  //
  // Arguments:
  //
  // solverParams [in/out] Parameters for the current solve.
  //
  // solverValidParams [in] Valid parameter list for the Belos solver.
  //   Result of calling the solver's getValidParameters() method.
  //
  // residualScalingType [in] String describing how the solver should
  //   scale residuals.  Valid values include "Norm of RHS" and "Norm
  //   of Initial Residual" (these are the only two options this file
  //   currently uses, though Belos offers other options).
  void
  setResidualScalingType (const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
			  const Teuchos::RCP<const Teuchos::ParameterList>& solverValidParams,
			  const std::string& residualScalingType)
  {
    // Many Belos solvers (especially the GMRES variants) define both
    // "Implicit Residual Scaling" and "Explicit Residual Scaling"
    // options.
    //
    // "Implicit" means "the left-preconditioned approximate
    // a.k.a. 'recursive' residual as computed by the Krylov method."
    // 
    // "Explicit" means ||B - A*X||, the unpreconditioned, "exact"
    // residual.
    //
    // Belos' GMRES implementations chain these two tests in sequence.
    // Implicit comes first, and explicit is not evaluated unless
    // implicit passes.  In some cases (e.g., no left preconditioner),
    // GMRES _only_ uses the implicit tests.  This means that only
    // setting "Explicit Residual Scaling" won't change the solver's
    // behavior.  Stratimikos tends to prefer using a right
    // preconditioner, in which case setting only the "Explicit
    // Residual Scaling" argument has no effect.  Furthermore, if
    // "Explicit Residual Scaling" is set to something other than the
    // default (initial residual norm), without "Implicit Residual
    // Scaling" getting the same setting, then the implicit residual
    // test will be using a radically different scaling factor than
    // the user wanted.
    // 
    // Not all Belos solvers support both options.  We check the
    // solver's valid parameter list first before attempting to set
    // the option.
    if (solverValidParams->isParameter ("Implicit Residual Scaling")) {
      solverParams->set ("Implicit Residual Scaling", residualScalingType);
    }
    if (solverValidParams->isParameter ("Explicit Residual Scaling")) {
      solverParams->set ("Explicit Residual Scaling", residualScalingType);
    }
  }

} // namespace (anonymous)


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
BelosLinearOpWithSolve<Scalar>::BelosLinearOpWithSolve()
  :convergenceTestFrequency_(-1),
  isExternalPrec_(false),
  supportSolveUse_(SUPPORT_SOLVE_UNSPECIFIED),
  defaultTol_(-1.0)
{}


template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::initialize(
  const RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> > &lp,
  const RCP<Teuchos::ParameterList> &solverPL,
  const RCP<Belos::SolverManager<Scalar,MV_t,LO_t> > &iterativeSolver,
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const PreconditionerBase<Scalar> > &prec,
  const bool isExternalPrec,
  const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
  const ESupportSolveUse &supportSolveUse,
  const int convergenceTestFrequency
  )
{
  this->setLinePrefix("BELOS/T");
  // ToDo: Validate input
  lp_ = lp;
  solverPL_ = solverPL;
  iterativeSolver_ = iterativeSolver;
  fwdOpSrc_ = fwdOpSrc;
  prec_ = prec;
  isExternalPrec_ = isExternalPrec;
  approxFwdOpSrc_ = approxFwdOpSrc;
  supportSolveUse_ = supportSolveUse;
  convergenceTestFrequency_ = convergenceTestFrequency;
  // Check if "Convergence Tolerance" is in the solver parameter list.  If
  // not, use the default from the solver.
  if ( !is_null(solverPL_) ) {
    if (solverPL_->isParameter("Convergence Tolerance")) {
      defaultTol_ = solverPL_->get<double>("Convergence Tolerance");
    }
  }
  else {
    RCP<const Teuchos::ParameterList> defaultPL =
      iterativeSolver->getValidParameters();
    defaultTol_ = defaultPL->get<double>("Convergence Tolerance");
  }
}


template<class Scalar>
RCP<const LinearOpSourceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::extract_fwdOpSrc()
{
  RCP<const LinearOpSourceBase<Scalar> >
    _fwdOpSrc = fwdOpSrc_;
  fwdOpSrc_ = Teuchos::null;
  return _fwdOpSrc;
}


template<class Scalar>
RCP<const PreconditionerBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::extract_prec()
{
  RCP<const PreconditionerBase<Scalar> >
    _prec = prec_;
  prec_ = Teuchos::null;
  return _prec;
}


template<class Scalar>
bool BelosLinearOpWithSolve<Scalar>::isExternalPrec() const
{
  return isExternalPrec_;
}


template<class Scalar>
RCP<const LinearOpSourceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::extract_approxFwdOpSrc()
{
  RCP<const LinearOpSourceBase<Scalar> >
    _approxFwdOpSrc = approxFwdOpSrc_;
  approxFwdOpSrc_ = Teuchos::null;
  return _approxFwdOpSrc;
}


template<class Scalar>
ESupportSolveUse BelosLinearOpWithSolve<Scalar>::supportSolveUse() const
{
  return supportSolveUse_;
}


template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::uninitialize(
  RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> > *lp,
  RCP<Teuchos::ParameterList> *solverPL,
  RCP<Belos::SolverManager<Scalar,MV_t,LO_t> > *iterativeSolver,
  RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
  RCP<const PreconditionerBase<Scalar> > *prec,
  bool *isExternalPrec,
  RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
  ESupportSolveUse *supportSolveUse
  )
{
  if (lp) *lp = lp_;
  if (solverPL) *solverPL = solverPL_;
  if (iterativeSolver) *iterativeSolver = iterativeSolver_;
  if (fwdOpSrc) *fwdOpSrc = fwdOpSrc_;
  if (prec) *prec = prec_;
  if (isExternalPrec) *isExternalPrec = isExternalPrec_;
  if (approxFwdOpSrc) *approxFwdOpSrc = approxFwdOpSrc_;
  if (supportSolveUse) *supportSolveUse = supportSolveUse_;

  lp_ = Teuchos::null;
  solverPL_ = Teuchos::null;
  iterativeSolver_ = Teuchos::null;
  fwdOpSrc_ = Teuchos::null;
  prec_ = Teuchos::null;
  isExternalPrec_ = false;
  approxFwdOpSrc_ = Teuchos::null;
  supportSolveUse_ = SUPPORT_SOLVE_UNSPECIFIED;
}


// Overridden from LinearOpBase


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::range() const
{
  if (!is_null(lp_))
    return lp_->getOperator()->range();
  return Teuchos::null;
}


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::domain() const
{
  if (!is_null(lp_))
    return lp_->getOperator()->domain();
  return Teuchos::null;
}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}


// Overridden from Teuchos::Describable


template<class Scalar>
std::string BelosLinearOpWithSolve<Scalar>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if ( !is_null(lp_) && !is_null(lp_->getOperator()) ) {
    oss << "{";
    oss << "iterativeSolver=\'"<<iterativeSolver_->description()<<"\'";
    oss << ",fwdOp=\'"<<lp_->getOperator()->description()<<"\'";
    if (lp_->getLeftPrec().get())
      oss << ",leftPrecOp=\'"<<lp_->getLeftPrec()->description()<<"\'";
    if (lp_->getRightPrec().get())
      oss << ",rightPrecOp=\'"<<lp_->getRightPrec()->description()<<"\'";
    oss << "}";
  }
  // ToDo: Make Belos::SolverManager derive from Teuchos::Describable so
  // that we can get better information.
  return oss.str();
}


template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::describe(
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch (verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << Teuchos::Describable::description()<< "{"
        << "rangeDim=" << this->range()->dim()
        << ",domainDim=" << this->domain()->dim() << "}\n";
      if (lp_->getOperator().get()) {
        OSTab tab(out);
        *out
          << "iterativeSolver = "<<describe(*iterativeSolver_,verbLevel)
          << "fwdOp = " << describe(*lp_->getOperator(),verbLevel);
        if (lp_->getLeftPrec().get())
          *out << "leftPrecOp = "<<describe(*lp_->getLeftPrec(),verbLevel);
        if (lp_->getRightPrec().get())
          *out << "rightPrecOp = "<<describe(*lp_->getRightPrec(),verbLevel);
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// protected


// Overridden from LinearOpBase


template<class Scalar>
bool BelosLinearOpWithSolve<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  return ::Thyra::opSupported(*lp_->getOperator(),M_trans);
}


template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  ::Thyra::apply<Scalar>(*lp_->getOperator(), M_trans, X, Y, alpha, beta);
}


// Overridden from LinearOpWithSolveBase


template<class Scalar>
bool
BelosLinearOpWithSolve<Scalar>::solveSupportsImpl(EOpTransp M_trans) const
{
  return solveSupportsNewImpl(M_trans, Teuchos::null);
}


template<class Scalar>
bool
BelosLinearOpWithSolve<Scalar>::solveSupportsNewImpl(EOpTransp transp,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria) const
{
  // Only support forward solve right now!
  if (real_trans(transp)==NOTRANS) return true;
  return false; // ToDo: Support adjoint solves!
  // Otherwise, Thyra/Belos now supports every solve criteria type that exists
  // because of the class Thyra::GeneralSolveCriteriaBelosStatusTest!
  return true;
/*
  if (real_trans(M_trans)==NOTRANS) {
    return (
      solveMeasureType.useDefault()
      ||
      solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS)
      ||
      solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_INIT_RESIDUAL)
      );
  }
*/
}


template<class Scalar>
bool
BelosLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureTypeImpl(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType) const
{
  SolveCriteria<Scalar> solveCriteria(solveMeasureType, SolveCriteria<Scalar>::unspecifiedTolerance());
  return solveSupportsNewImpl(M_trans, Teuchos::constOptInArg(solveCriteria));
}


template<class Scalar>
SolveStatus<Scalar>
BelosLinearOpWithSolve<Scalar>::solveImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &B,
  const Ptr<MultiVectorBase<Scalar> > &X,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria
  ) const
{

  THYRA_FUNC_TIME_MONITOR("Stratimikos: BelosLOWS");

  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::rcpFromPtr;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::describe;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  assertSolveSupports(*this, M_trans, solveCriteria);
  // 2010/08/22: rabartl: Bug 4915 ToDo: Move the above into the NIV function
  // solve(...).

  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  OSTab tab = this->getOSTab();
  if (out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_LOW)) {
    *out << "\nStarting iterations with Belos:\n";
    OSTab tab2(out);
    *out << "Using forward operator = " << describe(*fwdOpSrc_->getOp(),verbLevel);
    *out << "Using iterative solver = " << describe(*iterativeSolver_,verbLevel);
    *out << "With #Eqns="<<B.range()->dim()<<", #RHSs="<<B.domain()->dim()<<" ...\n";
  }

  //
  // Set RHS and LHS
  //

  bool ret = lp_->setProblem( rcpFromPtr(X), rcpFromRef(B) );
  TEUCHOS_TEST_FOR_EXCEPTION(
    ret == false, CatastrophicSolveFailure
    ,"Error, the Belos::LinearProblem could not be set for the current solve!"
    );

  //
  // Set the solution criteria
  //

  // Parameter list for the current solve.
  const RCP<ParameterList> tmpPL = Teuchos::parameterList();

  // The solver's valid parameter list.
  RCP<const ParameterList> validPL = iterativeSolver_->getValidParameters();

  SolveMeasureType solveMeasureType;
  RCP<GeneralSolveCriteriaBelosStatusTest<Scalar> > generalSolveCriteriaBelosStatusTest;
  if (nonnull(solveCriteria)) {
    solveMeasureType = solveCriteria->solveMeasureType;
    const ScalarMag requestedTol = solveCriteria->requestedTol;
    if (solveMeasureType.useDefault()) {
      tmpPL->set("Convergence Tolerance", defaultTol_);
    }
    else if (solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL, SOLVE_MEASURE_NORM_RHS)) {
      if (requestedTol != SolveCriteria<Scalar>::unspecifiedTolerance()) {
        tmpPL->set("Convergence Tolerance", requestedTol);
      }
      else {
        tmpPL->set("Convergence Tolerance", defaultTol_);
      }
      setResidualScalingType (tmpPL, validPL, "Norm of RHS");
    }
    else if (solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL, SOLVE_MEASURE_NORM_INIT_RESIDUAL)) {
      if (requestedTol != SolveCriteria<Scalar>::unspecifiedTolerance()) {
        tmpPL->set("Convergence Tolerance", requestedTol);
      }
      else {
        tmpPL->set("Convergence Tolerance", defaultTol_);
      }
      setResidualScalingType (tmpPL, validPL, "Norm of Initial Residual");
    }
    else {
      // Set the most generic (and inefficient) solve criteria
      generalSolveCriteriaBelosStatusTest = createGeneralSolveCriteriaBelosStatusTest(
        *solveCriteria, convergenceTestFrequency_);
      // Set the verbosity level (one level down)
      generalSolveCriteriaBelosStatusTest->setOStream(out);
      generalSolveCriteriaBelosStatusTest->setVerbLevel(incrVerbLevel(verbLevel, -1));
      // Set the default convergence tolerance to always converged to allow
      // the above status test to control things.
      tmpPL->set("Convergence Tolerance", 1.0);
    }
    // maximum iterations
    if (nonnull(solveCriteria->extraParameters)) {
      if (Teuchos::isParameterType<int>(*solveCriteria->extraParameters,"Maximum Iterations")) {
        tmpPL->set("Maximum Iterations", Teuchos::get<int>(*solveCriteria->extraParameters,"Maximum Iterations"));
      }
    }
  }
  else {
    // No solveCriteria was even passed in!
    tmpPL->set("Convergence Tolerance", defaultTol_);
  }

  //
  // Solve the linear system
  //

  Belos::ReturnType belosSolveStatus;
  {
    RCP<std::ostream>
      outUsed =
      ( static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_LOW)
        ? out
        : rcp(new FancyOStream(rcp(new Teuchos::oblackholestream())))
        );
    Teuchos::OSTab tab(outUsed,1,"BELOS");
    tmpPL->set("Output Stream", outUsed);
    iterativeSolver_->setParameters(tmpPL);
    if (nonnull(generalSolveCriteriaBelosStatusTest)) {
      iterativeSolver_->setUserConvStatusTest(generalSolveCriteriaBelosStatusTest);
    }
    belosSolveStatus = iterativeSolver_->solve();
  }

  //
  // Report the solve status
  //

  totalTimer.stop();

  SolveStatus<Scalar> solveStatus;

  switch (belosSolveStatus) {
    case Belos::Unconverged: {
      solveStatus.solveStatus = SOLVE_STATUS_UNCONVERGED;
      // Set achievedTol even if the solver did not converge.  This is
      // helpful for things like nonlinear solvers, which might be
      // able to use a partially converged result, and which would
      // like to know the achieved convergence tolerance for use in
      // computing bounds.  It's also helpful for estimating whether a
      // small increase in the maximum iteration count might be
      // helpful next time.
      try {
	// Some solvers might not have implemented achievedTol(). 
	// The default implementation throws std::runtime_error.
	solveStatus.achievedTol = iterativeSolver_->achievedTol();
      } catch (std::runtime_error&) {
	// Do nothing; use the default value of achievedTol.
      }
      break;
    }
    case Belos::Converged: {
      solveStatus.solveStatus = SOLVE_STATUS_CONVERGED;
      if (nonnull(generalSolveCriteriaBelosStatusTest)) {
	// The user set a custom status test.  This means that we
	// should ask the custom status test itself, rather than the
	// Belos solver, what the final achieved convergence tolerance
	// was.
        const ArrayView<const ScalarMag> achievedTol = 
          generalSolveCriteriaBelosStatusTest->achievedTol();
        solveStatus.achievedTol = ST::zero();
        for (Ordinal i = 0; i < achievedTol.size(); ++i) {
          solveStatus.achievedTol = std::max(solveStatus.achievedTol, achievedTol[i]);
        }
      }
      else {
	try {
	  // Some solvers might not have implemented achievedTol(). 
	  // The default implementation throws std::runtime_error.
	  solveStatus.achievedTol = iterativeSolver_->achievedTol();
	} catch (std::runtime_error&) {
	  // Use the default convergence tolerance.  This is a correct
	  // upper bound, since we did actually converge.
	  solveStatus.achievedTol = tmpPL->get("Convergence Tolerance", defaultTol_);
	}
      }
      break;
    }
    TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT();
  }

  std::ostringstream ossmessage;
  ossmessage
    << "The Belos solver of type \""<<iterativeSolver_->description()
    <<"\" returned a solve status of \""<< toString(solveStatus.solveStatus) << "\""
    << " in " << iterativeSolver_->getNumIters() << " iterations"
    << " with total CPU time of " << totalTimer.totalElapsedTime() << " sec" ;
  if (out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
    *out << "\n" << ossmessage.str() << "\n";

  solveStatus.message = ossmessage.str();

  // Dump the getNumIters() and the achieved convergence tolerance
  // into solveStatus.extraParameters, as the "Belos/Iteration Count"
  // resp. "Belos/Achieved Tolerance" parameters.
  if (solveStatus.extraParameters.is_null()) {
    solveStatus.extraParameters = parameterList ();
  }
  solveStatus.extraParameters->set ("Belos/Iteration Count", 
				    iterativeSolver_->getNumIters());\
  // package independent version of the same
  solveStatus.extraParameters->set ("Iteration Count", 
				    iterativeSolver_->getNumIters());\
  // NOTE (mfh 13 Dec 2011) Though the most commonly used Belos
  // solvers do implement achievedTol(), some Belos solvers currently
  // do not.  In the latter case, if the solver did not converge, the
  // reported achievedTol() value may just be the default "invalid"
  // value -1, and if the solver did converge, the reported value will
  // just be the convergence tolerance (a correct upper bound).
  solveStatus.extraParameters->set ("Belos/Achieved Tolerance", 
				    solveStatus.achievedTol);

//  This information is in the previous line, which is printed anytime the verbosity
//  is not set to Teuchos::VERB_NONE, so I'm commenting this out for now.
//  if (out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
//    *out << "\nTotal solve time in Belos = "<<totalTimer.totalElapsedTime()<<" sec\n";

  return solveStatus;

}


}	// end namespace Thyra


#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP
