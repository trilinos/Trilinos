// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperEPI_impl_hpp
#define Tempus_StepperEPI_impl_hpp

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"
#include "Tempus_StepperEPIModifierDefault.hpp"
#include "Tempus_StepperEPI_decl.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Tempus_StepperFactory.hpp"

// TODO: have to include this header to get LSP to work.
#include "Tempus_StepperEPI.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_VectorStdOps_decl.hpp"

namespace Tempus {

template <class Scalar>
StepperEPI<Scalar>::StepperEPI()
{
  this->setStepperName("EPI");
  this->setStepperType("EPI");
  this->setUseFSAL(false);
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);
  this->setZeroInitialGuess(false);

  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}


template<class Scalar>
StepperEPI<Scalar>::StepperEPI(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess,
  const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction)
{
  this->setStepperName("EPI");
  this->setStepperType("EPI");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setZeroInitialGuess(zeroInitialGuess);

  this->setAppAction(stepperEPIAppAction);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}


template<class Scalar>
void StepperEPI<Scalar>::setAppAction(
  Teuchos::RCP<StepperEPIAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperEPIAppAction_ =
      Teuchos::rcp(new StepperEPIModifierDefault<Scalar>());
  } else {
    stepperEPIAppAction_ = appAction;
  }
  this->isInitialized_ = false;
}


template<class Scalar>
void StepperEPI<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  StepperImplicit<Scalar>::setModel(appModel);

  phiEvaluator_->setModel(appModel);
  phiEvaluator_->initialize();

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperEPI<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();

  // Check if we need Stepper storage for xDot
  if (initialState->getXDot() == Teuchos::null)
    this->setStepperXDot(initialState->getX()->clone_v());
  else
    this->setStepperXDot(initialState->getXDot());

  StepperImplicit<Scalar>::setInitialConditions(solutionHistory);
}


template<class Scalar>
void StepperEPI<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperEPI::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION((solutionHistory->getNumStates() < 2 && order_ <= 2.0),
      std::logic_error,
      "Error - StepperEPI<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for EPI2.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\"\n"
      "  \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\" for EPI2.\n"
      "  or \"Storage Type\" = \"Unlimited\"\n");

    TEUCHOS_TEST_FOR_EXCEPTION((solutionHistory->getStorageLimit() < 3 && order_ > 2.0),
      std::logic_error,
      "Error - StepperEPI<Scalar>::takeStep(...)\n"
      "Need at least three SolutionStates for EPI3.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\"\n"
      "  \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"3\" for EPI3.\n"
      "  or \"Storage Type\" = \"Unlimited\"\n");
    bool useEPI3 = false;
    if (solutionHistory->getNumStates() >= 3 && order_ > 2.0) {
      useEPI3 = true;
    }

    Thyra::SolveStatus<Scalar> sStatus;

    RCP<StepperEPI<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    RCP<Thyra::VectorBase<Scalar> > xOld = currentState->getX();
    // TODO: Figure out why we need this hack:
    // We must have applied Dirichlet BCs to the workingState on first step
    if (workingState->getIndex() == 1)
      Thyra::assign(xOld.ptr(), *workingState->getX());

    RCP<Thyra::VectorBase<Scalar> > x = workingState->getX();
    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();

    const Scalar time = workingState->getTime();
    const Scalar t0   = currentState->getTime();
    const Scalar dt = workingState->getTimeStep();

    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXP);

    // TODO: Transition away from using implicit solver methods
    // and use ModelEvaluator directly

    // Setup TimeDerivative
    RCP<TimeDerivative<Scalar> > timeDer;
    timeDer = Teuchos::rcp(new StepperEPITimeDerivative<Scalar>(1/dt, xOld));
    auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(timeDer, dt, Scalar(0.0), Scalar(1.0)));

    // set the right hand side for the Phi_1 function
    RCP<Thyra::VectorBase<Scalar>> Mf = x->clone_v();
    this->evaluateImplicitODE(Mf, xOld, xDot, t0, p);

    Teuchos::ArrayRCP<RCP<const Thyra::VectorBase<Scalar>>> Mrhs_B(3);
    // Mrhs_B is default initialized with Teuchos::null
    Mrhs_B[1] = Mf;

    // Using the appModel
    RCP<const Thyra::ModelEvaluator<Scalar>> appModel = this->getModel();
    Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = appModel->createInArgs();
    // Model evaluator builds: alpha*u_dot + beta*F(u) = 0
    inArgs.set_x(xOld);
    inArgs.set_t(t0);
    // set x_dot == 0 to signal to some model evaluators that we want the implicit version
    inArgs.set_x_dot(xDot);

    // initialize vector for the update
    RCP<Thyra::VectorBase<Scalar>> vphi = x->clone_v();
    assign(vphi.ptr(), ST::zero());  // Must initialize to a guess before solve!

    // setup system Jacobian computation at the current time
    phiEvaluator_->setLinearizationPoint(inArgs);

    // if requested, compute the time derivative for the nonaotonomous correction
    RCP<Thyra::VectorBase<Scalar>> dt_Mf_deriv = Teuchos::null;
    if (temporal_finite_difference_eps_ > 0.0) {
      dt_Mf_deriv = Mf->clone_v();
      this->evaluateImplicitODE(
          dt_Mf_deriv, xOld, xDot, t0 + dt*temporal_finite_difference_eps_, p);
      // compute dt times the temporal finite difference of Mf:
      // dt_Mf_deriv = (Mf(t + eps*dt) - Mf(t)) / eps
      Scalar one_over_eps = Scalar(1. / temporal_finite_difference_eps_);
      Thyra::linear_combination<Scalar>(Teuchos::tuple(-one_over_eps),
                                        Teuchos::tuple(Mf.getConst().ptr()),
                                        one_over_eps,
                                        dt_Mf_deriv.ptr());
      // compute rhs for the phi_2 function.
      Mrhs_B[2] = dt_Mf_deriv;
    }

    // if requested, use the EPI3 3rd order update
    if (useEPI3) {
      // Retrieve x_{n-1} when available (used by EPI3 formula, added in next step).
      // Note this last solution is obtined with getStateTimeIndexNM2()
      // On the first step only two states exist so
      // method runs as pure EPI2.
      Scalar tOldOld = solutionHistory->getStateTimeIndexNM2()->getTime();
      RCP<const Thyra::VectorBase<Scalar>> xOldOld = solutionHistory->getStateTimeIndexNM2()->getX();
      RCP<Thyra::VectorBase<Scalar>> xDotOldOld = solutionHistory->getStateTimeIndexNM2()->getXDot();
      RCP<Thyra::VectorBase<Scalar>> Remf = computeRemf(
          xOldOld, xDotOldOld, tOldOld, xOld, t0, dt, Mf, dt_Mf_deriv);
      // Subtract (2/3)*R from phi_2 term
      Thyra::scale((-2.0 / 3.0), Remf.ptr());
      if (Mrhs_B[2] != Teuchos::null) {
        Thyra::Vp_V(Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar>>(Mrhs_B[2]).ptr(), *Remf);
      } else {
        Mrhs_B[2] = Remf;
      };
    }

    // if temporal_finite_difference_eps_ is positive or we use EPI3,
    // use the extension formula.
    if (temporal_finite_difference_eps_ > 0.0 || useEPI3) {
      sStatus = phiEvaluator_->computePhis(vphi.ptr(), dt, Mrhs_B());
    }
    else {
      sStatus = phiEvaluator_->computePhi(vphi.ptr(), 1, dt, Mf);
    }

    // TODO: make this configurable
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    // Teuchos::OSTab ostab(out, 1, "StepperEPI::takeStep");
    int current_iters = -1;
    if(!sStatus.extraParameters.is_null()) {
      current_iters = sStatus.extraParameters->get("Iteration Count", 0);
    }
    Scalar achieved_tol = sStatus.achievedTol;

    if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED && current_iters >=0) {
      *out << "Phi converged: iters: " << current_iters << " tol: " << achieved_tol << std::endl;
    }
    else {
      *out << sStatus.message << std::endl;
    }

    Thyra::V_VpStV(x.ptr(), *xOld, Scalar(-1.0)*dt, *vphi);

    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::AFTER_EXP);

    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::END_STEP);
  }
  return;
}


template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar>>
StepperEPI<Scalar>::computeRemf(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& xr,
    const Teuchos::RCP<Thyra::VectorBase<Scalar>>& xrDot,
    const Scalar tr,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x0,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& dt_Mf_deriv
)
{
  Teuchos::RCP<TimeDerivative<Scalar> > timeDer;
  timeDer = Teuchos::rcp(new StepperEPITimeDerivative<Scalar>(1/dt, xr));
  auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(timeDer, dt, Scalar(0.0), Scalar(1.0)));

  // Eval the rhs at (xr, tr)
  // Mf_old = -M*F(xr, tr) = F_impl(xDot = 0, xr, tr)
  Teuchos::RCP<Thyra::VectorBase<Scalar>> Mf_old = Thyra::createMember(xr->space());
  this->evaluateImplicitODE(Mf_old, Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar>>(xr), xrDot, tr, p);

  // xd = xr - x0
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xd = Thyra::createMember(xr->space());
  Thyra::V_VpStV(xd.ptr(), *xr, Scalar(-1.0), *x0);

  // J_xd = -M*J * (xr - x0)
  Teuchos::RCP<Thyra::VectorBase<Scalar>> MJ_xd = Thyra::createMember(x0->space());
  phiEvaluator_->applyJacobian(MJ_xd.ptr(), xd);

  // R = (Mf - Mf_old) + J_xd
  // Mf is rhs at the current time: -M*F(x0, t0)
  Thyra::Vp_V(MJ_xd.ptr(), *Mf);
  Thyra::Vp_StV(MJ_xd.ptr(), Scalar(-1.0), *Mf_old);

  // Time derivative remainder term only nonzero in nonautonomous case
  // add  (M * F') * (tr - t0) = (dt * M * F') * ((tr - t0) / dt)
  if (dt_Mf_deriv != Teuchos::null) {
    Thyra::Vp_StV(MJ_xd.ptr(), Scalar((tr - t0) / dt), *dt_Mf_deriv);
  }

  // TODO: let the outer function provide the memory.
  return MJ_xd;
}


/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperEPI<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperEPI<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperEPI ---\n";
  out << "  stepperEPIAppAction_                = "
      << stepperEPIAppAction_ << std::endl;
  out << "----------------------------" << std::endl;
}


template<class Scalar>
bool StepperEPI<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperImplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (stepperEPIAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The EPI AppAction is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperEPI<Scalar>::getValidParameters() const
{
  auto pl = this->getValidParametersBasicImplicit();
  return pl;
}


template <class Scalar>
void StepperEPI<Scalar>::setStepperExponentialValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  Teuchos::RCP<Teuchos::ParameterList> phiPL = Teuchos::null;

  temporal_finite_difference_eps_ = pl->get<double>("Epsilon for RHS finite difference", 1e-4);
  order_ = pl->get<Scalar>("EPI Order", 2.0);

  if (pl != Teuchos::null) {
    // TODO read in the pl for the exponential solver
    phiPL = sublist(pl, "PhiEvaluator");
  }

  // TODO: is this the right place to initialize the PhiEvaluator?
  auto phif = Teuchos::rcp(new PhiEvaluatorFactory<Scalar>());
  phiEvaluator_ = phif->createPhiEvaluator(phiPL);
}


// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template<class Scalar>
Teuchos::RCP<StepperEPI<Scalar> >
createStepperEPI(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperEPI<Scalar>());

  stepper->setStepperImplicitValues(pl);

  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}


} // namespace Tempus
#endif // Tempus_StepperEPI_impl_hpp
