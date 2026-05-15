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
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperEPI<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for EPI2, three for EPI3.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\" for EPI2,\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"3\" for EPI3.\n");

    Thyra::SolveStatus<Scalar> sStatus;

    RCP<StepperEPI<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    // TODO: Figure out why we need this hack:
    // RCP<const Thyra::VectorBase<Scalar> > xOld = currentState->getX();

    RCP<const Thyra::VectorBase<Scalar> > xOld = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > x = workingState->getX();
    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();

    const Scalar time = workingState->getTime();
    const Scalar dt = workingState->getTimeStep();

    // Retrieve x_{n-2} when available (used by EPI3 formula, added in next step).
    // On the first step only two states exist so xOldOld remains null and the
    // method runs as pure EPI2.
    RCP<const Thyra::VectorBase<Scalar> > xOldOld = Teuchos::null;
    if (solutionHistory->getNumStates() >= 3)
      xOldOld = solutionHistory->getStateTimeIndexNM2()->getX();
    // TODO: use xOldOld in EPI3 formula

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
    this->evaluateImplicitODE(Mf, x, xDot, time, p);

    Teuchos::ArrayRCP<RCP<const Thyra::VectorBase<Scalar>>> Mrhs_B(3);
    // Mrhs_B is default initialized with Teuchos::null
    Mrhs_B[1] = Mf;

    // if requested, compute the time derivative for the nonaotonomous correction
    if (temporal_finite_difference_eps_ > 0.0) {
      RCP<Thyra::VectorBase<Scalar>> dt_Mf_deriv = Mf->clone_v();
      this->evaluateImplicitODE(dt_Mf_deriv, x, xDot, time + dt*temporal_finite_difference_eps_, p);
      // compute -dt times the temporal finite difference of Mf, subtract dt_Mf_deriv / eps from Mf / eps
      Scalar one_over_eps = Scalar(1. / temporal_finite_difference_eps_);
      Thyra::linear_combination<Scalar>(Teuchos::tuple(one_over_eps),
                                        Teuchos::tuple(Mf.getConst().ptr()),
                                        -one_over_eps,
                                        dt_Mf_deriv.ptr());
      // compute rhs for the phi_2 function.
      Mrhs_B[2] = dt_Mf_deriv;
    }

    // Using the appModel
    RCP<const Thyra::ModelEvaluator<Scalar>> appModel = this->getModel();
    Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = appModel->createInArgs();
    // Model evaluator builds: alpha*u_dot + beta*F(u) = 0
    inArgs.set_x(x);
    inArgs.set_t(time);
    // set x_dot == 0 to signal to some model evaluators that we want the implicit version
    inArgs.set_x_dot(xDot);

    // initialize vector for the update
    RCP<Thyra::VectorBase<Scalar>> vphi = x->clone_v();
    assign(vphi.ptr(), ST::zero());  // Must initialize to a guess before solve!

    phiEvaluator_->setLinearizationPoint(inArgs);
    // TODO: Avoid using hard coded EPI2 (p=2) and adjust the logic for general p.

    // unless temporal_finite_difference_eps_ is not positive, add the nonautonomous correction term
    if (temporal_finite_difference_eps_ > 0.0) {
      sStatus = phiEvaluator_->computePhis(vphi.ptr(), dt, Mrhs_B());
    }
    else {
      sStatus = phiEvaluator_->computePhi(vphi.ptr(), 1, dt, Mf);
    }

    // std::cout << "xO[0,1] = " << Thyra::get_ele(*xOld, 0) << " " << Thyra::get_ele(*xOld, 1) << std::endl;
    // std::cout << "x[0,1]  = " << Thyra::get_ele(*x, 0) << " " << Thyra::get_ele(*x, 1) << std::endl;
    // std::cout << "vphi[0,1]  = " << Thyra::get_ele(*vphi, 0) << " " << Thyra::get_ele(*vphi, 1) << std::endl;
    Thyra::V_VpStV(x.ptr(), *xOld, Scalar(-1.0)*dt, *vphi);

    // std::cout << sStatus << std::endl;

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
