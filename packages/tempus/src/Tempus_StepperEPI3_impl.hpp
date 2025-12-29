// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperEPI3_impl_hpp
#define Tempus_StepperEPI3_impl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"
#include "Tempus_StepperEPI3ModifierDefault.hpp"
#include "Tempus_StepperEPI3_decl.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Tempus_StepperFactory.hpp"

// TODO: have to include this header to get LSP to work.
#include "Tempus_StepperEPI3.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_OperatorVectorTypes.hpp"

namespace Tempus {


template<class Scalar>
StepperEPI3<Scalar>::StepperEPI3()
{
  this->setStepperName("EPI3");
  this->setStepperType("EPI3");
  this->setUseFSAL(false);
  // TODO: think about the right setting here
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);
  this->setZeroInitialGuess(false);

  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}


template<class Scalar>
StepperEPI3<Scalar>::StepperEPI3(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess,
  const Teuchos::RCP<StepperEPI3AppAction<Scalar> >& stepperEEAppAction)
{
  this->setStepperName("EPI3");
  this->setStepperType("EPI3");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setZeroInitialGuess(zeroInitialGuess);

  this->setAppAction(stepperEEAppAction);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}


template<class Scalar>
void StepperEPI3<Scalar>::setAppAction(
  Teuchos::RCP<StepperEPI3AppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperEEAppAction_ =
      Teuchos::rcp(new StepperEPI3ModifierDefault<Scalar>());
  } else {
    stepperEEAppAction_ = appAction;
  }
  this->isInitialized_ = false;
}


template<class Scalar>
void StepperEPI3<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  StepperImplicit<Scalar>::setModel(appModel);

  phiEvaluator_->setModel(appModel);
  phiEvaluator_->initialize();

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperEPI3<Scalar>::setInitialConditions(
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
void StepperEPI3<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();
  phiEvaluator_->checkInitialized();

  using Teuchos::RCP;

  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperEPI3::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperEPI3<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for EPI3.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<StepperEPI3<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperEPI3AppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    RCP<const Thyra::VectorBase<Scalar> > xOld = currentState->getX();
    RCP<Thyra::VectorBase<Scalar> > x = workingState->getX();
    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();

    const Scalar time = workingState->getTime();
    const Scalar dt = workingState->getTimeStep();

    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperEPI3AppAction<Scalar>::ACTION_LOCATION::BEFORE_EXP);

    //{
    //  Teuchos::basic_FancyOStream<char> ostr(Teuchos::rcp(&std::cout, false));
    //  this->describe(ostr, Teuchos::VERB_EXTREME);
    //}
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer;

    bool exponential = true;

    Thyra::SolveStatus<Scalar> sStatus;
    if (!exponential)
    {
      // Setup TimeDerivative
      timeDer = Teuchos::rcp(new StepperEPI3TimeDerivative<Scalar>(
        Scalar(1.0)/dt,xOld));

      const Scalar alpha = Scalar(1.0)/dt;
      const Scalar beta  = Scalar(1.0);
      auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
      timeDer, dt, alpha, beta));

      sStatus =
	this->solveImplicitODE(x, xDot, time, p);
    }
    else
    {
      // Setup TimeDerivative
      timeDer = Teuchos::rcp(new StepperEPI3TimeDerivative<Scalar>(
        Scalar(0.0),xOld));
      auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
          timeDer, dt, Scalar(0.0), Scalar(1.0)));


      // TODO: Transition away from using implicit solver methods and use ModelEvaluator directly
      RCP<Thyra::VectorBase<Scalar> > f = x->clone_v();
      this->evaluateImplicitODE(f, x, xDot, time, p);

      // Using the appModel
      RCP<const Thyra::ModelEvaluator<Scalar>> appModel = this->getModel();
      Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = appModel->createInArgs();
      // Model evaluator builds: alpha*u_dot + beta*F(u) = 0
      inArgs.set_x(x);
      inArgs.set_t(time);
      inArgs.set_x_dot(xDot); // for what? xDot is zero at this point, updating of it not decided

      // initialize space for the update
      RCP<Thyra::VectorBase<Scalar>> vphi = x->clone_v();
      assign(vphi.ptr(), ST::zero());  // Must initialize to a guess before solve!

      bool use_phi_eval = true;
      Scalar factor = Scalar(-dt);

auto thyraModel = this->getModel();

auto inArgs_tyra  = thyraModel->createInArgs();
auto outArgs_tyra = thyraModel->createOutArgs();

auto curr = solutionHistory->getCurrentState();

inArgs_tyra.set_x(curr->getX());
inArgs_tyra.set_t(curr->getTime());

auto W = thyraModel->create_W_op();
outArgs_tyra.set_W_op(W);

thyraModel->evalModel(inArgs_tyra, outArgs_tyra);


      phiEvaluator_->setLinearizationPoint(inArgs);
      sStatus = phiEvaluator_->computePhi(vphi.ptr(), 1, dt, f);


      Thyra::V_VpStV(x.ptr(), *xOld, factor, *vphi);
    }

    std::cout << sStatus << std::endl;

    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperEPI3AppAction<Scalar>::ACTION_LOCATION::AFTER_EXP);

    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperEPI3AppAction<Scalar>::ACTION_LOCATION::END_STEP);
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
StepperEPI3<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperEPI3<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperEPI3 ---\n";
  out << "  stepperEEAppAction_                = "
      << stepperEEAppAction_ << std::endl;
  out << "----------------------------" << std::endl;
}


template<class Scalar>
bool StepperEPI3<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperImplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (stepperEEAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Backward Euler AppAction is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperEPI3<Scalar>::getValidParameters() const
{
  auto pl = this->getValidParametersBasicImplicit();
  return pl;
}


template <class Scalar>
void StepperEPI3<Scalar>::setStepperExponentialValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  Teuchos::RCP<Teuchos::ParameterList> phiPL = Teuchos::null;

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
Teuchos::RCP<StepperEPI3<Scalar> >
createStepperEPI3(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperEPI3<Scalar>());

  stepper->setStepperImplicitValues(pl);

  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}


} // namespace Tempus
#endif // Tempus_StepperEPI3_impl_hpp
