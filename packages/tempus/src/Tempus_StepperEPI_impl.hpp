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
  phiEvaluator_->checkInitialized();
//std::cout << "Inside takeStep 1." << std::endl;
  using Teuchos::RCP;

  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperEPI::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperEPI<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for EPI.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<StepperEPI<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);
//std::cout << "Inside takeStep 2." << std::endl;
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
//std::cout << "Inside takeStep 3." << std::endl;
    RCP<const Thyra::VectorBase<Scalar> > xOld = currentState->getX();
    RCP<Thyra::VectorBase<Scalar> > x = workingState->getX();
    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();
//std::cout << "Inside takeStep 4." << std::endl;
    const Scalar time = workingState->getTime();
    const Scalar dt = workingState->getTimeStep();
//std::cout << "Inside takeStep 5." << std::endl;
    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXP);
//std::cout << "Inside takeStep 6." << std::endl;
    //{
    //  Teuchos::basic_FancyOStream<char> ostr(Teuchos::rcp(&std::cout, false));
    //  this->describe(ostr, Teuchos::VERB_EXTREME);
    //}
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer;

    Thyra::SolveStatus<Scalar> sStatus;
//std::cout << "Inside takeStep 7." << std::endl;    
      // Setup TimeDerivative
      timeDer = Teuchos::rcp(new StepperEPITimeDerivative<Scalar>(
        Scalar(0.0),xOld));
      auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
          timeDer, dt, Scalar(0.0), Scalar(1.0)));
//std::cout << "Inside takeStep 8." << std::endl;

      // TODO: Transition away from using implicit solver methods and use ModelEvaluator directly
      RCP<Thyra::VectorBase<Scalar> > Mf = x->clone_v();
      // std::cout << "xO[0,1] = " << Thyra::get_ele(*xOld, 0) << " " << Thyra::get_ele(*xOld, 1) << std::endl;
      // std::cout << "x[0,1]  = " << Thyra::get_ele(*x, 0) << " " << Thyra::get_ele(*x, 1) << std::endl;
      // std::cout << "f[0,1]  = " << Thyra::get_ele(*Mf, 0) << " " << Thyra::get_ele(*Mf, 1) << std::endl;
      this->evaluateImplicitODE(Mf, x, xDot, time, p);
      // f = M*xDot in here
      // std::cout << "xO[0,1] = " << Thyra::get_ele(*xOld, 0) << " " << Thyra::get_ele(*xOld, 1) << std::endl;
      // std::cout << "x[0,1]  = " << Thyra::get_ele(*x, 0) << " " << Thyra::get_ele(*x, 1) << std::endl;
      // std::cout << "f[0,1]  = " << Thyra::get_ele(*f, 0) << " " << Thyra::get_ele(*f, 1) << std::endl;
//std::cout << "Inside takeStep 9." << std::endl;
      // Using the appModel
      RCP<const Thyra::ModelEvaluator<Scalar>> appModel = this->getModel();
      Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = appModel->createInArgs();
      // Model evaluator builds: alpha*u_dot + beta*F(u) = 0
      inArgs.set_x(x);
      inArgs.set_t(time);
      inArgs.set_x_dot(xDot); // for what? xDot is zero at this point, updating of it not decided
//std::cout << "Inside takeStep 10." << std::endl;
      // initialize space for the update
      RCP<Thyra::VectorBase<Scalar>> vphi = x->clone_v();
      assign(vphi.ptr(), ST::zero());  // Must initialize to a guess before solve!
//std::cout << "Inside takeStep 11." << std::endl;
      // auto thyraModel = this->getModel();
//std::cout << "Inside takeStep 12." << std::endl;
      // auto inArgs_tyra  = thyraModel->createInArgs();
      // auto outArgs_tyra = thyraModel->createOutArgs();
//std::cout << "Inside takeStep 13." << std::endl;
      // auto curr = solutionHistory->getCurrentState();

      // inArgs_tyra.set_x(curr->getX());
      // inArgs_tyra.set_t(curr->getTime());
//std::cout << "Inside takeStep 14." << std::endl;
      // auto W = thyraModel->create_W_op();
      // outArgs_tyra.set_W_op(W);
//std::cout << "Inside takeStep 15." << std::endl;
      // thyraModel->evalModel(inArgs_tyra, outArgs_tyra);

//std::cout << "Inside takeStep 16." << std::endl;
      phiEvaluator_->setLinearizationPoint(inArgs);
//std::cout << "Inside takeStep 16.5." << std::endl;
      // TODO: Avoid using hard coded EPI2 (p=2) and adjust the logic for general p.
      // TODO: right now, we have p=1 for EPI2, since we do not compute dF/dt yet
      // (correct only for autonomous problems)
      sStatus = phiEvaluator_->computePhi(vphi.ptr(), 1, dt, Mf);
//std::cout << "Inside takeStep 17." << std::endl;

      Thyra::V_VpStV(x.ptr(), *xOld, Scalar(-1.0)*dt, *vphi);
//std::cout << "Inside takeStep 18." << std::endl;

    // std::cout << sStatus << std::endl;

    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::AFTER_EXP);
//std::cout << "Inside takeStep 19." << std::endl;
    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);
//std::cout << "Inside takeStep 20." << std::endl;
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
