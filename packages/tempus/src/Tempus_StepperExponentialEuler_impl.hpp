// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExponentialEuler_impl_hpp
#define Tempus_StepperExponentialEuler_impl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"
#include "Tempus_StepperExponentialEulerModifierDefault.hpp"
#include "Tempus_StepperExponentialEuler_decl.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Tempus_StepperFactory.hpp"

// TODO: have to include this header to get LSP to work.
#include "Tempus_StepperExponentialEuler.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_OperatorVectorTypes.hpp"

namespace Tempus {


template<class Scalar>
StepperExponentialEuler<Scalar>::StepperExponentialEuler()
{
  this->setStepperName("Exponential Euler");
  this->setStepperType("Exponential Euler");
  this->setUseFSAL(false);
  // TODO: think about the right setting here
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);
  this->setZeroInitialGuess(false);

  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}


template<class Scalar>
StepperExponentialEuler<Scalar>::StepperExponentialEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess,
  const Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> >& stepperEEAppAction)
{
  this->setStepperName("Exponential Euler");
  this->setStepperType("Exponential Euler");
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
void StepperExponentialEuler<Scalar>::setAppAction(
  Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperEEAppAction_ =
      Teuchos::rcp(new StepperExponentialEulerModifierDefault<Scalar>());
  } else {
    stepperEEAppAction_ = appAction;
  }
  this->isInitialized_ = false;
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  StepperImplicit<Scalar>::setModel(appModel);

  phiEvaluator_->setModel(appModel);
  phiEvaluator_->initialize();

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::setInitialConditions(
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
void StepperExponentialEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();
  phiEvaluator_->checkInitialized();

  using Teuchos::RCP;

  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperExponentialEuler::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperExponentialEuler<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for Exponential Euler.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<StepperExponentialEuler<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

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
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXP);

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
      timeDer = Teuchos::rcp(new StepperExponentialEulerTimeDerivative<Scalar>(
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
      timeDer = Teuchos::rcp(new StepperExponentialEulerTimeDerivative<Scalar>(
        Scalar(0.0),xOld));
      auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
          timeDer, dt, Scalar(0.0), Scalar(1.0)));

      //std::cout << "x[0,1]  = " << Thyra::get_ele(*x, 0) << " " << Thyra::get_ele(*x, 1) << std::endl;

      // TODO: Transition away from using implicit solver methods and use ModelEvaluator directly
      RCP<Thyra::VectorBase<Scalar> > f = x->clone_v();
      this->evaluateImplicitODE(f, x, xDot, time, p);

      //std::cout << "xO[0,1] = " << Thyra::get_ele(*xOld, 0) << " " << Thyra::get_ele(*xOld, 1) << std::endl;
      //std::cout << "x[0,1]  = " << Thyra::get_ele(*x, 0) << " " << Thyra::get_ele(*x, 1) << std::endl;
      //std::cout << "f[0,1]  = " << Thyra::get_ele(*f, 0) << " " << Thyra::get_ele(*f, 1) << std::endl;

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
      if (use_phi_eval) {
	// use the PhiEvaluator to compute update
	//Teuchos::RCP<Teuchos::FancyOStream> out =
	//  Teuchos::VerboseObjectBase::getDefaultOStream();
        //out->setOutputToRootOnly(0);
	//phiEvaluator_->describe(*out, Teuchos::VERB_EXTREME);

	phiEvaluator_->setLinearizationPoint(inArgs);
        sStatus = phiEvaluator_->computePhi(vphi.ptr(), 1, dt, f);
      }
      else {
	const Scalar alpha = Scalar(1.0)/dt;
	const Scalar beta  = Scalar(0.5);
	inArgs.set_alpha(alpha);
	inArgs.set_beta(beta);

        Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs = appModel->createOutArgs();

	// first allocate space for the jacobian
	RCP<Thyra::LinearOpBase<Scalar>> jac = appModel->create_W_op();
	RCP<Thyra::PreconditionerBase<Scalar>> jac_p = Teuchos::null;
	if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
	  jac_p = appModel->create_W_prec();
	}
	// set only the jacobian matrix
	outArgs.set_W_op(jac);
	if (jac_p != Teuchos::null){
	  outArgs.set_W_prec(jac_p);
	}

	// this will fill the Jacobian operator
	appModel->evalModel(inArgs, outArgs);

	// TODO: const-cast why?
	RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar>> const_lowsFactory = appModel->get_W_factory();
	RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar>> lowsFactory =
          Teuchos::rcp_const_cast<Thyra::LinearOpWithSolveFactoryBase<Scalar>>(const_lowsFactory);

	RCP<Thyra::LinearOpWithSolveBase<Scalar>> LOWSB = Teuchos::null;
	if (jac_p == Teuchos::null){
	  // without preconditioner
	  RCP<const Thyra::LinearOpBase<Scalar>> const_jac = Teuchos::rcpFromRef(*jac);
	  LOWSB = Thyra::linearOpWithSolve(*lowsFactory, const_jac);
	}
	else {
	  // with preconditioner
	  LOWSB = lowsFactory->createOp();
	  Thyra::initializePreconditionedOp<Scalar>(*lowsFactory, jac, jac_p, LOWSB.ptr());
	}

	// compute an approximation to dt*phi_1(dt*J)*f and write it to x
	sStatus = LOWSB->solve(Thyra::NOTRANS, *f, vphi.ptr());
	factor = Scalar(-1.);
      }

      //std::cout << "ph[0,1] = " << Thyra::get_ele(*x, 0) << " " << Thyra::get_ele(*x, 1) << std::endl;
      //assign(x.ptr(), ST::zero());

      // x = xOld - dt*phi_1(dt*J)*f
      Thyra::V_VpStV(x.ptr(), *xOld, factor, *vphi);
    }

    //std::cout << sStatus << std::endl;

    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::AFTER_EXP);

    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::END_STEP);
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
StepperExponentialEuler<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperExponentialEuler ---\n";
  out << "  stepperEEAppAction_                = "
      << stepperEEAppAction_ << std::endl;
  out << "----------------------------" << std::endl;
}


template<class Scalar>
bool StepperExponentialEuler<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
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
StepperExponentialEuler<Scalar>::getValidParameters() const
{
  auto pl = this->getValidParametersBasicImplicit();
  return pl;
}


template <class Scalar>
void StepperExponentialEuler<Scalar>::setStepperExponentialValues(
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
Teuchos::RCP<StepperExponentialEuler<Scalar> >
createStepperExponentialEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperExponentialEuler<Scalar>());

  stepper->setStepperImplicitValues(pl);

  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}


} // namespace Tempus
#endif // Tempus_StepperExponentialEuler_impl_hpp
