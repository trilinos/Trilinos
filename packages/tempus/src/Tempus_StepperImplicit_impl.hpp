// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperImplicit_impl_hpp
#define Tempus_StepperImplicit_impl_hpp

// Tempus
//#include "Tempus_Stepper.hpp"
//#include "Tempus_TimeDerivative.hpp"

// Thrya
//#include "Thyra_VectorBase.hpp"
//#include "Thyra_VectorStdOps.hpp"
#include "NOX_Thyra.H"


namespace Tempus {


template<class Scalar>
void StepperImplicit<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  validImplicitODE_DAE(appModel);
  wrapperModel_ =
    Teuchos::rcp(new WrapperModelEvaluatorBasic<Scalar>(appModel));
}


template<class Scalar>
void StepperImplicit<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}


template<class Scalar>
void StepperImplicit<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(numStates < 1, std::logic_error,
    "Error - setInitialConditions() needs at least one SolutionState\n"
    "        to set the initial condition.  Number of States = " << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperImplicit::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"<<std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x = initialState->getX();

  auto inArgs = this->wrapperModel_->getNominalValues();
  if (this->getOrderODE() == SECOND_ORDER_ODE) {
    RCP<Thyra::VectorBase<Scalar> > xDot = initialState->getXDot();

    // If initialState has x and xDot set, treat them as the initial conditions.
    // Otherwise use the x and xDot from getNominalValues() as the ICs.
    TEUCHOS_TEST_FOR_EXCEPTION(
      !((x != Teuchos::null && xDot != Teuchos::null) ||
        (inArgs.get_x() != Teuchos::null &&
         inArgs.get_x_dot() != Teuchos::null)), std::logic_error,
      "Error - We need to set the initial conditions for x and xDot from\n"
      "        either initialState or appModel_->getNominalValues::InArgs\n"
      "        (but not from a mixture of the two).\n");
  }

  if (this->getOrderODE() == FIRST_ORDER_ODE) {
    // Use x from inArgs as ICs.
    if (x == Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION( (x == Teuchos::null) &&
        (inArgs.get_x() == Teuchos::null), std::logic_error,
        "Error - setInitialConditions needs the ICs from the SolutionHistory\n"
        "        or getNominalValues()!\n");

      x = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
      initialState->setX(x);
    }
  }
  else if (this->getOrderODE() == SECOND_ORDER_ODE) {
    // Use the x and xDot from getNominalValues() as the ICs.
    using Teuchos::rcp_const_cast;
    RCP<Thyra::VectorBase<Scalar> > xDot = initialState->getXDot();
    if ( x == Teuchos::null || xDot == Teuchos::null ) {
      TEUCHOS_TEST_FOR_EXCEPTION( (inArgs.get_x() == Teuchos::null) ||
        (inArgs.get_x_dot() == Teuchos::null), std::logic_error,
        "Error - setInitialConditions() needs the ICs from the initialState\n"
        "        or getNominalValues()!\n");
      x = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
      initialState->setX(x);
      xDot = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x_dot());
      initialState->setXDot(xDot);
    }
  }


  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      if (initialState->getXDot() == Teuchos::null) {
        RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,
          "StepperImplicit::setInitialConditions()");
        *out << "Warning -- Requested IC consistency of 'None' but\n"
             << "           initialState does not have an xDot.\n"
             << "           Setting a 'Zero' xDot!\n" << std::endl;

        Thyra::assign(this->getStepperXDot(initialState).ptr(), Scalar(0.0));
      }
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {
      if (initialState->getXDotDot() == Teuchos::null) {
        RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,
          "StepperImplicit::setInitialConditions()");
        *out << "Warning -- Requested IC consistency of 'None' but\n"
             << "           initialState does not have an xDotDot.\n"
             << "           Setting a 'Zero' xDotDot!\n" << std::endl;

        Thyra::assign(this->getStepperXDotDot(initialState).ptr(), Scalar(0.0));
      }
    }
  }
  else if (icConsistency == "Zero") {
    if (this->getOrderODE() == FIRST_ORDER_ODE)
      Thyra::assign(this->getStepperXDot(initialState).ptr(), Scalar(0.0));
    else if (this->getOrderODE() == SECOND_ORDER_ODE)
      Thyra::assign(this->getStepperXDotDot(initialState).ptr(), Scalar(0.0));
  }
  else if (icConsistency == "App") {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      auto xDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
                    inArgs.get_x_dot());
      TEUCHOS_TEST_FOR_EXCEPTION(xDot == Teuchos::null, std::logic_error,
        "Error - setInitialConditions() requested 'App' for IC consistency,\n"
        "        but 'App' returned a null pointer for xDot!\n");
      Thyra::assign(this->getStepperXDot(initialState).ptr(), *xDot);
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {
      auto xDotDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
                       inArgs.get_x_dot_dot());
      TEUCHOS_TEST_FOR_EXCEPTION(xDotDot == Teuchos::null, std::logic_error,
        "Error - setInitialConditions() requested 'App' for IC consistency,\n"
        "        but 'App' returned a null pointer for xDotDot!\n");
      Thyra::assign(getStepperXDotDot(initialState).ptr(), *xDotDot);
    }
  }
  else if (icConsistency == "Consistent") {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      // Solve f(x, xDot,t) = 0.
      const Scalar time = initialState->getTime();
      const Scalar dt   = initialState->getTimeStep();
      RCP<TimeDerivative<Scalar> > timeDer = Teuchos::null;
      const Scalar alpha = Scalar(1.0);    // d(xDot)/d(xDot)
      const Scalar beta  = Scalar(0.0);    // d(x   )/d(xDot)
      auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
        timeDer, dt, alpha, beta, SOLVE_FOR_XDOT_CONST_X));

      auto xDot = this->getStepperXDot(initialState);
      const Thyra::SolveStatus<Scalar> sStatus =
        this->solveImplicitODE(x, xDot, time, p);

      TEUCHOS_TEST_FOR_EXCEPTION(
        sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED, std::logic_error,
        "Error - Solver failed while determining the initial conditions.\n"
        "        Solver status is "<<Thyra::toString(sStatus.solveStatus)<<".\n");
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {

      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - setInitialConditions(): 'Consistent' for second-order ODE\n"
        "        has not been implemented.\n");
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Error - setInitialConditions() invalid IC consistency, "
      << icConsistency << ".\n");
  }

  // At this point, x and xDot are sync'ed or consistent
  // at the same time level for the initialState.
  initialState->setIsSynced(true);

  // Test for consistency.
  if (this->getICConsistencyCheck()) {
    auto f    = initialState->getX()->clone_v();
    auto xDot = this->getStepperXDot(initialState);

    const Scalar time = initialState->getTime();
    const Scalar dt   = initialState->getTimeStep();
    RCP<TimeDerivative<Scalar> > timeDer = Teuchos::null;
    const Scalar alpha = Scalar(0.0);
    const Scalar beta  = Scalar(0.0);
    auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
      timeDer, dt, alpha, beta, EVALUATE_RESIDUAL));

    this->evaluateImplicitODE(f, x, xDot, time, p);

    Scalar normX = Thyra::norm(*x);
    Scalar reldiff = Scalar(0.0);
    if (normX == Scalar(0.0)) reldiff = Thyra::norm(*f);
    else reldiff = Thyra::norm(*f)/normX;

    Scalar eps = Scalar(100.0)*std::abs(Teuchos::ScalarTraits<Scalar>::eps());
    if (reldiff > eps) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"StepperImplicit::setInitialConditions()");
      *out << "Warning -- Failed consistency check but continuing!\n"
         << "  ||f(x,xDot,t)||/||x|| > eps" << std::endl
         << "  ||f(x,xDot,t)||       = " << Thyra::norm(*f) << std::endl
         << "  ||x||                 = " << Thyra::norm(*x) << std::endl
         << "  ||f(x,xDot,t)||/||x|| = " << reldiff         << std::endl
         << "                    eps = " << eps             << std::endl;
    }
  }
}


template<class Scalar>
void StepperImplicit<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  if (solver == Teuchos::null) {
    solver = rcp(new Thyra::NOXNonlinearSolver());
    solver->setParameterList(defaultSolverParameters());
  }

  solver_ = solver;

  TEUCHOS_TEST_FOR_EXCEPTION(wrapperModel_ == Teuchos::null, std::logic_error,
       "Error - ModelEvaluator is unset!  Should call setModel() first.\n");

  solver_->setModel(wrapperModel_);

}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperImplicit<Scalar>::
getStepperXDot(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getXDot() != Teuchos::null) stepperXDot_ = state->getXDot();
  // Else use temporary storage stepperXDot_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperXDot_ == Teuchos::null, std::logic_error,
    "Error - stepperXDot_ has not been set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperXDot_;
}


template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperImplicit<Scalar>::
getStepperXDotDot(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getXDotDot() != Teuchos::null) stepperXDotDot_=state->getXDotDot();
  // Else use temporary storage stepperXDotDot_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperXDotDot_ == Teuchos::null,std::logic_error,
    "Error - stepperXDotDot_ has not been set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperXDotDot_;
}


template<class Scalar>
const Thyra::SolveStatus<Scalar>
StepperImplicit<Scalar>::solveImplicitODE(
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x)
{
  if (getZeroInitialGuess())
    Thyra::assign(x.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

  const Thyra::SolveStatus<Scalar> sStatus = (*solver_).solve(&*x);

  return sStatus;
}


template<class Scalar>
const Thyra::SolveStatus<Scalar>
StepperImplicit<Scalar>::solveImplicitODE(
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & xDot,
  const Scalar time,
  const Teuchos::RCP<ImplicitODEParameters<Scalar> > & p )
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar>  inArgs  = wrapperModel_->getInArgs();
  MEB::OutArgs<Scalar> outArgs = wrapperModel_->getOutArgs();
  inArgs.set_x(x);
  if (inArgs.supports(MEB::IN_ARG_x_dot    )) inArgs.set_x_dot    (xDot);
  if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (time);
  if (inArgs.supports(MEB::IN_ARG_step_size))
    inArgs.set_step_size(p->timeStepSize_);
  if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (p->alpha_);
  if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (p->beta_);
  if (inArgs.supports(MEB::IN_ARG_stage_number))
    inArgs.set_stage_number(p->stageNumber_);

  wrapperModel_->setForSolve(p->timeDer_, inArgs, outArgs, p->evaluationType_);

  Thyra::SolveStatus<Scalar> sStatus;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  switch (p->evaluationType_)
  {
    case SOLVE_FOR_X: {
      if (getZeroInitialGuess()) Thyra::assign(x.ptr(), ST::zero());
      sStatus = (*solver_).solve(&*x);
      break;
    }
    case SOLVE_FOR_XDOT_CONST_X: {
      //if (getZeroInitialGuess()) Thyra::assign(xDot.ptr(), ST::zero());
      sStatus = (*solver_).solve(&*xDot);
      break;
    }
    default: {
      TEUCHOS_TEST_FOR_EXCEPT("Invalid EVALUATION_TYPE!");
    }
  }

  return sStatus;
}


template<class Scalar>
void
StepperImplicit<Scalar>::evaluateImplicitODE(
        Teuchos::RCP<Thyra::VectorBase<Scalar> > & f,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & xDot,
  const Scalar time,
  const Teuchos::RCP<ImplicitODEParameters<Scalar> > & p )
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar>  inArgs  = wrapperModel_->getInArgs();
  inArgs.set_x(x);
  if (inArgs.supports(MEB::IN_ARG_x_dot    )) inArgs.set_x_dot    (xDot);
  if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (time);
  if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(p->timeStepSize_);
  if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (Scalar(0.0));
  if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (Scalar(0.0));

  MEB::OutArgs<Scalar> outArgs = wrapperModel_->getOutArgs();
  outArgs.set_f(f);

  wrapperModel_->setForSolve(Teuchos::null,inArgs,outArgs,p->evaluationType_);

  wrapperModel_->evalModel(inArgs, outArgs);
}


} // namespace Tempus
#endif // Tempus_StepperImplicit_impl_hpp
