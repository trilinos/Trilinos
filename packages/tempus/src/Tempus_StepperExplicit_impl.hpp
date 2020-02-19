// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicit_impl_hpp
#define Tempus_StepperExplicit_impl_hpp


namespace Tempus {


template<class Scalar>
void StepperExplicit<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  validExplicitODE(appModel);
  appModel_ = appModel;

  inArgs_  = appModel_->getNominalValues();
  outArgs_ = appModel_->createOutArgs();
}


template<class Scalar>
void StepperExplicit<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(numStates < 1, std::logic_error,
    "Error - setInitialConditions() needs at least one SolutionState\n"
    "        to set the initial condition.  Number of States = " << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1, "StepperExplicit::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"<<std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x = initialState->getX();

  inArgs_ = appModel_->getNominalValues();
  if (this->getOrderODE() == SECOND_ORDER_ODE) {
    RCP<Thyra::VectorBase<Scalar> > xDot = initialState->getXDot();

    // If initialState has x and xDot set, treat them as the initial conditions.
    // Otherwise use the x and xDot from getNominalValues() as the ICs.
    TEUCHOS_TEST_FOR_EXCEPTION(
      !((x != Teuchos::null && xDot != Teuchos::null) ||
        (inArgs_.get_x() != Teuchos::null &&
         inArgs_.get_x_dot() != Teuchos::null)), std::logic_error,
      "Error - We need to set the initial conditions for x and xDot from\n"
      "        either initialState or appModel_->getNominalValues::InArgs\n"
      "        (but not from a mixture of the two).\n");
  }

  if (this->getOrderODE() == FIRST_ORDER_ODE) {
    if (x == Teuchos::null) {
      // Use x from inArgs as ICs.
      TEUCHOS_TEST_FOR_EXCEPTION( (x == Teuchos::null) &&
        (inArgs_.get_x() == Teuchos::null), std::logic_error,
        "Error - setInitialConditions needs the ICs from the SolutionHistory\n"
        "        or getNominalValues()!\n");

      x = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs_.get_x());
      initialState->setX(x);
    }
  }
  else if (this->getOrderODE() == SECOND_ORDER_ODE) {
    using Teuchos::rcp_const_cast;
    // Use the x and xDot from getNominalValues() as the ICs.
    if ( initialState->getX() == Teuchos::null ||
         initialState->getXDot() == Teuchos::null ) {
      TEUCHOS_TEST_FOR_EXCEPTION( (inArgs_.get_x() == Teuchos::null) ||
        (inArgs_.get_x_dot() == Teuchos::null), std::logic_error,
        "Error - setInitialConditions() needs the ICs from the initialState\n"
        "        or getNominalValues()!\n");
      x = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs_.get_x());
      initialState->setX(x);
      RCP<Thyra::VectorBase<Scalar> > xDot
        = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs_.get_x_dot());
      initialState->setXDot(xDot);
    }
  }

  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      if (initialState->getXDot() == Teuchos::null) {
        RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"StepperExplicit::setInitialConditions()");
        *out << "Warning -- Requested IC consistency of 'None' but\n"
             << "           initialState does not have an xDot.\n"
             << "           Setting a 'Consistent' xDot!\n" << std::endl;
        auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(0.0));
        evaluateExplicitODE(getStepperXDot(initialState), x,
                            initialState->getTime(), p);
        initialState->setIsSynced(true);
      }
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {
      if (initialState->getXDotDot() == Teuchos::null) {
        RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"StepperExplicit::setInitialConditions()");
        *out << "Warning -- Requested IC consistency of 'None' but\n"
             << "           initialState does not have an xDotDot.\n"
             << "           Setting a 'Consistent' xDotDot!\n" << std::endl;
        auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(0.0));
        this->evaluateExplicitODE(getStepperXDotDot(initialState), x,
                                  initialState->getXDot(),
                                  initialState->getTime(), p);
        initialState->setIsSynced(true);
      }
    }
  }
  else if (icConsistency == "Zero") {
    if (this->getOrderODE() == FIRST_ORDER_ODE)
      Thyra::assign(getStepperXDot(initialState).ptr(), Scalar(0.0));
    else if (this->getOrderODE() == SECOND_ORDER_ODE)
      Thyra::assign(getStepperXDotDot(initialState).ptr(), Scalar(0.0));
  }
  else if (icConsistency == "App") {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      auto xDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
                    inArgs_.get_x_dot());
      TEUCHOS_TEST_FOR_EXCEPTION(xDot == Teuchos::null, std::logic_error,
        "Error - setInitialConditions() requested 'App' for IC consistency,\n"
        "        but 'App' returned a null pointer for xDot!\n");
      Thyra::assign(getStepperXDot(initialState).ptr(), *xDot);
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {
      auto xDotDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
                       inArgs_.get_x_dot_dot());
      TEUCHOS_TEST_FOR_EXCEPTION(xDotDot == Teuchos::null, std::logic_error,
        "Error - setInitialConditions() requested 'App' for IC consistency,\n"
        "        but 'App' returned a null pointer for xDotDot!\n");
      Thyra::assign(getStepperXDotDot(initialState).ptr(), *xDotDot);
    }
  }
  else if (icConsistency == "Consistent") {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      // Evaluate xDot = f(x,t).
      auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(0.0));
      evaluateExplicitODE(getStepperXDot(initialState), x,
                          initialState->getTime(), p);
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {
      // Evaluate xDotDot = f(x,xDot,t).
      auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(0.0));
      this->evaluateExplicitODE(initialState->getXDotDot(), x,
                                initialState->getXDot(),
                                initialState->getTime(), p);
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Error - setInitialConditions() invalid IC consistency, '"
      << icConsistency << "'.\n");
  }

  // At this point, x, and xDot (and xDotDot) sync'ed or consistent
  // at the same time level for the initialState.
  initialState->setIsSynced(true);

  // Test for consistency.
  if (this->getICConsistencyCheck()) {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      auto xDot = getStepperXDot(initialState);
      auto f    = initialState->getX()->clone_v();
      auto p    = Teuchos::rcp(new ExplicitODEParameters<Scalar>(0.0));
      evaluateExplicitODE(f, x, initialState->getTime(), p);
      Thyra::Vp_StV(f.ptr(), Scalar(-1.0), *(xDot));
      Scalar normX = Thyra::norm(*x);
      Scalar reldiff = Scalar(0.0);
      if (normX == Scalar(0.0)) reldiff = Thyra::norm(*f);
      else reldiff = Thyra::norm(*f)/normX;

      Scalar eps = Scalar(100.0)*std::abs(Teuchos::ScalarTraits<Scalar>::eps());
      if (reldiff > eps) {
        RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"StepperExplicit::setInitialConditions()");
        *out << "Warning -- Failed consistency check but continuing!\n"
           << "  ||xDot-f(x,t)||/||x|| > eps" << std::endl
           << "  ||xDot-f(x,t)||       = " << Thyra::norm(*f) << std::endl
           << "  ||x||                 = " << Thyra::norm(*x) << std::endl
           << "  ||xDot-f(x,t)||/||x|| = " << reldiff         << std::endl
           << "                    eps = " << eps             << std::endl;
      }
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {
      auto xDotDot = initialState->getXDotDot();
      auto f       = initialState->getX()->clone_v();
      auto p       = Teuchos::rcp(new ExplicitODEParameters<Scalar>(0.0));
      this->evaluateExplicitODE(f, x, initialState->getXDot(),
                                initialState->getTime(), p);
      Thyra::Vp_StV(f.ptr(), Scalar(-1.0), *(xDotDot));
      Scalar normX = Thyra::norm(*x);
      Scalar reldiff = Scalar(0.0);
      if (normX == Scalar(0.0)) reldiff = Thyra::norm(*f);
      else reldiff = Thyra::norm(*f)/normX;

      Scalar eps = Scalar(100.0)*std::abs(Teuchos::ScalarTraits<Scalar>::eps());
      if (reldiff > eps) {
        RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"StepperExplicit::setInitialConditions()");
        *out << "Warning -- Failed consistency check but continuing!\n"
           << "  ||xDotDot-f(x,xDot,t)||/||x|| > eps" << std::endl
           << "  ||xDotDot-f(x,xDot,t)||       = " << Thyra::norm(*f)
           << std::endl
           << "  ||x||                         = " << Thyra::norm(*x)
           << std::endl
           << "  ||xDotDot-f(x,xDot,t)||/||x|| = " << reldiff << std::endl
           << "                            eps = " << eps     << std::endl;
      }
    }
  }
}

template<class Scalar>
void StepperExplicit<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > /* solver */)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperExplicit::setSolver()");
  *out << "Warning -- No solver to set for StepperExplicit "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperExplicit<Scalar>::
getStepperX(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getX() != Teuchos::null) stepperX_ = state->getX();
  // Else use temporary storage stepperXp_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperX_ == Teuchos::null, std::logic_error,
    "Error - stepperX_ has not been set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperX_;
}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperExplicit<Scalar>::
getStepperXDot(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getXDot() != Teuchos::null) stepperXDot_ = state->getXDot();
  // Else use temporary storage stepperXDot_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperXDot_ == Teuchos::null, std::logic_error,
    "Error - stepperXDot_ has not set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperXDot_;
}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperExplicit<Scalar>::
getStepperXDotDot(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getXDotDot() != Teuchos::null) stepperXDotDot_=state->getXDotDot();
  // Else use temporary storage stepperXDotDot_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperXDotDot_==Teuchos::null, std::logic_error,
    "Error - stepperXDotDot_ has not set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperXDotDot_;
}

template<class Scalar>
void
StepperExplicit<Scalar>::
evaluateExplicitODE(Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
                    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
                    const Scalar time,
                    const Teuchos::RCP<ExplicitODEParameters<Scalar> > & p )
{
  typedef Thyra::ModelEvaluatorBase MEB;

  inArgs_.set_x(x);
  if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time);
  if (inArgs_.supports(MEB::IN_ARG_step_size))
    inArgs_.set_step_size(p->timeStepSize_);
  if (inArgs_.supports(MEB::IN_ARG_stage_number))
    inArgs_.set_stage_number(p->stageNumber_);

  // For model evaluators whose state function f(x, xDot, t) describes
  // an implicit ODE, and which accept an optional xDot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // xDot = f(x, t).
  if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);

  outArgs_.set_f(xDot);

  appModel_->evalModel(inArgs_, outArgs_);
}

template<class Scalar>
void
StepperExplicit<Scalar>::
evaluateExplicitODE(Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot,
                    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
                    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDot,
                    const Scalar time,
                    const Teuchos::RCP<ExplicitODEParameters<Scalar> > & p )
{
  typedef Thyra::ModelEvaluatorBase MEB;

  inArgs_.set_x(x);
  if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(xDot);
  if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time);
  if (inArgs_.supports(MEB::IN_ARG_step_size))
    inArgs_.set_step_size(p->timeStepSize_);
  if (inArgs_.supports(MEB::IN_ARG_stage_number))
    inArgs_.set_stage_number(p->stageNumber_);

  // For model evaluators whose state function f(x, xDot, xDotDot, t) describes
  // an implicit ODE, and which accept an optional xDotDot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // xDotDot = f(x, xDot, t).
  if (inArgs_.supports(MEB::IN_ARG_x_dot_dot))
    inArgs_.set_x_dot_dot(Teuchos::null);

  outArgs_.set_f(xDotDot);

  appModel_->evalModel(inArgs_, outArgs_);
}



template<class Scalar>
void StepperExplicit<Scalar>::describe(Teuchos::FancyOStream        & out,
                               const Teuchos::EVerbosityLevel verbLevel) const
{
  out << "--- StepperExplicit ---\n";
  out << "  appModel_         = " << appModel_ << std::endl;
  out << "  inArgs_           = " << inArgs_ << std::endl;
  out << "  outArgs_          = " << outArgs_ << std::endl;
  out << "  stepperX_         = " << stepperX_ << std::endl;
  out << "  stepperXDot_      = " << stepperXDot_ << std::endl;
  out << "  stepperXDotDot_   = " << stepperXDotDot_ << std::endl;
}


template<class Scalar>
bool StepperExplicit<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if (appModel_ == Teuchos::null) {
    isValidSetup = false;
    out << "The application ModelEvaluator is not set!\n";
  }

  return isValidSetup;
}


} // namespace Tempus
#endif // Tempus_StepperExplicit_impl_hpp
