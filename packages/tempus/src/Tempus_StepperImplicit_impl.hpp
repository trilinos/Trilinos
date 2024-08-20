//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperImplicit_impl_hpp
#define Tempus_StepperImplicit_impl_hpp

// Thrya
#include "NOX_Thyra.H"

namespace Tempus {

template <class Scalar>
void StepperImplicit<Scalar>::setModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  validImplicitODE_DAE(appModel);
  wrapperModel_ =
      Teuchos::rcp(new WrapperModelEvaluatorBasic<Scalar>(appModel));

  TEUCHOS_TEST_FOR_EXCEPTION(solver_ == Teuchos::null, std::logic_error,
                             "Error - Solver is not set!\n");
  solver_->setModel(wrapperModel_);

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperImplicit<Scalar>::setInitialConditions(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(
      numStates < 1, std::logic_error,
      "Error - setInitialConditions() needs at least one SolutionState\n"
      "        to set the initial condition.  Number of States = "
          << numStates);

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x        = initialState->getX();
  RCP<Thyra::VectorBase<Scalar> > xDot     = initialState->getXDot();
  if (xDot() == Teuchos::null) xDot = this->getStepperXDot();

  auto inArgs = this->wrapperModel_->getNominalValues();
  if (this->getOrderODE() == SECOND_ORDER_ODE) {
    RCP<Thyra::VectorBase<Scalar> > initialXDot = initialState->getXDot();
    // If initialState has x and xDot set, treat them as the initial conditions.
    // Otherwise use the x and xDot from getNominalValues() as the ICs.
    TEUCHOS_TEST_FOR_EXCEPTION(
        !((x != Teuchos::null && initialXDot != Teuchos::null) ||
          (inArgs.get_x() != Teuchos::null &&
           inArgs.get_x_dot() != Teuchos::null)),
        std::logic_error,
        "Error - We need to set the initial conditions for x and xDot from\n"
        "        either initialState or appModel_->getNominalValues::InArgs\n"
        "        (but not from a mixture of the two).\n");
  }

  if (this->getOrderODE() == FIRST_ORDER_ODE) {
    // Use x from inArgs as ICs.
    if (x == Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          (x == Teuchos::null) && (inArgs.get_x() == Teuchos::null),
          std::logic_error,
          "Error - setInitialConditions needs the ICs from the "
          "SolutionHistory\n"
          "        or getNominalValues()!\n");

      x = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
      initialState->setX(x);
    }
  }
  else if (this->getOrderODE() == SECOND_ORDER_ODE) {
    // Use the x and xDot from getNominalValues() as the ICs.
    using Teuchos::rcp_const_cast;
    if (x == Teuchos::null || xDot == Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          (inArgs.get_x() == Teuchos::null) ||
              (inArgs.get_x_dot() == Teuchos::null),
          std::logic_error,
          "Error - setInitialConditions() needs the ICs from the initialState\n"
          "        or getNominalValues()!\n");
      x = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
      initialState->setX(x);
      RCP<Thyra::VectorBase<Scalar> > x_dot =
          rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x_dot());
      initialState->setXDot(x_dot);
    }
  }

  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      if (initialState->getXDot() == Teuchos::null)
        Thyra::assign(xDot.ptr(), Scalar(0.0));
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {
      if (initialState->getXDotDot() == Teuchos::null)
        Thyra::assign(this->getStepperXDotDot(initialState).ptr(), Scalar(0.0));
    }
  }
  else if (icConsistency == "Zero") {
    if (this->getOrderODE() == FIRST_ORDER_ODE)
      Thyra::assign(xDot.ptr(), Scalar(0.0));
    else if (this->getOrderODE() == SECOND_ORDER_ODE)
      Thyra::assign(this->getStepperXDotDot(initialState).ptr(), Scalar(0.0));
  }
  else if (icConsistency == "App") {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      auto x_dot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
          inArgs.get_x_dot());
      TEUCHOS_TEST_FOR_EXCEPTION(
          x_dot == Teuchos::null, std::logic_error,
          "Error - setInitialConditions() requested 'App' for IC consistency,\n"
          "        but 'App' returned a null pointer for xDot!\n");
      Thyra::assign(xDot.ptr(), *x_dot);
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {
      auto xDotDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
          inArgs.get_x_dot_dot());
      TEUCHOS_TEST_FOR_EXCEPTION(
          xDotDot == Teuchos::null, std::logic_error,
          "Error - setInitialConditions() requested 'App' for IC consistency,\n"
          "        but 'App' returned a null pointer for xDotDot!\n");
      Thyra::assign(this->getStepperXDotDot(initialState).ptr(), *xDotDot);
    }
  }
  else if (icConsistency == "Consistent") {
    if (this->getOrderODE() == FIRST_ORDER_ODE) {
      // Solve f(x, xDot,t) = 0.
      const Scalar time                    = initialState->getTime();
      const Scalar dt                      = initialState->getTimeStep();
      RCP<TimeDerivative<Scalar> > timeDer = Teuchos::null;
      const Scalar alpha                   = Scalar(1.0);  // d(xDot)/d(xDot)
      const Scalar beta                    = Scalar(0.0);  // d(x   )/d(xDot)
      auto p                               = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
          timeDer, dt, alpha, beta, SOLVE_FOR_XDOT_CONST_X));

      const Thyra::SolveStatus<Scalar> sStatus =
          this->solveImplicitODE(x, xDot, time, p);

      TEUCHOS_TEST_FOR_EXCEPTION(
          sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED,
          std::logic_error,
          "Error - Solver failed while determining the initial conditions.\n"
          "        Solver status is "
              << Thyra::toString(sStatus.solveStatus) << ".\n");
    }
    else if (this->getOrderODE() == SECOND_ORDER_ODE) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Error - setInitialConditions(): 'Consistent' for second-order ODE\n"
          "        has not been implemented.\n");
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error - setInitialConditions() invalid IC consistency, "
            << icConsistency << ".\n");
  }

  // At this point, x and xDot are sync'ed or consistent
  // at the same time level for the initialState.
  initialState->setIsSynced(true);

  // Test for consistency.
  if (this->getICConsistencyCheck()) {
    auto f = initialState->getX()->clone_v();

    const Scalar time                    = initialState->getTime();
    const Scalar dt                      = initialState->getTimeStep();
    RCP<TimeDerivative<Scalar> > timeDer = Teuchos::null;
    const Scalar alpha                   = Scalar(0.0);
    const Scalar beta                    = Scalar(0.0);
    auto p                               = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
        timeDer, dt, alpha, beta, EVALUATE_RESIDUAL));

    this->evaluateImplicitODE(f, x, xDot, time, p);

    Scalar normX   = Thyra::norm(*x);
    Scalar reldiff = Scalar(0.0);
    if (normX == Scalar(0.0))
      reldiff = Thyra::norm(*f);
    else
      reldiff = Thyra::norm(*f) / normX;

    Scalar eps                     = Scalar(100.0) * std::abs(Teuchos::ScalarTraits<Scalar>::eps());
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out, 1, "StepperImplicit::setInitialConditions()");
    if (reldiff < eps) {
      *out << "\n---------------------------------------------------\n"
           << "Info -- Stepper = " << this->getStepperType() << "\n"
           << "  Initial condition PASSED consistency check!\n"
           << "  (||f(x,xDot,t)||/||x|| = " << reldiff << ") < "
           << "(eps = " << eps << ")" << std::endl
           << "---------------------------------------------------\n"
           << std::endl;
    }
    else {
      *out << "\n---------------------------------------------------\n"
           << "Info -- Stepper = " << this->getStepperType() << "\n"
           << "  Initial condition FAILED consistency check but continuing!\n"
           << "  (||f(x,xDot,t)||/||x|| = " << reldiff << ") > "
           << "(eps = " << eps << ")" << std::endl
           << "  ||f(x,xDot,t)|| = " << Thyra::norm(*f) << std::endl
           << "  ||x||           = " << Thyra::norm(*x) << std::endl
           << "---------------------------------------------------\n"
           << std::endl;
    }
  }
}

template <class Scalar>
void StepperImplicit<Scalar>::setDefaultSolver()
{
  solver_       = rcp(new Thyra::NOXNonlinearSolver());
  auto solverPL = defaultSolverParameters();
  this->setSolverName("Default Solver");
  auto subPL = sublist(solverPL, "NOX");
  solver_->setParameterList(subPL);

  if (wrapperModel_ != Teuchos::null) solver_->setModel(wrapperModel_);

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperImplicit<Scalar>::setSolver(
    Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      solver == Teuchos::null, std::logic_error,
      "Error - Can not set the solver to Teuchos::null.\n");

  solver_ = solver;
  if (wrapperModel_ != Teuchos::null) solver_->setModel(wrapperModel_);

  this->isInitialized_ = false;
}

template <class Scalar>
const Thyra::SolveStatus<Scalar> StepperImplicit<Scalar>::solveImplicitODE(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
    const Teuchos::RCP<ImplicitODEParameters<Scalar> >& p,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& y, const int index)
{
  wrapperModel_->setForSolve(x, xDot, time, p, y, index);

  Thyra::SolveStatus<Scalar> sStatus;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  switch (p->evaluationType_) {
    case SOLVE_FOR_X: {
      if (getZeroInitialGuess()) Thyra::assign(x.ptr(), ST::zero());
      sStatus = (*solver_).solve(&*x);
      break;
    }
    case SOLVE_FOR_XDOT_CONST_X: {
      // if (getZeroInitialGuess()) Thyra::assign(xDot.ptr(), ST::zero());
      sStatus = (*solver_).solve(&*xDot);
      break;
    }
    default: {
      TEUCHOS_TEST_FOR_EXCEPT("Invalid EVALUATION_TYPE!");
    }
  }

  return sStatus;
}

template <class Scalar>
void StepperImplicit<Scalar>::evaluateImplicitODE(
    Teuchos::RCP<Thyra::VectorBase<Scalar> >& f,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
    const Teuchos::RCP<ImplicitODEParameters<Scalar> >& p)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> inArgs = wrapperModel_->createInArgs();
  inArgs.set_x(x);
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(xDot);
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(time);
  if (inArgs.supports(MEB::IN_ARG_step_size))
    inArgs.set_step_size(p->timeStepSize_);
  if (inArgs.supports(MEB::IN_ARG_alpha)) inArgs.set_alpha(Scalar(0.0));
  if (inArgs.supports(MEB::IN_ARG_beta)) inArgs.set_beta(Scalar(0.0));

  MEB::OutArgs<Scalar> outArgs = wrapperModel_->createOutArgs();
  outArgs.set_f(f);

  wrapperModel_->setForSolve(x, xDot, time, p);

  wrapperModel_->evalModel(inArgs, outArgs);
}

template <class Scalar>
void StepperImplicit<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << "--- StepperImplicit ---\n";
  out << "  wrapperModel_     = " << wrapperModel_ << std::endl;
  out << "  solver_           = " << solver_ << std::endl;
  out << "  initialGuess_     = " << initialGuess_ << std::endl;
  out << "  zeroInitialGuess_ = " << Teuchos::toString(zeroInitialGuess_)
      << std::endl;
}

template <class Scalar>
bool StepperImplicit<Scalar>::isValidSetup(Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if (wrapperModel_->getAppModel() == Teuchos::null) {
    isValidSetup = false;
    out << "The application ModelEvaluator is not set!\n";
  }

  if (wrapperModel_ == Teuchos::null) {
    isValidSetup = false;
    out << "The wrapper ModelEvaluator is not set!\n";
  }

  if (solver_ == Teuchos::null) {
    isValidSetup = false;
    out << "The solver is not set!\n";
  }

  if (solver_ != Teuchos::null) {
    if (solver_->getModel() == Teuchos::null) {
      isValidSetup = false;
      out << "The solver's ModelEvaluator is not set!\n";
    }
  }

  return isValidSetup;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperImplicit<Scalar>::getValidParameters() const
{
  return this->getValidParametersBasicImplicit();
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperImplicit<Scalar>::getValidParametersBasicImplicit() const
{
  auto pl = this->getValidParametersBasic();
  pl->template set<std::string>("Solver Name", this->getSolverName());
  pl->template set<bool>("Zero Initial Guess", this->getZeroInitialGuess());
  auto noxSolverPL = this->getSolver()->getParameterList();
  auto solverPL    = Teuchos::parameterList(this->getSolverName());
  solverPL->set("NOX", *noxSolverPL);
  pl->set(this->getSolverName(), *solverPL);

  return pl;
}

template <class Scalar>
void StepperImplicit<Scalar>::setStepperImplicitValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  if (pl != Teuchos::null) {
    // Can not validate because of optional Parameters, e.g., 'Solver Name'.
    // pl->validateParametersAndSetDefaults(*this->getValidParameters());
    this->setStepperValues(pl);
    this->setZeroInitialGuess(pl->get<bool>("Zero Initial Guess", false));
  }
  this->setStepperSolverValues(pl);
}

template <class Scalar>
void StepperImplicit<Scalar>::setStepperSolverValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  if (pl != Teuchos::null) {
    setDefaultSolver();
    std::string solverName = pl->get<std::string>("Solver Name");
    if (pl->isSublist(solverName)) {
      auto solverPL = Teuchos::parameterList();
      solverPL      = Teuchos::sublist(pl, solverName);
      this->setSolverName(solverName);
      Teuchos::RCP<Teuchos::ParameterList> noxPL =
          Teuchos::sublist(solverPL, "NOX", true);
      getSolver()->setParameterList(noxPL);
    }
  }
}

}  // namespace Tempus
#endif  // Tempus_StepperImplicit_impl_hpp
