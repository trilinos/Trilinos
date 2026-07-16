//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperExponential_impl_hpp
#define Tempus_StepperExponential_impl_hpp

#include "Tempus_StepperExponential.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"
#include "Tempus_PhiEvaluator.hpp"

namespace Tempus {

template <class Scalar>
StepperExponential<Scalar>::StepperExponential():
  appModel_(),
  phiEvaluator_(),
  temporalFiniteDifferenceEps_{1.0e-5},
  operatorLinearizationInterval_{-1},
  adaptPhiEvaluatorInterval_{-1}
{
  this->setStepperType("Exponential");
  this->setUseFSAL(false);
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);
}

template <class Scalar>
void StepperExponential<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  Stepper<Scalar>::describe(out, verbLevel);

  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  //*l_out << "\n--- " << this->description() << " ---" << std::endl;

  if ((Teuchos::as<int>(verbLevel) ==
       Teuchos::as<int>(Teuchos::VERB_DEFAULT)) ||
      (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW))) {
    *l_out << "--- StepperExponential ---\n";
    *l_out << "  addModel_     = " << appModel_ << std::endl;
    *l_out << "  Epsilon for RHS finite difference = " << temporalFiniteDifferenceEps_ << std::endl;
    *l_out << "  Operator Linearization Interval   = " << operatorLinearizationInterval_ << std::endl;
    *l_out << "  Adapt PhiEvaluator Interval       = " << adaptPhiEvaluatorInterval_ << std::endl;
  }

  if (phiEvaluator_ != Teuchos::null)
    phiEvaluator_->describe(out, verbLevel);
}

template<class Scalar>
bool StepperExponential<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (this->getPhiEvaluator() == Teuchos::null) {
    isValidSetup = false;
    out << "The PhiEvaluator is not set!\n";
  }

  if (this->getModel() == Teuchos::null) {
    isValidSetup = false;
    out << "The application ModelEvaluator is not set!\n";
  }

  return isValidSetup;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperExponential<Scalar>::getValidParameters() const
{
  return this->getValidParametersBasicExponential();
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperExponential<Scalar>::getValidParametersBasicExponential() const
{
  auto pl = this->getValidParametersBasic();

  // default values and docstrings for parameters are provided here:
  // TODO: document and potentially rename these options.
  pl->template set<double>("Epsilon for RHS finite difference", temporalFiniteDifferenceEps_);
  pl->template set<int>("Operator Linearization Interval", operatorLinearizationInterval_);
  pl->template set<int>("Adapt PhiEvaluator Interval", adaptPhiEvaluatorInterval_);

  // add the PhiEvaluator sublist
  auto phiPL = Teuchos::parameterList("PhiEvaluator");
  if (getPhiEvaluator() != Teuchos::null)
  {
    // get and copy the valid parameters and defaults of the current PhiEvaluator
    auto validPhiPL = this->getPhiEvaluator()->getValidParameters();
    phiPL = Teuchos::rcp(new Teuchos::ParameterList(*validPhiPL));
  }
  pl->set("PhiEvaluator", *phiPL);

  // TODO: should we remove this?
  // add some dummy variables for compatibility with StepperImplicit
  // this will validate and default initialize some parameters used by Implicit Steppers
  // to avoid throwing a validation error when those are provided
  pl->template set<std::string>("Solver Name", "Demo Solver");
  auto noxSolverPL = Tempus::defaultSolverParameters();
  auto solverPL = Teuchos::parameterList("Demo Solver");
  solverPL->set("NOX", *sublist(noxSolverPL, "NOX"));
  pl->set("Demo Solver", *solverPL);
  pl->template set<std::string>("Predictor Stepper Type", "None");
  pl->template set<bool>("Zero Initial Guess", 0);

  //pl->print(*this->getOStream());

  return pl;
}

template <class Scalar>
void StepperExponential<Scalar>::setStepperExponentialValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  Teuchos::RCP<Teuchos::ParameterList> phiPL = Teuchos::null;

  if (pl != Teuchos::null) {
    if (pl->isSublist("PhiEvaluator")){
      phiPL = sublist(pl, "PhiEvaluator");

      // we always construct a PhiEvaluator at initialization, but the selected
      // PhiEvaluator depends on the pList, and we need to reinitialize it here
      auto phif = Teuchos::rcp(new PhiEvaluatorFactory<Scalar>());
      auto phiEvaluator = phif->createPhiEvaluator(phiPL);
      this->setPhiEvaluator(phiEvaluator);
    }

    // validate that the parameters from the user are valid, and set defaults for missing values
    pl->validateParametersAndSetDefaults(*this->getValidParameters());

    // set the validated values, or their defaults
    this->setStepperValues(pl);
    temporalFiniteDifferenceEps_ = pl->get<double>("Epsilon for RHS finite difference");
    operatorLinearizationInterval_ = pl->get<int>("Operator Linearization Interval");
    adaptPhiEvaluatorInterval_ = pl->get<int>("Adapt PhiEvaluator Interval");
  }

}

template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperExponential<Scalar>::getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperStateExponential<Scalar> > stepperState =
    rcp(new StepperStateExponential<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void StepperExponential<Scalar>::setDefaultPhiEvaluator()
{
  // get a default PhiEvaluator
  auto phif = Teuchos::rcp(new PhiEvaluatorFactory<Scalar>());
  auto phiEvaluator = phif->createPhiEvaluator();

  this->setPhiEvaluator(phiEvaluator);
}

template<class Scalar>
void StepperExponential<Scalar>::setPhiEvaluator(
  const Teuchos::RCP<Tempus::PhiEvaluator<Scalar> >& phiEvaluator)
{
  phiEvaluator_ = phiEvaluator;
  if (getModel() != Teuchos::null)
  {
    phiEvaluator_->setModel(getModel());
    phiEvaluator_->initialize();
  }
  this->isInitialized_ = false;
}

template<class Scalar>
Teuchos::RCP<Tempus::PhiEvaluator<Scalar>> StepperExponential<Scalar>::getPhiEvaluator() const
{
  return phiEvaluator_;
}

template <class Scalar>
void StepperExponential<Scalar>::setModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  validImplicitODE_DAE(appModel);
  appModel_ = appModel;

  if (phiEvaluator_ != Teuchos::null)
  {
    phiEvaluator_->setModel(appModel);
    phiEvaluator_->initialize();
  }

  this->isInitialized_ = false;
}

template<class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> StepperExponential<Scalar>::getModel() const
{
  return appModel_;
}

template<class Scalar>
void StepperExponential<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using magScalar = typename ST::magnitudeType;
  typedef Teuchos::ScalarTraits<magScalar> STM;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(
      numStates < 1, std::logic_error,
      "Error - setInitialConditions() needs at least one SolutionState\n"
      "        to set the initial condition.  Number of States = "
          << numStates);
  TEUCHOS_TEST_FOR_EXCEPTION(
      (this->getOrderODE() == SECOND_ORDER_ODE), std::logic_error,
      "Error - StepperEPI does not support SECOND_ORDER_ODE.\n");

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();

  RCP<Thyra::VectorBase<Scalar> > x = initialState->getX();
  RCP<Thyra::VectorBase<Scalar> > xDot = initialState->getXDot();

  auto inArgs = this->getModel()->getNominalValues();
  // Use x from inArgs as ICs if null.
  if (x == Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        (x == Teuchos::null) && (inArgs.get_x() == Teuchos::null),
        std::logic_error,
        "Error - setInitialConditions needs the ICs from the "
        "SolutionHistory\n"
        "        or getNominalValues()!\n");

    x = inArgs.get_x()->clone_v();
    initialState->setX(x);
  }

  // Reset the lastLinearizationPoint flag to ensure that a fresh Jacobian will be computed in the first step
  RCP<const StepperStateExponential<Scalar>> ss_exp =
    Teuchos::rcp_dynamic_cast<const StepperStateExponential<Scalar>>(initialState.getConst()->getStepperState());
  RCP<StepperState<Scalar>> ss_nc;
  if (ss_exp == Teuchos::null) {
    ss_nc = this->getDefaultStepperState();
  }
  else {
    // ensure that the stepper state is not constant, by resetting to a nonconst copy
    ss_nc = ss_exp->clone();
  }
  initialState->setStepperState(ss_nc);
  RCP<StepperStateExponential<Scalar>> ss_exp_nc =
    Teuchos::rcp_dynamic_cast<StepperStateExponential<Scalar>>(initialState->getStepperState());
  ss_exp_nc->lastLinearizationPoint_ = -1;

  // check and remember if the solutionHistory stores xdot
  bool xDotHistoryStored = !(xDot == Teuchos::null);

  // Check if we need to create local Stepper storage for xDot
  if (!xDotHistoryStored) {
    xDot = Thyra::createMember(x->space());
    Thyra::assign(xDot.ptr(), ST::zero());
    this->setStepperXDot(xDot);
  }
  else
    this->setStepperXDot(xDot);

  magScalar reldiff = STM::zero();
  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    // Mark xDot as not synced
    initialState->setIsSynced(false);
  }
  else if (icConsistency == "Zero") {
    // assign zero to xDot, even if something was provided in currentState
    Thyra::assign(xDot.ptr(), ST::zero());

    // Mark xDot as not synced
    initialState->setIsSynced(false);
  }
  else if (icConsistency == "App") {
    auto x_dot = inArgs.get_x_dot();
    TEUCHOS_TEST_FOR_EXCEPTION(
        x_dot == Teuchos::null, std::logic_error,
        "Error - setInitialConditions() requested 'App' for IC consistency,\n"
        "        but 'App' returned a null pointer for xDot!\n");
    Thyra::assign(xDot.ptr(), *x_dot);

    // Mark xDot as not synced
    initialState->setIsSynced(false);
  }
  else if (icConsistency == "Consistent") {
    TEUCHOS_TEST_FOR_EXCEPTION(
        !xDotHistoryStored, std::logic_error,
        "Error - icConsistency = 'Consistent' requested, but xDot not provided in history.\n");
    // TODO: since the user requested a "Consistent" initial condition,
    // should we make sure that xDot will be saved to the history, even if Teuchos::zero initially?
    // initialState->setXDot(xDot);
    // xDotHistoryStored = true;

    // xDot will be set to an all zero vector, to signal to the ModelEvaluator that we desire the implicit mode
    // it must remain zero until the last call to ModelEvaluator, until the end of this method
    Thyra::assign(xDot.ptr(), ST::zero());

    const Scalar t0 = initialState->getTime();
    const Scalar dt = initialState->getTimeStep();

    auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));
    // compute the right hand side for x
    // and potentially set the correct Dirichlet BC to x
    RCP<Thyra::VectorBase<Scalar>> Mf = Thyra::createMember(x->space());
    this->evaluateExponentialODE(Mf, x, xDot, t0, p);

    // setup system mass matrix
    phiEvaluator_->setLinearizationPoint(inArgs, PhiInitialization::ONLY_MASS);

    // overwrite xDot with f = (-M) \ Mf
    Thyra::scale(Scalar(-1.0), Mf.ptr());
    phiEvaluator_->solveMass(xDot.ptr(), Mf);
    reldiff = STM::zero();  // TODO: make solveMass return the tolerance.

    // At this point, x and xDot are sync'ed or consistent
    // at the same time level for the initialState.
    initialState->setIsSynced(true);

    // TODO: we do not track potential failures in mass solve
    // TEUCHOS_TEST_FOR_EXCEPTION(
    //    sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED,
    //    std::logic_error,
    //    "Error - Solver failed while determining the initial conditions.\n"
    //    "        Solver status is "
    //        << Thyra::toString(sStatus.solveStatus) << ".\n");
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error - setInitialConditions() invalid IC consistency, "
            << icConsistency << ".\n");
  }

  // Test for consistency.
  if (this->getICConsistencyCheck()) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        !xDotHistoryStored, std::logic_error,
        "Error - ICConsistenceCheck requested, but xDot not provided in history.\n");

    RCP<Thyra::VectorBase<Scalar>> xDot_temp = Thyra::createMember(x->space());
    Thyra::assign(xDot_temp.ptr(), ST::zero());

    const Scalar t0 = initialState->getTime();
    const Scalar dt = initialState->getTimeStep();

    auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));
    // compute the right hand side for x
    // and potentially set the correct Dirichlet BC to x
    RCP<Thyra::VectorBase<Scalar>> Mf = Thyra::createMember(x->space());
    this->evaluateExponentialODE(Mf, x, xDot_temp, t0, p);

    if (icConsistency != "Consistent")
      // setup system mass matrix
      phiEvaluator_->setLinearizationPoint(inArgs, PhiInitialization::ONLY_MASS);

    // overwrite xDot_temp with Mf - (-M) * xDot
    phiEvaluator_->applyMass(xDot_temp.ptr(), xDot);
    Thyra::Vp_V(xDot_temp.ptr(), *Mf);

    magScalar normX = Thyra::norm(*x);
    if (normX == STM::zero())
      reldiff = Thyra::norm(*xDot_temp);
    else
      reldiff = Thyra::norm(*xDot_temp) / normX;

    magScalar eps = magScalar(100.0) * STM::eps();
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
           << "  ||f(x,xDot,t)|| = " << Thyra::norm(*xDot_temp) << std::endl
           << "  ||x||           = " << Thyra::norm(*x) << std::endl
           << "---------------------------------------------------\n"
           << std::endl;
    }
  }
}


template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
StepperExponential<Scalar>::createInArgsExponentialODE(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
    const Teuchos::RCP<ExponentialODEParameters<Scalar> >& p)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel = this->getModel();
  MEB::InArgs<Scalar> inArgs = appModel->createInArgs();
  inArgs.set_x(x);
  TEUCHOS_ASSERT(inArgs.supports(MEB::IN_ARG_x_dot))
  inArgs.set_x_dot(xDot);

  if (inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(time);
  if (inArgs.supports(MEB::IN_ARG_step_size))
    inArgs.set_step_size(p->timeStepSize_);
  if (inArgs.supports(MEB::IN_ARG_alpha))
    inArgs.set_alpha(Scalar(0.0));
  if (inArgs.supports(MEB::IN_ARG_beta))
    inArgs.set_beta(Scalar(1.0));
  if (inArgs.supports(MEB::IN_ARG_stage_number))
    inArgs.set_stage_number(p->stageNumber_);

  return inArgs;
}


template <class Scalar>
void StepperExponential<Scalar>::evaluateExponentialODE(
    Teuchos::RCP<Thyra::VectorBase<Scalar> >& f,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
    const Teuchos::RCP<ExponentialODEParameters<Scalar> >& p)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel = this->getModel();

  // TODO: we could also rely on the user (takeStep) to set this, but double setting may be more robust
  Thyra::assign(xDot.ptr(), ST::zero());

  MEB::InArgs<Scalar> inArgs = createInArgsExponentialODE(x, xDot, time, p);
  MEB::OutArgs<Scalar> outArgs = appModel->createOutArgs();
  outArgs.set_f(f);
  appModel->evalModel(inArgs, outArgs);
}


template <class Scalar>
void StepperExponential<Scalar>::computeTemporalFD(
    Teuchos::RCP<Thyra::VectorBase<Scalar>>& dt_Mf_deriv,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf
)
{
  // evaluate Mf at t + dt * eps
  auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));
  Teuchos::RCP<Thyra::VectorBase<Scalar>> x_BC = x->clone_v();
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xDot = this->getStepperXDot();
  this->evaluateExponentialODE(
      dt_Mf_deriv, x_BC, xDot, t0 + dt * this->temporalFiniteDifferenceEps_, p);
  // this sets the Dirichlet BC in x_BC to the value at t0 + dt * eps
  // thus we need to use a temporary x_BC instead of x here

  // compute dt times the temporal finite difference of Mf:
  // dt_Mf_deriv = (Mf(t + eps*dt) - Mf(t)) / eps
  Scalar one_over_eps = Scalar(1. / this->temporalFiniteDifferenceEps_);
  Thyra::linear_combination<Scalar>(Teuchos::tuple(-one_over_eps),
                                    Teuchos::tuple(Mf.getConst().ptr()),
                                    one_over_eps, dt_Mf_deriv.ptr());

  // also compute the finite difference (x_BC - x) / eps, reusing storage
  Teuchos::RCP<Thyra::VectorBase<Scalar>> dt_x_BC = x_BC;
  Thyra::linear_combination<Scalar>(Teuchos::tuple(-one_over_eps),
                                    Teuchos::tuple(x.getConst().ptr()),
                                    one_over_eps, dt_x_BC.ptr());

  // add dt_x_BC to dt_Mf_deriv:
  // this simple logic should take care of _either_ inhomogeneous BC or RHS.
  //
  // TODO: In case both are inhomogeneous, dt_Mf_deriv contains two contributions:
  // dt_Mf_deriv = (Mf(t + eps*dt, x + eps*dt_x_BC) - Mf(t, x)) / eps
  //             = (Mf(t + eps*dt, x) - Mf(t, x)) / eps
  //             + (Mf(t + eps*dt, x + eps*dt_x_BC) - Mf(t + eps*dt, x)) / eps
  // The second term in the sum is approximately J(t + eps*dt, x) * dt_x_BC,
  // which we may have to correct for
  Thyra::Vp_V(dt_Mf_deriv.ptr(), *dt_x_BC);
}


template <class Scalar>
void StepperExponential<Scalar>::computeRemf(
    Teuchos::RCP<Thyra::VectorBase<Scalar>>& remf,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& xr,
    const Scalar tr,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x0,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& dt_Mf_deriv,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mfr
)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  // xd = xr - x0
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xd = Thyra::createMember(x0->space());
  Thyra::V_VpStV(xd.ptr(), *xr, Scalar(-1.0), *x0);

  // remf = J_xd = -M*J * (xr - x0)
  Thyra::assign(remf.ptr(), ST::zero());
  phiEvaluator_->applyJacobian(remf.ptr(), xd);
  // TODO: could save one temp vector by extending the interface of apply Jacobian

  // Evaluate the rhs -M*F at (xr, tr)
  //   Mfr = -M*F(xr, tr) = F_impl(xDot = 0, xr, tr)
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mfr_const;
  if (Mfr != Teuchos::null) {
    // use the provided Mfr vector
    Mfr_const = Mfr;
  }
  else {
    auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));
    Teuchos::RCP<Thyra::VectorBase<Scalar>> Mfr = Thyra::createMember(x0->space());
    Teuchos::RCP<Thyra::VectorBase<Scalar>> xr_temp = xd; // reuse xd memory
    Thyra::copy(*xr, xr_temp.ptr());
    // this will reset xDot to zero and re-set the BC in xr (that is why we need a non_const temporary)
    // if we save Mfr in the previous iteration correctly, we can reuse it.
    Teuchos::RCP<Thyra::VectorBase<Scalar>> xDot = this->getStepperXDot();
    this->evaluateExponentialODE(Mfr, xr_temp, xDot, tr, p);
    Mfr_const = Mfr;
  }

  // update remf = (Mf - Mfr) + J_xd
  // Mf is rhs at the current time: -M*F(x0, t0)
  Thyra::linear_combination<Scalar>(Teuchos::tuple(Scalar(1.0), Scalar(-1.0)),
                                    Teuchos::tuple(Mf.getConst().ptr(), Mfr.getConst().ptr()),
                                    Scalar(1.0), remf.ptr());

  // add time derivative remainder term; only nonzero in nonautonomous case
  // add  -(M * F') * (tr - t0) = -(dt * M * F') * ((tr - t0) / dt)
  if (dt_Mf_deriv != Teuchos::null) {
    Thyra::Vp_StV(remf.ptr(), Scalar((tr - t0) / dt), *dt_Mf_deriv);
  }
}


template <class Scalar>
bool StepperExponential<Scalar>::needsOperatorLinearization(
    const Teuchos::RCP<const Tempus::SolutionState<Scalar>>& currentState,
    const Teuchos::RCP<Tempus::SolutionState<Scalar>>& workingState
  )
{
  // get a const StepperStateExponential, if possible (nonconst version triggers assertion if not present)
  auto ss_exp =
    Teuchos::rcp_dynamic_cast<const StepperStateExponential<Scalar>>(workingState.getConst()->getStepperState());

  bool needsOpLin =
    !this->getPhiEvaluator()->checkLinearizationPoint(PhiInitialization::JACOBIAN_AND_MASS)
    || (ss_exp == Teuchos::null)
    || ss_exp->lastLinearizationPoint_ < 0
    || (currentState->getIndex() - ss_exp->lastLinearizationPoint_ >= getOperatorLinearizationInterval());

  //auto out = this->getOStream();
  //workingState->getStepperState()->describe(*out, Teuchos::VERB_EXTREME);
  //*out << "index: " << currentState->getIndex() << std::endl;

  if (needsOpLin)
  {
    Teuchos::RCP<StepperState<Scalar>> ss_nc;
    if (ss_exp == Teuchos::null) {
      ss_nc = this->getDefaultStepperState();
    }
    else {
      // ensure that the stepper state is not constant, by resetting to a nonconst copy
      ss_nc = ss_exp->clone();
    }
    workingState->setStepperState(ss_nc);
    auto ss_exp_nc =
      Teuchos::rcp_dynamic_cast<StepperStateExponential<Scalar>>(workingState->getStepperState());
    ss_exp_nc->lastLinearizationPoint_ = currentState->getIndex();
  }
  //workingState->getStepperState()->describe(*out, Teuchos::VERB_EXTREME);

  return needsOpLin;
}


template <class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar>> StepperStateExponential<Scalar>::clone() const
{
  Teuchos::RCP<StepperStateExponential<Scalar>> ss_out =
    Teuchos::rcp(new StepperStateExponential<Scalar>(this->stepperName_));
  ss_out->lastLinearizationPoint_ = this->lastLinearizationPoint_;
  return ss_out;
}

template <class Scalar>
void StepperStateExponential<Scalar>::copy(const Teuchos::RCP<const StepperState<Scalar> >& ss)
{
  this->stepperName_ = ss->stepperName_;
  auto ss_exp = Teuchos::rcp_dynamic_cast<const StepperStateExponential<Scalar>>(ss);
  if (ss_exp != Teuchos::null)
  {
    this->lastLinearizationPoint_ = ss_exp->lastLinearizationPoint_;
  }
}

template <class Scalar>
void StepperStateExponential<Scalar>::describe(
    Teuchos::FancyOStream& out,
    const Teuchos::EVerbosityLevel verbLevel) const
{
  Tempus::StepperState<Scalar>::describe(out, verbLevel);

  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << "last linearization point: " << lastLinearizationPoint_ << std::endl;
}

}  // namespace Tempus
#endif  // Tempus_StepperExponential_impl_hpp
