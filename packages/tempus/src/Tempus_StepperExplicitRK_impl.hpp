// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRK_impl_hpp
#define Tempus_StepperExplicitRK_impl_hpp

#include "Tempus_RKButcherTableau.hpp"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Tempus {


template<class Scalar>
void StepperExplicitRK<Scalar>::setupDefault()
{
  this->setUseEmbedded(false);
  this->setStageNumber(-1);
  this->setAppAction(Teuchos::null);
}


template<class Scalar>
void StepperExplicitRK<Scalar>::setup(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool useEmbedded,
  const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
{
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setUseEmbedded(        useEmbedded);
  this->setStageNumber(-1);

  this->setAppAction(stepperRKAppAction);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}


template<class Scalar>
Scalar StepperExplicitRK<Scalar>::getInitTimeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& sh) const
{

   Scalar dt = Scalar(1.0e+99);
   if (!this->getUseEmbedded()) return dt;

   Teuchos::RCP<SolutionState<Scalar> > currentState=sh->getCurrentState();
   const int order = currentState->getOrder();
   const Scalar time = currentState->getTime();
   const Scalar errorRel = currentState->getTolRel();
   const Scalar errorAbs = currentState->getTolAbs();

   Teuchos::RCP<Thyra::VectorBase<Scalar> > stageX, scratchX;
   stageX = Thyra::createMember(this->appModel_->get_f_space());
   scratchX = Thyra::createMember(this->appModel_->get_f_space());
   Thyra::assign(stageX.ptr(), *(currentState->getX()));

   std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageXDot(2);
   for (int i=0; i<2; ++i) {
      stageXDot[i] = Thyra::createMember(this->appModel_->get_f_space());
      assign(stageXDot[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
   }

   // A: one function evaluation at F(t_0, X_0)
   typedef Thyra::ModelEvaluatorBase MEB;
   MEB::InArgs<Scalar> inArgs = this->appModel_->getNominalValues();
   MEB::OutArgs<Scalar> outArgs = this->appModel_->createOutArgs();
   inArgs.set_x(stageX);
   if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(time);
   if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(Teuchos::null);
   outArgs.set_f(stageXDot[0]); // K1
   this->appModel_->evalModel(inArgs, outArgs);

   auto err_func = [] (Teuchos::RCP<Thyra::VectorBase<Scalar> > U,
         const Scalar rtol, const Scalar atol,
         Teuchos::RCP<Thyra::VectorBase<Scalar> > absU)
   {
      // compute err = Norm_{WRMS} with w = Atol + Rtol * | U |
      Thyra::assign(absU.ptr(), *U);
      Thyra::abs(*U, absU.ptr()); // absU = | X0 |
      Thyra::Vt_S(absU.ptr(), rtol); // absU *= Rtol
      Thyra::Vp_S(absU.ptr(), atol); // absU += Atol
      Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *U, *absU, absU.ptr());
      Scalar err = Thyra::norm_inf(*absU);
      return err;
   };

   Scalar d0 = err_func(stageX, errorRel, errorAbs, scratchX);
   Scalar d1 = err_func(stageXDot[0], errorRel, errorAbs, scratchX);

   // b) first guess for the step size
   dt = Teuchos::as<Scalar>(0.01)*(d0/d1);

   // c) perform one explicit Euler step (X_1)
   Thyra::Vp_StV(stageX.ptr(), dt, *(stageXDot[0]));

   // compute F(t_0 + dt, X_1)
   inArgs.set_x(stageX);
   if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(time + dt);
   if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(Teuchos::null);
   outArgs.set_f(stageXDot[1]); // K2
   this->appModel_->evalModel(inArgs, outArgs);

   // d) compute estimate of the second derivative of the solution
   // d2 = || f(t_0 + dt, X_1) - f(t_0, X_0) || / dt
   Teuchos::RCP<Thyra::VectorBase<Scalar> > errX;
   errX = Thyra::createMember(this->appModel_->get_f_space());
   assign(errX.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
   Thyra::V_VmV(errX.ptr(), *(stageXDot[1]), *(stageXDot[0]));
   Scalar d2 = err_func(errX, errorRel, errorAbs, scratchX) / dt;

   // e) compute step size h_1 (from m = 0 order Taylor series)
   Scalar max_d1_d2 = std::max(d1, d2);
   Scalar h1 = std::pow((0.01/max_d1_d2),(1.0/(order+1)));

   // f) propose starting step size
   dt = std::min(100*dt, h1);
   return dt;
}


template<class Scalar>
void StepperExplicitRK<Scalar>::getValidParametersBasicERK(
  Teuchos::RCP<Teuchos::ParameterList> pl) const
{
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<bool>("Use Embedded", false,
    "'Whether to use Embedded Stepper (if available) or not\n"
    "  'true' - Stepper will compute embedded solution and is adaptive.\n"
    "  'false' - Stepper is not embedded(adaptive).\n");
  pl->set<std::string>("Description", this->getDescription());
}


template<class Scalar>
void StepperExplicitRK<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( this->tableau_ == Teuchos::null, std::logic_error,
    "Error - Need to set the tableau, before calling "
    "StepperExplicitRK::initialize()\n");

  TEUCHOS_TEST_FOR_EXCEPTION( this->appModel_==Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperExplicitRK::initialize()\n");

  // Initialize the stage vectors
  int numStages = this->tableau_->numStages();
  stageXDot_.resize(numStages);
  for (int i=0; i<numStages; ++i) {
    stageXDot_[i] = Thyra::createMember(this->appModel_->get_f_space());
    assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  }

  if ( this->tableau_->isEmbedded() and this->getUseEmbedded() ){
     this->ee_ = Thyra::createMember(this->appModel_->get_f_space());
     this->abs_u0 = Thyra::createMember(this->appModel_->get_f_space());
     this->abs_u = Thyra::createMember(this->appModel_->get_f_space());
     this->sc = Thyra::createMember(this->appModel_->get_f_space());
  }

  Stepper<Scalar>::initialize();
}


template<class Scalar>
void StepperExplicitRK<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  if (this->getUseFSAL())
    this->setStepperXDot(stageXDot_.back());
  else
    this->setStepperXDot(stageXDot_[0]);

  StepperExplicit<Scalar>::setInitialConditions(solutionHistory);

  auto xDot = solutionHistory->getCurrentState()->getXDot();
  if (xDot != Teuchos::null && this->getUseFSAL())
    Thyra::assign(this->getStepperXDot().ptr(), *(xDot));
}


template<class Scalar>
void StepperExplicitRK<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperExplicitRK::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperExplicitRK<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for ExplicitRK.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    const int numStages = this->tableau_->numStages();
    Teuchos::SerialDenseMatrix<int,Scalar> A = this->tableau_->A();
    Teuchos::SerialDenseVector<int,Scalar> b = this->tableau_->b();
    Teuchos::SerialDenseVector<int,Scalar> c = this->tableau_->c();

    Thyra::assign(workingState->getX().ptr(), *(currentState->getX()));

    RCP<StepperExplicitRK<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
      StepperRKAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    // Compute stage solutions
    for (int i=0; i < numStages; ++i) {
      this->setStageNumber(i);
      Thyra::assign(workingState->getX().ptr(), *(currentState->getX()));
      for (int j=0; j < i; ++j) {
        if (A(i,j) != Teuchos::ScalarTraits<Scalar>::zero()) {
          Thyra::Vp_StV(workingState->getX().ptr(), dt*A(i,j), *stageXDot_[j]);
        }
      }
      this->setStepperXDot(stageXDot_[i]);

      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::BEGIN_STAGE);
      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::BEFORE_SOLVE);
      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::AFTER_SOLVE);
      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);

      if ( i == 0 && this->getUseFSAL() &&
           workingState->getNConsecutiveFailures() == 0 ) {
        RCP<Thyra::VectorBase<Scalar> > tmp = stageXDot_[0];
        stageXDot_[0] = stageXDot_.back();
        stageXDot_.back() = tmp;
        this->setStepperXDot(stageXDot_[0]);

      } else {
        const Scalar ts = time + c(i)*dt;
        auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));

        // Evaluate xDot = f(x,t).
        this->evaluateExplicitODE(stageXDot_[i], workingState->getX(), ts, p);
      }

      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::END_STAGE);
    }
    // reset the stage number
    this->setStageNumber(-1);

    // Sum for solution: x_n = x_n-1 + Sum{ b(i) * dt*f(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int i=0; i < numStages; ++i) {
      if (b(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
        Thyra::Vp_StV((workingState->getX()).ptr(), dt*b(i), *(stageXDot_[i]));
      }
    }

    if (this->getUseFSAL()) {
      if (numStages == 1) {
        const Scalar ts = time + dt;
        auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));
        // Evaluate xDot = f(x,t).
        this->evaluateExplicitODE(stageXDot_[0], workingState->getX(), ts, p);
      }
      if (workingState->getXDot() != Teuchos::null)
        Thyra::assign((workingState->getXDot()).ptr(), *(stageXDot_.back()));
    }


    // At this point, the stepper has passed.
    // But when using adaptive time stepping, the embedded method
    // can change the step status
    workingState->setSolutionStatus(Status::PASSED);

    if (this->tableau_->isEmbedded() and this->getUseEmbedded()) {

      const Scalar tolRel = workingState->getTolRel();
      const Scalar tolAbs = workingState->getTolAbs();

      // just compute the error weight vector
      // (all that is needed is the error, and not the embedded solution)
      Teuchos::SerialDenseVector<int,Scalar> errWght = b ;
      errWght -= this->tableau_->bstar();

      //compute local truncation error estimate: | u^{n+1} - \hat{u}^{n+1} |
      // Sum for solution: ee_n = Sum{ (b(i) - bstar(i)) * dt*f(i) }
      assign(this->ee_.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
      for (int i=0; i < numStages; ++i) {
         if (errWght(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
            Thyra::Vp_StV(this->ee_.ptr(), dt*errWght(i), *(stageXDot_[i]));
         }
      }

      // compute: Atol + max(|u^n|, |u^{n+1}| ) * Rtol
      Thyra::abs( *(currentState->getX()), this->abs_u0.ptr());
      Thyra::abs( *(workingState->getX()), this->abs_u.ptr());
      Thyra::pair_wise_max_update(tolRel, *this->abs_u0, this->abs_u.ptr());
      Thyra::add_scalar(tolAbs, this->abs_u.ptr());

      //compute: || ee / sc ||
      assign(this->sc.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
      Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *this->ee_, *this->abs_u,this->sc.ptr());

      const auto space_dim = this->ee_->space()->dim();
      Scalar err = std::abs(Thyra::norm(*this->sc)) / space_dim ;
      workingState->setErrorRel(err);

      // test if step should be rejected
      if (std::isinf(err) || std::isnan(err) || err > Teuchos::as<Scalar>(1.0))
        workingState->setSolutionStatus(Status::FAILED);
    }

    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
      StepperRKAppAction<Scalar>::ACTION_LOCATION::END_STEP);
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
Teuchos::RCP<Tempus::StepperState<Scalar> > StepperExplicitRK<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperExplicitRK<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperExplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperExplicitRK ---\n";
  if (this->tableau_ != Teuchos::null) this->tableau_->describe(out, verbLevel);
  out << "  tableau_           = " << this->tableau_ << std::endl;
  out << "  stepperRKAppAction_= " << this->stepperRKAppAction_ << std::endl;
  out << "  stageXDot_.size()  = " << stageXDot_.size() << std::endl;
  const int numStages = stageXDot_.size();
  for (int i=0; i<numStages; ++i)
    out << "    stageXDot_["<<i<<"] = " << stageXDot_[i] << std::endl;
  out << "  useEmbedded_       = "
      << Teuchos::toString(this->useEmbedded_) << std::endl;
  out << "  ee_                = " << this->ee_ << std::endl;
  out << "  abs_u0             = " << this->abs_u0 << std::endl;
  out << "  abs_u              = " << this->abs_u << std::endl;
  out << "  sc                 = " << this->sc << std::endl;
  out << "-------------------------" << std::endl;
}


template<class Scalar>
bool StepperExplicitRK<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperExplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (this->tableau_ == Teuchos::null) {
    isValidSetup = false;
    out << "The tableau is not set!\n";
  }

  if (this->stepperRKAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The AppAction is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperExplicitRK<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  this->getValidParametersBasicERK(pl);
  return pl;
}


} // namespace Tempus
#endif // Tempus_StepperExplicitRK_impl_hpp
