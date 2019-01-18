// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkExplicitAForm_impl_hpp
#define Tempus_StepperNewmarkExplicitAForm_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"

//#define DEBUG_OUTPUT

namespace Tempus {


template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::
predictVelocity(Thyra::VectorBase<Scalar>& vPred,
                const Thyra::VectorBase<Scalar>& v,
                const Thyra::VectorBase<Scalar>& a,
                const Scalar dt) const
{
  //vPred = v + dt*(1.0-gamma_)*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0, v, dt*(1.0-gamma_), a);
}

template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::
predictDisplacement(Thyra::VectorBase<Scalar>& dPred,
                    const Thyra::VectorBase<Scalar>& d,
                    const Thyra::VectorBase<Scalar>& v,
                    const Thyra::VectorBase<Scalar>& a,
                    const Scalar dt) const
{
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > tmp =
    Thyra::createMember<Scalar>(dPred.space());
  //dPred = dt*v + dt*dt/2.0*a
  Scalar aConst = dt*dt/2.0;
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), dt, v, aConst, a);
  //dPred += d;
  Thyra::Vp_V(Teuchos::ptrFromRef(dPred), d, 1.0);
}

template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::
correctVelocity(Thyra::VectorBase<Scalar>& v,
                const Thyra::VectorBase<Scalar>& vPred,
                const Thyra::VectorBase<Scalar>& a,
                const Scalar dt) const
{
  //v = vPred + dt*gamma_*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(v), 1.0, vPred, dt*gamma_, a);
}


template<class Scalar>
StepperNewmarkExplicitAForm<Scalar>::StepperNewmarkExplicitAForm(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList) :
  out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  // Set all the input parameters and call initialize
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}

template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->validExplicitODE(appModel);
  appModel_ = appModel;

  inArgs_  = appModel_->createInArgs();
  outArgs_ = appModel_->createOutArgs();
  inArgs_  = appModel_->getNominalValues();
}

template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}

template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::setSolver(std::string solverName)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperNewmarkExplicitAForm::setSolver()");
  *out << "Warning -- No solver to set for StepperNewmarkExplicitAForm "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperNewmarkExplicitAForm::setSolver()");
  *out << "Warning -- No solver to set for StepperNewmarkExplicitAForm "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperNewmarkExplicitAForm::setSolver()");
  *out << "Warning -- No solver to set for StepperNewmarkExplicitAForm "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->appModel_ == Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperNewmarkExplicitAForm::initialize()\n");

  this->setParameterList(this->stepperPL_);
}

template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperNewmarkExplicitAForm::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperNewmarkExplicitAForm<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for NewmarkExplicitAForm.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();

    //Get values of d, v and a from previous step
    RCP<const Thyra::VectorBase<Scalar> > d_old = currentState->getX();
    RCP<const Thyra::VectorBase<Scalar> > v_old = currentState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > a_old = currentState->getXDotDot();

    //Get dt and time
    const Scalar dt = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    typedef Thyra::ModelEvaluatorBase MEB;

#ifdef DEBUG_OUTPUT
    *out_ << "IKT d_old = " << Thyra::max(*d_old) << "\n";
    *out_ << "IKT v_old = " << Thyra::max(*v_old) << "\n";
    *out_ << "IKT a_old prescribed = " << Thyra::max(*a_old) << "\n";
#endif

    //Compute initial acceleration, a_old, using initial displacement (d_old) and initial
    //velocity (v_old) if in 1st time step
    //allocate a_init
    RCP<Thyra::VectorBase<Scalar> > a_init = Thyra::createMember(d_old->space());
    Thyra::put_scalar(0.0, a_init.ptr());
    if (time == solutionHistory->minTime()) {
      //Set x and x_dot in inArgs_ to be initial d and v, respectively
      inArgs_.set_x(d_old);
      inArgs_.set_x_dot(v_old);
      if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time);
      // For model evaluators whose state function f(x, x_dot, x_dotdot, t) describes
      // an implicit ODE, and which accept an optional x_dotdot input argument,
      // make sure the latter is set to null in order to request the evaluation
      // of a state function corresponding to the explicit ODE formulation
      // x_dotdot = f(x, x_dot, t)
      if (inArgs_.supports(MEB::IN_ARG_x_dot_dot)) inArgs_.set_x_dot_dot(Teuchos::null);
      outArgs_.set_f(a_init);
      appModel_->evalModel(inArgs_,outArgs_);
      Thyra::copy(*a_init, a_old.ptr());
#ifdef DEBUG_OUTPUT
      *out_ << "IKT a_init computed = " << Thyra::max(*a_old) << "\n";
#endif
    }


    //New d, v and a to be computed here
    RCP<Thyra::VectorBase<Scalar> > d_new = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > v_new = workingState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > a_new = workingState->getXDotDot();

    //allocate d and v predictors
    RCP<Thyra::VectorBase<Scalar> > d_pred =Thyra::createMember(d_old->space());
    RCP<Thyra::VectorBase<Scalar> > v_pred =Thyra::createMember(v_old->space());

    //compute displacement and velocity predictors
    predictDisplacement(*d_pred, *d_old, *v_old, *a_old, dt);
    predictVelocity(*v_pred, *v_old, *a_old, dt);

#ifdef DEBUG_OUTPUT
    *out_ << "IKT d_pred = " << Thyra::max(*d_pred) << "\n";
    *out_ << "IKT v_pred = " << Thyra::max(*v_pred) << "\n";
#endif

    //Set x and x_dot in inArgs_ to be d and v predictors, respectively
    inArgs_.set_x(d_pred);
    inArgs_.set_x_dot(v_pred);
    if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(currentState->getTime());

    // For model evaluators whose state function f(x, x_dot, x_dotdot, t) describes
    // an implicit ODE, and which accept an optional x_dotdot input argument,
    // make sure the latter is set to null in order to request the evaluation
    // of a state function corresponding to the explicit ODE formulation
    // x_dotdot = f(x, x_dot, t)
    if (inArgs_.supports(MEB::IN_ARG_x_dot_dot)) inArgs_.set_x_dot_dot(Teuchos::null);
    outArgs_.set_f(a_old);

    appModel_->evalModel(inArgs_,outArgs_);

    Thyra::copy(*(outArgs_.get_f()), a_new.ptr());
#ifdef DEBUG_OUTPUT
    *out_ << "IKT a_new = " << Thyra::max(*(workingState()->getXDotDot())) << "\n";
#endif

    //Set x in workingState to displacement predictor
    Thyra::copy(*d_pred, d_new.ptr());

    //set xdot in workingState to velocity corrector
    correctVelocity(*v_new, *v_pred, *a_new, dt);
#ifdef DEBUG_OUTPUT
    *out_ << "IKT d_new = " << Thyra::max(*(workingState()->getX())) << "\n";
    *out_ << "IKT v_new = " << Thyra::max(*(workingState()->getXDot())) << "\n";
#endif

    workingState->setSolutionStatus(Status::PASSED);
    workingState->setOrder(this->getOrder());
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
Teuchos::RCP<Tempus::StepperState<Scalar> > StepperNewmarkExplicitAForm<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperNewmarkExplicitAForm<Scalar>::description() const
{
  std::string name = "Newmark Explicit a-Form";
  return(name);
}


template<class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "appModel_ = " << appModel_->description() << std::endl;
}


template <class Scalar>
void StepperNewmarkExplicitAForm<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (stepperPL_ == Teuchos::null) stepperPL_ = this->getDefaultParameters();
  } else {
    stepperPL_ = pList;
  }
  stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string stepperType = stepperPL_->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Newmark Explicit a-Form",
    std::logic_error,
       "Error - Stepper Type is not 'Newmark Explicit a-Form'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");
  gamma_ = 0.5; //default value
  if (stepperPL_->isSublist("Newmark Explicit Parameters")) {
    Teuchos::ParameterList &newmarkPL =
      stepperPL_->sublist("Newmark Explicit Parameters", true);
    gamma_ = newmarkPL.get("Gamma", 0.5);
    *out_ << "\nSetting Gamma = " << gamma_ << " from input file.\n";
  }
  TEUCHOS_TEST_FOR_EXCEPTION( (gamma_ > 1.0) || (gamma_ < 0.0),
      std::logic_error,
      "\nError in 'Newmark Explicit a-Form' stepper: invalid value of Gamma = " <<gamma_ << ".  Please select Gamma >= 0 and <= 1. \n");

}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperNewmarkExplicitAForm<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set("Stepper Type", "Newmark Explicit a-Form",
          "'Stepper Type' must be 'Newmark Explicit a-Form'.");
  pl->sublist("Newmark Explicit Parameters", false, "");
  pl->sublist("Newmark Explicit Parameters", false, "").set("Gamma",
               0.5, "Newmark Explicit parameter");

  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperNewmarkExplicitAForm<Scalar>::getDefaultParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  *pl = *(this->getValidParameters());
  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperNewmarkExplicitAForm<Scalar>::getNonconstParameterList()
{
  return(stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperNewmarkExplicitAForm<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperNewmarkExplicitAForm_impl_hpp
