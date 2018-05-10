// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkImplicitAForm_impl_hpp
#define Tempus_StepperNewmarkImplicitAForm_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "NOX_Thyra.H"

//#define VERBOSE_DEBUG_OUTPUT
//#define DEBUG_OUTPUT

namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;

template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::
predictVelocity(Thyra::VectorBase<Scalar>& vPred,
                const Thyra::VectorBase<Scalar>& v,
                const Thyra::VectorBase<Scalar>& a,
                const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //vPred = v + dt*(1.0-gamma_)*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0, v, dt*(1.0-gamma_), a);
}

template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::
predictDisplacement(Thyra::VectorBase<Scalar>& dPred,
                    const Thyra::VectorBase<Scalar>& d,
                    const Thyra::VectorBase<Scalar>& v,
                    const Thyra::VectorBase<Scalar>& a,
                    const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > tmp =
    Thyra::createMember<Scalar>(dPred.space());
  //dPred = dt*v + dt*dt/2.0*(1.0-2.0*beta_)*a
  Scalar aConst = dt*dt/2.0*(1.0-2.0*beta_);
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), dt, v, aConst, a);
  //dPred += d;
  Thyra::Vp_V(Teuchos::ptrFromRef(dPred), d, 1.0);
}

template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::
correctVelocity(Thyra::VectorBase<Scalar>& v,
                const Thyra::VectorBase<Scalar>& vPred,
                const Thyra::VectorBase<Scalar>& a,
                const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //v = vPred + dt*gamma_*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(v), 1.0, vPred, dt*gamma_, a);
}

template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::
correctDisplacement(Thyra::VectorBase<Scalar>& d,
                    const Thyra::VectorBase<Scalar>& dPred,
                    const Thyra::VectorBase<Scalar>& a,
                    const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //d = dPred + beta_*dt*dt*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(d), 1.0, dPred, beta_*dt*dt, a);
}


template<class Scalar>
StepperNewmarkImplicitAForm<Scalar>::StepperNewmarkImplicitAForm(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList) :
  out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Set all the input parameters and call initialize
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->validSecondOrderODE_DAE(appModel);
  Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
    Teuchos::rcp(new WrapperModelEvaluatorSecondOrder<Scalar>(appModel,
                                              "Newmark Implicit a-Form"));
  this->wrapperModel_ = wrapperModel;
  inArgs_  = this->wrapperModel_->getNominalValues();
  outArgs_ = this->wrapperModel_->createOutArgs();
}


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( this->wrapperModel_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperNewmarkImplicitAForm::initialize()\n");

#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->setSolver();
}


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperNewmarkImplicitAForm::takeStep()");
  {
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
      Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar> >(
        this->wrapperModel_);

    //Get values of d, v and a from previous step
    RCP<const Thyra::VectorBase<Scalar> > d_old = currentState->getX();
    RCP<const Thyra::VectorBase<Scalar> > v_old = currentState->getXDot();
    //RCP<const Thyra::VectorBase<Scalar> > a_old = currentState->getXDotDot();
    RCP<Thyra::VectorBase<Scalar> > a_old = currentState->getXDotDot();

#ifdef DEBUG_OUTPUT
    //IKT, 3/21/17, debug output: pring d_old, v_old to check for
    // correctness.
    *out_ << "IKT d_old = " << Thyra::max(*d_old) << "\n";
    *out_ << "IKT v_old = " << Thyra::max(*v_old) << "\n";
#endif

    //Get new values of d, v and a from current workingState
    //(to be updated here)
    RCP<Thyra::VectorBase<Scalar> > d_new = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > v_new = workingState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > a_new = workingState->getXDotDot();

    //Get time and dt
    const Scalar time = currentState->getTime();
    const Scalar dt = workingState->getTimeStep();
    //Update time
    Scalar t = time+dt;

    //Compute initial acceleration, a_old, using initial displacement (d_old) and initial
    //velocity (v_old) if in 1st time step
    if (time == solutionHistory->minTime()) {
      RCP<Thyra::VectorBase<Scalar> > d_init = Thyra::createMember(d_old->space());
      RCP<Thyra::VectorBase<Scalar> > v_init = Thyra::createMember(v_old->space());
      RCP<Thyra::VectorBase<Scalar> > a_init = Thyra::createMember(a_old->space());
      Thyra::copy(*d_old, d_init.ptr());
      Thyra::copy(*v_old, v_init.ptr());
      Thyra::put_scalar(0.0, a_init.ptr());
      wrapperModel->initializeNewmark(a_init,v_init,d_init,0.0,time,beta_,gamma_);
      const Thyra::SolveStatus<Scalar> sStatus =
        this->solveImplicitODE(a_init);

      if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED )
        workingState->getStepperState()->stepperStatus_ = Status::PASSED;
      else
        workingState->getStepperState()->stepperStatus_ = Status::FAILED;
      Thyra::copy(*a_init, a_old.ptr());
    }
#ifdef DEBUG_OUTPUT
    //IKT, 3/30/17, debug output: pring a_old to check for correctness.
    *out_ << "IKT a_old = " << Thyra::max(*a_old) << "\n";
#endif


    //allocate d and v predictors
    RCP<Thyra::VectorBase<Scalar> > d_pred =Thyra::createMember(d_old->space());
    RCP<Thyra::VectorBase<Scalar> > v_pred =Thyra::createMember(v_old->space());

    //compute displacement and velocity predictors
    predictDisplacement(*d_pred, *d_old, *v_old, *a_old, dt);
    predictVelocity(*v_pred, *v_old, *a_old, dt);

    //inject d_pred, v_pred, a and other relevant data into wrapperModel
    wrapperModel->initializeNewmark(a_old,v_pred,d_pred,dt,t,beta_,gamma_);

    const Thyra::SolveStatus<Scalar> sStatus = this->solveImplicitODE(a_old);

    if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED )
      workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    else
      workingState->getStepperState()->stepperStatus_ = Status::FAILED;

    Thyra::copy(*a_old, a_new.ptr());
    correctVelocity(*v_new, *v_pred, *a_new, dt);
    correctDisplacement(*d_new, *d_pred, *a_new, dt);

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
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperNewmarkImplicitAForm<Scalar>::
getDefaultStepperState()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperNewmarkImplicitAForm<Scalar>::description() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  std::string name = "Newmark Implicit a-Form";
  return(name);
}


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  out << description() << "::describe:" << std::endl
      << "wrapperModel = " << this->wrapperModel_->description() << std::endl;
}


template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (this->stepperPL_ == Teuchos::null) this->stepperPL_ = this->getDefaultParameters();
  } else {
    this->stepperPL_ = pList;
  }
  // Can not validate because of optional Parameters.
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
  //Get beta and gamma from parameter list
  //IKT, FIXME: does parameter list get validated somewhere?  validateParameters above is commented out...

  Teuchos::RCP<Teuchos::ParameterList> stepperPL = this->stepperPL_;
  std::string stepperType = stepperPL->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Newmark Implicit a-Form",
    std::logic_error,
       "Error - Stepper Type is not 'Newmark Implicit a-Form'!\n"
       << "  Stepper Type = "<< stepperPL->get<std::string>("Stepper Type")
       << "\n");
  beta_ = 0.25; //default value
  gamma_ = 0.5; //default value
    Teuchos::VerboseObjectBase::getDefaultOStream();
  if (this->stepperPL_->isSublist("Newmark Parameters")) {
    Teuchos::ParameterList &newmarkPL =
      this->stepperPL_->sublist("Newmark Parameters", true);
    std::string scheme_name = newmarkPL.get("Scheme Name", "Not Specified");
    if (scheme_name == "Not Specified") {
      beta_ = newmarkPL.get("Beta", 0.25);
      gamma_ = newmarkPL.get("Gamma", 0.5);
      TEUCHOS_TEST_FOR_EXCEPTION( (beta_ > 1.0) || (beta_ < 0.0),
        std::logic_error,
        "\nError in 'Newmark Implicit a-Form' stepper: invalid value of Beta = "
        << beta_ << ".  Please select Beta >= 0 and <= 1. \n");
      TEUCHOS_TEST_FOR_EXCEPTION( (gamma_ > 1.0) || (gamma_ < 0.0),
        std::logic_error,
        "\nError in 'Newmark Implicit a-Form' stepper: invalid value of Gamma ="
        <<gamma_ << ".  Please select Gamma >= 0 and <= 1. \n");
      *out_ << "\nSetting Beta = " << beta_ << " and Gamma = " << gamma_
            << " from Newmark Parameters in input file.\n";
    }
    else {
      *out_ << "\nScheme Name = " << scheme_name << ".  Using values \n"
            << "of Beta and Gamma for this scheme (ignoring values of "
            << "Beta and Gamma \n"
            << "in input file, if provided).\n";
       if (scheme_name == "Average Acceleration") {
         beta_ = 0.25; gamma_ = 0.5;
       }
       else if (scheme_name == "Linear Acceleration") {
         beta_ = 0.25; gamma_ = 1.0/6.0;
       }
       else if (scheme_name == "Central Difference") {
         beta_ = 0.0; gamma_ = 0.5;
       }
       else {
         TEUCHOS_TEST_FOR_EXCEPTION(true,
            std::logic_error,
            "\nError in Tempus::StepperNewmarkImplicitAForm!  "
            <<"Invalid Scheme Name = " << scheme_name <<".  \n"
            <<"Valid Scheme Names are: 'Average Acceleration', "
            <<"'Linear Acceleration', \n"
            <<"'Central Difference' and 'Not Specified'.\n");
       }
       *out_ << "===> Beta = " << beta_ << ", Gamma = " << gamma_ << "\n";
    }
    if (beta_ == 0.0) {
      *out_ << "\nWARNING: Running (implicit implementation of) Newmark "
            << "Implicit a-Form Stepper with Beta = 0.0, which \n"
            << "specifies an explicit scheme.  Mass lumping is not possible, "
            << "so this will be slow!  To run explicit \n"
            << "implementation of Newmark Implicit a-Form Stepper, please "
            << "re-run with 'Stepper Type' = 'Newmark Explicit a-Form'.\n"
            << "This stepper allows for mass lumping when called through "
            << "Piro::TempusSolver.\n";
    }
  }
  else {
    *out_ << "\nNo Newmark Parameters sublist found in input file; using "
          << "default values of Beta = "
          << beta_ << " and Gamma = " << gamma_ << ".\n";
  }
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperNewmarkImplicitAForm<Scalar>::getValidParameters() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set("Stepper Type", this->description());
  pl->set("Zero Initial Guess", false);
  pl->set("Solver Name", "",
          "Name of ParameterList containing the solver specifications.");

  return pl;
}
template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperNewmarkImplicitAForm<Scalar>::getDefaultParameters() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set<std::string>("Stepper Type", this->description());
  pl->set<bool>       ("Zero Initial Guess", false);
  pl->set<std::string>("Solver Name", "Default Solver");

  RCP<ParameterList> solverPL = this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperNewmarkImplicitAForm<Scalar>::getNonconstParameterList()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  return(this->stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperNewmarkImplicitAForm<Scalar>::unsetParameterList()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = this->stepperPL_;
  this->stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperNewmarkImplicitAForm_impl_hpp
