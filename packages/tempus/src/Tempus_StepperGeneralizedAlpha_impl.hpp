// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperGeneralizedAlpha_impl_hpp
#define Tempus_StepperGeneralizedAlpha_impl_hpp

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
void StepperGeneralizedAlpha<Scalar>::
prediction(Thyra::VectorBase<Scalar>& dPred,
		   Thyra::VectorBase<Scalar>& vPred,
		   const Thyra::VectorBase<Scalar>& d,
           const Thyra::VectorBase<Scalar>& v,
           const Thyra::VectorBase<Scalar>& a,
           const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //dPred = d_n + dt*v_n + dt*dt/2.0*(1.0-2.0*\beta)*a_n
  Scalar aConst = dt*dt*(0.5-beta_);
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), dt, v, aConst, a);
  Thyra::Vp_V(Teuchos::ptrFromRef(dPred), d, 1.0);
  //dPred = (1-\alpha)*dPred + \alpha*d_n
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), 1.0-alpha_f_, dPred, alpha_f_, d);
  
  //vPred = v_n + dt*(1.0-\gamma)*a_n
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0, v, dt*(1.0-gamma_), a);
  //vPred = (1-\alpha)*vPred + \alpha*v_n
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0-alpha_f_, vPred, alpha_f_, v);
}

template<class Scalar>
void StepperGeneralizedAlpha<Scalar>::
correction(Thyra::VectorBase<Scalar>& d,
                    Thyra::VectorBase<Scalar>& v,
                    const Thyra::VectorBase<Scalar>& a_old,
					const Thyra::VectorBase<Scalar>& a,
                    const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //d_{n+1} = d_n + dt*v_n + dt*dt*( (0.5-\beta)*a_n + \beta*a_{n+1} )
  Thyra::Vp_StV(Teuchos::ptrFromRef(d), dt, v);
  Scalar c = dt*dt*(0.5-beta_);
  Thyra::Vp_StV(Teuchos::ptrFromRef(d), c, a_old);
  c = dt*dt*beta_;
  Thyra::Vp_StV(Teuchos::ptrFromRef(d), c, a);
  
  //v_{n+1} = v_n + dt*(1.0-\gamma)*a_n + dt*\gamma*a_{n+1}
  Thyra::Vp_StV(Teuchos::ptrFromRef(v), dt*(1.0-gamma_), a_old);
  Thyra::Vp_StV(Teuchos::ptrFromRef(v), dt*gamma_, a);
}




template<class Scalar>
StepperGeneralizedAlpha<Scalar>::StepperGeneralizedAlpha(
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
void StepperGeneralizedAlpha<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->validSecondOrderODE_DAE(appModel);
  Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
    Teuchos::rcp(new WrapperModelEvaluatorSecondOrder<Scalar>(appModel,
                                                      "Generalized-Alpha"));
  this->wrapperModel_ = wrapperModel;
}


template<class Scalar>
void StepperGeneralizedAlpha<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( this->wrapperModel_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperGeneralizedAlpha::initialize()\n");

#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->setSolver();
}


template<class Scalar>
void StepperGeneralizedAlpha<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperBackardEuler::takeStep()");
  {
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
      Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar> >(
        this->wrapperModel_);

    //Get values of d, v and a from previous step
    RCP<const Thyra::VectorBase<Scalar> > d_old = currentState->getX();
    RCP<const Thyra::VectorBase<Scalar> > v_old = currentState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > a_old = currentState->getXDotDot();

#ifdef DEBUG_OUTPUT
    //IKT, 3/21/17, debug output: pring d_old, v_old, a_old to check for
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
    const Scalar dt   = workingState->getTimeStep();

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
      const Thyra::SolveStatus<Scalar> sStatus=this->solveImplicitODE(a_init);

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
    prediction(*d_pred, *v_pred, *d_old, *v_old, *a_old, dt);
	//Update time. PA: It is supposed that wrapperModel would fetch external loads at this time point
    Scalar t = time+ (1.0-alpha_f_)*dt;

    //inject d_pred, v_pred, a and other relevant data into wrapperModel
    wrapperModel->initializeNewmark(a_old,v_pred,d_pred,dt,t,beta_,gamma_);

    //Solve for new acceleration
    const Thyra::SolveStatus<Scalar> sStatus = this->solveImplicitODE(a_new);

    //correct velocity and displacement
    correction(*d_new, *v_new, *a_old, *a_new, dt);

    if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED )
      workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    else
      workingState->getStepperState()->stepperStatus_ = Status::FAILED;
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
StepperGeneralizedAlpha<Scalar>::
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
std::string StepperGeneralizedAlpha<Scalar>::description() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  std::string name = "Generalized-Alpha";
  return(name);
}


template<class Scalar>
void StepperGeneralizedAlpha<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  out << description() << "::describe:" << std::endl
      << "wrapperModel_ = " << this->wrapperModel_->description() << std::endl;
}


template <class Scalar>
void StepperGeneralizedAlpha<Scalar>::setParameterList(
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
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Generalized-Alpha",
    std::logic_error,
       "\nError - Stepper Type is not 'Generalized-Alpha'!\n" << "Stepper Type = "
       << stepperPL->get<std::string>("Stepper Type") << "\n");
  alpha_f_ = 0.0; alpha_m_=0.0;  //default value. 
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  if (this->stepperPL_->isSublist("Generalized-Alpha Parameters")) {
    Teuchos::ParameterList &HHTalphaPL =
      this->stepperPL_->sublist("Generalized-Alpha Parameters", true);
    std::string scheme_name = HHTalphaPL.get("Scheme Name", "Not Specified");
    alpha_f_ = HHTalphaPL.get("Alpha_f", 0.0);
	alpha_m_ = HHTalphaPL.get("Alpha_m", 0.0);
    TEUCHOS_TEST_FOR_EXCEPTION( (alpha_f_ <0.0) || (alpha_m_ < 0.0) || 
	  (alpha_f_ >1.0) || (alpha_m_ >1.0) || alpha_f_<alpha_m_,
      std::logic_error,
         "\nError in 'Generalized-Alpha' stepper: invalid value of Alpha_f = "
         << alpha_f_ << "; Alpha_m = " << alpha_m_ 
		 << ".  Please select 0<=Alpha_f<=1, 0<=Alpha_m<=1 and Alpha_f>=Alpha_m \n");
    *out << "\n \nScheme Name = Generalized-Alpha.  Setting Alpha_f = "
           << alpha_f_ << "; Alpha_m= " << alpha_m_ 
           << "\n from Generalized-Alpha Parameters in input file.\n\n";
  }
  else {
    *out << "\n  \nNo Generalized-Alpha Parameters sublist found in input file; "
         << "using default values of Alpha_f=0 and Alpha_m=0 .\n\n";
  }
  gamma_ = 0.5 + alpha_f_ - alpha_m_;
  beta_ = 0.25*(0.5+gamma_)*(0.5+gamma_); 
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperGeneralizedAlpha<Scalar>::getValidParameters() const
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
StepperGeneralizedAlpha<Scalar>::getDefaultParameters() const
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
StepperGeneralizedAlpha<Scalar>::getNonconstParameterList()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  return(this->stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperGeneralizedAlpha<Scalar>::unsetParameterList()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = this->stepperPL_;
  this->stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperGeneralizedAlpha_impl_hpp
