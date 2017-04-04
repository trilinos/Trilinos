// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperHHTAlpha_impl_hpp
#define Tempus_StepperHHTAlpha_impl_hpp

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
void StepperHHTAlpha<Scalar>::
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
void StepperHHTAlpha<Scalar>::
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
void StepperHHTAlpha<Scalar>::
predictVelocity_alpha_f(Thyra::VectorBase<Scalar>& vPred,
                        const Thyra::VectorBase<Scalar>& v) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //vPred = (1-alpha_f)*vPred + alpha_f*v
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0-alpha_f_, vPred, alpha_f_, v);
}


template<class Scalar>
void StepperHHTAlpha<Scalar>::
predictDisplacement_alpha_f(Thyra::VectorBase<Scalar>& dPred,
                            const Thyra::VectorBase<Scalar>& d) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //dPred = (1-alpha_f)*dPred + alpha_f*d
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), 1.0-alpha_f_, dPred, alpha_f_, d);
}

template<class Scalar>
void StepperHHTAlpha<Scalar>::
correctAcceleration(Thyra::VectorBase<Scalar>& a_n_plus1,
                    const Thyra::VectorBase<Scalar>& a_n) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Scalar c = 1.0/(1.0-alpha_m_); 
  //a_n_plus1 = 1.0/(1.0-alpha_m_)*a_n_plus1 - alpha_m/(1.0-alpha_m)*a_n = (1-alpha_f)*vPred + alpha_f*v
  Thyra::V_StVpStV(Teuchos::ptrFromRef(a_n_plus1), c, a_n_plus1, -c*alpha_m_, a_n);
}



template<class Scalar>
void StepperHHTAlpha<Scalar>::
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
void StepperHHTAlpha<Scalar>::
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


// StepperHHTAlpha definitions:
template<class Scalar>
StepperHHTAlpha<Scalar>::StepperHHTAlpha(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel,
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
  this->setModel(transientModel);
  this->initialize();
}


template<class Scalar>
void StepperHHTAlpha<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->validSecondOrderODE_DAE(transientModel);
  if (residualModel_ != Teuchos::null) residualModel_ = Teuchos::null;
  residualModel_ =
    Teuchos::rcp(new SecondOrderResidualModelEvaluator<Scalar>(transientModel,
                                                      "HHT-Alpha"));
  inArgs_  = residualModel_->getNominalValues();
  outArgs_ = residualModel_->createOutArgs();
}


template<class Scalar>
void StepperHHTAlpha<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->setModel(transientModel);
}


/** \brief Set the solver to a pre-defined solver in the ParameterList.
 *  The solver is set to solverName sublist in the Stepper's ParameterList.
 *  The solverName sublist should already be defined in the Stepper's
 *  ParameterList.  Otherwise it will fail.
 */
template<class Scalar>
void StepperHHTAlpha<Scalar>::setSolver(std::string solverName)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> solverPL = Teuchos::sublist(stepperPL_, solverName, true);
  stepperPL_->set("Solver Name", solverName);
  solver_ = rcp(new Thyra::NOXNonlinearSolver());
  RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
  solver_->setParameterList(noxPL);
}


/** \brief Set the solver to the supplied Parameter sublist.
 *  This adds a new solver Parameter sublist to the Stepper's ParameterList.
 *  If the solver sublist is null, the solver is set to the solver name
 *  in the Stepper's ParameterList.
 */
template<class Scalar>
void StepperHHTAlpha<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  std::string solverName = stepperPL_->get<std::string>("Solver Name");
  if (is_null(solverPL)) {
    // Create default solver, otherwise keep current solver.
    if (solver_ == Teuchos::null) {
      solverPL = Teuchos::sublist(stepperPL_, solverName, true);
      solver_ = rcp(new Thyra::NOXNonlinearSolver());
      RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
      solver_->setParameterList(noxPL);
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( solverName == solverPL->name(),
      std::logic_error,
         "Error - Trying to add a solver that is already in ParameterList!\n"
      << "  Stepper Type = "<< stepperPL_->get<std::string>("Stepper Type")
      << "\n" << "  Solver Name  = "<<solverName<<"\n");
    solverName = solverPL->name();
    stepperPL_->set("Solver Name", solverName);
    stepperPL_->set(solverName, solverPL);      // Add sublist
    solver_ = rcp(new Thyra::NOXNonlinearSolver());
    RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
    solver_->setParameterList(noxPL);
  }
}


/** \brief Set the solver.
 *  This sets the solver to supplied solver and adds solver's ParameterList
 *  to the Stepper ParameterList.
 */
template<class Scalar>
void StepperHHTAlpha<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> solverPL = solver->getNonconstParameterList();
  std::string solverName = solverPL->name();
  stepperPL_->set("Solver Name", solverName);
  stepperPL_->set(solverName, solverPL);      // Add sublist
  solver_ = solver;
}


template<class Scalar>
void StepperHHTAlpha<Scalar>::initialize()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->setSolver();
}


template<class Scalar>
void StepperHHTAlpha<Scalar>::takeStep(
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
    const Scalar time = workingState->getTime();
    const Scalar dt   = workingState->getTimeStep();
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
      residualModel_->initializeNewmark(a_init,v_init,d_init,0.0,time,beta_,gamma_);
      const Thyra::SolveStatus<double> sStatus =
        this->solveNonLinear(residualModel_, *solver_, a_init, inArgs_);
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

    //compute second displacement and velocity predictors (those that are functions of alpha_f) 
    predictDisplacement_alpha_f(*d_pred, *d_old); 
    predictVelocity_alpha_f(*v_pred, *v_old);  

    //inject d_pred, v_pred, a and other relevant data into residualModel_
    residualModel_->initializeNewmark(a_old,v_pred,d_pred,dt,t,beta_,gamma_);

    //Solve for new acceleration
    const Thyra::SolveStatus<double> sStatus =
      this->solveNonLinear(residualModel_, *solver_, a_new, inArgs_);

    //correct acceleration (function of alpha_m) 
    correctAcceleration(*a_new, *a_old); 

    //correct velocity and displacement 
    correctVelocity(*v_new, *v_pred, *a_new, dt);
    correctDisplacement(*d_new, *d_pred, *a_new, dt);

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
StepperHHTAlpha<Scalar>::
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
std::string StepperHHTAlpha<Scalar>::description() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  std::string name = "HHTAlpha ";
  return(name);
}


template<class Scalar>
void StepperHHTAlpha<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  out << description() << "::describe:" << std::endl
      << "residualModel_ = " << residualModel_->description() << std::endl;
}


template <class Scalar>
void StepperHHTAlpha<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  if (pList == Teuchos::null) stepperPL_ = this->getDefaultParameters();
  else stepperPL_ = pList;
  // Can not validate because of optional Parameters.
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
  //Get beta and gamma from parameter list
  //IKT, FIXME: does parameter list get validated somewhere?  validateParameters above is commented out...

  std::string stepperType = stepperPL_->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "HHT-Alpha",
    std::logic_error,
       "\nError - Stepper Type is not 'HHT-Alpha'!\n"
    << "Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");
  beta_ = 0.25; //default value
  gamma_ = 0.5; //default value
  //IKT, FIXME: test this scheme for alpha_f and alpha_m != 0.0. 
  //Once that is done, logic should be changed to allow user to select these options from 
  //the parameter list. 
  alpha_f_ = 0.0; //default value.  Hard-coded for Newmark-Beta for now.
  alpha_m_ = 0.0; //default value.  Hard-coded for Newmark-Beta for now.
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  if (stepperPL_->isSublist("HHT-Alpha Parameters")) {
    Teuchos::ParameterList &HHTalphaPL =
      stepperPL_->sublist("HHT-Alpha Parameters", true);
    std::string scheme_name = HHTalphaPL.get("Scheme Name", "Not Specified");
    alpha_m_ = HHTalphaPL.get("Alpha_m", 0.0);
    alpha_f_ = HHTalphaPL.get("Alpha_f", 0.0);
    TEUCHOS_TEST_FOR_EXCEPTION( (alpha_m_ >= 1.0) || (alpha_m_ < 0.0),
      std::logic_error,
         "\nError in 'HHT-Alpha' stepper: invalid value of Alpha_m = " << alpha_m_ << ".  Please select Alpha_m >= 0 and < 1. \n"); 
    TEUCHOS_TEST_FOR_EXCEPTION( (alpha_f_ > 1.0) || (alpha_f_ < 0.0),
      std::logic_error,
         "\nError in 'HHT-Alpha' stepper: invalid value of Alpha_f = " << alpha_f_ << ".  Please select Alpha_f >= 0 and <= 1. \n"); 
    TEUCHOS_TEST_FOR_EXCEPTION( (alpha_m_ != 0.0) || (alpha_f_ != 0.0),
      std::logic_error,
         "\nError - 'HHT-Alpha' stepper has not been verified yet for Alpha_m, Alpha_f != 0! \n"
          << "You have specified Alpha_m = " << alpha_m_ << ", Alpha_f = " << alpha_f_ <<".\n");
    if (scheme_name == "Newmark Beta") {
      beta_ = HHTalphaPL.get("Beta", 0.25);
      gamma_ = HHTalphaPL.get("Gamma", 0.5);
      TEUCHOS_TEST_FOR_EXCEPTION( (beta_ > 1.0) || (beta_ < 0.0),
        std::logic_error,
           "\nError in 'HHT-Alpha' stepper: invalid value of Beta = " << beta_ << ".  Please select Beta >= 0 and <= 1. \n"); 
      TEUCHOS_TEST_FOR_EXCEPTION( (gamma_ > 1.0) || (gamma_ < 0.0),
        std::logic_error,
           "\nError in 'HHT-Alpha' stepper: invalid value of Gamma = " <<gamma_ << ".  Please select Gamma >= 0 and <= 1. \n"); 
      *out << "\n \nScheme Name = Newmark Beta.  Setting Alpha_f = Alpha_m = 0. Setting \n"
           << "Beta = " << beta_ << " and Gamma = " << gamma_ << " from HHT-Alpha Parameters in input file.\n\n";
    } 
    else {
      *out << "\n \nScheme Name = " << scheme_name << ".  Using values of Alpha_m, Alpha_f, \n"
           << "Beta and Gamma for this scheme (ignoring values of Alpha_m, Alpha_f, Beta and Gamma \n"
           << "in input file, if provided).\n"; 
       if (scheme_name == "Newmark Beta Average Acceleration") {
         beta_ = 0.25; gamma_ = 0.5; 
       }
       else if (scheme_name == "Newmark Beta Linear Acceleration") {
         beta_ = 0.25; gamma_ = 1.0/6.0; 
       }
       else if (scheme_name == "Newmark Beta Central Difference") {
         beta_ = 0.0; gamma_ = 0.5; 
       }
       else {
         TEUCHOS_TEST_FOR_EXCEPTION(true,
            std::logic_error,
            "\nError in Tempus::StepperHHTAlpha!  Invalid Scheme Name = " << scheme_name <<".  \n" 
            <<"Valid Scheme Names are: 'Newmark Beta', 'Newmark Beta Average Acceleration', \n"
            <<"'Newmark Beta Linear Acceleration', and 'Newmark Beta Central Difference'.\n"); 
       }
       *out << "===> Alpha_m = " << alpha_m_ << ", Alpha_f = " << alpha_f_ << ", Beta = " << beta_ << ", Gamma = " << gamma_ << "\n"; 
    }
    if (beta_ == 0.0) {
      *out << "\n \nRunning  HHT-Alpha Stepper with Beta = 0.0, which \n"
           << "specifies an explicit scheme.  WARNING: code has not been optimized \n"
           << "yet for this case (no mass lumping)\n"; 
    }
  }
  else {
    *out << "\n  \nNo HHT-Alpha Parameters sublist found in input file; using default values of Beta = "
         << beta_ << " and Gamma = " << gamma_ << ".\n\n";
  }
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperHHTAlpha<Scalar>::getValidParameters() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set("Stepper Type", this->description());
  pl->set("Solver Name", "",
          "Name of ParameterList containing the solver specifications.");

  return pl;
}
template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperHHTAlpha<Scalar>::getDefaultParameters() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set<std::string>("Stepper Type", this->description());
  pl->set<std::string>("Solver Name", "Default Solver");

  RCP<ParameterList> solverPL = this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperHHTAlpha<Scalar>::getNonconstParameterList()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  return(stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperHHTAlpha<Scalar>::unsetParameterList()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperHHTAlpha_impl_hpp
