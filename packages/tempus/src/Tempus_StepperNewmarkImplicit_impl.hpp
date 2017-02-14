// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkImplicit_impl_hpp
#define Tempus_StepperNewmarkImplicit_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


//IKT, FIXME: I don't think the following routine is needed here.
template <typename Scalar>
std::function<void (const Thyra::VectorBase<Scalar> &,
                          Thyra::VectorBase<Scalar> &)>
StepperNewmarkImplicit<Scalar>::xDotFunction(
  Scalar dt, Teuchos::RCP<const Thyra::VectorBase<Scalar> > x_old)
{
  return [=](const Thyra::VectorBase<Scalar> & x,
                   Thyra::VectorBase<Scalar> & x_dot)
    {
      // this is the Euler x dot vector
      Thyra::V_StVpStV(ptrFromRef(x_dot),1.0/dt,x,-1.0/dt,*x_old);
    };
}

template<class Scalar>
void StepperNewmarkImplicit<Scalar>::predictVelocity(Thyra::VectorBase<Scalar>& vPred,
                                                 const Thyra::VectorBase<Scalar>& v,
                                                 const Thyra::VectorBase<Scalar>& a,
                                                 const Scalar dt) const
{
  //vPred = v + dt*(1.0-gamma_)*a 
  Thyra::V_StVpStV(ptrFromRef(vPred), 1.0, v, dt*(1.0-gamma_), a);  
}
    
template<class Scalar>
void StepperNewmarkImplicit<Scalar>::predictDisplacement(Thyra::VectorBase<Scalar>& dPred,
                                                   const Thyra::VectorBase<Scalar>& d,
                                                   const Thyra::VectorBase<Scalar>& v,
                                                   const Thyra::VectorBase<Scalar>& a,
                                                   const Scalar dt) const
{
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > tmp = Thyra::createMember<Scalar>(dPred.space());
  //dPred = dt*v + dt*dt/2.0*(1.0-2.0*beta_)*a 
  Scalar aConst = dt*dt/2.0*(1.0-2.0*beta_); 
  Thyra::V_StVpStV(ptrFromRef(dPred), dt, v, aConst, a); 
  //dPred += d; 
  Thyra::Vp_V(ptrFromRef(dPred), d, 1.0); 
}
    
template<class Scalar>
void StepperNewmarkImplicit<Scalar>::correctVelocity(Thyra::VectorBase<Scalar>& v,
                                                 const Thyra::VectorBase<Scalar>& vPred,
                                                 const Thyra::VectorBase<Scalar>& a,
                                                 const Scalar dt) const
{
  //v = vPred + dt*gamma_*a 
  Thyra::V_StVpStV(ptrFromRef(v), 1.0, vPred, dt*gamma_, a);  
}

template<class Scalar>
void StepperNewmarkImplicit<Scalar>::correctDisplacement(Thyra::VectorBase<Scalar>& d,
                                                   const Thyra::VectorBase<Scalar>& dPred,
                                                   const Thyra::VectorBase<Scalar>& a,
                                                   const Scalar dt) const
{
  //d = dPred + beta_*dt*dt*a 
  Thyra::V_StVpStV(ptrFromRef(d), 1.0, dPred, beta_*dt*dt, a); 
}


// StepperNewmarkImplicit definitions:
template<class Scalar>
StepperNewmarkImplicit<Scalar>::StepperNewmarkImplicit(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel )
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Set all the input parameters and call initialize
  this->setParameterList(pList);
  this->setModel(transientModel);
  this->initialize();
}


template<class Scalar>
void StepperNewmarkImplicit<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel)
{
  this->validImplicitSecondOrderODE_DAE(transientModel);
  if (residualModel_ != Teuchos::null) residualModel_ = Teuchos::null;
  //IKT, 1/25/16: it is important to set 2nd argument of SecondOrderResidualModelEvaluator to true! 
  //This tells the code that it is a second order time integration scheme, and therefore 
  //to use xdotdot support.
  residualModel_ =
    Teuchos::rcp(new SecondOrderResidualModelEvaluator<Scalar>(transientModel));

  inArgs_  = residualModel_->getNominalValues();
  outArgs_ = residualModel_->createOutArgs();
}


template<class Scalar>
void StepperNewmarkImplicit<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel)
{
  this->setModel(transientModel);
}


/** \brief Set the solver to a pre-defined solver in the ParameterList.
 *  The solver is set to solverName sublist in the Stepper's ParameterList.
 *  The solverName sublist should already be defined in the Stepper's
 *  ParameterList.  Otherwise it will fail.
 */
template<class Scalar>
void StepperNewmarkImplicit<Scalar>::setSolver(std::string solverName)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> solverPL = Teuchos::sublist(pList_, solverName, true);
  pList_->set("Solver Name", solverName);
  if (solver_ != Teuchos::null) solver_ = Teuchos::null;
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
void StepperNewmarkImplicit<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  std::string solverName = pList_->get<std::string>("Solver Name");
  if (is_null(solverPL)) {
    solverPL = Teuchos::sublist(pList_, solverName, true);
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( solverName == solverPL->name(),
      std::logic_error,
         "Error - Trying to add a solver that is already in ParameterList!\n"
      << "  Stepper Type = "<< pList_->get<std::string>("Stepper Type") << "\n"
      << "  Solver Name  = "<<solverName<<"\n");
    solverName = solverPL->name();
    pList_->set("Solver Name", solverName);
    pList_->set(solverName, solverPL);      // Add sublist
  }
  if (solver_ != Teuchos::null) solver_ = Teuchos::null;
  solver_ = rcp(new Thyra::NOXNonlinearSolver());
  RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
  solver_->setParameterList(noxPL);
}


/** \brief Set the predictor to a pre-defined predictor in the ParameterList.
 *  The predictor is set to predictorName sublist in the Stepper's
 *  ParameterList.  The predictorName sublist should already be defined
 *  in the Stepper's ParameterList.  Otherwise it will fail.
 */
template<class Scalar>
void StepperNewmarkImplicit<Scalar>::setPredictor(std::string predictorName)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> predPL = Teuchos::sublist(pList_, predictorName, true);
  pList_->set("Predictor Name", predictorName);
  if (predictorStepper_ != Teuchos::null) predictorStepper_ = Teuchos::null;
  RCP<StepperFactory<Scalar> > sf = Teuchos::rcp(new StepperFactory<Scalar>());
}


/** \brief Set the predictor to the supplied Parameter sublist.
 *  This adds a new predictor Parameter sublist to the Stepper's ParameterList.
 *  If the predictor sublist is null, it tests if the predictor is set in
 *  the Stepper's ParameterList.
 */
template<class Scalar>
void StepperNewmarkImplicit<Scalar>::setPredictor(
  Teuchos::RCP<Teuchos::ParameterList> predPL)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  std::string predictorName = pList_->get<std::string>("Predictor Name","None");
  if (is_null(predPL)) {
    if (predictorName != "None") {
      RCP<ParameterList> predPL = Teuchos::sublist(pList_, predictorName, true);
      RCP<StepperFactory<Scalar> > sf =
        Teuchos::rcp(new StepperFactory<Scalar>());
      predictorStepper_ =
        sf->createStepper(predPL, residualModel_->getTransientModel());
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( predictorName == predPL->name(),
      std::logic_error,
         "Error - Trying to add a predictor that is already in ParameterList!\n"
      << "  Stepper Type = "<< pList_->get<std::string>("Stepper Type") << "\n"
      << "  Predictor Name  = "<<predictorName<<"\n");
    predictorName = predPL->name();
    pList_->set("Predictor Name", predictorName);
    pList_->set(predictorName, predPL);           // Add sublist
    if (predictorStepper_ != Teuchos::null) predictorStepper_ = Teuchos::null;
    RCP<StepperFactory<Scalar> > sf =
      Teuchos::rcp(new StepperFactory<Scalar>());
    predictorStepper_ =
      sf->createStepper(predPL, residualModel_->getTransientModel());
  }
}


template<class Scalar>
void StepperNewmarkImplicit<Scalar>::initialize()
{
  this->setSolver();
  this->setPredictor();
}


template<class Scalar>
void StepperNewmarkImplicit<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;
  
  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperBackardEuler::takeStep()");
  { 
    RCP<SolutionState<Scalar> > workingState = solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState = solutionHistory->getCurrentState();

    //IKT, FIXME, 1/25/17: replace the following code, cut/paste from Backward Euler, 
    //with code for Newmark Beta Implicit scheme.
   
    RCP<const Thyra::VectorBase<Scalar> > xOld = currentState->getX();
    RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > xDot = workingState->getXDot();
    
    computePredictor(solutionHistory);
    if (workingState->stepperState_->stepperStatus_ == Status::FAILED) return;
    
    //typedef Thyra::ModelEvaluatorBase MEB;
    const Scalar time = workingState->getTime();
    const Scalar dt   = workingState->getTimeStep();
    
    //IKT, compute predictor disp and vel using predictDisplacement
    //and predictVelocity routines 

    // constant variable capture of xOld pointer
    auto computeXDot = xDotFunction(dt, xOld);
    
    //IKT, FIXME: replace the following with values for NewmarkImplicit-Beta. 
    //See Piro::Epetra::NewmarkImplicitDecorator::evalModel 
    Scalar alpha = 1.0/dt;
    Scalar beta = 1.0;
    Scalar omega = 0.0; //FIXME set to correct value 
    Scalar t = time+dt;
    
    //IKT, FIXME: modify initialize to take in disp_pred and vel_pred 
    //(predictors), if possible.  We don't want to compute disp_pred 
    //and vel_pred in SecondOrderResidualModelEvaluator::evalModelImpl, but rather save it 
    //and use it there.  Will not need to pass in computeXDot. 

    residualModel_->initialize(computeXDot, t, alpha, beta, omega);
   
    //IKT, FIXME: x will need to be replaced below, I believe with acceleration. 
    const Thyra::SolveStatus<double> sStatus =
      this->solveNonLinear(residualModel_, *solver_, x, inArgs_);
   
    //IKT, FIXME: replace the following with computation of corrector 
    //disp and vel using correctDisplacement and correctVelocity
    //routines. 

    computeXDot(*x, *xDot);
    
    if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED )
      workingState->stepperState_->stepperStatus_ = Status::PASSED;
    else
      workingState->stepperState_->stepperStatus_ = Status::FAILED;
  }
  
  return;
}

template<class Scalar>
void StepperNewmarkImplicit<Scalar>::computePredictor(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  if (predictorStepper_ == Teuchos::null) return;
  predictorStepper_->takeStep(solutionHistory);

  Status & stepperStatus =
    solutionHistory->getWorkingState()->stepperState_->stepperStatus_;

  if (stepperStatus == Status::FAILED) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperNewmarkImplicit::computePredictor");
    *out << "Warning - predictorStepper has failed." << std::endl;
  } else {
    // Reset status to WORKING since this is the predictor
    stepperStatus = Status::WORKING;
  }
}


/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperNewmarkImplicit<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperNewmarkImplicit<Scalar>::description() const
{
  std::string name = "Backward Euler";
  return(name);
}


template<class Scalar>
void StepperNewmarkImplicit<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "residualModel_ = " << residualModel_->description() << std::endl;
}


template <class Scalar>
void StepperNewmarkImplicit<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pList));
  //pList->validateParameters(*this->getValidParameters());
  pList_ = pList;
  //Get beta and gamma from parameter list 
  //IKT, FIXME: does parameter list get validated somewhere?  validateParameters above is commented out...

  std::string stepperType = pList_->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Newmark Beta Implicit",
    std::logic_error,
       "Error - Stepper Type is not 'Newmark Beta Implicit'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");
  Teuchos::readVerboseObjectSublist(&*pList_,this);
  beta_ = 0.25; //default value
  gamma_ = 0.5; //default value
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream(); 
  if (pList_->isSublist("Newmark Beta Parameters")) {
    Teuchos::ParameterList &newmarkPL = pList_->sublist("Newmark Beta Parameters", true);
    beta_ = newmarkPL.get("Beta", 0.25); 
    gamma_ = newmarkPL.get("Gamma", 0.5);
    if (gamma_ == 0.0) {
      TEUCHOS_TEST_FOR_EXCEPTION( true,
        std::logic_error,
           "Error - Stepper Type = Newmark Beta Implicit is not value with Gamma = 0.0, as this \n"
            << " value specifies an explicit scheme.  Please run with Gamma > 0.0.  Explicit Newmark Beta \n"
            << " scheme is not yet available in Tempus. \n");
    }
    *out << "\n \nSetting Beta = " << beta_ << " and Gamma = " << gamma_ << " from Newmark Beta Parameters in input file.\n\n";  
  }
  else {
    *out << "\n  \nNo Newmark Beta Parameters sublist found in input file; using default values of Beta = " 
         << beta_ << " and Gamma = " << gamma_ << ".\n\n"; 
  }
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperNewmarkImplicit<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    std::ostringstream tmp;
    pl->set("Stepper Type", this->description());
    pl->set("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");
    validPL = pl;
  }
  return validPL;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperNewmarkImplicit<Scalar>::getNonconstParameterList()
{
  return(pList_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperNewmarkImplicit<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = pList_;
  pList_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperNewmarkImplicit_impl_hpp
