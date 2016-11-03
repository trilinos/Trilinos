#ifndef Tempus_StepperBackwardEuler_impl_hpp
#define Tempus_StepperBackwardEuler_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


template <typename Scalar>
std::function<void (const Thyra::VectorBase<Scalar> &,
                          Thyra::VectorBase<Scalar> &)>
StepperBackwardEuler<Scalar>::xDotFunction(
  Scalar dt, Teuchos::RCP<const Thyra::VectorBase<Scalar> > x_old)
{
  return [=](const Thyra::VectorBase<Scalar> & x,
                   Thyra::VectorBase<Scalar> & x_dot)
    {
      // this is the Euler x dot vector
      Thyra::V_StVpStV(ptrFromRef(x_dot),1.0/dt,x,-1.0/dt,*x_old);
    };
}


// StepperBackwardEuler definitions:
template<class Scalar>
StepperBackwardEuler<Scalar>::StepperBackwardEuler(
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
void StepperBackwardEuler<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel)
{
  this->validImplicitODE_DAE(transientModel);
  if (residualModel_ != Teuchos::null) residualModel_ = Teuchos::null;
  residualModel_ =
    Teuchos::rcp(new ResidualModelEvaluator<Scalar>(transientModel));

  inArgs_  = residualModel_->getNominalValues();
  outArgs_ = residualModel_->createOutArgs();
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::setNonConstModel(
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
void StepperBackwardEuler<Scalar>::setSolver(std::string solverName)
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
void StepperBackwardEuler<Scalar>::setSolver(
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
void StepperBackwardEuler<Scalar>::setPredictor(std::string predictorName)
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
void StepperBackwardEuler<Scalar>::setPredictor(
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
void StepperBackwardEuler<Scalar>::initialize()
{
  this->setSolver();
  this->setPredictor();
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  RCP<SolutionState<Scalar> > workingState = solutionHistory->getWorkingState();
  RCP<SolutionState<Scalar> > currentState = solutionHistory->getCurrentState();

  RCP<const Thyra::VectorBase<Scalar> > xOld = currentState->getX();
  RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();
  RCP<Thyra::VectorBase<Scalar> > xDot = workingState->getXDot();

  computePredictor(solutionHistory);
  if (workingState->stepperState_->stepperStatus_ == Status::FAILED) return;

  //typedef Thyra::ModelEvaluatorBase MEB;
  const Scalar time = workingState->getTime();
  const Scalar dt   = workingState->getTimeStep();

  // constant variable capture of xOld pointer
  auto computeXDot = xDotFunction(dt, xOld);

  Scalar alpha = 1.0/dt;
  Scalar beta = 1.0;
  Scalar t = time+dt;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> base_point =
    residualModel_->getNominalValues();
  residualModel_->initialize(computeXDot, t, alpha, beta, base_point);

  const Thyra::SolveStatus<double> solve_status =
    this->solveNonLinear(residualModel_, *solver_, x, inArgs_);

  computeXDot(*x, *xDot);

  if (solve_status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED )
    workingState->stepperState_->stepperStatus_ = Status::PASSED;
  else
    workingState->stepperState_->stepperStatus_ = Status::FAILED;

  return;
}

template<class Scalar>
void StepperBackwardEuler<Scalar>::computePredictor(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  if (predictorStepper_ == Teuchos::null) return;
  predictorStepper_->takeStep(solutionHistory);

  Status & stepperStatus =
    solutionHistory->getWorkingState()->stepperState_->stepperStatus_;

  if (stepperStatus == Status::FAILED) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperBackwardEuler::computePredictor");
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
StepperBackwardEuler<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperBackwardEuler<Scalar>::description() const
{
  std::string name = "Backward Euler";
  return(name);
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "residualModel_ = " << residualModel_->description() << std::endl;
}


template <class Scalar>
void StepperBackwardEuler<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pList));
  //pList->validateParameters(*this->getValidParameters());
  pList_ = pList;

  std::string stepperType = pList_->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Backward Euler",
    std::logic_error,
       "Error - Stepper Type is not 'Backward Euler'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");

  std::string solverName = pList_->get<std::string>("Solver Name");

  Teuchos::readVerboseObjectSublist(&*pList_,this);
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    std::ostringstream tmp;
    tmp << "'Stepper Type' must be 'Backward Euler'.";
    pl->set("Stepper Type", "Backward Euler", tmp.str());

    validPL = pl;
  }
  return validPL;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::getNonconstParameterList()
{
  return(pList_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = pList_;
  pList_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperBackwardEuler_impl_hpp
