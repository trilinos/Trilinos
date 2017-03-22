// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

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
      Thyra::V_StVpStV(Teuchos::ptrFromRef(x_dot),1.0/dt,x,-1.0/dt,*x_old);
    };
}


// StepperBackwardEuler definitions:
template<class Scalar>
StepperBackwardEuler<Scalar>::StepperBackwardEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel,
  Teuchos::RCP<Teuchos::ParameterList> pList)
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
    // Create default solver, otherwise keep current solver.
    if (solver_ == Teuchos::null) {
      solverPL = Teuchos::sublist(pList_, solverName, true);
      solver_ = rcp(new Thyra::NOXNonlinearSolver());
      RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
      solver_->setParameterList(noxPL);
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( solverName == solverPL->name(),
      std::logic_error,
         "Error - Trying to add a solver that is already in ParameterList!\n"
      << "  Stepper Type = "<< pList_->get<std::string>("Stepper Type") << "\n"
      << "  Solver Name  = "<<solverName<<"\n");
    solverName = solverPL->name();
    pList_->set("Solver Name", solverName);
    pList_->set(solverName, solverPL);      // Add sublist
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
void StepperBackwardEuler<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> solverPL = solver->getNonconstParameterList();
  std::string solverName = solverPL->name();
  pList_->set("Solver Name", solverName);
  pList_->set(solverName, solverPL);      // Add sublist
  solver_ = solver;
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
        sf->createStepper(residualModel_->getTransientModel(), predPL);
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
      sf->createStepper(residualModel_->getTransientModel(), predPL);
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

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperBackardEuler::takeStep()");
  {
    RCP<SolutionState<Scalar> > workingState = solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState = solutionHistory->getCurrentState();

    RCP<const Thyra::VectorBase<Scalar> > xOld = currentState->getX();
    RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > xDot = workingState->getXDot();

    computePredictor(solutionHistory);
    if (workingState->getStepperState()->stepperStatus_ == Status::FAILED)
      return;

    //typedef Thyra::ModelEvaluatorBase MEB;
    const Scalar time = workingState->getTime();
    const Scalar dt   = workingState->getTimeStep();

    // constant variable capture of xOld pointer
    auto computeXDot = xDotFunction(dt, xOld);

    Scalar alpha = 1.0/dt;
    Scalar beta = 1.0;
    Scalar t = time+dt;

    residualModel_->initialize(computeXDot, t, alpha, beta);

    const Thyra::SolveStatus<double> sStatus =
      this->solveNonLinear(residualModel_, *solver_, x, inArgs_);

    computeXDot(*x, *xDot);

    if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED )
      workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    else
      workingState->getStepperState()->stepperStatus_ = Status::FAILED;
  }

  return;
}

template<class Scalar>
void StepperBackwardEuler<Scalar>::computePredictor(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  if (predictorStepper_ == Teuchos::null) return;
  predictorStepper_->takeStep(solutionHistory);

  Status & stepperStatus =
    solutionHistory->getWorkingState()->getStepperState()->stepperStatus_;

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
  if (pList == Teuchos::null) pList_ = this->getDefaultParameters();
  else pList_ = pList;
  // Can not validate because of optional Parameters (e.g., Solver Name).
  //pList_->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string stepperType = pList_->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Backward Euler",
    std::logic_error,
       "Error - Stepper Type is not 'Backward Euler'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set("Stepper Type", this->description());
  pl->set("Solver Name", "",
    "Name of ParameterList containing the solver specifications.");

  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::getDefaultParameters() const
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set<std::string>("Stepper Type", this->description());
  pl->set<std::string>("Solver Name", "Default Solver");
  pl->set<std::string>("Predictor Name", "Default Predictor");

  RCP<ParameterList> solverPL = this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  // Predictor ParameterList
  RCP<ParameterList> predPL = Teuchos::parameterList();
  predPL->set("Stepper Type", "Forward Euler");
  pl->set("Default Predictor", *predPL);

  return pl;
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
