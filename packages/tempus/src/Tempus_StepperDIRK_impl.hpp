// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperDIRK_impl_hpp
#define Tempus_StepperDIRK_impl_hpp

#include "Tempus_RKButcherTableauBuilder.hpp"
#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;

template<class Scalar>
StepperDIRK<Scalar>::StepperDIRK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  std::string stepperType)
{
  this->setTableau(Teuchos::null, stepperType);
  this->setParameterList(Teuchos::null);
  this->setModel(appModel);
  this->setSolver();
  this->setObserver();
  this->initialize();
}

template<class Scalar>
StepperDIRK<Scalar>::StepperDIRK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList>                      pList)
{
  this->setTableau(pList, "SDIRK 2 Stage 2nd order");
  this->setParameterList(pList);
  this->setModel(appModel);
  this->setSolver();
  this->setObserver();
  this->initialize();
}

template<class Scalar>
StepperDIRK<Scalar>::StepperDIRK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  std::string stepperType,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  this->setTableau(pList, stepperType);
  this->setParameterList(pList);
  this->setModel(appModel);
  this->setSolver();
  this->setObserver();
  this->initialize();
}


template<class Scalar>
void StepperDIRK<Scalar>::setTableau(
  Teuchos::RCP<Teuchos::ParameterList> pList,
  std::string stepperType)
{
  if (stepperType == "") {
    if (pList == Teuchos::null)
      stepperType = "SDIRK 2 Stage 2nd order";
    else
      stepperType = pList->get<std::string>("Stepper Type",
                                            "SDIRK 2 Stage 2nd order");
  }

  DIRK_ButcherTableau_ = createRKBT<Scalar>(stepperType,pList);

  TEUCHOS_TEST_FOR_EXCEPTION( DIRK_ButcherTableau_->isDIRK() != true,
    std::logic_error,
       "Error - StepperDIRK did not receive a DIRK Butcher Tableau!\n"
    << "  Stepper Type = " << stepperType <<  "\n");
  description_ = DIRK_ButcherTableau_->description();
}


template<class Scalar>
void StepperDIRK<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (stepperObserver_ == Teuchos::null) {
      stepperDIRKObserver_ =
        Teuchos::rcp(new StepperDIRKObserver<Scalar>());
      stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >
          (stepperDIRKObserver_);
     }
  } else {
    stepperObserver_ = obs;
    stepperDIRKObserver_ =
      Teuchos::rcp_dynamic_cast<StepperDIRKObserver<Scalar> >(stepperObserver_);
  }
}


template<class Scalar>
void StepperDIRK<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( DIRK_ButcherTableau_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the Butcher Tableau, setTableau(), before calling "
    "StepperDIRK::initialize()\n");

  TEUCHOS_TEST_FOR_EXCEPTION( this->wrapperModel_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperDIRK::initialize()\n");

  // Initialize the stage vectors
  const int numStages = DIRK_ButcherTableau_->numStages();
  stageX_    = this->wrapperModel_->getNominalValues().get_x()->clone_v();
  stageXDot_.resize(numStages);
  for (int i=0; i<numStages; ++i) {
    stageXDot_[i] = Thyra::createMember(this->wrapperModel_->get_f_space());
    assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  }
  xTilde_    = Thyra::createMember(this->wrapperModel_->get_x_space());
  assign(xTilde_.ptr(),    Teuchos::ScalarTraits<Scalar>::zero());
}

template<class Scalar>
void StepperDIRK<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperDIRK::takeStep()");
  {
    stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    const int numStages = DIRK_ButcherTableau_->numStages();
    Teuchos::SerialDenseMatrix<int,Scalar> A = DIRK_ButcherTableau_->A();
    Teuchos::SerialDenseVector<int,Scalar> b = DIRK_ButcherTableau_->b();
    Teuchos::SerialDenseVector<int,Scalar> c = DIRK_ButcherTableau_->c();

    // Compute stage solutions
    bool pass = true;
    Thyra::SolveStatus<Scalar> sStatus;
    for (int i=0; i < numStages; ++i) {
      if (!Teuchos::is_null(stepperDIRKObserver_))
        stepperDIRKObserver_->observeBeginStage(solutionHistory, *this);
      Thyra::assign(xTilde_.ptr(), *(currentState->getX()));
      for (int j=0; j < i; ++j) {
        if (A(i,j) != Teuchos::ScalarTraits<Scalar>::zero()) {
          Thyra::Vp_StV(xTilde_.ptr(), dt*A(i,j), *(stageXDot_[j]));
        }
      }

      Scalar ts = time + c(i)*dt;
      if (A(i,i) == Teuchos::ScalarTraits<Scalar>::zero()) {
        // Explicit stage for the ImplicitODE_DAE
        bool isNeeded = false;
        for (int k=i+1; k<numStages; ++k) if (A(k,i) != 0.0) isNeeded = true;
        if (b(i) != 0.0) isNeeded = true;
        if (isNeeded == false) {
          // stageXDot_[i] is not needed.
          assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
        } else {
          typedef Thyra::ModelEvaluatorBase MEB;
          MEB::InArgs<Scalar>  inArgs  = this->wrapperModel_->getInArgs();
          MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
          inArgs.set_x(xTilde_);
          if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(ts);
          if (inArgs.supports(MEB::IN_ARG_x_dot))
            inArgs.set_x_dot(Teuchos::null);
          outArgs.set_f(stageXDot_[i]);

          if (!Teuchos::is_null(stepperDIRKObserver_))
            stepperDIRKObserver_->observeBeforeExplicit(solutionHistory,*this);
          this->wrapperModel_->getAppModel()->evalModel(inArgs,outArgs);
        }
      } else {
        // Implicit stage for the ImplicitODE_DAE
        Scalar alpha = 1.0/(dt*A(i,i));

        // Setup TimeDerivative
        Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
          Teuchos::rcp(new StepperDIRKTimeDerivative<Scalar>(
            alpha,xTilde_.getConst()));

        // Setup InArgs and OutArgs
        typedef Thyra::ModelEvaluatorBase MEB;
        MEB::InArgs<Scalar>  inArgs  = this->wrapperModel_->getInArgs();
        MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
        inArgs.set_x(stageX_);
        if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(stageXDot_[i]);
        if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (ts);
        if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(dt);
        if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (alpha);
        if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (1.0);

        this->wrapperModel_->setForSolve(timeDer, inArgs, outArgs);

        if (!Teuchos::is_null(stepperDIRKObserver_))
          stepperDIRKObserver_->observeBeforeSolve(solutionHistory, *this);

        sStatus = this->solveImplicitODE(stageX_);

        if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED ) pass=false;

        if (!Teuchos::is_null(stepperDIRKObserver_))
          stepperDIRKObserver_->observeAfterSolve(solutionHistory, *this);

        timeDer->compute(stageX_, stageXDot_[i]);
      }
      if (!Teuchos::is_null(stepperDIRKObserver_))
        stepperDIRKObserver_->observeEndStage(solutionHistory, *this);
    }

    // Sum for solution: x_n = x_n-1 + Sum{ dt*b(i) * f(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int i=0; i < numStages; ++i) {
      if (b(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
        Thyra::Vp_StV((workingState->getX()).ptr(), dt*b(i), *(stageXDot_[i]));
      }
    }

    if (DIRK_ButcherTableau_->isEmbedded() ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - Implicit RK embedded methods not implemented yet!.\n");
    }

    if (pass == true)
      workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    else
      workingState->getStepperState()->stepperStatus_ = Status::FAILED;
    workingState->setOrder(this->getOrder());
    stepperObserver_->observeEndTakeStep(solutionHistory, *this);
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
StepperDIRK<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperDIRK<Scalar>::description() const
{
  return(description_);
}


template<class Scalar>
void StepperDIRK<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "wrapperModel_ = " << this->wrapperModel_->description() << std::endl;
}


template <class Scalar>
void StepperDIRK<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (this->stepperPL_ == Teuchos::null) this->stepperPL_ = this->getDefaultParameters();
  } else {
    this->stepperPL_ = pList;
  }
  // Can not validate because of optional Parameters.
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperDIRK<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  *pl = *(DIRK_ButcherTableau_->getValidParameters());
  pl->set<bool>       ("Zero Initial Guess", false);
  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperDIRK<Scalar>::getDefaultParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  *pl = *(DIRK_ButcherTableau_->getValidParameters());
  pl->set<std::string>("Solver Name", "Default Solver");
  pl->set<bool>       ("Zero Initial Guess", false);
  Teuchos::RCP<Teuchos::ParameterList> solverPL=this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperDIRK<Scalar>::getNonconstParameterList()
{
  return(this->stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperDIRK<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = this->stepperPL_;
  this->stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperDIRK_impl_hpp
