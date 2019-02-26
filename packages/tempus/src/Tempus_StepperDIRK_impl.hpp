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
  this->setTableau(stepperType);
  this->setModel(appModel);
  this->initialize();
}

template<class Scalar>
StepperDIRK<Scalar>::StepperDIRK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList>                      pList)
{
  this->setTableau(pList);
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}

template<class Scalar>
StepperDIRK<Scalar>::StepperDIRK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  std::string stepperType,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  this->setTableau(stepperType);
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}


template<class Scalar>
void StepperDIRK<Scalar>::setTableau(std::string stepperType)
{
  if (stepperType == "") {
    this->setTableau();
  } else {
    DIRK_ButcherTableau_ = createRKBT<Scalar>(stepperType, this->stepperPL_);
  }

  TEUCHOS_TEST_FOR_EXCEPTION( DIRK_ButcherTableau_->isDIRK() != true,
    std::logic_error,
       "Error - StepperDIRK did not receive a DIRK Butcher Tableau!\n"
    << "  Stepper Type = " << stepperType <<  "\n");
  description_ = DIRK_ButcherTableau_->description();
}


template<class Scalar>
void StepperDIRK<Scalar>::setTableau(Teuchos::RCP<Teuchos::ParameterList> pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (this->stepperPL_ == Teuchos::null)
      this->stepperPL_ = this->getDefaultParameters();
  } else {
    this->stepperPL_ = pList;
  }

  std::string stepperType =
    this->stepperPL_->template get<std::string>("Stepper Type",
                                                "SDIRK 2 Stage 2nd order");
  DIRK_ButcherTableau_ = createRKBT<Scalar>(stepperType, this->stepperPL_);

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
    if (this->stepperObserver_ == Teuchos::null) {
      stepperDIRKObserver_ =
        Teuchos::rcp(new StepperDIRKObserver<Scalar>());
      this->stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >
          (stepperDIRKObserver_);
     }
  } else {
    this->stepperObserver_ = obs;
    stepperDIRKObserver_ =
      Teuchos::rcp_dynamic_cast<StepperDIRKObserver<Scalar> >(this->stepperObserver_);
  }
}


template<class Scalar>
void StepperDIRK<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( this->wrapperModel_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperDIRK::initialize()\n");

  this->setTableau(this->stepperPL_);
  this->setParameterList(this->stepperPL_);
  this->setSolver();
  this->setObserver();

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

    if (DIRK_ButcherTableau_->isEmbedded() and this->getEmbedded()) {
     ee_    = Thyra::createMember(this->wrapperModel_->get_f_space());
     abs_u0 = Thyra::createMember(this->wrapperModel_->get_f_space());
     abs_u  = Thyra::createMember(this->wrapperModel_->get_f_space());
     sc     = Thyra::createMember(this->wrapperModel_->get_f_space());
  }
}


template<class Scalar>
void StepperDIRK<Scalar>::setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(numStates < 1, std::logic_error,
    "Error - setInitialConditions() needs at least one SolutionState\n"
    "        to set the initial condition.  Number of States = " << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperDIRK::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"<<std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x = initialState->getX();

  // Use x from inArgs as ICs, if needed.
  auto inArgs = this->wrapperModel_->getNominalValues();
  if (x == Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION( (x == Teuchos::null) &&
      (inArgs.get_x() == Teuchos::null), std::logic_error,
      "Error - setInitialConditions() needs the ICs from the SolutionHistory\n"
      "        or getNominalValues()!\n");

    x = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
    initialState->setX(x);
  }

  // Perform IC Consistency.  Only needed if useFSAL=true or by SolutionState.
  // Put f in the last stage of stageXDot_.  This way
  // it would be available for the first step.
  if (this->getUseFSAL() || (initialState->getXDot() != Teuchos::null)) {
    std::string icConsistency = this->getICConsistency();
    if (icConsistency == "None") {
      if (initialState->getXDot() == Teuchos::null) {
        RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"StepperDIRK::setInitialConditions()");
        *out << "Warning -- Requested IC consistency of 'None' but\n"
             << "           initialState does not have an xDot.\n"
             << "           Setting a 'Zero' xDot!\n" << std::endl;
        Thyra::assign(stageXDot_.back().ptr(), Scalar(0.0));
      } else {
        Thyra::assign(stageXDot_.back().ptr(), *(initialState->getXDot()));
      }
    }
    else if (icConsistency == "Zero")
      Thyra::assign(stageXDot_.back().ptr(), Scalar(0.0));
    else if (icConsistency == "App") {
      auto xDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
                    inArgs.get_x_dot());
      TEUCHOS_TEST_FOR_EXCEPTION(xDot == Teuchos::null, std::logic_error,
        "Error - setInitialConditions() requested 'App' for IC consistency,\n"
        "        but 'App' returned a null pointer for xDot!\n");
      Thyra::assign(stageXDot_.back().ptr(), *xDot);
    }
    else if (icConsistency == "Consistent") {
    // Solve f(x, xDot,t) = 0.
    const Scalar time = initialState->getTime();
    const Scalar dt   = initialState->getTimeStep();
    RCP<TimeDerivative<Scalar> > timeDer = Teuchos::null;
    const Scalar alpha = 1.0;    // d(xDot)/d(xDot)
    const Scalar beta  = 0.0;    // d(x   )/d(xDot)
    RCP<ImplicitODEParameters<Scalar> > p =
      Teuchos::rcp(new ImplicitODEParameters<Scalar>(timeDer,dt,alpha,beta,
                                                     SOLVE_FOR_XDOT_CONST_X));

    auto xDot = this->getStepperXDot(initialState);
    const Thyra::SolveStatus<Scalar> sStatus =
      this->solveImplicitODE(x, xDot, time, p);

    TEUCHOS_TEST_FOR_EXCEPTION(
      sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED, std::logic_error,
      "Error - Solver failed while determining the initial conditions.\n"
      "        Solver status is "<<Thyra::toString(sStatus.solveStatus)<<".\n");    }
  }

  // Put f into initialState->getXDot(), if available.
  if (initialState->getXDot() != Teuchos::null) {
    Thyra::assign(initialState->getXDot().ptr(), *(stageXDot_.back()));

    // At this point, x, and xDot sync'ed or consistent
    // at the same time level for the initialState.
    initialState->setIsSynced(true);
  }

  // Test for consistency.
  if ( this->getICConsistencyCheck() &&
       (this->getUseFSAL() || (initialState->getXDot() != Teuchos::null)) ) {
    auto xDot = stageXDot_.back();
    auto f    = initialState->getX()->clone_v();

    const Scalar time = initialState->getTime();
    const Scalar dt   = initialState->getTimeStep();
    RCP<TimeDerivative<Scalar> > timeDer = Teuchos::null;
    const Scalar alpha = 0.0;
    const Scalar beta  = 0.0;
    RCP<ImplicitODEParameters<Scalar> > p =
      Teuchos::rcp(new ImplicitODEParameters<Scalar>(timeDer,dt,alpha,beta,
                                                     EVALUATE_RESIDUAL));

    this->evaluateImplicitODE(f, x, xDot, time, p);

    Scalar reldiff = Thyra::norm(*f)/Thyra::norm(*x);
    Scalar eps = Scalar(100.0)*std::abs(Teuchos::ScalarTraits<Scalar>::eps());
    if (reldiff > eps) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"StepperBackwardEuler::setInitialConditions()");
      *out << "Warning -- Failed consistency check but continuing!\n"
         << "  ||f(x,xDot,t)||/||x|| > eps" << std::endl
         << "  ||f(x,xDot,t)||       = " << Thyra::norm(*f) << std::endl
         << "  ||x||                 = " << Thyra::norm(*x) << std::endl
         << "  ||f(x,xDot,t)||/||x|| = " << reldiff         << std::endl
         << "                    eps = " << eps             << std::endl;
    }
  }
}


template<class Scalar>
void StepperDIRK<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperDIRK::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperDIRK<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for DIRK.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
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

      if ( i == 0 && this->getUseFSAL() &&
           workingState->getNConsecutiveFailures() == 0 ) {

        RCP<Thyra::VectorBase<Scalar> > tmp = stageXDot_[0];
        stageXDot_[0] = stageXDot_.back();
        stageXDot_.back() = tmp;

      } else {

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
          const Scalar alpha = 1.0/(dt*A(i,i));
          const Scalar beta  = 1.0;

          // Setup TimeDerivative
          Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
            Teuchos::rcp(new StepperDIRKTimeDerivative<Scalar>(
              alpha,xTilde_.getConst()));

          Teuchos::RCP<ImplicitODEParameters<Scalar> > p =
            Teuchos::rcp(new ImplicitODEParameters<Scalar>(
              timeDer, dt, alpha, beta));
          p->stageNumber_ = i;

          if (!Teuchos::is_null(stepperDIRKObserver_))
            stepperDIRKObserver_->observeBeforeSolve(solutionHistory, *this);

          sStatus = this->solveImplicitODE(stageX_, stageXDot_[i], ts, p);

          if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) pass=false;

          if (!Teuchos::is_null(stepperDIRKObserver_))
            stepperDIRKObserver_->observeAfterSolve(solutionHistory, *this);

          timeDer->compute(stageX_, stageXDot_[i]);
        }
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

    if (DIRK_ButcherTableau_->isEmbedded() and this->getEmbedded()) {
      RCP<SolutionStateMetaData<Scalar> > metaData=workingState->getMetaData();
      const Scalar tolAbs = metaData->getTolRel();
      const Scalar tolRel = metaData->getTolAbs();

      // just compute the error weight vector
      // (all that is needed is the error, and not the embedded solution)
      Teuchos::SerialDenseVector<int,Scalar> errWght = b ;
      errWght -= DIRK_ButcherTableau_->bstar();

      // compute local truncation error estimate: | u^{n+1} - \hat{u}^{n+1} |
      // Sum for solution: ee_n = Sum{ (b(i) - bstar(i)) * dt*f(i) }
      assign(ee_.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
      for (int i=0; i < numStages; ++i) {
         if (errWght(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
            Thyra::Vp_StV(ee_.ptr(), dt*errWght(i), *(stageXDot_[i]));
         }
      }

      // compute: Atol + max(|u^n|, |u^{n+1}| ) * Rtol
      Thyra::abs( *(currentState->getX()), abs_u0.ptr());
      Thyra::abs( *(workingState->getX()), abs_u.ptr());
      Thyra::pair_wise_max_update(tolRel, *abs_u0, abs_u.ptr());
      Thyra::add_scalar(tolAbs, abs_u.ptr());

      // compute: || ee / sc ||
      assign(sc.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
      Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *ee_, *abs_u, sc.ptr());
      Scalar err = std::abs(Thyra::norm_inf(*sc));
      metaData->setErrorRel(err);

      // test if step should be rejected
      if (std::isinf(err) || std::isnan(err) || err > Teuchos::as<Scalar>(1.0))
        pass = false;
    }

    if (pass) workingState->setSolutionStatus(Status::PASSED);
    else      workingState->setSolutionStatus(Status::FAILED);

    workingState->setOrder(this->getOrder());
    this->stepperObserver_->observeEndTakeStep(solutionHistory, *this);
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
  this->getValidParametersBasic(pl);
  pl->set<bool>("Initial Condition Consistency Check", false);
  pl->set<bool>("Zero Initial Guess", false);
  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperDIRK<Scalar>::getDefaultParameters() const
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::rcp_const_cast;

  RCP<ParameterList> pl =
    rcp_const_cast<ParameterList>(this->getValidParameters());

  if (DIRK_ButcherTableau_ == Teuchos::null) {
    auto DIRK_ButcherTableau =
      createRKBT<Scalar>("SDIRK 2 Stage 2nd order", Teuchos::null);
    pl->setParameters(*(DIRK_ButcherTableau->getValidParameters()));
  } else {
    pl->setParameters(*(DIRK_ButcherTableau_->getValidParameters()));
  }

  pl->set<std::string>("Solver Name", "Default Solver");
  RCP<ParameterList> solverPL = this->defaultSolverParameters();
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
