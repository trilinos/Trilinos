// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRK_impl_hpp
#define Tempus_StepperExplicitRK_impl_hpp

#include "Tempus_RKButcherTableauBuilder.hpp"
#include "Tempus_RKButcherTableau.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Tempus {

template<class Scalar>
StepperExplicitRK<Scalar>::StepperExplicitRK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  std::string stepperType)
{
  this->setTableau(Teuchos::null, stepperType);
  this->setParameterList(Teuchos::null);
  this->setModel(appModel);
  this->setObserver();
  this->initialize();
}

template<class Scalar>
StepperExplicitRK<Scalar>::StepperExplicitRK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList>                      pList)
{
  this->setTableau(pList, "RK Explicit 4 Stage");
  this->setParameterList(pList);
  this->setModel(appModel);
  this->setObserver();
  this->initialize();
}

template<class Scalar>
StepperExplicitRK<Scalar>::StepperExplicitRK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  std::string stepperType,
  Teuchos::RCP<Teuchos::ParameterList>                      pList)
{
  this->setTableau(pList, stepperType);
  this->setParameterList(pList);
  this->setModel(appModel);
  this->setObserver();
  this->initialize();
}

template<class Scalar>
void StepperExplicitRK<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->validExplicitODE(appModel);
  appModel_ = appModel;

  inArgs_  = appModel_->getNominalValues();
  outArgs_ = appModel_->createOutArgs();
}

template<class Scalar>
void StepperExplicitRK<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}

template<class Scalar>
void StepperExplicitRK<Scalar>::setSolver(std::string solverName)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperExplicitRK::setSolver()");
  *out << "Warning -- No solver to set for StepperExplicitRK "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperExplicitRK<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperExplicitRK::setSolver()");
  *out << "Warning -- No solver to set for StepperExplicitRK "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperExplicitRK<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperExplicitRK::setSolver()");
  *out << "Warning -- No solver to set for StepperExplicitRK "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperExplicitRK<Scalar>::setTableau(
  Teuchos::RCP<Teuchos::ParameterList> pList,
  std::string stepperType)
{
  if (stepperType == "") {
    if (pList == Teuchos::null)
      stepperType = "RK Explicit 4 Stage";
    else
      stepperType =
        pList->get<std::string>("Stepper Type", "RK Explicit 4 Stage");
  }

  ERK_ButcherTableau_ = createRKBT<Scalar>(stepperType, pList);

  TEUCHOS_TEST_FOR_EXCEPTION(ERK_ButcherTableau_->isImplicit() == true,
    std::logic_error,
       "Error - StepperExplicitRK received an implicit Butcher Tableau!\n"
    << "  Stepper Type = " << stepperType << "\n");
  description_ = ERK_ButcherTableau_->description();
}


template<class Scalar>
void StepperExplicitRK<Scalar>::setObserver(
  Teuchos::RCP<StepperExplicitRKObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (stepperExplicitRKObserver_ == Teuchos::null) {
      stepperExplicitRKObserver_ =
        Teuchos::rcp(new StepperExplicitRKObserver<Scalar>());
    }
  } else {
    stepperExplicitRKObserver_ = obs;
  }
}


template<class Scalar>
void StepperExplicitRK<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( ERK_ButcherTableau_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the Butcher Tableau, setTableau(), before calling "
    "StepperExplicitRK::initialize()\n");

  TEUCHOS_TEST_FOR_EXCEPTION( appModel_ == Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperExplicitRK::initialize()\n");

  // Initialize the stage vectors
  int numStages = ERK_ButcherTableau_->numStages();
  stageX_ = Thyra::createMember(appModel_->get_f_space());
  stageXDot_.resize(numStages);
  for (int i=0; i<numStages; ++i) {
    stageXDot_[i] = Thyra::createMember(appModel_->get_f_space());
    assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  }

  if (ERK_ButcherTableau_->isEmbedded() and stepperPL_->get<bool>("Use Embedded")){
     ee_ = Thyra::createMember(appModel_->get_f_space());
     abs_u0 = Thyra::createMember(appModel_->get_f_space());
     abs_u = Thyra::createMember(appModel_->get_f_space());
     sc = Thyra::createMember(appModel_->get_f_space());
  }
}

template<class Scalar>
void StepperExplicitRK<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperExplicitRK::takeStep()");
  {
    stepperExplicitRKObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    const int numStages = ERK_ButcherTableau_->numStages();
    Teuchos::SerialDenseMatrix<int,Scalar> A = ERK_ButcherTableau_->A();
    Teuchos::SerialDenseVector<int,Scalar> b = ERK_ButcherTableau_->b();
    Teuchos::SerialDenseVector<int,Scalar> c = ERK_ButcherTableau_->c();

    // Compute stage solutions
    for (int i=0; i < numStages; ++i) {
      stepperExplicitRKObserver_->observeBeginStage(solutionHistory, *this);
      Thyra::assign(stageX_.ptr(), *(currentState->getX()));
      for (int j=0; j < i; ++j) {
        if (A(i,j) != Teuchos::ScalarTraits<Scalar>::zero()) {
          Thyra::Vp_StV(stageX_.ptr(), dt*A(i,j), *stageXDot_[j]);
        }
      }
      const Scalar ts = time + c(i)*dt;

      // Evaluate model -----------------
      typedef Thyra::ModelEvaluatorBase MEB;
      inArgs_.set_x(stageX_);
      if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(ts);

      if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);
      outArgs_.set_f(stageXDot_[i]);

      stepperExplicitRKObserver_->observeBeforeExplicit(solutionHistory, *this);
      appModel_->evalModel(inArgs_,outArgs_);
      // --------------------------------
      stepperExplicitRKObserver_->observeEndStage(solutionHistory, *this);
    }

    // Sum for solution: x_n = x_n-1 + Sum{ b(i) * dt*f(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int i=0; i < numStages; ++i) {
      if (b(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
        Thyra::Vp_StV((workingState->getX()).ptr(), dt*b(i), *(stageXDot_[i]));
      }
    }


    // At this point, the stepper has passed.
    // but when using adaptive time stepping, the embedded method can change the step status
    workingState->getStepperState()->stepperStatus_ = Status::PASSED;

    if (ERK_ButcherTableau_->isEmbedded() and stepperPL_->get<bool>("Use Embedded")){

       RCP<SolutionStateMetaData<Scalar> > metaData = workingState->getMetaData();
       const Scalar tolAbs = metaData->getTolRel();
       const Scalar tolRel = metaData->getTolAbs();

       // just compute the error weight vector
       // (all that is needed is the error, and not the embedded solution)
       Teuchos::SerialDenseVector<int,Scalar> errWght = b ;
       errWght -= ERK_ButcherTableau_->bstar();

       //compute local truncation error estimate: | u^{n+1} - \hat{u}^{n+1} |
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

       //compute: || ee / sc ||
       assign(sc.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
       Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *ee_, *abs_u, sc.ptr());
       Scalar err = Thyra::norm_inf(*sc);
       metaData->setErrorRel(err);

       // test if step should be rejected
       if (err > 1.0){
          workingState->getStepperState()->stepperStatus_ = Status::FAILED;
       }
    }

    workingState->setOrder(this->getOrder());
    stepperExplicitRKObserver_->observeEndTakeStep(solutionHistory, *this);
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
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperExplicitRK<Scalar>::description() const
{
  return(description_);
}


template<class Scalar>
void StepperExplicitRK<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "appModel_ = " << appModel_->description() << std::endl;
}


template <class Scalar>
void StepperExplicitRK<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (stepperPL_ == Teuchos::null) stepperPL_ = this->getDefaultParameters();
  } else {
    stepperPL_ = pList;
  }
  // Can not validate because of optional Parameters.
  stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters(),0);
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperExplicitRK<Scalar>::getValidParameters() const
{
  //std::stringstream Description;
  //Description << "'Stepper Type' sets the stepper method.\n"
  //            << "For Explicit RK the following methods are valid:\n"
  //            << "  General ERK\n"
  //            << "  RK Forward Euler\n"
  //            << "  RK Explicit 4 Stage\n"
  //            << "  RK Explicit 3/8 Rule\n"
  //            << "  RK Explicit 4 Stage 3rd order by Runge\n"
  //            << "  RK Explicit 5 Stage 3rd order by Kinnmark and Gray\n"
  //            << "  RK Explicit 3 Stage 3rd order\n"
  //            << "  RK Explicit 3 Stage 3rd order TVD\n"
  //            << "  RK Explicit 3 Stage 3rd order by Heun\n"
  //            << "  RK Explicit 2 Stage 2nd order by Runge\n"
  //            << "  RK Explicit Trapezoidal\n";

  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->set<bool>("Use Embedded", false,
    "'Whether to use Embedded Stepper (if available) or not\n"
    "  'true' - Stepper will compute embedded solution and is adaptive.\n"
    "  'false' - Stepper is not embedded(adaptive).\n");
  pl->setParameters(*(ERK_ButcherTableau_->getValidParameters()));
  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperExplicitRK<Scalar>::getDefaultParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->set<bool>("Use Embedded", false,
    "'Whether to use Embedded Stepper (if available) or not\n"
    "  'true' - Stepper will compute embedded solution and is adaptive.\n"
    "  'false' - Stepper is not embedded(adaptive).\n");
  pl->setParameters(*(ERK_ButcherTableau_->getValidParameters()));
  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperExplicitRK<Scalar>::getNonconstParameterList()
{
  return(stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperExplicitRK<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperExplicitRK_impl_hpp
