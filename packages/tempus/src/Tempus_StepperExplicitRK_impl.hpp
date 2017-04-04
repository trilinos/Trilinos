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

// StepperExplicitRK definitions:
template<class Scalar>
StepperExplicitRK<Scalar>::StepperExplicitRK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel,
  std::string stepperType,
  Teuchos::RCP<Teuchos::ParameterList>                      pList)
{
  this->setTableau(pList, stepperType);
  this->setParameterList(pList);
  this->setModel(transientModel);
  this->initialize();
}

template<class Scalar>
void StepperExplicitRK<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel)
{
  this->validExplicitODE(transientModel);
  eODEModel_ = transientModel;

  inArgs_  = eODEModel_->createInArgs();
  outArgs_ = eODEModel_->createOutArgs();
  inArgs_  = eODEModel_->getNominalValues();
}

template<class Scalar>
void StepperExplicitRK<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel)
{
  this->setModel(transientModel);
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
      stepperType = pList->get<std::string>("Stepper Type");
  }

  ERK_ButcherTableau_ = createRKBT<Scalar>(stepperType,pList);

  TEUCHOS_TEST_FOR_EXCEPTION(ERK_ButcherTableau_->isImplicit() == true,
    std::logic_error,
       "Error - StepperExplicitRK received an implicit Butcher Tableau!\n"
    << "  Stepper Type = " << stepperType << "\n");
  description_ = ERK_ButcherTableau_->description();
}

template<class Scalar>
void StepperExplicitRK<Scalar>::initialize()
{
  stageX_ = Thyra::createMember(eODEModel_->get_f_space());
  // Initialize the stage vectors
  int numStages = ERK_ButcherTableau_->numStages();
  stagef_.reserve(numStages);
  for (int i=0; i<numStages; ++i) {
    stagef_.push_back(Thyra::createMember(eODEModel_->get_f_space()));
  }
}

template<class Scalar>
void StepperExplicitRK<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperExplicitRK::takeStep()");
  {
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    int stages = ERK_ButcherTableau_->numStages();
    Teuchos::SerialDenseMatrix<int,Scalar> A = ERK_ButcherTableau_->A();
    Teuchos::SerialDenseVector<int,Scalar> b = ERK_ButcherTableau_->b();
    Teuchos::SerialDenseVector<int,Scalar> c = ERK_ButcherTableau_->c();

    // Compute stage solutions
    for (int s=0 ; s < stages ; ++s) {
      Thyra::assign(stageX_.ptr(), *(currentState->getX()));
      for (int j=0 ; j < s ; ++j) {
        if (A(s,j) != Teuchos::ScalarTraits<Scalar>::zero()) {
          Thyra::Vp_StV(stageX_.ptr(), A(s,j), *stagef_[j]);
        }
      }
      typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag TScalarMag;
      TScalarMag ts = time + c(s)*dt;

      // Evaluate model -----------------
      //explicitEvalModel(currentState);
      typedef Thyra::ModelEvaluatorBase MEB;
      inArgs_.set_x(stageX_);
      if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(ts);

      if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);
      outArgs_.set_f(stagef_[s]);

      eODEModel_->evalModel(inArgs_,outArgs_);
      // --------------------------------

      Thyra::Vt_S(stagef_[s].ptr(),dt);
    }

    // Sum for solution: x_n = x_n-1 + Sum{ b(s) * dt*f(s) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int s=0 ; s < stages ; ++s) {
      if (b(s) != Teuchos::ScalarTraits<Scalar>::zero()) {
        Thyra::Vp_StV((workingState->getX()).ptr(), b(s), *(stagef_[s]));
      }
    }

    if (ERK_ButcherTableau_->isEmbedded() ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - Explicit RK embedded methods not implemented yet!.\n");
    }

    workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    workingState->setOrder(this->getOrder());
  }
  return;
}


template<class Scalar>
void StepperExplicitRK<Scalar>::
explicitEvalModel(Teuchos::RCP<SolutionState<Scalar> > state)
{
  // NOTE: on input state->getX() has the current solution x, and
  //       on output state->getXDot() has the current f(x,t) [=xDot].
  typedef Thyra::ModelEvaluatorBase MEB;
  inArgs_.set_x(state->getX());
  if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(state->getTime());

  // For model evaluators whose state function f(x, x_dot, t) describes
  // an implicit ODE, and which accept an optional x_dot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // x_dot = f(x, t)
  if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);
  outArgs_.set_f(state->getXDot());

  eODEModel_->evalModel(inArgs_,outArgs_);
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
      << "eODEModel_ = " << eODEModel_->description() << std::endl;
}


template <class Scalar>
void StepperExplicitRK<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  if (pList == Teuchos::null) {
    stepperPL_ = this->getDefaultParameters();
  } else {
    stepperPL_ = pList;
  }
  // Can not validate because of optional Parameters.
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
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

  return ERK_ButcherTableau_->getValidParameters();
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperExplicitRK<Scalar>::getDefaultParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  *pl = *(ERK_ButcherTableau_->getValidParameters());
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
