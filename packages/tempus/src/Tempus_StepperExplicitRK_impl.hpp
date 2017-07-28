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
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  std::string stepperType,
  Teuchos::RCP<Teuchos::ParameterList>                      pList)
{
  this->setTableau(pList, stepperType);
  this->setParameterList(pList);
  this->setModel(appModel);
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
  // Initialize the stage vectors
  int numStages = ERK_ButcherTableau_->numStages();
  stageX_ = Thyra::createMember(appModel_->get_f_space());
  stageXDot_.resize(numStages);
  for (int i=0; i<numStages; ++i) {
    stageXDot_[i] = Thyra::createMember(appModel_->get_f_space());
    assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
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

    const int numStages = ERK_ButcherTableau_->numStages();
    Teuchos::SerialDenseMatrix<int,Scalar> A = ERK_ButcherTableau_->A();
    Teuchos::SerialDenseVector<int,Scalar> b = ERK_ButcherTableau_->b();
    Teuchos::SerialDenseVector<int,Scalar> c = ERK_ButcherTableau_->c();

    // Compute stage solutions
    for (int i=0; i < numStages; ++i) {
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

      appModel_->evalModel(inArgs_,outArgs_);
      // --------------------------------
    }

    // Sum for solution: x_n = x_n-1 + Sum{ b(i) * dt*f(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int i=0; i < numStages; ++i) {
      if (b(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
        Thyra::Vp_StV((workingState->getX()).ptr(), dt*b(i), *(stageXDot_[i]));
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
