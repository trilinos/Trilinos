#ifndef Tempus_StepperExplicitRK_impl_hpp
#define Tempus_StepperExplicitRK_impl_hpp

#include "Tempus_RKButcherTableauBuilder.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Tempus {

// StepperExplicitRK definitions:
template<class Scalar>
StepperExplicitRK<Scalar>::StepperExplicitRK(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel )
{
  this->setParameterList(pList);
  std::string stepperType = pList_->get<std::string>("Stepper Type");

  ERK_ButcherTableau_ = createRKBT<Scalar>(stepperType);
  description_ = ERK_ButcherTableau_->description();
  TEUCHOS_TEST_FOR_EXCEPTION( ERK_ButcherTableau_->isImplicit() == true,
    std::logic_error,
       "Error - StepperExplicitRK received an implicit Butcher Tableau!\n"
    << "  Stepper Type = "<< pList_->get<std::string>("Stepper Type") << "\n");

  this->validExplicitODE(transientModel);
  eODEModel_ = transientModel;

  inArgs_  = eODEModel_->createInArgs();
  outArgs_ = eODEModel_->createOutArgs();
  inArgs_  = eODEModel_->getNominalValues();

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
  Teuchos::RCP<SolutionState<Scalar> > currentState =
    solutionHistory->getCurrentState();
  Teuchos::RCP<SolutionState<Scalar> > workingState =
    solutionHistory->getWorkingState();
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

  workingState->stepperState_->stepperStatus_ = Status::PASSED;
  workingState->setOrder(ERK_ButcherTableau_->order());
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
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pList));
  //pList->validateParameters(*this->getValidParameters());
  pList_ = pList;

  std::string stepperType = pList_->get<std::string>("Stepper Type");

  Teuchos::readVerboseObjectSublist(&*pList_,this);
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperExplicitRK<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    // pl->set("Stepper Type", "RK Forward Euler",
    //   "'Stepper Type' sets the stepper method.\n"
    //   "For Explicit RK the following methods are valid:\n"
    //   "  RK Forward Euler\n",
    //   "  RK Explicit 4 Stage\n",
    //   "  RK Explicit 3/8 Rule\n",
    //   "  RK Explicit 4 Stage 3rd order by Runge\n",
    //   "  RK Explicit 5 Stage 3rd order by Kinnmark and Gray\n",
    //   "  RK Explicit 3 Stage 3rd order\n",
    //   "  RK Explicit 3 Stage 3rd order TVD\n",
    //   "  RK Explicit 3 Stage 3rd order by Heun\n",
    //   "  RK Explicit 2 Stage 2nd order by Runge\n",
    //   "  RK Explicit Trapezoidal\n",
    //   StepperERKTypeValidator);

    validPL = pl;
  }
  return validPL;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperExplicitRK<Scalar>::getNonconstParameterList()
{
  return(pList_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperExplicitRK<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = pList_;
  pList_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperExplicitRK_impl_hpp
