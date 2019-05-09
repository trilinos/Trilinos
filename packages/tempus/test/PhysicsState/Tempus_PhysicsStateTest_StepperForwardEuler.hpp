// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_PhysicsStateTest_StepperForwardEuler_hpp
#define Tempus_PhysicsStateTest_StepperForwardEuler_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_PhysicsStateCounter.hpp"


namespace Tempus_Test {

/** \brief This is a Forward Euler time stepper to test the PhysicsState.
 *
 *  It is derived from StepperForwardEuler, and simply increments
 *  a physics counter.
 */
template<class Scalar>
class PhysicsStateTest_StepperForwardEuler
  : virtual public Tempus::StepperForwardEuler<Scalar>
{
public:

  /// Constructor
  PhysicsStateTest_StepperForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null)
    : Tempus::StepperForwardEuler<Scalar>(appModel, pList) {}

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
    const Teuchos::RCP<Tempus::SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperForwardEuler::takeStep()");
  {
    this->stepperFEObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<Tempus::SolutionState<Scalar> > currentState =
      solutionHistory->getCurrentState();

    typedef Thyra::ModelEvaluatorBase MEB;
    this->inArgs_.set_x(currentState->getX());
    if (this->inArgs_.supports(MEB::IN_ARG_t))
      this->inArgs_.set_t(currentState->getTime());

    // For model evaluators whose state function f(x, x_dot, t) describes
    // an implicit ODE, and which accept an optional x_dot input argument,
    // make sure the latter is set to null in order to request the evaluation
    // of a state function corresponding to the explicit ODE formulation
    // x_dot = f(x, t)
    if (this->inArgs_.supports(MEB::IN_ARG_x_dot))
      this->inArgs_.set_x_dot(Teuchos::null);
    this->outArgs_.set_f(currentState->getXDot());

    this->stepperFEObserver_->observeBeforeExplicit(solutionHistory, *this);

    this->appModel_->evalModel(this->inArgs_,this->outArgs_);

    // Forward Euler update, x = x + dt*xdot
    RCP<Tempus::SolutionState<Scalar> > workingState =
      solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
      *(currentState->getX()),dt,*(currentState->getXDot()));

    RCP<PhysicsStateCounter<Scalar> > pSC =
      Teuchos::rcp_dynamic_cast<PhysicsStateCounter<Scalar> >
        (workingState->getPhysicsState());
    int counter = pSC->getCounter();
    counter++;
    pSC->setCounter(counter);

    workingState->setSolutionStatus(Tempus::Status::PASSED);
    workingState->setOrder(this->getOrder());
    this->stepperFEObserver_->observeEndTakeStep(solutionHistory, *this);
  }
  return;
}

};

} // namespace Tempus_Test

#endif // Tempus_PhysicsStateTest_StepperForwardEuler_hpp
