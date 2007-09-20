//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_SIMPLE_INTEGRATOR_H
#define Rythmos_SIMPLE_INTEGRATOR_H


#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_InterpolationBufferHelpers.hpp"
#include "Rythmos_IntegrationControlStrategyBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_StepperHelpers.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief A very simple concrete subclass for <tt>IntegratorBase</tt> that
 * allows just for simple fixed steps or variable steps.
 */
template<class Scalar> 
class SimpleIntegrator
  : virtual public IntegratorBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:
  
  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors, Initializers, Misc */
  //@{
  
  /** \brief . */
  SimpleIntegrator();

  /** \brief . */
  void setIntegrationControlStrategy(
    const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy
    );

  /** \brief . */
  void setIntegrationObserver(
    const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
    );

  //@}

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from IntegratorBase */
  //@{

  /** \brief . */
  RCP<IntegratorBase<Scalar> > cloneIntegrator() const;
  
  /** \brief . */
  void setStepper(
    const RCP<StepperBase<Scalar> > &stepper,
    const Scalar &finalTime,
    const bool landOnFinalTime
    );

  /** \brief . */
  RCP<StepperBase<Scalar> > unSetStepper();

  /** \brief . */
  Teuchos::RCP<const StepperBase<Scalar> > getStepper() const;

  /** \brief . */
  void getFwdPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    );

  /** \brief . */
  TimeRange<Scalar> getFwdTimeRange() const;

  //@}

  /** \name Overridden from InterpolationBufferBase */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
    
  /** \brief . */
  void addPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    );

  /** \brief . */
  void getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const;

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;

  /** \brief . */
  void getNodes(Array<Scalar>* time_vec) const;

  /** \brief . */
  void removeNodes(Array<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}

private:

  // ////////////////////////
  // Private data members

  RCP<IntegrationControlStrategyBase<Scalar> > integrationControlStrategy_;
  RCP<IntegrationObserverBase<Scalar> > integrationObserver_;
  RCP<StepperBase<Scalar> > stepper_;
  TimeRange<Scalar> integrationTimeDomain_;
  bool landOnFinalTime_;

  int maxNumTimeSteps_;

  int currTimeStepIndex_;
  StepControlInfo<Scalar> stepCtrlInfoLast_;

  static const std::string maxNumTimeSteps_name_;
  static const int maxNumTimeSteps_default_;

  // /////////////////////////
  // Private member functions

  bool advanceStepperToTime( const Scalar& t );

};


/** \brief .
 *
 * \relates SimpleIntegrator
 */
template<class Scalar> 
RCP<SimpleIntegrator<Scalar> >
simpleIntegrator()
{
  RCP<SimpleIntegrator<Scalar> >
    integrator = Teuchos::rcp(new SimpleIntegrator<Scalar>());
  return integrator;
}


/** \brief .
 *
 * \relates SimpleIntegrator
 */
template<class Scalar> 
RCP<SimpleIntegrator<Scalar> >
controlledSimpleIntegrator(
  const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy
  )
{
  RCP<SimpleIntegrator<Scalar> >
    integrator = Teuchos::rcp(new SimpleIntegrator<Scalar>());
  integrator->setIntegrationControlStrategy(integrationControlStrategy);
  return integrator;
}


/** \brief .
 *
 * \relates SimpleIntegrator
 */
template<class Scalar> 
RCP<SimpleIntegrator<Scalar> >
observedSimpleIntegrator(
  const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
  )
{
  RCP<SimpleIntegrator<Scalar> >
    integrator = Teuchos::rcp(new SimpleIntegrator<Scalar>());
  integrator->setIntegrationObserver(integrationObserver);
  return integrator;
}


/** \brief .
 *
 * \relates SimpleIntegrator
 */
template<class Scalar> 
RCP<SimpleIntegrator<Scalar> >
simpleIntegrator(
  const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy,
  const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
  )
{
  RCP<SimpleIntegrator<Scalar> >
    integrator = Teuchos::rcp(new SimpleIntegrator<Scalar>());
  integrator->setIntegrationControlStrategy(integrationControlStrategy);
  integrator->setIntegrationObserver(integrationObserver);
  return integrator;
}


// 2007/08/30: rabartl: Above, note that I had to name the nonmember
// constructors taking an single RCP argument different names from each other
// in order to get around the classic ambiguity problem with implicit
// conversions of smart pointers.


// ////////////////////////////
// Defintions


// Static data members


template<class Scalar>
const std::string
SimpleIntegrator<Scalar>::maxNumTimeSteps_name_ = "Max Number Time Steps";

template<class Scalar>
const int
SimpleIntegrator<Scalar>::maxNumTimeSteps_default_ = 10000;


// Constructors, Initializers, Misc


template<class Scalar>
SimpleIntegrator<Scalar>::SimpleIntegrator()
  :landOnFinalTime_(true),
   maxNumTimeSteps_(maxNumTimeSteps_default_),
   currTimeStepIndex_(-1)
{}


template<class Scalar>
void SimpleIntegrator<Scalar>::setIntegrationControlStrategy(
  const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(integrationControlStrategy));
#endif
  integrationControlStrategy_ = integrationControlStrategy;
}


template<class Scalar>
void SimpleIntegrator<Scalar>::setIntegrationObserver(
  const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(integrationObserver));
#endif
  integrationObserver_ = integrationObserver;
}


// Overridden from ParameterListAcceptor


template<class Scalar> 
void SimpleIntegrator<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters());
  this->setMyParamList(paramList);
  maxNumTimeSteps_ = paramList->get(
    maxNumTimeSteps_name_, maxNumTimeSteps_default_);
  Teuchos::readVerboseObjectSublist(&*paramList,this);
}


template<class Scalar> 
RCP<const ParameterList>
SimpleIntegrator<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL) ) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(maxNumTimeSteps_name_, maxNumTimeSteps_default_,
      "Set the maximum number of integration time-steps allowed."
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from IntegratorBase


template<class Scalar>
RCP<IntegratorBase<Scalar> >
SimpleIntegrator<Scalar>::cloneIntegrator() const
{
  RCP<SimpleIntegrator<Scalar> >
    newIntegrator = Teuchos::rcp(new SimpleIntegrator<Scalar>());
  // Only copy control information, not the state of an existing integration!
  newIntegrator->stepper_ = Teuchos::null;
  const RCP<const ParameterList> paramList = this->getParameterList();
  if (!is_null(paramList))
    newIntegrator->setParameterList(Teuchos::parameterList(*paramList));
  if(!is_null(integrationControlStrategy_))
    newIntegrator->integrationControlStrategy_ =
      integrationControlStrategy_->cloneIntegrationControlStrategy().assert_not_null();
  if(!is_null(integrationObserver_))
    newIntegrator->integrationObserver_ =
      integrationObserver_->cloneIntegrationObserver().assert_not_null();
  return newIntegrator;
}


template<class Scalar> 
void SimpleIntegrator<Scalar>::setStepper(
  const RCP<StepperBase<Scalar> > &stepper,
  const Scalar &finalTime,
  const bool landOnFinalTime
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT(is_null(stepper));
  TEST_FOR_EXCEPT( finalTime <= stepper->getTimeRange().lower() );
  TEUCHOS_ASSERT( stepper->getTimeRange().length() == ST::zero() );
  // 2007/07/25: rabartl: ToDo: Validate state of the stepper!
  stepper_ = stepper;
  integrationTimeDomain_ = timeRange(stepper_->getTimeRange().lower(), finalTime);
  landOnFinalTime_ = landOnFinalTime;
  currTimeStepIndex_ = 0;
  stepCtrlInfoLast_ = StepControlInfo<Scalar>();
  if (!is_null(integrationControlStrategy_))
    integrationControlStrategy_->resetIntegrationControlStrategy(
      integrationTimeDomain_
      );
  if (!is_null(integrationObserver_))
    integrationObserver_->resetIntegrationObserver(
      integrationTimeDomain_
      );
}

template<class Scalar>
RCP<StepperBase<Scalar> > SimpleIntegrator<Scalar>::unSetStepper()
{
  RCP<StepperBase<Scalar> > stepper_temp = stepper_;
  stepper_ = Teuchos::null;
  return(stepper_temp);
}

template<class Scalar>
Teuchos::RCP<const StepperBase<Scalar> > SimpleIntegrator<Scalar>::getStepper() const
{
  return(stepper_);
}

template<class Scalar>
void SimpleIntegrator<Scalar>::getFwdPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  )
{

#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Rythmos:SimpleIntegrator::getFwdPoints");
#endif

  using Teuchos::incrVerbLevel;
  using Teuchos::Describable;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef InterpolationBufferBase<Scalar> IBB;
  typedef Teuchos::VerboseObjectTempState<IBB> VOTSIBB;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  VOTSIBB stepper_outputTempState(stepper_,out,incrVerbLevel(verbLevel,-1));

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nEntering " << this->Describable::description() << "::getFwdPoints(...) ...\n"
         << "\nStepper: " << Teuchos::describe(*stepper_,verbLevel);

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) )
    *out << "\nRequested time points: " << Teuchos::toString(time_vec) << "\n";

  //
  // 1) Initial setup
  //

  const int numTimePoints = time_vec.size();

  // Assert preconditions
  assertTimePointsAreSorted(time_vec);
  TEST_FOR_EXCEPT(accuracy_vec!=0); // ToDo: Remove accuracy_vec!

  // Resize the storage for the output arrays
  if (x_vec)
    x_vec->resize(numTimePoints);
  if (xdot_vec)
    xdot_vec->resize(numTimePoints);

  // This int records the next time point offset in time_vec[timePointIndex]
  // that needs to be handled.  This gets updated as the time points are
  // filled below.
  int nextTimePointIndex = 0;
  
  assertNoTimePointsBeforeCurrentTimeRange(*this,time_vec,nextTimePointIndex);

  //
  // 1) First, get all time points that fall within the current time range
  //

  {
#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("Rythmos:SimpleIntegrator::getFwdPoints: getPoints");
#endif
    getCurrentPoints(*stepper_,time_vec,x_vec,xdot_vec,&nextTimePointIndex);
  }

  //
  // 2) Advance the stepper to satisfy time points in time_vec that fall
  // before the current time.
  //

  while ( nextTimePointIndex < numTimePoints ) {
    
    // Use the time stepping algorithm to step up to or past the next
    // requested time but not so far as to step past the point entirely.
    const Scalar t = time_vec[nextTimePointIndex];
    bool advanceStepperToTimeSucceeded = false;
    {
#ifdef ENABLE_RYTHMOS_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:SimpleIntegrator::getFwdPoints: advanceStepperToTime");
#endif
      advanceStepperToTimeSucceeded= advanceStepperToTime(t);
    }
    TEST_FOR_EXCEPTION(
      !advanceStepperToTimeSucceeded, Exceptions::GetFwdPointsFailed,
      this->description() << "\n\n"
      "Error:  The integration failed to get to time " << t << " and only achieved\n"
      "getting to " << stepper_->getTimeRange().upper() << "!"
      );
    
    // Extract the next set of points (perhaps just one) from the stepper
    {
#ifdef ENABLE_RYTHMOS_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:SimpleIntegrator::getFwdPoints: getPoints (fwd)");
#endif
      getCurrentPoints(*stepper_,time_vec,x_vec,xdot_vec,&nextTimePointIndex);
    }
    
  }

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nLeaving " << this->Describable::description() << "::getFwdPoints(...) ...\n";
  
}


template<class Scalar> 
TimeRange<Scalar>
SimpleIntegrator<Scalar>::getFwdTimeRange() const
{
  return timeRange(
    stepper_->getTimeRange().lower(),
    integrationTimeDomain_.upper()
    );
}


// Overridden from InterpolationBufferBase

template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> > SimpleIntegrator<Scalar>::get_x_space() const
{
  if (stepper_ == Teuchos::null) {
    RCP<const Thyra::VectorSpaceBase<Scalar> > space;
    return(space);
  } else {
    return(stepper_->get_x_space());
  }
}

template<class Scalar> 
void SimpleIntegrator<Scalar>::addPoints(
  const Array<Scalar>& time_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  )
{
  stepper_->addPoints(time_vec,x_vec,xdot_vec);
}


template<class Scalar> 
void SimpleIntegrator<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{
  stepper_->getPoints(time_vec,x_vec,xdot_vec,accuracy_vec);
}


template<class Scalar> 
TimeRange<Scalar> SimpleIntegrator<Scalar>::getTimeRange() const
{
  return stepper_->getTimeRange();
}


template<class Scalar> 
void SimpleIntegrator<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  stepper_->getNodes(time_vec);
}


template<class Scalar> 
void SimpleIntegrator<Scalar>::removeNodes(Array<Scalar>& time_vec)
{
  stepper_->removeNodes(time_vec);
}


template<class Scalar> 
int SimpleIntegrator<Scalar>::getOrder() const
{
  return stepper_->getOrder();
}


// private


template<class Scalar> 
bool SimpleIntegrator<Scalar>::advanceStepperToTime( const Scalar& advance_to_t )
{

#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Rythmos:SimpleIntegrator::advanceStepperToTime");
#endif

  using std::endl;
  typedef std::numeric_limits<Scalar> NL;
  using Teuchos::incrVerbLevel;
  using Teuchos::Describable;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);

  if (!is_null(integrationControlStrategy_)) {
    integrationControlStrategy_->setOStream(out);
    integrationControlStrategy_->setVerbLevel(incrVerbLevel(verbLevel,-1));
  }

  if (!is_null(integrationObserver_)) {
    integrationObserver_->setOStream(out);
    integrationObserver_->setVerbLevel(incrVerbLevel(verbLevel,-1));
  }

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nEntering " << this->Describable::description()
         << "::advanceStepperToTime("<<advance_to_t<<") ...\n";

  // Remember what timestep index we are on so we can report it later
  const int initCurrTimeStepIndex = currTimeStepIndex_;

  // Take steps until we the requested time is reached (or passed)

  TimeRange<Scalar> currStepperTimeRange = stepper_->getTimeRange();

  // Start by assume we can reach the time advance_to_t
  bool return_val = true;
  
  while ( !currStepperTimeRange.isInRange(advance_to_t) ) {

    // Halt immediatly if exceeded max iterations
    if (currTimeStepIndex_ >= maxNumTimeSteps_) {
      if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
        *out
          << "\n***"
          << "\n*** NOTICE: currTimeStepIndex = "<<currTimeStepIndex_
          << " >= maxNumTimeSteps = "<<maxNumTimeSteps_<< ", halting time integration!"
          << "\n***\n";
      return_val = false;
      break; // Exit the loop immediately!
    }

    if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
      *out << "\nTake step:  current_stepper_t = " << currStepperTimeRange.upper()
           << ", currTimeStepIndex = " << currTimeStepIndex_ << endl;
    Teuchos::OSTab tab(out);

    //
    // A) Reinitialize if a hard breakpoint was reached on the last time step
    //

    if (stepCtrlInfoLast_.limitedByBreakPoint) {
      if ( stepCtrlInfoLast_.breakPointType == BREAK_POINT_TYPE_HARD ) {
#ifdef ENABLE_RYTHMOS_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("Rythmos:SimpleIntegrator::restart");
#endif
        if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
          *out << "\nAt a hard-breakpoint, restarting time integrator ...\n";
        restart(&*stepper_);
      }
      else  {
        if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
          *out << "\nAt a soft-breakpoint, NOT restarting time integrator ...\n";
      }
    }

    //
    // B) Get the trial step control info
    //

    StepControlInfo<Scalar> trialStepCtrlInfo;
    {
#ifdef ENABLE_RYTHMOS_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:SimpleIntegrator::advanceStepperToTime: getStepCtrl");
#endif
      if (!is_null(integrationControlStrategy_)) {
        // Let an external strategy object determine the step size and type.
        // Note that any breakpoint info is also related through this call.
        trialStepCtrlInfo = integrationControlStrategy_->getNextStepControlInfo(
          *stepper_, stepCtrlInfoLast_, currTimeStepIndex_
          );
      }
      else {
        // Take a variable step if we have no control strategy
        trialStepCtrlInfo.stepType = STEP_TYPE_VARIABLE;
        trialStepCtrlInfo.stepSize = NL::max();
      }
    }

    // Print the initial trial step
    if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) ) {
      *out << "\nTrial step:\n";
      OSTab tab(out);
      *out << trialStepCtrlInfo;
    }

    // Halt immediately if we where told to do so
    if (trialStepCtrlInfo.stepSize < ST::zero()) {
      if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) )
        *out
          << "\n***"
          << "\n*** NOTICE: The IntegrationControlStrategy object return stepSize < 0.0, halting time integration!"
          << "\n***\n";
      return_val = false;
      break; // Exit the loop immediately!
    }

    // Make sure we don't step past the final time if asked not to
    bool updatedTrialStepCtrlInfo = false;
    {
      const Scalar finalTime = integrationTimeDomain_.upper();
      if (landOnFinalTime_ && trialStepCtrlInfo.stepSize + currStepperTimeRange.upper() > finalTime) {
        if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
          *out << "\nCutting trial step to avoid stepping past final time ...\n";
        trialStepCtrlInfo.stepSize = finalTime - currStepperTimeRange.upper();
        updatedTrialStepCtrlInfo = true;
      }
    }
    
    // Print the modified trial step
    if ( updatedTrialStepCtrlInfo
      && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) )
    {
      *out << "\nUpdated trial step:\n";
      OSTab tab(out);
      *out << trialStepCtrlInfo;
    }

    //
    // C) Take the step
    //

    // Print step type and size
    if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) ) {
      if (trialStepCtrlInfo.stepType == STEP_TYPE_VARIABLE)
        *out << "\nTaking a variable time step with max step size = "
             << trialStepCtrlInfo.stepSize << " ....\n";
      else
        *out << "\nTaking a fixed time step of size = "
             << trialStepCtrlInfo.stepSize << " ....\n";
    }

    // Take step
    Scalar stepSizeTaken;
    {
#ifdef ENABLE_RYTHMOS_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:SimpleIntegrator::advanceStepperToTime: takeStep");
#endif
      stepSizeTaken = stepper_->takeStep(
        trialStepCtrlInfo.stepSize, trialStepCtrlInfo.stepType
        );
    }

    // Validate step taken
    if (trialStepCtrlInfo.stepType == STEP_TYPE_VARIABLE) {
      TEST_FOR_EXCEPTION(
        stepSizeTaken < ST::zero(), std::logic_error,
        "Error, stepper took negative step of dt = " << stepSizeTaken << "!\n"
        );
      TEST_FOR_EXCEPTION(
        stepSizeTaken > trialStepCtrlInfo.stepSize, std::logic_error,
        "Error, stepper took step of dt = " << stepSizeTaken
        << " > max step size of = " << trialStepCtrlInfo.stepSize << "!\n"
        );
    }
    else { // STEP_TYPE_FIXED
      TEST_FOR_EXCEPTION(
        stepSizeTaken != trialStepCtrlInfo.stepSize, std::logic_error,
        "Error, stepper took step of dt = " << stepSizeTaken 
        << " when asked to take step of dt = " << trialStepCtrlInfo.stepSize << "\n"
        );
    }

    // Update info about this step
    currStepperTimeRange = stepper_->getTimeRange();
    const StepControlInfo<Scalar> stepCtrlInfo =
      stepCtrlInfoTaken(trialStepCtrlInfo,stepSizeTaken);

    // Print the step actually taken 
    if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) ) {
      *out << "\nStep actually taken:\n";
      OSTab tab(out);
      *out << stepCtrlInfo;
    }

    //
    // D) Output info about step
    //

    {

#ifdef ENABLE_RYTHMOS_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:SimpleIntegrator::advanceStepperToTime: output");
#endif
      
      // Print our own brief output
      if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) ) {
        StepStatus<Scalar> stepStatus = stepper_->getStepStatus();
        *out << "\nTime point reached = " << stepStatus.time << endl;
        *out << "\nstepStatus:\n" << stepStatus;
        if ( includesVerbLevel(verbLevel,Teuchos::VERB_EXTREME) ) {
          RCP<const Thyra::VectorBase<Scalar> >
            solution = stepStatus.solution,
            solutionDot = stepStatus.solutionDot;
          if (!is_null(solution))
            *out << "\nsolution = \n" << Teuchos::describe(*solution,verbLevel);
          if (!is_null(solutionDot))
            *out << "\nsolutionDot = \n" << Teuchos::describe(*solutionDot,verbLevel);
        }
      }
      
      // Output to the observer
      if (!is_null(integrationObserver_))
        integrationObserver_->observeCompletedTimeStep(
          *stepper_, stepCtrlInfo, currTimeStepIndex_
          );

    }

    //
    // E) Update info for next time step
    //

    stepCtrlInfoLast_ = stepCtrlInfo;
    ++currTimeStepIndex_;
    
  }

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nNumber of steps taken in this call to advanceStepperToTime(...) = "
         << (currTimeStepIndex_ - initCurrTimeStepIndex) << endl
         << "\nLeaving" << this->Describable::description()
         << "::advanceStepperToTime("<<advance_to_t<<") ...\n";

  return return_val;
  
}


} // namespace Rythmos


#endif //Rythmos_SIMPLE_INTEGRATOR_H
