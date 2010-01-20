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

#ifndef RYTHMOS_DEFAULT_INTEGRATOR_DEF_HPP
#define RYTHMOS_DEFAULT_INTEGRATOR_DEF_HPP

#include "Rythmos_DefaultIntegrator_decl.hpp"
#include "Rythmos_InterpolationBufferHelpers.hpp"
#include "Rythmos_IntegrationControlStrategyBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_InterpolationBufferAppenderBase.hpp"
#include "Rythmos_PointwiseInterpolationBufferAppender.hpp"
#include "Rythmos_StepperHelpers.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

namespace Rythmos {

/** \brief .
 *
 * \relates DefaultIntegrator
 */
template<class Scalar> 
RCP<DefaultIntegrator<Scalar> >
defaultIntegrator()
{
  RCP<DefaultIntegrator<Scalar> >
    integrator = Teuchos::rcp(new DefaultIntegrator<Scalar>());
  return integrator;
}


/** \brief .
 *
 * \relates DefaultIntegrator
 */
template<class Scalar> 
RCP<DefaultIntegrator<Scalar> >
defaultIntegrator(
  const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy,
  const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
  )
{
  RCP<DefaultIntegrator<Scalar> >
    integrator = Teuchos::rcp(new DefaultIntegrator<Scalar>());
  integrator->setIntegrationControlStrategy(integrationControlStrategy);
  integrator->setIntegrationObserver(integrationObserver);
  return integrator;
}


/** \brief .
 *
 * \relates DefaultIntegrator
 */
template<class Scalar> 
RCP<DefaultIntegrator<Scalar> >
controlledDefaultIntegrator(
  const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy
  )
{
  RCP<DefaultIntegrator<Scalar> >
    integrator = Teuchos::rcp(new DefaultIntegrator<Scalar>());
  integrator->setIntegrationControlStrategy(integrationControlStrategy);
  return integrator;
}


/** \brief .
 *
 * \relates DefaultIntegrator
 */
template<class Scalar> 
RCP<DefaultIntegrator<Scalar> >
observedDefaultIntegrator(
  const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
  )
{
  RCP<DefaultIntegrator<Scalar> >
    integrator = Teuchos::rcp(new DefaultIntegrator<Scalar>());
  integrator->setIntegrationObserver(integrationObserver);
  return integrator;
}



// ////////////////////////////
// Defintions


// Static data members


template<class Scalar>
const std::string
DefaultIntegrator<Scalar>::maxNumTimeSteps_name_ = "Max Number Time Steps";

template<class Scalar>
const int
DefaultIntegrator<Scalar>::maxNumTimeSteps_default_ = 10000;


// Constructors, Initializers, Misc


template<class Scalar>
DefaultIntegrator<Scalar>::DefaultIntegrator()
  :landOnFinalTime_(true),
   maxNumTimeSteps_(maxNumTimeSteps_default_),
   currTimeStepIndex_(-1)
{}


template<class Scalar>
void DefaultIntegrator<Scalar>::setIntegrationControlStrategy(
  const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy
  )
{
#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT(is_null(integrationControlStrategy));
#endif
  integrationControlStrategy_ = integrationControlStrategy;
}

template<class Scalar>
RCP<IntegrationControlStrategyBase<Scalar> > 
  DefaultIntegrator<Scalar>::getNonconstIntegrationControlStrategy()
{
  return integrationControlStrategy_;
}

template<class Scalar>
RCP<const IntegrationControlStrategyBase<Scalar> > 
  DefaultIntegrator<Scalar>::getIntegrationControlStrategy() const
{
  return integrationControlStrategy_;
}


template<class Scalar>
void DefaultIntegrator<Scalar>::setIntegrationObserver(
  const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
  )
{
#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT(is_null(integrationObserver));
#endif
  integrationObserver_ = integrationObserver;
}


template<class Scalar> 
void DefaultIntegrator<Scalar>::setInterpolationBufferAppender(
  const RCP<InterpolationBufferAppenderBase<Scalar> > &interpBufferAppender
  )
{
  interpBufferAppender_ = interpBufferAppender.assert_not_null();
}


template<class Scalar> 
RCP<const InterpolationBufferAppenderBase<Scalar> >
DefaultIntegrator<Scalar>::getInterpolationBufferAppender()
{
  return interpBufferAppender_;
}

template<class Scalar> 
RCP<InterpolationBufferAppenderBase<Scalar> >
DefaultIntegrator<Scalar>::getNonconstInterpolationBufferAppender()
{
  return interpBufferAppender_;
}

template<class Scalar> 
RCP<InterpolationBufferAppenderBase<Scalar> >
DefaultIntegrator<Scalar>::unSetInterpolationBufferAppender()
{
  RCP<InterpolationBufferAppenderBase<Scalar> > InterpBufferAppender;
  std::swap( InterpBufferAppender, interpBufferAppender_ );
  return InterpBufferAppender;
}


// Overridden from ParameterListAcceptor


template<class Scalar> 
void DefaultIntegrator<Scalar>::setParameterList(
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
DefaultIntegrator<Scalar>::getValidParameters() const
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
DefaultIntegrator<Scalar>::cloneIntegrator() const
{

  using Teuchos::null;
  RCP<DefaultIntegrator<Scalar> >
    newIntegrator = Teuchos::rcp(new DefaultIntegrator<Scalar>());
  // Only copy control information, not the stepper or the model it contains!
  newIntegrator->stepper_ = Teuchos::null;
  const RCP<const ParameterList> paramList = this->getParameterList();
  if (!is_null(paramList)) {
    newIntegrator->setParameterList(Teuchos::parameterList(*paramList));
  }
  if (!is_null(integrationControlStrategy_)) {
    newIntegrator->integrationControlStrategy_ =
      integrationControlStrategy_->cloneIntegrationControlStrategy().assert_not_null();
  }
  if (!is_null(integrationObserver_)) {
    newIntegrator->integrationObserver_ =
      integrationObserver_->cloneIntegrationObserver().assert_not_null();
  }
  if (!is_null(trailingInterpBuffer_)) {
    // ToDo: implement the clone!
    newIntegrator->trailingInterpBuffer_ = null;
    //newIntegrator->trailingInterpBuffer_ =
    //  trailingInterpBuffer_->cloneInterploationBuffer().assert_not_null();
  }
  if (!is_null(interpBufferAppender_)) {
    // ToDo: implement the clone!
    newIntegrator->interpBufferAppender_ = null;
    //newIntegrator->interpBufferAppender_ =
    //  interpBufferAppender_->cloneInterpolationBufferAppender().assert_not_null();
  }
  return newIntegrator;
}


template<class Scalar> 
void DefaultIntegrator<Scalar>::setStepper(
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
RCP<StepperBase<Scalar> > DefaultIntegrator<Scalar>::unSetStepper()
{
  RCP<StepperBase<Scalar> > stepper_temp = stepper_;
  stepper_ = Teuchos::null;
  return(stepper_temp);
}


template<class Scalar>
RCP<const StepperBase<Scalar> > DefaultIntegrator<Scalar>::getStepper() const
{
  return(stepper_);
}


template<class Scalar>
RCP<StepperBase<Scalar> > DefaultIntegrator<Scalar>::getNonconstStepper() const
{
  return(stepper_);
}


template<class Scalar>
void DefaultIntegrator<Scalar>::setTrailingInterpolationBuffer(
  const RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
  )
{
  trailingInterpBuffer_ = trailingInterpBuffer.assert_not_null();
}


template<class Scalar>
RCP<InterpolationBufferBase<Scalar> >
DefaultIntegrator<Scalar>::getNonconstTrailingInterpolationBuffer()
{
  return trailingInterpBuffer_;
}


template<class Scalar>
RCP<const InterpolationBufferBase<Scalar> >
DefaultIntegrator<Scalar>::getTrailingInterpolationBuffer() const
{
  return trailingInterpBuffer_;
}

template<class Scalar>
RCP<InterpolationBufferBase<Scalar> >
DefaultIntegrator<Scalar>::unSetTrailingInterpolationBuffer()
{
  RCP<InterpolationBufferBase<Scalar> > trailingInterpBuffer;
  std::swap( trailingInterpBuffer, trailingInterpBuffer_ );
  return trailingInterpBuffer;
}


template<class Scalar>
void DefaultIntegrator<Scalar>::getFwdPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  )
{

#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR_DIFF("Rythmos:DefaultIntegrator::getFwdPoints",
    TopLevel);
#endif

  using Teuchos::incrVerbLevel;
#ifndef _MSC_VER
  using Teuchos::Describable;
#endif
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef InterpolationBufferBase<Scalar> IBB;
  typedef Teuchos::VerboseObjectTempState<IBB> VOTSIBB;

  finalizeSetup();

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  VOTSIBB stepper_outputTempState(stepper_,out,incrVerbLevel(verbLevel,-1));

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nEntering " << this->Describable::description() << "::getFwdPoints(...) ...\n"
         << "\nStepper: " << Teuchos::describe(*stepper_,verbLevel);

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) )
    *out << "\nRequested time points: " << Teuchos::toString(time_vec) << "\n";

  //
  // 0) Initial setup
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
    TEUCHOS_FUNC_TIME_MONITOR("Rythmos:DefaultIntegrator::getFwdPoints: getPoints");
#endif
    // 2007/10/05: rabartl: ToDo: Get points from trailingInterpBuffer_ first!
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
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:DefaultIntegrator::getFwdPoints: advanceStepperToTime");
#endif
      advanceStepperToTimeSucceeded= advanceStepperToTime(t);
    }
    if (!advanceStepperToTimeSucceeded) {
      bool reachedMaxNumTimeSteps = (currTimeStepIndex_ >= maxNumTimeSteps_);
      if (reachedMaxNumTimeSteps) {
        // Break out of the while loop and attempt to exit gracefully.
        break;
      }
      TEST_FOR_EXCEPTION(
          !advanceStepperToTimeSucceeded, Exceptions::GetFwdPointsFailed,
          this->description() << "\n\n"
          "Error:  The integration failed to get to time " << t << " and only achieved\n"
          "getting to " << stepper_->getTimeRange().upper() << "!"
          );
    }
    
    // Extract the next set of points (perhaps just one) from the stepper
    {
#ifdef ENABLE_RYTHMOS_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:DefaultIntegrator::getFwdPoints: getPoints (fwd)");
#endif
      getCurrentPoints(*stepper_,time_vec,x_vec,xdot_vec,&nextTimePointIndex);
    }
    
  }

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nLeaving " << this->Describable::description() << "::getFwdPoints(...) ...\n";
  
}


template<class Scalar> 
TimeRange<Scalar>
DefaultIntegrator<Scalar>::getFwdTimeRange() const
{
  return timeRange(
    stepper_->getTimeRange().lower(),
    integrationTimeDomain_.upper()
    );
}


// Overridden from InterpolationBufferBase


template<class Scalar> 
RCP<const Thyra::VectorSpaceBase<Scalar> >
DefaultIntegrator<Scalar>::get_x_space() const
{
  return stepper_->get_x_space();
}


template<class Scalar> 
void DefaultIntegrator<Scalar>::addPoints(
  const Array<Scalar>& time_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  )
{
  stepper_->addPoints(time_vec,x_vec,xdot_vec);
}


template<class Scalar> 
void DefaultIntegrator<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{
//  if (nonnull(trailingInterpBuffer_)) {
//    int nextTimePointIndex = 0;
//    getCurrentPoints(*trailingInterpBuffer_, time_vec, x_vec, xdot_vec, &nextTimePointIndex);
//    getCurrentPoints(*stepper_, time_vec, x_vec, xdot_vec, &nextTimePointIndex);
//    TEST_FOR_EXCEPTION( nextTimePointIndex < Teuchos::as<int>(time_vec.size()),
//      std::out_of_range,
//      "Error, the time point time_vec["<<nextTimePointIndex<<"] = "
//      << time_vec[nextTimePointIndex] << " falls outside of the time range "
//      << getTimeRange() << "!"
//      );
//  }
//  else {
  stepper_->getPoints(time_vec, x_vec, xdot_vec, accuracy_vec);
//  }
}


template<class Scalar> 
TimeRange<Scalar> DefaultIntegrator<Scalar>::getTimeRange() const
{
  if (nonnull(trailingInterpBuffer_)) {
    return timeRange(trailingInterpBuffer_->getTimeRange().lower(),
      stepper_->getTimeRange().upper());
  }
  return stepper_->getTimeRange();
}


template<class Scalar> 
void DefaultIntegrator<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  stepper_->getNodes(time_vec);
}


template<class Scalar> 
void DefaultIntegrator<Scalar>::removeNodes(Array<Scalar>& time_vec)
{
  stepper_->removeNodes(time_vec);
}


template<class Scalar> 
int DefaultIntegrator<Scalar>::getOrder() const
{
  return stepper_->getOrder();
}


// private


template<class Scalar> 
void DefaultIntegrator<Scalar>::finalizeSetup()
{
  if (!is_null(trailingInterpBuffer_) && is_null(interpBufferAppender_))
    interpBufferAppender_ = pointwiseInterpolationBufferAppender<Scalar>();
  // ToDo: Do other setup?
}


template<class Scalar> 
bool DefaultIntegrator<Scalar>::advanceStepperToTime( const Scalar& advance_to_t )
{

#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR_DIFF("Rythmos:DefaultIntegrator::advanceStepperToTime",
    TopLevel);
#endif

  using std::endl;
  typedef std::numeric_limits<Scalar> NL;
  using Teuchos::incrVerbLevel; 
#ifndef _MSC_VER
  using Teuchos::Describable;
#endif
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

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

    //
    // A) Reinitialize if a hard breakpoint was reached on the last time step
    //

    if (stepCtrlInfoLast_.limitedByBreakPoint) {
      if ( stepCtrlInfoLast_.breakPointType == BREAK_POINT_TYPE_HARD ) {
#ifdef ENABLE_RYTHMOS_TIMERS
        TEUCHOS_FUNC_TIME_MONITOR("Rythmos:DefaultIntegrator::restart");
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
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:DefaultIntegrator::advanceStepperToTime: getStepCtrl");
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
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:DefaultIntegrator::advanceStepperToTime: takeStep");
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

    // Append the trailing interpolation buffer (if defined)
    if (!is_null(trailingInterpBuffer_)) {
      interpBufferAppender_->append(*stepper_,currStepperTimeRange,
        trailingInterpBuffer_.ptr() );
    }

    //
    // D) Output info about step
    //

    {

#ifdef ENABLE_RYTHMOS_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("Rythmos:DefaultIntegrator::advanceStepperToTime: output");
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

// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_DEFAULT_INTEGRATOR_INSTANT(SCALAR) \
  \
  template class DefaultIntegrator< SCALAR >; \
  \
  template RCP<DefaultIntegrator< SCALAR > > \
  defaultIntegrator(); \
  \
  template RCP<DefaultIntegrator< SCALAR > > \
  defaultIntegrator( \
    const RCP<IntegrationControlStrategyBase< SCALAR > > &integrationControlStrategy, \
    const RCP<IntegrationObserverBase< SCALAR > > &integrationObserver \
    ); \
  \
  template RCP<DefaultIntegrator< SCALAR > > \
  controlledDefaultIntegrator( \
    const RCP<IntegrationControlStrategyBase< SCALAR > > &integrationControlStrategy \
    ); 


} // namespace Rythmos


#endif //RYTHMOS_DEFAULT_INTEGRATOR_DEF_HPP
