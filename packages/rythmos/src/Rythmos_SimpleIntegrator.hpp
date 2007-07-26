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
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief A very simple concrete subclass for <tt>IntegratorBase</tt> that
 * allows just for simple fixed steps or variable steps.
 */
template<class Scalar> 
class SimpleIntegrator : virtual public IntegratorBase<Scalar>
{
public:
  
  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors, Initializers, Misc */
  //@{
  
  /** \brief . */
  SimpleIntegrator();

  //@}

  /** \name Overridden from IntegratorBase */
  //@{
  
  /** \brief . */
  void setInterpolationBuffer(
    const RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
    );

  /** \brief . */
  void setStepper(
    const RCP<StepperBase<Scalar> > &stepper,
    const Scalar &finalTime
    );

  /** \brief . */
  bool getFwdPoints(
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
  bool setPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec,
    const Array<ScalarMag> & accuracy_vec 
    );

  /** \brief . */
  bool getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const;

  /** \brief . */
  bool setRange(
    const TimeRange<Scalar>& range,
    const InterpolationBufferBase<Scalar>& interpBuffer
    );

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;

  /** \brief . */
  bool getNodes(Array<Scalar>* time_vec) const;

  /** \brief . */
  bool removeNodes(Array<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<ParameterList> getParameterList();

  /** \brief . */
  RCP<ParameterList> unsetParameterList();

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

private:

  // ////////////////////////
  // Private data members

  RCP<StepperBase<Scalar> > stepper_;
  RCP<ParameterList> paramList_;
  bool takeVariableSteps_;
  Scalar fixed_dt_;
  Scalar finalTime_;

  static const std::string takeVariableSteps_name_;
  static const bool takeVariableSteps_default_;

  static const std::string fixed_dt_name_;
  static const double fixed_dt_default_;

  // /////////////////////////
  // Private member functions

  void advanceStepperToTime( const Scalar& t );

};


/** \brief .
 *
 * \relates SimpleIntegrator
 */
template<class Scalar> 
RCP<SimpleIntegrator<Scalar> > simpleIntegrator()
{
  return Teuchos::rcp(new SimpleIntegrator<Scalar>());
}


/** \brief .
 *
 * \relates SimpleIntegrator
 */
template<class Scalar> 
RCP<SimpleIntegrator<Scalar> >
simpleIntegrator( const RCP<ParameterList> &paramList )
{
  RCP<SimpleIntegrator<Scalar> >
    integrator = Teuchos::rcp(new SimpleIntegrator<Scalar>());
  integrator->setParameterList(paramList);
  return integrator;
}


// ////////////////////////////
// Defintions


// Static members


template<class Scalar> 
const std::string
SimpleIntegrator<Scalar>::takeVariableSteps_name_ = "Take Variable Steps";

template<class Scalar> 
const bool
SimpleIntegrator<Scalar>::takeVariableSteps_default_ = true;


template<class Scalar> 
const std::string
SimpleIntegrator<Scalar>::fixed_dt_name_ = "Fixed dt";

template<class Scalar> 
const double
SimpleIntegrator<Scalar>::fixed_dt_default_ = -1.0;


// Constructors, Initializers, Misc


template<class Scalar>
SimpleIntegrator<Scalar>::SimpleIntegrator()
  :takeVariableSteps_(takeVariableSteps_default_),
   fixed_dt_(fixed_dt_default_),
   finalTime_(-std::numeric_limits<Scalar>::max())
{}


// Overridden from IntegratorBase


template<class Scalar>
void SimpleIntegrator<Scalar>::setInterpolationBuffer(
    const RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
    )
{
  TEST_FOR_EXCEPT(true); // We will never implement this function in this class!
}


template<class Scalar> 
void SimpleIntegrator<Scalar>::setStepper(
  const RCP<StepperBase<Scalar> > &stepper,
  const Scalar &finalTime
  )
{
  TEST_FOR_EXCEPT(is_null(stepper));
  TEST_FOR_EXCEPT( finalTime <= stepper->getTimeRange().lower() );
  // 2007/07/25: rabartl: ToDo: Validate state of the stepper!
  stepper_ = stepper;
  finalTime_ = finalTime;
}


template<class Scalar>
bool SimpleIntegrator<Scalar>::getFwdPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  )
{

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
    *out << "\nEntering " << this->Describable::description() << "::getFwdPoints(...) ...\n";

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) )
    *out << "\nRequested time points: " << Teuchos::toString(time_vec) << "\n";

  //
  // 1) Initial setup
  //

  const int numTimePoints = time_vec.size();

  // Assert preconditions
  assertTimePointsAreSorted(time_vec);
  TEST_FOR_EXCEPT(accuracy_vec!=0); // ToDo: Remove this!

  // Resize the storage for the output arrays
  if (x_vec)
    x_vec->resize(numTimePoints);
  if (xdot_vec)
    xdot_vec->resize(numTimePoints);

  // This index records the next time point in time_vec[timePointIndex] that
  // needs to be handled.  This gets updated in each step below
  int nextTimePointIndex = 0;
  
  assertNoTimePointsBeforeCurrentTimeRange(*this,time_vec,nextTimePointIndex);

  //
  // 1) Get all time points that fall within the current time range
  //

  getCurrentPoints(*stepper_,time_vec,x_vec,xdot_vec,&nextTimePointIndex);

  //
  // 2) Advance the stepper to satisfy time points in time_vec that fall
  // before the current time.
  //

  while ( nextTimePointIndex < numTimePoints ) {
    
    // Use the time stepping algorithm to step up to or past the next
    // requested time.
    const Scalar t = time_vec[nextTimePointIndex];
    advanceStepperToTime(t);
    
    // Extract the next set of points (perhaps just one) from the stepper
    getCurrentPoints(*stepper_,time_vec,x_vec,xdot_vec,&nextTimePointIndex);
    
  }

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nLeaving " << this->Describable::description() << "::getFwdPoints(...) ...\n";
  
  return true; // ToDo: Remove this!
  
}


template<class Scalar> 
TimeRange<Scalar>
SimpleIntegrator<Scalar>::getFwdTimeRange() const
{
  return timeRange(stepper_->getTimeRange().lower(),finalTime_);
}


// Overridden from InterpolationBufferBase


template<class Scalar> 
bool SimpleIntegrator<Scalar>::setPoints(
  const Array<Scalar>& time_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec,
  const Array<ScalarMag> & accuracy_vec 
  )
{
  return stepper_->setPoints(time_vec,x_vec,xdot_vec,accuracy_vec);
}


template<class Scalar> 
bool SimpleIntegrator<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{
  return stepper_->getPoints(time_vec,x_vec,xdot_vec,accuracy_vec);
}


template<class Scalar> 
bool SimpleIntegrator<Scalar>::setRange(
  const TimeRange<Scalar>& range,
  const InterpolationBufferBase<Scalar>& interpBuffer
  )
{
  return stepper_->setRange(range,interpBuffer);
}


template<class Scalar> 
TimeRange<Scalar> SimpleIntegrator<Scalar>::getTimeRange() const
{
  return stepper_->getTimeRange();
}


template<class Scalar> 
bool SimpleIntegrator<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  return stepper_->getNodes(time_vec);
}


template<class Scalar> 
bool SimpleIntegrator<Scalar>::removeNodes(Array<Scalar>& time_vec)
{
  return stepper_->removeNodes(time_vec);
}


template<class Scalar> 
int SimpleIntegrator<Scalar>::getOrder() const
{
  return stepper_->getOrder();
}


// Overridden from Teuchos::Describable


template<class Scalar> 
void SimpleIntegrator<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  TEST_FOR_EXCEPT(true);
}


// Overridden from ParameterListAcceptor


template<class Scalar> 
void SimpleIntegrator<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  using Teuchos::as;
  using Teuchos::get;
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters());
  paramList_ = paramList;
  takeVariableSteps_ = paramList_->get(
    takeVariableSteps_name_, takeVariableSteps_ );
  if (!takeVariableSteps_) {
    // The fixed_dt parameter must exsist!
    fixed_dt_ = get<double>(*paramList_,fixed_dt_name_);
    // ToDo: Generate a better error messae
    //Teuchos::ParameterEntry
    //  fixed_dt_entry = paramList_->getPtr(fixed_dt_name_);
  }
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
}


template<class Scalar> 
RCP<ParameterList> SimpleIntegrator<Scalar>::getParameterList()
{
  return paramList_;
}


template<class Scalar> 
RCP<ParameterList> SimpleIntegrator<Scalar>::unsetParameterList()
{
  RCP<ParameterList> tempParamList = paramList_;
  paramList_ = Teuchos::null;
  return tempParamList;
}


template<class Scalar> 
RCP<const ParameterList> SimpleIntegrator<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL) ) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(
      takeVariableSteps_name_, takeVariableSteps_default_,
      "Take variable time steps or fixed time steps.\n"
      "If set to false, then the parameter \"" + fixed_dt_name_ + "\"\n"
      "must be specified to give the size of the time steps."
      );
    pl->set(
      fixed_dt_name_, fixed_dt_default_,
      "Gives the size of the fixed time steps.  This is only read and used if\n"
      "\"" + takeVariableSteps_name_ + "\" is set to true."
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// private


template<class Scalar> 
void SimpleIntegrator<Scalar>::advanceStepperToTime( const Scalar& t )
{

  using std::endl;
  using Teuchos::incrVerbLevel;
  using Teuchos::Describable;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nEntering " << this->Describable::description()
         << "::advanceStepperToTime("<<t<<") ...\n";

  // Take steps until we the requested time is reached (or passed)

  TimeRange<Scalar> currentTimeRange = stepper_->getTimeRange();
  
  while ( !currentTimeRange.isInRange(t) ) {

    const Scalar current_t = currentTimeRange.upper();

    if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
      *out << "\nTime before time step is taken current_t = " << current_t << endl;
    Teuchos::OSTab tab(out);

    // Take a a variable or a fixed time step
    
    if (takeVariableSteps_) {

      const Scalar max_dt = finalTime_ - current_t;

      if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) )
        *out << "\nTaking a variable time step with max_dt = " << max_dt << " ....\n";

      const Scalar dt_taken = stepper_->takeStep(max_dt,VARIABLE_STEP);

      TEST_FOR_EXCEPTION(
        dt_taken < ST::zero(), std::logic_error,
        "Error, stepper took negative step of dt = " << dt_taken << "!\n"
        );
      TEST_FOR_EXCEPTION(
        dt_taken > max_dt, std::logic_error,
        "Error, stepper took step of dt = " << dt_taken
        << " > max_dt = " << max_dt << "!\n"
        );

    }
    else {

      if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) )
        *out << "\nTaking a fixed time step of size fixed_dt = " << fixed_dt_ << " ....\n";
      
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT( fixed_dt_ >= ST::zero() );
#endif

      const Scalar dt_taken = stepper_->takeStep(fixed_dt_,FIXED_STEP);

      TEST_FOR_EXCEPTION(
        dt_taken != fixed_dt_, std::logic_error,
        "Error, stepper took step of dt = " << dt_taken 
        << " when asked to take step of dt = " << fixed_dt_ << "\n"
        );

    }

    // Update the current time range that will be used below and for the next
    // loop control check!
    currentTimeRange = stepper_->getTimeRange();

    if ( includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM) )
      *out << "\nTime point reached = " << currentTimeRange.upper() << endl;
      
    if ( includesVerbLevel(verbLevel,Teuchos::VERB_EXTREME) ) {
      
      StepStatus<Scalar>
        stepStatus = stepper_->getStepStatus();
      
      RCP<const Thyra::VectorBase<Scalar> >
        solution = stepStatus.solution,
        solutionDot = stepStatus.solutionDot;

      *out << "\nsolution = \n" << Teuchos::describe(*solution,verbLevel);
      *out << "\nsolutionDot = \n" << Teuchos::describe(*solutionDot,verbLevel);

    }
    
  }

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nLeaving " << this->Describable::description()
         << "::advanceStepperToTime("<<t<<") ...\n";
  
}


} // namespace Rythmos


#endif //Rythmos_SIMPLE_INTEGRATOR_H
