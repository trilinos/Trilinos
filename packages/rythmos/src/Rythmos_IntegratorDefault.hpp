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

#ifndef Rythmos_INTEGRATOR_DEFAULT_H
#define Rythmos_INTEGRATOR_DEFAULT_H

#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_InterpolationBufferBase.hpp"
#include "Rythmos_InterpolationBufferAppender.hpp"
#include "Rythmos_BreakPointInformerBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Rythmos_InterpolationBufferHelpers.hpp"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief Concrete subclass of <tt>InterpolationBufferBase</tt> implemented in terms of
 * a <tt>StepperBase</tt> object and a trailing <tt>InterpolationBufferBase</tt> object.
 */
template<class Scalar> 
class IntegratorDefault : virtual public IntegratorBase<Scalar>
{
public:
  
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors, Initializers, Misc */
  //@{
  
  /** \brief . */
  ~IntegratorDefault() {};
 
  /** \brief . */
  IntegratorDefault();
  
  /** \brief . */
  IntegratorDefault(
    const Teuchos::RCP<StepperBase<Scalar> > &stepper,
    const Teuchos::RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer,
    const Teuchos::RCP<Teuchos::ParameterList> &paramList = Teuchos::null
    );
  
  /** \brief . */
  void setTrailingInterpolationBuffer(
    const Teuchos::RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
    );

  /** \brief . */
  bool acceptsTrailingInterpolationBuffer() const;

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

  /** \brief . */
  Teuchos::RCP<const InterpolationBufferBase<Scalar> > getInterpolationBuffer();

  /** \brief . */
  Teuchos::RCP<InterpolationBufferBase<Scalar> > unSetInterpolationBuffer();

  /** \brief . */
  void setInterpolationBufferAppender(
    const Teuchos::RCP<InterpolationBufferAppenderBase<Scalar> > &interpolationBufferAppender
    );

  /** \brief . */
  Teuchos::RCP<const InterpolationBufferAppenderBase<Scalar> > getInterpolationBufferAppender();

  /** \brief . */
  Teuchos::RCP<InterpolationBufferAppenderBase<Scalar> > unSetInterpolationBufferAppender();

  /** \brief . */
  void setStepper(
    const Teuchos::RCP<StepperBase<Scalar> > &stepper,
    const Scalar &finalTime,
    const bool landOnFinalTime
    );

  /** \brief . */
  Teuchos::RCP<const Rythmos::StepperBase<Scalar> > getStepper() const;

  /** \brief . */
  Teuchos::RCP<StepperBase<Scalar> > unSetStepper();
  
  /** \brief.  */
  void setBreakPointInformer(Teuchos::RCP<BreakPointInformerBase<Scalar> >& breakPointInformer);

  /** \brief.  */
  Teuchos::RCP<const BreakPointInformerBase<Scalar> > getBreakPointInformer();

  /** \brief.  */
  Teuchos::RCP<BreakPointInformerBase<Scalar> > unSetBreakPointInformer();
  
  /** \brief. */
  void setObserver(Teuchos::RCP<IntegrationObserverBase<Scalar> >& observer);

  /** \brief. */
  Teuchos::RCP<const IntegrationObserverBase<Scalar> > getObserver();

  /** \brief. */
  Teuchos::RCP<IntegrationObserverBase<Scalar> > unSetObserver();
  
  //@}
  
  void getFwdPoints(
    const Array<Scalar>& time_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    );

  /** \brief . */
  TimeRange<Scalar> getFwdTimeRange() const;
  

  /** \name Overridden from InterpolationBufferBase */
  //@{
    
  /** \brief . */
  void addPoints(
    const Array<Scalar>& time_vec,
    const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    );

  /** \brief . */
  void getPoints(
    const Array<Scalar>& time_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const;

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;

  /** \brief . */
  void getNodes(Array<Scalar>* time_vec) const;

  /** \brief . */
  int getOrder() const;

  //@}

  /** \brief . */
  virtual void removeNodes(Array<Scalar>& time_vec);

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getParameterList();

  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

private:


  // Interpolation Buffer used to store past results
  Teuchos::RCP<InterpolationBufferBase<Scalar> > trailingInterpBuffer_;

  // Interpolation Buffer Appender strategy object for transfering data from
  // the stepper to the traling interpolation buffer
  Teuchos::RCP<InterpolationBufferAppenderBase<Scalar> > interpolationBufferAppender_;

  // Stepper used to fill interpolation buffer.
  Teuchos::RCP<StepperBase<Scalar> > stepper_;

  // Stepper Observer 
  Teuchos::RCP<IntegrationObserverBase<Scalar> > observer_;

  // ParameterList to control behavior
  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  // BreakPoint informer strategy object
  // If this is null, then it has no effect
  Teuchos::RCP<BreakPointInformerBase<Scalar> > breakPointInformer_;

  // Take variable steps or not
  bool takeVariableSteps_;

  // The fixed step size to take
  Scalar fixed_dt_;

  // The final time to integrate to
  Scalar finalTime_;

  // Flag for whether the integrator is initialized or not.
  bool isInitialized_;

  // Private member functions:

  /** \brief . */
  void initialize_();

};


// ////////////////////////////
// Defintions


// Constructors, Initializers, Misc


template<class Scalar>
IntegratorDefault<Scalar>::IntegratorDefault()
  : takeVariableSteps_(true), fixed_dt_(-1.0), isInitialized_(false)
{}


template<class Scalar>
IntegratorDefault<Scalar>::IntegratorDefault(
  const Teuchos::RCP<StepperBase<Scalar> > &stepper,
  const Teuchos::RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer,
  const Teuchos::RCP<Teuchos::ParameterList> &paramList 
  )
  : takeVariableSteps_(true), fixed_dt_(-1.0), isInitialized_(false)
{
  using Teuchos::as;
  if (!is_null(paramList))
    setParameterList(paramList);
  const Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"IntegratorDefault::constructor");
  *out << "Initializing IntegratorDefault" << std::endl;
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH)) {
    *out << "Calling setStepper..." << std::endl;
  }
  setStepper(stepper,0.0,BREAK_POINT_TYPE_HARD);
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH)) {
    *out << "Calling setInterpolationBuffer..." << std::endl;
  }
  setTrailingInterpolationBuffer(trailingInterpBuffer);
  initialize_();
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > IntegratorDefault<Scalar>::get_x_space() const
{
  if (!isInitialized_) { 
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > space;
    return(space);
  }
  return(stepper_->get_x_space());
}

template<class Scalar>
void IntegratorDefault<Scalar>::initialize_()
{
  TEST_FOR_EXCEPTION(
      stepper_ == Teuchos::null, std::logic_error,
      "Error, Attempting to intialize integrator without setting a stepper!\n"
      );
  TEST_FOR_EXCEPTION(
      trailingInterpBuffer_ == Teuchos::null, std::logic_error,
      "Error, Attempting to intialize integrator without setting an interpolation buffer!\n"
      );
  if (paramList_ == Teuchos::null) {
    Teuchos::RCP<Teuchos::ParameterList> params;
    setParameterList(params);
  }
  if (interpolationBufferAppender_ == Teuchos::null) {
    Teuchos::RCP<InterpolationBufferAppenderDefault<Scalar> > 
      defaultInterpolationBufferAppender = Teuchos::rcp(new InterpolationBufferAppenderDefault<Scalar>);
    interpolationBufferAppender_ = defaultInterpolationBufferAppender;
  }
  isInitialized_ = true;
}

template<class Scalar>
void IntegratorDefault<Scalar>::setTrailingInterpolationBuffer(
  const Teuchos::RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
  )
{
#ifdef TEUCHOS_DEBUG
  // Check preconditions:
  TEST_FOR_EXCEPTION(trailingInterpBuffer == Teuchos::null,std::logic_error,"Error, trailingInterpBuffer == Teuchos::null!\n");
  if (stepper_ != Teuchos::null) {
    TEST_FOR_EXCEPTION(
        trailingInterpBuffer->getTimeRange().upper() != stepper_->getTimeRange().lower(),
        std::logic_error,
        "Error, specified interpolation buffer's upper time range = " <<
          trailingInterpBuffer->getTimeRange().upper() << 
          " does not match internal stepper's lower time range = " <<
          stepper_->getTimeRange().lower() << "!\n"
        );
  }
#endif // TEUCHOS_DEBUG
  using Teuchos::as;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  trailingInterpBuffer_ = trailingInterpBuffer;
  Teuchos::OSTab ostab(out,1,"IntegratorDefault::setInterpolationBuffer");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "trailingInterpBuffer_ = " << trailingInterpBuffer_->description() << std::endl;
  }
  if (stepper_ != Teuchos::null) {
    initialize_();
  }
}

template<class Scalar>
bool IntegratorDefault<Scalar>::acceptsTrailingInterpolationBuffer() const
{
  return(true);
}

template<class Scalar>
void IntegratorDefault<Scalar>::setStepper(
  const Teuchos::RCP<StepperBase<Scalar> > &stepper,
  const Scalar &finalTime,
  const bool landOnFinalTime // Not using this yet!
  )
{
#ifdef TEUCHOS_DEBUG
  // Check preconditions:
  TEST_FOR_EXCEPTION(stepper == Teuchos::null, std::logic_error, "Error, steppeer == Teuchos::null!\n");
  if (trailingInterpBuffer_ != Teuchos::null) {
    TEST_FOR_EXCEPTION(
        trailingInterpBuffer_->getTimeRange().upper() != stepper->getTimeRange().lower(),
        std::logic_error,
        "Error, specified stepper's lower time range = " <<
          stepper->getTimeRange().lower() <<
          " does not match internal trailing interpolation buffer's upper time range = " <<
          trailingInterpBuffer_->getTimeRange().upper() << "!\n"
        );
  }
#endif // TEUCHOS_DEBUG
  using Teuchos::as;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  stepper_ = stepper;
  Teuchos::OSTab ostab(out,1,"IntegratorDefault::setStepper");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "stepper_ = " << stepper_->description() << std::endl;
  }
  if (trailingInterpBuffer_ != Teuchos::null) {
    initialize_();
  }
}

template<class Scalar>
Teuchos::RCP<StepperBase<Scalar> > IntegratorDefault<Scalar>::unSetStepper()
{
  Teuchos::RCP<StepperBase<Scalar> > stepper_old = stepper_;
  stepper_ = Teuchos::null;
  isInitialized_ = false;
  return(stepper_old);
}

template<class Scalar>
void IntegratorDefault<Scalar>::setBreakPointInformer(
    Teuchos::RCP<BreakPointInformerBase<Scalar> >& breakPointInformer
    )
{
  breakPointInformer_ = breakPointInformer;
}


template<class Scalar>
Teuchos::RCP<const BreakPointInformerBase<Scalar> > IntegratorDefault<Scalar>::getBreakPointInformer()
{
  return(breakPointInformer_);
}


template<class Scalar>
Teuchos::RCP<BreakPointInformerBase<Scalar> > IntegratorDefault<Scalar>::unSetBreakPointInformer()
{
  Teuchos::RCP<BreakPointInformerBase<Scalar> > breakPointInformer_old = breakPointInformer_;
  breakPointInformer_ = Teuchos::null;
  return(breakPointInformer_old);
}

template<class Scalar>
Teuchos::RCP<const Rythmos::StepperBase<Scalar> > IntegratorDefault<Scalar>::getStepper() const
{
  return(stepper_);
}

template<class Scalar>
Teuchos::RCP<const InterpolationBufferBase<Scalar> > IntegratorDefault<Scalar>::getInterpolationBuffer()
{
  return(trailingInterpBuffer_);
}

template<class Scalar>
Teuchos::RCP<InterpolationBufferBase<Scalar> > IntegratorDefault<Scalar>::unSetInterpolationBuffer()
{
  Teuchos::RCP<InterpolationBufferBase<Scalar> > interpolationBuffer_old = trailingInterpBuffer_;
  trailingInterpBuffer_ = Teuchos::null;
  isInitialized_ = false;
  return(interpolationBuffer_old);
}

template<class Scalar>
void IntegratorDefault<Scalar>::setInterpolationBufferAppender(
    const Teuchos::RCP<InterpolationBufferAppenderBase<Scalar> > &interpolationBufferAppender
    )
{
  using Teuchos::as;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  interpolationBufferAppender_ = interpolationBufferAppender;
  Teuchos::OSTab ostab(out,1,"IntegratorDefault::setiInterpolationBufferAppender");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "interpolationBufferAppender_ = " << interpolationBufferAppender_->description() << std::endl;
  }
}

template<class Scalar>
Teuchos::RCP<const InterpolationBufferAppenderBase<Scalar> > IntegratorDefault<Scalar>::getInterpolationBufferAppender()
{
  return(interpolationBufferAppender_);
}

template<class Scalar>
Teuchos::RCP<InterpolationBufferAppenderBase<Scalar> > IntegratorDefault<Scalar>::unSetInterpolationBufferAppender()
{
  Teuchos::RCP<InterpolationBufferAppenderBase<Scalar> > interpolationBufferAppender_old = interpolationBufferAppender_;
  interpolationBufferAppender_ = Teuchos::null;
  isInitialized_ = false;
  return(interpolationBufferAppender_old);
}

template<class Scalar>
void IntegratorDefault<Scalar>::setObserver(Teuchos::RCP<IntegrationObserverBase<Scalar> >& observer)
{
  observer_ = observer;
}

template<class Scalar>
Teuchos::RCP<const IntegrationObserverBase<Scalar> > IntegratorDefault<Scalar>::getObserver()
{
  return(observer_);
}

template<class Scalar>
Teuchos::RCP<IntegrationObserverBase<Scalar> > IntegratorDefault<Scalar>::unSetObserver()
{
  Teuchos::RCP<IntegrationObserverBase<Scalar> > observer_old = observer_;
  observer_ = Teuchos::null;
  return(observer_old);
}
  

// Overridden from InterpolationBufferBase


template<class Scalar>
void IntegratorDefault<Scalar>::addPoints(
  const Array<Scalar>& time_vec
  ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec
  ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  ) 
{
  TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, addPoints is not defined for IntegratorDefault.\n"
      );
}


template<class Scalar>
void IntegratorDefault<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec_in,
  Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec_in,
  Array<ScalarMag>* accuracy_vec_in
  ) const
{
  TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,"Error, attemping to call getPoints (const) before initialization.\n");
#ifdef TEUCHOS_DEBUG
  // check preconditions:
  for (int i=0; i<Teuchos::as<int>(time_vec.size()) ; ++i) {
    TEST_FOR_EXCEPTION(
        ~(trailingInterpBuffer_->getTimeRange().isInRange(time_vec[i]) || stepper_->getTimeRange().isInRange(time_vec[i])),
        std::logic_error,
        "Error, time_vec[" << i << "] is not in TimeRange of trailing interpolation buffer = [" <<
          trailingInterpBuffer_->getTimeRange().lower() << "," <<
          trailingInterpBuffer_->getTimeRange().upper() << "] nor in TimeRange of stepper = [" << 
          stepper_->getTimeRange().lower() << "," << stepper_->getTimeRange().upper() << "]!\n"
        );
  }
#endif // TEUCHOS_DEBUG

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"IntegratorDefault::getPoints");

  // Dump the input

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "time_vec = " << std::endl;
    for (int i=0 ; i<Teuchos::as<int>(time_vec.size()) ; ++i) {
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
    }
    if (x_vec_in == NULL) {
      *out << "x_vec_in = NULL" << std::endl;
    }
    else if (x_vec_in->size() == 0) {
      *out << "x_vec_in = ptr to empty vector" << std::endl;
    }
    else {
      *out << "x_vec_in = " << std::endl;
      for (int i=0 ; i<Teuchos::as<int>(x_vec_in->size()) ; ++i) {
        *out << "x_vec[" << i << "] = " << std::endl;
        (*x_vec_in)[i]->describe(*out,Teuchos::VERB_EXTREME);
      }
    }
    if (xdot_vec_in == NULL) {
      *out << "xdot_vec_in = NULL" << std::endl;
    }
    else if (xdot_vec_in->size() == 0) {
      *out << "xdot_vec_in = ptr to empty vector" << std::endl;
    }
    else {
      *out << "xdot_vec = " << std::endl;
      for (int i=0 ; i<Teuchos::as<int>(xdot_vec_in->size()) ; ++i) {
        if ((*xdot_vec_in)[i] == Teuchos::null) {
          *out << "xdot_vec[" << i << "] = Teuchos::null" << std::endl;
        } else {
          *out << "xdot_vec[" << i << "] = " << std::endl;
          (*xdot_vec_in)[i]->describe(*out,Teuchos::VERB_EXTREME);
        }
      }
    }
    if (accuracy_vec_in == NULL) {
      *out << "accuracy_vec_in = NULL" << std::endl;
    }
    else if (accuracy_vec_in->size() == 0) {
      *out << "accuracy_vec_in = ptr to empty vector" << std::endl;
    }
    else {
      *out << "accuracy_vec = " << std::endl;
      for (int i=0 ; i<Teuchos::as<int>(accuracy_vec_in->size()) ; ++i) {
        *out << "accuracy_vec[" << i << "] = " << (*accuracy_vec_in)[i] << std::endl;
      }
    }
  }

  Array<Scalar> time_vec_local = time_vec;

  // First we get points from the stepper_.
  Array<Scalar> stepper_time_vec;
  Array<RCP<const Thyra::VectorBase<Scalar> > > stepper_x_vec, stepper_xdot_vec;
  Array<ScalarMag> stepper_accuracy_vec;

  selectPointsInTimeRange<Scalar>(&stepper_time_vec,time_vec_local,stepper_->getTimeRange());
  removePointsInTimeRange<Scalar>(&time_vec_local,stepper_->getTimeRange());
  stepper_->getPoints(
      stepper_time_vec, 
      &stepper_x_vec, 
      &stepper_xdot_vec, 
      &stepper_accuracy_vec
      );
  typename DataStore<Scalar>::DataStoreVector_t stepper_dsv;
  vectorToDataStoreVector(
      stepper_time_vec, 
      stepper_x_vec, 
      stepper_xdot_vec, 
      stepper_accuracy_vec, 
      &stepper_dsv
      );
  
  // Next we get points from the trailingInterpBuffer_
  Array<Scalar> IB_time_vec;
  Array<RCP<const Thyra::VectorBase<Scalar> > > IB_x_vec, IB_xdot_vec;
  Array<ScalarMag> IB_accuracy_vec;

  selectPointsInTimeRange<Scalar>(&IB_time_vec,time_vec_local,trailingInterpBuffer_->getTimeRange());
  removePointsInTimeRange<Scalar>(&time_vec_local,trailingInterpBuffer_->getTimeRange());
  trailingInterpBuffer_->getPoints(
      IB_time_vec, 
      &IB_x_vec, 
      &IB_xdot_vec, 
      &IB_accuracy_vec
      );
  typename DataStore<Scalar>::DataStoreVector_t IB_dsv;
  vectorToDataStoreVector(
      IB_time_vec, 
      IB_x_vec, 
      IB_xdot_vec, 
      IB_accuracy_vec, 
      &IB_dsv
      );

  TEST_FOR_EXCEPTION(
      time_vec_local.size() != 0, std::logic_error,
      "Error, there are " << time_vec_local.size() << 
      " points in time_vec that were not found in either the stepper_ or the trailingInterpBuffer_!\n"
      );

  int IB_N = IB_dsv.size();
  for (int i=0 ; i < IB_N ; ++i) {
    stepper_dsv.push_back(IB_dsv[i]);
  }

  std::sort(stepper_dsv.begin(),stepper_dsv.end());

  Array<Scalar> time_vec_out;
  dataStoreVectorToVector(
      stepper_dsv,
      &time_vec_out,
      x_vec_in,
      xdot_vec_in,
      accuracy_vec_in
      );

  TEST_FOR_EXCEPTION(
     time_vec.size() != time_vec_out.size(),
     std::logic_error,
     "Error, number of output points = " << 
       time_vec_out.size()  << " != " << 
       time_vec.size() << " = number of time point requested!\n"
     ); 
}

template<class Scalar>
TimeRange<Scalar> IntegratorDefault<Scalar>::getTimeRange() const
{
  if (!isInitialized_) {
    TimeRange<Scalar> range; // return invalid time range
    return(range);
  }
  TimeRange<Scalar> timerange(trailingInterpBuffer_->getTimeRange().lower(),stepper_->getTimeRange().upper());
  return(timerange);
}


template<class Scalar>
void IntegratorDefault<Scalar>::getNodes(
  Array<Scalar>* time_vec
  ) const
{
  TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,"Error, attempting to call getNodes before initialized!\n");
  using Teuchos::as;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"IntegratorDefault::getNodes");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << this->description() << std::endl;
  }
  return(trailingInterpBuffer_->getNodes(time_vec));
}


template<class Scalar>
void IntegratorDefault<Scalar>::removeNodes(
  Array<Scalar>& time_vec
  ) 
{
  TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,"Error, attempting to call removeNodes before initialized!\n");
  return(trailingInterpBuffer_->removeNodes(time_vec));
}


template<class Scalar>
int IntegratorDefault<Scalar>::getOrder() const
{
  TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,"Error, attempting to call getOrder before initialized!\n");
  return(trailingInterpBuffer_->getOrder());
}


// Overridden from Teuchos::Describable


template<class Scalar>
void IntegratorDefault<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{

  using Teuchos::as;
  using Teuchos::describe;

  const bool printSomething =
    (
      as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT)
      ||
      as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)
      );

  Teuchos::OSTab tab(out);

  if (printSomething) {
    out << this->description() << endl;
  }

  Teuchos::OSTab tab2(out);

  if (printSomething) {
    out << "trailing interpolation buffer = ";
    if (!is_null(trailingInterpBuffer_)) {
      out << describe(*trailingInterpBuffer_,verbLevel);
    } else {
      out << "NULL\n";
    }
    out << "stepper = ";
    if (!is_null(stepper_)) {
      out << describe(*stepper_,verbLevel);
    } else {
      out << "NULL\n";
    }
    out << "interpolationBufferAppender = ";
    if (!is_null(interpolationBufferAppender_)) {
      out << describe(*interpolationBufferAppender_,verbLevel);
    } else {
      out << "NULL\n";
    }
  }
  out << "paramList = ";
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM)) {
    if (!is_null(paramList_)) {
      out << paramList_->print(out);
    } else {
      out << "NULL\n";
    }
  }
  
}


// Overridden from Teuchos::ParameterListAcceptor


template <class Scalar>
void IntegratorDefault<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  using Teuchos::as;
  using Teuchos::get;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEST_FOR_EXCEPT(is_null(paramList));

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"IntegratorDefault::setParameterList");

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "paramList = " << paramList->print(*out) << std::endl;
  }

  paramList->validateParametersAndSetDefaults(*getValidParameters());
  paramList_ = paramList;

  takeVariableSteps_ = get<bool>(*paramList_,"Take Variable Steps");
  fixed_dt_ = get<Scalar>(*paramList_,"fixed_dt");
  finalTime_ = get<Scalar>(*paramList_,"Final Time");

  TEST_FOR_EXCEPTION(
    !takeVariableSteps_ && fixed_dt_ <= ST::zero(),
    Teuchos::Exceptions::InvalidParameterValue,
    "Error, if we are taking fixed steps then \"fixed_dt\" must be set to"
    " a positive number!"
    );

  Teuchos::readVerboseObjectSublist(&*paramList_,this);

#ifdef TEUCHOS_DEBUG
  paramList->validateParameters(*getValidParameters());
#endif
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorDefault<Scalar>::getParameterList()
{
  return(paramList_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorDefault<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = paramList_;
  paramList_ = Teuchos::null;
  return(temp_param_list);
}


template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
IntegratorDefault<Scalar>::getValidParameters() const
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(
      "Take Variable Steps", bool(true),
      "If set to true, then the stepper will take variable steps.\n"
      "If set to false, then fixed stepps will be taken of size\n"
      "\"fixed_dt\".  Warning, you must set \"fixed_dt\" in this case!"
      );
    pl->set(
      "fixed_dt", Scalar(-1.0),
      "If fixed steps are being taken, then this is the step size."
      );
    pl->set(
      "Final Time", Scalar(0.0),
      "This specifies the final time to which the integrator will integrate.\n"
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return (validPL);
}

/* 08/13/07 tscoffe:  To simplify getPoints and allow the time points to be
 * unsorted & to come from current values and forward values, I can write a
 * helper function to do all of this.
 */
template <class Scalar>
void IntegratorDefault<Scalar>::getFwdPoints(
    const Array<Scalar>& time_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  if (!isInitialized_) {
    initialize_();
  }
#ifdef TEUCHOS_DEBUG
  // Check preconditions:
  for (int i=0 ; i<Teuchos::as<int>(time_vec.size()) ; ++i) { 
    TEST_FOR_EXCEPTION(
        isInRange_oc(this->getFwdTimeRange(),time_vec[i]),
        std::logic_error,
        "Error, time_vec[" << i << "] = " << time_vec[i] << 
          " is not in this->getFwdTimeRange:  (" << 
          this->getFwdTimeRange().lower() << "," <<
          this->getFwdTimeRange().upper() << "]!\n"
        );
  }
  // Check that time_vec is sorted:
  assertTimePointsAreSorted(time_vec);
#endif // TEUCHOS_DEBUG

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"IntegratorDefault::getFwdPoints");

  Array<Scalar> time_vec_local = time_vec;
  Array<Scalar> temp_time_vec;
  
  // First we use this->getPoints to recover points in the trailingInterpBuffer_ and the stepper_.
  selectPointsInTimeRange<Scalar>(&temp_time_vec,time_vec_local,this->getTimeRange());
  removePointsInTimeRange<Scalar>(&time_vec_local,this->getTimeRange());
  this->getPoints(
      temp_time_vec, 
      x_vec, 
      xdot_vec, 
      accuracy_vec
      );
   
  // Then we integrate forward for new points
  
  // Assumptions:
  // stepper_ has initial condition and is ready to take steps.
  // trailingInterpBuffer_ is initialized and is ready for points to be imported to it from stepper_.
  //
  // 08/13/07 tscoffe:  Question:  We're getting RCPs from the Stepper here and
  // stuffing them into the trailing interpolation buffer, what if the stepper
  // decides to re-use the vector sitting bethind the RCP?  It will corrupt the
  // data in the trailing interpolation buffer.  Oof.
  // This is okay for ImplicitBDFStepper, but we need to define this behavior
  // at the StepperBase level.  TODO

  RCP<Array<RCP<const Thyra::VectorBase<Scalar> > > > local_x_vec, local_xdot_vec;
  RCP<Array<ScalarMag> > local_accuracy_vec;

  if (x_vec == 0) { 
    local_x_vec = Teuchos::null; 
  } else {
    local_x_vec = rcp(new Array<RCP<const Thyra::VectorBase<Scalar> > >());
  }
  if (xdot_vec == 0) { 
    local_xdot_vec = Teuchos::null; 
  } else {
    local_xdot_vec = rcp(new Array<RCP<const Thyra::VectorBase<Scalar> > >());
  }
  if (accuracy_vec == 0) { 
    local_accuracy_vec = Teuchos::null; 
  } else {
    local_accuracy_vec = rcp(new Array<ScalarMag>());
  }

  for (int i=0 ; i<Teuchos::as<int>(time_vec_local.size()) ; ++i) {
    // 08/13/07 tscoffe:  we're going to grab points one at a time from the stepper.

    while (!(stepper_->getTimeRange().isInRange(time_vec_local[i]))) {
      if (takeVariableSteps_) {
        Scalar step_taken = stepper_->takeStep(ST::zero(),STEP_TYPE_VARIABLE);
        if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH)) {
          *out << "step_taken = " << step_taken << std::endl;
        }
      } 
      else {
        Scalar step_taken = stepper_->takeStep(fixed_dt_,STEP_TYPE_FIXED);
        if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH)) {
          *out << "step_taken = " << step_taken << std::endl;
        }
      }
      interpolationBufferAppender_->import(&*trailingInterpBuffer_,*stepper_,stepper_->getTimeRange());
      // After we've stepped forward & transfered the data into the trailing interpolation buffer, call the observer.
      if (observer_ != Teuchos::null)
      {
        TEST_FOR_EXCEPT(true); // ToDo: Call the new interface!
        //observer_->notify(*stepper_);
      } 
    }
    Array<Scalar> temp_time_vec;
    temp_time_vec.push_back(time_vec_local[i]);

    stepper_->getPoints(temp_time_vec,&*local_x_vec,&*local_xdot_vec,&*local_accuracy_vec);

    if (x_vec) {
      x_vec->push_back((*local_x_vec)[0]);
    }
    if (xdot_vec) {
      xdot_vec->push_back((*local_xdot_vec)[0]);
    }
    if (accuracy_vec) {
      accuracy_vec->push_back((*local_accuracy_vec)[0]);
    }
  }

}

template <class Scalar>
TimeRange<Scalar> IntegratorDefault<Scalar>::getFwdTimeRange() const
{
  if (!isInitialized_) {
    TimeRange<Scalar> range; // invalid time range.
    return(range);
  }
  TimeRange<Scalar> timerange(trailingInterpBuffer_->getTimeRange().lower(),finalTime_);
  return(timerange);
}

} // namespace Rythmos


#endif //Rythmos_INTEGRATOR_DEFAULT_H
