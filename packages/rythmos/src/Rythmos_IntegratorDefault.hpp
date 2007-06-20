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
#include "Rythmos_StepperBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief Concrete subclass of <tt>InterpolationBufferBase</tt> implemented in terms of
 * a <tt>StepperBase</tt> object and a trailing <tt>InterpolationBufferBase</tt> object.
 *
 * This class is really the beginnings of a 
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
  void setInterpolationBuffer(
    const Teuchos::RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
    );

  /** \brief . */
  void setStepper(
    const Teuchos::RCP<StepperBase<Scalar> > &stepper_
    );

  //@}
  
  /** \brief . */
  bool getFwdPoints(
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
  bool setPoints(
    const Array<Scalar>& time_vec,
    const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec,
    const Array<ScalarMag> & accuracy_vec 
    );

  /** \brief . */
  bool getPoints(
    const Array<Scalar>& time_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
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
  virtual bool removeNodes(Array<Scalar>& time_vec);

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

  // Stepper used to fill interpolation buffer.
  Teuchos::RCP<StepperBase<Scalar> > stepper_;

  // ParameterList to control behavior
  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  // Take variable steps or not
  bool takeVariableSteps_;

  // The fixed step size to take
  Scalar fixed_dt_;

};


// ////////////////////////////
// Defintions


// Constructors, Initializers, Misc


template<class Scalar>
IntegratorDefault<Scalar>::IntegratorDefault()
  : takeVariableSteps_(true), fixed_dt_(-1.0)
{}


template<class Scalar>
IntegratorDefault<Scalar>::IntegratorDefault(
  const Teuchos::RCP<StepperBase<Scalar> > &stepper,
  const Teuchos::RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer,
  const Teuchos::RCP<Teuchos::ParameterList> &paramList 
  )
  : takeVariableSteps_(true), fixed_dt_(-1.0)
{
  using Teuchos::as;
  if (!is_null(paramList))
    setParameterList(paramList);
  const Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  out->precision(20);
  out->setMaxLenLinePrefix(40);
  //out->pushLinePrefix("IntegratorDefault");
  //out->setShowLinePrefix(true);
  //out->setTabIndentStr("    ");
  Teuchos::OSTab ostab(out,1,"IBAS::constructor");
  *out << "Initializing IntegratorDefault" << std::endl;
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH)) {
    *out << "Calling setStepper..." << std::endl;
  }
  setStepper(stepper);
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH)) {
    *out << "Calling setInterpolationBuffer..." << std::endl;
  }
  setInterpolationBuffer(trailingInterpBuffer);
}


template<class Scalar>
void IntegratorDefault<Scalar>::setInterpolationBuffer(
  const Teuchos::RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
  )
{
  using Teuchos::as;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  // 10/9/06 tscoffe: What should we do if this is called after
  // initialization?  Basically, you're swapping out the history for a new
  // one.  This could be the result of upgrading or downgrading the accuracy
  // of the buffer, or something I haven't thought of yet.  Since
  // trailingInterpBuffer_'s node_vec is checked each time getPoints is
  // called, this should be fine.  And the time values in
  // trailingInterpBuffer_ need not be synchronized with stepper_.  Note also:
  // this functionality is important for checkpointing.
  trailingInterpBuffer_ = trailingInterpBuffer;
  Teuchos::OSTab ostab(out,1,"IBAS::setInterpolationBuffer");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "trailingInterpBuffer_ = " << trailingInterpBuffer_->description() << std::endl;
  }
}


template<class Scalar>
void IntegratorDefault<Scalar>::setStepper(
  const Teuchos::RCP<StepperBase<Scalar> > &stepper
  )
{
  using Teuchos::as;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  // 10/9/06 tscoffe: What should we do if this is called after
  // initialization?  Basically, you're swapping out the stepper for a new
  // one.  Since we're copying the data into trailingInterpBuffer_ after each
  // stepper_->takeStep() call, this should be fine, and it will essentially
  // result in changing the stepper_ mid-stream.  If the new stepper_ has a
  // time value before the end of the data in trailingInterpBuffer_, then you
  // will not get new stepper_ data until you ask for time values after the
  // end of trailingInterpBuffer_'s data.  And then the stepper_ will walk
  // forward inserting new (potentially inconsistent) data into
  // trailingInterpBuffer_ until it can give you the time values you asked
  // for.  Then trailingInterpBuffer_ will potentially have old and new data
  // in it.  On the other hand, if you swap out the stepper_ and the time
  // value is synchronized with the old stepper_, then you will essentially
  // change the integrator mid-stream and everything should proceed without
  // problems.  Note also: this functionality is important for checkpointing.
  stepper_ = stepper;
  Teuchos::OSTab ostab(out,1,"IBAS::setStepper");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "stepper_ = " << stepper_->description() << std::endl;
  }
}


// Overridden from InterpolationBufferBase


template<class Scalar>
bool IntegratorDefault<Scalar>::setPoints(
  const Array<Scalar>& time_vec
  ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec
  ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  ,const Array<ScalarMag> & accuracy_vec 
  ) 
{
  return(trailingInterpBuffer_->setPoints(time_vec,x_vec,xdot_vec,accuracy_vec));
}


template<class Scalar>
bool IntegratorDefault<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec_in,
  Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec_in,
  Array<ScalarMag>* accuracy_vec_in
  ) const
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"IBAS::getPoints");

  // Dump the input

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "time_vec = " << std::endl;
    for (unsigned int i=0 ; i<time_vec.size() ; ++i) {
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
      for (unsigned int i=0 ; i<x_vec_in->size() ; ++i) {
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
      for (unsigned int i=0 ; i<xdot_vec_in->size() ; ++i) {
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
      for (unsigned int i=0 ; i<accuracy_vec_in->size() ; ++i) {
        *out << "accuracy_vec[" << i << "] = " << (*accuracy_vec_in)[i] << std::endl;
      }
    }
  }

  // See if all of the time points are in the trailing interpolation buffer
  // and return them if they are

  bool status = trailingInterpBuffer_->getPoints(
    time_vec,x_vec_in,xdot_vec_in,accuracy_vec_in );
  if (status) {
    return(status); 
  }
  status = true;

  // All of the time points where not in the trailing interpolation buffer so
  // we have to get some from perhaps the stepper too.

  x_vec_in->clear();
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "trailingInterpBuffer_->getPoints unsuccessful" << std::endl;
  }

  TEST_FOR_EXCEPT(0 == x_vec_in);
  TEST_FOR_EXCEPT(0 == xdot_vec_in);
  TEST_FOR_EXCEPT(0 == accuracy_vec_in);
  // 2007/06/08: rabartl: ToDo: This code should be updated to allow any and
  // all of these pointers to be null in which case the will be ignored.

  Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &x_vec = *x_vec_in;
  Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &xdot_vec = *xdot_vec_in;
  Array<ScalarMag> &accuracy_vec = *accuracy_vec_in;

  // Sort time_vec
  Array<Scalar> local_time_vec = time_vec;
  std::sort(local_time_vec.begin(),local_time_vec.end());
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "sorted local time_vec:" << std::endl;
    for (unsigned int i=0 ; i<local_time_vec.size() ; ++i) {
      *out << "local_time_vec[" << i << "] = " << local_time_vec[i] << std::endl;
    }
  }

  // Get nodes out of trailingInterpBuffer_:

  Array<Scalar> node_vec; 
  status = trailingInterpBuffer_->getNodes(&node_vec); 
  if (!status) {
    return(status);
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    if (node_vec.size() == 0) {
      *out << "trailingInterpBuffer_->getNodes: node_vec = empty vector" << std::endl;
    }
    else {
      *out << "trailingInterpBuffer_->getNodes:" << std::endl;
      for (unsigned int i=0 ; i<node_vec.size() ; ++i) {
        *out << "node_vec[" << i << "] = " << node_vec[i] << std::endl;
      }
    }
  }
  if (node_vec.size() == 0) {
    // Initialization case for an empty InterpolationBuffer
    Array<Scalar> stepper_vec;
    status = stepper_->getNodes(&stepper_vec);
    if (!status) {
      return(status);
    }
    if (stepper_vec.size() < 2) {
      // Stepper and trailingInterpBuffer_ are basically empty
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "Initializing empty Stepper and InterpolationBuffer" << std::endl;
      }
      Scalar step_taken;
      if (takeVariableSteps_) {
        step_taken = stepper_->takeStep(ST::zero(),VARIABLE_STEP);
      } 
      else {
        step_taken = stepper_->takeStep(fixed_dt_,FIXED_STEP);
      }
      // Pass information from stepper_ to trailingInterpBuffer_:
      status = stepper_->getNodes(&stepper_vec);
      if (!status) {
        return(status);
      }
      Scalar stepper_begin = stepper_vec.front();
      Scalar stepper_end = stepper_vec.back();
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "stepper_begin = " << stepper_begin << std::endl;
        *out << "stepper_end = " << stepper_end << std::endl;
      }
      status = trailingInterpBuffer_->setRange(timeRange(stepper_begin,stepper_end),*stepper_);
      if (!status) {
        return(status);
      }
      status = trailingInterpBuffer_->getNodes(&node_vec);
      if (!status) {
        return(status);
      }
    }
    else {
      // Just the trailingInterpBuffer_ is empty
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "Initializing empty InterpolationBuffer" << std::endl;
      }
      Array<Scalar> stepper_vec;
      status = stepper_->getNodes(&stepper_vec);
      if (!status) {
        return(status);
      }
      Scalar stepper_begin = stepper_vec.front();
      Scalar stepper_end = stepper_vec.back();
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "stepper_begin = " << stepper_begin << std::endl;
        *out << "stepper_end = " << stepper_end << std::endl;
      }
      status = trailingInterpBuffer_->setRange(timeRange(stepper_begin,stepper_end),*stepper_);
      if (!status) {
        return(status);
      }
      status = trailingInterpBuffer_->getNodes(&node_vec);
      if (!status) {
        return(status);
      }
    }
  }
  Scalar node_begin = node_vec.front();
  Scalar node_end = node_vec.back();
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "node_begin = " << node_begin << std::endl;
    *out << "node_end = " << node_end << std::endl;
  }

  // Check for valid input of local_time_vec:  (check initialization conditions)
  if ((*(local_time_vec.end()) < node_vec[0]) || (*(local_time_vec.begin()) < node_vec[0])) {
    return(false);
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Requested times are valid." << std::endl;
  }

  // Get time out of stepper_:
  Array<Scalar> stepper_vec;
  status = stepper_->getNodes(&stepper_vec);
  if (!status) {
    return(status);
  }
  if (stepper_vec.size() < 2) {
    // Initialization case for an empty Stepper
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "Initializing empty Stepper" << std::endl;
    }
    Scalar stepper_begin = *(stepper_vec.begin()); // There must be an IC in the stepper_.
    Scalar step_taken;
    if (takeVariableSteps_) {
      step_taken = stepper_->takeStep(ST::zero(),VARIABLE_STEP);
    }
    else {
      step_taken = stepper_->takeStep(fixed_dt_,FIXED_STEP);
    }
    status = stepper_->getNodes(&stepper_vec);
    if (!status) {
      return(status);
    }
    // Pass information from stepper_ to trailingInterpBuffer_:
    status = trailingInterpBuffer_->setRange(timeRange(stepper_begin,stepper_begin+step_taken),*stepper_);
    if (!status) { 
      return(status);
    }
    status = trailingInterpBuffer_->getNodes(&node_vec); 
    if (!status) { 
      return(status);
    }
    node_begin = node_vec[0];
    node_end = node_vec[node_vec.size()-1];
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "node_begin = " << node_begin << std::endl;
      *out << "node_end = " << node_end << std::endl;
    }
  }
  Scalar stepper_begin = stepper_vec[0];
  Scalar stepper_end = stepper_vec[stepper_vec.size()-1];
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "stepper_begin = " << stepper_begin << std::endl;
    *out << "stepper_end = " << stepper_end << std::endl;
  }
  int num = local_time_vec.size();
  for (int i=0; i<num ; ++i) {
    if ( ( node_begin <= local_time_vec[i] ) && ( local_time_vec[i] <= node_end ) ) {
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "Requested local_time_vec[" << i << "] = " << 
          local_time_vec[i] << " available in trailingInterpBuffer_[" << 
          node_begin << "," << node_end << "]" << std::endl;
      }
      Array<Scalar> tmp_time_vec; 
      tmp_time_vec.clear();
      tmp_time_vec.push_back(local_time_vec[i]);
      Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > tmp_x_vec, tmp_xdot_vec;
      Array<ScalarMag> tmp_accuracy_vec;
      status = trailingInterpBuffer_->getPoints(tmp_time_vec, &tmp_x_vec, &tmp_xdot_vec, &tmp_accuracy_vec); 
      x_vec.push_back(tmp_x_vec[0]);
      xdot_vec.push_back(tmp_xdot_vec[0]);
      accuracy_vec.push_back(tmp_accuracy_vec[0]);
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "trailingInterpBuffer_->getPoints returned:" << std::endl;
        *out << "tmp_x_vec = " << std::endl;
        tmp_x_vec[0]->describe(*out,Teuchos::VERB_EXTREME);
        if (tmp_xdot_vec[0] == Teuchos::null) {
          *out << "tmp_xdot_vec = Teuchos::null" << std::endl;
        } else {
          *out << "tmp_xdot_vec = " << std::endl;
          tmp_xdot_vec[0]->describe(*out,Teuchos::VERB_EXTREME);
        }
        *out << "tmp_accuracy_vec = " << accuracy_vec[0] << std::endl;
      }
      if (!status) { 
        return(status);
      }
    }
    else if ( ( stepper_begin <= local_time_vec[i] ) && ( local_time_vec[i] <= stepper_end ) ) {
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "Requested local_time_vec[" << i << "] = " << 
          local_time_vec[i] << " available in stepper_[" << 
          stepper_begin << "," << stepper_end << "]" << std::endl;
      }
      Array<Scalar> tmp_time_vec; 
      tmp_time_vec.push_back(local_time_vec[i]);
      Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > tmp_x_vec, tmp_xdot_vec;
      Array<ScalarMag> tmp_accuracy_vec;
      status = stepper_->getPoints(tmp_time_vec, &tmp_x_vec, &tmp_xdot_vec, &tmp_accuracy_vec); 
      x_vec.push_back(tmp_x_vec[0]);
      xdot_vec.push_back(tmp_xdot_vec[0]);
      accuracy_vec.push_back(tmp_accuracy_vec[0]);
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "stepper_->getPoints returned:" << std::endl;
        *out << "tmp_x_vec = " << std::endl;
        tmp_x_vec[0]->describe(*out,Teuchos::VERB_EXTREME);
        if (tmp_xdot_vec[0] == Teuchos::null) {
          *out << "tmp_xdot_vec = Teuchos::null" << std::endl;
        } else {
          *out << "tmp_xdot_vec = " << std::endl;
          tmp_xdot_vec[0]->describe(*out,Teuchos::VERB_EXTREME);
        }
        *out << "tmp_accuracy_vec = " << tmp_accuracy_vec[0] << std::endl;
      }
      if (!status) { 
        return(status);
      }
    }
    else {
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "Requested local_time_vec[" << i << "] = " <<
          local_time_vec[i] << " not available in trailingInterpBuffer_ or stepper_" << std::endl;
        *out << "Integrating forward with stepper_" << std::endl;
      }
      int num_local_steps_taken = 0;
      while (stepper_end < local_time_vec[i]) {
        // integrate forward with stepper_ 
        Scalar step_taken;
        if (takeVariableSteps_) {
          step_taken = stepper_->takeStep(ST::zero(),VARIABLE_STEP);
        }
        else {
          step_taken = stepper_->takeStep(fixed_dt_,FIXED_STEP);
        }
        num_local_steps_taken++;
        if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) ) {
          *out << "Took step of size " << step_taken << " with stepper_" << std::endl;
        }
        // Pass information from stepper_ to trailingInterpBuffer_:
        status = trailingInterpBuffer_->setRange(timeRange(stepper_end+step_taken,stepper_end+step_taken),*stepper_);
        if (!status) { 
          return(status);
        }
        // Check to see if we're past the current requested local_time_vec[i] point:
        if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
          *out << "Checking to see if we're past the current requested time" << std::endl;
          *out << "local_time_vec[" << i << "] = " << local_time_vec[i] 
               << " <= " << stepper_end+step_taken << " stepper_end+step_taken" << std::endl;
        }
        if (local_time_vec[i] <= stepper_end+step_taken) {
          if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
            *out << "We've integrated past local_time_vec[" << i << "] = " 
                 << local_time_vec[i] << "!" << std::endl;
          }
          Array<Scalar> tmp_time_vec;
          tmp_time_vec.push_back(local_time_vec[i]);
          Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > tmp_x_vec, tmp_xdot_vec;
          Array<ScalarMag> tmp_accuracy_vec;
          status = stepper_->getPoints(tmp_time_vec, &tmp_x_vec, &tmp_xdot_vec, &tmp_accuracy_vec); 
          x_vec.push_back(tmp_x_vec[0]);
          xdot_vec.push_back(tmp_xdot_vec[0]);
          accuracy_vec.push_back(tmp_accuracy_vec[0]);
          if (!status) { 
            return(status);
          }
        }
        // Update end-points:
        node_end += step_taken;
        stepper_begin = stepper_end;
        stepper_end += step_taken;
      }
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
        *out << "Took " << num_local_steps_taken << " steps with " << stepper_->description() << std::endl;
      }
    }
  }
  return(status);
}

// 2007/06/08: rabartl: ToDo: The above function getPoints(...) is way to long
// and too complicated.  I have not looked at it carefully enough but there
// must be some opportunities to resuse some intermediate code and to make
// this look cleaner.  In my opinion, this is how this routine should be
// written:
//
// 1) Sort the time points (or require the user to sort the time points and
// then check for this).
//
// 2) If any time points fall before the storage given in
// trailingInterpBuffer, then thrown an exception.
//
// 3) Any points that fall in the time range of the trailingInterpBuffer are
// first gotten.
//
// 4) An time points that fall in the current range of the stepper are gotten
//
// 5) Any points that fall before the current storage of the stepper are
// gotten, from lowest to largest, as the stepper is advanced.
//
// Given the above algorithm, we can satisfy any time_vec request when the
// stepper is just starting out.


template<class Scalar>
bool IntegratorDefault<Scalar>::setRange(
  const TimeRange<Scalar>& range,
  const InterpolationBufferBase<Scalar> & interpBuffer
  )
{
  using Teuchos::as;
  const Scalar time_lower = range.lower();
  const Scalar time_upper = range.upper();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"IBAS::setRange");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "time_lower = " << time_lower << std::endl;
    *out << "time_upper = " << time_upper << std::endl;
    *out << "interpBuffer_ = " << interpBuffer.description() << std::endl;
  }
  return(trailingInterpBuffer_->setRange(timeRange(time_lower,time_upper),interpBuffer));
}


template<class Scalar>
TimeRange<Scalar> IntegratorDefault<Scalar>::getTimeRange() const
{
  TEST_FOR_EXCEPT("ToDo: Implement this!");
  return invalidTimeRange<Scalar>(); // ToDo: Fill in this range, I know you have one!
}


template<class Scalar>
bool IntegratorDefault<Scalar>::getNodes(
  Array<Scalar>* time_vec
  ) const
{
  using Teuchos::as;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"IBAS::getNodes");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << this->description() << std::endl;
  }
  return(trailingInterpBuffer_->getNodes(time_vec));
}


template<class Scalar>
bool IntegratorDefault<Scalar>::removeNodes(
  Array<Scalar>& time_vec
  ) 
{
  return(trailingInterpBuffer_->removeNodes(time_vec));
}


template<class Scalar>
int IntegratorDefault<Scalar>::getOrder() const
{
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
    if (!is_null(trailingInterpBuffer_))
      out << describe(*trailingInterpBuffer_,verbLevel);
    else
      out << "NULL\n";
    out << "stepper = ";
    if (!is_null(stepper_))
      out << describe(*stepper_,verbLevel);
    else
      out << "NULL\n";
  }
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM)) {
    out << "paramList_ = " << paramList_->print(out);
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
  Teuchos::OSTab ostab(out,1,"IBAS::setParameterList");

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "paramList = " << paramList->print(*out) << std::endl;
  }

  paramList->validateParametersAndSetDefaults(*getValidParameters());
  paramList_ = paramList;

  takeVariableSteps_ = get<bool>(*paramList_,"Take Variable Steps");
  fixed_dt_ = get<double>(*paramList_,"fixed_dt");

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
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return (validPL);
}

template <class Scalar>
bool IntegratorDefault<Scalar>::getFwdPoints(
    const Array<Scalar>& time_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    )
{
  return(false);
}

template <class Scalar>
TimeRange<Scalar> IntegratorDefault<Scalar>::getFwdTimeRange() const
{
  TEST_FOR_EXCEPT("ToDo: Implement this!");
  return invalidTimeRange<Scalar>(); // ToDo: Fill in this range, I know you have one!
}

} // namespace Rythmos


#endif //Rythmos_INTEGRATOR_DEFAULT_H
