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

#ifndef Rythmos_INTERPOLATION_BUFFER_AS_STEPPER_H
#define Rythmos_INTERPOLATION_BUFFER_AS_STEPPER_H

#include "Rythmos_InterpolationBufferBase.hpp"
#include "Rythmos_StepperBase.hpp"

namespace Rythmos {

/** \brief Base class for defining interpolation buffer functionality. */
template<class Scalar> 
class InterpolationBufferAsStepper : virtual public Rythmos::InterpolationBufferBase<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /// Destructor
    ~InterpolationBufferAsStepper() {};

    /// Constructors
    InterpolationBufferAsStepper();
    InterpolationBufferAsStepper(
      const Teuchos::RefCountPtr<Rythmos::StepperBase<Scalar> > &stepper_
      ,const Teuchos::RefCountPtr<Rythmos::InterpolationBufferBase<Scalar> > &IB_
      ,const Teuchos::RefCountPtr<Teuchos::ParameterList> &parameterList_ = Teuchos::null
      );

    /// Set InterpolationBufferBase:
    void setInterpolationBuffer(const Teuchos::RefCountPtr<Rythmos::InterpolationBufferBase<Scalar> > &IB_);

    /// Set Stepper:
    void setStepper(const Teuchos::RefCountPtr<Rythmos::StepperBase<Scalar> > &stepper_);
    
    /// Redefined from InterpolationBufferBase
    /// This is a pass-through to the underlying InterpolationBufferBase:
    bool SetPoints(
      const std::vector<Scalar>& time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
      ,const std::vector<ScalarMag> & accuracy_vec 
      );

    // This is not a pass-through.
    bool GetPoints(
      const std::vector<Scalar>& time_vec_
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec_
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec_
      ,std::vector<ScalarMag>* accuracy_vec_) const;

    /// This is a pass-through to the underlying InterpolationBufferBase:
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar>& IB_);

    /// This is a pass-through to the underlying InterpolationBufferBase:
    bool GetNodes(std::vector<Scalar>* time_vec) const;

    /// This is a pass-through to the underlying InterpolationBufferBase:
    virtual bool RemoveNodes(std::vector<Scalar>& time_vec);

    /// This is a pass-through to the underlying InterpolationBufferBase:
    int GetOrder() const;

    /// Redefined from Teuchos::Describable
    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream       &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();

    /** \brief . */
    Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  private:

    // Interpolation Buffer used to store past results
    Teuchos::RefCountPtr<Rythmos::InterpolationBufferBase<Scalar> > IB;

    // Stepper used to fill interpolation buffer.
    Teuchos::RefCountPtr<Rythmos::StepperBase<Scalar> > stepper;

    // ParameterList to control behavior
    Teuchos::RefCountPtr<Teuchos::ParameterList> parameterList;

};

// ////////////////////////////
// Defintions
template<class Scalar>
InterpolationBufferAsStepper<Scalar>::InterpolationBufferAsStepper(
    const Teuchos::RefCountPtr<Rythmos::StepperBase<Scalar> > &stepper_
    ,const Teuchos::RefCountPtr<Rythmos::InterpolationBufferBase<Scalar> > &IB_
    ,const Teuchos::RefCountPtr<Teuchos::ParameterList> &parameterList_ 
    )
{
  setParameterList(parameterList_);
  int outputLevel = parameterList_->get( "outputLevel", int(-1) );
  outputLevel = min(max(outputLevel,-1),4);
  this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(20);
  out->setMaxLenLinePrefix(40);
  out->pushLinePrefix("Rythmos::InterpolationBufferAsStepper");
  out->setShowLinePrefix(true);
  out->setTabIndentStr("    ");
  *out << "Initializing InterpolationBufferAsStepper" << std::endl;
  if (outputLevel >= 3)
    *out << "Calling setStepper..." << std::endl;
  setStepper(stepper_);
  if (outputLevel >= 3)
    *out << "Calling setInterpolationBuffer..." << std::endl;
  setInterpolationBuffer(IB_);
}

template<class Scalar>
void InterpolationBufferAsStepper<Scalar>::setStepper(
    const Teuchos::RefCountPtr<Rythmos::StepperBase<Scalar> > &stepper_
    )
{
  // 10/9/06 tscoffe:  What should we do if this is called after initialization?
  //                   Basically, you're swapping out the stepper for a new one.
  //                   Since we're copying the data into IB after each
  //                   stepper->TakeStep() call, this should be fine, and it
  //                   will essentially result in changing the stepper
  //                   mid-stream.  If the new stepper has a time value before
  //                   the end of the data in IB, then you will not get new
  //                   stepper data until you ask for time values after the end
  //                   of IB's data.  And then the stepper will walk forward
  //                   inserting new (potentially inconsistent) data into IB
  //                   until it can give you the time values you asked for.
  //                   Then IB will potentially have old and new data in it.
  //                   On the other hand, if you swap out the stepper and the
  //                   time value is synchronized with the old stepper, then
  //                   you will essentially change the integrator mid-stream
  //                   and everything should proceed without problems.
  //                   Note also:  this functionality is important for checkpointing.
  stepper = stepper_;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IBAS::setStepper");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "stepper = " << stepper->description() << std::endl;
}

template<class Scalar>
void InterpolationBufferAsStepper<Scalar>::setInterpolationBuffer(
    const Teuchos::RefCountPtr<Rythmos::InterpolationBufferBase<Scalar> > &IB_
    )
{
  // 10/9/06 tscoffe:  What should we do if this is called after initialization?
  //                   Basically, you're swapping out the history for a new
  //                   one.  This could be the result of upgrading or
  //                   downgrading the accuracy of the buffer, or something I
  //                   haven't thought of yet.  Since IB's node_vec is checked
  //                   each time GetPoints is called, this should be fine.  And
  //                   the time values in IB need not be synchronized with
  //                   stepper.
  //                   Note also:  this functionality is important for checkpointing.
  IB = IB_;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IBAS::setInterpolationBuffer");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "IB = " << IB->description() << std::endl;
}


template<class Scalar>
bool InterpolationBufferAsStepper<Scalar>::SetPoints(
      const std::vector<Scalar>& time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
      ,const std::vector<ScalarMag> & accuracy_vec 
      ) 
{
  return(IB->SetPoints(time_vec,x_vec,xdot_vec,accuracy_vec));
}

template<class Scalar>
bool InterpolationBufferAsStepper<Scalar>::GetPoints(
      const std::vector<Scalar>& time_vec_
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec_ptr_
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec_ptr_
      ,std::vector<ScalarMag>* accuracy_vec_ptr_
      ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IBAS::GetPoints");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "time_vec_ = " << std::endl;
    for (unsigned int i=0 ; i<time_vec_.size() ; ++i)
      *out << "time_vec_[" << i << "] = " << time_vec_[i] << std::endl;
    if (x_vec_ptr_ == NULL)
      *out << "x_vec_ptr_ = NULL" << std::endl;
    else if (x_vec_ptr_->size() == 0)
      *out << "x_vec_ptr_ = ptr to empty vector" << std::endl;
    else
    {
      *out << "x_vec_ptr_ = " << std::endl;
      for (unsigned int i=0 ; i<x_vec_ptr_->size() ; ++i)
      {
        *out << "x_vec[" << i << "] = " << std::endl;
        (*x_vec_ptr_)[i]->describe(*out,Teuchos::VERB_EXTREME);
      }
    }
    if (xdot_vec_ptr_ == NULL)
      *out << "xdot_vec_ptr_ = NULL" << std::endl;
    else if (xdot_vec_ptr_->size() == 0)
      *out << "xdot_vec_ptr_ = ptr to empty vector" << std::endl;
    else
    {
      *out << "xdot_vec = " << std::endl;
      for (unsigned int i=0 ; i<xdot_vec_ptr_->size() ; ++i)
      {
        if ((*xdot_vec_ptr_)[i] == Teuchos::null)
          *out << "xdot_vec[" << i << "] = Teuchos::null" << std::endl;
        else
        {
          *out << "xdot_vec[" << i << "] = " << std::endl;
          (*xdot_vec_ptr_)[i]->describe(*out,Teuchos::VERB_EXTREME);
        }
      }
    }
    if (accuracy_vec_ptr_ == NULL)
      *out << "accuracy_vec_ptr_ = NULL" << std::endl;
    else if (accuracy_vec_ptr_->size() == 0)
      *out << "accuracy_vec_ptr_ = ptr to empty vector" << std::endl;
    else
    {
      *out << "accuracy_vec = " << std::endl;
      for (unsigned int i=0 ; i<accuracy_vec_ptr_->size() ; ++i)
        *out << "accuracy_vec[" << i << "] = " << (*accuracy_vec_ptr_)[i] << std::endl;
    }
  }
  bool status = IB->GetPoints(time_vec_,x_vec_ptr_,xdot_vec_ptr_,accuracy_vec_ptr_);
  if (status) return(status);
  x_vec_ptr_->clear();
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "IB->GetPoints unsuccessful" << std::endl;
  status = true;
  std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec = *x_vec_ptr_;
  std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec = *xdot_vec_ptr_;
  std::vector<ScalarMag> &accuracy_vec = *accuracy_vec_ptr_;
  // Sort time_vec_
  std::vector<Scalar> local_time_vec = time_vec_;
  std::sort(local_time_vec.begin(),local_time_vec.end());
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "sorted local time_vec:" << std::endl;
    for (unsigned int i=0 ; i<local_time_vec.size() ; ++i)
      *out << "local_time_vec[" << i << "] = " << local_time_vec[i] << std::endl;
  }
  // Get nodes out of IB:
  std::vector<Scalar> node_vec; 
  status = IB->GetNodes(&node_vec); 
  if (!status) return(status);
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    if (node_vec.size() == 0)
      *out << "IB->GetNodes: node_vec = empty vector" << std::endl;
    else
    {
      *out << "IB->GetNodes:" << std::endl;
      for (unsigned int i=0 ; i<node_vec.size() ; ++i)
        *out << "node_vec[" << i << "] = " << node_vec[i] << std::endl;
    }
  }
  if (node_vec.size() == 0)
  {
    // Initialization case for an empty InterpolationBuffer
    std::vector<Scalar> stepper_vec;
    status = stepper->GetNodes(&stepper_vec);
    if (!status) return(status);
    if (stepper_vec.size() < 2)
    {
      // Stepper and IB are basically empty
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
        *out << "Initializing empty Stepper and InterpolationBuffer" << std::endl;
      Scalar step_taken;
      if (parameterList->isParameter("fixed_dt"))
        step_taken = stepper->TakeStep(parameterList->get<Scalar>("fixed_dt"),FIXED_STEP);
      else
        step_taken = stepper->TakeStep(ST::zero(),VARIABLE_STEP);
      // Pass information from stepper to IB:
      status = stepper->GetNodes(&stepper_vec);
      if (!status) return(status);
      Scalar stepper_begin = stepper_vec.front();
      Scalar stepper_end = stepper_vec.back();
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "stepper_begin = " << stepper_begin << std::endl;
        *out << "stepper_end = " << stepper_end << std::endl;
      }
      status = IB->SetRange(stepper_begin,stepper_end,*stepper);
      if (!status) return(status);
      status = IB->GetNodes(&node_vec);
      if (!status) return(status);
    }
    else 
    {
      // Just the IB is empty
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
        *out << "Initializing empty InterpolationBuffer" << std::endl;
      std::vector<Scalar> stepper_vec;
      status = stepper->GetNodes(&stepper_vec);
      if (!status) return(status);
      Scalar stepper_begin = stepper_vec.front();
      Scalar stepper_end = stepper_vec.back();
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "stepper_begin = " << stepper_begin << std::endl;
        *out << "stepper_end = " << stepper_end << std::endl;
      }
      status = IB->SetRange(stepper_begin,stepper_end,*stepper);
      if (!status) return(status);
      status = IB->GetNodes(&node_vec);
      if (!status) return(status);
    }
  }
  Scalar node_begin = node_vec.front();
  Scalar node_end = node_vec.back();
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "node_begin = " << node_begin << std::endl;
    *out << "node_end = " << node_end << std::endl;
  }
  // Check for valid input of local_time_vec:  (check initialization conditions)
  if ((*(local_time_vec.end()) < node_vec[0]) || (*(local_time_vec.begin()) < node_vec[0]))
    return(false);
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "Requested times are valid." << std::endl;
  // Get time out of stepper:
  std::vector<Scalar> stepper_vec;
  status = stepper->GetNodes(&stepper_vec);
  if (!status) return(status);
  if (stepper_vec.size() < 2)
  {
    // Initialization case for an empty Stepper
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      *out << "Initializing empty Stepper" << std::endl;
    Scalar stepper_begin = *(stepper_vec.begin()); // There must be an IC in the stepper.
    Scalar step_taken;
    if (parameterList->isParameter("fixed_dt"))
      step_taken = stepper->TakeStep(parameterList->get<Scalar>("fixed_dt"),FIXED_STEP);
    else
      step_taken = stepper->TakeStep(ST::zero(),VARIABLE_STEP);
    status = stepper->GetNodes(&stepper_vec);
    if (!status) return(status);
    // Pass information from stepper to IB:
    status = IB->SetRange(stepper_begin,stepper_begin+step_taken,*stepper);
    if (!status) return(status);
    status = IB->GetNodes(&node_vec); 
    if (!status) return(status);
    node_begin = node_vec[0];
    node_end = node_vec[node_vec.size()-1];
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    {
      *out << "node_begin = " << node_begin << std::endl;
      *out << "node_end = " << node_end << std::endl;
    }
  }
  Scalar stepper_begin = stepper_vec[0];
  Scalar stepper_end = stepper_vec[stepper_vec.size()-1];
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "stepper_begin = " << stepper_begin << std::endl;
    *out << "stepper_end = " << stepper_end << std::endl;
  }
  int num = local_time_vec.size();
  for (int i=0; i<num ; ++i)
  {
    if ( ( node_begin <= local_time_vec[i] ) && ( local_time_vec[i] <= node_end ) )
    {
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "Requested local_time_vec[" << i << "] = " << 
          local_time_vec[i] << " available in IB[" << 
          node_begin << "," << node_end << "]" << std::endl;
      }
      std::vector<Scalar> tmp_time_vec; 
      tmp_time_vec.clear();
      tmp_time_vec.push_back(local_time_vec[i]);
      std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > tmp_x_vec, tmp_xdot_vec;
      std::vector<ScalarMag> tmp_accuracy_vec;
      status = IB->GetPoints(tmp_time_vec, &tmp_x_vec, &tmp_xdot_vec, &tmp_accuracy_vec); 
      x_vec.push_back(tmp_x_vec[0]);
      xdot_vec.push_back(tmp_xdot_vec[0]);
      accuracy_vec.push_back(tmp_accuracy_vec[0]);
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "IB->GetPoints returned:" << std::endl;
        *out << "tmp_x_vec = " << std::endl;
        tmp_x_vec[0]->describe(*out,Teuchos::VERB_EXTREME);
        if (tmp_xdot_vec[0] == Teuchos::null)
          *out << "tmp_xdot_vec = Teuchos::null" << std::endl;
        else
        {
          *out << "tmp_xdot_vec = " << std::endl;
          tmp_xdot_vec[0]->describe(*out,Teuchos::VERB_EXTREME);
        }
        *out << "tmp_accuracy_vec = " << accuracy_vec[0] << std::endl;
      }
      if (!status) return(status);
    }
    else if ( ( stepper_begin <= local_time_vec[i] ) && ( local_time_vec[i] <= stepper_end ) )
    {
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "Requested local_time_vec[" << i << "] = " << 
          local_time_vec[i] << " available in stepper[" << 
          stepper_begin << "," << stepper_end << "]" << std::endl;
      }
      std::vector<Scalar> tmp_time_vec; 
      tmp_time_vec.push_back(local_time_vec[i]);
      std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > tmp_x_vec, tmp_xdot_vec;
      std::vector<ScalarMag> tmp_accuracy_vec;
      status = stepper->GetPoints(tmp_time_vec, &tmp_x_vec, &tmp_xdot_vec, &tmp_accuracy_vec); 
      x_vec.push_back(tmp_x_vec[0]);
      xdot_vec.push_back(tmp_xdot_vec[0]);
      accuracy_vec.push_back(tmp_accuracy_vec[0]);
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "stepper->GetPoints returned:" << std::endl;
        *out << "tmp_x_vec = " << std::endl;
        tmp_x_vec[0]->describe(*out,Teuchos::VERB_EXTREME);
        if (tmp_xdot_vec[0] == Teuchos::null)
          *out << "tmp_xdot_vec = Teuchos::null" << std::endl;
        else
        {
          *out << "tmp_xdot_vec = " << std::endl;
          tmp_xdot_vec[0]->describe(*out,Teuchos::VERB_EXTREME);
        }
        *out << "tmp_accuracy_vec = " << tmp_accuracy_vec[0] << std::endl;
      }
      if (!status) return(status);
    }
    else
    {
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "Requested local_time_vec[" << i << "] = " <<
          local_time_vec[i] << " not available in IB or stepper" << std::endl;
        *out << "Integrating forward with stepper" << std::endl;
      }
      int num_local_steps_taken = 0;
      while (stepper_end < local_time_vec[i])
      {
        // integrate forward with stepper 
        Scalar step_taken;
        if (parameterList->isParameter("fixed_dt"))
          step_taken = stepper->TakeStep(parameterList->get<Scalar>("fixed_dt"),FIXED_STEP);
        else
          step_taken = stepper->TakeStep(ST::zero(),VARIABLE_STEP);
        num_local_steps_taken++;
        if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_MEDIUM) )
          *out << "Took step of size " << step_taken << " with stepper" << std::endl;
        // Pass information from stepper to IB:
        status = IB->SetRange(stepper_end+step_taken,stepper_end+step_taken,*stepper);
        if (!status) return(status);
        // Check to see if we're past the current requested local_time_vec[i] point:
        if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
        {
          *out << "Checking to see if we're past the current requested time" << std::endl;
          *out << "local_time_vec[" << i << "] = " << local_time_vec[i] 
            << " <= " << stepper_end+step_taken << " stepper_end+step_taken" << std::endl;
        }
        if (local_time_vec[i] <= stepper_end+step_taken)
        {
          if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
          {
            *out << "We've integrated past local_time_vec[" << i << "] = " 
              << local_time_vec[i] << "!" << std::endl;
          }
          std::vector<Scalar> tmp_time_vec;
          tmp_time_vec.push_back(local_time_vec[i]);
          std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > tmp_x_vec, tmp_xdot_vec;
          std::vector<ScalarMag> tmp_accuracy_vec;
          status = stepper->GetPoints(tmp_time_vec, &tmp_x_vec, &tmp_xdot_vec, &tmp_accuracy_vec); 
          x_vec.push_back(tmp_x_vec[0]);
          xdot_vec.push_back(tmp_xdot_vec[0]);
          accuracy_vec.push_back(tmp_accuracy_vec[0]);
          if (!status) return(status);
        }
        // Update end-points:
        node_end += step_taken;
        stepper_begin = stepper_end;
        stepper_end += step_taken;
      }
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_LOW) )
        *out << "Took " << num_local_steps_taken << " steps with " << stepper->description() << std::endl;
    }
  }
  return(status);
}

template<class Scalar>
bool InterpolationBufferAsStepper<Scalar>::SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar> & IB_
      )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IBAS::SetRange");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "time_lower = " << time_lower << std::endl;
    *out << "time_upper = " << time_upper << std::endl;
    *out << "IB = " << IB_.description() << std::endl;
  }
  return(IB->SetRange(time_lower,time_upper,IB_));
}

template<class Scalar>
bool InterpolationBufferAsStepper<Scalar>::GetNodes(
    std::vector<Scalar>* time_vec
    ) const
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IBAS::GetNodes");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << this->description() << std::endl;
  return(IB->GetNodes(time_vec));
}

template<class Scalar>
bool InterpolationBufferAsStepper<Scalar>::RemoveNodes(
    std::vector<Scalar>& time_vec
    ) 
{
  return(IB->RemoveNodes(time_vec));
}

template<class Scalar>
int InterpolationBufferAsStepper<Scalar>::GetOrder() const
{
  return(IB->GetOrder());
}

template<class Scalar>
std::string InterpolationBufferAsStepper<Scalar>::description() const
{
  std::string name = "Rythmos::InterpolationBufferAsStepper";
  return(name);
}

template<class Scalar>
void InterpolationBufferAsStepper<Scalar>::describe(
      Teuchos::FancyOStream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const
{
  if ( (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_DEFAULT) ) ||
       (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)     )
     )
  {
    out << description() << "::describe" << std::endl;
    out << "interpolation buffer = " << IB->description() << std::endl;
    out << "stepper = " << stepper->description() << std::endl;
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM))
  {
    out << "parameterList = " << parameterList->print(out) << std::endl;
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH))
  {
  }
}

template <class Scalar>
void InterpolationBufferAsStepper<Scalar>::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& parameterList_)
{
  if (parameterList_ == Teuchos::null)
    parameterList = Teuchos::rcp(new Teuchos::ParameterList);
  else
    parameterList = parameterList_;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IBAS::setParameterList");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "parameterList = " << parameterList->print(*out) << std::endl;
  }
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> InterpolationBufferAsStepper<Scalar>::getParameterList()
{
  return(parameterList);
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> InterpolationBufferAsStepper<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> temp_param_list = parameterList;
  parameterList = Teuchos::null;
  return(temp_param_list);
}

template <class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList> InterpolationBufferAsStepper<Scalar>::getValidParameters() const
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> temp_param_list = Teuchos::rcp(new Teuchos::ParameterList);
  temp_param_list->set<Scalar>( "fixed_dt", Scalar(0.1) );
  return(temp_param_list);
}

} // namespace Rythmos

#endif //Rythmos_INTERPOLATION_BUFFER_AS_STEPPER_H
