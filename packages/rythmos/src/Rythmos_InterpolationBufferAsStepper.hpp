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

#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_Stepper.hpp"

namespace Rythmos {

/** \brief Base class for defining interpolation buffer functionality. */
template<class Scalar> 
class InterpolationBufferAsStepper : virtual public Rythmos::InterpolationBuffer<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /// Destructor
    ~InterpolationBufferAsStepper() {};

    /// Constructor
    InterpolationBufferAsStepper();

    /// Add InterpolationBufferAsStepper:
    setInterpolationBuffer(const Teuchos::RefCountPtr<const Rythmos::InterpolationBuffer<Scalar> > &IB_);

    /// Add Stepper:
    setStepper(const Teuchos::RefCountPtr<const Rythmos::Stepper<Scalar> > &stepper_);

    /// Redefined from InterpolationBuffer
    /// This is a pass-through to the underlying InterpolationBuffer:
    bool SetPoints(
      const std::vector<Scalar>& time_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& x_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& xdot_list);

    // This is not a pass-through.
    bool GetPoints(
      const std::vector<Scalar>& time_list_
      ,std::vector<Thyra::VectorBase<Scalar> >* x_list_
      ,std::vector<Thyra::VectorBase<Scalar> >* xdot_list_
      ,std::vector<ScalarMag>* accuracy_list_) const;

    /// This is a pass-through to the underlying InterpolationBuffer:
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBuffer<Scalar> & IB_);

    /// This is a pass-through to the underlying InterpolationBuffer:
    bool GetNodes(std::vector<Scalar>* time_list) const;

    /// This is a pass-through to the underlying InterpolationBuffer:
    int GetOrder() const;

  private:

    // Interpolation Buffer used to store past results
    Teuchos::RefCountPtr<Rythmos::InterpolationBuffer<Scalar> > IB;

    // Stepper used to fill interpolation buffer.
    Teuchos::RefCountPtr<Rythmos::Stepper<Scalar> > stepper;

};

// ////////////////////////////
// Defintions
template<class Scalar>
InterpolationBufferAsStepper<Scalar>::InterpolationBufferAsStepper(
    const Teuchos::RefCountPtr<const Rythmos::Stepper<Scalar> > &stepper_
    ,const Teuchos::RefCountPtr<const Rythmos::InterpolationBuffer<Scalar> > &IB_
    )
{
  setStepper(stepper_);
  setInterpolationBuffer(IB_);
}

template<class Scalar>
InterpolationBufferAsStepper<Scalar>::setStepper(
    const Teuchos::RefCountPtr<const Rythmos::Stepper<Scalar> > &stepper_
    )
{
  stepper = stepper_;
}

template<class Scalar>
InterpolationBufferAsStepper<Scalar>::setInterpolationBuffer(
    const Teuchos::RefCountPtr<const Rythmos::InterpolationBuffer<Scalar> > &IB_
    )
{
  IB = IB_;
}

template<class Scalar>
bool InterpolationBufferAsStepper<Scalar>::SetPoints(
      const std::vector<Scalar>& time_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& x_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& xdot_list
      ) 
{
  return(IB->SetPoints(time_list,x_list,xdot_list));
}

template<class Scalar>
bool InterpolationBufferAsStepper<Scalar>::GetPoints(
      const std::vector<Scalar>& time_list_
      ,std::vector<Thyra::VectorBase<Scalar> >* x_list_
      ,std::vector<Thyra::VectorBase<Scalar> >* xdot_list_
      ,std::vector<ScalarMag>* accuracy_list_
      ) const
{
  bool status = IB->GetPoints(time_list_,x_list,xdot_list_,accuracy_list_);
  if (status) return(status);
  status = true;
  // Sort time_list_
  std::sort(time_list_.begin(),time_list_.end());
  // Check for valid input of time_list_:  (check initialization conditions)
  if (time_list_.end() < node_begin)
    return(false);
  // Get nodes out of IB:
  std::vector<Scalar> node_list; 
  status = IB.GetNodes(&node_list); 
  if (!status) return(status);
  Scalar node_begin = node_list.begin();
  Scalar node_end = node_list.end();
  // Get time out of stepper:
  std::vector<Scalar> stepper_list;
  stepper.GetNodes(&stepper_list);
  Scalar stepper_begin = stepper_list.begin();
  Scalar stepper_end = stepper_list.end();
  int num = time_list_->size();
  for (int i=0; i<num ; ++i)
  {
    if ( ( node_begin < time_list_[i] ) && ( time_list_[i] < node_end ) )
    {
      std::vector<Scalar> tmp_time_list = time_list_[i];
      std::vector<Scalar> tmp_x_list, tmp_xdot_list, tmp_accuracy_list;
      status = IB->GetPoints(tmp_time_list, &tmp_x_list, &tmp_xdot_list, &tmp_accuracy_list); 
      x_list_[i] = tmp_x_list[0];
      xdot_list_[i] = tmp_xdot_list[0];
      accuracy_list_[i] = tmp_accuracy_list[0];
      if (!status) return(status);
    }
    else if ( ( stepper_begin < time_list_[i] ) && ( time_list_[i] < stepper_end ) )
    {
      std::vector<Scalar> tmp_time_list = time_list_[i];
      std::vector<Scalar> tmp_x_list, tmp_xdot_list, tmp_accuracy_list;
      status = stepper->GetPoints(tmp_time_list, &tmp_x_list, &tmp_xdot_list, &tmp_accuracy_list); 
      x_list_[i] = tmp_x_list[0];
      xdot_list_[i] = tmp_xdot_list[0];
      accuracy_list_[i] = tmp_accuracy_list[0];
      if (!status) return(status);
    }
    else
    {
      while (stepper_end < time_list_[i])
      {
        // integrate forward with stepper 
        Scalar step_taken = stepper.TakeStep();
        // Pass information from stepper to IB:
        status = IB.SetRange(stepper_end,stepper_end+step_taken,stepper);
        if (!status) return(status);
        // Check to see if we're past the current requested time_list_[i] point:
        if (time_list_[i] <= stepper_end+step_taken)
        {
          std::vector<Scalar> tmp_time_list = time_list_[i];
          std::vector<Scalar> tmp_x_list, tmp_xdot_list, tmp_accuracy_list;
          status = stepper->GetPoints(tmp_time_list, &tmp_x_list, &tmp_xdot_list, &tmp_accuracy_list); 
          x_list_[i] = tmp_x_list[0];
          xdot_list_[i] = tmp_xdot_list[0];
          accuracy_list_[i] = tmp_accuracy_list[0];
          if (!status) return(status);
        }
        // Update end-points:
        node_end += step_taken;
        stepper_begin = stepper_end;
        stepper_end += step_taken;
      }
    }
  }
  return(status);
}

template<class Scalar>
bool InterpolationBufferAsStepper<Scalar>::SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBuffer<Scalar> & IB_
      )
{
  return(IB->SetRange(time_lower,time_upper,IB_));
}

template<class Scalar>
bool InterpolationBufferAsStepper<Scalar>::GetNodes(
    std::vector<Scalar>* time_list
    ) const
{
  return(IB->GetNodes(time_list));
}

template<class Scalar>
int InterpolationBufferAsStepper<Scalar>::GetOrder() const
{
  return(IB->GetOrder());
}


} // namespace Rythmos

#endif //Rythmos_INTERPOLATION_BUFFER_AS_STEPPER_H
