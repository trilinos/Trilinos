//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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

#ifndef Rythmos_LINEAR_INTERPOLATOR_H
#define Rythmos_LINEAR_INTERPOLATOR_H

#include "Rythmos_Interpolator.hpp"

namespace Rythmos {

template<class Scalar>
class LinearInterpolator : virtual public Interpolator<Scalar>
{
  public:

    /// Destructor
    ~LinearInterpolator() {};

    /// Interpolation:
    bool interpolate(
        const std::vector<DataStore<Scalar> > &data_in
        ,const std::vector<Scalar> &t_values
        ,std::vector<DataStore<Scalar> > *data_out
        ) const;

  private:

    Thyra::VectorBase<Scalar> tmp_vec;

};

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::GetPoints(
    const std::vector<Scalar>& time_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
    ,std::vector<Tuechos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
    ,std::vector<ScalarMag>* accuracy_vec) const
{
  // Copy the const time_vec to a local sorted time_vec
  std::vector<Scalar> local_time_vec = time_vec;
  std::sort(local_time_vec.begin(),local_time_vec.end());
  // If there are fewer than 2 points in node_list, then return failure
  if (node_list.size() < 2)
    return(false);
  // If time is outside range of t_values, then return failure
  if ( (*local_time_vec.begin() < node_list.begin().t) || (*local_time_vec.end() > node_list.end().t) )
    return(false);
  // Find t values on either side of time
  std::vector<Scalar>::iterator input_it = local_time_vec.begin();
  std::vector<DataStore<Scalar> >::iterator node_it = node_list.begin();
  for ( ; node_it != node_list.end() ; node_it++ )
  {
    while ((*input_it >= node_it->t) && (*input_it <= (node_it+1)->t))
    {
      Scalar& t = *input_it;
      Scalar& ti = node_it->t;
      Scalar& tip1 = (node_it+1)->t;
      Thyra::VectorBase<Scalar>& xi = *(node_it->x);
      Thyra::VectorBase<Scalar>& xip1 = *((node_it+1)->x);
      Thyra::VectorBase<Scalar>& xdoti = *(node_it->xdot);
      Thyra::VectorBase<Scalar>& xdotip1 = *((node_it+1)->xdot);

      // interpolate this point
      Scalar h = tip1-ti;
      // First we work on x.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),xi,Scalar(-ST::one()/h),xip1);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x = (node_it->x).clone_v();
      V_StVpStV(&*x, ST::one(), xi, t-ti, tmp_vec);
      x_vec.pushback(x);
      // Then we work on xdot.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),xdoti,Scalar(-ST::one()/h),xdotip1);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot = (node_it->xdot).clone_v();
      V_StVpStV(&*xdot, ST::one(), xdoti, t-ti, tmp_vec);
      xdot_vec.pushback(xdot);
      // And finally we estimate our order of accuracy
      accuracy_vec.pushback(h); 
      // Now we increment iterator for local_time_vec
      input_it++;
    }
  }
  return(true);
}

} // namespace Rythmos

#endif // Rythmos_LINEAR_INTERPOLATOR_H
