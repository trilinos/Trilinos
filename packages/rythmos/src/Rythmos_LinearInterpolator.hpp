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
        const std::vector<Teuchos::RefCountPtr<DataStore<Scalar> > > &data_in
        ,const std::vector<Scalar> &t_values
        ,std::vector<Teuchos::RefCountPtr<DataStore<Scalar> > > *data_out
        ) const;

    /// Order of interpolation:
    int order() const; 

  private:

    Thyra::VectorBase<Scalar> tmp_vec;

};

template<class Scalar>
bool LinearInterpolator<Scalar>::interpolate(
    const std::vector<Teuchos::RefCountPtr<DataStore<Scalar> > > &data_in
    ,const std::vector<Scalar> &t_values
    ,std::vector<Teuchos::RefCountPtr<DataStore<Scalar> > > *data_out ) const
{
  // If there are fewer than 2 points in data_in, then return failure
  if (data_in.size() < 2)
    return(false);
  // Sort data_in: 
  std::vector<Scalar> local_data_in = data_in;
  std::sort(local_data_in.begin(),local_data_in.end());
  // Sort t_values:
  std::vector<Scalar> local_time_vec = t_values;
  std::sort(local_time_vec.begin(),local_time_vec.end());
  // If time is outside range of t_values, then return failure
  if ( (*local_time_vec.begin() < local_data_in.begin().time) || (*local_time_vec.end() > local_data_in.end().time) )
    return(false);
  // Find t values on either side of time
  std::vector<Scalar>::iterator input_it = local_time_vec.begin();
  std::vector<DataStore<Scalar> >::iterator node_it = local_data_in.begin();
  for ( ; node_it != local_data_in.end() ; node_it++ )
  {
    while ((*input_it >= node_it->time) && (*input_it <= (node_it+1)->time))
    {
      Scalar& t = *input_it;
      Scalar& ti = node_it->time;
      Scalar& tip1 = (node_it+1)->time;
      Thyra::VectorBase<Scalar>& xi = *(node_it->x);
      Thyra::VectorBase<Scalar>& xip1 = *((node_it+1)->x);
      Thyra::VectorBase<Scalar>& xdoti = *(node_it->xdot);
      Thyra::VectorBase<Scalar>& xdotip1 = *((node_it+1)->xdot);

      Thyra::VectorBase<Scalar> tmp_vec;

      // interpolate this point
      Teuchos::RefCountPtr<DataStore<Scalar> > DS = Teuchos::rcp(new DataStore<Scalar>);
      DS->time = t;
      Scalar h = tip1-ti;
      // First we work on x.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),xi,Scalar(-ST::one()/h),xip1);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x = (node_it->x).clone_v();
      V_StVpStV(&*x, ST::one(), xi, t-ti, tmp_vec);
      DS->x = x;
      // Then we work on xdot.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),xdoti,Scalar(-ST::one()/h),xdotip1);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot = (node_it->xdot).clone_v();
      V_StVpStV(&*xdot, ST::one(), xdoti, t-ti, tmp_vec);
      DS->xdot = xdot;
      // And finally we estimate our order of accuracy
      DS->accuracy = h;
      // Now we increment iterator for local_time_vec
      input_it++;
      data_out->push_back(DS);
    }
  }
  return(true);
}

template<class Scalar>
int LinearInterpolator<Scalar>::order()
{
  return(1);
}

} // namespace Rythmos

#endif // Rythmos_LINEAR_INTERPOLATOR_H
