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

#ifndef Rythmos_INTERPOLATION_BUFFER_LINEAR_H
#define Rythmos_INTERPOLATION_BUFFER_LINEAR_H

#include "Rythmos_InterpolationBuffer.hpp"
#include "Thyra_VectorBase.hpp"

namespace Rythmos {

/** \brief class for defining linear interpolation buffer functionality. */
template<class Scalar> 
class LinearInterpolationBuffer : virtual public InterpolationBuffer<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    
    /** \brief. */
    LinearInterpolationBuffer();
    LinearInterpolationBuffer( ScalarMag accuracy, int storage );

    // 06/29/06 tscoffe:  The idea of accuracy is not something that the buffer
    // can really control.  It is related to the intrinsic accuracy of
    // interpolation accuracy and the local spacing of the node points.  So the
    // accuracy of linear interpolation is first order and hence O(h) with h
    // being the node spacing (which will change for different intervals).
    SetAccuracy( ScalarMag accuracy );
    SetStorage( int storage );
        
    /// Destructor
    ~LinearInterpolationBuffer() {};

    /// Add point to buffer
    bool SetPoints(
      const std::vector<Scalar>& time_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& x_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& xdot_list);

    /// Get value from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_list
      ,std::vector<Thyra::VectorBase<Scalar> >* x_list
      ,std::vector<Thyra::VectorBase<Scalar> >* xdot_list
      ,std::vector<ScalarMag>* accuracy_list) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBuffer<Scalar>& IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_list) const;

    /// Get order of interpolation
    int GetOrder() const;

  private:

    std::vector<Scalar> indep_values
    std::vector<std::vector<Thyra::VectorBase<Scalar> > > dep_values;
    ScalarMag accuracy;
    int storage_limit;
};

// ////////////////////////////
// Defintions

template<class Scalar>
LinearInterpolationBuffer<Scalar>::LinearInterpolationBuffer()
{
  SetAccuracy(Scalar(1.0e-4));
  SetStorage(2);
}

template<class Scalar>
LinearInterpolationBuffer<Scalar>::LinearInterpolationBuffer( ScalarMag accuracy_, int storage_ )
{
  SetAccuracy(accuracy_);
  SetStorage(storate_);
}

template<class Scalar>
LinearInterpolationBuffer<Scalar>::SetAccuracy( ScalarMag accuracy_ )
{
  accuracy = accuracy_;
}

template<class Scalar>
LinearInterpolationBuffer<Scalar>::SetStorage( int storage_ )
{
  storage_limit = min(2,storage_); // Minimum of two points so interpolation is possible
}

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::SetPoints( 
    const std::vector<Scalar>& time_list
    ,const std::vector<Thyra::VectorBase<Scalar> >& x
    ,const std::vector<Thyra::VectorBase<Scalar> >& xdot );
{
  // Determine if time is already in list and if so, replace existing data with new data.
  // If we're already at our max_storage limit, then report failure.
  // Determine where time should be in list and insert it along with x and xdot.
  return(false);
}

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::GetPoints(
    const std::vector<Scalar>& time_list
    ,std::vector<Thyra::VectorBase<Scalar> >* x
    ,std::vector<Thyra::VectorBase<Scalar> >* xdot
    ,std::vector<ScalarMag>* accuracy_list) const
{
  // If time is outside range of indep_values, then return failure
  // Find indep_values on either side of time
  // Do interpolation, storing data in provided x and xdot pointers and return success
  return(false);
}

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::SetRange(
    const ScalarMagRange& time_range
    ,const InterpolationBuffer<Scalar>& IB )
{
  // If IB has a sense of accuracy, its lower than our sense, and its node
  // points are too far apart, then it will be impossible for us to maintain
  // our accuracy level while importing this IB, so return failure.
  // Otherwise, grab node values from IB and copy those over if they're in the interval of time_range.
  // If these points are too far apart for our accuracy level, then ask IB to interpolate more.
  // Use SetPoint and check return value to make sure we observe storage_limit.
  return(false);
}

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::GetNodes( std::vector<Scalar>* time_list ) const
{
  time_list = indep_values; // std::copy of data (this may be incorrect)
  return(true);
}

template<class Scalar>
int LinearInterpolationBuffer<Scalar>::GetOrder() const
{
  return(1);
}

} // namespace Rythmos

#endif //Rythmos_INTERPOLATION_BUFFER_LINEAR_H
