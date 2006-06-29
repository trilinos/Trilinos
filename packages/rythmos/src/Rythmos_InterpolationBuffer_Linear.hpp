//
// @HEADER
// ***********************************************************************
// 
//                           Rythmos Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

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
    bool SetPoint(const ScalarMag& time, const Thyra::VectorBase<Scalar>& x, const Thyra::VectorBase<Scalar>& xdot);

    /// Get value from buffer
    bool GetPoint(const ScalarMag& time, Thyra::VectorBase<Scalar>* x, Thyra::VectorBase<Scalar>* xdot) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(const ScalarMagRange& time_range, const InterpolationBuffer<Scalar>& IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<ScalarMag>* time_list) const;

  private:

    std::vector<ScalarMag> indep_values
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
bool LinearInterpolationBuffer<Scalar>::SetPoint( 
    const ScalarMag& time
    ,const Thyra::VectorBase<Scalar>& x
    ,const Thyra::VectorBase<Scalar>& xdot );
{
  // Determine if time is already in list and if so, replace existing data with new data.
  // If we're already at our max_storage limit, then report failure.
  // Determine where time should be in list and insert it along with x and xdot.
  return(false);
}

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::GetPoint(
    const ScalarMag& time
    ,Thyra::VectorBase<Scalar>* x
    ,Thyra::VectorBase<Scalar>* xdot) const
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
bool LinearInterpolationBuffer<Scalar>::GetNodes( std::vector<ScalarMag>* time_list ) const
{
  time_list = indep_values; // std::copy of data (this may be incorrect)
  return(true);
}

} // namespace Rythmos

#endif //Rythmos_INTERPOLATION_BUFFER_LINEAR_H
