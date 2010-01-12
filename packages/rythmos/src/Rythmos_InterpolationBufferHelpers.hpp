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

#ifndef Rythmos_INTERPOLATION_BUFFER_HELPERS_HPP
#define Rythmos_INTERPOLATION_BUFFER_HELPERS_HPP


#include "Rythmos_InterpolationBufferBase.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief Assert that a time point vector is sorted.
 *
 * \relates InterpolationBufferBase
 */
template<class Scalar>
void assertTimePointsAreSorted(const Array<Scalar>& time_vec);


/** \brief Assert that none of the time points fall before the current
 * time range for an interpolation buffer object.
 *
 * \param interpBuffer [in] The interpolation buffer defining the time range.
 * The reason that we accept the interpolation buffer and not just the time
 * range is so we can create a better error message.
 *
 * \param time_vec [in] The array of time points
 *
 * \param startingTimePointIndex [in] The time index in <tt>time_vec</tt> to
 * begin asserting time points.  The default is 0.
 *
 * This function with throw an exception if any of the time points
 * are found to be before the current time range.
 *
 * \relates InterpolationBufferBase
 */
template<class Scalar>
void assertNoTimePointsBeforeCurrentTimeRange(
  const InterpolationBufferBase<Scalar> &interpBuffer,
  const Array<Scalar>& time_vec,
  const int &startingTimePointIndex = 0
  );

  
/** \brief Assert that none of the time points fall inside the current
 * time range for an interpolation buffer object.
 *
 * \param interpBuffer [in] The interpolation buffer defining the time range.
 * The reason that we accept the interpolation buffer and not just the time
 * range is so we can create a better error message.
 *
 * \param time_vec [in] The array of time points
 *
 * This function with throw an exception if any of the time points
 * are found to be inside the current time range.
 *
 * \relates InterpolationBufferBase
 */
template<class Scalar>
void assertNoTimePointsInsideCurrentTimeRange(
  const InterpolationBufferBase<Scalar> &interpBuffer,
  const Array<Scalar>& time_vec
  );


/** \brief Select points from an Array that sit in a TimeRange
 *
 * \relates InterpolationBufferBase
 */
template<class TimeType>
void selectPointsInTimeRange(
  const Array<TimeType>& points_in,
  const TimeRange<TimeType>& range,
  const Ptr<Array<TimeType> >& points_out 
  );


/** \brief Remove points from an Array that sit in a TimeRange
 *
 * \relates InterpolationBufferBase
 */
template<class TimeType>
void removePointsInTimeRange(
  Array<TimeType>* points_in, 
  const TimeRange<TimeType>& range 
  );


/** \brief Get time points in the current range of an interpolation buffer
 * object.
 *
 * \param interpBuffer [in] The interpolation buffer object that will be used
 * to get state values at different time points
 *
 * \param time_vec [in] The whole time vector, including time points that fall
 * outside of the current time range in <tt>interpBuffer</tt>.
 *
 * \param x_vec [in/out] This argument is optional and it is allowed for
 * <tt>x_vec==0</tt>.  However, if <tt>x_vec!=0</tt>, then on input
 * <tt>x_vec->size()==time_vec.size()</tt> must be true.  On output,
 * <tt>*x_vec</tt> will be filled with the state vectors for the current time
 * points given in <tt>time_vec</tt>.
 *
 * \param xdot_vec [in/out] This argument is optional and it is allowed for
 * <tt>xdot_vec==0</tt>.  However, if <tt>xdot_vec!=0</tt>, then on input
 * <tt>xdot_vec->size()==time_vec.size()</tt> must be true.  On output,
 * <tt>*xdot_vec</tt> will be filled with the state derivative vectors for the
 * current time points given in <tt>time_vec</tt>.
 *
 * \param nextTimePointIndex [in/out] On input, <tt>*nextTimePointIndex</tt>
 * gives first time point <tt>time_vec[*nextTimePointIndex]</tt> to extract
 * state values for.  On output, <tt>*nextTimePointIndex</tt> will be
 * incremented such that there are no more time points left (i.e.
 * <tt>nextTimePointIndex==time_vec.size()</tt> or such that
 * <tt>time_vec[*nextTimePointIndex] >
 * interpBuffer.getTimeRange().upper()</tt>
 *
 * <b>Preconditions:</b><ul>
 * <tt> <tt>nextTimePointIndex!=0</tt>
 * <tt> <tt>0 <= *nextTimePointIndex < time_vec.size()</tt>
 * <li> <tt>time_vec[*nextTimePointIndex] >= interpBuffer.getTimeRange().lower()</tt>
 * </ul>
 *
 * <b>Preconditions:</b><ul>
 * <li> [<tt>returnVal==true</tt>] <tt>*nextTimePointIndex</tt> is greater on output
 *      than on input
 * <li> [<tt>returnVal==false</tt>] <tt>*nextTimePointIndex</tt> is unchanged on output
 * </ul>
 *
 * \returns Returns <tt>true</tt> if one or more time points where extracted
 * from <tt>interpBuffer</tt> and returns <tt>false</tt> otherwise.
 *
 * \relates InterpolationBufferBase
 */
template<class Scalar>
bool getCurrentPoints(
  const InterpolationBufferBase<Scalar> &interpBuffer,
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  int *nextTimePointIndex
  );


} // namespace Rythmos


//
// Implementations
//


template<class Scalar>
void Rythmos::assertTimePointsAreSorted(const Array<Scalar>& time_vec)
{
  const int numTimePoints = time_vec.size();
  for ( int i = 0; i < numTimePoints-1; ++ i ) {
    TEST_FOR_EXCEPTION(
      time_vec[i] >= time_vec[i+1], std::logic_error,
      "Error, the time vector points time_vec["<<i<<"] = " << time_vec[i]
      << " >= time_vec["<<i+1<<"] = " << time_vec[i+1] << " are not [unique|sorted]!"
      );
  }
}


template<class Scalar>
void Rythmos::assertNoTimePointsBeforeCurrentTimeRange(
  const InterpolationBufferBase<Scalar> &interpBuffer,
  const Array<Scalar>& time_vec,
  const int &startingTimePointIndex
  )
{
  typedef ScalarTraits<Scalar> ST;
  const int numTimePoints = time_vec.size();
  const TimeRange<Scalar> currentTimeRange = interpBuffer.getTimeRange();
  if (currentTimeRange.length() >= ST::zero()) {
    for ( int i = 0; i < numTimePoints; ++i ) {
      TEST_FOR_EXCEPTION(
        time_vec[i] < currentTimeRange.lower(), std::out_of_range,
        "Error, time_vec["<<i<<"] = " << time_vec[i] << " < currentTimeRange.lower() = "
        << currentTimeRange.lower() << " for " << interpBuffer.description() << "!"
        );
    }
  }
}


template<class Scalar>
void Rythmos::assertNoTimePointsInsideCurrentTimeRange(
  const InterpolationBufferBase<Scalar>& interpBuffer,
  const Array<Scalar>& time_vec
  )
{
  typedef ScalarTraits<Scalar> ST;
  const int numTimePoints = time_vec.size();
  const TimeRange<Scalar> currentTimeRange = interpBuffer.getTimeRange();
  if (currentTimeRange.length() >= ST::zero()) {
    for ( int i = 0; i < numTimePoints; ++i ) {
      TEST_FOR_EXCEPTION(
        currentTimeRange.isInRange(time_vec[i]), std::out_of_range,
        "Error, time_vec["<<i<<"] = " << time_vec[i] << " is in TimeRange of " 
        << interpBuffer.description() << " = ["
        << currentTimeRange.lower() << "," << currentTimeRange.upper() << "]!"
        );
    }
  }
}


template<class TimeType>
void Rythmos::selectPointsInTimeRange(
    const Array<TimeType>& points_in,
    const TimeRange<TimeType>& range,
    const Ptr<Array<TimeType> >& points_out 
    )
{
  points_out->clear();
  int Nt = Teuchos::as<int>(points_in.size());
  for (int i=0; i < Nt ; ++i) {
    if (range.isInRange(points_in[i])) {
      points_out->push_back(points_in[i]);
    }
  }
}


template<class TimeType>
void Rythmos::removePointsInTimeRange(
    Array<TimeType>* points_in, 
    const TimeRange<TimeType>& range 
    )
{
  Array<TimeType> values_to_remove;
  for (int i=0 ; i<Teuchos::as<int>(points_in->size()) ; ++i) {
    if (range.isInRange((*points_in)[i])) {
      values_to_remove.push_back((*points_in)[i]);
    }
  }
  typename Array<TimeType>::iterator point_it;
  for (int i=0 ; i< Teuchos::as<int>(values_to_remove.size()) ; ++i) {
    point_it = std::find(points_in->begin(),points_in->end(),values_to_remove[i]);
    TEST_FOR_EXCEPTION(
        point_it == points_in->end(), std::logic_error,
        "Error, point to remove = " << values_to_remove[i] << " not found with std:find!\n"
        );
    points_in->erase(point_it);
  }
}


template<class Scalar>
bool Rythmos::getCurrentPoints(
  const InterpolationBufferBase<Scalar> &interpBuffer,
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  int *nextTimePointIndex_inout
  )
{

  typedef ScalarTraits<Scalar> ST;
  using Teuchos::as;

  const int numTotalTimePoints = time_vec.size();

  // Validate input
#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT(nextTimePointIndex_inout==0);
  TEUCHOS_ASSERT( 0 <= *nextTimePointIndex_inout && *nextTimePointIndex_inout < numTotalTimePoints );
  TEUCHOS_ASSERT( x_vec == 0 || as<int>(x_vec->size()) == numTotalTimePoints );
  TEUCHOS_ASSERT( xdot_vec == 0 || as<int>(xdot_vec->size()) == numTotalTimePoints );
#endif // RYTHMOS_DEBUG

  int &nextTimePointIndex = *nextTimePointIndex_inout;
  const int initNextTimePointIndex = nextTimePointIndex;

  const TimeRange<Scalar> currentTimeRange = interpBuffer.getTimeRange();
  
  if (currentTimeRange.length() >= ST::zero()) {

    // Load a temp array with all of the current time points that fall in the
    // current time range.
    Array<Scalar> current_time_vec;
    { // scope for i to remove shadow warning.
      int i;
      for ( i = 0; i < numTotalTimePoints-nextTimePointIndex; ++i ) {
        const Scalar t = time_vec[nextTimePointIndex];
#ifdef RYTHMOS_DEBUG
        TEUCHOS_ASSERT( t >= currentTimeRange.lower() );
#endif // RYTHMOS_DEBUG
        if ( currentTimeRange.isInRange(t) ) {
          ++nextTimePointIndex;
          current_time_vec.push_back(t);
        }
        else {
          break;
        }
      }
#ifdef RYTHMOS_DEBUG
      // Here I am just checking that the loop worked as expected with the data
      // in the current time range all comming first.
      TEUCHOS_ASSERT( nextTimePointIndex-initNextTimePointIndex == i );
#endif
    }

    // Get points in current time range if any such points exist

    const int numCurrentTimePoints = current_time_vec.size();

    if ( numCurrentTimePoints > 0 ) {

      // Get the state(s) for current time points from the stepper and put
      // them into temp arrays
      Array<RCP<const Thyra::VectorBase<Scalar> > > current_x_vec;
      Array<RCP<const Thyra::VectorBase<Scalar> > > current_xdot_vec;
      if (x_vec || xdot_vec) {
        interpBuffer.getPoints(
          current_time_vec,
          x_vec ? &current_x_vec : 0,
          xdot_vec ? &current_xdot_vec : 0,
          0 // accuracy_vec
          );
      }

      // Copy the gotten x and xdot vectors from the temp arrays to the output
      // arrays.
      for ( int i = initNextTimePointIndex; i < nextTimePointIndex; ++i ) {
        if (x_vec)
          (*x_vec)[i] = current_x_vec[i-initNextTimePointIndex];
        if (xdot_vec)
          (*xdot_vec)[i] = current_xdot_vec[i-initNextTimePointIndex];
      }

    }

  }

  return ( nextTimePointIndex == initNextTimePointIndex ? false : true );

}


#endif //Rythmos_INTERPOLATION_BUFFER_HELPERS_HPP
