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

#ifndef RYTHMOS_TIME_RANGE_H
#define RYTHMOS_TIME_RANGE_H

#include "Rythmos_ConfigDefs.h"
#include "Teuchos_Assert.hpp"

namespace Rythmos {


/** \brief Compare two times taking into account floating point errors.
 *
 * \returns Return <tt> -1 </tt> if <tt> v1 < v2  </tt>, 
 *                 <tt>  0 </tt> if <tt> v1 == v2 </tt> and 
 *                 <tt> +1 </tt> if <tt> v1 > v2  </tt>.
 *
 * Note that this function compares t1 to t2 and not the other way around.
 *
 * This function is designed to solve the problem where in:
 
 \code

 const TimeType timeStep = endTime - time;
 const TimeType updatedEndTime = time + timeStep;
 const bool isAtOrPastEndTime = ( updatedEndTime >= endTime );

 \endcode

 * the bool <tt>isAtOrPastEndTime</tt> may actually be <tt>false</tt>!  This
 * is especially a danger if IEEE floating point numerics are not being used
 * (and there are many systems that do not for speed and other reasons).
 */
template<class TimeType>
int compareTimeValues( const TimeType &t1, const TimeType &t2 )
{
  // Here we will do the comparison based on the magnitude of t1
  typedef Teuchos::ScalarTraits<TimeType> ST;
  const TimeType epsMore = 10.0*std::numeric_limits<TimeType>::epsilon();
  const TimeType t1Mag = ST::magnitude(t1);
  const TimeType t1Tol = t1Mag*epsMore;
  if ( t2 - t1Tol <= t1 && t1 <= t2 + t1Tol )
    return 0;
  else if ( t1 > t2 + t1Tol )
    return +1;
  // t1 < t2 - t1Tol
  return -1;
}


/** \brief Represent a time range.
 *
 * The compiler-generated default constructor, copy constructor, and
 * assignment operators are allowed and perform correctly.
 *
 * ToDo: Put in checks for the range if needed.
 */
template<class TimeType>
class TimeRange {
public:
  /** \brief Construct an invalid range. */
  TimeRange()
    : lower_(0.0), upper_(-1.0)
    {}
  /** \brief Construct a valid range. */
  TimeRange( const TimeType &lower, const TimeType &upper )
    : lower_(lower), upper_(upper)
    {
    }
  /** \brief . */
  bool isValid() const { return (lower_ <= upper_); }
  /** \brief . */
  TimeType lower() const { return lower_; }
  /** \brief . */
  TimeType upper() const { return upper_; }
  /** \brief . */
  TimeType length() const { return (upper_ - lower_); }
  /** \brief . */
  bool isInRange ( const TimeType &t ) const
    {
      return (
        compareTimeValues(lower_,t) <= 0
        && compareTimeValues(upper_,t) >= 0
        );
    }
  /** \brief . */
  TimeRange<TimeType> copyAndScale( const TimeType &scale ) const
    {
      TimeRange<TimeType> newRange = *this;
      if (!newRange.isValid())
        return newRange;
      newRange.lower_ *= scale;
      newRange.upper_ *= scale;
      return newRange;
    }

private:
  TimeType lower_;
  TimeType upper_;
};


/** \brief Nonmember constructor.
 *
 * \relates TimeRange
 */
template<class TimeType>
inline
TimeRange<TimeType> timeRange(const TimeType lower, const TimeType upper)
{
  return TimeRange<TimeType>(lower,upper);
}


/** \brief Nonmember constructor.
 *
 * \relates TimeRange
 */
template<class TimeType>
inline
TimeRange<TimeType> invalidTimeRange()
{
  return TimeRange<TimeType>();
}


/** \brief Output operator.
 *
 * \relates TimeRange
 */
template<class TimeType>
std::ostream& operator<<( std::ostream& out, const TimeRange<TimeType>& range )
{
  out << "[";
  if (range.isValid()) {
    out << range.lower() << "," << range.upper();
  }
  else {
    out <<"INVALID";
  }
  out << "]";
  return out;
}


/** \brief Nonmember isInRange function [closed, closed].
 *
 * <tt>tr.lower() <= p <= tr.upper()</tt>
 *
 * \relates TimeRange
 */
template<class TimeType>
bool isInRange_cc(const TimeRange<TimeType> &tr, const TimeType &p)
{
  return (
    compareTimeValues(tr.lower(),p) <= 0
    && compareTimeValues(tr.upper(),p) >= 0
    );
}


/** \brief Nonmember isInRange function (open, closed].
 *
 * <tt>tr.lower() < p <= tr.upper()</tt>
 *
 * \relates TimeRange
 */
template<class TimeType>
bool isInRange_oc(const TimeRange<TimeType> &tr, const TimeType &p)
{
  return (
    compareTimeValues(tr.lower(),p) < 0
    && compareTimeValues(tr.upper(),p) >= 0
    );
}


/** \brief Nonmember isInRange function [closed, open).
 *
 * <tt>tr.lower() <= p < tr.upper()</tt>
 *
 * \relates TimeRange
 */
template<class TimeType>
bool isInRange_co(const TimeRange<TimeType> &tr, const TimeType &p)
{
  return (
    compareTimeValues(tr.lower(),p) <= 0
    && compareTimeValues(tr.upper(),p) > 0
    );
}


/** \brief Nonmember isInRange function (open, open).
 *
 * <tt>tr.lower() < p < tr.upper()</tt>
 *
 * \relates TimeRange
 */
template<class TimeType>
bool isInRange_oo(const TimeRange<TimeType> &tr, const TimeType &p)
{
  return (
    compareTimeValues(tr.lower(),p) < 0
    && compareTimeValues(tr.upper(),p) > 0
    );
}


} // namespace Rythmos


#endif //RYTHMOS_TIME_RANGE_H
