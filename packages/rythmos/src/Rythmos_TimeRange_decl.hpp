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

#ifndef RYTHMOS_TIME_RANGE_DECL_H
#define RYTHMOS_TIME_RANGE_DECL_H

#include "Rythmos_ConfigDefs.h"

namespace Rythmos {


/** \brief Compare two times taking into account floating point errors.
 *
 * \returns Return <tt> -1 </tt> if <tt> t1 < t2  </tt>, 
 *                 <tt>  0 </tt> if <tt> t1 == t2 </tt> and 
 *                 <tt> +1 </tt> if <tt> t1 > t2  </tt>.
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
 *
 * \relates TimeRange
 */
template<class TimeType>
int compareTimeValues( const TimeType &t1, const TimeType &t2 );

/** \brief Represent a time range.
 *
 * The compiler-generated default constructor, copy constructor, and
 * assignment operators are allowed and perform correctly.
 *
 * ToDo: Put in checks for the range if needed.
 */
template<class TimeType>
class TimeRange 
{
public:
  /** \brief Construct an invalid range. */
  TimeRange()
    : lower_(0.0), upper_(-1.0)
    {}
  /** \brief Construct a valid range. */
  TimeRange( const TimeType &my_lower, const TimeType &my_upper )
    : lower_(my_lower), upper_(my_upper)
    {
    }
  /** \brief Copy constructor. */
  TimeRange( const TimeRange<TimeType>& tr )
    : lower_(tr.lower()), upper_(tr.upper())
    {
    }
  /** \brief . */
  virtual ~TimeRange() {}
  /** \brief . */
  bool isValid() const { return (lower_ <= upper_); }
  /** \brief . */
  TimeType lower() const { return lower_; }
  /** \brief . */
  TimeType upper() const { return upper_; }
  /** \brief . */
  TimeType length() const { return (upper_ - lower_); }
  /** \brief . */
  virtual bool isInRange ( const TimeType &t ) const
    {
      return (
        compareTimeValues(t,lower_) >= 0
        && compareTimeValues(t,upper_) <= 0
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
TimeRange<TimeType> timeRange(const TimeType my_lower, const TimeType my_upper);


/** \brief Nonmember constructor.
 *
 * \relates TimeRange
 */
template<class TimeType>
TimeRange<TimeType> invalidTimeRange();


/** \brief Output operator.
 *
 * \relates TimeRange
 */
template<class TimeType>
std::ostream& operator<<( std::ostream& out, const TimeRange<TimeType>& range );


/** \brief Assert point is in the time range and if it is not, throw an
 * exception with a very good error message.
 */
template<class TimeType>
void asssertInTimeRange( const TimeRange<TimeType> &timeRange,
  const TimeType &time );


/** \brief Nonmember isInRange function [closed, closed].
 *
 * <tt>tr.lower() <= p <= tr.upper()</tt>
 *
 * \relates TimeRange
 */
template<class TimeType>
bool isInRange_cc(const TimeRange<TimeType> &tr, const TimeType &p);


/** \brief Nonmember isInRange function (open, closed].
 *
 * <tt>tr.lower() < p <= tr.upper()</tt>
 *
 * \relates TimeRange
 */
template<class TimeType>
bool isInRange_oc(const TimeRange<TimeType> &tr, const TimeType &p);


/** \brief Nonmember isInRange function [closed, open).
 *
 * <tt>tr.lower() <= p < tr.upper()</tt>
 *
 * \relates TimeRange
 */
template<class TimeType>
bool isInRange_co(const TimeRange<TimeType> &tr, const TimeType &p);


/** \brief Nonmember isInRange function (open, open).
 *
 * <tt>tr.lower() < p < tr.upper()</tt>
 *
 * \relates TimeRange
 */
template<class TimeType>
bool isInRange_oo(const TimeRange<TimeType> &tr, const TimeType &p);

template<class TimeType>
class TimeRange_cc : virtual public TimeRange<TimeType> 
{
public:
  TimeRange_cc(const TimeRange<TimeType>& tr)
    :TimeRange<TimeType>(tr)
    {
    }
  TimeRange_cc( const TimeType &my_lower, const TimeType &my_upper )
    :TimeRange<TimeType>(my_lower,my_upper) 
    {
    }
  bool isInRange ( const TimeType &t ) const
    {
      return ( isInRange_cc<TimeType>(*this,t) );
    }
};

template<class TimeType>
class TimeRange_co : virtual public TimeRange<TimeType> 
{
public:
  TimeRange_co(const TimeRange<TimeType>& tr)
    :TimeRange<TimeType>(tr)
    {
    }
  TimeRange_co( const TimeType &my_lower, const TimeType &my_upper )
    :TimeRange<TimeType>(my_lower,my_upper) 
    {
    }
  bool isInRange ( const TimeType &t ) const
    {
      return ( isInRange_co<TimeType>(*this,t) );
    }
};

template<class TimeType>
class TimeRange_oo : virtual public TimeRange<TimeType> 
{
public:
  TimeRange_oo(const TimeRange<TimeType>& tr)
    :TimeRange<TimeType>(tr)
    {
    }
  TimeRange_oo( const TimeType &my_lower, const TimeType &my_upper )
    :TimeRange<TimeType>(my_lower,my_upper) 
    {
    }
  bool isInRange ( const TimeType &t ) const
    {
      return ( isInRange_oo<TimeType>(*this,t) );
    }
};

template<class TimeType>
class TimeRange_oc : virtual public TimeRange<TimeType> 
{
public:
  TimeRange_oc(const TimeRange<TimeType>& tr)
    :TimeRange<TimeType>(tr)
    {
    }
  TimeRange_oc( const TimeType &my_lower, const TimeType &my_upper )
    :TimeRange<TimeType>(my_lower,my_upper) 
    {
    }
  bool isInRange ( const TimeType &t ) const
    {
      return ( isInRange_oc<TimeType>(*this,t) );
    }
};


} // namespace Rythmos


#endif //RYTHMOS_TIME_RANGE_DECL_H
