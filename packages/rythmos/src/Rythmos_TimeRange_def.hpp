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

#ifndef RYTHMOS_TIME_RANGE_DEF_H
#define RYTHMOS_TIME_RANGE_DEF_H

#include "Rythmos_TimeRange_decl.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_ScalarTraits.hpp"


template<class TimeType>
int Rythmos::compareTimeValues( const TimeType &t1, const TimeType &t2 )
{
  // Here we will do the comparison based on the magnitude of t1
  const TimeType epsMore = 10.0*std::numeric_limits<TimeType>::epsilon();
  const TimeType t1Mag = Teuchos::ScalarTraits<TimeType>::magnitude(t1);
  const TimeType t1Tol = t1Mag*epsMore;
  if ( t2 - t1Tol <= t1 && t1 <= t2 + t1Tol )
    return 0;
  else if ( t1 > t2 + t1Tol )
    return +1;
  // t1 < t2 - t1Tol
  return -1;
}


template<class TimeType>
Rythmos::TimeRange<TimeType>
Rythmos::timeRange(const TimeType lower, const TimeType upper)
{
  return TimeRange<TimeType>(lower,upper);
}


template<class TimeType>
Rythmos::TimeRange<TimeType>
Rythmos::invalidTimeRange()
{
  return TimeRange<TimeType>();
}


template<class TimeType>
std::ostream&
Rythmos::operator<<( std::ostream& out, const TimeRange<TimeType>& range )
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


template<class TimeType>
void Rythmos::asssertInTimeRange( const TimeRange<TimeType> &timeRange,
  const TimeType &time )
{
  TEST_FOR_EXCEPTION( !timeRange.isInRange(time), std::out_of_range,
    "Error, the time = " << time
    << " is out of the range = " << timeRange << "!"
    );
}


template<class TimeType>
bool Rythmos::isInRange_cc(const TimeRange<TimeType> &tr, const TimeType &p)
{
  return (
    compareTimeValues(p,tr.lower()) >= 0
    && compareTimeValues(p,tr.upper()) <= 0
    );
}


template<class TimeType>
bool Rythmos::isInRange_oc(const TimeRange<TimeType> &tr, const TimeType &p)
{
  return (
    compareTimeValues(p,tr.lower()) > 0
    && compareTimeValues(p,tr.upper()) <= 0
    );
}


template<class TimeType>
bool Rythmos::isInRange_co(const TimeRange<TimeType> &tr, const TimeType &p)
{
  return (
    compareTimeValues(p,tr.lower()) >= 0
    && compareTimeValues(p,tr.upper()) < 0
    );
}


template<class TimeType>
bool Rythmos::isInRange_oo(const TimeRange<TimeType> &tr, const TimeType &p)
{
  return (
    compareTimeValues(p,tr.lower()) > 0
    && compareTimeValues(p,tr.upper()) < 0
    );
}


#define RYTHMOS_TIME_RANGE_INSTANT(SCALAR) \
  \
  template class TimeRange< SCALAR >; \
  \
  template int compareTimeValues( const  SCALAR  &t1, const  SCALAR  &t2 ); \
  template TimeRange< SCALAR > timeRange(const  SCALAR  lower, const  SCALAR  upper); \
  template TimeRange< SCALAR > invalidTimeRange(); \
  template std::ostream& operator<<( std::ostream& out, const TimeRange< SCALAR >& range ); \
  template void asssertInTimeRange( const TimeRange<SCALAR > &timeRange, const SCALAR &time ); \
  template bool isInRange_cc(const TimeRange< SCALAR > &tr, const  SCALAR  &p); \
  template bool isInRange_oc(const TimeRange< SCALAR > &tr, const  SCALAR  &p); \
  template bool isInRange_co(const TimeRange< SCALAR > &tr, const  SCALAR  &p); \
  template bool isInRange_oo(const TimeRange< SCALAR > &tr, const  SCALAR  &p); \
  template class TimeRange_cc< SCALAR >;  \
  template class TimeRange_co< SCALAR >;  \
  template class TimeRange_oo< SCALAR >;  \
  template class TimeRange_oc< SCALAR >;  


#endif //RYTHMOS_TIME_RANGE_DEF_H
