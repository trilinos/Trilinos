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

#ifndef RYTHMOS_TIME_RANGE_UNITTEST_H
#define RYTHMOS_TIME_RANGE_UNITTEST_H

#include "../test/CppUnitLite/TestHarness.h"

#include "Rythmos_TimeRange.hpp"

namespace Rythmos {

TEST( TimeRange, newTimeRange ) {
  TimeRange<double> tr;
  // it should be initialized as [0,-1]
  CHECK( !(tr.isValid()) );
  CHECK( tr.lower() > tr.upper() );
  CHECK( !(tr.isInRange(0.5)) );
  CHECK( !tr.isInRange(0.0) );
  CHECK( !tr.isInRange(-1.0) );
  CHECK( tr.length() == -1.0 );
}

TEST( TimeRange, isValid ) {
  TimeRange<double> tr(0.0,1.0);
  CHECK( tr.isValid() );
  CHECK( tr.isInRange(0.0) );
  CHECK( tr.isInRange(0.5) );
  CHECK( tr.isInRange(1.0) );
  CHECK( tr.length() == 1.0 );
}

TEST( TimeRange, copyAndScale ) {
  TimeRange<double> tr(1.0,2.0);
  TimeRange<double> newTr = tr.copyAndScale(5.0);
  CHECK( newTr.isValid() );
  CHECK( newTr.lower() == 5.0 );
  CHECK( newTr.upper() == 10.0 );
  CHECK( newTr.length() == 5.0 );
}

TEST( TimeRange, copyAndScaleInvalid ) {
  TimeRange<double> tr;
  TimeRange<double> newTr = tr.copyAndScale(5.0);
  CHECK( !(newTr.isValid()) );
  CHECK( newTr.lower() == tr.lower() );
  CHECK( newTr.upper() == tr.upper() );
  CHECK( newTr.length() == tr.length() );
}

TEST( TimeRange, nonMemberConstructor ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  CHECK( tr.isValid() );
  CHECK( tr.lower() == 1.25 );
  CHECK( tr.upper() == 3.45 );
}

TEST( TimeRange, invalidTimeRange ) {
  TimeRange<double> tr = invalidTimeRange<double>();
  CHECK( !tr.isValid() );
  CHECK( tr.lower() > tr.upper() );
  CHECK( !(tr.isInRange(0.5)) );
}

// How do I check the << stream operator?

TEST( TimeRange, cc ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  CHECK( isInRange_cc(tr, 1.25) );
  CHECK( isInRange_cc(tr, 2.0) );
  CHECK( isInRange_cc(tr, 3.45) );
}

TEST( TimeRange, oc ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  CHECK( !isInRange_oc(tr, 1.25) );
  CHECK( isInRange_oc(tr, 2.0) );
  CHECK( isInRange_oc(tr, 3.45) );
}

TEST( TimeRange, co ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  CHECK( isInRange_co(tr, 1.25) );
  CHECK( isInRange_co(tr, 2.0) );
  CHECK( !isInRange_co(tr, 3.45) );
}

TEST( TimeRange, oo ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  CHECK( !isInRange_oo(tr, 1.25) );
  CHECK( isInRange_oo(tr, 2.0) );
  CHECK( !isInRange_oo(tr, 3.45) );
}

} // namespace Rythmos


#endif // RYTHMOS_TIME_RANGE_UNITTEST_H

