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

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_TimeRange.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, newTimeRange ) {
  TimeRange<double> tr;
  // it should be initialized as [0,-1]
  TEST_EQUALITY_CONST( tr.isValid(), false );
  TEST_COMPARE( tr.lower(), >, tr.upper() );
  TEST_EQUALITY_CONST( tr.isInRange(0.5), false );
  TEST_EQUALITY_CONST( tr.isInRange(0.0), false );
  TEST_EQUALITY_CONST( tr.isInRange(-1.0), false );
  TEST_EQUALITY_CONST( tr.length(), -1.0 );
}

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, isValid ) {
  TimeRange<double> tr(0.0,1.0);
  TEST_EQUALITY_CONST( tr.isValid(), true );
  TEST_EQUALITY_CONST( tr.isInRange(0.0), true );
  TEST_EQUALITY_CONST( tr.isInRange(0.5), true );
  TEST_EQUALITY_CONST( tr.isInRange(1.0), true );
  TEST_EQUALITY_CONST( tr.length(), 1.0 );
}

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, copyAndScale ) {
  TimeRange<double> tr(1.0,2.0);
  TimeRange<double> newTr = tr.copyAndScale(5.0);
  TEST_EQUALITY_CONST( newTr.isValid(), true );
  TEST_EQUALITY_CONST( newTr.lower(), 5.0 );
  TEST_EQUALITY_CONST( newTr.upper(), 10.0 );
  TEST_EQUALITY_CONST( newTr.length(), 5.0 );
}

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, copyAndScaleInvalid ) {
  TimeRange<double> tr;
  TimeRange<double> newTr = tr.copyAndScale(5.0);
  TEST_EQUALITY_CONST( newTr.isValid(), false );
  TEST_EQUALITY( newTr.lower(), tr.lower() );
  TEST_EQUALITY( newTr.upper(), tr.upper() );
  TEST_EQUALITY( newTr.length(), tr.length() );
}

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, nonMemberConstructor ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  TEST_EQUALITY_CONST( tr.isValid(), true );
  TEST_EQUALITY_CONST( tr.lower(), 1.25 );
  TEST_EQUALITY_CONST( tr.upper(), 3.45 );
}

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, invalidTimeRange ) {
  TimeRange<double> tr = invalidTimeRange<double>();
  TEST_EQUALITY_CONST( tr.isValid(), false );
  TEST_COMPARE( tr.lower(), >, tr.upper() );
  TEST_EQUALITY_CONST( tr.isInRange(0.5), false );
}

// How do I check the << stream operator?

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, cc ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  TEST_EQUALITY_CONST( isInRange_cc(tr, 1.25), true );
  TEST_EQUALITY_CONST( isInRange_cc(tr, 2.0), true );
  TEST_EQUALITY_CONST( isInRange_cc(tr, 3.45), true );
}

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, oc ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  TEST_EQUALITY_CONST( isInRange_oc(tr, 1.25), false );
  TEST_EQUALITY_CONST( isInRange_oc(tr, 2.0), true );
  TEST_EQUALITY_CONST( isInRange_oc(tr, 3.45), true );
}

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, co ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  TEST_EQUALITY_CONST( isInRange_co(tr, 1.25), true );
  TEST_EQUALITY_CONST( isInRange_co(tr, 2.0), true );
  TEST_EQUALITY_CONST( isInRange_co(tr, 3.45), false );
}

TEUCHOS_UNIT_TEST( Rythmos_TimeRange, oo ) {
  TimeRange<double> tr = timeRange(1.25,3.45);
  TEST_EQUALITY_CONST( isInRange_oo(tr, 1.25), false );
  TEST_EQUALITY_CONST( isInRange_oo(tr, 2.0), true );
  TEST_EQUALITY_CONST( isInRange_oo(tr, 3.45), false );
}

} // namespace Rythmos



