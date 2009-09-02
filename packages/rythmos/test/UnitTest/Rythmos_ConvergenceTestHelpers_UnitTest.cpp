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

#include "Teuchos_UnitTestHarness.hpp"

#include "../ConvergenceTest/Rythmos_ConvergenceTestHelpers.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

namespace Rythmos {


TEUCHOS_UNIT_TEST( Rythmos_ConvergenceTestHelpers, create ) {
  LinearRegression<double> lr;
  TEST_THROW(lr.getSlope(), std::logic_error);
  TEST_THROW(lr.getYIntercept(), std::logic_error);
}
TEUCHOS_UNIT_TEST( Rythmos_ConvergenceTestHelpers, invalidData ) {
  LinearRegression<double> lr;
  Array<double> x,y;
  TEST_THROW(lr.setData(x,y), std::logic_error);
  x.push_back(0.0);
  y.push_back(0.0);
  TEST_THROW(lr.setData(x,y), std::logic_error);
  x.push_back(1.0);
  TEST_THROW(lr.setData(x,y), std::logic_error);
}
TEUCHOS_UNIT_TEST( Rythmos_ConvergenceTestHelpers, trivialData ) {
  LinearRegression<double> lr;
  Array<double> x,y;
  x.push_back(0.0);
  x.push_back(1.0);
  y.push_back(0.0);
  y.push_back(1.0);
  lr.setData(x,y);
  TEST_EQUALITY_CONST( lr.getSlope(), 1.0 );
  TEST_EQUALITY_CONST( lr.getYIntercept(), 0.0 );
}
TEUCHOS_UNIT_TEST( Rythmos_ConvergenceTestHelpers, lessTrivialData ) {
  LinearRegression<double> lr;
  Array<double> x,y;
  x.push_back(0.0);
  x.push_back(1.0);
  x.push_back(2.0);
  y.push_back(0.0);
  y.push_back(1.0);
  y.push_back(2.0);
  lr.setData(x,y);
  TEST_EQUALITY_CONST( lr.getSlope(), 1.0 );
  TEST_EQUALITY_CONST( lr.getYIntercept(), 0.0 );
}
TEUCHOS_UNIT_TEST( Rythmos_ConvergenceTestHelpers, zeroData ) {
  LinearRegression<double> lr;
  Array<double> x,y;
  x.push_back(0.0);
  x.push_back(0.0);
  y.push_back(1.0);
  y.push_back(2.0);
  TEST_THROW(lr.setData(x,y), std::logic_error);
}
TEUCHOS_UNIT_TEST( Rythmos_ConvergenceTestHelpers, nonUniqueXData ) {
  LinearRegression<double> lr;
  Array<double> x,y;
  x.push_back(0.0);
  x.push_back(1.0);
  x.push_back(1.0);
  y.push_back(0.0);
  y.push_back(0.5);
  y.push_back(1.5);
  lr.setData(x,y);
  double tol = 1.0e-10;
  TEST_FLOATING_EQUALITY( lr.getSlope(), 1.0, tol );
  TEST_FLOATING_EQUALITY( lr.getYIntercept(), 0.0, tol );
}
TEUCHOS_UNIT_TEST( Rythmos_ConvergenceTestHelpers, slightlyOffData ) {
  LinearRegression<double> lr;
  Array<double> x,y;
  x.push_back(0.0);
  x.push_back(1.0);
  x.push_back(2.0);
  y.push_back(0.0);
  y.push_back(0.5);
  y.push_back(2.5);
  lr.setData(x,y);
  TEST_EQUALITY_CONST( lr.getSlope(), 1.25 );
  TEST_EQUALITY_CONST( lr.getYIntercept(), -0.25 );
}
TEUCHOS_UNIT_TEST( Rythmos_ConvergenceTestHelpers, SineData ) {
  LinearRegression<double> lr;
  Array<double> x,y;
  int N = 11;
  double dt = 0.1;
  for (int i=0 ; i<N ; ++i) {
    double xval = dt*i;
    x.push_back( xval );
    y.push_back( sin(xval) );
  }
  lr.setData(x,y);
  // 1.0e-14 works on rancilio but not on gabriel, exetazo, or s858352
  double tol = 1.0e-13;
  TEST_FLOATING_EQUALITY( lr.getSlope(), 8.518189335013251e-01, tol ); // These came from matlab
  TEST_FLOATING_EQUALITY( lr.getYIntercept(), 2.989789515694744e-02, tol );
}
TEUCHOS_UNIT_TEST( Rythmos_ConvergenceTestHelpers, CoSineData ) {
  LinearRegression<double> lr;
  Array<double> x,y;
  int N = 11;
  double dt = 0.1;
  for (int i=0 ; i<N ; ++i) {
    double xval = dt*i;
    x.push_back( xval );
    y.push_back( cos(xval) );
  }
  lr.setData(x,y);
  // 1.0e-14 works on rancilio but not on gabriel, exetazo, or s858352
  double tol = 1.0e-13;
  TEST_FLOATING_EQUALITY( lr.getSlope(), -4.653508042678562e-01, tol ); // These came from matlab
  TEST_FLOATING_EQUALITY( lr.getYIntercept(), 1.067025181571952, tol );
}




} // namespace Rythmos


