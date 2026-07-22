// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Teuchos_Polynomial.hpp"
#include "Teuchos_Array.hpp"

using Teuchos::Polynomial;
using Teuchos::as;
using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;


TEUCHOS_UNIT_TEST( Teuchos_Polynomial, create ) {
  Polynomial<double> P(0,1.0);
  TEST_EQUALITY_CONST( P.degree(), 0 );
}

TEUCHOS_UNIT_TEST( Teuchos_Polynomial, degrees ) {
  for (unsigned int degree=0 ; degree<10 ; ++degree) {
    Polynomial<double> P(degree,1.0);
    TEST_EQUALITY_CONST( P.degree(), degree );
  }
}

TEUCHOS_UNIT_TEST( Teuchos_Polynomial, coeffs ) {
  unsigned int degree = 10;
  Polynomial<double> P(degree,20.0);
  for (unsigned int d=0 ; d <= degree ; ++d) {
    P.setCoefficient(d,d*1.0);
  }
  for (unsigned int d=0 ; d <= degree ; ++d) {
    TEST_EQUALITY_CONST( *(P.getCoefficient(d)), d*1.0 );
  }
}

TEUCHOS_UNIT_TEST( Teuchos_Polynomial, coeffsPtr ) {
  unsigned int degree = 10;
  Polynomial<double> P(degree);
  for (unsigned int d=0 ; d <= degree ; ++d) {
    RCP<double> coeffPtr = rcp(new double(d*1.0));
    P.setCoefficientPtr(d,coeffPtr);
  }
  for (unsigned int d=0 ; d <= degree ; ++d) {
    TEST_EQUALITY_CONST( *(P.getCoefficient(d)), d*1.0 );
  }
}

TEUCHOS_UNIT_TEST( Teuchos_Polynomial, RCPcoeffs ) {
  int degree = 10;
  Polynomial<double> P(degree,20.0);
  for (int d=0 ; d <= degree ; ++d) {
    P.setCoefficient(d,d*1.0);
  }
  RCP<const double> constCoeff = rcp(new double);
  constCoeff = P.getCoefficient(8);
  TEST_EQUALITY_CONST( *constCoeff, 8.0 );

  RCP<double> coeff = rcp(new double);
  coeff = P.getCoefficient(4);
  TEST_EQUALITY_CONST( *coeff, 4.0 );

}

TEUCHOS_UNIT_TEST( Teuchos_Polynomial, evaluate ) {
  int degree = 2;
  Polynomial<double> P(degree,0.0);
  P.setCoefficient(0,-1.0);
  P.setCoefficient(1, 0.0);
  P.setCoefficient(2, 1.0);
  int numTests = 11;
  Array<double> testValues(numTests);
  for (int i=0 ; i<numTests ; ++i) {
    testValues[i] = (i-5);
  }
  Array<double> polyValues(numTests);
  Array<double> polyDotValues(numTests);
  for (int i=0 ; i<numTests ; ++i) {
    polyValues[i] = pow(testValues[i],2.0)-1.0;
    polyDotValues[i] = 2*testValues[i];
  }
  for (int i=0 ; i<numTests  ; ++i ) {
    double x_out;
    double x_dot_out;
    P.evaluate(testValues[i], &x_out, &x_dot_out );
    TEST_EQUALITY( x_out, polyValues[i] );
    TEST_EQUALITY( x_dot_out, polyDotValues[i] );
  }
}

#ifdef TEUCHOS_DEBUG
TEUCHOS_UNIT_TEST( Teuchos_Polynomial, errors ) {
  {
    unsigned int degree = 2;
    const Polynomial<double> constP(degree,20.0);
    RCP<const double> constCoeff = rcp(new double);
    TEST_THROW( constCoeff = constP.getCoefficient(3), std::out_of_range );
  }
  {
    unsigned int degree = 2;
    Polynomial<double> P(degree,20.0);
    RCP<double> coeff = rcp(new double);
    TEST_THROW( coeff = P.getCoefficient(3), std::out_of_range );
  }
  {
    unsigned int degree = 2;
    Polynomial<double> P(degree,20.0);
    unsigned int i = 3;
    const double coeff = 5.0;
    TEST_THROW( P.setCoefficient(i,coeff), std::out_of_range );
  }
  {
    unsigned int degree = 2;
    Polynomial<double> P(degree);
    unsigned int i = 0;
    const double coeff = 5.0;
    TEST_THROW( P.setCoefficient(i,coeff), std::runtime_error );
  }
  {
    unsigned int degree = 2;
    Polynomial<double> P(degree);
    unsigned int i = 3;
    RCP<double> coeff = rcp(new double(5.0));
    TEST_THROW( P.setCoefficientPtr(i,coeff), std::out_of_range );
  }
  {
    unsigned int degree = 2;
    Polynomial<double> P(degree);
    double t = 2.5;
    double x;
    double x_dot;
    TEST_THROW( P.evaluate(t,&x,&x_dot), std::runtime_error );
  }
}
#endif // TEUCHOS_DEBUG

