// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef CACHEDFADUNITTESTS_HPP
#define CACHEDFADUNITTESTS_HPP

// Sacado includes
#include "Sacado.hpp"
#include "Sacado_Random.hpp"

typedef Sacado::CacheFad::DFad<double> DFadType;

// Fad includes
#include "Fad/fad.h"

// Cppunit includes
#include <cppunit/extensions/HelperMacros.h>

#define BINARY_OP_TEST(TESTNAME,OP) \
  void TESTNAME () {		    \
    c_dfad = a_dfad OP b_dfad;	    \
    c_fad = a_fad OP b_fad;	    \
    compareFads(c_dfad, c_fad);	    \
				    \
    double val = urand.number();    \
    c_dfad = a_dfad OP val;	    \
    c_fad = a_fad OP val;	    \
    compareFads(c_dfad, c_fad);	    \
				    \
    c_dfad = val OP b_dfad;	    \
    c_fad = val OP b_fad;	    \
    compareFads(c_dfad, c_fad);	    \
  }

#define BINARY_FUNC_TEST(TESTNAME,FUNC) \
  void TESTNAME () {			\
    c_dfad = FUNC (a_dfad,b_dfad);	\
    c_fad = FUNC (a_fad,b_fad);		\
    compareFads(c_dfad, c_fad);		\
    					\
    double val = urand.number();	\
    c_dfad = FUNC (a_dfad,val);		\
    c_fad = FUNC (a_fad,val);		\
    compareFads(c_dfad, c_fad);		\
    					\
    c_dfad = FUNC (val,b_dfad);		\
    c_fad = FUNC (val,b_fad);		\
    compareFads(c_dfad, c_fad);		\
  }

#define UNARY_OP_TEST(TESTNAME,OP)	    \
  void TESTNAME () {			    \
    c_dfad = OP a_dfad;			    \
    c_fad = OP a_fad;			    \
    compareFads(c_dfad, c_fad);		    \
  }

#define UNARY_FUNC_TEST(TESTNAME,FUNC)	    \
  void TESTNAME () {			    \
    c_dfad = FUNC (a_dfad);		    \
    c_fad = FUNC (a_fad);		    \
    compareFads(c_dfad, c_fad);		    \
  }

#define UNARY_ASSIGNOP_TEST(TESTNAME,OP)    \
  void TESTNAME () {			    \
    c_dfad OP a_dfad;			    \
    c_fad OP a_fad;			    \
    compareFads(c_dfad, c_fad);		    \
					    \
    double val = urand.number();	    \
    c_dfad OP val;			    \
    c_fad OP val;			    \
    compareFads(c_dfad, c_fad);		    \
  }

// A class for testing each DFad operation
class CacheDFadOpsUnitTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE( CacheDFadOpsUnitTest );
  
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST(testSubtraction);
  CPPUNIT_TEST(testMultiplication);
  CPPUNIT_TEST(testDivision);

  CPPUNIT_TEST(testPow);
  CPPUNIT_TEST(testMax);
  CPPUNIT_TEST(testMin);

  CPPUNIT_TEST(testUnaryPlus);
  CPPUNIT_TEST(testUnaryMinus);
  
  CPPUNIT_TEST(testExp);
  CPPUNIT_TEST(testLog);
  CPPUNIT_TEST(testLog10);
  CPPUNIT_TEST(testSqrt);
  CPPUNIT_TEST(testCos);
  CPPUNIT_TEST(testSin);
  CPPUNIT_TEST(testTan);
  CPPUNIT_TEST(testACos);
  CPPUNIT_TEST(testASin);
  CPPUNIT_TEST(testATan);
  CPPUNIT_TEST(testCosh);
  CPPUNIT_TEST(testSinh);
  CPPUNIT_TEST(testTanh);
  CPPUNIT_TEST(testAbs);
  CPPUNIT_TEST(testFAbs);

  CPPUNIT_TEST(testPlusEquals);
  CPPUNIT_TEST(testMinusEquals);
  CPPUNIT_TEST(testTimesEquals);
  CPPUNIT_TEST(testDivideEquals);

  CPPUNIT_TEST(testComposite1);

  CPPUNIT_TEST_SUITE_END();

public:

  CacheDFadOpsUnitTest();

  CacheDFadOpsUnitTest(int numComponents, double absolute_tolerance, 
		       double relative_tolerance);

  void setUp();

  void tearDown();

  // Assert to Fad objects are the same
  void compareFads(const DFadType& x_dfad,
		   const FAD::Fad<double>& x_fad);

  // Assert to doubles are the same to relative precision
  void compareDoubles(double a, double b);

  BINARY_OP_TEST(testAddition, +);
  BINARY_OP_TEST(testSubtraction, -);
  BINARY_OP_TEST(testMultiplication, *);
  BINARY_OP_TEST(testDivision, /);

  BINARY_FUNC_TEST(testPow, pow);

  UNARY_OP_TEST(testUnaryPlus, +);
  UNARY_OP_TEST(testUnaryMinus, -);

  UNARY_FUNC_TEST(testExp, exp);
  UNARY_FUNC_TEST(testLog, log);
  UNARY_FUNC_TEST(testLog10, log10);
  UNARY_FUNC_TEST(testSqrt, sqrt);
  UNARY_FUNC_TEST(testCos, cos);
  UNARY_FUNC_TEST(testSin, sin);
  UNARY_FUNC_TEST(testTan, tan);
  UNARY_FUNC_TEST(testACos, acos);
  UNARY_FUNC_TEST(testASin, asin);
  UNARY_FUNC_TEST(testATan, atan);
  UNARY_FUNC_TEST(testCosh, cosh);
  UNARY_FUNC_TEST(testSinh, sinh);
  UNARY_FUNC_TEST(testTanh, tanh);
  UNARY_FUNC_TEST(testAbs, abs);
  UNARY_FUNC_TEST(testFAbs, fabs);

  UNARY_ASSIGNOP_TEST(testPlusEquals, +=);
  UNARY_ASSIGNOP_TEST(testMinusEquals, -=);
  UNARY_ASSIGNOP_TEST(testTimesEquals, *=);
  UNARY_ASSIGNOP_TEST(testDivideEquals, /=);

  void testMax();
  void testMin();

  template <typename ScalarT>
  ScalarT composite1(const ScalarT& a, const ScalarT& b) {
    ScalarT t1 = 3. * a + sin(b) / log(fabs(a - b * 7.));
    ScalarT t2 = 1.0e3;
    ScalarT t3 = 5.7e4;
    ScalarT t4 = 3.2e5;
    t1 *= cos(a + exp(t1)) / 6. - tan(t1*sqrt(abs(a * log10(abs(b)))));
    t1 -= acos((6.+asin(pow(fabs(a),b)/t2))/t3) * asin(pow(fabs(b),2.)*1.0/t4) * atan((b*pow(2.,log(abs(a))))/(t3*t4));
    t1 /= cosh(b - 0.7) + 7.*sinh(t1 + 0.8)*tanh(9./a) - 9.;
    t1 += pow(abs(a*4.),b-8.)/cos(a*b*a);
    
  return t1;
}

  void testComposite1() {
    c_dfad = composite1(a_dfad, b_dfad);
    c_fad = composite1(a_fad, b_fad);
    compareFads(c_dfad, c_fad);
  }

protected:

  // DFad variables
  DFadType a_dfad, b_dfad, c_dfad;

  // Fad variables
  FAD::Fad<double> a_fad, b_fad, c_fad;

  // Random number generator
  Sacado::Random urand;

  // Number of derivative components
  int n;

  // Tolerances to which fad objects should be the same
  double tol_a, tol_r;

}; // class CacheDFadOpsUnitTest

#endif // CACHEDFADUNITTESTS_HPP
