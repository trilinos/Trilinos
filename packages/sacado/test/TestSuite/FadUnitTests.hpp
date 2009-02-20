// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#ifndef FADUNITTESTS_HPP
#define FADUNITTESTS_HPP

// Sacado includes
#include "Sacado.hpp"
#include "Sacado_Random.hpp"

// Fad includes
#include "Fad/fad.h"

// Cppunit includes
#include <cppunit/extensions/HelperMacros.h>

#define COMPARE_VALUES(a, b) \
  CPPUNIT_ASSERT( std::abs(a-b) < this->tol_a + this->tol_r*std::abs(a) );

#define COMPARE_FADS(a, b)                              \
CPPUNIT_ASSERT(a.size() == b.size());			\
CPPUNIT_ASSERT(a.hasFastAccess() == b.hasFastAccess()); \
COMPARE_VALUES(a.val(), b.val());			\
for (int i=0; i<a.size(); i++) {			\
  COMPARE_VALUES(a.dx(i), b.dx(i));			\
  COMPARE_VALUES(a.fastAccessDx(i), b.fastAccessDx(i)); \
 }							\
 ;

#define BINARY_OP_TEST(TESTNAME,OP) \
  void TESTNAME () {		    \
    c_dfad = a_dfad OP b_dfad;	    \
    c_fad = a_fad OP b_fad;	    \
    COMPARE_FADS(c_dfad, c_fad);	    \
				    \
    double val = urand.number();    \
    c_dfad = a_dfad OP val;	    \
    c_fad = a_fad OP val;	    \
    COMPARE_FADS(c_dfad, c_fad);	    \
				    \
    c_dfad = val OP b_dfad;	    \
    c_fad = val OP b_fad;	    \
    COMPARE_FADS(c_dfad, c_fad);	    \
  }

#define RELOP_TEST(TESTNAME,OP)     \
  void TESTNAME () {		    \
    bool r1 = a_dfad OP b_dfad;	    \
    bool r2 = a_fad OP b_fad;	    \
    CPPUNIT_ASSERT(r1 == r2);	    \
				    \
    double val = urand.number();    \
    r1 = a_dfad OP val;	            \
    r2 = a_fad OP val;	            \
    CPPUNIT_ASSERT(r1 == r2);	    \
				    \
    r1 = val OP b_dfad;	            \
    r2 = val OP b_fad;	            \
    CPPUNIT_ASSERT(r1 == r2);	    \
  }

#define BINARY_FUNC_TEST(TESTNAME,FUNC) \
  void TESTNAME () {			\
    c_dfad = FUNC (a_dfad,b_dfad);	\
    c_fad = FUNC (a_fad,b_fad);		\
    COMPARE_FADS(c_dfad, c_fad);		\
    					\
    double val = urand.number();	\
    c_dfad = FUNC (a_dfad,val);		\
    c_fad = FUNC (a_fad,val);		\
    COMPARE_FADS(c_dfad, c_fad);		\
    					\
    c_dfad = FUNC (val,b_dfad);		\
    c_fad = FUNC (val,b_fad);		\
    COMPARE_FADS(c_dfad, c_fad);		\
  }

#define UNARY_OP_TEST(TESTNAME,OP)	    \
  void TESTNAME () {			    \
    c_dfad = OP a_dfad;			    \
    c_fad = OP a_fad;			    \
    COMPARE_FADS(c_dfad, c_fad);		    \
  }

#define UNARY_FUNC_TEST(TESTNAME,FUNC)	    \
  void TESTNAME () {			    \
    c_dfad = FUNC (a_dfad);		    \
    c_fad = FUNC (a_fad);		    \
    COMPARE_FADS(c_dfad, c_fad);		    \
  }

#define UNARY_ASSIGNOP_TEST(TESTNAME,OP)    \
  void TESTNAME () {			    \
    c_dfad OP a_dfad;			    \
    c_fad OP a_fad;			    \
    COMPARE_FADS(c_dfad, c_fad);		    \
					    \
    double val = urand.number();	    \
    c_dfad OP val;			    \
    c_fad OP val;			    \
    COMPARE_FADS(c_dfad, c_fad);		    \
  }

// A class for testing each Fad operation
template <class FadType, class ScalarType>
class FadOpsUnitTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE( FadOpsUnitTest );
  
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST(testSubtraction);
  CPPUNIT_TEST(testMultiplication);
  CPPUNIT_TEST(testDivision);

  CPPUNIT_TEST(testEquals);
  CPPUNIT_TEST(testNotEquals);
  CPPUNIT_TEST(testLessThanOrEquals);
  CPPUNIT_TEST(testGreaterThanOrEquals);
  CPPUNIT_TEST(testLessThan);
  CPPUNIT_TEST(testGreaterThan);

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

  CPPUNIT_TEST(testPlusLR);
  CPPUNIT_TEST(testMinusLR);
  CPPUNIT_TEST(testTimesLR);
  CPPUNIT_TEST(testDivideLR);

  CPPUNIT_TEST_SUITE_END();

public:

  FadOpsUnitTest();

  FadOpsUnitTest(int numComponents, ScalarType absolute_tolerance, 
		 ScalarType relative_tolerance);

  void setUp();

  void tearDown();

  BINARY_OP_TEST(testAddition, +);
  BINARY_OP_TEST(testSubtraction, -);
  BINARY_OP_TEST(testMultiplication, *);
  BINARY_OP_TEST(testDivision, /);

  RELOP_TEST(testEquals, ==);
  RELOP_TEST(testNotEquals, !=);
  RELOP_TEST(testLessThanOrEquals, <=);
  RELOP_TEST(testGreaterThanOrEquals, >=);
  RELOP_TEST(testLessThan, <);
  RELOP_TEST(testGreaterThan, >);

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
    COMPARE_FADS(c_dfad, c_fad);
  }

  void testPlusLR() {
    FadType aa_dfad = a_dfad;
    FAD::Fad<ScalarType> aa_fad = a_fad;
    aa_dfad = 1.0;
    aa_fad = 1.0;
    aa_dfad = aa_dfad + b_dfad;
    aa_fad = aa_fad + b_fad;
    COMPARE_FADS(aa_dfad, aa_fad);
  }

  void testMinusLR() {
    FadType aa_dfad = a_dfad;
    FAD::Fad<ScalarType> aa_fad = a_fad;
    aa_dfad = 1.0;
    aa_fad = 1.0;
    aa_dfad = aa_dfad - b_dfad;
    aa_fad = aa_fad - b_fad;
    COMPARE_FADS(aa_dfad, aa_fad);
  }

  void testTimesLR() {
    FadType aa_dfad = a_dfad;
    FAD::Fad<ScalarType> aa_fad = a_fad;
    aa_dfad = 2.0;
    aa_fad = 2.0;
    aa_dfad = aa_dfad * b_dfad;
    aa_fad = aa_fad * b_fad;
    COMPARE_FADS(aa_dfad, aa_fad);
  }

  void testDivideLR() {
    FadType aa_dfad = a_dfad;
    FAD::Fad<ScalarType> aa_fad = a_fad;
    aa_dfad = 2.0;
    aa_fad = 2.0;
    aa_dfad = aa_dfad / b_dfad;
    aa_fad = aa_fad / b_fad;
    COMPARE_FADS(aa_dfad, aa_fad);
  }

protected:

  // DFad variables
  FadType a_dfad, b_dfad, c_dfad;

  // Fad variables
  FAD::Fad<ScalarType> a_fad, b_fad, c_fad;

  // Random number generator
  Sacado::Random urand;

  // Number of derivative components
  int n;

  // Tolerances to which fad objects should be the same
  ScalarType tol_a, tol_r;

}; // class FadOpsUnitTest

template <class FadType, class ScalarType>
FadOpsUnitTest<FadType,ScalarType>::
FadOpsUnitTest() :
  urand(0.0, 1.0), n(5), tol_a(1.0e-15), tol_r(1.0e-14) {}

template <class FadType, class ScalarType>
FadOpsUnitTest<FadType,ScalarType>::
FadOpsUnitTest(int numComponents, ScalarType absolute_tolerance, 
	       ScalarType relative_tolerance) :
  urand(0.0, 1.0), 
  n(numComponents), 
  tol_a(absolute_tolerance), 
  tol_r(relative_tolerance) {}

template <class FadType, class ScalarType>
void FadOpsUnitTest<FadType,ScalarType>::setUp() {
  ScalarType val;

  val = urand.number();
  a_dfad = FadType(n,val);
  a_fad = FAD::Fad<ScalarType>(n,val);
  
  val = urand.number();
  b_dfad = FadType(n,val);
  b_fad = FAD::Fad<ScalarType>(n,val);

  for (int i=0; i<n; i++) {
    val = urand.number();
    a_dfad.fastAccessDx(i) = val;
    a_fad.fastAccessDx(i) = val;

    val = urand.number();
    b_dfad.fastAccessDx(i) = val;
    b_fad.fastAccessDx(i) = val;
  }
}

template <class FadType, class ScalarType>
void FadOpsUnitTest<FadType,ScalarType>::
tearDown() {}

template <class FadType, class ScalarType>
void FadOpsUnitTest<FadType,ScalarType>::
testMax() {
  ScalarType val;

  FadType aa_dfad = a_dfad + 1.0;
  c_dfad = max(aa_dfad, a_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }
  
  c_dfad = max(a_dfad, aa_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = max(a_dfad+1.0, a_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }
  
  c_dfad = max(a_dfad, a_dfad+1.0);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }
  
  val = a_dfad.val() + 1;
  c_dfad = max(a_dfad, val);
  COMPARE_VALUES(c_dfad.val(), val);
  for (int i=0; i<n; i++)
    COMPARE_VALUES(c_dfad.dx(i), 0.0);
  
  val = a_dfad.val() - 1;
  c_dfad = max(a_dfad, val);
  COMPARE_VALUES(c_dfad.val(), a_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), a_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), a_dfad.fastAccessDx(i));
  }

  val = b_dfad.val() + 1;
  c_dfad = max(val, b_dfad);
  COMPARE_VALUES(c_dfad.val(), val);
  for (int i=0; i<n; i++)
    COMPARE_VALUES(c_dfad.dx(i), 0.0);
  
  val = b_dfad.val() - 1;
  c_dfad = max(val, b_dfad);
  COMPARE_VALUES(c_dfad.val(), b_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), b_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), b_dfad.fastAccessDx(i));
  }
}

template <class FadType, class ScalarType>
void FadOpsUnitTest<FadType,ScalarType>::
testMin() {
  ScalarType val;

  FadType aa_dfad = a_dfad - 1.0;
  c_dfad = min(aa_dfad, a_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = min(a_dfad, aa_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  val = a_dfad.val() - 1;
  c_dfad = min(a_dfad, val);
  COMPARE_VALUES(c_dfad.val(), val);
  for (int i=0; i<n; i++)
    COMPARE_VALUES(c_dfad.dx(i), 0.0);
  
  val = a_dfad.val() + 1;
  c_dfad = min(a_dfad, val);
  COMPARE_VALUES(c_dfad.val(), a_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), a_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), a_dfad.fastAccessDx(i));
  }

  val = b_dfad.val() - 1;
  c_dfad = min(val, b_dfad);
  COMPARE_VALUES(c_dfad.val(), val);
  for (int i=0; i<n; i++)
    COMPARE_VALUES(c_dfad.dx(i), 0.0);
  
  val = b_dfad.val() + 1;
  c_dfad = min(val, b_dfad);
  COMPARE_VALUES(c_dfad.val(), b_dfad.val());
  for (int i=0; i<n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), b_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), b_dfad.fastAccessDx(i));
  }
}

#undef COMPARE_VALUES(a, b)
#undef COMPARE_FADS(a, b)

#endif // FADUNITTESTS_HPP
