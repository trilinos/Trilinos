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

#ifndef FADUNITTESTS2_HPP
#define FADUNITTESTS2_HPP

// Sacado includes
#include "Sacado.hpp"
#include "Sacado_Random.hpp"

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

// A class for testing each Fad operation based on hand-coded derivatives.
// Only those operations that are defined for both real and complex
// types are tested here.  The others are tested below
template <class FadType, class ScalarType>
class FadOpsUnitTest2 : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE( FadOpsUnitTest2 );
  
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST(testSubtraction);
  CPPUNIT_TEST(testMultiplication);
  CPPUNIT_TEST(testDivision);

  CPPUNIT_TEST(testEquals);
  CPPUNIT_TEST(testNotEquals);

  CPPUNIT_TEST(testPow);

  CPPUNIT_TEST(testUnaryPlus);
  CPPUNIT_TEST(testUnaryMinus);
  
  CPPUNIT_TEST(testExp);
  CPPUNIT_TEST(testLog);
  CPPUNIT_TEST(testLog10);
  CPPUNIT_TEST(testSqrt);
  CPPUNIT_TEST(testCos);
  CPPUNIT_TEST(testSin);
  CPPUNIT_TEST(testTan);
  CPPUNIT_TEST(testCosh);
  CPPUNIT_TEST(testSinh);
  CPPUNIT_TEST(testTanh);

  CPPUNIT_TEST(testPlusEquals);
  CPPUNIT_TEST(testMinusEquals);
  CPPUNIT_TEST(testTimesEquals);
  CPPUNIT_TEST(testDivideEquals);

  CPPUNIT_TEST(testPlusLR);
  CPPUNIT_TEST(testMinusLR);
  CPPUNIT_TEST(testTimesLR);
  CPPUNIT_TEST(testDivideLR);

  CPPUNIT_TEST_SUITE_END();

public:

  FadOpsUnitTest2();

  FadOpsUnitTest2(int numComponents, double absolute_tolerance, 
		 double relative_tolerance);

  void setUp();

  void tearDown();

  void testAddition();
  void testSubtraction();
  void testMultiplication ();
  void testDivision();

  void testEquals();
  void testNotEquals();

  void testPow();

  void testUnaryPlus();
  void testUnaryMinus();

  void testExp();
  void testLog();
  void testLog10();
  void testSqrt();
  void testCos();
  void testSin();
  void testTan();
  void testCosh();
  void testSinh();
  void testTanh();

  void testPlusEquals();
  void testMinusEquals();
  void testTimesEquals();
  void testDivideEquals();

  void testPlusLR();
  void testMinusLR();
  void testTimesLR();
  void testDivideLR();

protected:

  // DFad variables
  FadType a_fad, b_fad, c_fad;

  // Random number generator
  Sacado::Random<ScalarType> urand;

  // Number of derivative components
  int n;

  // Tolerances to which fad objects should be the same
  double tol_a, tol_r;

}; // class FadOpsUnitTest2

template <class FadType, class ScalarType>
FadOpsUnitTest2<FadType,ScalarType>::
FadOpsUnitTest2() :
  urand(), n(5), tol_a(1.0e-15), tol_r(1.0e-14) {}

template <class FadType, class ScalarType>
FadOpsUnitTest2<FadType,ScalarType>::
FadOpsUnitTest2(int numComponents, double absolute_tolerance, 
	       double relative_tolerance) :
  urand(), 
  n(numComponents), 
  tol_a(absolute_tolerance), 
  tol_r(relative_tolerance) {}

template <class FadType, class ScalarType>
void FadOpsUnitTest2<FadType,ScalarType>::setUp() {
  ScalarType val;

  val = urand.number();
  a_fad = FadType(n,val);
  
  val = urand.number();
  b_fad = FadType(n,val);

  for (int i=0; i<n; i++) {
    val = urand.number();
    a_fad.fastAccessDx(i) = val;

    val = urand.number();
    b_fad.fastAccessDx(i) = val;
  }
}

template <class FadType, class ScalarType>
void FadOpsUnitTest2<FadType,ScalarType>::
tearDown() {}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testAddition() {
  c_fad = a_fad + b_fad;
  FadType t1(n, a_fad.val()+b_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i) + b_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
  
  ScalarType val = urand.number();
  c_fad = a_fad + val;
  FadType t2(n, a_fad.val()+val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = a_fad.dx(i);
  COMPARE_FADS(c_fad, t2);

  c_fad = val + b_fad;
  FadType t3(n, val+b_fad.val());
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = b_fad.dx(i);
  COMPARE_FADS(c_fad, t3);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testSubtraction() {
  c_fad = a_fad - b_fad;
  FadType t1(n, a_fad.val()-b_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i) - b_fad.dx(i);
  COMPARE_FADS(c_fad, t1);

  ScalarType val = urand.number();
  c_fad = a_fad - val;
  FadType t2(n, a_fad.val()-val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = a_fad.dx(i);
  COMPARE_FADS(c_fad, t2);

  c_fad = val - b_fad;
  FadType t3(n, val-b_fad.val());
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = -b_fad.dx(i);
  COMPARE_FADS(c_fad, t3);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testMultiplication() {
  c_fad = a_fad * b_fad;
  FadType t1(n, a_fad.val()*b_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)*b_fad.val() + a_fad.val()*b_fad.dx(i);
  COMPARE_FADS(c_fad, t1);

  ScalarType val = urand.number();
  c_fad = a_fad * val;
  FadType t2(n, a_fad.val()*val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = a_fad.dx(i)*val;
  COMPARE_FADS(c_fad, t2);

  c_fad = val * b_fad;
  FadType t3(n, val*b_fad.val());
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = val*b_fad.dx(i);
  COMPARE_FADS(c_fad, t3);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testDivision() {
  c_fad = a_fad / b_fad;
  FadType t1(n, a_fad.val()/b_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = 
      (a_fad.dx(i)*b_fad.val() - a_fad.val()*b_fad.dx(i)) / 
      (b_fad.val()*b_fad.val());
  COMPARE_FADS(c_fad, t1);

  ScalarType val = urand.number();
  c_fad = a_fad / val;
  FadType t2(n, a_fad.val()/val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = a_fad.dx(i)/val;
  COMPARE_FADS(c_fad, t2);

  c_fad = val / b_fad;
  FadType t3(n, val/b_fad.val());
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = -val*b_fad.dx(i)/(b_fad.val()*b_fad.val());
  COMPARE_FADS(c_fad, t3);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testEquals() {
  bool r1 = a_fad == b_fad;
  bool r2 = a_fad.val() == b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);

  ScalarType val = urand.number();
  r1 = a_fad == val;
  r2 = a_fad.val() == val;
  CPPUNIT_ASSERT(r1 == r2);

  r1 = val == b_fad;
  r2 = val == b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testNotEquals() {
  bool r1 = a_fad != b_fad;
  bool r2 = a_fad.val() != b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);

  ScalarType val = urand.number();
  r1 = a_fad != val;
  r2 = a_fad.val() != val;
  CPPUNIT_ASSERT(r1 == r2);

  r1 = val != b_fad;
  r2 = val != b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testUnaryPlus() {
  c_fad = +(a_fad);
  FadType t1(n, a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testUnaryMinus() {
  c_fad = -(a_fad);
  FadType t1(n, -a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = -a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testExp() {
  c_fad = std::exp(a_fad);
  FadType t1(n, std::exp(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = std::exp(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testLog() {
  c_fad = std::log(a_fad);
  FadType t1(n, std::log(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/a_fad.val();
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testLog10() {
  c_fad = std::log10(a_fad);
  FadType t1(n, std::log10(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/(a_fad.val()*std::log(10));
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testSqrt() {
  c_fad = std::sqrt(a_fad);
  FadType t1(n, std::sqrt(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/(2.*std::sqrt(a_fad.val()));
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testCos() {
  c_fad = std::cos(a_fad);
  FadType t1(n, std::cos(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = -std::sin(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testSin() {
  c_fad = std::sin(a_fad);
  FadType t1(n, std::sin(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = std::cos(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testTan() {
  c_fad = std::tan(a_fad);
  FadType t1(n, std::tan(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = 
      a_fad.dx(i)/(std::cos(a_fad.val())*std::cos(a_fad.val()));;
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testCosh() {
  c_fad = std::cosh(a_fad);
  FadType t1(n, std::cosh(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = std::sinh(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testSinh() {
  c_fad = std::sinh(a_fad);
  FadType t1(n, std::sinh(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = std::cosh(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testTanh() {
  c_fad = std::tanh(a_fad);
  FadType t1(n, std::tanh(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = 
      a_fad.dx(i)/(std::cosh(a_fad.val())*std::cosh(a_fad.val()));
  COMPARE_FADS(c_fad, t1);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testPlusEquals() {
  FadType t1(n, c_fad.val()+a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = c_fad.dx(i) + a_fad.dx(i);
  c_fad += a_fad;
  COMPARE_FADS(c_fad, t1);
  
  ScalarType val = urand.number();
  FadType t2(n, c_fad.val()+val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = c_fad.dx(i);
  c_fad += val;
  COMPARE_FADS(c_fad, t2);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testMinusEquals() {
  FadType t1(n, c_fad.val()-a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = c_fad.dx(i) - a_fad.dx(i);
  c_fad -= a_fad;
  COMPARE_FADS(c_fad, t1);
  
  ScalarType val = urand.number();
  FadType t2(n, c_fad.val()-val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = c_fad.dx(i);
  c_fad -= val;
  COMPARE_FADS(c_fad, t2);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testTimesEquals() {
  FadType t1(n, c_fad.val()*a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = c_fad.dx(i)*a_fad.val() + a_fad.dx(i)*c_fad.val();
  c_fad *= a_fad;
  COMPARE_FADS(c_fad, t1);
  
  ScalarType val = urand.number();
  FadType t2(n, c_fad.val()*val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = c_fad.dx(i)*val;
  c_fad *= val;
  COMPARE_FADS(c_fad, t2);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testDivideEquals() {
  FadType t1(n, c_fad.val()/a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = 
      (a_fad.dx(i)*c_fad.val() - c_fad.dx(i)*a_fad.val()) /
      (a_fad.val()*a_fad.val());
  c_fad /= a_fad;
  COMPARE_FADS(c_fad, t1);
  
  ScalarType val = urand.number();
  FadType t2(n, c_fad.val()/val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = c_fad.dx(i)/val;
  c_fad /= val;
  COMPARE_FADS(c_fad, t2);
}

template <class FadType, class ScalarType>
void
FadOpsUnitTest2<FadType,ScalarType>::
testPow() {
  c_fad = std::pow(a_fad, b_fad);
  FadType t1(n, std::pow(a_fad.val(),b_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = 
      std::pow(a_fad.val(),b_fad.val())*(b_fad.val()*a_fad.dx(i)/a_fad.val() + 
					 std::log(a_fad.val())*b_fad.dx(i));
  COMPARE_FADS(c_fad, t1);
  
  ScalarType val = urand.number();
  c_fad = std::pow(a_fad, val);
  FadType t2(n, std::pow(a_fad.val(), val));
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = 
      std::pow(a_fad.val(), val)*(val*a_fad.dx(i)/a_fad.val());
  COMPARE_FADS(c_fad, t2);

  c_fad = std::pow(val, b_fad);
  FadType t3(n, std::pow(val, b_fad.val()));
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = 
      std::pow(val, b_fad.val())*std::log(val)*b_fad.dx(i);
  COMPARE_FADS(c_fad, t3);
}

template <class FadType, class ScalarType>
void 
FadOpsUnitTest2<FadType,ScalarType>::
testPlusLR() {
  FadType aa_fad = a_fad;
  aa_fad = 1.0;
  aa_fad = aa_fad + b_fad;
  c_fad = 1.0 + b_fad;
  COMPARE_FADS(aa_fad, c_fad);
}

template <class FadType, class ScalarType>
void 
FadOpsUnitTest2<FadType,ScalarType>::
testMinusLR() {
  FadType aa_fad = a_fad;
  aa_fad = 1.0;
  aa_fad = aa_fad - b_fad;
  c_fad = 1.0 - b_fad;
  COMPARE_FADS(aa_fad, c_fad);
}

template <class FadType, class ScalarType>
void 
FadOpsUnitTest2<FadType,ScalarType>::
testTimesLR() {
  FadType aa_fad = a_fad;
  aa_fad = 2.0;
  aa_fad = aa_fad * b_fad;
  c_fad = 2.0 * b_fad;
  COMPARE_FADS(aa_fad, c_fad);
}

template <class FadType, class ScalarType>
void 
FadOpsUnitTest2<FadType,ScalarType>::
testDivideLR() {
  FadType aa_fad = a_fad;
  aa_fad = 2.0;
  aa_fad = aa_fad / b_fad;
  c_fad = 2.0 / b_fad;
  COMPARE_FADS(aa_fad, c_fad);
}

// A class for testing each real Fad operation
// This class tests additional functions that aren't define for complex
// types
template <class FadType, class ScalarType>
class RealFadOpsUnitTest2 : public FadOpsUnitTest2<FadType,ScalarType> {

  CPPUNIT_TEST_SUITE( RealFadOpsUnitTest2 );
  
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
  CPPUNIT_TEST(testATan2);
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
  CPPUNIT_TEST(testACosh);
  CPPUNIT_TEST(testASinh);
  CPPUNIT_TEST(testATanh);
  CPPUNIT_TEST(testAbs);
  CPPUNIT_TEST(testFAbs);

  CPPUNIT_TEST(testPlusEquals);
  CPPUNIT_TEST(testMinusEquals);
  CPPUNIT_TEST(testTimesEquals);
  CPPUNIT_TEST(testDivideEquals);

  CPPUNIT_TEST(testPlusLR);
  CPPUNIT_TEST(testMinusLR);
  CPPUNIT_TEST(testTimesLR);
  CPPUNIT_TEST(testDivideLR);

  CPPUNIT_TEST_SUITE_END();

public:

  RealFadOpsUnitTest2() {}

  RealFadOpsUnitTest2(int numComponents, double absolute_tolerance, 
		      double relative_tolerance) :
    FadOpsUnitTest2<FadType,ScalarType>(numComponents, absolute_tolerance, relative_tolerance) {}

  void testLessThanOrEquals();
  void testGreaterThanOrEquals();
  void testLessThan();
  void testGreaterThan();

  void testACos();
  void testASin();
  void testATan();
  void testACosh();
  void testASinh();
  void testATanh();
  void testAbs();
  void testFAbs();

  void testATan2();
  void testMax();
  void testMin();
};

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testLessThanOrEquals() {
  bool r1 = this->a_fad <= this->b_fad;
  bool r2 = this->a_fad.val() <= this->b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = this->a_fad <= val;
  r2 = this->a_fad.val() <= val;
  CPPUNIT_ASSERT(r1 == r2);

  r1 = val <= this->b_fad;
  r2 = val <= this->b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testGreaterThanOrEquals() {
  bool r1 = this->a_fad >= this->b_fad;
  bool r2 = this->a_fad.val() >= this->b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = this->a_fad >= val;
  r2 = this->a_fad.val() >= val;
  CPPUNIT_ASSERT(r1 == r2);

  r1 = val >= this->b_fad;
  r2 = val >= this->b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testLessThan() {
  bool r1 = this->a_fad < this->b_fad;
  bool r2 = this->a_fad.val() < this->b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = this->a_fad < val;
  r2 = this->a_fad.val() < val;
  CPPUNIT_ASSERT(r1 == r2);

  r1 = val < this->b_fad;
  r2 = val < this->b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testGreaterThan() {
  bool r1 = this->a_fad > this->b_fad;
  bool r2 = this->a_fad.val() > this->b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = this->a_fad > val;
  r2 = this->a_fad.val() > val;
  CPPUNIT_ASSERT(r1 == r2);

  r1 = val > this->b_fad;
  r2 = val > this->b_fad.val();
  CPPUNIT_ASSERT(r1 == r2);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testACos() {
  this->c_fad = std::acos(this->a_fad);
  FadType t1(this->n, std::acos(this->a_fad.val()));
  for (int i=0; i<this->n; i++)
    t1.fastAccessDx(i) = -this->a_fad.dx(i)/std::sqrt(1.0 - this->a_fad.val()*this->a_fad.val());
  COMPARE_FADS(this->c_fad, t1);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testASin() {
  this->c_fad = std::asin(this->a_fad);
  FadType t1(this->n, std::asin(this->a_fad.val()));
  for (int i=0; i<this->n; i++)
    t1.fastAccessDx(i) = this->a_fad.dx(i)/std::sqrt(1.0 - this->a_fad.val()*this->a_fad.val());
  COMPARE_FADS(this->c_fad, t1);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testATan() {
  this->c_fad = std::atan(this->a_fad);
  FadType t1(this->n, std::atan(this->a_fad.val()));
  for (int i=0; i<this->n; i++)
    t1.fastAccessDx(i) = this->a_fad.dx(i)/(1.0 + this->a_fad.val()*this->a_fad.val());
  COMPARE_FADS(this->c_fad, t1);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testACosh() {
  FadType aa_fad = this->a_fad;
  if (this->a_fad.val() < 1.0)
    aa_fad.val() = 1.0 / this->a_fad.val();
  this->c_fad = std::acosh(aa_fad);
  FadType t1(this->n, std::acosh(aa_fad.val()));
  for (int i=0; i<this->n; i++)
    t1.fastAccessDx(i) = aa_fad.dx(i)/std::sqrt(aa_fad.val()*aa_fad.val()-1.0);
  COMPARE_FADS(this->c_fad, t1);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testASinh() {
  this->c_fad = std::asinh(this->a_fad);
  FadType t1(this->n, std::asinh(this->a_fad.val()));
  for (int i=0; i<this->n; i++)
    t1.fastAccessDx(i) = this->a_fad.dx(i)/std::sqrt(this->a_fad.val()*this->a_fad.val()+1.0);
  COMPARE_FADS(this->c_fad, t1);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testATanh() {
  this->c_fad = std::atanh(this->a_fad);
  FadType t1(this->n, std::atanh(this->a_fad.val()));
  for (int i=0; i<this->n; i++)
    t1.fastAccessDx(i) = this->a_fad.dx(i)/(1.0 - this->a_fad.val()*this->a_fad.val());
  COMPARE_FADS(this->c_fad, t1);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testAbs() {
  this->c_fad = std::abs(this->a_fad);
  FadType t1(this->n, std::abs(this->a_fad.val()));
  for (int i=0; i<this->n; i++) {
    if (this->a_fad.val() >= 0)
      t1.fastAccessDx(i) = this->a_fad.dx(i);
    else
      t1.fastAccessDx(i) = -this->a_fad.dx(i);
  }
  COMPARE_FADS(this->c_fad, t1);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testFAbs() {
  this->c_fad = std::fabs(this->a_fad);
  FadType t1(this->n, std::fabs(this->a_fad.val()));
  for (int i=0; i<this->n; i++) {
    if (this->a_fad.val() >= 0)
      t1.fastAccessDx(i) = this->a_fad.dx(i);
    else
      t1.fastAccessDx(i) = -this->a_fad.dx(i);
  }
  COMPARE_FADS(this->c_fad, t1);
}

template <class FadType, class ScalarType>
void
RealFadOpsUnitTest2<FadType,ScalarType>::
testATan2() {
  this->c_fad = std::atan2(this->a_fad, this->b_fad);
  FadType t1(this->n, std::atan2(this->a_fad.val(),this->b_fad.val()));
  ScalarType t = this->a_fad.val()*this->a_fad.val() + 
    this->b_fad.val()*this->b_fad.val();
  for (int i=0; i<this->n; i++)
    t1.fastAccessDx(i) = (this->b_fad.val()*this->a_fad.dx(i) - 
			  this->a_fad.val()*this->b_fad.dx(i))/t; 
  COMPARE_FADS(this->c_fad, t1);
  
  ScalarType val = this->urand.number();
  this->c_fad = std::atan2(this->a_fad, val);
  FadType t2(this->n, std::atan2(this->a_fad.val(), val));
  t = this->a_fad.val()*this->a_fad.val() + val*val;
  for (int i=0; i<this->n; i++)
    t2.fastAccessDx(i) = val*this->a_fad.dx(i)/t;
  COMPARE_FADS(this->c_fad, t2);

  this->c_fad = std::atan2(val, this->b_fad);
  FadType t3(this->n, std::atan2(val, this->b_fad.val()));
  t = val*val + this->b_fad.val()*this->b_fad.val();
  for (int i=0; i<this->n; i++)
    t3.fastAccessDx(i) = -val*this->b_fad.dx(i)/t;
  COMPARE_FADS(this->c_fad, t3);
}

template <class FadType, class ScalarType>
void 
RealFadOpsUnitTest2<FadType,ScalarType>::
testMax() {
  ScalarType val;

  // Fad, Fad
  FadType aa_fad = this->a_fad + 1.0;
  this->c_fad = max(aa_fad, this->a_fad);
  COMPARE_FADS(this->c_fad, aa_fad);
  this->c_fad = max(this->a_fad, aa_fad);
  COMPARE_FADS(this->c_fad, aa_fad);

  // Expr, Fad
  this->c_fad = max(this->a_fad+1.0, this->a_fad);
  COMPARE_FADS(this->c_fad, aa_fad);
  this->c_fad = max(this->a_fad, this->a_fad+1.0);
  COMPARE_FADS(this->c_fad, aa_fad);

  // Expr, Expr (same)
  this->c_fad = max(this->a_fad+1.0, this->a_fad+1.0);
  COMPARE_FADS(this->c_fad, aa_fad);

  // Expr, Expr (different)
  this->c_fad = max(this->a_fad+1.0, this->a_fad-1.0);
  COMPARE_FADS(this->c_fad, aa_fad);
  this->c_fad = max(this->a_fad-1.0, this->a_fad+1.0);
  COMPARE_FADS(this->c_fad, aa_fad);
  
  // Fad, const
  val = this->a_fad.val() + 1;
  this->c_fad = max(this->a_fad, val);
  COMPARE_VALUES(this->c_fad.val(), val);
  for (int i=0; i<this->n; i++)
    COMPARE_VALUES(this->c_fad.dx(i), 0.0);
  val = this->a_fad.val() - 1;
  this->c_fad = max(this->a_fad, val);
  COMPARE_FADS(this->c_fad, this->a_fad);
  val = this->b_fad.val() + 1;
  this->c_fad = max(val, this->b_fad);
  COMPARE_VALUES(this->c_fad.val(), val);
  for (int i=0; i<this->n; i++)
    COMPARE_VALUES(this->c_fad.dx(i), 0.0);
  val = this->b_fad.val() - 1;
  this->c_fad = max(val, this->b_fad);
  COMPARE_FADS(this->c_fad, this->b_fad);

  // Expr, const
  val = this->a_fad.val();
  this->c_fad = max(this->a_fad+1.0, val);
  COMPARE_FADS(this->c_fad, aa_fad);
  this->c_fad = max(val, this->a_fad+1.0);
  COMPARE_FADS(this->c_fad, aa_fad);
}

template <class FadType, class ScalarType>
void 
RealFadOpsUnitTest2<FadType,ScalarType>::
testMin() {
  ScalarType val;

  // Fad, Fad
  FadType aa_fad = this->a_fad - 1.0;
  this->c_fad = min(aa_fad, this->a_fad);
  COMPARE_FADS(this->c_fad, aa_fad);
  this->c_fad = min(this->a_fad, aa_fad);
  COMPARE_FADS(this->c_fad, aa_fad);

  // Expr, Fad
  this->c_fad = min(this->a_fad-1.0, this->a_fad);
  COMPARE_FADS(this->c_fad, aa_fad);
  this->c_fad = min(this->a_fad, this->a_fad-1.0);
  COMPARE_FADS(this->c_fad, aa_fad);

  // Expr, Expr (same)
  this->c_fad = min(this->a_fad-1.0, this->a_fad-1.0);
  COMPARE_FADS(this->c_fad, aa_fad);

  // Expr, Expr (different)
  this->c_fad = min(this->a_fad+1.0, this->a_fad-1.0);
  COMPARE_FADS(this->c_fad, aa_fad);
  this->c_fad = min(this->a_fad-1.0, this->a_fad+1.0);
  COMPARE_FADS(this->c_fad, aa_fad);

  // Fad, const
  val = this->a_fad.val() - 1;
  this->c_fad = min(this->a_fad, val);
  COMPARE_VALUES(this->c_fad.val(), val);
  for (int i=0; i<this->n; i++)
    COMPARE_VALUES(this->c_fad.dx(i), 0.0);
  val = this->a_fad.val() + 1;
  this->c_fad = min(this->a_fad, val);
  COMPARE_FADS(this->c_fad, this->a_fad);
  val = this->b_fad.val() - 1;
  this->c_fad = min(val, this->b_fad);
  COMPARE_VALUES(this->c_fad.val(), val);
  for (int i=0; i<this->n; i++)
    COMPARE_VALUES(this->c_fad.dx(i), 0.0);
  val = this->b_fad.val() + 1;
  this->c_fad = min(val, this->b_fad);
  COMPARE_FADS(this->c_fad, this->b_fad);

  // Expr, const
  val = this->a_fad.val();
  this->c_fad = min(this->a_fad-1.0, val);
  COMPARE_FADS(this->c_fad, aa_fad);
  this->c_fad = min(val, this->a_fad-1.0);
  COMPARE_FADS(this->c_fad, aa_fad);
}

#undef COMPARE_VALUES
#undef COMPARE_FADS

#endif // FADUNITTESTS2_HPP
