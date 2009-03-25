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

#ifndef LOGICALSPARSEUNITTESTS_HPP
#define LOGICALSPARSEUNITTESTS_HPP

// Sacado includes
#include "Sacado.hpp"
#include "Sacado_Random.hpp"

typedef Sacado::Fad::DFad<double> DFadType;
typedef Sacado::LFad::LogicalSparse<double,bool> LSType;

// Cppunit includes
#include <cppunit/extensions/HelperMacros.h>

#define BINARY_OP_TEST(TESTNAME,OP) \
  void TESTNAME () {		    \
    c_dfad = a_dfad OP b_dfad;	    \
    c_ls = a_ls OP b_ls;	    \
    compareFads(c_dfad, c_ls);	    \
				    \
    double val = urand.number();    \
    c_dfad = a_dfad OP val;	    \
    c_ls = a_ls OP val;		    \
    compareFads(c_dfad, c_ls);	    \
				    \
    c_dfad = val OP b_dfad;	    \
    c_ls = val OP b_ls;		    \
    compareFads(c_dfad, c_ls);	    \
  }

#define RELOP_TEST(TESTNAME,OP)     \
  void TESTNAME () {		    \
    bool r1 = a_dfad OP b_dfad;	    \
    bool r2 = a_ls OP b_ls;	    \
    CPPUNIT_ASSERT(r1 == r2);	    \
				    \
    double val = urand.number();    \
    r1 = a_dfad OP val;	            \
    r2 = a_ls OP val;	            \
    CPPUNIT_ASSERT(r1 == r2);	    \
				    \
    r1 = val OP b_dfad;	            \
    r2 = val OP b_ls;	            \
    CPPUNIT_ASSERT(r1 == r2);	    \
  }

#define BINARY_FUNC_TEST(TESTNAME,FUNC) \
  void TESTNAME () {			\
    c_dfad = FUNC (a_dfad,b_dfad);	\
    c_ls = FUNC (a_ls,b_ls);		\
    compareFads(c_dfad, c_ls);		\
    					\
    double val = urand.number();	\
    c_dfad = FUNC (a_dfad,val);		\
    c_ls = FUNC (a_ls,val);		\
    compareFads(c_dfad, c_ls);		\
    					\
    c_dfad = FUNC (val,b_dfad);		\
    c_ls = FUNC (val,b_ls);		\
    compareFads(c_dfad, c_ls);		\
  }

#define UNARY_OP_TEST(TESTNAME,OP)	    \
  void TESTNAME () {			    \
    c_dfad = OP a_dfad;			    \
    c_ls = OP a_ls;			    \
    compareFads(c_dfad, c_ls);		    \
  }

#define UNARY_FUNC_TEST(TESTNAME,FUNC)	    \
  void TESTNAME () {			    \
    c_dfad = FUNC (a_dfad);		    \
    c_ls = FUNC (a_ls);			    \
    compareFads(c_dfad, c_ls);		    \
  }

#define UNARY_ASSIGNOP_TEST(TESTNAME,OP)    \
  void TESTNAME () {			    \
    c_dfad OP a_dfad;			    \
    c_ls OP a_ls;			    \
    compareFads(c_dfad, c_ls);		    \
					    \
    double val = urand.number();	    \
    c_dfad OP val;			    \
    c_ls OP val;			    \
    compareFads(c_dfad, c_ls);		    \
  }

// A class for testing each DFad operation
class LogicalSparseOpsUnitTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE( LogicalSparseOpsUnitTest );
  
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

  LogicalSparseOpsUnitTest();

  LogicalSparseOpsUnitTest(int numComponents, double absolute_tolerance, 
			   double relative_tolerance);

  void setUp();

  void tearDown();

  // Assert to Fad objects are the same
  void compareFads(const DFadType& x_dfad, const LSType& x_ls);

  // Assert two doubles are the same to relative precision
  void compareDoubles(double a, double b);
  
  // Assert two bools are the same
  void compareBools(bool a, bool b);

  // Assert a double and bool are same (logically)
  void compareDx(double a, bool b);

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
    c_ls = composite1(a_ls, b_ls);
    compareFads(c_dfad, c_ls);
  }

  void testPlusLR() {
    DFadType aa_dfad = a_dfad;
    LSType aa_ls = a_ls;
    aa_dfad = 1.0;
    aa_ls = 1.0;
    aa_dfad = aa_dfad + b_dfad;
    aa_ls = aa_ls + b_ls;
    compareFads(aa_dfad, aa_ls);
  }

  void testMinusLR() {
    DFadType aa_dfad = a_dfad;
    LSType aa_ls = a_ls;
    aa_dfad = 1.0;
    aa_ls = 1.0;
    aa_dfad = aa_dfad - b_dfad;
    aa_ls = aa_ls - b_ls;
    compareFads(aa_dfad, aa_ls);
  }

  void testTimesLR() {
    DFadType aa_dfad = a_dfad;
    LSType aa_ls = a_ls;
    aa_dfad = 2.0;
    aa_ls = 2.0;
    aa_dfad = aa_dfad * b_dfad;
    aa_ls = aa_ls * b_ls;
    compareFads(aa_dfad, aa_ls);
  }

  void testDivideLR() {
    DFadType aa_dfad = a_dfad;
    LSType aa_ls = a_ls;
    aa_dfad = 2.0;
    aa_ls = 2.0;
    aa_dfad = aa_dfad / b_dfad;
    aa_ls = aa_ls / b_ls;
    compareFads(aa_dfad, aa_ls);
  }

protected:

  // DFad variables
  DFadType a_dfad, b_dfad, c_dfad;

  // Logical sparse variables
  LSType a_ls, b_ls, c_ls;

  // Random number generator
  Sacado::Random<double> urand;

  // Number of derivative components
  int n;

  // Tolerances to which fad objects should be the same
  double tol_a, tol_r;

}; // class LogicalSparseOpsUnitTest

#endif // LOGICALSPARSEUNITTESTS_HPP
