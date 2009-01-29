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

#ifndef HERMITEUNITTESTS_HPP
#define HERMITEUNITTESTS_HPP

// Sacado includes
#include "Sacado.hpp"
#include "Sacado_Random.hpp"
#include "Sacado_PCE_OrthogPoly.hpp"
#include "Stokhos_HermiteEBasis.hpp"

typedef Stokhos::HermiteEBasis<double> basis_type;
typedef Sacado::PCE::OrthogPoly<double>::expansion_type exp_type;
typedef Sacado::PCE::OrthogPoly<double> pce_type;

// Cppunit includes
#include <cppunit/extensions/HelperMacros.h>

#define BINARY_OP_TEST(TESTNAME,OP) \
  void TESTNAME () {		    \
    cc = ac OP bc;		    \
    c = a OP b;			    \
    comparePCEs(cc, c);		    \
				    \
    cc = ac OP b;		    \
    c = a OP b;			    \
    comparePCEs(cc, c);		    \
				    \
    cc = a OP bc;		    \
    c = a OP b;			    \
    comparePCEs(cc, c);		    \
  }

#define RELOP_TEST(TESTNAME,OP)     \
  void TESTNAME () {		    \
    bool r1 = ac OP bc;		    \
    bool r2 = a OP b;		    \
    CPPUNIT_ASSERT(r1 == r2);	    \
				    \
    r1 = ac OP b;	            \
    r2 = a OP b;	            \
    CPPUNIT_ASSERT(r1 == r2);	    \
				    \
    r1 = a OP bc;	            \
    r2 = a OP b;	            \
    CPPUNIT_ASSERT(r1 == r2);	    \
  }

#define BINARY_FUNC_TEST(TESTNAME,FUNC) \
  void TESTNAME () {			\
    cc = FUNC (ac, bc);			\
    c = FUNC (a, b);			\
    comparePCEs(cc, c);			\
    					\
    cc = FUNC (ac, b);			\
    c = FUNC (a, b);			\
    comparePCEs(cc, c);			\
    					\
    cc = FUNC (a, bc);			\
    c = FUNC (a, b);			\
    comparePCEs(cc, c);			\
  }

#define UNARY_OP_TEST(TESTNAME,OP)	    \
  void TESTNAME () {			    \
    cc = OP ac;				    \
    c = OP a;				    \
    comparePCEs(cc, c);			    \
  }

#define UNARY_FUNC_TEST(TESTNAME,FUNC)	    \
  void TESTNAME () {			    \
    cc = FUNC (ac);			    \
    c = FUNC (a);			    \
    comparePCEs(cc, c);			    \
  }

#define UNARY_ASSIGNOP_TEST(TESTNAME,OP)    \
  void TESTNAME () {			    \
    cc OP ac;				    \
    c OP a;				    \
    comparePCEs(cc, c);			    \
					    \
    cc OP a;				    \
    c OP a;				    \
    comparePCEs(cc, c);			    \
  }

// A class for testing each DFad operation
class HermiteUnitTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE( HermiteUnitTest );
  
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

  HermiteUnitTest();

  HermiteUnitTest(double absolute_tolerance, 
		  double relative_tolerance);

  void setUp();

  void tearDown();

  // Assert to PCE objects are the same
  void comparePCEs(const pce_type& xc, double x);

  // Assert to doubles are the same to relative precision
  void compareDoubles(double a, double b);

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
  UNARY_FUNC_TEST(testAbs, fabs);
  UNARY_FUNC_TEST(testFAbs, fabs);

  UNARY_ASSIGNOP_TEST(testPlusEquals, +=);
  UNARY_ASSIGNOP_TEST(testMinusEquals, -=);
  UNARY_ASSIGNOP_TEST(testTimesEquals, *=);
  UNARY_ASSIGNOP_TEST(testDivideEquals, /=);

  template <typename ScalarT>
  ScalarT composite1(const ScalarT& x, const ScalarT& y) {
    ScalarT t1 = 3. * x + sin(y) / log(fabs(x - y * 7.));
    ScalarT t2 = 1.0e3;
    ScalarT t3 = 5.7e4;
    ScalarT t4 = 3.2e5;
    t1 *= cos(x + exp(t1)) / 6. - tan(t1*sqrt(fabs(x * log10(fabs(y)))));
    t1 -= acos((6.+asin(pow(fabs(x),y)/t2))/t3) * asin(pow(fabs(y),2.)*1.0/t4) * atan((y*pow(2.,log(fabs(x))))/(t3*t4));
    t1 /= cosh(y - 0.7) + 7.*sinh(t1 + 0.8)*tanh(9./x) - 9.;
    t1 += pow(fabs(x*4.),y-8.)/cos(x*y*x);
    
  return t1;
}

  void testComposite1() {
    cc = composite1(ac, bc);
    c = composite1(a, b);
    comparePCEs(cc, c);
  }

protected:

  // PCE variables
  pce_type ac, bc, cc;

  // Scalars
  double a, b, c;

  // Random number generator
  Sacado::Random urand;

  // Tolerances to which expansions should be the same
  double tol_a, tol_r;

}; // class HermiteUnitTest

#endif // HERMITEUNITTESTS_HPP
