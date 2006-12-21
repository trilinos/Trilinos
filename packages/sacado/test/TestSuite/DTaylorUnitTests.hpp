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

#ifndef DTAYLORUNITTESTS_HPP
#define DTAYLORUNITTESTS_HPP

// Sacado includes
#include "Sacado.hpp"
#include "Sacado_Random.hpp"

typedef Sacado::Taylor::DTaylor<double> DTaylorType;

// ADOL-C includes
#include "adouble.h"
#include "interfaces.h"

// Cppunit includes
#include <cppunit/extensions/HelperMacros.h>

#define BINARY_OP2_TEST(TESTNAME,OP)	    \
  void TESTNAME () {			    \
    c_dtay = a_dtay OP b_dtay;		    \
    trace_on(0);			    \
    adouble aa, ab, ac;			    \
    aa <<= X[0][0];			    \
    ab <<= X[1][0];			    \
    ac = aa OP ab;			    \
    ac >>= Y[0][0];			    \
    trace_off();			    \
    forward(0,1,2,d,0,X,Y);		    \
    comparePolys(c_dtay,Y[0]);		    \
  }

#define BINARY_OPRC_TEST(TESTNAME,OP)	    \
  void TESTNAME () {			    \
    double val = urand.number();	    \
    c_dtay = a_dtay OP val;		    \
    trace_on(0);			    \
    adouble aa, ac;			    \
    aa <<= X[0][0];			    \
    ac = aa OP val;			    \
    ac >>= Y[0][0];			    \
    trace_off();			    \
    forward(0,1,1,d,0,X,Y);		    \
    comparePolys(c_dtay,Y[0]);		    \
  }

#define BINARY_OPLC_TEST(TESTNAME,OP)	    \
  void TESTNAME () {			    \
    double val = urand.number();	    \
    c_dtay = val OP a_dtay;		    \
    trace_on(0);			    \
    adouble aa, ac;			    \
    aa <<= X[0][0];			    \
    ac = val OP aa;			    \
    ac >>= Y[0][0];			    \
    trace_off();			    \
    forward(0,1,1,d,0,X,Y);		    \
    comparePolys(c_dtay,Y[0]);		    \
  }

#define BINARY_OP_TEST(TESTNAME,OP)			\
  BINARY_OP2_TEST(TESTNAME,OP);				\
  BINARY_OPLC_TEST(TESTNAME ## LeftConstant,OP);	\
  BINARY_OPRC_TEST(TESTNAME ## RightConstant,OP)

#define CPPUNIT_BINARY_OP_TEST(TESTNAME)	\
  CPPUNIT_TEST(TESTNAME);			\
  CPPUNIT_TEST(TESTNAME ## LeftConstant);	\
  CPPUNIT_TEST(TESTNAME ## RightConstant)

#define BINARY_FUNC2_TEST(TESTNAME,FUNC)    \
  void TESTNAME () {			    \
    c_dtay = FUNC (a_dtay, b_dtay);	    \
    trace_on(0);			    \
    adouble aa, ab, ac;			    \
    aa <<= X[0][0];			    \
    ab <<= X[1][0];			    \
    ac = FUNC (aa, ab);			    \
    ac >>= Y[0][0];			    \
    trace_off();			    \
    forward(0,1,2,d,0,X,Y);		    \
    comparePolys(c_dtay,Y[0]);		    \
  }

#define BINARY_FUNCRC_TEST(TESTNAME,FUNC)   \
  void TESTNAME () {			    \
    double val = urand.number();	    \
    c_dtay = FUNC (a_dtay, val);	    \
    trace_on(0);			    \
    adouble aa, ac;			    \
    aa <<= X[0][0];			    \
    ac = FUNC (aa, val);		    \
    ac >>= Y[0][0];			    \
    trace_off();			    \
    forward(0,1,1,d,0,X,Y);		    \
    comparePolys(c_dtay,Y[0]);		    \
  }

#define BINARY_FUNCLC_TEST(TESTNAME,FUNC)   \
  void TESTNAME () {			    \
    double val = urand.number();	    \
    c_dtay = FUNC (val, a_dtay);	    \
    trace_on(0);			    \
    adouble aa, ac;			    \
    aa <<= X[0][0];			    \
    ac = FUNC (val, aa);		    \
    ac >>= Y[0][0];			    \
    trace_off();			    \
    forward(0,1,1,d,0,X,Y);		    \
    comparePolys(c_dtay,Y[0]);		    \
  }

#define BINARY_FUNC_TEST(TESTNAME,FUNC)				\
  BINARY_FUNC2_TEST(TESTNAME,FUNC);				\
  BINARY_FUNCLC_TEST(TESTNAME ## LeftConstant,FUNC);		\
  BINARY_FUNCRC_TEST(TESTNAME ## RightConstant,FUNC)

#define CPPUNIT_BINARY_FUNC_TEST(TESTNAME)	\
  CPPUNIT_TEST(TESTNAME);			\
  CPPUNIT_TEST(TESTNAME ## LeftConstant);	\
  CPPUNIT_TEST(TESTNAME ## RightConstant)

#define UNARY_OP_TEST(TESTNAME,OP)		    \
  void TESTNAME () {				    \
    c_dtay = OP a_dtay;				    \
    trace_on(0);				    \
    adouble aa, ac;				    \
    aa <<= X[0][0];				    \
    ac = OP aa;					    \
    ac >>= Y[0][0];				    \
    trace_off();				    \
    forward(0,1,1,d,0,X,Y);			    \
    comparePolys(c_dtay,Y[0]);			    \
  }

#define UNARY_FUNC_TEST(TESTNAME,FUNC)		    \
  void TESTNAME () {				    \
    c_dtay = FUNC (a_dtay);			    \
    trace_on(0);				    \
    adouble aa, ac;				    \
    aa <<= X[0][0];				    \
    ac = FUNC (aa);				    \
    ac >>= Y[0][0];				    \
    trace_off();				    \
    forward(0,1,1,d,0,X,Y);			    \
    comparePolys(c_dtay,Y[0]);			    \
  }

#define UNARY_ASSIGNOP2_TEST(TESTNAME,OP)	    \
  void TESTNAME () {				    \
    c_dtay = a_dtay;				    \
    c_dtay OP b_dtay;				    \
    trace_on(0);				    \
    adouble aa, ab, ac;				    \
    aa <<= X[0][0];				    \
    ab <<= X[1][0];				    \
    ac = aa;					    \
    ac OP ab;					    \
    ac >>= Y[0][0];				    \
    trace_off();				    \
    forward(0,1,2,d,0,X,Y);			    \
    comparePolys(c_dtay,Y[0]);			    \
  }

#define UNARY_ASSIGNOPRC_TEST(TESTNAME,OP)	    \
  void TESTNAME () {				    \
    double val = urand.number();		    \
    c_dtay = a_dtay;				    \
    c_dtay OP val;				    \
    trace_on(0);				    \
    adouble aa, ac;				    \
    aa <<= X[0][0];				    \
    ac = aa;					    \
    ac OP val;					    \
    ac >>= Y[0][0];				    \
    trace_off();				    \
    forward(0,1,1,d,0,X,Y);			    \
    comparePolys(c_dtay,Y[0]);			    \
  }

#define UNARY_ASSIGNOPLC_TEST(TESTNAME,OP)	    \
  void TESTNAME () {				    \
    double val = urand.number();		    \
    c_dtay = val;				    \
    c_dtay OP a_dtay;				    \
    trace_on(0);				    \
    adouble aa, ac;				    \
    aa <<= X[0][0];				    \
    ac = val;					    \
    ac OP aa;					    \
    ac >>= Y[0][0];				    \
    trace_off();				    \
    forward(0,1,1,d,0,X,Y);			    \
    comparePolys(c_dtay,Y[0]);			    \
  }

#define UNARY_ASSIGNOP_TEST(TESTNAME,OP)		\
  UNARY_ASSIGNOP2_TEST(TESTNAME,OP);			\
  UNARY_ASSIGNOPLC_TEST(TESTNAME ## LeftConstant,OP);	\
  UNARY_ASSIGNOPRC_TEST(TESTNAME ## RightConstant,OP)

#define CPPUNIT_UNARY_ASSIGNOP_TEST(TESTNAME)	\
  CPPUNIT_TEST(TESTNAME);			\
  CPPUNIT_TEST(TESTNAME ## LeftConstant);	\
  CPPUNIT_TEST(TESTNAME ## RightConstant)

// A class for testing each DTaylor operation
class DTaylorOpsUnitTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE( DTaylorOpsUnitTest );

  CPPUNIT_BINARY_OP_TEST(testAddition);
  CPPUNIT_BINARY_OP_TEST(testSubtraction);
  CPPUNIT_BINARY_OP_TEST(testMultiplication);
  CPPUNIT_BINARY_OP_TEST(testDivision);

  CPPUNIT_BINARY_FUNC_TEST(testPow);

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
  CPPUNIT_TEST(testFAbs);

  CPPUNIT_UNARY_ASSIGNOP_TEST(testPlusEquals);
  CPPUNIT_UNARY_ASSIGNOP_TEST(testMinusEquals);
  CPPUNIT_UNARY_ASSIGNOP_TEST(testTimesEquals);
  CPPUNIT_UNARY_ASSIGNOP_TEST(testDivideEquals);

  CPPUNIT_TEST(testComposite1);

  CPPUNIT_TEST_SUITE_END();

public:

  DTaylorOpsUnitTest();

  DTaylorOpsUnitTest(unsigned int degree, double absolute_tolerance, 
		     double relative_tolerance);

  ~DTaylorOpsUnitTest();

  void setUp();

  void tearDown();

  // Assert to Fad objects are the same
  void comparePolys(const DTaylorType& x_dtay,
		    double* x_adolc);

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
  UNARY_FUNC_TEST(testFAbs, fabs);

  UNARY_ASSIGNOP_TEST(testPlusEquals, +=);
  UNARY_ASSIGNOP_TEST(testMinusEquals, -=);
  UNARY_ASSIGNOP_TEST(testTimesEquals, *=);
  UNARY_ASSIGNOP_TEST(testDivideEquals, /=);

  template <typename ScalarT>
  ScalarT composite1(const ScalarT& a, const ScalarT& b) {
    ScalarT t1 = 3. * a + sin(b) / log(fabs(a - b * 7.));
    ScalarT t2 = 1.0e3;
    ScalarT t3 = 5.7e4;
    ScalarT t4 = 3.2e5;
    t1 *= cos(a + exp(t1)) / 6. - tan(t1*sqrt(fabs(a * log10(fabs(b)))));
    t1 -= acos((6.+asin(pow(fabs(a),b)/t2))/t3) * asin(pow(fabs(b),2.)*1.0/t4) * atan((b*pow(2.,log(fabs(a))))/(t3*t4));
    t1 /= cosh(b - 0.7) + 7.*sinh(t1 + 0.8)*tanh(9./a) - 9.;
    t1 += pow(fabs(a*4.),b-8.)/cos(a*b*a);
    
  return t1;
}

  void testComposite1() {
    c_dtay = composite1(a_dtay, b_dtay);
    trace_on(0);
    adouble aa, ab, ac;
    aa <<= X[0][0];
    ab <<= X[1][0];
    ac = composite1(aa,ab);
    ac >>= Y[0][0];
    trace_off();
    forward(0,1,2,d,0,X,Y);
    comparePolys(c_dtay,Y[0]);
  }

  void print_poly(double *x);

  void print_diff(const DTaylorType& x_dtay, double* x_adolc);

protected:

  // DTaylor variables
  DTaylorType a_dtay, b_dtay, c_dtay;

  // ADOL-C arrays
  double **X, **Y;

  // Random number generator
  Sacado::Random urand;

  // Degree of polynomials
  unsigned int d;

  // Tolerances to which fad objects should be the same
  double tol_a, tol_r;

}; // class DTaylorOpsUnitTest

#endif // DTAYUNITTESTS_HPP
