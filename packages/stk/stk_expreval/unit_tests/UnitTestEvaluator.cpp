/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <math.h>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_expreval/Evaluator.hpp>

class UnitTestEvaluator
{
public:
  void testEvaluator();
};

using namespace stk::expreval;

//  expr_eval.bind("x", x);				\								x
//      expr_eval.setValue("x", x);			\								x
//      std::cout << std::setprecision(20) << x << "," << std::setprecision(20) << y << std::endl;		\	x
namespace {

STKUNIT_UNIT_TEST( UnitTestEvaluator, testEvaluator)
{
  UnitTestEvaluator unit;

  unit.testEvaluator();
}


bool
syntax(
  const char *	expr)
{
  std::cout << "Syntax " << expr << " ... ";
  try {
    stk::expreval::Eval expr_eval(stk::expreval::VariableMap::getDefaultResolver(), expr);
  }
  catch (std::runtime_error &x) {
    std::cout << "fail, " << x.what() << std::endl;
    return false;
  }
  std::cout << "pass" << std::endl;
  return true;
}


bool
fail_syntax(
  const char *	expr)
{
  std::cout << "Invalid syntax " << expr << " ...  ";
  try {
    stk::expreval::Eval expr_eval(stk::expreval::VariableMap::getDefaultResolver(), expr);
  }
  catch (std::runtime_error &x) {
    std::cout << "pass" << std::endl;
    return true;
  }
  std::cout << "fail, should have parse error" << std::endl;
  return false;
}

bool
vectest(
  const char *	expr)
{
  std::cout << " syntax " << expr << " ...  ";
  try {
    stk::expreval::Eval expr_eval(stk::expreval::VariableMap::getDefaultResolver(), expr);
  }
  catch (std::runtime_error &x) {
    std::cout << "pass" << std::endl;
    return true;
  }
  std::cout << "fail, should have parse error" << std::endl;
  return false;
}

typedef double (TestFunc)(double);

bool
test(
  const char *	expr,
  TestFunc	c_expr)
{
  bool failed = false;
  std::cout << "Evaluate " << expr << " ... ";
  std::string by_expr = std::string("by=") + expr + ";";
  stk::expreval::Eval expr_eval(stk::expreval::VariableMap::getDefaultResolver(), by_expr.c_str());

  double x, y, by, result = 0.0;
  double v[2];

  expr_eval.bindVariable("x", x);
  expr_eval.bindVariable("by", by);
  expr_eval.bindVariable("v", *v);
  for (int i = 1; i < 100; ++i) {
    x = v[1] = i*0.01;
    y = (*c_expr)(x);
    try {
      result = expr_eval.evaluate();
    }
    catch (std::runtime_error &exc) {
      std::cout << expr << " at "
		<< std::setprecision(20) << x << " is "
		<< std::setprecision(20) << result
		<< " should be " << y
		<< "(" << std::setprecision(20) << by
		<< ") threw exception " << exc.what()
		<< std::endl;
      failed = true;
    }
    double absolute_error = fabs(result - y);
    if (absolute_error > fabs(1.0e-14*result)) {
      std::cout << expr << " at "
		<< std::setprecision(2) << x << " is "
		<< std::setprecision(20) << result
		<< " should be " << y
		<< " error is " << absolute_error
		<< std::endl;
      failed = true;
    }
    else if (by != result) {
      std::cout << expr << " at "
		<< std::setprecision(2) << x << " is "
		<< std::setprecision(20) << result
		<< " does not match bound value "
		<< std::setprecision(20) << by
		<< std::endl;
      failed = true;
    }
  }

  std::cout << (failed ? "fail" : "pass") << std::endl;
  return !failed;
}

#define EXPREVAL_DEFINE_TEST(name,expr1,expr2)			\
double name(double x) {return expr2;}			        \
const char *name##_expr = #expr1

#define EXPREVAL_DEFINE_TEST1(name,expr) EXPREVAL_DEFINE_TEST(name,expr,expr)

#define EXPREVAL_TEST(name) test(name##_expr, name)

  // Hiearchy tests
EXPREVAL_DEFINE_TEST1(h1, x*1.0/2.0*3.0);
EXPREVAL_DEFINE_TEST1(h2, x*1.0/2.0*3.0);
EXPREVAL_DEFINE_TEST1(h3, x*(4.0+5.0)/6.0);
EXPREVAL_DEFINE_TEST1(h4, x==0.5);
EXPREVAL_DEFINE_TEST1(h5, x>=0.5);
EXPREVAL_DEFINE_TEST1(h6, x<0.25 || x>0.75);
EXPREVAL_DEFINE_TEST1(h7, x>0.25 && x<0.75);
EXPREVAL_DEFINE_TEST1(h8, x*2>0.25 && x<0.75);
EXPREVAL_DEFINE_TEST1(h9, !(x - 0.5));
EXPREVAL_DEFINE_TEST1(h10, x > 0.5 ? 5.0 : 7.0);
EXPREVAL_DEFINE_TEST1(h11, x*(x+1.0));
EXPREVAL_DEFINE_TEST1(h12, x*1.0+2.0);
EXPREVAL_DEFINE_TEST1(h13, x*x+1.0);
EXPREVAL_DEFINE_TEST1(h14, x+x*1.0+2.0);
EXPREVAL_DEFINE_TEST1(h15, x > 0.5 ? x/2.0 + 1.0 : x*2.0 - 1.0);
EXPREVAL_DEFINE_TEST1(h16, x > 0.5 ? x > 0.75 ? x/2.0 + 1.0 : x*2.0 - 1.0 : x*5.0/2.0);
EXPREVAL_DEFINE_TEST(h17, v[1]=x*0.5;y=v[1],x*0.5);
EXPREVAL_DEFINE_TEST1(h18, x - -7);
EXPREVAL_DEFINE_TEST1(h19, x - -x);
EXPREVAL_DEFINE_TEST1(h20, x - - - 7);

  // Function tests
EXPREVAL_DEFINE_TEST(f1, abs(x), fabs(x));
EXPREVAL_DEFINE_TEST(f2, mod(x,10.0),fmod(x,10.0));
EXPREVAL_DEFINE_TEST1(f3, fabs(x));
EXPREVAL_DEFINE_TEST1(f4, fmod(x,10.0));
EXPREVAL_DEFINE_TEST1(f5, acos(x));
EXPREVAL_DEFINE_TEST1(f6, asin(x));
EXPREVAL_DEFINE_TEST1(f7, atan(x));
EXPREVAL_DEFINE_TEST1(f8, ceil(x));
EXPREVAL_DEFINE_TEST1(f9, cos(x));
EXPREVAL_DEFINE_TEST1(f10, cosh(x));
EXPREVAL_DEFINE_TEST1(f11, exp(x));
EXPREVAL_DEFINE_TEST1(f12, floor(x));
EXPREVAL_DEFINE_TEST1(f13, log(x));
EXPREVAL_DEFINE_TEST1(f14, pow(x, 10.0));
// EXPREVAL_DEFINE_TEST(f15, pow10(x),pow(x,10.0));
EXPREVAL_DEFINE_TEST1(f16, sin(x));
EXPREVAL_DEFINE_TEST1(f17, sinh(x));
EXPREVAL_DEFINE_TEST1(f18, sqrt(x));
EXPREVAL_DEFINE_TEST1(f19, tan(x));
EXPREVAL_DEFINE_TEST1(f20, tanh(x));
EXPREVAL_DEFINE_TEST(f21, atan2(x, PI),atan2(x, stk::expreval::s_pi));
EXPREVAL_DEFINE_TEST(f22, ln(x),log(x));
EXPREVAL_DEFINE_TEST(f23, deg(x),(180.0 / stk::expreval::s_pi) * x);
EXPREVAL_DEFINE_TEST(f24, rad(x),(stk::expreval::s_pi / 180.0) * x);
EXPREVAL_DEFINE_TEST(f25, max(x,1.0),std::max(x,1.0));
EXPREVAL_DEFINE_TEST(f26, min(x,1.0),std::min(x,1.0));
EXPREVAL_DEFINE_TEST(f27, recttopolr(x,1.0),sqrt(x*x+1.0*1.0));
EXPREVAL_DEFINE_TEST(f28, recttopola(x,1.0),atan2(1.0, x));
EXPREVAL_DEFINE_TEST(f29, poltorectx(x,PI/4.0),x*cos(stk::expreval::s_pi/4.0));
EXPREVAL_DEFINE_TEST(f30, poltorecty(x,PI/4.0),x*sin(stk::expreval::s_pi/4.0));
EXPREVAL_DEFINE_TEST1(f31, 0.4209+4.5e-4*x);

// EXPREVAL_DEFINE_TEST2(f1, fpart(x),modf(x,&y));
// EXPREVAL_DEFINE_TEST2(f1, rand(x),rand());
// EXPREVAL_DEFINE_TEST2(f1, random(x), random());
// EXPREVAL_DEFINE_TEST1(f1, randomize(x));
// EXPREVAL_DEFINE_TEST1(f1, srand(x));

// Bova tests
EXPREVAL_DEFINE_TEST1(b1, sin(x*.5));
EXPREVAL_DEFINE_TEST1(b2, .5*.2*sin(.5*x));
EXPREVAL_DEFINE_TEST1(b3, .5*sin(x));

} // namespace <unnamed>


void
UnitTestEvaluator::testEvaluator()
{
  STKUNIT_EXPECT_TRUE(syntax(""));
  STKUNIT_EXPECT_TRUE(syntax(""));
  STKUNIT_EXPECT_TRUE(syntax(";"));
  STKUNIT_EXPECT_TRUE(syntax(";;"));
  STKUNIT_EXPECT_TRUE(syntax(";;;"));
  STKUNIT_EXPECT_TRUE(syntax("x*0.1"));
  STKUNIT_EXPECT_TRUE(syntax("x*-0.1"));
  STKUNIT_EXPECT_TRUE(syntax("x*+0.1"));
  STKUNIT_EXPECT_TRUE(syntax("x--7.0"));
  STKUNIT_EXPECT_TRUE(syntax("x*-x"));
  STKUNIT_EXPECT_TRUE(syntax("x*+x"));
  STKUNIT_EXPECT_TRUE(syntax("v[0]=v[1]*0.1"));
  STKUNIT_EXPECT_TRUE(syntax("x--x"));
  STKUNIT_EXPECT_TRUE(syntax("x---x"));
  STKUNIT_EXPECT_TRUE(fail_syntax("0.01.02"));
  STKUNIT_EXPECT_TRUE(fail_syntax("5*.e+10"));
  STKUNIT_EXPECT_TRUE(fail_syntax("x y"));
  STKUNIT_EXPECT_TRUE(fail_syntax("x(y"));
  STKUNIT_EXPECT_TRUE(fail_syntax("x*"));
  STKUNIT_EXPECT_TRUE(fail_syntax("x*(y+1"));
  STKUNIT_EXPECT_TRUE(fail_syntax("x*y)"));
  STKUNIT_EXPECT_TRUE(fail_syntax("cos(x"));
  STKUNIT_EXPECT_TRUE(fail_syntax("(x)y"));
  STKUNIT_EXPECT_TRUE(fail_syntax("()"));
  STKUNIT_EXPECT_TRUE(syntax("rand()"));

  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h1));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h2));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h3));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h4));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h5));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h6));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h7));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h8));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h9));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h10));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h11));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h12));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h13));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h14));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h15));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h16));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h17));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h18));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h19));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(h20));

  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f1));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f2));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f3));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f4));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f5));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f6));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f7));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f8));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f9));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f10));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f11));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f12));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f13));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f14));
//  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f15));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f16));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f17));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f18));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f19));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f20));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f21));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f22));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f23));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f24));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f25));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f26));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f27));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f28));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f29));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f30));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(f31));

  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(b1));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(b2));
  STKUNIT_EXPECT_TRUE(EXPREVAL_TEST(b3));
}
