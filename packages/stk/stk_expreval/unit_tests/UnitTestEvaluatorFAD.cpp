#include <iostream>
#include <iomanip>
#include <math.h>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_expreval/ExprEvalFAD.hpp>

class UnitTestEvaluatorFAD
{
public:
  void testEvaluator();
};

using namespace stk::expreval;

//  expr_eval.bind("x", x);				\								x
//      expr_eval.setValue("x", x);			\								x
//      std::cout << std::setprecision(20) << x << "," << std::setprecision(20) << y << std::endl;		\	x
namespace {

STKUNIT_UNIT_TEST( UnitTestEvaluatorFAD, testEvaluator)
{
  UnitTestEvaluatorFAD unit;

  unit.testEvaluator();
}



bool
syntax(
  const char *	expr)
{
  std::cout << "Syntax " << expr << " ... ";
  try {
    stk::expreval::fad::Eval expr_eval(stk::expreval::fad::Eval::VariableMap::getDefaultResolver(), expr);
    expr_eval.parse();
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
    stk::expreval::fad::Eval expr_eval(stk::expreval::fad::Eval::VariableMap::getDefaultResolver(), expr);
    expr_eval.parse();
  }
  catch (std::runtime_error &x) {
    std::cout << "pass" << std::endl;
    return true;
  }
  std::cout << "fail, should have parse error" << std::endl;
  return false;
}

/*
bool
vectest(
  const char *	expr)
{
  std::cout << " syntax " << expr << " ...  ";
  try {
    stk::expreval::fad::Eval expr_eval(stk::expreval::fad::Eval::VariableMap::getDefaultResolver(), expr);
    expr_eval.parse();
  }
  catch (std::runtime_error &x) {
    std::cout << "pass" << std::endl;
    return true;
  }
  std::cout << "fail, should have parse error" << std::endl;
  return false;
}
*/

typedef stk::expreval::fad::FADDouble FADDouble;
typedef FADDouble (TestFunc)(stk::expreval::fad::FADDouble);

bool
test(
  const char *	expr,
  TestFunc	c_expr)
{
  bool failed = false;
  std::cout << "Evaluate " << expr << " ... ";
  std::string by_expr = std::string("by=") + expr + ";";
  stk::expreval::fad::Eval expr_eval(stk::expreval::fad::Eval::VariableMap::getDefaultResolver(), by_expr.c_str());
  expr_eval.parse();

  // result and by are the ExprEvalFAD results
  // y is the c++ result
  
  FADDouble x, y, by, result;
  FADDouble v[2];

  expr_eval.bindVariable("x", x);
  expr_eval.bindVariable("by", by);
  expr_eval.bindVariable("v", *v);
  for (int i = 1; i < 100; ++i) {
    x = v[1] = FADDouble(1, 0, i*0.01);
    
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
    FADDouble absolute_error = fabs(result - y);
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
    if (absolute_error.dx(0) > fabs(1.0e-14*result.dx(0))) {
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

#define EXPREVALFAD_DEFINE_TEST(name,expr1,expr2)       \
FADDouble name(FADDouble x) {return expr2;}		\
const char *name##_expr = #expr1

#define EXPREVALFAD_DEFINE_TEST_1(name,expr) EXPREVALFAD_DEFINE_TEST(name,expr,expr)

  // Hiearchy tests
EXPREVALFAD_DEFINE_TEST_1(h1, x*1.0/2.0*3.0);
EXPREVALFAD_DEFINE_TEST_1(h2, x*1.0/2.0*3.0);
EXPREVALFAD_DEFINE_TEST_1(h3, x*(4.0+5.0)/6.0);
// EXPREVALFAD_DEFINE_TEST_1(h4, x==0.5);
// EXPREVALFAD_DEFINE_TEST_1(h5, x>=0.5);
// EXPREVALFAD_DEFINE_TEST_1(h6, x<0.25 || x>0.75);
// EXPREVALFAD_DEFINE_TEST_1(h7, x>0.25 && x<0.75);
// EXPREVALFAD_DEFINE_TEST_1(h8, x*2>0.25 && x<0.75);
// EXPREVALFAD_DEFINE_TEST_1(h9, !(x - 0.5));
// EXPREVALFAD_DEFINE_TEST_1(h10, x > 0.5 ? 5.0 : 7.0);
EXPREVALFAD_DEFINE_TEST_1(h11, x*(x+1.0));
EXPREVALFAD_DEFINE_TEST_1(h12, x*1.0+2.0);
EXPREVALFAD_DEFINE_TEST_1(h13, x*x+1.0);
EXPREVALFAD_DEFINE_TEST_1(h14, x+x*1.0+2.0);
EXPREVALFAD_DEFINE_TEST(h17, v[1]=x*0.5;y=v[1],x*0.5);
EXPREVALFAD_DEFINE_TEST_1(h18, x - -7);
EXPREVALFAD_DEFINE_TEST_1(h19, x - -x);
EXPREVALFAD_DEFINE_TEST_1(h20, x - - - 7);

  // Function tests
// EXPREVALFAD_DEFINE_TEST(f1, abs(x), fabs(x));
// EXPREVALFAD_DEFINE_TEST(f2, mod(x,10.0),fmod(x,10.0));
// EXPREVALFAD_DEFINE_TEST_1(f3, fabs(x));
// EXPREVALFAD_DEFINE_TEST_1(f4, fmod(x,10.0));
EXPREVALFAD_DEFINE_TEST_1(f5, acos(x));
EXPREVALFAD_DEFINE_TEST_1(f6, asin(x));
EXPREVALFAD_DEFINE_TEST_1(f7, atan(x));
EXPREVALFAD_DEFINE_TEST_1(f8, cos(x));
EXPREVALFAD_DEFINE_TEST_1(f9, cosh(x));
EXPREVALFAD_DEFINE_TEST_1(f10, exp(x));
// EXPREVALFAD_DEFINE_TEST_1(f11, ceil(x));
// EXPREVALFAD_DEFINE_TEST_1(f12, floor(x));
EXPREVALFAD_DEFINE_TEST_1(f13, log(x));
EXPREVALFAD_DEFINE_TEST_1(f14, pow(x, 10.0));
// EXPREVALFAD_DEFINE_TEST(f15, pow10(x),pow(x,10.0));
EXPREVALFAD_DEFINE_TEST_1(f16, sin(x));
EXPREVALFAD_DEFINE_TEST_1(f17, sinh(x));
EXPREVALFAD_DEFINE_TEST_1(f18, sqrt(x));
EXPREVALFAD_DEFINE_TEST_1(f19, tan(x));
EXPREVALFAD_DEFINE_TEST_1(f20, tanh(x));
// EXPREVALFAD_DEFINE_TEST(f21, atan2(x, PI),atan2(x, stk::expreval::fad::Eval::s_pi));
EXPREVALFAD_DEFINE_TEST(f22, ln(x), log(x));
// EXPREVALFAD_DEFINE_TEST(f23, deg(x),(180.0 / Expr::Eval::s_pi) * x);
// EXPREVALFAD_DEFINE_TEST(f24, rad(x),(Expr::Eval::s_pi / 180.0) * x);
// EXPREVALFAD_DEFINE_TEST(f25, max(x,1.0),std::max(x,1.0));
// EXPREVALFAD_DEFINE_TEST(f26, min(x,1.0),std::min(x,1.0));
// EXPREVALFAD_DEFINE_TEST(f27, recttopolr(x,1.0),sqrt(x*x+1.0*1.0));
// EXPREVALFAD_DEFINE_TEST(f28, recttopola(x,1.0),atan2(1.0, x));
// EXPREVALFAD_DEFINE_TEST(f29, poltorectx(x,PI/4.0),x*cos(Expr::Eval::s_pi/4.0));
// EXPREVALFAD_DEFINE_TEST(f30, poltorecty(x,PI/4.0),x*sin(Expr::Eval::s_pi/4.0));
EXPREVALFAD_DEFINE_TEST_1(f31, 0.4209+4.5e-4*x);

// Bova tests
EXPREVALFAD_DEFINE_TEST_1(b1, sin(x*.5));
EXPREVALFAD_DEFINE_TEST_1(b2, .5*.2*sin(.5*x));
EXPREVALFAD_DEFINE_TEST_1(b3, .5*sin(x));

} // namespace <unnamed>

void
UnitTestEvaluatorFAD::testEvaluator()
{
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

#define EXPREVALFAD_TEST(name) test(name##_expr, name)

  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h1));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h2));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h3));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h4));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h5));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h6));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h7));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h8));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h9));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h10));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h11));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h12));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h13));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h14));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h16));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h17));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h18));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h19));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(h20));

//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f1));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f2));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f3));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f4));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f5));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f6));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f7));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f8));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f9));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f10));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f11));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f12));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f13));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f14));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f15));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f16));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f17));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f18));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f19));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f20));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f21));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f22));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f23));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f24));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f25));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f26));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f27));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f28));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f29));
//   STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f30));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(f31));

  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(b1));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(b2));
  STKUNIT_EXPECT_TRUE(EXPREVALFAD_TEST(b3));

#undef EXPREVALFAD_TEST
}
