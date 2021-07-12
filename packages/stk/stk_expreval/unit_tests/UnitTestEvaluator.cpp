// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <fstream>
#include <iostream>
#include <limits>
#include <iomanip>
#include <cmath>

#include <Kokkos_Core.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gtest/gtest.h>

#include <stk_expreval/Evaluator.hpp>
#include <stk_util/util/ThreadLocalData.hpp>

using ExecutionSpace = Kokkos::Serial;

class UnitTestEvaluator
{
public:
  void testEvaluator();
};

using namespace stk::expreval;

namespace {

TEST( UnitTestEvaluator, testEvaluator)
{
  UnitTestEvaluator unit;

  unit.testEvaluator();
}

std::string generate_expression(std::string expression) {
  return std::string("by=") + expression + ";";
}

double
kernel_evaluate(stk::expreval::Eval& inputEval) 
{
  double result;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(0,1), KOKKOS_LAMBDA (const int& i, double& localResult) {
    localResult = inputEval.evaluate();
  }, result);

  return result;
}

bool
has_variable(const std::vector<std::string>& variableNames, const std::string& variableName) 
{
  return (std::find(variableNames.begin(), variableNames.end(), variableName) != variableNames.end());
}

TEST(UnitTestEvaluator, isConstantExpression_empty)
{
  stk::expreval::Eval eval;
  eval.parse();
  EXPECT_EQ(eval.is_constant_expression(), true);
}

TEST(UnitTestEvaluator, isConstantExpression_constant)
{
  stk::expreval::Eval eval("2");
  eval.parse();
  EXPECT_EQ(eval.is_constant_expression(), true);
}

TEST(UnitTestEvaluator, isConstantExpression_variable)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_EQ(eval.is_constant_expression(), false);
}

TEST(UnitTestEvaluator, isVariable_no)
{
  stk::expreval::Eval eval;
  eval.parse();
  EXPECT_EQ(eval.is_variable("x"), false);
}

TEST(UnitTestEvaluator, isVariable_yes)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_EQ(eval.is_variable("x"), true);
}

TEST(UnitTestEvaluator, isVariable_twoVariables_yes)
{
  stk::expreval::Eval eval("x + y");
  eval.parse();
  EXPECT_EQ(eval.is_variable("x"), true);
  EXPECT_EQ(eval.is_variable("y"), true);
}

TEST(UnitTestEvaluator, isVariable_twoVariables_no)
{
  stk::expreval::Eval eval("x + y");
  eval.parse();
  EXPECT_EQ(eval.is_variable("x"), true);
  EXPECT_EQ(eval.is_variable("z"), false);
}

TEST(UnitTestEvaluator, getAllVariables_noVariables)
{
  stk::expreval::Eval eval;
  eval.parse();
  EXPECT_EQ(eval.get_variable_names().size(), 0u);
}

TEST(UnitTestEvaluator, getAllVariables_noAssign)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_variable_names();
  EXPECT_EQ(variableNames.size(), 1u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
}

TEST(UnitTestEvaluator, getAllVariables_constant)
{
  stk::expreval::Eval eval("x = 2");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_variable_names();
  EXPECT_EQ(variableNames.size(), 1u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
}

TEST(UnitTestEvaluator, getAllVariables_oneDependent)
{
  stk::expreval::Eval eval("x = sin(y)");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_variable_names();
  EXPECT_EQ(variableNames.size(), 2u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
  EXPECT_TRUE(has_variable(variableNames, "y"));
}

TEST(UnitTestEvaluator, getAllVariables_constantAssign)
{
  stk::expreval::Eval eval("x = 2; y = x");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_variable_names();
  EXPECT_EQ(variableNames.size(), 2u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
  EXPECT_TRUE(has_variable(variableNames, "y"));
}

TEST(UnitTestEvaluator, getAllVariables_twoIdenticalVariables)
{
  stk::expreval::Eval eval("x = y; z = y");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_variable_names();
  EXPECT_EQ(variableNames.size(), 3u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
  EXPECT_TRUE(has_variable(variableNames, "y"));
  EXPECT_TRUE(has_variable(variableNames, "z"));
}

TEST(UnitTestEvaluator, getAllVariables_twoVariables)
{
  stk::expreval::Eval eval("y = x; w = z");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_variable_names();
  EXPECT_EQ(variableNames.size(), 4u);
  EXPECT_TRUE(has_variable(variableNames, "y"));
  EXPECT_TRUE(has_variable(variableNames, "x"));
  EXPECT_TRUE(has_variable(variableNames, "w"));
  EXPECT_TRUE(has_variable(variableNames, "z"));
}

TEST(UnitTestEvaluator, getDependentVariables_noVariables)
{
  stk::expreval::Eval eval;
  eval.parse();
  EXPECT_EQ(eval.get_dependent_variable_names().size(), 0u);
}

TEST(UnitTestEvaluator, getDependentVariables_noAssign)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_EQ(eval.get_dependent_variable_names().size(), 0u);
}

TEST(UnitTestEvaluator, getDependentVariables_constant)
{
  stk::expreval::Eval eval("x = 2");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_dependent_variable_names();
  EXPECT_EQ(variableNames.size(), 1u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
}

TEST(UnitTestEvaluator, getDependentVariables_oneDependent)
{
  stk::expreval::Eval eval("x = sin(y)");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_dependent_variable_names();
  EXPECT_EQ(variableNames.size(), 1u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
}

TEST(UnitTestEvaluator, getDependentVariables_constantAssign)
{
  stk::expreval::Eval eval("x = 2; y = x");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_dependent_variable_names();
  EXPECT_EQ(variableNames.size(), 2u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
  EXPECT_TRUE(has_variable(variableNames, "y"));
}

TEST(UnitTestEvaluator, getDependentVariables_twoIdenticalVariables)
{
  stk::expreval::Eval eval("x = y; z = y");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_dependent_variable_names();
  EXPECT_EQ(variableNames.size(), 2u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
  EXPECT_TRUE(has_variable(variableNames, "z"));
}

TEST(UnitTestEvaluator, getDependentVariables_twoVariables)
{
  stk::expreval::Eval eval("y = x; w = z");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_dependent_variable_names();
  EXPECT_EQ(variableNames.size(), 2u);
  EXPECT_TRUE(has_variable(variableNames, "y"));
  EXPECT_TRUE(has_variable(variableNames, "w"));
}

TEST(UnitTestEvaluator, getIndependentVariables_noVariables)
{
  stk::expreval::Eval eval;
  eval.parse();
  EXPECT_EQ(eval.get_independent_variable_names().size(), 0u);
}

TEST(UnitTestEvaluator, getIndependentVariables_noAssign)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_independent_variable_names();
  EXPECT_EQ(variableNames.size(), 1u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
}

TEST(UnitTestEvaluator, getIndependentVariables_constant)
{
  stk::expreval::Eval eval("x = 2");
  eval.parse();
  EXPECT_EQ(eval.get_independent_variable_names().size(), 0u);
}

TEST(UnitTestEvaluator, getIndependentVariables_oneDependent)
{
  stk::expreval::Eval eval("x = sin(y)");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_independent_variable_names();
  EXPECT_EQ(variableNames.size(), 1u);
  EXPECT_TRUE(has_variable(variableNames, "y"));
}

TEST(UnitTestEvaluator, getIndependentVariables_constantAssign)
{
  stk::expreval::Eval eval("x = 2; y = x");
  eval.parse();
  EXPECT_EQ(eval.get_independent_variable_names().size(), 0u);
}

TEST(UnitTestEvaluator, getIndependentVariables_twoIdenticalVariables)
{
  stk::expreval::Eval eval("x = y; z = y");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_independent_variable_names();
  EXPECT_EQ(variableNames.size(), 1u);
  EXPECT_TRUE(has_variable(variableNames, "y"));
}

TEST(UnitTestEvaluator, getIndependentVariables_twoVariables)
{
  stk::expreval::Eval eval("y = x; w = z");
  eval.parse();
  std::vector<std::string> variableNames = eval.get_independent_variable_names();
  EXPECT_EQ(variableNames.size(), 2u);
  EXPECT_TRUE(has_variable(variableNames, "x"));
  EXPECT_TRUE(has_variable(variableNames, "z"));
}

TEST(UnitTestEvaluator, testThreadedEvaluation)
{
  // complicated expression that simplifies into 'x'
  std::string expression = " k = 7; h = -7; 1.0*x + 2.0 - (1.0/1.0) - 1.0*1.0 + 0.0*x + (x - x) + (k+h)";
  stk::expreval::Eval expr_eval(generate_expression(expression));

  expr_eval.parse();

  stk::ThreadLocalData<double> x;
  expr_eval.bindVariable("x", x);

  const int N = 10000;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < N; ++i) {
    x.getMyThreadEntry() = i;
    double y = expr_eval.evaluate();
    EXPECT_EQ(x.getMyThreadEntry(), y);
  }
}

TEST( UnitTestEvaluator, testEvaluateEmptyString)
{
    std::string emptyExpression = "";
    stk::expreval::Eval expr_eval(emptyExpression);
    expr_eval.parse();
    double result = kernel_evaluate(expr_eval);
    EXPECT_EQ(0.0, result);
}

TEST( UnitTestEvaluator, testIndexing)
{
  double y[3] = {1.1, 2.2, 3.3};

  {
    //  
    //  Check basic indexing with default zero offset
    //
    std::string expr0 = "y[0] + y[1] + y[2]";
    stk::expreval::Eval expr_eval0(expr0);
    expr_eval0.parse();
    expr_eval0.bindVariable("y", *y, 3);
    double result = kernel_evaluate(expr_eval0);
    EXPECT_EQ(6.6, result);
  }

  {
    //  
    //  Check basic indexing with default one offset
    //
    std::string expr1 = "y[1] + y[2] + y[3]";
    stk::expreval::Eval expr_eval1(expr1, stk::expreval::Variable::ONE_BASED_INDEX);
    expr_eval1.parse();
    expr_eval1.bindVariable("y", *y, 3);
    double result = kernel_evaluate(expr_eval1);
    EXPECT_EQ(6.6, result);
  }

  {
    //
    //  Error check for index above valid range
    //
    std::string expr1 = "y[0] + y[1] + y[3]";
    stk::expreval::Eval expr_eval1(expr1, stk::expreval::Variable::ZERO_BASED_INDEX);
    expr_eval1.parse();
    expr_eval1.bindVariable("y", *y, 3);
    try {
      expr_eval1.evaluate();
    } 
    catch (std::runtime_error &x) {
      std::string errMess = x.what();
      EXPECT_EQ(0, strcmp(errMess.c_str(),  "In analytic expression evaluator, processing variable 'y'.  Attempting to access invalid component '3' in analytic function.  Valid components are 0 to '2'.  "));
    }
  }

  {
    //
    //  Error check for index below valid range
    //
    std::string expr1 = "y[0] + y[1] + y[2]";
    stk::expreval::Eval expr_eval1(expr1, stk::expreval::Variable::ONE_BASED_INDEX);
    expr_eval1.parse();
    expr_eval1.bindVariable("y", *y, 3);
    try {
      expr_eval1.evaluate();
    } 
    catch (std::runtime_error &x) {
      std::string errMess = x.what();

      EXPECT_EQ(0, strcmp(errMess.c_str(),  "In analytic expression evaluator, processing variable 'y'.  Attempting to access invalid component '0' in analytic function.  Valid components are 1 to '3'.  "));
    }
  }

  {
    //
    //  Error check for lack of index
    //
    std::string expr1 = "y";
    stk::expreval::Eval expr_eval1(expr1, stk::expreval::Variable::ONE_BASED_INDEX);
    expr_eval1.parse();
    expr_eval1.bindVariable("y", *y, 3);
    try {
      expr_eval1.evaluate();
    } 
    catch (std::runtime_error &x) {
      std::string errMess = x.what();
      EXPECT_EQ(0, strcmp(errMess.c_str(),  "In analytic expression evaluator, processing variable 'y'.  Invalid direct access of array variable, must access by index"));
    }
  }

  {
    //  
    //  Check for variable index
    //
    std::string expr0 = "z=1; y[0] + y[z] + y[2]";
    stk::expreval::Eval expr_eval0(expr0);
    expr_eval0.parse();
    expr_eval0.bindVariable("y", *y, 3);
    double result = kernel_evaluate(expr_eval0);
    EXPECT_EQ(6.6, result);
  }

  {
    //  
    //  Check for compound index
    //
    std::string expr0 = "y[z[2]] + y[z[1]] + y[2]";
    stk::expreval::Eval expr_eval0(expr0);
    expr_eval0.parse();

    double z[3] = {2.0, 1.0, 0};

    expr_eval0.bindVariable("y", *y, 3);
    expr_eval0.bindVariable("z", *z, 3);
    double result = kernel_evaluate(expr_eval0);
    EXPECT_EQ(6.6, result);
  }

}

bool
checkUndefinedFunction(
    const char *	expr)
{
  stk::expreval::Eval expr_eval(expr);
  try {
    expr_eval.parse();
  }
  catch (std::runtime_error &x) {
    // parse throws on undefined function(s) and lists them
    std::cerr << x.what();
  }
  if (expr_eval.undefinedFunction()) {
    return true;
  } else {
    return false;
  }
}

bool
notUndefinedFunction(
    const char *	expr)
{
  std::cout << "Not an undefined function " << expr << " ... ";
  return !checkUndefinedFunction(expr);
}

bool
undefinedFunction(
    const char *	expr)
{
  std::cout << "Undefined function " << expr << " ... ";
  return checkUndefinedFunction(expr);
}

bool
syntax(
  const char *	expr)
{
  std::cout << "Syntax " << expr << " ... ";
  try {
    stk::expreval::Eval expr_eval(expr);
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
    stk::expreval::Eval expr_eval(expr);
    expr_eval.parse();
  }
  catch (std::runtime_error &x) {
    std::cout << "pass" << std::endl;
    return true;
  }
  std::cout << "fail, should have parse error" << std::endl;
  return false;
}

bool
test_one_value(const char *expression, double gold_value)
{
  bool failed = false;
  std::cout << "Evaluate " << expression << " ... ";
  stk::expreval::Eval expr_eval(generate_expression(expression));
  expr_eval.parse();

  double result = kernel_evaluate(expr_eval);

  double absolute_error = std::fabs(result - gold_value);
  if (absolute_error > std::fabs(1.0e-14*result)) 
  {
    std::cout << expression << " = " << std::setprecision(20) << result << " should be " << gold_value << " error= " << absolute_error << std::endl;
    failed = true;
  } else {
    std::cout << "Expression= " << expression << " == " << result << "\n";
  }
  std::cout << (failed ? "fail" : "pass") << std::endl;
  return !failed;
}

bool
test_one_value_is_equal(const char *expression, const char* gold_expression)
{
  stk::expreval::Eval gold_expr_eval(generate_expression(gold_expression));
  gold_expr_eval.parse();

  double gold_value = kernel_evaluate(gold_expr_eval); 

  return test_one_value(expression, gold_value);
}

typedef double (TestFunc)(double);

bool
test(const char* expr, TestFunc c_expr)
{
  bool failed = false;
  std::cout << "Evaluate " << expr << " ... ";
  stk::expreval::Eval expr_eval(generate_expression(expr));
  expr_eval.parse();

  double x=0, y=0, by=0, result = 0.0;
  double v[2];

  expr_eval.bindVariable("x", x);
  expr_eval.bindVariable("by", by);
  expr_eval.bindVariable("v", *v);
  for (int i = 1; i < 100; ++i) {
    x = v[1] = i*0.01;
    y = (*c_expr)(x);
    try {
      result = kernel_evaluate(expr_eval);
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
    double absolute_error = std::fabs(result - y);
    if (absolute_error > std::fabs(1.0e-14*result)) {
      std::cout << expr << " at " << std::setprecision(2) << x << " is " << std::setprecision(20) << result << " should be " << y << ", error= " << absolute_error
		<< std::endl;
      failed = true;
    }
    else if (by != result) {
      std::cout << expr << " at " << std::setprecision(2) << x << " is " << std::setprecision(20) << result << " does not match bound value " << std::setprecision(20) << by << std::endl;
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

// Hierarchy tests
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
EXPREVAL_DEFINE_TEST(f2, mod(x,10.0), fmod(x,10.0));
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
EXPREVAL_DEFINE_TEST(f15, x^2, pow(x, 2.0));

EXPREVAL_DEFINE_TEST1(f16, sin(x));
EXPREVAL_DEFINE_TEST1(f17, sinh(x));
EXPREVAL_DEFINE_TEST1(f18, sqrt(x));
EXPREVAL_DEFINE_TEST1(f19, tan(x));
EXPREVAL_DEFINE_TEST1(f20, tanh(x));
EXPREVAL_DEFINE_TEST(f21, atan2(x, PI),atan2(x, stk::expreval::pi() ));
EXPREVAL_DEFINE_TEST(f22, ln(x),log(x));
EXPREVAL_DEFINE_TEST(f23, deg(x),(180.0 / stk::expreval::pi() ) * x);
EXPREVAL_DEFINE_TEST(f24, rad(x),(stk::expreval::pi()  / 180.0) * x);
EXPREVAL_DEFINE_TEST(f25, max(x,1.0),std::max(x,1.0));
EXPREVAL_DEFINE_TEST(f26, min(x,1.0),std::min(x,1.0));
EXPREVAL_DEFINE_TEST(f27, recttopolr(x,1.0),sqrt(x*x+1.0*1.0));
EXPREVAL_DEFINE_TEST(f28, recttopola(x,1.0),atan2(1.0, x));
EXPREVAL_DEFINE_TEST(f29, poltorectx(x,PI/4.0),x*cos(stk::expreval::pi() /4.0));
EXPREVAL_DEFINE_TEST(f30, poltorecty(x,PI/4.0),x*sin(stk::expreval::pi() /4.0));
EXPREVAL_DEFINE_TEST1(f32, atanh(x));
EXPREVAL_DEFINE_TEST1(f33, 0.4209+4.5e-4*x);

// Bova tests
EXPREVAL_DEFINE_TEST1(b1, sin(x*.5));
EXPREVAL_DEFINE_TEST1(b2, .5*.2*sin(.5*x));
EXPREVAL_DEFINE_TEST1(b3, .5*sin(x));

// Pierson tests
EXPREVAL_DEFINE_TEST(k1, x^2, x*x);
EXPREVAL_DEFINE_TEST(k2, cosine_ramp(x),           (1.0-cos(x*pi() ))/2);
EXPREVAL_DEFINE_TEST(k3, cosine_ramp(x, 1.0),      (1.0-cos(x*pi() /1.0))/2);
EXPREVAL_DEFINE_TEST(k4, cosine_ramp(x, 0.0, 1.0), (1.0-cos(x*pi() /1.0))/2);

EXPREVAL_DEFINE_TEST(k5, haversine_pulse(x, 0.0, 1.0), std::pow(std::sin(pi() *x),2)   );
EXPREVAL_DEFINE_TEST(k6,  cycloidal_ramp(x, 0.0, 1.0), x-1/(two_pi())*sin(two_pi()*x) );

#undef EXPREVAL_DEFINE_TEST1

} // namespace <unnamed>

void
UnitTestEvaluator::testEvaluator()
{
  EXPECT_TRUE(syntax("3^2"));
  EXPECT_TRUE(test_one_value("a=1;b=2;a+b",3));
  EXPECT_TRUE(test_one_value("9%3" ,0));
  EXPECT_TRUE(test_one_value("9 % 4" ,1));
  EXPECT_TRUE(test_one_value("15%(1+1+1)",0));

  EXPECT_TRUE(test_one_value("max(1,2)",2));
  EXPECT_TRUE(test_one_value("max(1,2,3)",3));
  EXPECT_TRUE(test_one_value("max(1,2,3,4)",4));

  //
  //  Test overloaded function with wrong number of arguments
  //
  //EXPECT_THROW(test_one_value("min(4)",4), std::exception);

  EXPECT_TRUE(test_one_value("min(4,3)",3));
  EXPECT_TRUE(test_one_value("min(4,3,2)",2));
  EXPECT_TRUE(test_one_value("min(4,3,2,1)",1));

  EXPECT_TRUE(test_one_value("a=3;a^2",9.));
  EXPECT_TRUE(test_one_value("(1+2+3)^2",36.));
  EXPECT_TRUE(test_one_value("(1+2+3+4)^(1+1)",100.));
  EXPECT_TRUE(test_one_value("a=1;b=2;c=(a!=b);",1));
  EXPECT_TRUE(test_one_value("a=1;b=2;c=(a==b);",0));
  EXPECT_TRUE(test_one_value("a=1;b=!a;",0));
  EXPECT_TRUE(test_one_value("sign(7)", 1));
  EXPECT_TRUE(test_one_value("sign(-5)", -1));

  EXPECT_TRUE(test_one_value("ipart(5.1)", 5));
  EXPECT_TRUE(test_one_value_is_equal("ipart(5.2)", "ipart(5.1)")); 
  EXPECT_TRUE(test_one_value_is_equal("fpart(5.1234) + ipart(5.1234)", "5.1234")); 

  EXPECT_TRUE(test_one_value_is_equal("pow(5,2)", "5*5")); 
  EXPECT_TRUE(test_one_value_is_equal("min(1,2)", "min(2,1)")); 
  EXPECT_TRUE(test_one_value_is_equal("min(-1,0,1)", "min(0,1,-1)")); 

  EXPECT_TRUE(test_one_value_is_equal("x=-1; cycloidal_ramp(x,0,1)", "0.0")); 
  EXPECT_TRUE(test_one_value_is_equal("x= 2; cycloidal_ramp(x,0,1)", "1.0")); 
  EXPECT_TRUE(test_one_value_is_equal("x=-1; haversine_pulse(x,0,1)", "0.0")); 
  EXPECT_TRUE(test_one_value_is_equal("x= 2; haversine_pulse(x,0,1)", "0.0")); 

  //
  //  Test some composite functions
  //
  EXPECT_TRUE(test_one_value("pow(pow(2,2), pow(2,2))", 256));
  EXPECT_TRUE(test_one_value("pow(2^2, 2^2)", 256));
  EXPECT_TRUE(test_one_value("(2^2)^(2^2)", 256));
  EXPECT_TRUE(test_one_value("(2^2)^pow(2,2)", 256));
  EXPECT_TRUE(test_one_value("sin(pow(2,2))", sin(4)));
  EXPECT_TRUE(test_one_value("sin(sin(sin(sin(100))))",sin(sin(sin(sin(100.0))))));
  EXPECT_TRUE(test_one_value("sin(sin(sin(sin(pow(2,2)))))",sin(sin(sin(sin(4.0))))));

  const double weibull_gold_value = 3.6787944117144233402;
  EXPECT_TRUE(test_one_value("shape=10;scale=1;weibull_pdf(1.0,shape,scale)", weibull_gold_value));

  // Need a better test for distributions, perhaps something that computes the
  // mean and standard deviation of the distribution.
  const double normal_gold_value = 0.79788456080286540573;
  EXPECT_TRUE(test_one_value("mean=0;standard_deviation=0.5;normal_pdf(0.0,mean,standard_deviation)", normal_gold_value));

  // Simple test of ts_random
  const double tsrand_gold_value = 0.64737415343866411277;
  EXPECT_TRUE(test_one_value("ts_random(1e-3,0,1e-6,1e-4)", tsrand_gold_value));

  // Better test of ts_random
  {
    double t = 1e-4;
    double x = 0;
    double y = 0;
    double z = 0;
    stk::expreval::Eval expr_eval("ts_random(t,x,y,z)");
    expr_eval.parse();
    expr_eval.bindVariable("t",t);
    expr_eval.bindVariable("x",x);
    expr_eval.bindVariable("y",y);
    expr_eval.bindVariable("z",z);

    stk::expreval::Eval expr_eval_norm("ts_normal(t,x,y,z,1,1.2,-10,10)");
    expr_eval_norm.parse();
    expr_eval_norm.bindVariable("t",t);
    expr_eval_norm.bindVariable("x",x);
    expr_eval_norm.bindVariable("y",y);
    expr_eval_norm.bindVariable("z",z);

    // Ensure uniform distribution with different spatial length scales
    for(auto&& L : {1e-6,1e3}) {

      std::vector<int> bins(10,0);
      double mean = 0.0;
      double den = 0.0;
      double sigma = 0.0;

      for(int i = 0; i < 100; ++i) {
        x = L*(i-50)/100.0;
        for(int j = 0; j < 100; ++j) {
          y = L*(j-50)/100.0;
          double result = kernel_evaluate(expr_eval);
          int idx = static_cast<int>(result*10);
          bins[idx] += 1;

          result = kernel_evaluate(expr_eval_norm);
          mean += result;
          sigma += (result-1)*(result-1);
          den += 1;
        }
      }

      mean /= den;
      sigma = std::sqrt(sigma/den);

      std::cout << "Mu = " << mean << ", s = " << sigma << std::endl;

      EXPECT_TRUE(std::abs(mean - 1) < 0.02);
      EXPECT_TRUE(std::abs(sigma - 1.2) < 0.02);

      int maxN = *std::max_element(bins.begin(), bins.end());
      int minN = *std::min_element(bins.begin(), bins.end());

      std::cout << "Min/Max = " << minN << "/" << maxN << std::endl;
      EXPECT_TRUE(maxN < 1100);
      EXPECT_TRUE(minN > 900);
    }

    // test uniformity in time
    x = 1e-4;
    y = 2.2e-3;
    z = 1e-5;

    {
      std::vector<int> bins(10,0);
      for(int i = 0; i < 10000; ++i) {
        double result = kernel_evaluate(expr_eval);
        int idx = static_cast<int>(result*10);
        bins[idx] += 1;
        t += 1e-4;
      }

      int maxN = *std::max_element(bins.begin(), bins.end());
      int minN = *std::min_element(bins.begin(), bins.end());
      std::cout << "Min/Max = " << minN << "/" << maxN << std::endl;
      EXPECT_TRUE(maxN < 1100);
      EXPECT_TRUE(minN > 900);
    }

  }

  // Test point mask functions
  EXPECT_TRUE(test_one_value("point2d(0,0,1,0.1)", 1.0));
  EXPECT_TRUE(test_one_value("point2d(0,1,1,0.1)", 0.5));
  EXPECT_TRUE(test_one_value("point2d(0,2,1,0.1)", 0.0));

  EXPECT_TRUE(test_one_value("point3d(0,0,0,1,0.1)", 1.0));
  EXPECT_TRUE(test_one_value("point3d(0,1,0,1,0.1)", 0.5));
  EXPECT_TRUE(test_one_value("point3d(0,2,0,1,0.1)", 0.0));

  EXPECT_TRUE(syntax("2*2"));
  EXPECT_TRUE(syntax(""));
  EXPECT_TRUE(syntax(";"));
  EXPECT_TRUE(syntax(";;"));
  EXPECT_TRUE(syntax(";;;"));
  EXPECT_TRUE(syntax("x*0.1"));
  EXPECT_TRUE(syntax("x*-0.1"));
  EXPECT_TRUE(syntax("x*+0.1"));
  EXPECT_TRUE(syntax("x--7.0"));
  EXPECT_TRUE(syntax("x*-x"));
  EXPECT_TRUE(syntax("x*+x"));
  EXPECT_TRUE(syntax("v[0]=v[1]*0.1"));
  EXPECT_TRUE(syntax("x--x"));
  EXPECT_TRUE(syntax("x---x"));

  EXPECT_TRUE(fail_syntax("0.01.02"));
  EXPECT_TRUE(fail_syntax("5*.e+10"));
  EXPECT_TRUE(fail_syntax("x y"));
  EXPECT_TRUE(fail_syntax("x(y"));
  EXPECT_TRUE(fail_syntax("x*"));
  EXPECT_TRUE(fail_syntax("x*(y+1"));
  EXPECT_TRUE(fail_syntax("x*y)"));
  EXPECT_TRUE(fail_syntax("cos(x"));
  EXPECT_TRUE(fail_syntax("(x)y"));
  EXPECT_TRUE(fail_syntax("()"));

  EXPECT_TRUE(undefinedFunction("stress(1)"));
  EXPECT_TRUE(notUndefinedFunction("sin(1)"));
  EXPECT_TRUE(notUndefinedFunction("0.01.02"));
  EXPECT_TRUE(syntax("rand()"));
  EXPECT_TRUE(syntax("srand()"));
  EXPECT_TRUE(syntax("time()"));
  EXPECT_TRUE(syntax("random()"));
  EXPECT_TRUE(syntax("random(1)"));
  EXPECT_TRUE(syntax("random(time())"));
  EXPECT_TRUE(syntax("cosine_ramp(x,y)"));
  EXPECT_TRUE(syntax("sign(x)"));

  EXPECT_TRUE(syntax("weibull_pdf(x, alpha, beta)"));
  EXPECT_TRUE(syntax( "normal_pdf(x, alpha, beta)"));

#define EXPREVAL_TEST(name) test(name##_expr, name)

  EXPECT_TRUE(EXPREVAL_TEST(h1));
  EXPECT_TRUE(EXPREVAL_TEST(h2));
  EXPECT_TRUE(EXPREVAL_TEST(h3));
  EXPECT_TRUE(EXPREVAL_TEST(h4));
  EXPECT_TRUE(EXPREVAL_TEST(h5));
  EXPECT_TRUE(EXPREVAL_TEST(h6));
  EXPECT_TRUE(EXPREVAL_TEST(h7));
  EXPECT_TRUE(EXPREVAL_TEST(h8));
  EXPECT_TRUE(EXPREVAL_TEST(h9));
  EXPECT_TRUE(EXPREVAL_TEST(h10));
  EXPECT_TRUE(EXPREVAL_TEST(h11));
  EXPECT_TRUE(EXPREVAL_TEST(h12));
  EXPECT_TRUE(EXPREVAL_TEST(h13));
  EXPECT_TRUE(EXPREVAL_TEST(h14));
  EXPECT_TRUE(EXPREVAL_TEST(h15));
  EXPECT_TRUE(EXPREVAL_TEST(h16));
  EXPECT_TRUE(EXPREVAL_TEST(h17));
  EXPECT_TRUE(EXPREVAL_TEST(h18));
  EXPECT_TRUE(EXPREVAL_TEST(h19));
  EXPECT_TRUE(EXPREVAL_TEST(h20));

  EXPECT_TRUE(EXPREVAL_TEST(f1));
  EXPECT_TRUE(EXPREVAL_TEST(f2));
  EXPECT_TRUE(EXPREVAL_TEST(f3));
  EXPECT_TRUE(EXPREVAL_TEST(f4));
  EXPECT_TRUE(EXPREVAL_TEST(f5));
  EXPECT_TRUE(EXPREVAL_TEST(f6));
  EXPECT_TRUE(EXPREVAL_TEST(f7));
  EXPECT_TRUE(EXPREVAL_TEST(f8));
  EXPECT_TRUE(EXPREVAL_TEST(f9));
  EXPECT_TRUE(EXPREVAL_TEST(f10));
  EXPECT_TRUE(EXPREVAL_TEST(f11));
  EXPECT_TRUE(EXPREVAL_TEST(f12));
  EXPECT_TRUE(EXPREVAL_TEST(f13));
  EXPECT_TRUE(EXPREVAL_TEST(f14));
  EXPECT_TRUE(EXPREVAL_TEST(f15));
  EXPECT_TRUE(EXPREVAL_TEST(f16));
  EXPECT_TRUE(EXPREVAL_TEST(f17));
  EXPECT_TRUE(EXPREVAL_TEST(f18));
  EXPECT_TRUE(EXPREVAL_TEST(f19));
  EXPECT_TRUE(EXPREVAL_TEST(f20));
  EXPECT_TRUE(EXPREVAL_TEST(f21));
  EXPECT_TRUE(EXPREVAL_TEST(f22));
  EXPECT_TRUE(EXPREVAL_TEST(f23));
  EXPECT_TRUE(EXPREVAL_TEST(f24));
  EXPECT_TRUE(EXPREVAL_TEST(f25));
  EXPECT_TRUE(EXPREVAL_TEST(f26));
  EXPECT_TRUE(EXPREVAL_TEST(f27));
  EXPECT_TRUE(EXPREVAL_TEST(f28));
  EXPECT_TRUE(EXPREVAL_TEST(f29));
  EXPECT_TRUE(EXPREVAL_TEST(f30));
  EXPECT_TRUE(EXPREVAL_TEST(f32));
  EXPECT_TRUE(EXPREVAL_TEST(f33));

  EXPECT_TRUE(EXPREVAL_TEST(b1));
  EXPECT_TRUE(EXPREVAL_TEST(b2));
  EXPECT_TRUE(EXPREVAL_TEST(b3));

  EXPECT_TRUE(EXPREVAL_TEST(k1));
  EXPECT_TRUE(EXPREVAL_TEST(k2));
  EXPECT_TRUE(EXPREVAL_TEST(k3));
  EXPECT_TRUE(EXPREVAL_TEST(k4));
  EXPECT_TRUE(EXPREVAL_TEST(k5));
  EXPECT_TRUE(EXPREVAL_TEST(k6));

#undef EXPREVAL_TEST
}
