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

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>
#include <stk_expreval/Evaluator.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <iomanip>
#include <cmath>

using ExecutionSpace = Kokkos::Serial;

namespace {

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

bool
test_is_scalar(const stk::expreval::Eval& eval, const std::string& variableName)
{
  return eval.is_scalar(variableName);
}

TEST(UnitTestEvaluator, isConstantExpression_empty)
{
  stk::expreval::Eval eval;
  eval.parse();
  EXPECT_TRUE(eval.is_constant_expression());
}

TEST(UnitTestEvaluator, isConstantExpression_constant)
{
  stk::expreval::Eval eval("2");
  eval.parse();
  EXPECT_TRUE(eval.is_constant_expression());
}

TEST(UnitTestEvaluator, isConstantExpression_variable)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_FALSE(eval.is_constant_expression());
}

TEST(UnitTestEvaluator, isVariable_no)
{
  stk::expreval::Eval eval;
  eval.parse();
  EXPECT_FALSE(eval.is_variable("x"));
}

TEST(UnitTestEvaluator, isVariable_yes)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_TRUE(eval.is_variable("x"));
}

TEST(UnitTestEvaluator, isVariable_twoVariables_yes)
{
  stk::expreval::Eval eval("x + y");
  eval.parse();
  EXPECT_TRUE(eval.is_variable("x"));
  EXPECT_TRUE(eval.is_variable("y"));
}

TEST(UnitTestEvaluator, isVariable_twoVariables_no)
{
  stk::expreval::Eval eval("x + y");
  eval.parse();
  EXPECT_TRUE(eval.is_variable("x"));
  EXPECT_FALSE(eval.is_variable("z"));
}

TEST(UnitTestEvaluator, isScalar_default)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_TRUE(test_is_scalar(eval, "x"));
}

TEST(UnitTestEvaluator, isScalar_notPresent)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_FALSE(test_is_scalar(eval, "y"));
}

TEST(UnitTestEvaluator, isScalar_assignYes)
{
  stk::expreval::Eval eval("x = 2.0");
  eval.parse();
  EXPECT_TRUE(test_is_scalar(eval, "x"));
}

TEST(UnitTestEvaluator, isScalar_bindDefaultYes)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  double x = 3.0;
  eval.bindVariable("x", x);
  EXPECT_TRUE(test_is_scalar(eval, "x"));
}

TEST(UnitTestEvaluator, isScalar_bindYes)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  double x = 3.0;
  eval.bindVariable("x", x, 1);
  EXPECT_TRUE(test_is_scalar(eval, "x"));
}

TEST(UnitTestEvaluator, isScalar_bindNo)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  double x[3] = {4.0, 5.0, 6.0};
  eval.bindVariable("x", *x, 3);
  EXPECT_FALSE(test_is_scalar(eval, "x"));
}

TEST(UnitTestEvaluator, isScalar_bindYesAndNo)
{
  stk::expreval::Eval eval("z = y[1]");
  eval.parse();
  double y[3] = {4.0, 5.0, 6.0};
  eval.bindVariable("y", *y, 3);
  EXPECT_FALSE(test_is_scalar(eval, "y"));
  EXPECT_TRUE(test_is_scalar(eval, "z"));
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
isValidParse(const char *expr)
{
  stk::expreval::Eval expr_eval(expr);
  EXPECT_NO_THROW(expr_eval.parse());
  return expr_eval.getSyntaxStatus();
}

bool
isInvalidParse(const char *expr)
{
  stk::expreval::Eval expr_eval(expr);
  try {
    expr_eval.parse();
  }
  catch (std::runtime_error& ) {
    return !expr_eval.getSyntaxStatus();
  }

  return false;
}

TEST(UnitTestEvaluator, testAlgebraicSyntax)
{
  EXPECT_TRUE(isValidParse(""));
  EXPECT_TRUE(isValidParse(";"));
  EXPECT_TRUE(isValidParse(";;"));
  EXPECT_TRUE(isValidParse(";;;"));
  EXPECT_TRUE(isValidParse("2*2"));
  EXPECT_TRUE(isValidParse("3^2"));
  EXPECT_TRUE(isValidParse("x*0.1"));
  EXPECT_TRUE(isValidParse("x*-0.1"));
  EXPECT_TRUE(isValidParse("x*+0.1"));
  EXPECT_TRUE(isValidParse("x--7.0"));
  EXPECT_TRUE(isValidParse("x*-x"));
  EXPECT_TRUE(isValidParse("x*+x"));
  EXPECT_TRUE(isValidParse("v[0]=v[1]*0.1"));
  EXPECT_TRUE(isValidParse("x--x"));
  EXPECT_TRUE(isValidParse("x---x"));

  EXPECT_FALSE(isInvalidParse("2*2"));
  EXPECT_TRUE(isInvalidParse("0.01.02"));
  EXPECT_TRUE(isInvalidParse("5*.e+10"));
  EXPECT_TRUE(isInvalidParse("x y"));
  EXPECT_TRUE(isInvalidParse("x(y"));
  EXPECT_TRUE(isInvalidParse("x*"));
  EXPECT_TRUE(isInvalidParse("x*(y+1"));
  EXPECT_TRUE(isInvalidParse("x*y)"));
  EXPECT_TRUE(isInvalidParse("cos(x"));
  EXPECT_TRUE(isInvalidParse("(x)y"));
  EXPECT_TRUE(isInvalidParse("()"));
}

bool
isValidFunction(const char *expr)
{
  stk::expreval::Eval expr_eval(expr);
  EXPECT_NO_THROW(expr_eval.parse());
  return !expr_eval.undefinedFunction();
}

bool
isInvalidFunction(const char *expr)
{
  stk::expreval::Eval expr_eval(expr);
  try {
    expr_eval.parse();
  }
  catch (std::runtime_error& ) {
    return expr_eval.undefinedFunction();
  }

  return false;
}

TEST(UnitTestEvaluator, testFunctionSyntax)
{
  EXPECT_TRUE(isValidFunction("sin(1)"));
  EXPECT_TRUE(isValidFunction("SIN(1)"));
  EXPECT_TRUE(isValidFunction("rand()"));
  EXPECT_TRUE(isValidFunction("srand()"));
  EXPECT_TRUE(isValidFunction("time()"));
  EXPECT_TRUE(isValidFunction("random()"));
  EXPECT_TRUE(isValidFunction("random(1)"));
  EXPECT_TRUE(isValidFunction("random(time())"));
  EXPECT_TRUE(isValidFunction("cosine_ramp(x,y)"));
  EXPECT_TRUE(isValidFunction("sign(x)"));
  EXPECT_TRUE(isValidFunction("weibull_pdf(x, alpha, beta)"));
  EXPECT_TRUE(isValidFunction("normal_pdf(x, alpha, beta)"));

  EXPECT_FALSE(isInvalidFunction("sin(1)"));
  EXPECT_TRUE(isInvalidFunction("stress(1)"));
  EXPECT_TRUE(isInvalidFunction("gamma(1)"));
}

struct ScalarBinding {
  std::string varName;
  double varValue;
};

struct VectorBinding {
  std::string varName;
  std::vector<double> varValues;
};

double evaluate(const std::string & expression,
                std::vector<ScalarBinding> boundScalars = std::vector<ScalarBinding>(),
                std::vector<VectorBinding> boundVectors = std::vector<VectorBinding>(),
                const stk::expreval::Variable::ArrayOffset arrayOffsetType = stk::expreval::Variable::ZERO_BASED_INDEX)
{
  stk::expreval::Eval eval(expression, arrayOffsetType);
  eval.parse();

  for (ScalarBinding & scalar : boundScalars) {
    eval.bindVariable(scalar.varName, scalar.varValue, 1);
  }

  for (VectorBinding & vector : boundVectors) {
    eval.bindVariable(vector.varName, *vector.varValues.data(), vector.varValues.size());
  }

  return eval.evaluate();
}

TEST(UnitTestEvaluator, testOpcode_CONSTANT)
{
  EXPECT_DOUBLE_EQ(evaluate("-1000"),  -1000);
  EXPECT_DOUBLE_EQ(evaluate("-1.333"), -1.333);
  EXPECT_DOUBLE_EQ(evaluate("-1"),     -1);
  EXPECT_DOUBLE_EQ(evaluate("0"),       0);
  EXPECT_DOUBLE_EQ(evaluate("1.01"),    1.01);
  EXPECT_DOUBLE_EQ(evaluate("999"),     999);
  EXPECT_DOUBLE_EQ(evaluate("3.e-2"),   0.03);
 }

TEST(UnitTestEvaluator, testOpcode_ADD)
{
  EXPECT_DOUBLE_EQ(evaluate("1+2"),          3);
  EXPECT_DOUBLE_EQ(evaluate("1+4+9"),        14);
  EXPECT_DOUBLE_EQ(evaluate("1+4+9+16"),     30);
  EXPECT_DOUBLE_EQ(evaluate("(1+4)+(9+16)"), 30);
  EXPECT_DOUBLE_EQ(evaluate("1.25+2.5"),     3.75);
}

TEST(UnitTestEvaluator, testOpcode_SUBTRACT)
{
  EXPECT_DOUBLE_EQ(evaluate("-1-2"),         -3);
  EXPECT_DOUBLE_EQ(evaluate("-1-4-9"),       -14);
  EXPECT_DOUBLE_EQ(evaluate("-1-4-9-16"),    -30);
  EXPECT_DOUBLE_EQ(evaluate("(-1-4)-(9-16)"), 2);
  EXPECT_DOUBLE_EQ(evaluate("1.25-2.5"),     -1.25);
}

TEST(UnitTestEvaluator, testOpcode_MULTIPLY)
{
  EXPECT_DOUBLE_EQ(evaluate("2*3"),         6);
  EXPECT_DOUBLE_EQ(evaluate("2*3*4"),       24);
  EXPECT_DOUBLE_EQ(evaluate("2*3*4*5"),     120);
  EXPECT_DOUBLE_EQ(evaluate("(2*3)*(4*5)"), 120);
  EXPECT_DOUBLE_EQ(evaluate("0.25*4"),      1);
  EXPECT_DOUBLE_EQ(evaluate("-0.25*4"),    -1);
  EXPECT_DOUBLE_EQ(evaluate("0.25*-4"),    -1);
}

TEST(UnitTestEvaluator, testOpcode_DIVIDE)
{
  EXPECT_DOUBLE_EQ(evaluate("6/3"),           2);
  EXPECT_DOUBLE_EQ(evaluate("24/4/3"),        2);
  EXPECT_DOUBLE_EQ(evaluate("120/5/4/3"),     2);
  EXPECT_DOUBLE_EQ(evaluate("(120/5)/(4/3)"), 18);
  EXPECT_DOUBLE_EQ(evaluate("1/4"),           0.25);
  EXPECT_DOUBLE_EQ(evaluate("-1/4"),         -0.25);
  EXPECT_DOUBLE_EQ(evaluate("1/-4"),         -0.25);
}

TEST(UnitTestEvaluator, testOpcode_UNARY_MINUS)
{
  EXPECT_DOUBLE_EQ(evaluate("-2"),       -2);
  EXPECT_DOUBLE_EQ(evaluate("--2"),       2);
  EXPECT_DOUBLE_EQ(evaluate("---2"),     -2);
  EXPECT_DOUBLE_EQ(evaluate("-(1+2)"),   -3);
  EXPECT_DOUBLE_EQ(evaluate("-(-(1+2))"), 3);
}

TEST(UnitTestEvaluator, testOpcode_MODULUS)
{
  EXPECT_DOUBLE_EQ(evaluate("-9 % 3"),   0);
  EXPECT_DOUBLE_EQ(evaluate("-9 % -4"), -1);
  EXPECT_DOUBLE_EQ(evaluate("-9 % 4"),  -1);
  EXPECT_DOUBLE_EQ(evaluate("9 % -4"),   1);
  EXPECT_DOUBLE_EQ(evaluate("9 % 3"),    0);
  EXPECT_DOUBLE_EQ(evaluate("9 % 4"),    1);
  EXPECT_DOUBLE_EQ(evaluate("9.1 % 4"),  1.1);
  EXPECT_DOUBLE_EQ(evaluate("9 % 5"),    4);
  EXPECT_DOUBLE_EQ(evaluate("9 % 9"),    0);
}

TEST(UnitTestEvaluator, testOpcode_EXPONENTIATION)
{
  EXPECT_DOUBLE_EQ(evaluate("3^2"),         9);
  EXPECT_DOUBLE_EQ(evaluate("3^2.5"),       std::pow(3, 2.5));
  EXPECT_DOUBLE_EQ(evaluate("3^2^2"),       81);
  EXPECT_DOUBLE_EQ(evaluate("(2^2)^(2^2)"), 256);
  EXPECT_DOUBLE_EQ(evaluate("2^-1"),        0.5);
}

TEST(UnitTestEvaluator, testOpcode_compoundSimpleMath)
{
  EXPECT_DOUBLE_EQ(evaluate("1+2*3"),           7);
  EXPECT_DOUBLE_EQ(evaluate("1+2*-3"),         -5);
  EXPECT_DOUBLE_EQ(evaluate("2*3+1"),           7);
  EXPECT_DOUBLE_EQ(evaluate("12*3/3*5/5"),      12);
  EXPECT_DOUBLE_EQ(evaluate("(4+5)/3"),         3);
  EXPECT_DOUBLE_EQ(evaluate("(1+2+3)^2"),       36);
  EXPECT_DOUBLE_EQ(evaluate("(1+2+3+4)^(1+1)"), 100);
  EXPECT_DOUBLE_EQ(evaluate("15%(1+1+1)"),      0);
  EXPECT_DOUBLE_EQ(evaluate("1-7"),            -6);
  EXPECT_DOUBLE_EQ(evaluate("1--7"),            8);
  EXPECT_DOUBLE_EQ(evaluate("1---7"),          -6);
}

TEST(UnitTestEvaluator, testOpcode_EQUAL)
{
  EXPECT_DOUBLE_EQ(evaluate("2==2"),          1);
  EXPECT_DOUBLE_EQ(evaluate("2==(1+1)"),      1);
  EXPECT_DOUBLE_EQ(evaluate("2==1"),          0);
  EXPECT_DOUBLE_EQ(evaluate("0.1==0.999999"), 0);
}

TEST(UnitTestEvaluator, testOpcode_NOT_EQUAL)
{
  EXPECT_DOUBLE_EQ(evaluate("2!=2"),          0);
  EXPECT_DOUBLE_EQ(evaluate("2!=(1+1)"),      0);
  EXPECT_DOUBLE_EQ(evaluate("2!=1"),          1);
  EXPECT_DOUBLE_EQ(evaluate("0.1!=0.999999"), 1);
}

TEST(UnitTestEvaluator, testOpcode_LESS)
{
  EXPECT_DOUBLE_EQ(evaluate("1<2"),        1);
  EXPECT_DOUBLE_EQ(evaluate("2<1"),        0);
  EXPECT_DOUBLE_EQ(evaluate("-2<-1"),      1);
  EXPECT_DOUBLE_EQ(evaluate("1<1"),        0);
  EXPECT_DOUBLE_EQ(evaluate("1<1.000001"), 1);
  EXPECT_DOUBLE_EQ(evaluate("2<(1+2)"),    1);
}

TEST(UnitTestEvaluator, testOpcode_GREATER)
{
  EXPECT_DOUBLE_EQ(evaluate("1>2"),        0);
  EXPECT_DOUBLE_EQ(evaluate("2>1"),        1);
  EXPECT_DOUBLE_EQ(evaluate("-2>-1"),      0);
  EXPECT_DOUBLE_EQ(evaluate("1>1"),        0);
  EXPECT_DOUBLE_EQ(evaluate("1>1.000001"), 0);
  EXPECT_DOUBLE_EQ(evaluate("2>(1+2)"),    0);
}

TEST(UnitTestEvaluator, testOpcode_LESS_EQUAL)
{
  EXPECT_DOUBLE_EQ(evaluate("1<=2"),        1);
  EXPECT_DOUBLE_EQ(evaluate("2<=1"),        0);
  EXPECT_DOUBLE_EQ(evaluate("-2<=-1"),      1);
  EXPECT_DOUBLE_EQ(evaluate("1<=1"),        1);
  EXPECT_DOUBLE_EQ(evaluate("1<=1.000001"), 1);
  EXPECT_DOUBLE_EQ(evaluate("2<=(1+2)"),    1);
}

TEST(UnitTestEvaluator, testOpcode_GREATER_EQUAL)
{
  EXPECT_DOUBLE_EQ(evaluate("1>=2"),        0);
  EXPECT_DOUBLE_EQ(evaluate("2>=1"),        1);
  EXPECT_DOUBLE_EQ(evaluate("-2>=-1"),      0);
  EXPECT_DOUBLE_EQ(evaluate("1>=1"),        1);
  EXPECT_DOUBLE_EQ(evaluate("1>=1.000001"), 0);
  EXPECT_DOUBLE_EQ(evaluate("2>=(1+2)"),    0);
}

TEST(UnitTestEvaluator, testOpcode_UNARY_NOT)
{
  EXPECT_DOUBLE_EQ(evaluate("!0"),        1);
  EXPECT_DOUBLE_EQ(evaluate("!0.000001"), 0);
  EXPECT_DOUBLE_EQ(evaluate("!1"),        0);
  EXPECT_DOUBLE_EQ(evaluate("!10"),       0);
  EXPECT_DOUBLE_EQ(evaluate("!-1"),       0);
  EXPECT_DOUBLE_EQ(evaluate("!-10"),      0);
}

TEST(UnitTestEvaluator, testOpcode_LOGICAL_AND)
{
  EXPECT_DOUBLE_EQ(evaluate("0 && 0"),        0);
  EXPECT_DOUBLE_EQ(evaluate("0 && 1"),        0);
  EXPECT_DOUBLE_EQ(evaluate("1 && 0"),        0);
  EXPECT_DOUBLE_EQ(evaluate("1 && 1"),        1);
  EXPECT_DOUBLE_EQ(evaluate("1 && -1"),       1);
  EXPECT_DOUBLE_EQ(evaluate("1 && 0.000001"), 1);
 }

TEST(UnitTestEvaluator, testOpcode_LOGICAL_OR)
{
  EXPECT_DOUBLE_EQ(evaluate("0 || 0"),        0);
  EXPECT_DOUBLE_EQ(evaluate("0 || 1"),        1);
  EXPECT_DOUBLE_EQ(evaluate("1 || 0"),        1);
  EXPECT_DOUBLE_EQ(evaluate("1 || 1"),        1);
  EXPECT_DOUBLE_EQ(evaluate("0 || -1"),       1);
  EXPECT_DOUBLE_EQ(evaluate("0 || 0.000001"), 1);
}

TEST(UnitTestEvaluator, testOpcode_TERNARY)
{
  EXPECT_DOUBLE_EQ(evaluate("1 ? 1 : 2"),           1);
  EXPECT_DOUBLE_EQ(evaluate("0 ? 1 : 2"),           2);
  EXPECT_DOUBLE_EQ(evaluate("1 ? (1+1) : 1"),       2);
  EXPECT_DOUBLE_EQ(evaluate("0 ? 1 : (1+1)"),       2);
  EXPECT_DOUBLE_EQ(evaluate("0 ? 1 : 2"),           2);
  EXPECT_DOUBLE_EQ(evaluate("0.000001 ? 1 : 2"),    1);
  EXPECT_DOUBLE_EQ(evaluate("(1 ? 0 : 1) ? 2 : 3"), 3);
}

TEST(UnitTestEvaluator, testOpcode_ASSIGN)
{
  EXPECT_DOUBLE_EQ(evaluate("x=1"),               1);
  EXPECT_DOUBLE_EQ(evaluate("x=1; y=2"),          2);
  EXPECT_DOUBLE_EQ(evaluate("x=1; y=2; z=3"),     3);
  EXPECT_DOUBLE_EQ(evaluate("x=1; y=2; z=(1+2)"), 3);
}

TEST(UnitTestEvaluator, testOpcode_RVALUE)
{
  EXPECT_DOUBLE_EQ(evaluate("x=1; !x"),                  0);
  EXPECT_DOUBLE_EQ(evaluate("x=1; y=!x"),                0);
  EXPECT_DOUBLE_EQ(evaluate("x=0.5; !(x - 0.5)"),        1);
  EXPECT_DOUBLE_EQ(evaluate("x=0.6; !(x - 0.5)"),        0);
  EXPECT_DOUBLE_EQ(evaluate("x=2; x^2"),                 4);
  EXPECT_DOUBLE_EQ(evaluate("x=1; y=2; x+y"),            3);
  EXPECT_DOUBLE_EQ(evaluate("x=2; y=3; x*y"),            6);
  EXPECT_DOUBLE_EQ(evaluate("x=3; y=2; z=x/y; w=z-1"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("x=1; y=2; z=(x != y)"),     1);
  EXPECT_DOUBLE_EQ(evaluate("x=1; y=2; z=(x == y)"),     0);
  EXPECT_DOUBLE_EQ(evaluate("x=0.1; x<0.25 || x>0.75"),  1);
  EXPECT_DOUBLE_EQ(evaluate("x=0.3; x<0.25 || x>0.75"),  0);
  EXPECT_DOUBLE_EQ(evaluate("x=0.1; x>0.25 && x<0.75"),  0);
  EXPECT_DOUBLE_EQ(evaluate("x=0.3; x>0.25 && x<0.75"),  1);
  EXPECT_DOUBLE_EQ(evaluate("x=0.4; (x > 0.5) ? 2 : 3"), 3);
  EXPECT_DOUBLE_EQ(evaluate("x=0.6; (x > 0.5) ? 2 : 3"), 2);
  EXPECT_DOUBLE_EQ(evaluate("x=0.5; x*(x + 1)"),         0.75);
}

TEST(UnitTestEvaluator, bindScalar)
{
  EXPECT_DOUBLE_EQ(evaluate("x",               {{"x", 2}}),                     2);
  EXPECT_DOUBLE_EQ(evaluate("x*x",             {{"x", 2}}),                     4);
  EXPECT_DOUBLE_EQ(evaluate("x*y",             {{"x", 2}, {"y", 3}}),           6);
  EXPECT_DOUBLE_EQ(evaluate("x+y*z",           {{"x", 2}, {"y", 3}, {"z", 4}}), 14);
  EXPECT_DOUBLE_EQ(evaluate("x=5; y=y+x; y+z", {{"x", 2}, {"y", 3}, {"z", 4}}), 12);
}

TEST(UnitTestEvaluator, bindVector)
{
  EXPECT_DOUBLE_EQ(evaluate("a[0]+a[1]+a[2]",           {},         {{"a", {1, 2, 3}}}), 6);
  EXPECT_DOUBLE_EQ(evaluate("a[0]=x*a[1]+a[2]",         {{"x", 4}}, {{"a", {1, 2, 3}}}), 11);
  EXPECT_DOUBLE_EQ(evaluate("a[1]=x*0.5+a[0]; a[1]*2",  {{"x", 2}}, {{"a", {3, 4}}}),    8);
  EXPECT_DOUBLE_EQ(evaluate("(a[0]*b[0] + a[1]*b[1] + a[2]*b[2])^0.5",
                            {}, {{"a", {1, 2, 3}}, {"b", {5, 4, 4}}}),                   5);
  EXPECT_DOUBLE_EQ(evaluate("(a[1]*b[1] + a[2]*b[2] + a[3]*b[3])^0.5",
                            {}, {{"a", {1, 2, 3}}, {"b", {5, 4, 4}}},
                            stk::expreval::Variable::ONE_BASED_INDEX),                   5);
}

TEST(UnitTestEvaluator, testFunction_abs)
{
  EXPECT_DOUBLE_EQ(evaluate("abs(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("abs(1)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("abs(1.5)"),  1.5);
  EXPECT_DOUBLE_EQ(evaluate("abs(-1)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("abs(-2*3)"), 6);
}

TEST(UnitTestEvaluator, testFunction_fabs)
{
  EXPECT_DOUBLE_EQ(evaluate("fabs(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("fabs(1)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("fabs(1.5)"),  1.5);
  EXPECT_DOUBLE_EQ(evaluate("fabs(-1)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("fabs(-2*3)"), 6);
}

TEST(UnitTestEvaluator, testFunction_max)
{
  EXPECT_DOUBLE_EQ(evaluate("max(1,2)"),          2);
  EXPECT_DOUBLE_EQ(evaluate("max(2,1)"),          2);
  EXPECT_DOUBLE_EQ(evaluate("max(1,2,3)"),        3);
  EXPECT_DOUBLE_EQ(evaluate("max(1,2,3,4)"),      4);
  EXPECT_DOUBLE_EQ(evaluate("max(-1,-2)"),       -1);
  EXPECT_DOUBLE_EQ(evaluate("max(-1,-2,-3)"),    -1);
  EXPECT_DOUBLE_EQ(evaluate("max(-1,-2,-3,-4)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("max(3+2,2+1)"),      5);
}

TEST(UnitTestEvaluator, testFunction_min)
{
  EXPECT_DOUBLE_EQ(evaluate("min(4,3)"),          3);
  EXPECT_DOUBLE_EQ(evaluate("min(3,4)"),          3);
  EXPECT_DOUBLE_EQ(evaluate("min(4,3,2)"),        2);
  EXPECT_DOUBLE_EQ(evaluate("min(4,3,2,1)"),      1);
  EXPECT_DOUBLE_EQ(evaluate("min(-1,-2)"),       -2);
  EXPECT_DOUBLE_EQ(evaluate("min(-1,-2,-3)"),    -3);
  EXPECT_DOUBLE_EQ(evaluate("min(-1,-2,-3,-4)"), -4);
  EXPECT_DOUBLE_EQ(evaluate("min(3+2,2+1)"),      3);
}

TEST(UnitTestEvaluator, testFunction_sign)
{
  EXPECT_DOUBLE_EQ(evaluate("sign(-10)"),       -1);
  EXPECT_DOUBLE_EQ(evaluate("sign(-1)"),        -1);
  EXPECT_DOUBLE_EQ(evaluate("sign(-0.5)"),      -1);
  EXPECT_DOUBLE_EQ(evaluate("sign(-0.000001)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("sign(-0)"),         1);
  EXPECT_DOUBLE_EQ(evaluate("sign(0)"),          1);
  EXPECT_DOUBLE_EQ(evaluate("sign(0.000001)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("sign(0.5)"),        1);
  EXPECT_DOUBLE_EQ(evaluate("sign(1)"),          1);
  EXPECT_DOUBLE_EQ(evaluate("sign(10)"),         1);
}

TEST(UnitTestEvaluator, testFunction_ipart)
{
  EXPECT_DOUBLE_EQ(evaluate("ipart(-9.25)"),     -9);
  EXPECT_DOUBLE_EQ(evaluate("ipart(-2.5)"),      -2);
  EXPECT_DOUBLE_EQ(evaluate("ipart(-1.0)"),      -1);
  EXPECT_DOUBLE_EQ(evaluate("ipart(-0.000001)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("ipart(0.0)"),        0);
  EXPECT_DOUBLE_EQ(evaluate("ipart(0.000001)"),   0);
  EXPECT_DOUBLE_EQ(evaluate("ipart(1.0)"),        1);
  EXPECT_DOUBLE_EQ(evaluate("ipart(2.5)"),        2);
  EXPECT_DOUBLE_EQ(evaluate("ipart(9.25)"),       9);
}

TEST(UnitTestEvaluator, testFunction_fpart)
{
  EXPECT_DOUBLE_EQ(evaluate("fpart(-9.25)"),     -0.25);
  EXPECT_DOUBLE_EQ(evaluate("fpart(-2.5)"),      -0.5);
  EXPECT_DOUBLE_EQ(evaluate("fpart(-1.0)"),       0.0);
  EXPECT_DOUBLE_EQ(evaluate("fpart(-0.000001)"), -0.000001);
  EXPECT_DOUBLE_EQ(evaluate("fpart(0.0)"),        0.0);
  EXPECT_DOUBLE_EQ(evaluate("fpart(0.000001)"),   0.000001);
  EXPECT_DOUBLE_EQ(evaluate("fpart(1.0)"),        0.0);
  EXPECT_DOUBLE_EQ(evaluate("fpart(2.5)"),        0.5);
  EXPECT_DOUBLE_EQ(evaluate("fpart(9.25)"),       0.25);
}

TEST(UnitTestEvaluator, testFunction_ceil)
{
  EXPECT_DOUBLE_EQ(evaluate("ceil(-1.000001)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("ceil(-1)"),        -1);
  EXPECT_DOUBLE_EQ(evaluate("ceil(-0.999999)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("ceil(-0.5)"),       0);
  EXPECT_DOUBLE_EQ(evaluate("ceil(-0.000001)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("ceil(-0)"),        -0);
  EXPECT_DOUBLE_EQ(evaluate("ceil(0)"),          0);
  EXPECT_DOUBLE_EQ(evaluate("ceil(0.000001)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("ceil(0.5)"),        1);
  EXPECT_DOUBLE_EQ(evaluate("ceil(0.999999)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("ceil(1)"),          1);
  EXPECT_DOUBLE_EQ(evaluate("ceil(1.000001)"),   2);
}

TEST(UnitTestEvaluator, testFunction_floor)
{
  EXPECT_DOUBLE_EQ(evaluate("floor(-1.000001)"), -2);
  EXPECT_DOUBLE_EQ(evaluate("floor(-1)"),        -1);
  EXPECT_DOUBLE_EQ(evaluate("floor(-0.999999)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("floor(-0.5)"),      -1);
  EXPECT_DOUBLE_EQ(evaluate("floor(-0.000001)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("floor(-0)"),        -0);
  EXPECT_DOUBLE_EQ(evaluate("floor(0)"),          0);
  EXPECT_DOUBLE_EQ(evaluate("floor(0.000001)"),   0);
  EXPECT_DOUBLE_EQ(evaluate("floor(0.5)"),        0);
  EXPECT_DOUBLE_EQ(evaluate("floor(0.999999)"),   0);
  EXPECT_DOUBLE_EQ(evaluate("floor(1)"),          1);
  EXPECT_DOUBLE_EQ(evaluate("floor(1.000001)"),   1);
}

TEST(UnitTestEvaluator, testFunction_mod)
{
  EXPECT_DOUBLE_EQ(evaluate("mod(-9, 3)"),   0);
  EXPECT_DOUBLE_EQ(evaluate("mod(-9, -4)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("mod(-9, 4)"),  -1);
  EXPECT_DOUBLE_EQ(evaluate("mod(9, -4)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("mod(9, 3)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("mod(9, 4)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("mod(9.1, 4)"),  1.1);
  EXPECT_DOUBLE_EQ(evaluate("mod(9, 4.5)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("mod(9, 3.5)"),  2);
  EXPECT_DOUBLE_EQ(evaluate("mod(9, 5)"),    4);
  EXPECT_DOUBLE_EQ(evaluate("mod(9, 9)"),    0);
}

TEST(UnitTestEvaluator, testFunction_fmod)
{
  EXPECT_DOUBLE_EQ(evaluate("fmod(-9, 3)"),   0);
  EXPECT_DOUBLE_EQ(evaluate("fmod(-9, -4)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("fmod(-9, 4)"),  -1);
  EXPECT_DOUBLE_EQ(evaluate("fmod(9, -4)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("fmod(9, 3)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("fmod(9, 4)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("fmod(9.1, 4)"),  1.1);
  EXPECT_DOUBLE_EQ(evaluate("fmod(9, 4.5)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("fmod(9, 3.5)"),  2);
  EXPECT_DOUBLE_EQ(evaluate("fmod(9, 5)"),    4);
  EXPECT_DOUBLE_EQ(evaluate("fmod(9, 9)"),    0);
}

TEST(UnitTestEvaluator, testFunction_pow)
{
  EXPECT_DOUBLE_EQ(evaluate("pow(0, 2)"),                 0);
  EXPECT_DOUBLE_EQ(evaluate("pow(2, 0)"),                 1);
  EXPECT_DOUBLE_EQ(evaluate("pow(3, 2)"),                 9);
  EXPECT_DOUBLE_EQ(evaluate("pow(3, 2.5)"),               std::pow(3, 2.5));
  EXPECT_DOUBLE_EQ(evaluate("pow(pow(3, 2), 2)"),         81);
  EXPECT_DOUBLE_EQ(evaluate("pow(pow(2, 2), pow(2, 2))"), 256);
  EXPECT_DOUBLE_EQ(evaluate("pow(2^2, 2^2)"),             256);
  EXPECT_DOUBLE_EQ(evaluate("pow(2, -1)"),                0.5);
}

TEST(UnitTestEvaluator, testFunction_sqrt)
{
  EXPECT_DOUBLE_EQ(evaluate("sqrt(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("sqrt(1)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("sqrt(4)"),    2);
  EXPECT_DOUBLE_EQ(evaluate("sqrt(9)"),    3);
  EXPECT_DOUBLE_EQ(evaluate("sqrt(2)"),    std::sqrt(2));
  EXPECT_DOUBLE_EQ(evaluate("sqrt(1.21)"), 1.1);
}

TEST(UnitTestEvaluator, testFunction_exp)
{
  EXPECT_DOUBLE_EQ(evaluate("exp(-2)"),  std::exp(-2));
  EXPECT_DOUBLE_EQ(evaluate("exp(-1)"),  std::exp(-1));
  EXPECT_DOUBLE_EQ(evaluate("exp(0)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("exp(1)"),   std::exp(1));
  EXPECT_DOUBLE_EQ(evaluate("exp(1.5)"), std::exp(1.5));
  EXPECT_DOUBLE_EQ(evaluate("exp(2)"),   std::exp(2));
}

TEST(UnitTestEvaluator, testFunction_ln)
{
  EXPECT_DOUBLE_EQ(evaluate("ln(1)"),        0);
  EXPECT_DOUBLE_EQ(evaluate("ln(0.5)"),      std::log(0.5));
  EXPECT_DOUBLE_EQ(evaluate("ln(exp(1))"),   1);
  EXPECT_DOUBLE_EQ(evaluate("ln(exp(1.5))"), 1.5);
  EXPECT_DOUBLE_EQ(evaluate("ln(exp(2))"),   2);
}

TEST(UnitTestEvaluator, testFunction_log)
{
  EXPECT_DOUBLE_EQ(evaluate("log(1)"),        0);
  EXPECT_DOUBLE_EQ(evaluate("log(0.5)"),      std::log(0.5));
  EXPECT_DOUBLE_EQ(evaluate("log(exp(1))"),   1);
  EXPECT_DOUBLE_EQ(evaluate("log(exp(1.5))"), 1.5);
  EXPECT_DOUBLE_EQ(evaluate("log(exp(2))"),   2);
}

TEST(UnitTestEvaluator, testFunction_log10)
{
  EXPECT_DOUBLE_EQ(evaluate("log10(0.001)"),   -3);
  EXPECT_DOUBLE_EQ(evaluate("log10(1)"),        0);
  EXPECT_DOUBLE_EQ(evaluate("log10(10)"),       1);
  EXPECT_DOUBLE_EQ(evaluate("log10(12)"),       std::log10(12));
  EXPECT_DOUBLE_EQ(evaluate("log10(1000)"),     3);
}

TEST(UnitTestEvaluator, testFunction_deg)
{
  EXPECT_DOUBLE_EQ(evaluate("deg(-TWO_PI)"), -360);
  EXPECT_DOUBLE_EQ(evaluate("deg(-PI)"),     -180);
  EXPECT_DOUBLE_EQ(evaluate("deg(-PI/2)"),   -90);
  EXPECT_DOUBLE_EQ(evaluate("deg(-PI/4)"),   -45);
  EXPECT_DOUBLE_EQ(evaluate("deg(0)"),        0);
  EXPECT_DOUBLE_EQ(evaluate("deg(PI/4)"),     45);
  EXPECT_DOUBLE_EQ(evaluate("deg(PI/2)"),     90);
  EXPECT_DOUBLE_EQ(evaluate("deg(PI)"),       180);
  EXPECT_DOUBLE_EQ(evaluate("deg(TWO_PI)"),   360);
}

TEST(UnitTestEvaluator, testFunction_rad)
{
  EXPECT_DOUBLE_EQ(evaluate("rad(-360)"), -stk::expreval::two_pi());
  EXPECT_DOUBLE_EQ(evaluate("rad(-180)"), -stk::expreval::pi());
  EXPECT_DOUBLE_EQ(evaluate("rad(-90)"),  -stk::expreval::pi()/2);
  EXPECT_DOUBLE_EQ(evaluate("rad(-45)"),  -stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("rad(0)"),     0);
  EXPECT_DOUBLE_EQ(evaluate("rad(45)"),    stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("rad(90)"),    stk::expreval::pi()/2);
  EXPECT_DOUBLE_EQ(evaluate("rad(180)"),   stk::expreval::pi());
  EXPECT_DOUBLE_EQ(evaluate("rad(360)"),   stk::expreval::two_pi());
}

TEST(UnitTestEvaluator, testFunction_sin)
{
  EXPECT_DOUBLE_EQ(evaluate("sin(0)"),       0);
  EXPECT_DOUBLE_EQ(evaluate("sin(PI/6)"),    0.5);
  EXPECT_DOUBLE_EQ(evaluate("sin(PI/4)"),    std::sqrt(2)/2);
  EXPECT_DOUBLE_EQ(evaluate("sin(PI/2)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("sin(PI)"),      std::sin(stk::expreval::pi()));
  EXPECT_DOUBLE_EQ(evaluate("sin(3*PI/2)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("sin(TWO_PI)"),  std::sin(stk::expreval::two_pi()));
}

TEST(UnitTestEvaluator, testFunction_cos)
{
  EXPECT_DOUBLE_EQ(evaluate("cos(0)"),       1);
  EXPECT_DOUBLE_EQ(evaluate("cos(PI/4)"),    std::sqrt(2)/2);
  EXPECT_DOUBLE_EQ(evaluate("cos(PI/3)"),    0.5);
  EXPECT_DOUBLE_EQ(evaluate("cos(PI/2)"),    std::cos(stk::expreval::pi()/2));
  EXPECT_DOUBLE_EQ(evaluate("cos(PI)"),     -1);
  EXPECT_DOUBLE_EQ(evaluate("cos(3*PI/2)"),  std::cos(3*stk::expreval::pi()/2));
  EXPECT_DOUBLE_EQ(evaluate("cos(TWO_PI)"),  1);
}

TEST(UnitTestEvaluator, testFunction_tan)
{
  EXPECT_DOUBLE_EQ(evaluate("tan(0)"),       0);
  EXPECT_DOUBLE_EQ(evaluate("tan(PI/3)"),    std::sqrt(3));
  EXPECT_DOUBLE_EQ(evaluate("tan(PI/4)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("tan(PI)"),      std::tan(stk::expreval::pi()));
  EXPECT_DOUBLE_EQ(evaluate("tan(3*PI/4)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("tan(TWO_PI)"),  std::tan(stk::expreval::two_pi()));
}

TEST(UnitTestEvaluator, testFunction_asin)
{
  EXPECT_DOUBLE_EQ(evaluate("asin(0)"),         0);
  EXPECT_DOUBLE_EQ(evaluate("asin(0.5)"),       stk::expreval::pi()/6);
  EXPECT_DOUBLE_EQ(evaluate("asin(sqrt(2)/2)"), stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("asin(1)"),         stk::expreval::pi()/2);
}

TEST(UnitTestEvaluator, testFunction_acos)
{
  EXPECT_DOUBLE_EQ(evaluate("acos(1)"),         0);
  EXPECT_DOUBLE_EQ(evaluate("acos(sqrt(2)/2)"), stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("acos(0.5)"),       stk::expreval::pi()/3);
  EXPECT_DOUBLE_EQ(evaluate("acos(0)"),         stk::expreval::pi()/2);
}

TEST(UnitTestEvaluator, testFunction_atan)
{
  EXPECT_DOUBLE_EQ(evaluate("atan(0)"),       0);
  EXPECT_DOUBLE_EQ(evaluate("atan(sqrt(3))"), stk::expreval::pi()/3);
  EXPECT_DOUBLE_EQ(evaluate("atan(1)"),       stk::expreval::pi()/4);
}

TEST(UnitTestEvaluator, testFunction_atan2)
{
  EXPECT_DOUBLE_EQ(evaluate("atan2(0, 1)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("atan2(1, 1)"),    stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("atan2(-1, 1)"),  -stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("atan2(1, -1)"),   3*stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("atan2(-1, -1)"), -3*stk::expreval::pi()/4);
}

TEST(UnitTestEvaluator, testFunction_sinh)
{
  EXPECT_DOUBLE_EQ(evaluate("sinh(-0.1)"), std::sinh(-0.1));
  EXPECT_DOUBLE_EQ(evaluate("sinh(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("sinh(0.1)"),  std::sinh(0.1));
  EXPECT_DOUBLE_EQ(evaluate("sinh(0.5)"),  std::sinh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("sinh(1)"),    std::sinh(1));
  EXPECT_DOUBLE_EQ(evaluate("sinh(2)"),    std::sinh(2));
  EXPECT_DOUBLE_EQ(evaluate("sinh(10)"),   std::sinh(10));
}

TEST(UnitTestEvaluator, testFunction_cosh)
{
  EXPECT_DOUBLE_EQ(evaluate("cosh(-0.1)"), std::cosh(-0.1));
  EXPECT_DOUBLE_EQ(evaluate("cosh(0)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("cosh(0.1)"),  std::cosh(0.1));
  EXPECT_DOUBLE_EQ(evaluate("cosh(0.5)"),  std::cosh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("cosh(1)"),    std::cosh(1));
  EXPECT_DOUBLE_EQ(evaluate("cosh(2)"),    std::cosh(2));
  EXPECT_DOUBLE_EQ(evaluate("cosh(10)"),   std::cosh(10));
}

TEST(UnitTestEvaluator, testFunction_tanh)
{
  EXPECT_DOUBLE_EQ(evaluate("tanh(-0.1)"), std::tanh(-0.1));
  EXPECT_DOUBLE_EQ(evaluate("tanh(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("tanh(0.1)"),  std::tanh(0.1));
  EXPECT_DOUBLE_EQ(evaluate("tanh(0.5)"),  std::tanh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("tanh(1)"),    std::tanh(1));
  EXPECT_DOUBLE_EQ(evaluate("tanh(2)"),    std::tanh(2));
  EXPECT_DOUBLE_EQ(evaluate("tanh(10)"),   std::tanh(10));
}

TEST(UnitTestEvaluator, testFunction_asinh)
{
  EXPECT_DOUBLE_EQ(evaluate("asinh(-0.1)"), std::asinh(-0.1));
  EXPECT_DOUBLE_EQ(evaluate("asinh(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("asinh(0.1)"),  std::asinh(0.1));
  EXPECT_DOUBLE_EQ(evaluate("asinh(0.5)"),  std::asinh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("asinh(1)"),    std::asinh(1));
  EXPECT_DOUBLE_EQ(evaluate("asinh(2)"),    std::asinh(2));
  EXPECT_DOUBLE_EQ(evaluate("asinh(10)"),   std::asinh(10));
}

TEST(UnitTestEvaluator, testFunction_acosh)
{
  EXPECT_DOUBLE_EQ(evaluate("acosh(1)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("acosh(2)"),    std::acosh(2));
  EXPECT_DOUBLE_EQ(evaluate("acosh(10)"),   std::acosh(10));
}

TEST(UnitTestEvaluator, testFunction_atanh)
{
  EXPECT_DOUBLE_EQ(evaluate("atanh(-0.1)"), std::atanh(-0.1));
  EXPECT_DOUBLE_EQ(evaluate("atanh(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("atanh(0.1)"),  std::atanh(0.1));
  EXPECT_DOUBLE_EQ(evaluate("atanh(0.5)"),  std::atanh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("atanh(1)"),    std::atanh(1));
}

TEST(UnitTestEvaluator, testFunction_erf)
{
  EXPECT_DOUBLE_EQ(evaluate("erf(-2)"),   std::erf(-2));
  EXPECT_DOUBLE_EQ(evaluate("erf(-1.5)"), std::erf(-1.5));
  EXPECT_DOUBLE_EQ(evaluate("erf(-1)"),   std::erf(-1));
  EXPECT_DOUBLE_EQ(evaluate("erf(-0)"),  -0);
  EXPECT_DOUBLE_EQ(evaluate("erf(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("erf(1)"),    std::erf(1));
  EXPECT_DOUBLE_EQ(evaluate("erf(1.5)"),  std::erf(1.5));
  EXPECT_DOUBLE_EQ(evaluate("erf(2)"),    std::erf(2));
}

TEST(UnitTestEvaluator, testFunction_erfc)
{
  EXPECT_DOUBLE_EQ(evaluate("erfc(-2)"),   std::erfc(-2));
  EXPECT_DOUBLE_EQ(evaluate("erfc(-1.5)"), std::erfc(-1.5));
  EXPECT_DOUBLE_EQ(evaluate("erfc(-1)"),   std::erfc(-1));
  EXPECT_DOUBLE_EQ(evaluate("erfc(0)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("erfc(1)"),    std::erfc(1));
  EXPECT_DOUBLE_EQ(evaluate("erfc(1.5)"),  std::erfc(1.5));
  EXPECT_DOUBLE_EQ(evaluate("erfc(2)"),    std::erfc(2));
}

TEST(UnitTestEvaluator, testFunction_poltorectx)
{
  EXPECT_DOUBLE_EQ(evaluate("poltorectx(-5, 0)"),    -5);
  EXPECT_DOUBLE_EQ(evaluate("poltorectx(5, 0)"),      5);
  EXPECT_DOUBLE_EQ(evaluate("poltorectx(5, PI/2)"),   5*std::cos(stk::expreval::pi()/2));
  EXPECT_DOUBLE_EQ(evaluate("poltorectx(5, PI)"),    -5);
  EXPECT_DOUBLE_EQ(evaluate("poltorectx(5, 3*PI/2)"), 5*std::cos(3*stk::expreval::pi()/2));
  EXPECT_DOUBLE_EQ(evaluate("poltorectx(0, 0)"),      0);
  EXPECT_DOUBLE_EQ(evaluate("poltorectx(0, PI)"),     0);
}

TEST(UnitTestEvaluator, testFunction_poltorecty)
{
  EXPECT_DOUBLE_EQ(evaluate("poltorecty(-5, PI/2)"),  -5);
  EXPECT_DOUBLE_EQ(evaluate("poltorecty(5, 0)"),       0);
  EXPECT_DOUBLE_EQ(evaluate("poltorecty(5, PI/2)"),    5);
  EXPECT_DOUBLE_EQ(evaluate("poltorecty(5, PI)"),      5*std::sin(stk::expreval::pi()));
  EXPECT_DOUBLE_EQ(evaluate("poltorecty(5, 3*PI/2)"), -5);
  EXPECT_DOUBLE_EQ(evaluate("poltorecty(0, PI/2)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("poltorecty(0, 3*PI/2)"),  0);
}

TEST(UnitTestEvaluator, testFunction_recttopolr)
{
  EXPECT_DOUBLE_EQ(evaluate("recttopolr(0, 0)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("recttopolr(1, 0)"),  1);
  EXPECT_DOUBLE_EQ(evaluate("recttopolr(0, 1)"),  1);
  EXPECT_DOUBLE_EQ(evaluate("recttopolr(-1, 0)"), 1);
  EXPECT_DOUBLE_EQ(evaluate("recttopolr(0, -1)"), 1);
  EXPECT_DOUBLE_EQ(evaluate("recttopolr(1, 1)"),  std::sqrt(2));
}

TEST(UnitTestEvaluator, testFunction_recttopola)
{
  EXPECT_DOUBLE_EQ(evaluate("recttopola(0, 0)"),   0);
  EXPECT_DOUBLE_EQ(evaluate("recttopola(1, 0)"),   0);
  EXPECT_DOUBLE_EQ(evaluate("recttopola(0, 1)"),   stk::expreval::pi()/2);
  EXPECT_DOUBLE_EQ(evaluate("recttopola(-1, 0)"),  stk::expreval::pi());
  EXPECT_DOUBLE_EQ(evaluate("recttopola(0, -1)"),  3*stk::expreval::pi()/2);
  EXPECT_DOUBLE_EQ(evaluate("recttopola(1, 1)"),   stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("recttopola(-1, 1)"),  3*stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("recttopola(-1, -1)"), 5*stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("recttopola(1, -1)"),  7*stk::expreval::pi()/4);
}

TEST(UnitTestEvaluator, testFunction_unit_step)
{
  EXPECT_DOUBLE_EQ(evaluate("unit_step(-0.5, 0, 1)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("unit_step(0, 0, 1)"),     1);
  EXPECT_DOUBLE_EQ(evaluate("unit_step(0.5, 0, 1)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("unit_step(1, 0, 1)"),     1);
  EXPECT_DOUBLE_EQ(evaluate("unit_step(1.5, 0, 1)"),   0);
}

TEST(UnitTestEvaluator, testFunction_cycloidal_ramp)
{
  EXPECT_DOUBLE_EQ(evaluate("cycloidal_ramp(-0.5, 0, 1)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("cycloidal_ramp(0, 0, 1)"),     0);
  EXPECT_DOUBLE_EQ(evaluate("cycloidal_ramp(0.25, 0, 1)"),  0.25-1/(2*stk::expreval::pi()));
  EXPECT_DOUBLE_EQ(evaluate("cycloidal_ramp(0.5, 0, 1)"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("cycloidal_ramp(0.75, 0, 1)"),  0.75+1/(2*stk::expreval::pi()));
  EXPECT_DOUBLE_EQ(evaluate("cycloidal_ramp(1, 0, 1)"),     1);
  EXPECT_DOUBLE_EQ(evaluate("cycloidal_ramp(1.5, 0, 1)"),   1);
}

TEST(UnitTestEvaluator, testFunction_cos_ramp3)
{
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(-0.5, 0, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(0, 0, 1)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1/3, 0, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(0.5, 0, 1)"),  0.5);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(2/3, 0, 1)"),  0.75);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1, 0, 1)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1.5, 0, 1)"),  1);
}

TEST(UnitTestEvaluator, testFunction_cos_ramp2)
{
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(-0.5, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(0, 1)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1/3, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(0.5, 1)"),  0.5);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(2/3, 1)"),  0.75);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1, 1)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1.5, 1)"),  1);
}

TEST(UnitTestEvaluator, testFunction_cos_ramp1)
{
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(-0.5)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1/3)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(0.5)"),  0.5);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(2/3)"),  0.75);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1.5)"),  1);
}

TEST(UnitTestEvaluator, testFunction_cosine_ramp3)
{
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(-0.5, 0, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(0, 0, 1)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1/3, 0, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(0.5, 0, 1)"),  0.5);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(2/3, 0, 1)"),  0.75);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1, 0, 1)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1.5, 0, 1)"),  1);
}

TEST(UnitTestEvaluator, testFunction_cosine_ramp2)
{
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(-0.5, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(0, 1)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1/3, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(0.5, 1)"),  0.5);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(2/3, 1)"),  0.75);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1, 1)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1.5, 1)"),  1);
}

TEST(UnitTestEvaluator, testFunction_cosine_ramp1)
{
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(-0.5)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(0)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1/3)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(0.5)"),  0.5);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(2/3)"),  0.75);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1.5)"),  1);
}

TEST(UnitTestEvaluator, testFunction_haversine_pulse)
{
  EXPECT_DOUBLE_EQ(evaluate("haversine_pulse(-0.5, 0, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("haversine_pulse(0, 0, 1)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("haversine_pulse(1/6, 0, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("haversine_pulse(0.5, 0, 1)"),  1);
  EXPECT_DOUBLE_EQ(evaluate("haversine_pulse(5/6, 0, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("haversine_pulse(1, 0, 1)"),    0);
  EXPECT_DOUBLE_EQ(evaluate("haversine_pulse(1.5, 0, 1)"),  0);
}

TEST(UnitTestEvaluator, testFunction_point2d)
{
  EXPECT_DOUBLE_EQ(evaluate("point2d(0, 0, 1, 1)"),   1);

  EXPECT_DOUBLE_EQ(evaluate("point2d(0.5, 0, 1, 1)"), 1);
  EXPECT_DOUBLE_EQ(evaluate("point2d(5/6, 0, 1, 1)"), 0.75);
  EXPECT_DOUBLE_EQ(evaluate("point2d(1, 0, 1, 1)"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("point2d(7/6, 0, 1, 1)"), 0.25);
  EXPECT_DOUBLE_EQ(evaluate("point2d(1.5, 0, 1, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("point2d(2, 0, 1, 1)"),   0);

  EXPECT_DOUBLE_EQ(evaluate("point2d(0, -0.5, 1, 1)"), 1);
  EXPECT_DOUBLE_EQ(evaluate("point2d(0, -5/6, 1, 1)"), 0.75);
  EXPECT_DOUBLE_EQ(evaluate("point2d(0, -1, 1, 1)"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("point2d(0, -7/6, 1, 1)"), 0.25);
  EXPECT_DOUBLE_EQ(evaluate("point2d(0, -1.5, 1, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("point2d(0, -2, 1, 1)"),   0);
}

TEST(UnitTestEvaluator, testFunction_point3d)
{
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, 0, 0, 1, 1)"),   1);

  EXPECT_DOUBLE_EQ(evaluate("point3d(0.5, 0, 0, 1, 1)"), 1);
  EXPECT_DOUBLE_EQ(evaluate("point3d(5/6, 0, 0, 1, 1)"), 0.75);
  EXPECT_DOUBLE_EQ(evaluate("point3d(1, 0, 0, 1, 1)"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("point3d(7/6, 0, 0, 1, 1)"), 0.25);
  EXPECT_DOUBLE_EQ(evaluate("point3d(1.5, 0, 0, 1, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("point3d(2, 0, 0, 1, 1)"),   0);

  EXPECT_DOUBLE_EQ(evaluate("point3d(0, -0.5, 0, 1, 1)"), 1);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, -5/6, 0, 1, 1)"), 0.75);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, -1, 0, 1, 1)"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, -7/6, 0, 1, 1)"), 0.25);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, -1.5, 0, 1, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, -2, 0, 1, 1)"),   0);

  EXPECT_DOUBLE_EQ(evaluate("point3d(0, 0, -0.5, 1, 1)"), 1);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, 0, -5/6, 1, 1)"), 0.75);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, 0, -1, 1, 1)"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, 0, -7/6, 1, 1)"), 0.25);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, 0, -1.5, 1, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, 0, -2, 1, 1)"),   0);
}

TEST(UnitTestEvaluator, testFunction_exponential_pdf)
{
  EXPECT_DOUBLE_EQ(evaluate("exponential_pdf(-1, 1)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("exponential_pdf(0, 1)"),  1);
  EXPECT_DOUBLE_EQ(evaluate("exponential_pdf(1, 1)"),  1/std::exp(1));
}

TEST(UnitTestEvaluator, testFunction_log_uniform)
{
  EXPECT_DOUBLE_EQ(evaluate("log_uniform_pdf(0.5, 1, E)"), 0);
  EXPECT_DOUBLE_EQ(evaluate("log_uniform_pdf(1, 1, E)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("log_uniform_pdf(2, 1, E)"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("log_uniform_pdf(E, 1, E)"),   1/std::exp(1));
  EXPECT_DOUBLE_EQ(evaluate("log_uniform_pdf(3, 1, E)"),   0);
}

double reference_normal_pdf(double x, double mu, double sigma) {
  return std::exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))/std::sqrt(2.0*stk::expreval::pi()*sigma*sigma);
}

TEST(UnitTestEvaluator, testFunction_normal_pdf)
{
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(-0.5, 1, 0.5)"), reference_normal_pdf(-0.5, 1, 0.5));
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(0, 1, 0.5)"),    reference_normal_pdf(0, 1, 0.5));
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(0.5, 1, 0.5)"),  reference_normal_pdf(0.5, 1, 0.5));
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(0.75, 1, 0.5)"), reference_normal_pdf(0.75, 1, 0.5));
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(1, 1, 0.5)"),    reference_normal_pdf(1, 1, 0.5));
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(1.25, 1, 0.5)"), reference_normal_pdf(1.25, 1, 0.5));
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(1.5, 1, 0.5)"),  reference_normal_pdf(1.5, 1, 0.5));
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(2, 1, 0.5)"),    reference_normal_pdf(2, 1, 0.5));
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(2.5, 1, 0.5)"),  reference_normal_pdf(2.5, 1, 0.5));
}

double reference_weibull_pdf(double x, double k, double lambda) {
  return (x >= 0) ? (k/lambda)*std::pow(x/lambda, k-1)*std::exp(-std::pow(x/lambda, k)) : 0;
}

TEST(UnitTestEvaluator, testFunction_weibull_pdf)
{
  EXPECT_DOUBLE_EQ(evaluate("weibull_pdf(-1, 5, 1)"),   0);
  EXPECT_DOUBLE_EQ(evaluate("weibull_pdf(0, 5, 1)"),    reference_weibull_pdf(0, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("weibull_pdf(0.5, 5, 1)"),  reference_weibull_pdf(0.5, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("weibull_pdf(0.75, 5, 1)"), reference_weibull_pdf(0.75, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("weibull_pdf(1, 5, 1)"),    reference_weibull_pdf(1, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("weibull_pdf(1.25, 5, 1)"), reference_weibull_pdf(1.25, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("weibull_pdf(1.5, 5, 1)"),  reference_weibull_pdf(1.5, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("weibull_pdf(2, 5, 1)"),    reference_weibull_pdf(2, 5, 1));
}

double reference_gamma_pdf(double x, double k, double theta) {
  return (x >= 0) ? 1/(std::tgamma(k)*std::pow(theta, k))*std::pow(x, k-1)*std::exp(-x/theta) : 0;
}

TEST(UnitTestEvaluator, testFunction_gamma_pdf)
{
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(-1, 5, 1)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(0, 5, 1)"),   reference_gamma_pdf(0, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(0.5, 5, 1)"), reference_gamma_pdf(0.5, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(1, 5, 1)"),   reference_gamma_pdf(1, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(2, 5, 1)"),   reference_gamma_pdf(2, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(3, 5, 1)"),   reference_gamma_pdf(3, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(3.5, 5, 1)"), reference_gamma_pdf(3.5, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(4, 5, 1)"),   reference_gamma_pdf(4, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(4.5, 5, 1)"), reference_gamma_pdf(4.5, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(5, 5, 1)"),   reference_gamma_pdf(5, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(6, 5, 1)"),   reference_gamma_pdf(6, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(8, 5, 1)"),   reference_gamma_pdf(8, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(10, 5, 1)"),  reference_gamma_pdf(10, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(12, 5, 1)"),  reference_gamma_pdf(12, 5, 1));
}

void testRandom(const char * expression)
{
  const int NUM_SAMPLES = 10000;
  const double EXPECTED_MEAN = 0.5;
  const double EXPECTED_SIGMA = 1./std::sqrt(12.);
  std::vector<int> bins(10, 0);
  double mean = 0;
  double sigma = 0;

  for (int i = 0; i < NUM_SAMPLES; ++i) {
    const double result = evaluate(expression);
    const int binIndex = static_cast<int>(result*10);
    bins[binIndex] += 1;

    mean += result;
    sigma += (result-EXPECTED_MEAN)*(result-EXPECTED_MEAN);
  }

  mean /= NUM_SAMPLES;
  sigma = std::sqrt(sigma/NUM_SAMPLES);

  EXPECT_NEAR(mean, EXPECTED_MEAN, 0.01);
  EXPECT_NEAR(sigma, EXPECTED_SIGMA, 0.005);

  const int maxN = *std::max_element(bins.begin(), bins.end());
  const int minN = *std::min_element(bins.begin(), bins.end());

  EXPECT_NEAR(maxN, NUM_SAMPLES/10, 100);
  EXPECT_NEAR(minN, NUM_SAMPLES/10, 100);
}

TEST(UnitTestEvaluator, testFunction_rand)
{
  testRandom("rand()");
}

TEST(UnitTestEvaluator, testFunction_srand_repeatability)
{
  std::vector<double> result(10);
  evaluate("srand(123.)");
  for (unsigned i = 0; i < result.size(); ++i) {
    result[i] = evaluate("rand()");
  }

  evaluate("srand(123.)");
  for (unsigned i = 0; i < result.size(); ++i) {
    EXPECT_DOUBLE_EQ(evaluate("rand()"), result[i]);
  }
}

TEST(UnitTestEvaluator, testFunction_random)
{
  testRandom("random()");
}

TEST(UnitTestEvaluator, testFunction_random1_repeatability)
{
  std::vector<double> result(10);
  evaluate("random(123.)");
  for (unsigned i = 0; i < result.size(); ++i) {
    result[i] = evaluate("random()");
  }

  evaluate("random(123.)");
  for (unsigned i = 0; i < result.size(); ++i) {
    EXPECT_DOUBLE_EQ(evaluate("random()"), result[i]);
  }
}

TEST(UnitTestEvaluator, testFunction_ts_random_repeatability)
{
  std::vector<double> result(10);
  double time = 0;
  for (unsigned i = 0; i < result.size(); ++i) {
    result[i] = evaluate("ts_random(t, 1.0, 2.0, 3.0)", {{"t", time}});
    time += 0.1;
  }

  time = 0;
  for (unsigned i = 0; i < result.size(); ++i) {
    EXPECT_DOUBLE_EQ(evaluate("ts_random(t, 1.0, 2.0, 3.0)", {{"t", time}}), result[i]);
    time += 0.1;
  }
}

TEST(UnitTestEvaluator, testFunction_ts_normal_repeatability)
{
  std::vector<double> result(10);
  double time = 0;
  for (unsigned i = 0; i < result.size(); ++i) {
    result[i] = evaluate("ts_normal(t, 1.0, 2.0, 3.0, 1.0, 0.5, 0, 2)", {{"t", time}});
    time += 0.1;
  }

  time = 0;
  for (unsigned i = 0; i < result.size(); ++i) {
    EXPECT_DOUBLE_EQ(evaluate("ts_normal(t, 1.0, 2.0, 3.0, 1.0, 0.5, 0, 2)", {{"t", time}}), result[i]);
    time += 0.1;
  }
}

TEST(UnitTestEvaluator, testFunction_ts_normal_clipping)
{
  const size_t NUM_SAMPLES = 10000;
  const double lowerBound = 0;
  const double upperBound = 2;

  double time = 0;
  for (unsigned i = 0; i < NUM_SAMPLES; ++i) {
    const double result = evaluate("ts_normal(t, 1.0, 2.0, 3.0, 1.0, 0.5, " + std::to_string(lowerBound) + ", " +
                                   std::to_string(upperBound) + ")", {{"t", time}});
    time += 0.1;
    EXPECT_GE(result, lowerBound);
    EXPECT_LE(result, upperBound);
  }
}

TEST(UnitTestEvaluator, testFunction_time)
{
  EXPECT_NEAR(evaluate("time()"), std::time(nullptr), 1.1);
}


} // namespace <unnamed>
