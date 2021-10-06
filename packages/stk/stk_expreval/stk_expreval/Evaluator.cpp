/*
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
 */

/*
  Loosely based on:  (Less and less so each checkin, practically none)
  File: Eval.c
  Auth: Brian Allen Vanderburg II
  Date: Wednesday, April 30, 2003
  Desc: Evaluation routines for the Eval library
*/


#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <stdexcept>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <stk_util/util/ReportHandler.hpp>

#include <math.h>
#include <time.h>

#include <stk_expreval/Evaluator.hpp>
#include <stk_expreval/Lexer.hpp>

using std::cout;
using std::endl;

namespace {
struct expression_evaluation_exception : public virtual std::exception
{
  virtual const char* what() const throw() {
    return "Error evaluating expressions";
  }

};

struct expression_undefined_exception : public virtual std::exception
{
  virtual const char* what() const throw() {
    static std::string rtnMsg;
    rtnMsg = "Found undefined function with name: " + m_msg + " and " + std::to_string(m_numArgs)  + " argument(s)";
    return rtnMsg.c_str();
  }

  expression_undefined_exception(const char *msg, std::size_t numArgs) : m_msg(msg), m_numArgs(numArgs) {}

  std::string m_msg;
  std::size_t m_numArgs = 0;
};
}

namespace stk {
namespace expreval {

class Node
{
  
public:
  enum { MAXIMUM_NUMBER_OF_OVERLOADED_FUNCTION_NAMES = 5 };
  enum { MAXIMUM_FUNCTION_NAME_LENGTH = 32 };

  explicit Node(Opcode opcode, Eval* owner)
    : m_opcode(opcode),
      m_left(nullptr),
      m_right(nullptr),
      m_other(nullptr),
      m_owner(owner)
  {
      m_data.function.undefinedFunction = false;
      for(unsigned i=0; i<MAXIMUM_NUMBER_OF_OVERLOADED_FUNCTION_NAMES; ++i)
          m_data.function.function[i] = nullptr;
  }

private:
  explicit Node(const Node &);
  Node &operator=(const Node &);

public:
  ~Node()
  {}

  double eval() const;

  const Opcode m_opcode;

  union _data
  {
    struct _constant
    {
      double value;
    } constant;

    struct _variable
    {
      Variable* variable;
    } variable;

    struct _function
    {
      CFunctionBase* function[MAXIMUM_NUMBER_OF_OVERLOADED_FUNCTION_NAMES];
      bool undefinedFunction;
      char functionName[MAXIMUM_FUNCTION_NAME_LENGTH];
      FunctionType functionType;
    } function;
  } m_data;

  int m_currentNodeIndex;
  Node* m_left;
  Node* m_right;
  Node* m_other;
  Eval* m_owner;
};

NgpNode::NgpNode()
  : m_opcode(OPCODE_UNDEFINED),
    m_data{{0.0}},
    m_currentNodeIndex(-1),
    m_leftNodeIndex(-1),
    m_rightNodeIndex(-1),
    m_otherNodeIndex(-1)
{}

NgpNode::NgpNode(const Node& node)
  : m_opcode(node.m_opcode),
    m_currentNodeIndex(node.m_currentNodeIndex),
    m_leftNodeIndex((node.m_left != nullptr) ? node.m_left->m_currentNodeIndex : -1),
    m_rightNodeIndex((node.m_right != nullptr) ? node.m_right->m_currentNodeIndex : -1),
    m_otherNodeIndex((node.m_other != nullptr) ? node.m_other->m_currentNodeIndex : -1)
{
  if (m_opcode == OPCODE_CONSTANT) {
    m_data.constant.value = node.m_data.constant.value;
  } 
  else if (m_opcode == OPCODE_RVALUE) {
    m_data.variable.variableIndex = node.m_data.variable.variable->get_index(); 
  }
  else if (m_opcode == OPCODE_ASSIGN) {
    m_data.variable.variableIndex = node.m_data.variable.variable->get_index(); 
    m_data.variable.variableType = node.m_data.variable.variable->get_type();
    m_data.variable.variableSize = node.m_data.variable.variable->getLength(); 
  }
  else if (m_opcode == OPCODE_FUNCTION) {
    m_data.function.functionType = node.m_data.function.functionType;
  }
}

KOKKOS_FUNCTION
double
NgpNode::evaluate_function(int argumentCount, double* arguments) const
{
  switch (m_data.function.functionType) {
    case FunctionType::ABS :
      if (argumentCount == 1) {
        return std::fabs(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for abs function");  
      break;
    case FunctionType::MAX :
      if (argumentCount == 2) {
        return max_2(arguments[0], arguments[1]);
      }
      else if (argumentCount == 3) {
        return max_3(arguments[0], arguments[1], arguments[2]);
      }
      else if (argumentCount == 4) {
        return max_4(arguments[0], arguments[1], arguments[2], arguments[3]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for max function");
      break;
    case FunctionType::MIN :
      if (argumentCount == 2) {
        return min_2(arguments[0], arguments[1]);
      }
      else if (argumentCount == 3) {
        return min_3(arguments[0], arguments[1], arguments[2]);
      }
      else if (argumentCount == 4) {
        return min_4(arguments[0], arguments[1], arguments[2], arguments[3]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for min function");
      break;
    case FunctionType::SIGN :
      if (argumentCount == 1) {
        return sign(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for sign function");
      break;
    case FunctionType::IPART :
      if (argumentCount == 1) {
        return ipart(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for ipart function");
      break;
    case FunctionType::FPART :
      if (argumentCount == 1) {
        return fpart(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for fpart function");
      break;
    case FunctionType::CEIL :
      if (argumentCount == 1) {
        return std::ceil(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for ceil function");
      break;
    case FunctionType::FLOOR :
      if (argumentCount == 1) {
        return std::floor(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for floor function");
      break;
    case FunctionType::MOD :
      if (argumentCount == 2) {
        return std::fmod(arguments[0], arguments[1]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for fmod or mod function");
      break;
    case FunctionType::POW :
      if (argumentCount == 2) {
        return std::pow(arguments[0], arguments[1]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for pow function");
      break;
    case FunctionType::SQRT :
      if (argumentCount == 1) {
        return std::sqrt(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for sqrt function");
      break;
    case FunctionType::EXP :
      if (argumentCount == 1) {
        return std::exp(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for exp function");
      break;
    case FunctionType::LN :
      if (argumentCount == 1) {
        return std::log(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for ln or log function");
      break;
    case FunctionType::LOG10 :
      if (argumentCount == 1) {
        return std::log10(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for log10 function");
      break;
    case FunctionType::DEG :
      if (argumentCount == 1) {
        return radian_to_degree()*arguments[0];
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for deg function");
      break;
    case FunctionType::RAD :
      if (argumentCount == 1) {
        return degree_to_radian()*arguments[0];
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for rad function");
      break;
    case FunctionType::SIN :
      if (argumentCount == 1) {
        return std::sin(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for sin function");
      break;
    case FunctionType::COS :
      if (argumentCount == 1) {
        return std::cos(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for cos function");
      break;
    case FunctionType::TAN :
      if (argumentCount == 1) {
        return std::tan(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for tan function");
      break;
    case FunctionType::ASIN :
      if (argumentCount == 1) {
        return std::asin(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for asin function");
      break;
    case FunctionType::ACOS :
      if (argumentCount == 1) {
        return std::acos(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for acos function");
      break;
    case FunctionType::ATAN :
      if (argumentCount == 1) {
        return std::atan(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for atan function");
      break;
    case FunctionType::ATAN2 :
      if (argumentCount == 2) {
        return std::atan2(arguments[0], arguments[1]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for atan2 function");
      break;
    case FunctionType::SINH :
      if (argumentCount == 1) {
        return std::sinh(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for sinh function");
      break;
    case FunctionType::COSH :
      if (argumentCount == 1) {
        return std::cosh(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arugments for cosh function");
      break;
    case FunctionType::TANH :
      if (argumentCount == 1) {
        return std::tanh(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arugments for tanh function");
      break;
    case FunctionType::ASINH :
      if (argumentCount == 1) {
        return std::asinh(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for asinh function");
      break;
    case FunctionType::ACOSH :
      if (argumentCount == 1) {
        return std::acosh(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for acosh function");
      break;
    case FunctionType::ATANH :
      if (argumentCount == 1) {
        return std::atanh(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for atanh function");
      break;
    case FunctionType::ERF :
      if (argumentCount == 1) {
        return std::erf(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for erf function");
      break;
    case FunctionType::ERFC :
      if (argumentCount == 1) {
        return std::erfc(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for erfc function");
      break;
    case FunctionType::POLTORECTX :
      if (argumentCount == 2) {
        return poltorectx(arguments[0], arguments[1]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for poltorectx function");
      break;
    case FunctionType::POLTORECTY :
      if (argumentCount == 2) {
        return poltorecty(arguments[0], arguments[1]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for poltorecty function");
      break;
    case FunctionType::RECTTOPOLR :
      if (argumentCount == 2) {
        return recttopolr(arguments[0], arguments[1]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for recttopolr function");
      break;
    case FunctionType::RECTTOPOLA :
      if (argumentCount == 2) {
        return recttopola(arguments[0], arguments[1]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for recttopola function");
      break;
    case FunctionType::UNIT_STEP :
      if (argumentCount == 3) {
        return unit_step3(arguments[0], arguments[1], arguments[2]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for unit_step function");
      break;
    case FunctionType::CYCLOIDAL_RAMP :
      if (argumentCount == 3) {
        return cycloidal_ramp(arguments[0], arguments[1], arguments[2]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for cycloidal_ramp function");
      break;
    case FunctionType::COS_RAMP :
      if (argumentCount == 1) {
        return cosine_ramp1(arguments[0]);
      }
      else if (argumentCount == 2) {
        return cosine_ramp2(arguments[0], arguments[1]);
      }
      else if (argumentCount == 3) {
        return cosine_ramp3(arguments[0], arguments[1], arguments[2]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for cos_ramp or cosine_ramp function");
      break;
    case FunctionType::HAVERSINE_PULSE :
      if (argumentCount == 3) {
        return haversine_pulse(arguments[0], arguments[1], arguments[2]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for haversine_pulse function");
      break;
    case FunctionType::POINT2D :
      if (argumentCount == 4) {
        return point_2(arguments[0], arguments[1], arguments[2], arguments[3]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for pulse_2 function");
      break;
    case FunctionType::POINT3D :
      if (argumentCount == 5) {
        return point_3(arguments[0], arguments[1], arguments[2], arguments[3], arguments[4]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for pulse_3 function");
      break;
    case FunctionType::EXPONENTIAL_PDF :
      if (argumentCount == 2) {
        return exponential_pdf(arguments[0], arguments[1]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for exponential_pdf function");
      break;
    case FunctionType::LOG_UNIFORM_PDF :
      if (argumentCount == 3) {
        return log_uniform_pdf(arguments[0], arguments[1], arguments[2]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for log_uniform_pdf function");
      break;
    case FunctionType::NORMAL_PDF :
      if (argumentCount == 3) {
        return normal_pdf(arguments[0], arguments[1], arguments[2]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for normal_pdf function");
      break;
    case FunctionType::WEIBULL_PDF :
      if (argumentCount == 3) {
        return weibull_pdf(arguments[0], arguments[1], arguments[2]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for weibull_pdf function");
      break;
    case FunctionType::GAMMA_PDF :
      if (argumentCount == 3) {
        return gamma_pdf(arguments[0], arguments[1], arguments[2]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for gamma_pdf function");
      break;
    case FunctionType::RAND :
      if (argumentCount == 0) {
        return real_rand();
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for rand function");
      break;
    case FunctionType::SRAND :
      if (argumentCount == 1) {
        return real_srand(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for srand function");
      break;
    case FunctionType::RANDOM :
      if (argumentCount == 0) {
        return random0();
      }
      else if (argumentCount == 1) {
        return random1(arguments[0]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for random function");
      break;
    case FunctionType::TS_RANDOM :
      if (argumentCount == 4) {
        return time_space_random(arguments[0], arguments[1], arguments[2], arguments[3]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for ts_random function");
      break;
    case FunctionType::TS_NORMAL :
      if (argumentCount == 8) {
        return time_space_normal(arguments[0], arguments[1], arguments[2], arguments[3], arguments[4],
                                 arguments[5], arguments[6], arguments[7]);
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for ts_normal function");
      break;
    case FunctionType::TIME :
      if (argumentCount == 0) {
        return current_time();
      }
      NGP_ThrowErrorMsg("Incorrect number of arguments for time function");
      break;
    case FunctionType::UNDEFINED :
      NGP_ThrowErrorMsg("Undefined function type");
      break;
    default :
      break;
  }
  return 0.0;
}

double
Node::eval() const
{
  switch (m_opcode) {
  case OPCODE_STATEMENT:
    {
      double value = 0.0;
      for (const Node *statement = this; statement; statement = statement->m_right)
        value = statement->m_left->eval();
      return value;
    }

  case OPCODE_CONSTANT:
    return m_data.constant.value;

  case OPCODE_RVALUE:
    /* Directly access the variable */
    if (m_left) {
      return m_data.variable.variable->getArrayValue(m_left->eval(), m_owner->getArrayOffsetType());
    } else {
      return m_data.variable.variable->getValue();
    }

  case OPCODE_MULTIPLY:
    return m_left->eval()*m_right->eval();

  case OPCODE_EXPONENIATION:
    return std::pow(m_left->eval(),m_right->eval());

  case OPCODE_DIVIDE:
    return m_left->eval()/m_right->eval();

  case OPCODE_MODULUS:
    return std::fmod(m_left->eval(), m_right->eval());

  case OPCODE_ADD:
    return m_left->eval() + m_right->eval();

  case OPCODE_SUBTRACT:
    return m_left->eval() - m_right->eval();

  case OPCODE_EQUAL:
    return m_left->eval() == m_right->eval() ? s_true : s_false;

  case OPCODE_NOT_EQUAL:
    return m_left->eval() != m_right->eval() ? s_true : s_false;

  case OPCODE_LESS:
    return m_left->eval() < m_right->eval() ? s_true : s_false;

  case OPCODE_GREATER:
    return m_left->eval() > m_right->eval() ? s_true : s_false;

  case OPCODE_LESS_EQUAL:
    return m_left->eval() <= m_right->eval() ? s_true : s_false;

  case OPCODE_GREATER_EQUAL:
    return m_left->eval() >= m_right->eval() ? s_true : s_false;

  case OPCODE_LOGICAL_AND: {
    double left = m_left->eval();
    double right = m_right->eval();
    return  (left != s_false) && (right != s_false) ? s_true : s_false;
  }
  case OPCODE_LOGICAL_OR: {
    double left = m_left->eval();
    double right = m_right->eval();
    return (left != s_false) || (right != s_false) ? s_true : s_false;
  }
  case OPCODE_TERNARY:
    return m_left->eval() != s_false ? m_right->eval() : m_other->eval();

  case OPCODE_UNARY_MINUS:
    return -m_right->eval();

  case OPCODE_UNARY_NOT:
    return m_right->eval() == s_false ? s_true : s_false;

  case OPCODE_ASSIGN:
    if (m_left)
      return m_data.variable.variable->getArrayValue(m_left->eval(),  m_owner->getArrayOffsetType()) = m_right->eval();
    else {
      *m_data.variable.variable = m_right->eval();
      return m_data.variable.variable->getValue();
    }

  case OPCODE_FUNCTION:
    {
      double argv[20];

      int argc = 0;
      for (Node *arg = m_right; arg; arg = arg->m_right)
      {
        argv[argc++] = arg->m_left->eval();
      }

      for(unsigned int i=0; i<Node::MAXIMUM_NUMBER_OF_OVERLOADED_FUNCTION_NAMES; ++i)
      {
         if(nullptr == m_data.function.function[i])
         {
             throw expression_undefined_exception(m_data.function.functionName, argc);
         }

        // Linear search to match the function name and number of arguments.
        if( m_data.function.function[i]->getArgCount() == argc) {
          return (*m_data.function.function[i])(argc, argv);
        }
      }
      throw expression_evaluation_exception();
    }

  default: // Unknown opcode
    throw expression_evaluation_exception();
  }
}

namespace Parser {

Node *parseStatements(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseStatement(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseExpression(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseAssign(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator assign, LexemVector::const_iterator to);
Node *parseTerm(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator term, LexemVector::const_iterator to);
Node *parseFactor(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator factor, LexemVector::const_iterator to);
Node *parseRelation(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator factor, LexemVector::const_iterator to);
Node *parseLogical(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator factor, LexemVector::const_iterator to);
Node *parseUnary(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator unary, LexemVector::const_iterator to);
Node *parseTernary(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator question, LexemVector::const_iterator colon, LexemVector::const_iterator to);
Node *parseFunction(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator lparen, LexemVector::const_iterator rparen, LexemVector::const_iterator to);
Node *parseFunctionArg(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseRValue(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseIndex(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator lbrack, LexemVector::const_iterator rbrack, LexemVector::const_iterator to);

// Parser productions

Node *
parseStatements(
  Eval &eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator to)
{
  if ((*from).getToken() == TOKEN_END)
    return nullptr;

  if ((*from).getToken() == TOKEN_SEMI)
    return parseStatements(eval, from + 1, to);

  // Technically, there should be no check for TOKEN_END, but we allow a missing final TOKEN_SEMI
  LexemVector::const_iterator it;
  for (it = from; (*it).getToken() != TOKEN_SEMI && (*it).getToken() != TOKEN_END; ++it)
    ;

  Node *statement = eval.newNode(OPCODE_STATEMENT);
  statement->m_left = parseStatement(eval, from, it);

  // Technically, there should be no check for TOKEN_END, but we allow a missing final TOKEN_SEMI
  if ((*it).getToken() != TOKEN_END)
    statement->m_right = parseStatements(eval, it + 1, to);

  return statement;
}


Node *
parseStatement(
  Eval &eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator to)
{
  return parseExpression(eval, from, to);
}


Node *
parseExpression(
  Eval &eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator to)
{
  int paren_level = 0;                                  // Paren level
  int brack_level = 0;                                  // Brack level
  LexemVector::const_iterator lparen_open_it = to;      // First open paren
  LexemVector::const_iterator lparen_close_it = to;     // Corresponding close paren
  LexemVector::const_iterator lbrack_open_it = to;      // First open bracket
  LexemVector::const_iterator lbrack_close_it = to;     // Corresponding close brack
  LexemVector::const_iterator assign_it = to;           // First = at paren_level 0 for assignment
  LexemVector::const_iterator term_it = to;             // Last + or - at paren_level 0 for adding or subtracting
  LexemVector::const_iterator factor_it = to;           // Last * or / at paren_level 0 for multiplying or dividing
  LexemVector::const_iterator expon_it = to;            // Last ^ for exponenation


  LexemVector::const_iterator relation_it = to;         // Last relational at paren_level 0 for relational operator
  LexemVector::const_iterator logical_it = to;          // Last logical at paren_level 0 for logical operator
  LexemVector::const_iterator question_it = to;         // Last ternary at paren_level 0 for ternary operator
  LexemVector::const_iterator colon_it = to;
  LexemVector::const_iterator unary_it = to;            // First +,- at plevel 0 for positive,negative
  LexemVector::const_iterator last_unary_it = to;       // Last +,- found at plevel for for positive,negative

  // Scan the expression for the instances of the above tokens
  for (LexemVector::const_iterator it = from; it != to; ++it) {
    switch((*it).getToken()) {
    case TOKEN_LPAREN:
      if (paren_level == 0 && lparen_open_it == to
         && brack_level == 0 && lbrack_open_it == to)
      lparen_open_it = it;
      paren_level++;
      break;

    case TOKEN_RPAREN:
      paren_level--;

      if (paren_level == 0 && lparen_close_it == to
       && brack_level == 0 && lbrack_close_it == to)
        lparen_close_it = it;

      if (paren_level < 0) {
        throw std::runtime_error("mismatched parenthesis");
      }
      break;

    case TOKEN_LBRACK:
      if (paren_level == 0 && lparen_open_it == to
          && brack_level == 0 && lbrack_open_it == to)
        lbrack_open_it = it;
      brack_level++;
      break;

    case TOKEN_RBRACK:
      brack_level--;

      if (paren_level == 0 && lparen_close_it == to
          && brack_level == 0 && lbrack_close_it == to)
        lbrack_close_it = it;

      if (brack_level < 0)
        throw std::runtime_error("mismatched bracket");
      break;

    case TOKEN_ASSIGN:
      if (paren_level == 0 && assign_it == to)
        assign_it = it;
      break;

    case TOKEN_QUESTION:
      if (paren_level == 0 && question_it == to)
        question_it = it;
      break;

    case TOKEN_COLON:
      if (paren_level == 0) // && colon_it == to)
        colon_it = it;
      break;

    case TOKEN_EXPONENTIATION:
      if (paren_level == 0) // && expon_it == to)
        expon_it = it;
      break;


    case TOKEN_MULTIPLY:
    case TOKEN_DIVIDE:
    case TOKEN_PERCENT:
      if (paren_level == 0) // && factor_it == to)
        factor_it = it;
      break;

    case TOKEN_EQUAL:
    case TOKEN_NOT_EQUAL:
    case TOKEN_LESS:
    case TOKEN_GREATER:
    case TOKEN_LESS_EQUAL:
    case TOKEN_GREATER_EQUAL:
      if (paren_level == 0 && relation_it == to)
        relation_it = it;
      break;

    case TOKEN_LOGICAL_AND:
    case TOKEN_LOGICAL_OR:
      if (paren_level == 0 && logical_it == to)
        logical_it = it;
      break;

    case TOKEN_PLUS:
    case TOKEN_MINUS:
      if (paren_level == 0) {
        // After any of these, we are a unary operator, not a term
        if (it == from || it == assign_it + 1
            || it == term_it + 1 || it == factor_it + 1
            || it == last_unary_it + 1 || it == expon_it + 1)
        { // Unary operator
          if (unary_it == to) // First unary operator?
            unary_it = it;
          last_unary_it = it;
        }
        else { // Term
          term_it = it;
        }
      }
      break;

    case TOKEN_NOT:
      if (paren_level == 0) {
        if (unary_it == to) /// First unary operator
          unary_it = it;
        last_unary_it = it;
      }
      break;

    default:
      break;
    }
  }

  if (paren_level != 0) { // paren_level should now be zero */
    throw std::runtime_error("mismatched parenthesis");
  }

  // This implement the operator hiearchy
  // Assignment
  if (assign_it != to)
    return parseAssign(eval, from, assign_it, to);

  // Ternary operator
  if (question_it != to || colon_it != to)
    return parseTernary(eval, from, question_it, colon_it, to);

  // Logical
  if (logical_it != to)
    return parseLogical(eval, from, logical_it, to);

  // Relational
  if (relation_it != to)
    return parseRelation(eval, from, relation_it, to);

  // Term
  if (term_it != to)
    return parseTerm(eval, from, term_it, to);

  // Factor
  if (factor_it != to)
    return parseFactor(eval, from, factor_it, to);

  // Unary
  if (unary_it != to)
    return parseUnary(eval, from, unary_it, to);

  if (expon_it != to)
    return parseFactor(eval, from, expon_it, to);



  // Parenthetical
  if (lparen_open_it != to) {
    if (lparen_open_it == from) {
      if (lparen_close_it == to - 1 && lparen_close_it - lparen_open_it > 1) {
        return parseExpression(eval, lparen_open_it + 1, lparen_close_it);
      }
      else
        throw std::runtime_error("syntax error parsing parentheses");
    }

    // Function
    if (lparen_open_it == from + 1) {
      if (lparen_close_it == to - 1)
        return parseFunction(eval, from, lparen_open_it, lparen_close_it, to);
      else // Closing paren not at to
        throw std::runtime_error("syntax error 2");
    }

    throw std::runtime_error("syntax error 3");
  }

  // Bracket
  if (lbrack_open_it != to) {

    // Array index
    if (lbrack_open_it == from + 1) {
      if (lbrack_close_it == to - 1)
        return parseIndex(eval, from, lbrack_open_it, lbrack_close_it, to);
      else // Closing brack not at to
        throw std::runtime_error("syntax error 2");
    }

    throw std::runtime_error("syntax error 3");
  }

  // R-Value
  return parseRValue(eval, from, to);
}


Node *
parseAssign(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator assign_it,
  LexemVector::const_iterator to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER)
    throw std::runtime_error("stk::expreval::parseAssign: expected identifier");

  Node *assign;

  if ((*(from + 1)).getToken() == TOKEN_ASSIGN || (*(from + 1)).getToken() == TOKEN_LBRACK) {
    assign = eval.newNode(OPCODE_ASSIGN);

    assign->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
    assign->m_data.variable.variable->setDependent();

    assign->m_right = parseExpression(eval, assign_it + 1, to);

    if ((*(from + 1)).getToken() == TOKEN_LBRACK) {
      assign->m_left = parseExpression(eval, from + 2, assign_it - 1);
    }
  }
  else
    throw std::runtime_error("syntax error");

  return assign;
}


Node *
parseTerm(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator term_it,
  LexemVector::const_iterator to)
{
  Node *term = eval.newNode((*term_it).getToken() == TOKEN_PLUS ? OPCODE_ADD : OPCODE_SUBTRACT);

  term->m_left = parseExpression(eval, from, term_it);
  term->m_right = parseExpression(eval, term_it + 1, to);

  return term;
}


Node *
parseFactor(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator factor_it,
  LexemVector::const_iterator to)
{
  Node *factor = eval.newNode
    (
     (*factor_it).getToken() == TOKEN_MULTIPLY ? OPCODE_MULTIPLY :
     (
      (*factor_it).getToken() == TOKEN_DIVIDE ? OPCODE_DIVIDE :
      (
       (*factor_it).getToken() == TOKEN_EXPONENTIATION ? OPCODE_EXPONENIATION : OPCODE_MODULUS
      )
     )
    );
  factor->m_left = parseExpression(eval, from, factor_it);
  factor->m_right = parseExpression(eval, factor_it + 1, to);

  return factor;
}


Node *
parseRelation(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator relation_it,
  LexemVector::const_iterator to)
{
  Opcode relation_opcode = OPCODE_UNDEFINED;

  switch ((*relation_it).getToken()) {
  case TOKEN_EQUAL:
    relation_opcode = OPCODE_EQUAL;
    break;

  case TOKEN_NOT_EQUAL:
    relation_opcode = OPCODE_NOT_EQUAL;
    break;

  case TOKEN_LESS:
    relation_opcode = OPCODE_LESS;
    break;

  case TOKEN_GREATER:
    relation_opcode = OPCODE_GREATER;
    break;

  case TOKEN_LESS_EQUAL:
    relation_opcode = OPCODE_LESS_EQUAL;
    break;

  case TOKEN_GREATER_EQUAL:
    relation_opcode = OPCODE_GREATER_EQUAL;
    break;

  default:
    break;
  }

  Node *relation = eval.newNode(relation_opcode);
  relation->m_left = parseExpression(eval, from, relation_it);
  relation->m_right = parseExpression(eval, relation_it + 1, to);

  return relation;
}


Node *
parseLogical(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator logical_it,
  LexemVector::const_iterator to)
{
  Node *logical = eval.newNode(((*logical_it).getToken() == TOKEN_LOGICAL_AND ? OPCODE_LOGICAL_AND : OPCODE_LOGICAL_OR));

  logical->m_left = parseExpression(eval, from, logical_it);
  logical->m_right = parseExpression(eval, logical_it + 1, to);

  return logical;
}


Node *
parseTernary(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator question_it,
  LexemVector::const_iterator colon_it,
  LexemVector::const_iterator to)
{
  if (question_it == to || colon_it == to)
    throw std::runtime_error("syntax error parsing ?: operator");

  Node *ternary = eval.newNode(OPCODE_TERNARY);

  ternary->m_left = parseExpression(eval, from, question_it);
  ternary->m_right = parseExpression(eval, question_it + 1, colon_it);
  ternary->m_other = parseExpression(eval, colon_it + 1, to);

  return ternary;
}


Node *
parseUnary(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator unary_it,
  LexemVector::const_iterator to)
{

  /* If it is a positive, just parse the internal of it */
  if ((*unary_it).getToken() == TOKEN_PLUS)
    return parseExpression(eval, unary_it + 1, to);
  else if ((*unary_it).getToken() == TOKEN_MINUS) {
    Node *unary = eval.newNode(OPCODE_UNARY_MINUS);
    unary->m_right = parseExpression(eval, unary_it + 1, to);
    return unary;
  }
  else if ((*unary_it).getToken() == TOKEN_NOT) {
    Node *unary = eval.newNode(OPCODE_UNARY_NOT);
    unary->m_right = parseExpression(eval, unary_it + 1, to);
    return unary;
  }
  else
    throw std::runtime_error("syntax error parsing unary operator");
}


Node *
parseFunction(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator lparen,
  LexemVector::const_iterator rparen,
  LexemVector::const_iterator to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER)
    throw std::runtime_error("syntax error parsing function");

  const std::string &function_name = (*from).getString();

  CFunctionBase *c_function = nullptr;
  CFunctionMap::iterator it = getCFunctionMap().find(function_name);
  Node *function = eval.newNode(OPCODE_FUNCTION);

  FunctionType functionType = eval.get_function_type(function_name);
  function->m_data.function.functionType = functionType;

  // only 1 function found with that function name.
  if ( getCFunctionMap().count(function_name) == 1 ) {
    if (it != getCFunctionMap().end()) {
      c_function = (*it).second;
    }
    function->m_data.function.function[0] = c_function;
    std::strncpy(function->m_data.function.functionName,
                 function_name.c_str(),
                 function_name.length() < Node::MAXIMUM_FUNCTION_NAME_LENGTH-1 ? function_name.length() : Node::MAXIMUM_FUNCTION_NAME_LENGTH-1);

    if (!c_function)
      eval.getUndefinedFunctionSet().insert(function_name);
  } 
  else {
    std::pair<CFunctionMap::iterator, CFunctionMap::iterator> ppp;
    ppp = getCFunctionMap().equal_range(function_name);
    int iCount=0;
    for (CFunctionMap::iterator it2 = ppp.first;
        it2 != ppp.second;
        ++it2,
        ++iCount)
    {
      c_function = (*it2).second;
      function->m_data.function.function[iCount] = c_function;
      std::strncpy(function->m_data.function.functionName,
                   function_name.c_str(),
                   function_name.length() < Node::MAXIMUM_FUNCTION_NAME_LENGTH-1 ? function_name.length() : Node::MAXIMUM_FUNCTION_NAME_LENGTH-1);


      if (!c_function)
        eval.getUndefinedFunctionSet().insert(function_name);
    }
    if( iCount == Node::MAXIMUM_NUMBER_OF_OVERLOADED_FUNCTION_NAMES) {
      throw std::runtime_error(std::string("Exceeded maximum number of overloaded function names for function named= ") + function_name);
    }
  }

  function->m_right = parseFunctionArg(eval, lparen + 1, rparen);

  if (!c_function) {
    function->m_data.function.undefinedFunction = true;
    eval.getUndefinedFunctionSet().insert(function_name);
  }

  return function;
}


Node *
parseIndex(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator lbrack,
  LexemVector::const_iterator rbrack,
  LexemVector::const_iterator to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER)
    throw std::runtime_error("syntax error parsing array");

  Node *index = eval.newNode(OPCODE_RVALUE);
  index->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
  index->m_left = parseExpression(eval, lbrack + 1, rbrack);

  return index;
}


Node *
parseFunctionArg(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator to)
{
  if (from == to) {
    return nullptr;
  }

  //
  //  Identify the first argument by looking for the first argument.  However, keep anything inside a bracket together and
  //  ignore commas inside those brackets
  //
  int paren_level=0;
  int brack_level=0;

  LexemVector::const_iterator endIterator = to;

  LexemVector::const_iterator it = from;
  while(it != to) {
    Token curToken = (*it).getToken();
    if(curToken == TOKEN_LPAREN) {
      paren_level++;
    } else if (curToken == TOKEN_RPAREN) {
      paren_level--;
      if (paren_level < 0) {
        throw std::runtime_error("mismatched parenthesis");
      }
    } else if (curToken ==TOKEN_LBRACK) {
      brack_level++;
    } else if(curToken == TOKEN_RBRACK) {
      brack_level--;
      if (brack_level < 0) {
        throw std::runtime_error("mismatched bracket");
      }
    } else if (curToken == TOKEN_COMMA) {
      if(paren_level == 0 && brack_level == 0) {
        endIterator = it;
        break;
      }
    } else {
    }
    ++it;
  }

  Node *argument = eval.newNode(OPCODE_ARGUMENT);
  argument->m_left = parseExpression(eval, from, endIterator);
  if (endIterator != to) {
    argument->m_right = parseFunctionArg(eval, endIterator + 1, to);
  }

  return argument;
}


Node *
parseRValue(
  Eval & eval,
  LexemVector::const_iterator from,
  LexemVector::const_iterator to)
{
  if (from + 1 != to)
    throw std::runtime_error(std::string("r-value not allowed following ") + (*from).getString());

  switch ((*from).getToken()) {
  case TOKEN_IDENTIFIER:
    {
      ConstantMap::iterator it = getConstantMap().find((*from).getString());
      if (it != getConstantMap().end()) {
        Node *constant = eval.newNode(OPCODE_CONSTANT);
        constant->m_data.constant.value = (*it).second;
        return constant;
      }
      // Define a variable
      else {
        Node *variable = eval.newNode(OPCODE_RVALUE);
        variable->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
        return variable;
      }
    }

  case TOKEN_REAL_CONSTANT:
  case TOKEN_INTEGER_CONSTANT:
    {
      Node *constant = eval.newNode(OPCODE_CONSTANT);
      constant->m_data.constant.value = (*from).getValue<double>();
      return constant;
    }

  default:
    throw std::runtime_error("invalid rvalue");
  }
}

} // namespace Parser

Eval::Eval(VariableMap::Resolver & resolver, const std::string & expression, Variable::ArrayOffset arrayOffsetType)
  : m_variableMap(resolver),
    m_expression(expression),
    m_syntaxStatus(false),
    m_parseStatus(false),
    m_headNode(nullptr),
    m_arrayOffsetType(arrayOffsetType),
    m_parsedEval(nullptr)
{
  initialize_function_map();
}

Eval::Eval(const std::string & expression, Variable::ArrayOffset arrayOffsetType)
  : m_variableMap(VariableMap::getDefaultResolver()),
    m_expression(expression),
    m_syntaxStatus(false),
    m_parseStatus(false),
    m_headNode(nullptr),
    m_arrayOffsetType(arrayOffsetType),
    m_parsedEval(nullptr)
{
  initialize_function_map();
}

Eval::Eval(const Eval& otherEval) 
  : m_variableMap(otherEval.m_variableMap),
    m_undefinedFunctionSet(otherEval.m_undefinedFunctionSet),
    m_functionMap(otherEval.m_functionMap),
    m_expression(otherEval.m_expression),
    m_syntaxStatus(otherEval.m_syntaxStatus),
    m_parseStatus(otherEval.m_parseStatus),
    m_headNode(otherEval.m_headNode),
    m_nodes(otherEval.m_nodes),
    m_arrayOffsetType(otherEval.m_arrayOffsetType),
    m_parsedEval(nullptr)
{} 

Eval::~Eval() 
{
  delete m_parsedEval;
}

Node*
Eval::newNode(int opcode)
{
  m_nodes.push_back(std::make_shared<Node>(static_cast<Opcode>(opcode), this));
  m_nodes.back().get()->m_currentNodeIndex = m_nodes.size() - 1;
  return m_nodes.back().get();
}

void
Eval::syntax()
{
  m_syntaxStatus = false;
  m_parseStatus = false;

  try {
    // Validate the characters
    LexemVector lex_vector = tokenize(m_expression);

    // Call the multiparse routine to parse subexpressions
    m_headNode = Parser::parseStatements(*this, lex_vector.begin(), lex_vector.end());

    m_syntaxStatus = true;
  }
  catch (std::runtime_error &) {
  }
}

void
Eval::parse()
{
  try {
    syntax();

    if (m_syntaxStatus) {
      if (!m_undefinedFunctionSet.empty()) {
        std::ostringstream strout;
        strout << "In expression '" << m_expression << "', the following functions are not defined:" << std::endl;
        for (const auto& it : m_undefinedFunctionSet)
          strout << it << std::endl;
        throw std::runtime_error(strout.str());
      }

      resolve();

      m_parseStatus = true;
    } else {
      throw std::runtime_error("");
    }
  }
  catch (std::runtime_error & ) {
    throw;
  }

  int index = 0;
  for (VariableMap::iterator it = m_variableMap.begin(); it != m_variableMap.end(); ++it) {
    it->second->set_index(index);
    ++index;
  }
  
}

void
Eval::resolve()
{
  for (VariableMap::iterator it = m_variableMap.begin(); it != m_variableMap.end(); ++it) {
    m_variableMap.getResolver().resolve(it);
  }
}

double
Eval::evaluate() const
{
  /* Make sure it was parsed successfully */
  if (!m_parseStatus) {
    throw std::runtime_error(std::string("Expression '") + m_expression + "' did not parse successfully");
  }

  double returnValue = 0.0;
  try
  {
    if(m_headNode) {
      returnValue = m_headNode->eval();
    }
  }
  catch(expression_evaluation_exception &)
  {
    throw std::runtime_error(std::string("Expression '") + m_expression + "' did not evaluate successfully");
  }
  return returnValue;
}

bool
Eval::undefinedFunction() const
{
  /* Check for an undefined function in any allocated node */
  for (const auto& node : m_nodes) {
    if (node->m_data.function.undefinedFunction) return true;
  }
  return false;
}

bool 
Eval::is_constant_expression() const
{
  return m_variableMap.empty();
}

bool 
Eval::is_variable(const std::string& variableName) const
{
  return (m_variableMap.count(variableName) > 0);
}

bool
Eval::is_scalar(const std::string& variableName) const
{
  auto variableIterator = m_variableMap.find(variableName);

  if (variableIterator == m_variableMap.end()) { 
    return false; 
  }
  
  int variableLength = variableIterator->second->getLength();
  return variableLength == 1 || variableLength == std::numeric_limits<int>::max();  
}

std::vector<std::string> 
Eval::get_variable_names() const
{
  std::vector<std::string> variableList;
  for(auto& currentVariable : m_variableMap) {
    std::string variableName = currentVariable.first;
    variableList.push_back(variableName);
  }
  
  return variableList;
}

std::vector<std::string> 
Eval::get_dependent_variable_names() const
{
  std::vector<std::string> dependentVariableList;
  for(auto& currentVariable : m_variableMap) {
    std::string variableName = currentVariable.first;
    stk::expreval::Variable* variable = currentVariable.second.get();
    if (variable->isDependent()) {
      dependentVariableList.push_back(variableName);
    }
  }
  
  return dependentVariableList;
}

std::vector<std::string> 
Eval::get_independent_variable_names() const
{
  std::vector<std::string> independentVariableList;
  for(auto& currentVariable : m_variableMap) {
    std::string variableName = currentVariable.first;
    stk::expreval::Variable* variable = currentVariable.second.get();
    if (!(variable->isDependent())) {
      independentVariableList.push_back(variableName);
    }
  }
  
  return independentVariableList;
}

int
Eval::get_variable_index(const std::string& variable) const
{
  const auto variableIter = m_variableMap.find(variable);
  ThrowRequireMsg (variableIter != m_variableMap.end(), "Variable " + variable + " Not Found in VariableMap"); 
  return variableIter->second->get_index();
}

int 
Eval::get_head_node_index() const 
{ 
  return (m_headNode) ? m_headNode->m_currentNodeIndex : -1; 

}

FunctionType 
Eval::get_function_type(const std::string& functionName) const
{
  const auto functionIt = m_functionMap.find(functionName);
  if (functionIt != m_functionMap.end()) {
    return functionIt->second;
  }
  else {
    return FunctionType::UNDEFINED;
  }
}

void
Eval::initialize_function_map()
{
  m_functionMap["abs"] = FunctionType::ABS;
  m_functionMap["fabs"] = FunctionType::ABS;
  m_functionMap["max"] = FunctionType::MAX;
  m_functionMap["min"] = FunctionType::MIN;
  m_functionMap["sign"] = FunctionType::SIGN;
  m_functionMap["ipart"] = FunctionType::IPART;
  m_functionMap["fpart"] = FunctionType::FPART;
  m_functionMap["ceil"] = FunctionType::CEIL;
  m_functionMap["floor"] = FunctionType::FLOOR;
  m_functionMap["mod"] = FunctionType::MOD;
  m_functionMap["fmod"] = FunctionType::MOD;
  m_functionMap["pow"] = FunctionType::POW;
  m_functionMap["sqrt"] = FunctionType::SQRT;
  m_functionMap["exp"] = FunctionType::EXP;
  m_functionMap["ln"] = FunctionType::LN;
  m_functionMap["log"] = FunctionType::LN;
  m_functionMap["log10"] = FunctionType::LOG10;


  m_functionMap["deg"] = FunctionType::DEG;
  m_functionMap["rad"] = FunctionType::RAD;
  m_functionMap["sin"] = FunctionType::SIN;
  m_functionMap["cos"] = FunctionType::COS;
  m_functionMap["tan"] = FunctionType::TAN;
  m_functionMap["asin"] = FunctionType::ASIN;
  m_functionMap["acos"] = FunctionType::ACOS;
  m_functionMap["atan"] = FunctionType::ATAN;
  m_functionMap["atan2"] = FunctionType::ATAN2;
  m_functionMap["sinh"] = FunctionType::SINH;
  m_functionMap["cosh"] = FunctionType::COSH;
  m_functionMap["tanh"] = FunctionType::TANH;
  m_functionMap["asinh"] = FunctionType::ASINH;
  m_functionMap["acosh"] = FunctionType::ACOSH;
  m_functionMap["atanh"] = FunctionType::ATANH;
  m_functionMap["erf"] = FunctionType::ERF;
  m_functionMap["erfc"] = FunctionType::ERFC;
  m_functionMap["poltorectx"] = FunctionType::POLTORECTX;
  m_functionMap["poltorecty"] = FunctionType::POLTORECTY;
  m_functionMap["recttopolr"] = FunctionType::RECTTOPOLR;
  m_functionMap["recttopola"] = FunctionType::RECTTOPOLA;

  m_functionMap["unit_step"] = FunctionType::UNIT_STEP;
  m_functionMap["cycloidal_ramp"] = FunctionType::CYCLOIDAL_RAMP;
  m_functionMap["cos_ramp"] = FunctionType::COS_RAMP;
  m_functionMap["cosine_ramp"] = FunctionType::COS_RAMP;
  m_functionMap["haversine_pulse"] = FunctionType::HAVERSINE_PULSE;
  m_functionMap["point2d"] = FunctionType::POINT2D;
  m_functionMap["point3d"] = FunctionType::POINT3D;

  m_functionMap["exponential_pdf"] = FunctionType::EXPONENTIAL_PDF;
  m_functionMap["log_uniform_pdf"] = FunctionType::LOG_UNIFORM_PDF;
  m_functionMap["normal_pdf"] = FunctionType::NORMAL_PDF;
  m_functionMap["weibull_pdf"] = FunctionType::WEIBULL_PDF;
  m_functionMap["gamma_pdf"] = FunctionType::GAMMA_PDF;

  m_functionMap["rand"] = FunctionType::RAND;
  m_functionMap["srand"] = FunctionType::SRAND;
  m_functionMap["random"] = FunctionType::RANDOM;
  m_functionMap["ts_random"] = FunctionType::TS_RANDOM;
  m_functionMap["ts_normal"] = FunctionType::TS_NORMAL;
  m_functionMap["time"] = FunctionType::TIME;

}

ParsedEval& 
Eval::get_parsed_eval() {
  if (m_parsedEval == nullptr) {
    m_parsedEval = new ParsedEval(*this);
  }
  return *m_parsedEval;
}

ParsedEval::ParsedEval(Eval& eval)
  : m_numVariables(eval.get_variable_count()),
    m_arrayOffsetType(eval.getArrayOffsetType()),
    m_headNodeIndex(eval.get_head_node_index()),
    m_deviceNodes(Kokkos::view_alloc(Kokkos::WithoutInitializing, "deviceNodeTree"), eval.get_node_count()),
    m_hostNodes(Kokkos::view_alloc(Kokkos::WithoutInitializing, "hostNodeTree"), eval.get_node_count())
{
  ThrowRequireMsg (eval.getParseStatus(), std::string("Expression '") + eval.getExpression() + "' did not parse successfully");  

  for (int i = 0; i < eval.get_node_count(); ++i) {
      m_hostNodes(i) = NgpNode(*(eval.get_node(i)));
  } 
  Kokkos::deep_copy(m_deviceNodes, m_hostNodes);
}



} // namespace expreval
} // namespace stk
