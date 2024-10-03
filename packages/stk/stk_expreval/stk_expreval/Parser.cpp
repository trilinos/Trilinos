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

#include "stk_expreval/Parser.hpp"
#include "stk_expreval/Eval.hpp"
#include "stk_expreval/Node.hpp"

namespace stk {
namespace expreval {
namespace Parser {

Node *parseStatements(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseStatement(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseExpression(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseAssign(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator assign, LexemVector::const_iterator to);
Node *parseTerm(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator term, LexemVector::const_iterator to);
Node *parseFactor(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator factor, LexemVector::const_iterator to);
Node *parseRelation(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator factor, LexemVector::const_iterator to);
Node *parseLogical(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator factor, LexemVector::const_iterator to);
Node *parseUnary(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator unary, LexemVector::const_iterator to);
Node *parseTernary(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator question, LexemVector::const_iterator colon, LexemVector::const_iterator to);
Node *parseFunction(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator lparen, LexemVector::const_iterator rparen, LexemVector::const_iterator to);
Node *parseFunctionArg(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseRValue(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseIndex(Eval & eval, LexemVector::const_iterator from, LexemVector::const_iterator lbrack, LexemVector::const_iterator rbrack, LexemVector::const_iterator to);
int countFunctionArgs(Node* function);

Node *
parseStatements(Eval &eval,
                LexemVector::const_iterator from,
                LexemVector::const_iterator to)
{
  if ((*from).getToken() == TOKEN_END) {
    return nullptr;
  }

  if ((*from).getToken() == TOKEN_SEMI) {
    return parseStatements(eval, from + 1, to);
  }

  // Technically, there should be no check for TOKEN_END, but we allow a missing final TOKEN_SEMI
  LexemVector::const_iterator it;
  for (it = from; (*it).getToken() != TOKEN_SEMI && (*it).getToken() != TOKEN_END; ++it) {}

  Node *statement = eval.newNode(OPCODE_STATEMENT);
  statement->m_left = parseStatement(eval, from, it);

  // Technically, there should be no check for TOKEN_END, but we allow a missing final TOKEN_SEMI
  if ((*it).getToken() != TOKEN_END) {
    statement->m_right = parseStatements(eval, it + 1, to);
  }

  return statement;
}


Node *
parseStatement(Eval & eval,
               LexemVector::const_iterator from,
               LexemVector::const_iterator to)
{
  return parseExpression(eval, from, to);
}


Node *
parseExpression(Eval & eval,
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
    case TOKEN_LPAREN: {
      if (paren_level == 0 && lparen_open_it == to && brack_level == 0 && lbrack_open_it == to) {
        lparen_open_it = it;
      }
      paren_level++;
      break;
    }

    case TOKEN_RPAREN: {
      paren_level--;

      if (paren_level == 0 && lparen_close_it == to && brack_level == 0 && lbrack_close_it == to) {
        lparen_close_it = it;
      }

      if (paren_level < 0) {
        throw std::runtime_error("mismatched parenthesis");
      }
      break;
    }

    case TOKEN_LBRACK: {
      if (paren_level == 0 && lparen_open_it == to && brack_level == 0 && lbrack_open_it == to) {
        lbrack_open_it = it;
      }
      brack_level++;
      break;
    }

    case TOKEN_RBRACK: {
      brack_level--;

      if (paren_level == 0 && lparen_close_it == to && brack_level == 0 && lbrack_close_it == to) {
        lbrack_close_it = it;
      }

      if (brack_level < 0) {
        throw std::runtime_error("mismatched bracket");
      }
      break;
    }

    case TOKEN_ASSIGN: {
      if (paren_level == 0 && assign_it == to) {
        assign_it = it;
      }
      break;
    }

    case TOKEN_QUESTION: {
      if (paren_level == 0 && question_it == to) {
        question_it = it;
      }
      break;
    }

    case TOKEN_COLON: {
      if (paren_level == 0) {
        colon_it = it;
      }
      break;
    }

    case TOKEN_EXPONENTIATION: {
      if (paren_level == 0) {
        expon_it = it;
      }
      break;
    }

    case TOKEN_MULTIPLY:
    case TOKEN_DIVIDE:
    case TOKEN_PERCENT: {
      if (paren_level == 0) {
        factor_it = it;
      }
      break;
    }

    case TOKEN_EQUAL:
    case TOKEN_NOT_EQUAL:
    case TOKEN_LESS:
    case TOKEN_GREATER:
    case TOKEN_LESS_EQUAL:
    case TOKEN_GREATER_EQUAL: {
      if (paren_level == 0 && relation_it == to) {
        relation_it = it;
      }
      break;
    }

    case TOKEN_LOGICAL_AND:
    case TOKEN_LOGICAL_OR: {
      if (paren_level == 0 && logical_it == to) {
        logical_it = it;
      }
      break;
    }

    case TOKEN_PLUS:
    case TOKEN_MINUS: {
      if (paren_level == 0) {
        // After any of these, we are a unary operator, not a term
        if (it == from || it == assign_it + 1 || it == term_it + 1 || it == factor_it + 1 ||
            it == last_unary_it + 1 || it == expon_it + 1)
        { // Unary operator
          if (unary_it == to) { // First unary operator?
            unary_it = it;
          }
          last_unary_it = it;
        }
        else { // Term
          term_it = it;
        }
      }
      break;
    }

    case TOKEN_NOT: {
      if (paren_level == 0) {
        if (unary_it == to) { /// First unary operator
          unary_it = it;
        }
        last_unary_it = it;
      }
      break;
    }

    default: {
      break;
    }
    }
  }

  if (paren_level != 0) { // paren_level should now be zero
    throw std::runtime_error("mismatched parenthesis");
  }

  // This implement the operator hiearchy
  // Assignment
  if (assign_it != to) {
    return parseAssign(eval, from, assign_it, to);
  }

  // Ternary operator
  if (question_it != to || colon_it != to) {
    return parseTernary(eval, from, question_it, colon_it, to);
  }

  // Logical
  if (logical_it != to) {
    return parseLogical(eval, from, logical_it, to);
  }

  // Relational
  if (relation_it != to) {
    return parseRelation(eval, from, relation_it, to);
  }

  // Term
  if (term_it != to) {
    return parseTerm(eval, from, term_it, to);
  }

  // Factor
  if (factor_it != to) {
    return parseFactor(eval, from, factor_it, to);
  }

  // Unary
  if (unary_it != to) {
    return parseUnary(eval, from, unary_it, to);
  }

  if (expon_it != to) {
    return parseFactor(eval, from, expon_it, to);
  }


  // Parenthetical
  if (lparen_open_it != to) {
    if (lparen_open_it == from) {
      if (lparen_close_it == to - 1 && lparen_close_it - lparen_open_it > 1) {
        return parseExpression(eval, lparen_open_it + 1, lparen_close_it);
      }
      else {
        throw std::runtime_error("syntax error parsing parentheses");
      }
    }

    // Function
    if (lparen_open_it == from + 1) {
      if (lparen_close_it == to - 1) {
        return parseFunction(eval, from, lparen_open_it, lparen_close_it, to);
      }
      else { // Closing paren not at to
        throw std::runtime_error("syntax error 2");
      }
    }

    throw std::runtime_error("syntax error 3");
  }

  // Bracket
  if (lbrack_open_it != to) {

    // Array index
    if (lbrack_open_it == from + 1) {
      if (lbrack_close_it == to - 1) {
        return parseIndex(eval, from, lbrack_open_it, lbrack_close_it, to);
      }
      else { // Closing brack not at to
        throw std::runtime_error("syntax error 2");
      }
    }

    throw std::runtime_error("syntax error 3");
  }

  // R-Value
  return parseRValue(eval, from, to);
}

Node *
parseAssign(Eval & eval,
            LexemVector::const_iterator from,
            LexemVector::const_iterator assign_it,
            LexemVector::const_iterator to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER) {
    throw std::runtime_error("stk::expreval::parseAssign: expected identifier");
  }

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
  else {
    throw std::runtime_error("syntax error");
  }

  return assign;
}

Node *
parseTerm(Eval & eval,
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
parseFactor(Eval & eval,
            LexemVector::const_iterator from,
            LexemVector::const_iterator factor_it,
            LexemVector::const_iterator to)
{
  Node *factor = eval.newNode((*factor_it).getToken() == TOKEN_MULTIPLY ? OPCODE_MULTIPLY :
                              ((*factor_it).getToken() == TOKEN_DIVIDE ? OPCODE_DIVIDE :
                               ((*factor_it).getToken() == TOKEN_EXPONENTIATION ? OPCODE_EXPONENIATION : OPCODE_MODULUS)));

  factor->m_left = parseExpression(eval, from, factor_it);
  factor->m_right = parseExpression(eval, factor_it + 1, to);

  return factor;
}

bool isRelation(Node* node)
{
  switch (node->m_opcode) {
    case OPCODE_EQUAL:
    case OPCODE_NOT_EQUAL:
    case OPCODE_LESS:
    case OPCODE_GREATER:
    case OPCODE_LESS_EQUAL:
    case OPCODE_GREATER_EQUAL:
      return true;
    default:
      return false;
  }
}

Node *
parseRelation(Eval & eval,
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

  if (isRelation(relation->m_left) || isRelation(relation->m_right)) {
    throw std::runtime_error("stk::expreval::parseRelation: stk_expreval does not support chained comparisons");
  }

  return relation;
}

Node *
parseLogical(Eval & eval,
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
parseTernary(Eval & eval,
             LexemVector::const_iterator from,
             LexemVector::const_iterator question_it,
             LexemVector::const_iterator colon_it,
             LexemVector::const_iterator to)
{
  if (question_it == to || colon_it == to) {
    throw std::runtime_error("syntax error parsing ?: operator");
  }

  Node *predicate = eval.newNode(OPCODE_TERNARY_PREDICATE);
  predicate->m_left = parseExpression(eval, from, question_it);

  Node *join = eval.newNode(OPCODE_TERNARY_JOIN);
  join->m_left = parseExpression(eval, question_it + 1, colon_it);
  join->m_right = parseExpression(eval, colon_it + 1, to);
  join->m_ternaryOther = predicate;

  predicate->m_ternaryOther = join;

  return predicate;
}


Node *
parseUnary(Eval & eval,
           LexemVector::const_iterator from,
           LexemVector::const_iterator unary_it,
           LexemVector::const_iterator to)
{
  /* If it is a positive, just parse the internal of it */
  if ((*unary_it).getToken() == TOKEN_PLUS) {
    return parseExpression(eval, unary_it + 1, to);
  }
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
  else {
    throw std::runtime_error("syntax error parsing unary operator");
  }
}


Node *
parseFunction(Eval & eval,
              LexemVector::const_iterator from,
              LexemVector::const_iterator lparen,
              LexemVector::const_iterator rparen,
              LexemVector::const_iterator to)
{
  using CFunctionIterPair = std::pair<CFunctionMap::iterator, CFunctionMap::iterator>;

  if ((*from).getToken() != TOKEN_IDENTIFIER) {
    throw std::runtime_error("syntax error parsing function");
  }

  const std::string &function_name = (*from).getString();

  Node *function = eval.newNode(OPCODE_FUNCTION);

  function->m_right = parseFunctionArg(eval, lparen + 1, rparen);
  int argc = countFunctionArgs(function);
  if (argc > Node::MAXIMUM_NUMBER_OF_FUNCTION_ARGS)
  {
    throw std::runtime_error(std::string("too many arguments to function, limit is ") +
                                         std::to_string(Node::MAXIMUM_NUMBER_OF_FUNCTION_ARGS));
  }

  FunctionType functionType = eval.get_function_type(function_name);

  CFunctionIterPair iteratorPair = getCFunctionMap().equal_range(function_name);
  function->m_data.function.function = nullptr;
  for (auto it = iteratorPair.first; it != iteratorPair.second; ++it)
  {
    CFunctionBase* c_function = it->second;
    if (c_function->getArgCount() == argc)
    {
      function->m_data.function.function = c_function;
      function->m_data.function.functionType = functionType;
      std::strncpy(function->m_data.function.functionName,
                   function_name.c_str(),
                   function_name.length() < Node::MAXIMUM_FUNCTION_NAME_LENGTH-1 ? function_name.length() : Node::MAXIMUM_FUNCTION_NAME_LENGTH-1);
      break;
    }
  }

  if (!(function->m_data.function.function)) {
    eval.getUndefinedFunctionSet().insert(function_name);
  }

  return function;
}


Node *
parseIndex(Eval & eval,
           LexemVector::const_iterator from,
           LexemVector::const_iterator lbrack,
           LexemVector::const_iterator rbrack,
           LexemVector::const_iterator to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER) {
    throw std::runtime_error("syntax error parsing array");
  }

  Node *index = eval.newNode(OPCODE_RVALUE);
  index->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
  index->m_left = parseExpression(eval, lbrack + 1, rbrack);

  return index;
}

Node *
parseFunctionArg(Eval & eval,
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
  while (it != to) {
    Token curToken = (*it).getToken();
    if (curToken == TOKEN_LPAREN) {
      paren_level++;
    }
    else if (curToken == TOKEN_RPAREN) {
      paren_level--;
      if (paren_level < 0) {
        throw std::runtime_error("mismatched parenthesis");
      }
    }
    else if (curToken ==TOKEN_LBRACK) {
      brack_level++;
    }
    else if(curToken == TOKEN_RBRACK) {
      brack_level--;
      if (brack_level < 0) {
        throw std::runtime_error("mismatched bracket");
      }
    }
    else if (curToken == TOKEN_COMMA) {
      if (paren_level == 0 && brack_level == 0) {
        endIterator = it;
        break;
      }
    }
    else {
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
parseRValue(Eval & eval,
            LexemVector::const_iterator from,
            LexemVector::const_iterator to)
{
  if (from + 1 != to) {
    throw std::runtime_error(std::string("r-value not allowed following ") + (*from).getString());
  }

  switch ((*from).getToken()) {
  case TOKEN_IDENTIFIER: {
    ConstantMap::iterator it = getConstantMap().find((*from).getString());
    if (it != getConstantMap().end()) {
      Node *constant = eval.newNode(OPCODE_CONSTANT);
      constant->m_data.constant.value = (*it).second;
      return constant;
    }
    else {
      // Define a variable
      Node *variable = eval.newNode(OPCODE_RVALUE);
      variable->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
      return variable;
    }
  }

  case TOKEN_REAL_CONSTANT:
  case TOKEN_INTEGER_CONSTANT: {
    Node *constant = eval.newNode(OPCODE_CONSTANT);
    constant->m_data.constant.value = (*from).getValue<double>();
    return constant;
  }

  default: {
    throw std::runtime_error("invalid rvalue");
  }
  }
}

int countFunctionArgs(Node* function)
{
  if (function->m_opcode != OPCODE_FUNCTION)
  {
    throw std::runtime_error("Node must be a function to count its arguments");
  }

  int argc = 0;
  for (Node *arg = function->m_right; arg; arg = arg->m_right) {
    argc++;
  }

  return argc;
}

}
}
}
