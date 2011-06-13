/**   ------------------------------------------------------------
 *    Copyright 2003 - 2011 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
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
#include <stdexcept>
#include <cstdlib>
#include <cctype>
#include <cmath>

#include <math.h>
#include <time.h>

#include <stk_expreval/Evaluator.hpp>
#include <stk_expreval/Lexer.hpp>

namespace stk {
namespace expreval {

/**
 * @brief Enumeration <b>Opcode</b> lists the operation codes which can be
 * executed using the execution virtual machine.
 *
 */
enum Opcode {
  OPCODE_UNDEFINED,
  OPCODE_CONSTANT,
  OPCODE_RVALUE,
  OPCODE_STATEMENT,
  OPCODE_ARGUMENT,

  OPCODE_TIERNARY,
  OPCODE_MULTIPLY,
  OPCODE_DIVIDE,
  OPCODE_MODULUS,
  OPCODE_ADD,
  OPCODE_SUBTRACT,
  OPCODE_UNARY_MINUS,
  OPCODE_FUNCTION,

  OPCODE_EQUAL,
  OPCODE_NOT_EQUAL,
  OPCODE_LESS,
  OPCODE_GREATER,
  OPCODE_LESS_EQUAL,
  OPCODE_GREATER_EQUAL,

  OPCODE_UNARY_NOT,
  OPCODE_LOGICAL_AND,
  OPCODE_LOGICAL_OR,

  OPCODE_ASSIGN
};


class Node
{
public:
  explicit Node(Opcode opcode)
    : m_opcode(opcode),
      m_left(0),
      m_right(0),
      m_other(0)
  {}

private:
  explicit Node(const Node &);
  Node &operator=(const Node &);

public:
  ~Node() {
    delete m_left;
    delete m_right;
    delete m_other;
  }

  double eval() const;

  const Opcode	m_opcode;

  union _data
  {
    struct _constant
    {
      double value;
    } constant;

    struct _variable
    {
      Variable *variable;
    } variable;

    struct _function
    {
      CFunctionBase *	function;
    } function;
  } m_data;

  Node *		m_left;
  Node *		m_right;
  Node *		m_other;
};


double
Node::eval() const
{
  /* %TRACE[NONE]% */  /* %TRACE% */
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
    if (m_left)
      return (*m_data.variable.variable)[m_left->eval()];
    else
      return m_data.variable.variable->getValue();

  case OPCODE_MULTIPLY:
    return m_left->eval()*m_right->eval();

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
  case OPCODE_TIERNARY:
    return m_left->eval() != s_false ? m_right->eval() : m_other->eval();

  case OPCODE_UNARY_MINUS:
    return -m_right->eval();

  case OPCODE_UNARY_NOT:
    return m_right->eval() == s_false ? s_true : s_false;

  case OPCODE_ASSIGN:
    if (m_left)
      return (*m_data.variable.variable)[m_left->eval()] = m_right->eval();
    else {
      *m_data.variable.variable = m_right->eval();
      return m_data.variable.variable->getValue();
    }

  case OPCODE_FUNCTION:
    {
      double argv[20];

      int argc = 0;
      for (Node *arg = m_right; arg; arg = arg->m_right)
	argv[argc++] = arg->m_left->eval();

      return (*m_data.function.function)(argc, argv);
    }

  default: // Unknown opcode
    throw std::runtime_error("Evaluation error");
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
Node *parseTiernary(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator question, LexemVector::const_iterator colon, LexemVector::const_iterator to);
Node *parseFunction(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator lparen, LexemVector::const_iterator rparen, LexemVector::const_iterator to);
Node *parseFunctionArg(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseRValue(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseIndex(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator lbrack, LexemVector::const_iterator rbrack, LexemVector::const_iterator to);

// Parser productions

Node *
parseStatements(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  if ((*from).getToken() == TOKEN_END)
    return NULL;

  if ((*from).getToken() == TOKEN_SEMI)
    return parseStatements(eval, from + 1, to);

  // Technically, there should be no check for TOKEN_END, but we allow a missing final TOKEN_SEMI
  LexemVector::const_iterator it;
  for (it = from; (*it).getToken() != TOKEN_SEMI && (*it).getToken() != TOKEN_END; ++it)
    ;

  Node *statement = new Node(OPCODE_STATEMENT);
  statement->m_left = parseStatement(eval, from, it);

  // Technically, there should be no check for TOKEN_END, but we allow a missing final TOKEN_SEMI
  if ((*it).getToken() != TOKEN_END)
    statement->m_right = parseStatements(eval, it + 1, to);

  return statement;
}


Node *
parseStatement(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  return parseExpression(eval, from, to);
}


Node *
parseExpression(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  int paren_level = 0;					// Paren level
  int brack_level = 0;					// Brack level
  LexemVector::const_iterator lparen_open_it = to;	// First open paren
  LexemVector::const_iterator lparen_close_it = to;	// Corresponding close paren
  LexemVector::const_iterator lbrack_open_it = to;	// First open bracket
  LexemVector::const_iterator lbrack_close_it = to;	// Corresponding close brack
  LexemVector::const_iterator assign_it = to;		// First = at paren_level 0 for assignment
  LexemVector::const_iterator term_it = to;		// Last + or - at paren_level 0 for adding or subtracting
  LexemVector::const_iterator factor_it = to;		// Last * or / at paren_level 0 for multiplying or dividing
  LexemVector::const_iterator relation_it = to;		// Last relational at paren_level 0 for relational operator
  LexemVector::const_iterator logical_it = to;		// Last logical at paren_level 0 for logical operator
  LexemVector::const_iterator question_it = to;		// Last tiernary at paren_level 0 for tiernary operator
  LexemVector::const_iterator colon_it = to;
  LexemVector::const_iterator unary_it = to;		// First +,- at plevel 0 for positive,negative
  LexemVector::const_iterator last_unary_it = to;	// Last +,- found at plevel for for positive,negative

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

      if (paren_level < 0)
	throw std::runtime_error("mismatched parenthesis");
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
	    || it == last_unary_it + 1)
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

  if (paren_level != 0) // paren_level should now be zero */
    throw std::runtime_error("mismatched parenthesis");

  // This implement the operator hiearchy
  // Assignment
  if (assign_it != to)
    return parseAssign(eval, from, assign_it, to);

  // Tiernary operator
  if (question_it != to || colon_it != to)
    return parseTiernary(eval, from, question_it, colon_it, to);

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

  // Parenthetical
  if (lparen_open_it != to) {
    if (lparen_open_it == from) {
      if (lparen_close_it == to - 1 && lparen_close_it - lparen_open_it > 1)
	return parseExpression(eval, lparen_open_it + 1, lparen_close_it);
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
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	assign_it,
  LexemVector::const_iterator	to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER) //  || from + 1 != assign_it) {
    throw std::runtime_error("stk::expreval::parseAssign: expected identifier");

  Node *assign;

  if ((*(from + 1)).getToken() == TOKEN_ASSIGN || (*(from + 1)).getToken() == TOKEN_LBRACK) {
    assign = new Node(OPCODE_ASSIGN);

    assign->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
    assign->m_data.variable.variable->setDependent();
    assign->m_right = parseExpression(eval, assign_it + 1, to);

    if ((*(from + 1)).getToken() == TOKEN_LBRACK)
      assign->m_left = parseExpression(eval, from + 2, assign_it - 1);
  }
  else
    throw std::runtime_error("syntax error");

  return assign;
}


Node *
parseTerm(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	term_it,
  LexemVector::const_iterator	to)
{
  Node *term = new Node((*term_it).getToken() == TOKEN_PLUS ? OPCODE_ADD : OPCODE_SUBTRACT);

  term->m_left = parseExpression(eval, from, term_it);
  term->m_right = parseExpression(eval, term_it + 1, to);

  return term;
}


Node *
parseFactor(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	factor_it,
  LexemVector::const_iterator	to)
{
  Node *factor = new Node((*factor_it).getToken() == TOKEN_MULTIPLY ? OPCODE_MULTIPLY : ((*factor_it).getToken() == TOKEN_DIVIDE ? OPCODE_DIVIDE : OPCODE_MODULUS));

  factor->m_left = parseExpression(eval, from, factor_it);
  factor->m_right = parseExpression(eval, factor_it + 1, to);

  return factor;
}


Node *
parseRelation(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	relation_it,
  LexemVector::const_iterator	to)
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

  Node *relation = new Node(relation_opcode);

  relation->m_left = parseExpression(eval, from, relation_it);
  relation->m_right = parseExpression(eval, relation_it + 1, to);

  return relation;
}


Node *
parseLogical(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	logical_it,
  LexemVector::const_iterator	to)
{
  Node *logical = new Node(((*logical_it).getToken() == TOKEN_LOGICAL_AND ? OPCODE_LOGICAL_AND : OPCODE_LOGICAL_OR));

  logical->m_left = parseExpression(eval, from, logical_it);
  logical->m_right = parseExpression(eval, logical_it + 1, to);

  return logical;
}


Node *
parseTiernary(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	question_it,
  LexemVector::const_iterator	colon_it,
  LexemVector::const_iterator	to)
{
  if (question_it == to || colon_it == to)
    throw std::runtime_error("syntax error parsing ?: operator");

  Node *tiernary = new Node(OPCODE_TIERNARY);

  tiernary->m_left = parseExpression(eval, from, question_it);
  tiernary->m_right = parseExpression(eval, question_it + 1, colon_it);
  tiernary->m_other = parseExpression(eval, colon_it + 1, to);

  return tiernary;
}


Node *
parseUnary(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	unary_it,
  LexemVector::const_iterator	to)
{
  /* If it is a positive, just parse the internal of it */
  if ((*unary_it).getToken() == TOKEN_PLUS)
    return parseExpression(eval, unary_it + 1, to);
  else if ((*unary_it).getToken() == TOKEN_MINUS) {
    Node *unary = new Node(OPCODE_UNARY_MINUS);
    unary->m_right = parseExpression(eval, unary_it + 1, to);
    return unary;
  }
  else if ((*unary_it).getToken() == TOKEN_NOT) {
    Node *unary = new Node(OPCODE_UNARY_NOT);
    unary->m_right = parseExpression(eval, unary_it + 1, to);
    return unary;
  }
  else
    throw std::runtime_error("syntax error parsing unary operator");
}


Node *
parseFunction(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	lparen,
  LexemVector::const_iterator	rparen,
  LexemVector::const_iterator	to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER)
    throw std::runtime_error("syntax error parsing function");

  const std::string &function_name = (*from).getString();

  CFunctionBase *c_function = NULL;
  CFunctionMap::iterator it = getCFunctionMap().find(function_name);
  if (it != getCFunctionMap().end())
    c_function = (*it).second;

//   if (!c_function)
//     throw std::runtime_error(std::string("Undefined function ") + function_name);

  Node *function = new Node(OPCODE_FUNCTION);
  function->m_data.function.function = c_function;

  if (!c_function)
    eval.getUndefinedFunctionSet().insert(function_name);

  function->m_right = parseFunctionArg(eval, lparen + 1, rparen);

  return function;
}


Node *
parseIndex(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	lbrack,
  LexemVector::const_iterator	rbrack,
  LexemVector::const_iterator	to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER)
    throw std::runtime_error("syntax error parsing array");

  Node *index = new Node(OPCODE_RVALUE);
  index->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
  index->m_left = parseExpression(eval, lbrack + 1, rbrack);

  return index;
}


Node *
parseFunctionArg(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  if (from == to)
    return NULL;

  LexemVector::const_iterator it;
  for (it = from; it != to && (*it).getToken() != TOKEN_COMMA; ++it)
    ;

  Node *argument = new Node(OPCODE_ARGUMENT);
  argument->m_left = parseExpression(eval, from, it);
  if (it != to)
    argument->m_right = parseFunctionArg(eval, it + 1, to);
  return argument;
}


Node *
parseRValue(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  if (from + 1 != to)
    throw std::runtime_error(std::string("r-value not allowed following ") + (*from).getString());

  switch ((*from).getToken()) {
  case TOKEN_IDENTIFIER:
    {
      ConstantMap::iterator it = getConstantMap().find((*from).getString());
      if (it != getConstantMap().end()) {
	Node *constant = new Node(OPCODE_CONSTANT);
	constant->m_data.constant.value = (*it).second;
	return constant;
      }

      // Define a variable
      else {
	Node *variable = new Node(OPCODE_RVALUE);
	variable->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
	return variable;
      }
    }

  case TOKEN_REAL_CONSTANT:
  case TOKEN_INTEGER_CONSTANT:
    {
      Node *constant = new Node(OPCODE_CONSTANT);
      constant->m_data.constant.value = (*from).getValue<double>();
      return constant;
    }

  default:
    throw std::runtime_error("invalid rvalue");
  }
}

} // namespace Parser

Eval::Eval(
  VariableMap::Resolver &		resolver,
  const std::string &			expression)
  : m_variableMap(resolver),
    m_expression(expression),
    m_syntaxStatus(false),
    m_parseStatus(false),
    m_headNode(0)
{
  if (!m_expression.empty())
    parse();
}


Eval::~Eval() {
  delete m_headNode;
}


void
Eval::syntax()
{
  m_syntaxStatus = false;
  m_parseStatus = false;
  delete m_headNode;

  try {
    /* Validate the characters */
    LexemVector lex_vector = tokenize(m_expression);

    /* Call the multiparse routine to parse subexpressions */
    m_headNode = Parser::parseStatements(*this, lex_vector.begin(), lex_vector.end());

    m_syntaxStatus = true;
  }
  catch (std::runtime_error & /* x */) {
//     x << " while parsing expression: " << m_expression;
//     RuntimeDoomed() << x.what();
    throw;
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
//	for (iteration<UndefinedFunctionSet> it(m_undefinedFunctionSet); it; ++it)
	for (UndefinedFunctionSet::iterator it = m_undefinedFunctionSet.begin(); it != m_undefinedFunctionSet.end(); ++it)
	  strout << (*it) << std::endl;
	throw std::runtime_error(strout.str());
      }

      resolve();

      m_parseStatus = true;
    }
  }
  catch (std::runtime_error & /* x */) {
//     x << " while parsing expression: " << m_expression;
//     RuntimeDoomed() << x.what();
    throw;
  }
}


void
Eval::resolve()
{
  for (VariableMap::iterator it = m_variableMap.begin(); it != m_variableMap.end(); ++it)
    m_variableMap.getResolver().resolve(it);
}


double
Eval::evaluate() const
{
  /* Make sure it was parsed successfully */
  if (!m_parseStatus)
    throw std::runtime_error(std::string("Expression '") + m_expression + "' did not parse successfully");

  return m_headNode->eval();
}

} // namespace expreval
} // namespace stk
