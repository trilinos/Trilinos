// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef SIERRA_ExprEval_h
#define SIERRA_ExprEval_h

#include <string>
#include <limits>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cctype>

#include <stk_expreval/Function.hpp>
#include <stk_expreval/Constants.hpp>
#include <stk_expreval/Variable.hpp>
#include <stk_util/util/ThreadLocalData.hpp>

namespace stk {
namespace expreval {

class Node;

/**
 * Class <b>Eval</b> parses and evaluates mathematical expressions.
 */
class Eval
{
public:
  typedef std::set<std::string> UndefinedFunctionSet;

  /**
   * Creates a new <b>Eval</b> instance.
   *
   * @param expr		a <b>std::string</b> const reference to the
   *				expression to be parsed.
   */
  Eval(VariableMap::Resolver &resolver = VariableMap::getDefaultResolver(), const std::string &expr = "", const Variable::ArrayOffset arrayOffsetType = Variable::ZERO_BASED_INDEX);
  explicit Eval(const std::string &expr, const Variable::ArrayOffset arrayOffsetType = Variable::ZERO_BASED_INDEX);

private:
  explicit Eval(const Eval &);
  Eval &operator=(const Eval &);

public:
  /**
   * Destroys a <b>Eval</b> instance.
   */
  ~Eval();

  /**
   * @brief Member function <b>getExpression</b> returns the original text of the expression.
   *
   * @return			a <b>std::string</b> const reference to the
   *				expression.
   */
  const std::string &getExpression() const {
    return m_expression;
  }

  /**
   * @brief Member function <b>setExpression</b> gives the evaluator a new
   * expression to parse.
   *
   * @param expression		a <b>std::string</b> const refernce to the
   *				expression to be parsed.
   *
   * @return			an <b>Eval</b> reference to this expression
   *				evaluator.
   */
  Eval &setExpression(const std::string &expression) {
    m_expression = expression;
    m_parseStatus = false;
    return *this;
  }

  /**
   * @brief Member function <b>getVariableMap</b> returns a reference to the
   * variable map of this expression.
   *
   * @return			a <b>VariableMap</b> reference to the variable map.
   */
  VariableMap &getVariableMap() {
    return m_variableMap.getMyThreadEntry();
  }

  /**
   * @brief Member function <b>getUndefinedFunctionSet</b> returns a reference to the
   * variable map of this expression.
   *
   * @return			a <b>VariableMap</b> reference to the variable map.
   */
  UndefinedFunctionSet &getUndefinedFunctionSet() {
    return m_undefinedFunctionSet;
  }

  /**
   * @brief Member function <b>getSyntaxStatus</b> returns true if the expression has
   * been syntaxd successfully.
   *
   * @return			a <b>bool</b> value of true if the expression has
   *				been syntaxd successfully.
   */
  bool getSyntaxStatus() const {
    return m_syntaxStatus;
  }

  /**
   * @brief Member function <b>getParseStatus</b> returns true if the expression has
   * been parsed successfully.
   *
   * @return			a <b>bool</b> value of true if the expression has
   *				been parsed successfully.
   */
  bool getParseStatus() const {
    return m_parseStatus;
  }

  Eval &setValue(const std::string &name, double* value, int definedLength) {
    auto& variableMap = m_variableMap.getMyThreadEntry();
    VariableMap::iterator it = variableMap.find(name);
    if (it != variableMap.end()) {
      (*it).second->bind(*value, definedLength);
    }
    return *this;
  }

  /**
   * @brief Member function <b>getValue</b> returns the value of the variable
   * specified by <b>name</b>.
   *
   * @param name		a <b>std::string</b> const reference to the variable's
   *				name.
   *
   * @return			a <b>double</b> value of the variable.
   */
  double getValue(const std::string &name) {
    auto& variableMap = m_variableMap.getMyThreadEntry();
    VariableMap::iterator it = variableMap.find(name);
    if (it == variableMap.end())
      throw std::runtime_error(std::string("Variable ") + name  + " not defined");
    return (*it).second->getValue();
  }

  /**
   * @brief Member function <b>newNode</b> allocates a new node.  The
   * new node is allocated on a node list so that it may be
   * deallocated properly on exception.
   *
   * @param op		a <b>int</b> value of the opcode for the node.
   *
   * @return            a <b>Node</b> pointer to the newly allocated node.
   */
  Node *newNode(int op);
  
  /**
   * @brief Member function <b>bindVariable</b> binds the variable to the address of
   * the specified value.  This address must remain in scope during the lifetime of the
   * variable or until the variable is rebound to a new address.
   *
   * @param name		a <b>std::string</b> const reference to the variable's name.
   *
   * @param value_ref		a <b>double</b> reference to be used for this variable.
   *
   * @return			an <b>Eval</b> reference to this expression evaluator.
   */
  Eval &bindVariable(const std::string &name, double &value_ref, int definedLength=std::numeric_limits<int>::max()) {
    auto& variableMap = m_variableMap.getMyThreadEntry();
    VariableMap::iterator it = variableMap.find(name);
    if (it != variableMap.end()) {
      (*it).second->bind(value_ref, definedLength);
    }
    return *this;
  }

  Eval &bindVariable(const std::string &name, ThreadLocalData<double> &value_ref, int definedLength=std::numeric_limits<int>::max()) {
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      auto& variableMap = m_variableMap.getMyThreadEntry();
      VariableMap::iterator it = variableMap.find(name);
      if (it != variableMap.end()) {
        (*it).second->bind(value_ref.getMyThreadEntry(), definedLength);
      }
    }
    return *this;
  }

  /**
   * @brief Member function <b>getVariable</b> returns a reference to the variable
   * specified by <b>name</b>.
   *
   * @param name		a <b>std::string</b> const reference to the variable's
   *				name.
   *
   * @return			a <b>Variable</b> reference to the specified
   *				variable.
   */
  Variable &getVariable(const std::string &name) {
    auto& variableMap = m_variableMap.getMyThreadEntry();
    VariableMap::iterator it = variableMap.find(name);
    if (it == variableMap.end())
      throw std::runtime_error(std::string("Variable ") + name  + " not defined");
    return *(*it).second;
  }

  /**
   * @brief Member function <b>parse</b> parses the expression.  If successful, the
   * parse status is set to true.
   *
   * @param expr		a <b>std::string</b> const reference to the
   *				expression to parse.
   *
   */
  void syntaxCheck(const std::string &expr) {
    setExpression(expr);
    syntax();
  }

  /**
   * @brief Member function <b>syntax</b> performs a syntax check on the current
   * expression.  If successful, the syntax status is set to true.
   */
  void syntax();

  /**
   * @brief Member function <b>parse</b> parses the expression.  If successful, the
   * parse status is set to true.
   *
   * @param expr		a <b>std::string</b> const reference to the expression to parse.
   *
   */
  void parse(const std::string &expr) {
    setExpression(expr);
    parse();
  }

  /**
   * @brief Member function <b>parse</b> parses the current expression.  If
   * successful, the parse status is set to true.
   */
  void parse();

  /**
   * @brief Member function <b>resolve</b> calls the variable name resolver for each
   * variable in the variable map.
   */
  void resolve();

  /**
   * @brief Member function <b>evaluate</b> evaluates the expression.
   *
   * @return			a <b>double</b> value of the result on the
   *				expression evaluation.
   */
  double evaluate() const;

  /**
   * @brief Member function <b>undefinedFunction</b> checks if any allocated node
   * represents an undefined (i.e. unknown at this point) function.
   *
   * @return			The returned <b>bool</b> is true if any allocated node
   *                            is an undefined function, as indicated by bool in the 
   *                            member union, which is set whenever a function is parsed.
   */
  bool undefinedFunction() const;

  Variable::ArrayOffset getArrayOffsetType() {return m_arrayOffsetType;}

private:
  stk::ThreadLocalData<VariableMap>		m_variableMap;		///< Variable map
  UndefinedFunctionSet	m_undefinedFunctionSet;	///< Vector of undefined functions

  std::string		m_expression;		///< Expression which was parsed.
  bool			m_syntaxStatus;		///< True if syntax is correct
  bool			m_parseStatus;		///< True if parsed successfully

  stk::ThreadLocalData<Node *>		m_headNode;		///< Head of compiled expression
  stk::ThreadLocalData<std::vector<Node *>> m_nodes;                ///< Allocated nodes
  Variable::ArrayOffset m_arrayOffsetType;      ///< Zero or one based array indexing
};

} // namespace expreval
} // namespace stk

#endif // SIERRA_ExprEval_h
