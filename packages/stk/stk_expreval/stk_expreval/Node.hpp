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

#ifndef NODE_HPP
#define NODE_HPP

#include "stk_expreval/Function.hpp"
#include <ostream>
#include <map>
#include <vector>
#include <exception>

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

namespace stk {
namespace expreval {

class Node;
class Eval;
class Variable;
class CFunctionBase;

enum Opcode {
  OPCODE_UNDEFINED,
  OPCODE_CONSTANT,
  OPCODE_RVALUE,
  OPCODE_STATEMENT,
  OPCODE_ARGUMENT,

  OPCODE_TERNARY_PREDICATE,
  OPCODE_TERNARY_JOIN,

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

  OPCODE_EXPONENIATION,

  OPCODE_ASSIGN
};

inline std::ostream & operator<<(std::ostream & stream, Opcode opcode) {
  static std::vector<std::string> opcodeNames {"OPCODE_UNDEFINED",
                                               "OPCODE_CONSTANT",
                                               "OPCODE_RVALUE",
                                               "OPCODE_STATEMENT",
                                               "OPCODE_ARGUMENT",

                                               "OPCODE_TERNARY_PREDICATE",
                                               "OPCODE_TERNARY_JOIN",

                                               "OPCODE_MULTIPLY",
                                               "OPCODE_DIVIDE",
                                               "OPCODE_MODULUS",
                                               "OPCODE_ADD",
                                               "OPCODE_SUBTRACT",
                                               "OPCODE_UNARY_MINUS",
                                               "OPCODE_FUNCTION",

                                               "OPCODE_EQUAL",
                                               "OPCODE_NOT_EQUAL",
                                               "OPCODE_LESS",
                                               "OPCODE_GREATER",
                                               "OPCODE_LESS_EQUAL",
                                               "OPCODE_GREATER_EQUAL",

                                               "OPCODE_UNARY_NOT",
                                               "OPCODE_LOGICAL_AND",
                                               "OPCODE_LOGICAL_OR",

                                               "OPCODE_EXPONENIATION",

                                               "OPCODE_ASSIGN"};
  return stream << opcodeNames[opcode];
}

using EvalNodesType = std::vector<Node*>;
using NodeWeightMap = std::map<Node*, int>;

class Node
{
public:
  static constexpr int MAXIMUM_FUNCTION_NAME_LENGTH = 32;
  static constexpr int MAXIMUM_NUMBER_OF_FUNCTION_ARGS = 20;

  explicit Node(Opcode opcode, Eval* owner);

private:
  explicit Node(const Node &);
  Node &operator=(const Node &);

public:
  ~Node()
  {}

  void eval();

  void computeNodeWeight(NodeWeightMap & nodeWeights);
  void evalTrace(const NodeWeightMap & nodeWeights, EvalNodesType & evaluationNodes);

  int countNumFunctionArgs();

  int getNextNodeIndex();

  double getResult() const;

  double& setResult();

  const Opcode m_opcode;

  union _data
  {
    struct _constant
    {
      double value;
    } constant;

    struct _variable
    {
      Variable * variable;
    } variable;

    struct _function
    {
      CFunctionBase * function;
      char functionName[MAXIMUM_FUNCTION_NAME_LENGTH];
      FunctionType functionType;
    } function;
  } m_data;

  int m_resultIdx;
  int m_currentNodeIndex;
  int m_nextNodeIndex;
  int m_ternaryTrueNextNodeIndex;
  int m_ternaryFalseNextNodeIndex;
  Node * m_left;
  Node * m_right;
  Node * m_ternaryOther;
  Eval * m_owner;
  bool m_hasBeenEvaluated;
};

}
}

#endif // NODE_HPP
