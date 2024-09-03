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

#include "stk_expreval/Node.hpp"
#include "stk_expreval/Eval.hpp"
#include "stk_expreval/Constants.hpp"
#include <string>

namespace stk {
namespace expreval {

Node::Node(Opcode opcode, Eval* owner)
  : m_opcode(opcode),
    m_resultIdx(-1),
    m_currentNodeIndex(-1),
    m_nextNodeIndex(-1),
    m_ternaryTrueNextNodeIndex(-1),
    m_ternaryFalseNextNodeIndex(-1),
    m_left(nullptr),
    m_right(nullptr),
    m_ternaryOther(nullptr),
    m_owner(owner),
    m_hasBeenEvaluated(false)
{
  m_data.function.function = nullptr;
}

int
Node::getNextNodeIndex() {
  if (m_opcode == OPCODE_TERNARY_PREDICATE) {
    if (getResult() == s_false) {
      return m_ternaryFalseNextNodeIndex;
    }
    else {
      return m_ternaryTrueNextNodeIndex;
    }
  }
  return m_nextNodeIndex;
}

double Node::getResult() const {
  STK_ThrowAssertMsg(m_hasBeenEvaluated, "Requesting node result before it has been computed.");
  return m_owner->get_result_buffer_value(m_resultIdx);
}

double& Node::setResult() {
  return m_owner->get_result_buffer_value(m_resultIdx);
}

void
Node::eval()
{
  switch (m_opcode) {
  case OPCODE_STATEMENT: {
    setResult() = m_left->getResult();
    break;
  }
  case OPCODE_CONSTANT: {
    setResult() = m_data.constant.value;
    break;
  }
  case OPCODE_RVALUE: {
    /* Directly access the variable */
    if (m_left) {
      setResult() = m_data.variable.variable->getArrayValue(m_left->getResult(), m_owner->getArrayOffsetType());
    } else {
      setResult() = m_data.variable.variable->getValue();
    }
    break;
  }
  case OPCODE_MULTIPLY: {
    setResult() = m_left->getResult()*m_right->getResult();
    break;
  }
  case OPCODE_EXPONENIATION: {
    setResult() = std::pow(m_left->getResult(),m_right->getResult());
    break;
  }
  case OPCODE_DIVIDE: {
    setResult() = m_left->getResult()/m_right->getResult();
    break;
  }
  case OPCODE_MODULUS: {
    setResult() = std::fmod(m_left->getResult(), m_right->getResult());
    break;
  }
  case OPCODE_ADD: {
    setResult() = m_left->getResult() + m_right->getResult();
    break;
  }
  case OPCODE_SUBTRACT: {
    setResult() = m_left->getResult() - m_right->getResult();
    break;
  }
  case OPCODE_EQUAL: {
    setResult() = m_left->getResult() == m_right->getResult() ? s_true : s_false;
    break;
  }
  case OPCODE_NOT_EQUAL: {
    setResult() = m_left->getResult() != m_right->getResult() ? s_true : s_false;
    break;
  }
  case OPCODE_LESS: {
    setResult() = m_left->getResult() < m_right->getResult() ? s_true : s_false;
    break;
  }
  case OPCODE_GREATER: {
    setResult() = m_left->getResult() > m_right->getResult() ? s_true : s_false;
    break;
  }
  case OPCODE_LESS_EQUAL: {
    setResult() = m_left->getResult() <= m_right->getResult() ? s_true : s_false;
    break;
  }
  case OPCODE_GREATER_EQUAL: {
    setResult() = m_left->getResult() >= m_right->getResult() ? s_true : s_false;
    break;
  }
  case OPCODE_LOGICAL_AND: {
    double left = m_left->getResult();
    double right = m_right->getResult();
    setResult() = (left != s_false) && (right != s_false) ? s_true : s_false;
    break;
  }
  case OPCODE_LOGICAL_OR: {
    double left = m_left->getResult();
    double right = m_right->getResult();
    setResult() = (left != s_false) || (right != s_false) ? s_true : s_false;
    break;
  }
  case OPCODE_TERNARY_PREDICATE: {
    setResult() = m_left->getResult();
    break;
  }
  case OPCODE_TERNARY_JOIN: {
    if (m_ternaryOther->getResult() == s_false) {
      setResult() = m_right->getResult();
    }
    else {
      setResult() = m_left->getResult();
    }
    break;
  }
  case OPCODE_UNARY_MINUS: {
    setResult() = -m_right->getResult();
    break;
  }
  case OPCODE_UNARY_NOT: {
    setResult() = m_right->getResult() == s_false ? s_true : s_false;
    break;
  }
  case OPCODE_ASSIGN: {
    if (m_left) {
      m_data.variable.variable->assignArrayValue(m_left->getResult(),
                                                 m_owner->getArrayOffsetType(),
                                                 m_right->getResult());
    }
    else {
      *m_data.variable.variable = m_right->getResult();
    }
    setResult() = m_right->getResult();
    break;
  }
  case OPCODE_ARGUMENT: {
    setResult() = m_left->getResult();
    break;
  }
  case OPCODE_FUNCTION: {
    double argv[MAXIMUM_NUMBER_OF_FUNCTION_ARGS];

    int argc = 0;
    for (Node *arg = m_right; arg; arg = arg->m_right) {
      argv[argc++] = arg->getResult();
    }

    setResult() = (*m_data.function.function)(argc, argv);
    break;
  }
  default: {
    STK_ThrowErrorMsg("Unknown OpCode (" + std::to_string(m_opcode) + ")");
  }
  }

  m_hasBeenEvaluated = true;
}

void
Node::computeNodeWeight(NodeWeightMap & nodeWeights)
{
  switch (m_opcode) {
  case OPCODE_STATEMENT: {
    for (Node *statement = this; statement; statement = statement->m_right) {
      statement->m_left->computeNodeWeight(nodeWeights);
      nodeWeights[statement] = nodeWeights[m_left];
    }
    break;
  }

  case OPCODE_CONSTANT: {
    nodeWeights[this] = 1;
    break;
  }

  case OPCODE_RVALUE: {
    nodeWeights[this] = 1;
    if (m_left) {
      m_left->computeNodeWeight(nodeWeights);
      nodeWeights[this] += nodeWeights.at(m_left);
    }
    break;
  }

  case OPCODE_MULTIPLY:
  case OPCODE_EXPONENIATION:
  case OPCODE_DIVIDE:
  case OPCODE_MODULUS:
  case OPCODE_ADD:
  case OPCODE_SUBTRACT:
  case OPCODE_EQUAL:
  case OPCODE_NOT_EQUAL:
  case OPCODE_LESS:
  case OPCODE_GREATER:
  case OPCODE_LESS_EQUAL:
  case OPCODE_GREATER_EQUAL:
  case OPCODE_LOGICAL_AND:
  case OPCODE_LOGICAL_OR: {
    m_left->computeNodeWeight(nodeWeights);
    m_right->computeNodeWeight(nodeWeights);
    nodeWeights[this] = nodeWeights.at(m_left) + nodeWeights.at(m_right);
    break;
  }

  case OPCODE_TERNARY_PREDICATE: {
    m_left->computeNodeWeight(nodeWeights);
    m_ternaryOther->m_left->computeNodeWeight(nodeWeights);
    m_ternaryOther->m_right->computeNodeWeight(nodeWeights);

    nodeWeights[this] = nodeWeights.at(m_left) + std::max(nodeWeights.at(m_ternaryOther->m_left),
                                                          nodeWeights.at(m_ternaryOther->m_right));
    break;
  }

  case OPCODE_TERNARY_JOIN: {
    nodeWeights[this] = std::max(nodeWeights.at(m_left), nodeWeights.at(m_right));
    break;
  }

  case OPCODE_UNARY_MINUS:
  case OPCODE_UNARY_NOT: {
    m_right->computeNodeWeight(nodeWeights);
    nodeWeights[this] = nodeWeights[m_right];
    break;
  }

  case OPCODE_ASSIGN: {
    m_right->computeNodeWeight(nodeWeights);
    nodeWeights[this] = nodeWeights.at(m_right);
    if (m_left) {
      m_left->computeNodeWeight(nodeWeights);
      nodeWeights[this] += nodeWeights.at(m_left);
    }
    break;
  }

  case OPCODE_FUNCTION: {
    int functionWeight = 1;
    for (Node *arg = m_right; arg; arg = arg->m_right) {
      arg->m_left->computeNodeWeight(nodeWeights);
      nodeWeights[arg] = nodeWeights.at(arg->m_left);
      functionWeight = std::max(functionWeight, nodeWeights.at(arg->m_left));
    }
    nodeWeights[this] = functionWeight;
    break;
  }

  default: { // Unknown opcode
    throw expression_evaluation_exception();
  }
  }
}

void
Node::evalTrace(const NodeWeightMap & nodeWeights, EvalNodesType & evaluationNodes)
{
  switch (m_opcode) {
  case OPCODE_STATEMENT: {
    for (Node *statement = this; statement; statement = statement->m_right) {
      statement->m_left->evalTrace(nodeWeights, evaluationNodes);
      evaluationNodes.back()->m_nextNodeIndex = statement->m_currentNodeIndex;
      evaluationNodes.push_back(statement);
    }
    break;
  }

  case OPCODE_CONSTANT: {
    break;
  }

  case OPCODE_RVALUE: {
    if (m_left) {
      m_left->evalTrace(nodeWeights, evaluationNodes);
    }
    break;
  }

  case OPCODE_MULTIPLY:
  case OPCODE_EXPONENIATION:
  case OPCODE_DIVIDE:
  case OPCODE_MODULUS:
  case OPCODE_ADD:
  case OPCODE_SUBTRACT:
  case OPCODE_EQUAL:
  case OPCODE_NOT_EQUAL:
  case OPCODE_LESS:
  case OPCODE_GREATER:
  case OPCODE_LESS_EQUAL:
  case OPCODE_GREATER_EQUAL:
  case OPCODE_LOGICAL_AND:
  case OPCODE_LOGICAL_OR: {
    if (nodeWeights.at(m_left) >= nodeWeights.at(m_right)) {
      m_left->evalTrace(nodeWeights, evaluationNodes);
      m_right->evalTrace(nodeWeights, evaluationNodes);
    }
    else {
      m_right->evalTrace(nodeWeights, evaluationNodes);
      m_left->evalTrace(nodeWeights, evaluationNodes);
    }
    break;
  }

  case OPCODE_TERNARY_PREDICATE: {
    m_left->evalTrace(nodeWeights, evaluationNodes);
    evaluationNodes.back()->m_nextNodeIndex = m_currentNodeIndex;

    evaluationNodes.push_back(this);

    const size_t firstNodeTrueBranch = evaluationNodes.size();
    m_ternaryOther->m_left->evalTrace(nodeWeights, evaluationNodes);
    m_ternaryTrueNextNodeIndex = evaluationNodes[firstNodeTrueBranch]->m_currentNodeIndex;

    const size_t firstNodeFalseBranch = evaluationNodes.size();
    m_ternaryOther->m_right->evalTrace(nodeWeights, evaluationNodes);
    m_ternaryFalseNextNodeIndex = evaluationNodes[firstNodeFalseBranch]->m_currentNodeIndex;

    const size_t lastNodeTrueBranch = firstNodeFalseBranch - 1;
    evaluationNodes[lastNodeTrueBranch]->m_nextNodeIndex = m_ternaryOther->m_currentNodeIndex;
    break;
  }

  case OPCODE_TERNARY_JOIN:
    break;

  case OPCODE_UNARY_MINUS:
  case OPCODE_UNARY_NOT: {
    m_right->evalTrace(nodeWeights, evaluationNodes);
    break;
  }

  case OPCODE_ASSIGN: {
    if (m_left) {
      if (nodeWeights.at(m_left) >= nodeWeights.at(m_right)) {
        m_left->evalTrace(nodeWeights, evaluationNodes);
        m_right->evalTrace(nodeWeights, evaluationNodes);
      }
      else {
        m_right->evalTrace(nodeWeights, evaluationNodes);
        m_left->evalTrace(nodeWeights, evaluationNodes);
      }
    }
    else {
      m_right->evalTrace(nodeWeights, evaluationNodes);
    }
    break;
  }

  case OPCODE_FUNCTION: {
    for (Node *arg = m_right; arg; arg = arg->m_right) {
      arg->m_left->evalTrace(nodeWeights, evaluationNodes);
      evaluationNodes.back()->m_nextNodeIndex = arg->m_currentNodeIndex;
      evaluationNodes.push_back(arg);
    }
    break;
  }

  default: { // Unknown opcode
    throw expression_evaluation_exception();
  }
  }

  if (m_opcode == OPCODE_TERNARY_PREDICATE) {
    m_ternaryOther->evalTrace(nodeWeights, evaluationNodes);
  }
  else if (m_opcode != OPCODE_STATEMENT) {
    if (!evaluationNodes.empty()) {
      evaluationNodes.back()->m_nextNodeIndex = m_currentNodeIndex;
    }
    evaluationNodes.push_back(this);
  }
}

}
}
