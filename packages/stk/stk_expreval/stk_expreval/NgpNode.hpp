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

#ifndef NGPNODE_HPP
#define NGPNODE_HPP

#include "Kokkos_Core.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_expreval/Variable.hpp"
#include "stk_expreval/Node.hpp"

namespace stk {
namespace expreval {

static constexpr int DEFAULT_MAX_BOUND_VARIABLES = 16;

class Node;

template <int MAX_BOUND_VARIABLES=DEFAULT_MAX_BOUND_VARIABLES>
class DeviceVariableMap;

class NgpNode
{
public:
  enum { MAXIMUM_NUMBER_OF_OVERLOADED_FUNCTION_NAMES = 5 };
  enum { MAXIMUM_FUNCTION_NAME_LENGTH = 32 };

  KOKKOS_FUNCTION
  NgpNode();

  explicit NgpNode(const Node& node);

  KOKKOS_DEFAULTED_FUNCTION
  NgpNode(const NgpNode &) = default;

  KOKKOS_DEFAULTED_FUNCTION
  NgpNode &operator=(const NgpNode &) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~NgpNode() = default;

  template <int MAX_BOUND_VARIABLES>
  KOKKOS_FUNCTION
  double
  eval(DeviceVariableMap<MAX_BOUND_VARIABLES> & deviceVariableMap, double resultBuffer[])
  {
    switch (m_opcode) {
    case OPCODE_STATEMENT: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_CONSTANT: {
      setResult(resultBuffer) = m_data.constant.value;
      break;
    }
    case OPCODE_RVALUE: {
      if (get_left_node()) {
        setResult(resultBuffer) = deviceVariableMap[m_data.variable.variableIndex].getArrayValue(get_left_node()->getResult(resultBuffer), deviceVariableMap.get_array_offset_type());
      }
      else {
        setResult(resultBuffer) = deviceVariableMap[m_data.variable.variableIndex].getValue();
      }
      break;
    }
    case OPCODE_MULTIPLY: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer)*get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_EXPONENIATION: {
      setResult(resultBuffer) = std::pow(get_left_node()->getResult(resultBuffer),get_right_node()->getResult(resultBuffer));
      break;
    }
    case OPCODE_DIVIDE: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer)/get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_MODULUS: {
      setResult(resultBuffer) = std::fmod(get_left_node()->getResult(resultBuffer), get_right_node()->getResult(resultBuffer));
      break;
    }
    case OPCODE_ADD: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) + get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_SUBTRACT: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) - get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_EQUAL: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) == get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_NOT_EQUAL: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) != get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_LESS: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) < get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_GREATER: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) > get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_LESS_EQUAL: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) <= get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_GREATER_EQUAL: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer) >= get_right_node()->getResult(resultBuffer) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_LOGICAL_AND: {
      double left = get_left_node()->getResult(resultBuffer);
      double right = get_right_node()->getResult(resultBuffer);
      setResult(resultBuffer) = (left != 0.0) && (right != 0.0) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_LOGICAL_OR: {
      double left = get_left_node()->getResult(resultBuffer);
      double right = get_right_node()->getResult(resultBuffer);
      setResult(resultBuffer) = (left != 0.0) || (right != 0.0) ? 1.0 : 0.0;
      break;
    }
    case OPCODE_TERNARY_PREDICATE: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_TERNARY_JOIN: {
      if (get_ternary_other_node()->getResult(resultBuffer) == 0.0) {
        setResult(resultBuffer) = get_right_node()->getResult(resultBuffer);
      }
      else {
        setResult(resultBuffer) = get_left_node()->getResult(resultBuffer);
      }
      break;
    }
    case OPCODE_UNARY_MINUS: {
      setResult(resultBuffer) = -get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_UNARY_NOT: {
      setResult(resultBuffer) = get_right_node()->getResult(resultBuffer) == 0.0 ? 1.0 : 0.0;
      break;
    }
    case OPCODE_ASSIGN: {
      if (get_left_node()) {
        deviceVariableMap[m_data.variable.variableIndex].getArrayValue(get_left_node()->getResult(resultBuffer),
                                                                       deviceVariableMap.get_array_offset_type()) = get_right_node()->getResult(resultBuffer);
      }
      else {
        deviceVariableMap[m_data.variable.variableIndex] = get_right_node()->getResult(resultBuffer);
      }
      setResult(resultBuffer) = get_right_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_ARGUMENT: {
      setResult(resultBuffer) = get_left_node()->getResult(resultBuffer);
      break;
    }
    case OPCODE_FUNCTION: {
      double arguments[20];

      int argumentCount = 0;
      for (const NgpNode* arg = get_right_node(); arg; arg = arg->get_right_node()) {
        arguments[argumentCount++] = arg->getResult(resultBuffer);
      }

      setResult(resultBuffer) = evaluate_function(argumentCount, arguments);
      break;
    }
    default: { //  Unknown opcode
      setResult(resultBuffer) = 0.0;
    }
    }

    return 0.0;
  }

  KOKKOS_INLINE_FUNCTION
  int getNextNodeIndex(double resultBuffer[]) {
    if (m_opcode == OPCODE_TERNARY_PREDICATE) {
      if (setResult(resultBuffer) == 0.0) {
        return m_ternaryFalseNextNodeIndex;
      }
      else {
        return m_ternaryTrueNextNodeIndex;
      }
    }
    return m_nextNodeIndex;
  }

  KOKKOS_INLINE_FUNCTION
  double getResult(double resultBuffer[]) const { return resultBuffer[m_resultIdx]; }

  KOKKOS_INLINE_FUNCTION
  double& setResult(double resultBuffer[]) const { return resultBuffer[m_resultIdx]; }

  KOKKOS_INLINE_FUNCTION
  const NgpNode* get_node(int nodeIndex) const {
    return nodeIndex != -1 ? (this + (nodeIndex - m_currentNodeIndex)) : nullptr;
  }

  KOKKOS_INLINE_FUNCTION const NgpNode* get_left_node() const { return get_node(m_leftNodeIndex); }
  KOKKOS_INLINE_FUNCTION const NgpNode* get_right_node() const { return get_node(m_rightNodeIndex); }
  KOKKOS_INLINE_FUNCTION const NgpNode* get_ternary_other_node() const { return get_node(m_ternaryOtherNodeIndex); }

  Opcode m_opcode;

  union _data
  {
    struct _constant
    {
      double value;
    } constant;

    struct _variable
    {
      int variableIndex;
      Variable::Type variableType;
      int variableSize;
    } variable;

    struct _function
    {
      FunctionType functionType;
    } function;
  } m_data;

  int m_resultIdx;
  int m_currentNodeIndex;
  int m_nextNodeIndex;
  int m_ternaryTrueNextNodeIndex;
  int m_ternaryFalseNextNodeIndex;
  int m_leftNodeIndex;
  int m_rightNodeIndex;
  int m_ternaryOtherNodeIndex;

  KOKKOS_FUNCTION
  double evaluate_function(int argumentCount, double* arguments) const;
};

}
}

#endif // NGPNODE_HPP
