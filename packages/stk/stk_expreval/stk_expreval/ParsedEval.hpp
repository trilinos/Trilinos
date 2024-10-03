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

#ifndef PARSEDEVAL_HPP
#define PARSEDEVAL_HPP

#include "Kokkos_Core.hpp"
#include "stk_expreval/Function.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_expreval/Eval.hpp"
#include "stk_expreval/ParsedEvalBase.hpp"
#include "stk_expreval/NgpNode.hpp"

namespace stk {
namespace expreval {

template <int RESULT_BUFFER_SIZE>
class ParsedEval : public ParsedEvalBase
{
  using NodeView = Kokkos::View<NgpNode*, stk::ngp::MemSpace>;

public:
  ParsedEval() = default;

  ParsedEval(const Eval& eval)
    : ParsedEvalBase(),
      m_numVariables(eval.get_variable_count()),
      m_requiredResultBufferSize(eval.get_result_buffer_size()),
      m_arrayOffsetType(eval.getArrayOffsetType()),
      m_firstNodeIndex(eval.get_first_node_index()),
      m_lastNodeIndex(eval.get_last_node_index()),
      m_deviceNodes(Kokkos::view_alloc(Kokkos::WithoutInitializing, "deviceNodeTree"), eval.get_node_count()),
      m_hostNodes(Kokkos::view_alloc(Kokkos::WithoutInitializing, "hostNodeTree"), eval.get_node_count())
  {
    STK_ThrowRequireMsg(eval.getParseStatus(), "Expression '" << eval.getExpression() << "' did not parse successfully");

    STK_ThrowRequireMsg(RESULT_BUFFER_SIZE >= eval.get_result_buffer_size(), "ParsedEval result buffer size ("
                    << RESULT_BUFFER_SIZE << ") must be at least " << eval.get_result_buffer_size()
                    << " to support expression");

    for (int i = 0; i < eval.get_node_count(); ++i) {
      m_hostNodes(i) = NgpNode(*(eval.get_node(i)));
    }
    Kokkos::deep_copy(m_deviceNodes, m_hostNodes);

    check_for_errors();
  }

  KOKKOS_DEFAULTED_FUNCTION ParsedEval(const ParsedEval&) = default;

  KOKKOS_DEFAULTED_FUNCTION virtual ~ParsedEval() override = default;

  virtual int get_result_buffer_size() override { return RESULT_BUFFER_SIZE; }

  STK_DEPRECATED_MSG("check_for_errors is now called by the constructor.  No need to call it yourself")
  void check_for_errors(bool /*will_run_on_device*/) const override
  {
    check_for_errors();
  }


  KOKKOS_INLINE_FUNCTION
  int get_num_variables() const { return m_numVariables; }

  KOKKOS_INLINE_FUNCTION
  int get_required_result_buffer_size() const { return m_requiredResultBufferSize; }

  KOKKOS_INLINE_FUNCTION
  Variable::ArrayOffset get_array_offset_type() const { return m_arrayOffsetType; }

  template <int MAX_BOUND_VARIABLES>
  KOKKOS_INLINE_FUNCTION
  double
  evaluate(DeviceVariableMap<MAX_BOUND_VARIABLES>& deviceVariableMap) const
  {
    double nodeResultBuffer[RESULT_BUFFER_SIZE];

    int nodeIndex = m_firstNodeIndex;

    KOKKOS_IF_ON_DEVICE((
      while (nodeIndex >= 0) {
        NgpNode & node = m_deviceNodes[nodeIndex];
        node.eval(deviceVariableMap, nodeResultBuffer);
        nodeIndex = node.getNextNodeIndex(nodeResultBuffer);
      }
      return (m_lastNodeIndex >= 0) ? m_deviceNodes[m_lastNodeIndex].getResult(nodeResultBuffer) : 0.0;
    ))

    KOKKOS_IF_ON_HOST((
      while (nodeIndex >= 0) {
        NgpNode & node = m_hostNodes[nodeIndex];
        node.eval(deviceVariableMap, nodeResultBuffer);
        nodeIndex = node.getNextNodeIndex(nodeResultBuffer);
      }
      return (m_lastNodeIndex >= 0) ? m_hostNodes[m_lastNodeIndex].getResult(nodeResultBuffer) : 0.0;
    ))
  }

private:
  template <int MAX_BOUND_VARIABLES>
  friend class DeviceVariableMap;

  void check_for_errors() const
  {
    for (size_t i=0; i < m_hostNodes.size(); ++i)
    {
      const NgpNode& node = m_hostNodes(i);
      if (node.m_opcode == OPCODE_FUNCTION)
      {
        FunctionType funcType = node.m_data.function.functionType;
        STK_ThrowRequireMsg(funcType != FunctionType::UNDEFINED, "user defined functions and system functions (rand(), time() etc.) are not supported by ParsedEval");
      }
    }
  }

  int m_numVariables;
  int m_requiredResultBufferSize;
  Variable::ArrayOffset m_arrayOffsetType;
  int m_firstNodeIndex;
  int m_lastNodeIndex;
  NodeView m_deviceNodes;
  NodeView::HostMirror m_hostNodes;
};

}
}

#endif // PARSEDEVAL_HPP
