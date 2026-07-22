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

#include "stk_expreval/NgpNode.hpp"
#include "stk_expreval/Function.hpp"
#include "stk_expreval/Eval.hpp"

namespace stk {
namespace expreval {

NgpNode::NgpNode(const Node& node)
  : m_opcode(node.m_opcode),
    m_resultIdx(node.m_resultIdx),
    m_currentNodeIndex(node.m_currentNodeIndex),
    m_nextNodeIndex(node.m_nextNodeIndex),
    m_ternaryTrueNextNodeIndex(node.m_ternaryTrueNextNodeIndex),
    m_ternaryFalseNextNodeIndex(node.m_ternaryFalseNextNodeIndex),
    m_leftNodeIndex((node.m_left != nullptr) ? node.m_left->m_currentNodeIndex : -1),
    m_rightNodeIndex((node.m_right != nullptr) ? node.m_right->m_currentNodeIndex : -1),
    m_ternaryOtherNodeIndex((node.m_ternaryOther != nullptr) ? node.m_ternaryOther->m_currentNodeIndex : -1),
    m_fpErrorBehavior(node.m_owner->get_fp_error_behavior())
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

}
}
