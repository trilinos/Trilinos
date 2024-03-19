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

#ifndef DEVICEVARIABLEMAP_HPP
#define DEVICEVARIABLEMAP_HPP

#include "Kokkos_Core.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_expreval/DeviceVariable.hpp"
#include "stk_expreval/ParsedEval.hpp"

namespace stk {
namespace expreval {

template <int MAX_BOUND_VARIABLES>
class DeviceVariableMap
{
public:
  KOKKOS_DEFAULTED_FUNCTION ~DeviceVariableMap() = default;

  template <int RESULT_BUFFER_SIZE>
  KOKKOS_INLINE_FUNCTION
  explicit DeviceVariableMap(const ParsedEval<RESULT_BUFFER_SIZE> & parsedEval)
    : m_arrayOffsetType(parsedEval.m_arrayOffsetType)
  {
    STK_NGP_ThrowRequireMsg(parsedEval.get_num_variables() <= MAX_BOUND_VARIABLES,
                        "stk::expreval::DeviceVariableMap is not large enough to hold all variables");

    KOKKOS_IF_ON_DEVICE((
      const auto & deviceNodes = parsedEval.m_deviceNodes;
      for (unsigned nodeIndex = 0u; nodeIndex < deviceNodes.extent(0); ++nodeIndex) {
        if (deviceNodes(nodeIndex).m_opcode == OPCODE_ASSIGN) {
          m_deviceVariableMap[deviceNodes(nodeIndex).m_data.variable.variableIndex] = DeviceVariable(deviceNodes(nodeIndex).m_data.variable.variableType,
                                                                                                    deviceNodes(nodeIndex).m_data.variable.variableSize);
        }
      }
    ))

    KOKKOS_IF_ON_HOST((
      const auto & hostNodes = parsedEval.m_hostNodes;
      for (unsigned nodeIndex = 0u; nodeIndex < hostNodes.extent(0); ++nodeIndex) {
        if (hostNodes(nodeIndex).m_opcode == OPCODE_ASSIGN) {
          m_deviceVariableMap[hostNodes(nodeIndex).m_data.variable.variableIndex] = DeviceVariable(hostNodes(nodeIndex).m_data.variable.variableType,
                                                                                                   hostNodes(nodeIndex).m_data.variable.variableSize);
        }
      }
    ))
  }

  KOKKOS_INLINE_FUNCTION
  void bind(int variableIndex, const double& value_ref, int definedLength=1, int strideLength=1) {
    m_deviceVariableMap[variableIndex].bind(value_ref, definedLength, strideLength);
  }

  KOKKOS_INLINE_FUNCTION
  void bind(int variableIndex, double& value_ref, int definedLength=1, int strideLength=1) {
    m_deviceVariableMap[variableIndex].bind(value_ref, definedLength, strideLength);
  }

  KOKKOS_INLINE_FUNCTION
  void bind(int variableIndex, const int& value_ref, int definedLength=1, int strideLength=1) {
    m_deviceVariableMap[variableIndex].bind(value_ref, definedLength, strideLength);
  }

  KOKKOS_INLINE_FUNCTION
  void bind(int variableIndex, int& value_ref, int definedLength=1, int strideLength=1) {
    m_deviceVariableMap[variableIndex].bind(value_ref, definedLength, strideLength);
  }

  KOKKOS_INLINE_FUNCTION
  DeviceVariable & operator[](int index) { return m_deviceVariableMap[index]; }

  KOKKOS_INLINE_FUNCTION
  Variable::ArrayOffset get_array_offset_type() { return m_arrayOffsetType; }

private:
  DeviceVariable m_deviceVariableMap[MAX_BOUND_VARIABLES];
  Variable::ArrayOffset m_arrayOffsetType;
};

}
}

#endif // DEVICEVARIABLEMAP_HPP
