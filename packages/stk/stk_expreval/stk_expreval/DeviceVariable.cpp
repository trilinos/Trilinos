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

#include "DeviceVariable.hpp"

namespace stk {
namespace expreval {

KOKKOS_FUNCTION
DeviceVariable::DeviceVariable()
  : m_type(Variable::Type::DOUBLE),
    m_size(1),
    m_stride(1),
    m_doublePtr(&m_doubleValue),
    m_doubleValue(0.0)
{
}

KOKKOS_FUNCTION
DeviceVariable::DeviceVariable(const Variable::Type variableType, int variableSize, int variableStride)
  : m_type(variableType),
    m_size(variableSize),
    m_stride(variableStride)
{
  switch (variableType) {
  case Variable::Type::DOUBLE:
    m_doublePtr = &m_doubleValue;
    m_doubleValue = 0.0;
    break;
  case Variable::Type::INTEGER:
    m_intPtr = &m_intValue;
    m_intValue = 0;
    break;
  }
}

KOKKOS_FUNCTION
DeviceVariable &
DeviceVariable::operator=(const DeviceVariable & deviceVariable)
{
  m_type = deviceVariable.m_type;
  m_size = deviceVariable.m_size;
  m_stride = deviceVariable.m_stride;

  switch (m_type) {
  case Variable::Type::DOUBLE:
    m_doublePtr = &m_doubleValue;
    m_doubleValue = deviceVariable.m_doubleValue;
    break;
  case Variable::Type::INTEGER:
    m_intPtr = &m_intValue;
    m_intValue = deviceVariable.m_intValue;
    break;
  }
  return *this;
}

KOKKOS_FUNCTION
double &
DeviceVariable::getArrayValue(int index, Variable::ArrayOffset arrayOffsetType) const
{
  NGP_ThrowRequireMsg(m_type == Variable::DOUBLE, "Only double arrays are allowed.");
  NGP_ThrowRequireMsg(m_doublePtr != nullptr, "Unbound array variable.");

  if (arrayOffsetType == Variable::ArrayOffset::ZERO_BASED_INDEX) {
    NGP_ThrowRequireMsg(index >= 0, "Provided variable array index is less than 0.");
    NGP_ThrowRequireMsg(index < m_size, "Provided variable array index exceeds array upper bound.");

    return m_doublePtr[index*m_stride];
  }
  else if (arrayOffsetType == Variable::ArrayOffset::ONE_BASED_INDEX) {
    NGP_ThrowRequireMsg(index >= 1, "Provided variable array index is less than 1.");
    NGP_ThrowRequireMsg(index <= m_size, "Provided variable array index exceeds array upper bound.");

    return m_doublePtr[(index-1)*m_stride];
  }
  else {
    NGP_ThrowErrorMsg("Invalid ArrayOffsetType.")
    return m_doublePtr[0];
  }
}

KOKKOS_FUNCTION
double
DeviceVariable::getValue() const
{
  NGP_ThrowRequireMsg(m_size == 1, "getValue Cannot access vector variable as a scalar.");

  switch (m_type) {
  case Variable::DOUBLE: {
    NGP_ThrowRequireMsg(m_doublePtr != nullptr, "Unbound double variable.");
    return *m_doublePtr;
  }
  case Variable::INTEGER: {
    NGP_ThrowRequireMsg(m_intPtr != nullptr, "Unbound integer variable.");
    return *m_intPtr;
  }
  }

  NGP_ThrowErrorMsg("Invalid variable type.");
  return *m_doublePtr;
}

KOKKOS_FUNCTION
void
DeviceVariable::bind(double& value_ref, int definedLength, int strideLength)
{
  m_type = Variable::DOUBLE;
  m_doublePtr = &value_ref;
  m_size = definedLength;
  m_stride = strideLength;
}

KOKKOS_FUNCTION
void
DeviceVariable::bind(int& value_ref, int definedLength, int strideLength)
{
  m_type = Variable::INTEGER;
  m_intPtr = &value_ref;
  m_size = definedLength;
  m_stride = strideLength;
}

KOKKOS_FUNCTION
DeviceVariable &
DeviceVariable::operator=(const double& value)
{
  NGP_ThrowRequireMsg(m_size == 1, "double = Cannot access vector variable as a scalar.");

  if (m_type == Variable::INTEGER) {
    *m_intPtr = static_cast<int>(value);
  }
  else if (m_type == Variable::DOUBLE) {
    *m_doublePtr = value;
  }
  return *this;
}

KOKKOS_FUNCTION
DeviceVariable&
DeviceVariable::operator=(const int& value)
{
  NGP_ThrowRequireMsg(m_size == 1, "int = Cannot access vector variable as a scalar.");

  if (m_type == Variable::INTEGER) {
    *m_intPtr = value;
  }
  else if (m_type == Variable::DOUBLE) {
    *m_doublePtr = static_cast<double>(value);
  }
  return *this;
}

}
}
