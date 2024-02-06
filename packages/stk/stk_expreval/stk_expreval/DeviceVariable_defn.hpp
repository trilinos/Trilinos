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

#ifndef stk_expreval_DeviceVariable_defn_hpp
#define stk_expreval_DeviceVariable_defn_hpp

#include "DeviceVariable_decl.hpp"

namespace stk {
namespace expreval {

KOKKOS_INLINE_FUNCTION
DeviceVariable::DeviceVariable()
  : m_type(Variable::Type::DOUBLE),
    m_size(1),
    m_stride(1),
    m_doublePtr(&m_doubleValue),
    m_doubleValue(0.0),
    m_isModifiable(true)
{
}

KOKKOS_INLINE_FUNCTION
DeviceVariable::DeviceVariable(const Variable::Type variableType, int variableSize, int variableStride)
  : m_type(variableType),
    m_size(variableSize),
    m_stride(variableStride),
    m_isModifiable(true)
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

KOKKOS_INLINE_FUNCTION
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

KOKKOS_INLINE_FUNCTION
double
DeviceVariable::getArrayValue(int index, Variable::ArrayOffset arrayOffsetType) const
{
  if (m_type == Variable::Type::DOUBLE) {
    STK_NGP_ThrowRequireMsg(m_doublePtr != nullptr, "Invalid array variable.");

    if (arrayOffsetType == Variable::ArrayOffset::ZERO_BASED_INDEX) {
      STK_NGP_ThrowRequireMsg(index >= 0, "Provided variable array index is less than 0.");
      STK_NGP_ThrowRequireMsg(index < m_size, "Provided variable array index exceeds array upper bound.");

      return m_doublePtr[index*m_stride];
    }
    else if (arrayOffsetType == Variable::ArrayOffset::ONE_BASED_INDEX) {
      STK_NGP_ThrowRequireMsg(index >= 1, "Provided variable array index is less than 1.");
      STK_NGP_ThrowRequireMsg(index <= m_size, "Provided variable array index exceeds array upper bound.");

      return m_doublePtr[(index-1)*m_stride];
    }
    else {
      STK_NGP_ThrowErrorMsg("Invalid ArrayOffsetType.");
      return m_doublePtr[0];
    }
  }
  else {
    STK_NGP_ThrowRequireMsg(m_intPtr != nullptr, "Invalid array variable.");

    if (arrayOffsetType == Variable::ArrayOffset::ZERO_BASED_INDEX) {
      STK_NGP_ThrowRequireMsg(index >= 0, "Provided variable array index is less than 0.");
      STK_NGP_ThrowRequireMsg(index < m_size, "Provided variable array index exceeds array upper bound.");

      return static_cast<double>(m_intPtr[index*m_stride]);
    }
    else if (arrayOffsetType == Variable::ArrayOffset::ONE_BASED_INDEX) {
      STK_NGP_ThrowRequireMsg(index >= 1, "Provided variable array index is less than 1.");
      STK_NGP_ThrowRequireMsg(index <= m_size, "Provided variable array index exceeds array upper bound.");

      return static_cast<double>(m_intPtr[(index-1)*m_stride]);
    }
    else {
      STK_NGP_ThrowErrorMsg("Invalid ArrayOffsetType.");
      return static_cast<double>(m_intPtr[0]);
    }
  }
}

KOKKOS_INLINE_FUNCTION
void
DeviceVariable::assignArrayValue(int index, Variable::ArrayOffset arrayOffsetType, double value) const
{
  STK_NGP_ThrowRequireMsg(m_isModifiable, "Cannot modify const array.");

  if (m_type == Variable::Type::DOUBLE) {
    STK_NGP_ThrowRequireMsg(m_doublePtr != nullptr, "Invalid array variable.");

    if (arrayOffsetType == Variable::ArrayOffset::ZERO_BASED_INDEX) {
      STK_NGP_ThrowRequireMsg(index >= 0, "Provided variable array index is less than 0.");
      STK_NGP_ThrowRequireMsg(index < m_size, "Provided variable array index exceeds array upper bound.");

      const_cast<double&>(m_doublePtr[index*m_stride]) = value;
    }
    else if (arrayOffsetType == Variable::ArrayOffset::ONE_BASED_INDEX) {
      STK_NGP_ThrowRequireMsg(index >= 1, "Provided variable array index is less than 1.");
      STK_NGP_ThrowRequireMsg(index <= m_size, "Provided variable array index exceeds array upper bound.");

      const_cast<double&>(m_doublePtr[(index-1)*m_stride]) = value;
    }
    else {
      STK_NGP_ThrowErrorMsg("Invalid ArrayOffsetType.");
    }
  }
  else {
    STK_NGP_ThrowRequireMsg(m_intPtr != nullptr, "Invalid array variable.");

    if (arrayOffsetType == Variable::ArrayOffset::ZERO_BASED_INDEX) {
      STK_NGP_ThrowRequireMsg(index >= 0, "Provided variable array index is less than 0.");
      STK_NGP_ThrowRequireMsg(index < m_size, "Provided variable array index exceeds array upper bound.");

      const_cast<int&>(m_intPtr[index*m_stride]) = static_cast<int>(value);
    }
    else if (arrayOffsetType == Variable::ArrayOffset::ONE_BASED_INDEX) {
      STK_NGP_ThrowRequireMsg(index >= 1, "Provided variable array index is less than 1.");
      STK_NGP_ThrowRequireMsg(index <= m_size, "Provided variable array index exceeds array upper bound.");

      const_cast<int&>(m_intPtr[(index-1)*m_stride]) = static_cast<int>(value);
    }
    else {
      STK_NGP_ThrowErrorMsg("Invalid ArrayOffsetType.");
    }
  }
}

KOKKOS_INLINE_FUNCTION
double
DeviceVariable::getValue() const
{
  STK_NGP_ThrowRequireMsg(m_size == 1, "getValue Cannot access vector variable as a scalar.");

  switch (m_type) {
  case Variable::DOUBLE: {
    STK_NGP_ThrowRequireMsg(m_doublePtr != nullptr, "Double variable does not have a valid value.");
    return *m_doublePtr;
  }
  case Variable::INTEGER: {
    STK_NGP_ThrowRequireMsg(m_intPtr != nullptr, "Integer variable does not have a valid value.");
    return *m_intPtr;
  }
  }

  STK_NGP_ThrowErrorMsg("Invalid variable type.");
  return *m_doublePtr;
}

KOKKOS_INLINE_FUNCTION
void
DeviceVariable::bind(const double& value_ref, int definedLength, int strideLength)
{
  m_type = Variable::DOUBLE;
  m_doublePtr = &value_ref;
  m_size = definedLength;
  m_stride = strideLength;
  m_isModifiable = false;
}

KOKKOS_INLINE_FUNCTION
void
DeviceVariable::bind(double& value_ref, int definedLength, int strideLength)
{
  m_type = Variable::DOUBLE;
  m_doublePtr = &value_ref;
  m_size = definedLength;
  m_stride = strideLength;
  m_isModifiable = true;
}

KOKKOS_INLINE_FUNCTION
void
DeviceVariable::bind(const int& value_ref, int definedLength, int strideLength)
{
  m_type = Variable::INTEGER;
  m_intPtr = &value_ref;
  m_size = definedLength;
  m_stride = strideLength;
  m_isModifiable = false;
}

KOKKOS_INLINE_FUNCTION
void
DeviceVariable::bind(int& value_ref, int definedLength, int strideLength)
{
  m_type = Variable::INTEGER;
  m_intPtr = &value_ref;
  m_size = definedLength;
  m_stride = strideLength;
  m_isModifiable = true;
}

KOKKOS_INLINE_FUNCTION
DeviceVariable &
DeviceVariable::operator=(double value)
{
  STK_NGP_ThrowRequireMsg(m_size == 1, "double assignment cannot access vector variable as a scalar.");
  STK_NGP_ThrowRequireMsg(m_isModifiable, "double assignment to a const variable.");

  if (m_type == Variable::INTEGER) {
    *const_cast<int*>(m_intPtr) = static_cast<int>(value);
  }
  else if (m_type == Variable::DOUBLE) {
    *const_cast<double*>(m_doublePtr) = value;
  }
  return *this;
}

}
}

#endif

