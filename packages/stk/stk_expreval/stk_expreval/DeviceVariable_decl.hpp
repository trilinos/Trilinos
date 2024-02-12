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

#ifndef stk_expreval_DeviceVariable_decl_hpp
#define stk_expreval_DeviceVariable_decl_hpp

#include <Kokkos_Core.hpp>
#include <stk_expreval/Variable.hpp>

namespace stk {
namespace expreval {

class DeviceVariable
{
public:
  KOKKOS_FUNCTION
  DeviceVariable();

  KOKKOS_FUNCTION
  DeviceVariable(const Variable::Type variableType, int variableSize, int variableStride=1);

  KOKKOS_DEFAULTED_FUNCTION
  DeviceVariable(const DeviceVariable& deviceVariable) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~DeviceVariable() = default;

  KOKKOS_FUNCTION
  DeviceVariable& operator=(const DeviceVariable& deviceVariable);

  KOKKOS_FUNCTION
  double getArrayValue(int index, Variable::ArrayOffset arrayOffsetType) const;

  KOKKOS_FUNCTION
  void assignArrayValue(int index, Variable::ArrayOffset arrayOffsetType, double value) const;

  KOKKOS_FUNCTION
  double getValue() const;

  KOKKOS_FUNCTION
  void bind(const double& value_ref, int definedLength, int strideLength);

  KOKKOS_FUNCTION
  void bind(double& value_ref, int definedLength, int strideLength);

  KOKKOS_FUNCTION
  void bind(const int& value_ref, int definedLength, int strideLength);

  KOKKOS_FUNCTION
  void bind(int& value_ref, int definedLength, int strideLength);

  KOKKOS_FUNCTION
  DeviceVariable& operator=(double value);

private:
  Variable::Type m_type;
  int m_size;
  int m_stride;

  union {
    const double * m_doublePtr;
    const int * m_intPtr;
  };

  union {
    double m_doubleValue;
    int m_intValue;
  };

  bool m_isModifiable;
};

}
}

#endif // stk_expreval_DeviceVariable_decl_hpp

