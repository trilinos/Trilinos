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

#ifndef STK_MESH_ENTITY_FIELD_DATA_HPP
#define STK_MESH_ENTITY_FIELD_DATA_HPP

#include "Kokkos_Macros.hpp"
#include "stk_util/stk_config.h"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace mesh {

template<typename T>
class EntityFieldData {
public:
  KOKKOS_FUNCTION EntityFieldData(T* dataPtr, unsigned length, unsigned componentStride)
    : m_dataPtr(dataPtr),
      m_dataLength(length),
      m_componentStride(componentStride)
  {
  }

  KOKKOS_DEFAULTED_FUNCTION ~EntityFieldData() = default;

  KOKKOS_FUNCTION unsigned size() const { return m_dataLength; }

  KOKKOS_FUNCTION T& operator*() {
#ifdef STK_ENABLE_GPU
    STK_NGP_ThrowAssertMsg(m_dataLength == 1u,
                           "Invalid scalar access with EntityFieldData::operator*() for an Entity with multi-"
                           "component Field data.");
#else
    STK_ThrowAssertMsg(m_dataLength == 1u,
                       "Invalid scalar access with EntityFieldData::operator*() for an Entity with Field length "
                       << m_dataLength);
#endif

    return *m_dataPtr;
  }

  KOKKOS_FUNCTION T& operator[](unsigned component) {
#ifdef STK_ENABLE_GPU
    STK_NGP_ThrowAssertMsg(component < m_dataLength,
                           "Out-of-bounds access to EntityFieldData::operator[]");
#else
    STK_ThrowAssertMsg(component < m_dataLength,
                       "Out-of-bounds access to EntityFieldData::operator[] with component " << component
                       << " for an Entity with Field length " << m_dataLength);
#endif

    return m_dataPtr[component*m_componentStride];
  }

  KOKKOS_FUNCTION const std::remove_const_t<T>& operator[](unsigned component) const {
#ifdef STK_ENABLE_GPU
    STK_NGP_ThrowAssertMsg(component < m_dataLength,
                           "Out-of-bounds access to EntityFieldData::operator[]");
#else
    STK_ThrowAssertMsg(component < m_dataLength,
                       "Out-of-bounds access to EntityFieldData::operator[] with component " << component
                       << " for an Entity with Field length " << m_dataLength);
#endif

    return m_dataPtr[component*m_componentStride];
  }

private:
  T* m_dataPtr;
  unsigned m_dataLength;
  unsigned m_componentStride;
};

}
}

#endif
