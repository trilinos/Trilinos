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

#ifndef STK_MESH_IMPL_DEVICE_MESH_VIEW_VECTOR_HPP
#define STK_MESH_IMPL_DEVICE_MESH_VIEW_VECTOR_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/ViewVector.hpp>
#include "Kokkos_Core.hpp"
#include "Kokkos_Macros.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

namespace stk {
namespace mesh {
namespace impl {

template <typename T, typename MemSpace = stk::ngp::HostMemSpace, typename SizeT = unsigned>
class ImplDeviceMeshViewVector : public ViewVector<T, MemSpace, SizeT> {
 private:
  using base_t = ViewVector<T, MemSpace, SizeT>;
  using active_entries_size_view = Kokkos::View<SizeT, stk::ngp::UVMMemSpace>;

 public:
   ImplDeviceMeshViewVector(const std::string& name)
    : base_t(name),
      m_numActiveEntries("numActiveEntries")
  {
    clear_active_entries();
  }

  ImplDeviceMeshViewVector(const std::string& name, SizeT size)
    : base_t(name, size),
      m_numActiveEntries("numActiveEntries")
  {
    clear_active_entries();
  }

  KOKKOS_DEFAULTED_FUNCTION ImplDeviceMeshViewVector() = default;
  KOKKOS_DEFAULTED_FUNCTION ImplDeviceMeshViewVector(const ImplDeviceMeshViewVector& other) = default;
  KOKKOS_DEFAULTED_FUNCTION ImplDeviceMeshViewVector(ImplDeviceMeshViewVector&& other) = default;
  KOKKOS_DEFAULTED_FUNCTION ImplDeviceMeshViewVector& operator=(const ImplDeviceMeshViewVector& rhs) = default;
  KOKKOS_DEFAULTED_FUNCTION ImplDeviceMeshViewVector& operator=(ImplDeviceMeshViewVector&& rhs) = default;
  KOKKOS_DEFAULTED_FUNCTION ~ImplDeviceMeshViewVector() = default;

  void clear() {
    base_t::clear();
    m_numActiveEntries() = 0;
  }

  void push_back(const T& value) {
    base_t::push_back(value);
    increment_num_active_entries();
  }

  void push_back(T&& value) {
    base_t::push_back(value);
    increment_num_active_entries();
  }

  template <typename... Args>
  void emplace_back(Args&&... args) {
    base_t::emplace_back(args...);
    increment_num_active_entries();
  }

  KOKKOS_INLINE_FUNCTION SizeT num_active_entries() const { return m_numActiveEntries(); }
  KOKKOS_INLINE_FUNCTION void set_active_entries(SizeT size) { m_numActiveEntries() = size; }
  KOKKOS_INLINE_FUNCTION void clear_active_entries() { m_numActiveEntries() = 0; }
  KOKKOS_INLINE_FUNCTION void increment_num_active_entries() { ++m_numActiveEntries(); }
  KOKKOS_INLINE_FUNCTION void decrement_num_active_entries() { --m_numActiveEntries(); }

 private:
  active_entries_size_view m_numActiveEntries;
};

} } }

#endif
