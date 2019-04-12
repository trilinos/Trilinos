// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_STATICVIEW_HPP
#define TPETRA_DETAILS_STATICVIEW_HPP

#include "TpetraCore_config.h"
//#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_DualView.hpp"

namespace Tpetra {
namespace Details {
namespace Impl {

template<class MemorySpace>
class StaticKokkosAllocation {
public:
  StaticKokkosAllocation () = delete;
  ~StaticKokkosAllocation () = delete;
  StaticKokkosAllocation (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation& operator= (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation (StaticKokkosAllocation&&) = delete;
  StaticKokkosAllocation& operator= (StaticKokkosAllocation&&) = delete;

  // Allocation automatically registers deallocation to happen at
  // Kokkos::finalize.  Reallocation only happens if needed.
  static void* resize (MemorySpace space, const size_t size);
};

#ifdef KOKKOS_ENABLE_CUDA
template<>
class StaticKokkosAllocation<Kokkos::CudaSpace> {
public:
  StaticKokkosAllocation () = delete;
  ~StaticKokkosAllocation () = delete;
  StaticKokkosAllocation (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation& operator= (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation (StaticKokkosAllocation&&) = delete;
  StaticKokkosAllocation& operator= (StaticKokkosAllocation&&) = delete;

  static void* resize (Kokkos::CudaSpace space, const size_t size);
};

template<>
class StaticKokkosAllocation<Kokkos::CudaUVMSpace> {
public:
  StaticKokkosAllocation () = delete;
  ~StaticKokkosAllocation () = delete;
  StaticKokkosAllocation (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation& operator= (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation (StaticKokkosAllocation&&) = delete;
  StaticKokkosAllocation& operator= (StaticKokkosAllocation&&) = delete;

  static void* resize (Kokkos::CudaUVMSpace space, const size_t size);
};

template<>
class StaticKokkosAllocation<Kokkos::CudaHostPinnedSpace> {
public:
  StaticKokkosAllocation () = delete;
  ~StaticKokkosAllocation () = delete;
  StaticKokkosAllocation (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation& operator= (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation (StaticKokkosAllocation&&) = delete;
  StaticKokkosAllocation& operator= (StaticKokkosAllocation&&) = delete;

  static void* resize (Kokkos::CudaHostPinnedSpace space, const size_t size);
};
#endif // KOKKOS_ENABLE_CUDA

template<>
class StaticKokkosAllocation<Kokkos::HostSpace> {
public:
  StaticKokkosAllocation () = delete;
  ~StaticKokkosAllocation () = delete;
  StaticKokkosAllocation (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation& operator= (const StaticKokkosAllocation&) = delete;
  StaticKokkosAllocation (StaticKokkosAllocation&&) = delete;
  StaticKokkosAllocation& operator= (StaticKokkosAllocation&&) = delete;

  static void* resize (Kokkos::HostSpace space, const size_t size);
};

template<class ValueType, class MemorySpace>
ValueType*
getStaticKokkosMemory (MemorySpace space,
                       const size_t num_entries,
                       const size_t value_size = sizeof (ValueType))
{
  void* ptr = StaticKokkosAllocation<MemorySpace>::resize
    (space, num_entries * value_size);
  return reinterpret_cast<ValueType*> (ptr);
}

} // namespace Impl

template<class ValueType, class DeviceType>
Kokkos::View<ValueType*, DeviceType>
getStatic1dView (const size_t size)
{
  using Impl::getStaticKokkosMemory;
  using mem_space = typename DeviceType::memory_space;
  using view_type = Kokkos::View<ValueType*, DeviceType>;

  ValueType* ptr = getStaticKokkosMemory<ValueType> (mem_space (), size);
  return view_type (ptr, size);
}

template<class ValueType, class DeviceType>
Kokkos::View<ValueType**, Kokkos::LayoutLeft, DeviceType>
getStatic2dView (const size_t num_rows, const size_t num_cols)
{
  using Impl::getStaticKokkosMemory;
  using mem_space = typename DeviceType::memory_space;
  using view_type = Kokkos::View<ValueType**, Kokkos::LayoutLeft, DeviceType>;

  const size_t size = num_rows * num_cols;
  ValueType* ptr = getStaticKokkosMemory<ValueType> (mem_space (), size);
  return view_type (ptr, num_rows, num_cols);
}

template<class ValueType, class DeviceType>
Kokkos::DualView<ValueType**, Kokkos::LayoutLeft, DeviceType>
getStatic2dDualView (const size_t num_rows, const size_t num_cols)
{
  using dual_view_type =
    Kokkos::DualView<ValueType**, Kokkos::LayoutLeft, DeviceType>;
  using d_view_type = typename dual_view_type::t_dev;
  using h_view_type = typename dual_view_type::t_host;

  auto d_view = getStatic2dView<ValueType, DeviceType> (num_rows, num_cols);
  // Preserve the invariant that Kokkos::create_mirror_view returns
  // the input View here, if and only if it would have returned it if
  // the allocating View constructor were called.
  h_view_type h_view;
  if (std::is_same<typename d_view_type::memory_space,
                   typename h_view_type::memory_space>::value) {
    h_view = Kokkos::create_mirror_view (d_view);
  }
  else {
    h_view = getStatic2dView<ValueType,
      typename h_view_type::device_type> (num_rows, num_cols);
  }

  return dual_view_type (d_view, h_view);
}


} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_STATICVIEW_HPP
