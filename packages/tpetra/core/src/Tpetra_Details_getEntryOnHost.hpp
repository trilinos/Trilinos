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

#ifndef TPETRA_DETAILS_GETENTRYONHOST_HPP
#define TPETRA_DETAILS_GETENTRYONHOST_HPP

/// \file Tpetra_Details_getEntryOnHost.hpp
/// \brief Declaration and definition of
///   Tpetra::Details::getEntryOnHost.
/// \warning The contents of this file are implementation details of
///   Tpetra.  We make no promises of backwards compatibility.

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {
namespace Impl {

template<class ViewType,
         class IndexType = typename ViewType::size_type,
         const bool isHostSpace =
           std::is_same<typename ViewType::memory_space,
                        Kokkos::HostSpace>::value>
struct GetEntryOnHost {
  static typename ViewType::non_const_value_type
  getEntryOnHost (const ViewType& x,
                  const IndexType ind);
};

template<class ViewType,
         class IndexType>
struct GetEntryOnHost<ViewType, IndexType, true> {
  static typename ViewType::non_const_value_type
  getEntryOnHost (const ViewType& x,
                  const IndexType ind)
  {
    static_assert (ViewType::Rank == 1, "x must be a rank-1 Kokkos::View.");
    return x(ind);
  }
};

template<class ViewType,
         class IndexType>
struct GetEntryOnHost<ViewType, IndexType, false> {
  static typename ViewType::non_const_value_type
  getEntryOnHost (const ViewType& x,
                  const IndexType ind)
  {
    // Don't assume UVM.  Carefully get a 0-D subview of the entry of
    // the array, and copy to device.  Do not use host mirror, because
    // that could just be a UVM View if using UVM.
    static_assert (ViewType::Rank == 1, "x must be a rank-1 Kokkos::View.");
#ifdef KOKKOS_ENABLE_CUDA
    // Do not use Kokkos::create_mirror_view, because that could just
    // be a UVM View if using UVM.
    typedef typename ViewType::device_type device_type;
    // Hide this in ifdef, to avoid unused typedef warning.
    typedef typename device_type::execution_space dev_exec_space;
#endif // KOKKOS_ENABLE_CUDA
    typedef typename ViewType::HostMirror::execution_space host_exec_space;
    typedef Kokkos::Device<host_exec_space, Kokkos::HostSpace> host_device_type;
    typedef typename ViewType::non_const_value_type value_type;

    value_type val;
    Kokkos::View<value_type, host_device_type> view_h (&val);
    auto view_d = Kokkos::subview (x, ind); // 0-D View
#ifdef KOKKOS_ENABLE_CUDA
    typedef typename device_type::memory_space dev_memory_space;
    if (std::is_same<dev_memory_space, Kokkos::CudaUVMSpace>::value) {
      dev_exec_space::fence (); // for UVM's sake.
    }
#endif // KOKKOS_ENABLE_CUDA
    Kokkos::deep_copy (view_h, view_d);
#ifdef KOKKOS_ENABLE_CUDA
    if (std::is_same<dev_memory_space, Kokkos::CudaUVMSpace>::value) {
      dev_exec_space::fence (); // for UVM's sake.
    }
#endif // KOKKOS_ENABLE_CUDA
    return val;
  }
};

} // namespace Impl

template<class ViewType,
         class IndexType>
typename ViewType::non_const_value_type
getEntryOnHost (const ViewType& x,
                const IndexType ind)
{
  return Impl::GetEntryOnHost<ViewType, IndexType>::getEntryOnHost (x, ind);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GETENTRYONHOST_HPP
