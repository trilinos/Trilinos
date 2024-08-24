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

#ifndef NGPTYPES_HPP
#define NGPTYPES_HPP

#include <stk_util/ngp/NgpSpaces.hpp>
#include "stk_mesh/base/Types.hpp"
#include <Kokkos_Core.hpp>

namespace stk {
namespace mesh {

using DeviceCommMapIndices        = Kokkos::View<FastMeshIndex*, stk::ngp::MemSpace>;
using HostCommMapIndices          = DeviceCommMapIndices::HostMirror;
using EntityKeyViewType           = Kokkos::View<EntityKey*, stk::ngp::MemSpace>;
using EntityViewType              = Kokkos::View<Entity*, stk::ngp::MemSpace>;
using HostEntityViewType          = Kokkos::View<const Entity*, stk::ngp::HostExecSpace::memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using BucketConnectivityType      = Kokkos::View<Entity*, stk::ngp::MemSpace>;
using UnsignedViewType            = Kokkos::View<unsigned*, stk::ngp::MemSpace>;
using Unsigned2dViewType          = Kokkos::View<unsigned**, stk::ngp::MemSpace>;
using BoolViewType                = Kokkos::View<bool*, stk::ngp::MemSpace>;
using OrdinalViewType             = Kokkos::View<ConnectivityOrdinal*, stk::ngp::MemSpace>;
using PartOrdinalViewType         = Kokkos::View<PartOrdinal*, stk::ngp::MemSpace>;
using HostPartOrdinalViewType     = Kokkos::View<const PartOrdinal*, stk::ngp::HostExecSpace::memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using PermutationViewType         = Kokkos::View<Permutation*, stk::ngp::MemSpace>;
using FastSharedCommMapViewType   = DeviceCommMapIndices;
using HostMeshIndexType           = Kokkos::View<FastMeshIndex*>::HostMirror;
using MeshIndexType               = Kokkos::View<const FastMeshIndex*, stk::ngp::MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
using BucketEntityOffsetsViewType = Kokkos::View<int*, stk::ngp::MemSpace>;

template <typename T> using FieldDataDeviceViewType = Kokkos::View<T***, Kokkos::LayoutRight, stk::ngp::MemSpace>;
template <typename T> using FieldDataHostViewType   = Kokkos::View<T***, Kokkos::LayoutRight, stk::ngp::HostPinnedSpace>;

using FieldDataPointerHostViewType = Kokkos::View<uintptr_t*, Kokkos::LayoutRight, stk::ngp::HostPinnedSpace>;
using FieldDataPointerDeviceViewType = Kokkos::View<uintptr_t*, Kokkos::LayoutRight, stk::ngp::MemSpace>;

template <typename T> using UnmanagedHostInnerView = Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
template <typename T> using UnmanagedDevInnerView = Kokkos::View<T**, Kokkos::LayoutRight, stk::ngp::MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

#ifdef STK_USE_DEVICE_MESH
#define ORDER_INDICES(i,j) j,i
#else
#define ORDER_INDICES(i,j) i,j
#endif

}
}

#endif // NGPTYPES_HPP
