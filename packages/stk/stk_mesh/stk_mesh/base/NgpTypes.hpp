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

#include <stk_mesh/base/NgpSpaces.hpp>
#include "stk_mesh/base/Types.hpp"
#include <Kokkos_Core.hpp>

namespace stk {
namespace mesh {

using DeviceCommMapIndices      = Kokkos::View<FastMeshIndex*, MemSpace>;
using EntityKeyViewType         = Kokkos::View<EntityKey*, MemSpace>;
using EntityViewType            = Kokkos::View<Entity*, MemSpace>;
using BucketConnectivityType    = Kokkos::View<Entity**, MemSpace>;
using UnsignedViewType          = Kokkos::View<unsigned*, MemSpace>;
using BoolViewType              = Kokkos::View<bool*, MemSpace>;
using OrdinalViewType           = Kokkos::View<ConnectivityOrdinal*, MemSpace>;
using PartOrdinalViewType       = Kokkos::View<PartOrdinal*, MemSpace>;
using PermutationViewType       = Kokkos::View<Permutation*, MemSpace>;
using FastSharedCommMapViewType = Kokkos::View<FastMeshIndex*, MemSpace>;
using HostMeshIndexType         = Kokkos::View<FastMeshIndex*>::HostMirror;
using MeshIndexType             = Kokkos::View<const FastMeshIndex*, MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

template <typename T> using FieldDataDeviceViewType = Kokkos::View<T***, Kokkos::LayoutRight, MemSpace>;
template <typename T> using FieldDataHostViewType   = Kokkos::View<T***, Kokkos::LayoutRight, HostPinnedSpace>;

template <typename T> using UnmanagedHostInnerView = Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
template <typename T> using UnmanagedDevInnerView = Kokkos::View<T**, Kokkos::LayoutRight, stk::mesh::MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

#ifdef KOKKOS_ENABLE_CUDA
#define ORDER_INDICES(i,j) j,i
#else
#define ORDER_INDICES(i,j) i,j
#endif

}
}

#endif // NGPTYPES_HPP
