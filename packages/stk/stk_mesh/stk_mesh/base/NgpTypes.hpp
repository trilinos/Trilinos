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

#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/ViewVector.hpp"
#include "Kokkos_Core.hpp"

namespace stk {
namespace mesh {

template <typename MemSpace> using DeviceCommMapIndices           = Kokkos::View<FastMeshIndex*, MemSpace>;
template <typename MemSpace> using HostCommMapIndices             = typename DeviceCommMapIndices<MemSpace>::host_mirror_type;
template <typename MemSpace> using NgpCommMapIndices              = Kokkos::View<FastMeshIndex*, MemSpace>;
template <typename MemSpace> using NgpCommMapIndicesHostMirror    = typename NgpCommMapIndices<MemSpace>::host_mirror_type;

template <typename MemSpace> using EntityKeyViewType              = Kokkos::View<EntityKey*, MemSpace>;
template <typename MemSpace> using EntityViewType                 = Kokkos::View<Entity*, MemSpace>;
using HostEntityViewType = Kokkos::View<const Entity*, stk::ngp::HostExecSpace::memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
template <typename MemSpace> using BucketConnectivityType         = Kokkos::View<Entity*, MemSpace>;

template <typename MemSpace> using UnsignedViewType               = Kokkos::View<unsigned*, MemSpace>;
template <typename MemSpace> using Unsigned2dViewType             = Kokkos::View<unsigned**, MemSpace>;
template <typename MemSpace> using OrdinalViewType                = Kokkos::View<ConnectivityOrdinal*, MemSpace>;
template <typename MemSpace> using PartOrdinalViewType            = Kokkos::View<PartOrdinal*, MemSpace>;
using HostPartOrdinalViewType = Kokkos::View<const PartOrdinal*, stk::ngp::HostExecSpace::memory_space,
                                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
template <typename MemSpace> using PermutationViewType            = Kokkos::View<Permutation*, MemSpace>;
template <typename MemSpace> using FastSharedCommMapViewType      = DeviceCommMapIndices<MemSpace>;
template <typename MemSpace> using MeshIndexType                  = Kokkos::View<FastMeshIndex*, MemSpace,
                                                                                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
template <typename MemSpace> using HostMeshIndexType              = typename MeshIndexType<MemSpace>::host_mirror_type;

template <typename MemSpace> using UnsignedPairViewType           = Kokkos::View<Kokkos::pair<unsigned, unsigned>*, MemSpace>;
template <typename MemSpace> using EntityRankViewType             = Kokkos::View<EntityRank*, MemSpace>;

using DeviceStringType = Kokkos::View<char*, stk::ngp::HostPinnedSpace>;
using HostStringType = Kokkos::View<char*, stk::ngp::HostMemSpace>;

template <typename MemSpace> using DeviceFieldMetaDataArrayType   = Kokkos::View<DeviceFieldMetaData*, MemSpace>;
template <typename MemSpace> using HostFieldMetaDataArrayType     = typename DeviceFieldMetaDataArrayType<MemSpace>::host_mirror_type;
template <typename MemSpace> using DeviceBucketsModifiedCollectionType = Kokkos::View<int**, Kokkos::LayoutRight, MemSpace>;

using FieldMetaDataArrayType = impl::ViewVector<FieldMetaData, stk::ngp::HostMemSpace>;

using FieldMetaDataModCountType = Kokkos::View<unsigned, stk::ngp::HostPinnedSpace>;

template <typename Space>
struct DefaultLayoutSelector {
  static constexpr Layout layout = DefaultDeviceLayout;
};

template <>
struct DefaultLayoutSelector<stk::ngp::HostSpace> {
  static constexpr Layout layout = DefaultHostLayout;
};

}
}

#endif // NGPTYPES_HPP
