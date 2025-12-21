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

#include <gtest/gtest.h>
#include "stk_mesh/base/BucketConnectivity.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/NgpFieldParallel.hpp"
#include "stk_mesh/baseImpl/DeviceBucketRepository.hpp"
#include "stk_topology/topology.hpp"
#include "stk_io/FillMesh.hpp"

namespace stk {
namespace unit_test_util {

using Memspace = stk::mesh::NgpMeshDefaultMemSpace;
using DevicePartition = stk::mesh::impl::DevicePartition<Memspace>;
using DeviceBucket = stk::mesh::DeviceBucketT<Memspace>;
using DeviceBucketWrapper = stk::mesh::impl::DeviceBucketPtrWrapper<Memspace>;

Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace>
get_entities(const DeviceBucket deviceBucket);

void check_entities(const stk::mesh::Bucket& hostBucket, const DeviceBucket& deviceBucket);


Kokkos::View<stk::mesh::Entity*, stk::ngp::HostMemSpace>
get_connected_entities(const DeviceBucket& deviceBucket, stk::mesh::EntityRank rank,
                       unsigned entityOrdinal);

Kokkos::View<stk::mesh::Permutation*, stk::ngp::HostMemSpace>
get_connected_permutations(const DeviceBucket deviceBucket, stk::mesh::EntityRank rank,
                           unsigned entityOrdinal);

void check_connected_entities(const stk::mesh::Bucket& hostBucket, const DeviceBucket& deviceBucket);

void check_repo_counts(stk::mesh::impl::BucketRepository& hostBucketRepo,
                       stk::mesh::impl::DeviceBucketRepository<stk::mesh::NgpMeshDefaultMemSpace> deviceBucketRepo);

void check_partition_attributes(const stk::mesh::impl::Partition& hostPartition, const DevicePartition& devicePartition);

void check_partition_ordinals(const stk::mesh::impl::Partition& hostPartition, const DevicePartition& devicePartition);

void check_partition_bucket_attributes(const stk::mesh::impl::Partition& hostPartition, const DevicePartition& devicePartition);

}
}