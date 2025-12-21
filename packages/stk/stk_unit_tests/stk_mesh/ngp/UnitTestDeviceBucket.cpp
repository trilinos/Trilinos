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

#include "stk_mesh/base/Ngp.hpp"

#ifdef STK_USE_DEVICE_MESH

#include <gtest/gtest.h>
#include "stk_mesh/base/BucketConnectivity.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/NgpFieldParallel.hpp"
#include "stk_mesh/baseImpl/DeviceBucketRepository.hpp"
#include "stk_topology/topology.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_unit_test_utils/DeviceBucketTestUtils.hpp"

namespace {

using stk::unit_test_util::DeviceBucket;

void remove_entity(DeviceBucket* bucket, unsigned ordinal)
{
  stk::ngp::RangePolicy<stk::ngp::ExecSpace> policy(0, 1);
  auto func = KOKKOS_LAMBDA(int) { bucket->remove_entity(ordinal); };
  Kokkos::parallel_for("removing_entity_from_bucket", policy, func);
}

}

TEST(DeviceBucket, CopyToHost)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  unsigned maximumBucketCapacity = 2;
  stk::mesh::MeshBuilder builder(stk::parallel_machine_world());
  builder.set_spatial_dimension(3);
  builder.set_initial_bucket_capacity(maximumBucketCapacity);
  builder.set_maximum_bucket_capacity(maximumBucketCapacity);
  std::unique_ptr<stk::mesh::BulkData> bulk = builder.create();

  stk::io::fill_mesh("generated:1x1x1", *bulk);

  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(*bulk);
  stk::mesh::impl::DeviceBucketRepository<stk::mesh::NgpMeshDefaultMemSpace>& deviceBucketRepo = deviceMesh.get_device_bucket_repository();

  DeviceBucket* node_bucket = deviceBucketRepo.get_bucket(stk::topology::NODE_RANK, 0);
  remove_entity(node_bucket, 1);

  for (stk::mesh::EntityRank rank = stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; ++rank)
  {
    const stk::mesh::BucketVector& hostBuckets = bulk->buckets(rank);
    for (size_t i=0; i < hostBuckets.size(); ++i)
    {
      stk::mesh::Bucket& hostBucket = *(hostBuckets[i]);
      DeviceBucket&       deviceBucket = *(deviceBucketRepo.get_bucket(rank, i));
      deviceBucket.sync_to_host(hostBucket);
      stk::unit_test_util::check_entities(hostBucket, deviceBucket);
      stk::unit_test_util::check_connected_entities(hostBucket, deviceBucket);
    }
  }
}

#endif