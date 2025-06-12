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
#include <stk_util/stk_config.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

namespace {

class FieldDataAccess : public stk::unit_test_util::MeshFixture
{
public:
  FieldDataAccess()
    : batchTimer(get_comm())
  { }

protected:
  void declare_centroid_field()
  {
    m_centroidField = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "centroid");
    stk::mesh::put_field_on_mesh(*m_centroidField, get_meta().universal_part(), 3, nullptr);
  }

  void declare_centroid_partial_mesh(unsigned numBlocks)
  {
    m_centroidField = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "centroid");
    stk::mesh::PartVector elemBlocks;
    stk::mesh::fill_element_block_parts(get_meta(), stk::topology::HEX_8, elemBlocks);
    for(unsigned i = 0; i < numBlocks; ++i) {
      const stk::mesh::Part& part = *elemBlocks[i];
      stk::mesh::put_field_on_mesh(*m_centroidField, part, 3, nullptr);
    }
  }

  void fill_multi_block_mesh(unsigned numElemsPerDim)
  {
    stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(numElemsPerDim), get_bulk());
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemsPerDim);
  }

  template <typename CENTROID_FUNCTOR, typename VERIFIER_FUNCTOR>
  void run_single_block_test(int numIters, const CENTROID_FUNCTOR& centroidFunctor,
                             const VERIFIER_FUNCTOR& verifierFunctor)
  {
    const unsigned NUM_RUNS = 5;
    const int ELEMS_PER_DIM = 100;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    declare_centroid_field();
    stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(ELEMS_PER_DIM), get_bulk());

    stk::mesh::Selector selector(get_meta().locally_owned_part());
    const stk::mesh::Field<double>& coordsField = *dynamic_cast<const stk::mesh::Field<double>*>(get_meta().coordinate_field());

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        centroidFunctor(selector, *m_centroidField, coordsField);
      }
      batchTimer.stop_batch_timer();
    }

    verifierFunctor(ELEMS_PER_DIM, 1, 1, *m_centroidField);

    batchTimer.print_batch_timing(numIters);
  }

  template <typename CENTROID_FUNCTOR, typename VERIFIER_FUNCTOR>
  void run_multiple_block_test(int numIters, const CENTROID_FUNCTOR& centroidFunctor,
                               const VERIFIER_FUNCTOR& verifierFunctor)
  {
    const unsigned NUM_RUNS = 5;
    const int ELEMS_PER_DIM = 100;
    const int NUM_BLOCKS = 100;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::performance_tests::setup_multiple_blocks(get_meta(), NUM_BLOCKS);
    declare_centroid_field();
    fill_multi_block_mesh(ELEMS_PER_DIM);

    stk::mesh::Selector selector(get_meta().locally_owned_part());
    const stk::mesh::Field<double>& coordsField = *dynamic_cast<const stk::mesh::Field<double>*>(get_meta().coordinate_field());

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        centroidFunctor(selector, *m_centroidField, coordsField);
      }
      batchTimer.stop_batch_timer();
    }

    verifierFunctor(ELEMS_PER_DIM, NUM_BLOCKS, NUM_BLOCKS, *m_centroidField);

    batchTimer.print_batch_timing(numIters);
  }

  template <typename CENTROID_FUNCTOR, typename VERIFIER_FUNCTOR>
  void run_partial_block_test(int numIters, const CENTROID_FUNCTOR& centroidFunctor,
                              const VERIFIER_FUNCTOR& verifierFunctor)
  {
    const unsigned NUM_RUNS = 5;
    const int ELEMS_PER_DIM = 100;
    const int NUM_BLOCKS = 100;
    const int USED_BLOCKS = 50;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::performance_tests::setup_multiple_blocks(get_meta(), NUM_BLOCKS);
    declare_centroid_partial_mesh(USED_BLOCKS);
    fill_multi_block_mesh(ELEMS_PER_DIM);

    stk::mesh::Selector selector(get_meta().locally_owned_part());
    const stk::mesh::Field<double>& coordsField = *dynamic_cast<const stk::mesh::Field<double>*>(get_meta().coordinate_field());

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        centroidFunctor(selector, *m_centroidField, coordsField);
      }
      batchTimer.stop_batch_timer();
    }

    verifierFunctor(ELEMS_PER_DIM, NUM_BLOCKS, USED_BLOCKS, *m_centroidField);

    batchTimer.print_batch_timing(numIters);
  }

  stk::unit_test_util::BatchTimer batchTimer;
  stk::mesh::Field<double> *m_centroidField;
};

class LegacyFieldDataAccess : public FieldDataAccess {};

auto host_verify_averaged_centroids_are_center_of_mesh = [](int elemsPerDim, int numTotalBlocks, int numUsedBlocks,
                                                            stk::mesh::Field<double>& centroidField)
{
  stk::mesh::BulkData& bulk = centroidField.get_mesh();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  std::vector<double> average = stk::performance_tests::get_centroid_average_from_host(bulk, centroidField,
                                                                                       stk::mesh::Selector(meta.universal_part()));
  double meshCenterX = elemsPerDim * ((double)numUsedBlocks/numTotalBlocks) / 2.0;
  double meshCenterY = elemsPerDim / 2.0;
  double meshCenterZ = elemsPerDim / 2.0;

  EXPECT_DOUBLE_EQ(meshCenterX, average[0]);
  EXPECT_DOUBLE_EQ(meshCenterY, average[1]);
  EXPECT_DOUBLE_EQ(meshCenterZ, average[2]);
};

auto device_verify_averaged_centroids_are_center_of_mesh = [](int elemsPerDim, int numTotalBlocks, int numUsedBlocks,
                                                              stk::mesh::Field<double>& centroidField)
{
  stk::mesh::BulkData& bulk = centroidField.get_mesh();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  std::vector<double> average = stk::performance_tests::get_centroid_average_from_device(bulk, centroidField,
                                                                                         stk::mesh::Selector(meta.universal_part()));
  double meshCenterX = elemsPerDim * ((double)numUsedBlocks/numTotalBlocks) / 2.0;
  double meshCenterY = elemsPerDim / 2.0;
  double meshCenterZ = elemsPerDim / 2.0;

  EXPECT_DOUBLE_EQ(meshCenterX, average[0]);
  EXPECT_DOUBLE_EQ(meshCenterY, average[1]);
  EXPECT_DOUBLE_EQ(meshCenterZ, average[2]);
};


auto legacy_compute_centroid_entity_access = [](const stk::mesh::Selector& selector,
                                                stk::mesh::Field<double>& centroidField,
                                                const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  for (const stk::mesh::Bucket* bucket : elemBuckets) {
    const unsigned numNodes = bucket->topology().num_nodes();

    for (stk::mesh::Entity elem : *bucket) {
      double* centroid = stk::mesh::field_data(centroidField, elem);

      const unsigned numComponents = stk::mesh::field_scalars_per_entity(centroidField, elem);
      if (numComponents == 0u) {
        continue;
      }

      centroid[0] = 0.0;
      centroid[1] = 0.0;
      centroid[2] = 0.0;

      const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
      for (unsigned n = 0; n < numNodes; ++n) {
        const double* nodeCoords = stk::mesh::field_data(coordsField, nodes[n]);
        centroid[0] += nodeCoords[0];
        centroid[1] += nodeCoords[1];
        centroid[2] += nodeCoords[2];
      }

      centroid[0] /= numNodes;
      centroid[1] /= numNodes;
      centroid[2] /= numNodes;
    }
  }
};

auto compute_centroid_entity_access = [](const stk::mesh::Selector& selector,
                                         stk::mesh::Field<double>& centroidField,
                                         const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite>();
  auto coordsData = coordsField.data<stk::mesh::ReadOnly>();

  for (const stk::mesh::Bucket* bucket : elemBuckets) {
    const unsigned numNodes = bucket->topology().num_nodes();

    for (stk::mesh::Entity elem : *bucket) {
      auto centroid = centroidData.entity_values(elem);

      if (not centroid.is_field_defined()) {
        continue;
      }

      centroid(0_comp) = 0.0;
      centroid(1_comp) = 0.0;
      centroid(2_comp) = 0.0;

      const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
      for (unsigned n = 0; n < numNodes; ++n) {
        auto nodeCoords = coordsData.entity_values(nodes[n]);
        centroid(0_comp) += nodeCoords(0_comp);
        centroid(1_comp) += nodeCoords(1_comp);
        centroid(2_comp) += nodeCoords(2_comp);
      }

      centroid(0_comp) /= numNodes;
      centroid(1_comp) /= numNodes;
      centroid(2_comp) /= numNodes;
    }
  }
};

auto legacy_compute_centroid_bucket_access = [](const stk::mesh::Selector& selector,
                                                stk::mesh::Field<double>& centroidField,
                                                const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  for (const stk::mesh::Bucket* bucket : elemBuckets) {
    double* centroid = stk::mesh::field_data(centroidField, *bucket);

    const unsigned numComponents = stk::mesh::field_scalars_per_entity(centroidField, *bucket);
    if (numComponents == 0u) {
      continue;
    }

    const unsigned numNodes = bucket->topology().num_nodes();
    for (unsigned elem = 0; elem < bucket->size(); ++elem) {
      centroid[elem*3  ] = 0.0;
      centroid[elem*3+1] = 0.0;
      centroid[elem*3+2] = 0.0;

      const stk::mesh::Entity* nodes = bucket->begin_nodes(elem);
      for (unsigned n = 0; n < numNodes; ++n) {
        const double* nodeCoords = stk::mesh::field_data(coordsField, nodes[n]);
        centroid[elem*3  ] += nodeCoords[0];
        centroid[elem*3+1] += nodeCoords[1];
        centroid[elem*3+2] += nodeCoords[2];
      }

      centroid[elem*3  ] /= numNodes;
      centroid[elem*3+1] /= numNodes;
      centroid[elem*3+2] /= numNodes;
    }
  }
};

auto compute_centroid_bucket_access = [](const stk::mesh::Selector& selector,
                                         stk::mesh::Field<double>& centroidField,
                                         const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite>();
  auto coordsData = coordsField.data<stk::mesh::ReadOnly>();

  for (const stk::mesh::Bucket* bucket : elemBuckets) {
    auto centroid = centroidData.bucket_values(*bucket);

    if (not centroid.is_field_defined()) {
      continue;
    }

    const int numNodes = bucket->topology().num_nodes();
    for (stk::mesh::EntityIdx elem : bucket->entities()) {
      centroid(elem, 0_comp) = 0.0;
      centroid(elem, 1_comp) = 0.0;
      centroid(elem, 2_comp) = 0.0;

      const stk::mesh::Entity* nodes = bucket->begin_nodes(elem);
      for (int n = 0; n < numNodes; ++n) {
        auto nodeCoords = coordsData.entity_values(nodes[n]);
        centroid(elem, 0_comp) += nodeCoords(0_comp);
        centroid(elem, 1_comp) += nodeCoords(1_comp);
        centroid(elem, 2_comp) += nodeCoords(2_comp);
      }

      centroid(elem, 0_comp) /= numNodes;
      centroid(elem, 1_comp) /= numNodes;
      centroid(elem, 2_comp) /= numNodes;
    }
  }
};

#ifdef NDEBUG
  constexpr int numHostIters = 100;
#else
  constexpr int numHostIters = 1;
#endif


//------------------------------------------------------------------------------
TEST_F(LegacyFieldDataAccess, entity_SingleBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test(numHostIters, legacy_compute_centroid_entity_access,
                        host_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(LegacyFieldDataAccess, entity_MultiBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_multiple_block_test(numHostIters, legacy_compute_centroid_entity_access,
                          host_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(LegacyFieldDataAccess, entity_PartialBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_partial_block_test(numHostIters, legacy_compute_centroid_entity_access,
                         host_verify_averaged_centroids_are_center_of_mesh);
}


TEST_F(FieldDataAccess, entity_SingleBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test(numHostIters, compute_centroid_entity_access,
                        host_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(FieldDataAccess, entity_MultiBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_multiple_block_test(numHostIters, compute_centroid_entity_access,
                          host_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(FieldDataAccess, entity_PartialBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_partial_block_test(numHostIters, compute_centroid_entity_access,
                         host_verify_averaged_centroids_are_center_of_mesh);
}

//------------------------------------------------------------------------------
TEST_F(LegacyFieldDataAccess, bucket_SingleBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test(numHostIters, legacy_compute_centroid_bucket_access,
                        host_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(LegacyFieldDataAccess, bucket_MultiBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_multiple_block_test(numHostIters, legacy_compute_centroid_bucket_access,
                          host_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(LegacyFieldDataAccess, bucket_PartialBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_partial_block_test(numHostIters, legacy_compute_centroid_bucket_access,
                         host_verify_averaged_centroids_are_center_of_mesh);
}


TEST_F(FieldDataAccess, bucket_SingleBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test(numHostIters, compute_centroid_bucket_access,
                        host_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(FieldDataAccess, bucket_MultiBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_multiple_block_test(numHostIters, compute_centroid_bucket_access,
                          host_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(FieldDataAccess, bucket_PartialBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_partial_block_test(numHostIters, compute_centroid_bucket_access,
                         host_verify_averaged_centroids_are_center_of_mesh);
}



class LegacyDeviceFieldDataAccess : public FieldDataAccess {};
class DeviceFieldDataAccess : public FieldDataAccess {};

//------------------------------------------------------------------------------
// Can't put a device lambda in another lambda, so add a lambda->function->device_lambda layer
//
void legacy_device_compute_centroid_entity_component_access_function(const stk::mesh::Selector& selector,
                                                                     stk::mesh::Field<double>& centroidField,
                                                                     const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::mesh::NgpField<double>& ngpCentroidField = stk::mesh::get_updated_ngp_field<double>(centroidField);
  stk::mesh::NgpField<double>& ngpCoordsField = stk::mesh::get_updated_ngp_field<double>(coordsField);

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
      const unsigned numComponents = ngpCentroidField.get_num_components_per_entity(elem);
      if (numComponents == 0) {
        return;
      }

      ngpCentroidField(elem, 0) = 0.0;
      ngpCentroidField(elem, 1) = 0.0;
      ngpCentroidField(elem, 2) = 0.0;

      stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elem);
      for (size_t i = 0; i < nodes.size(); ++i) {
        stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[i]);

        ngpCentroidField(elem, 0) += ngpCoordsField(nodeIndex, 0);
        ngpCentroidField(elem, 1) += ngpCoordsField(nodeIndex, 1);
        ngpCentroidField(elem, 2) += ngpCoordsField(nodeIndex, 2);
      }

      ngpCentroidField(elem, 0) /= nodes.size();
      ngpCentroidField(elem, 1) /= nodes.size();
      ngpCentroidField(elem, 2) /= nodes.size();
    }
  );
};

void legacy_device_compute_centroid_entity_access_function(const stk::mesh::Selector& selector,
                                                           stk::mesh::Field<double>& centroidField,
                                                           const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::mesh::NgpField<double>& ngpCentroidField = stk::mesh::get_updated_ngp_field<double>(centroidField);
  stk::mesh::NgpField<double>& ngpCoordsField = stk::mesh::get_updated_ngp_field<double>(coordsField);

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
      const unsigned numComponents = ngpCentroidField.get_num_components_per_entity(elem);
      if (numComponents == 0) {
        return;
      }

      stk::mesh::EntityFieldData<double> elemCentroid = ngpCentroidField(elem);
      elemCentroid[0] = 0.0;
      elemCentroid[1] = 0.0;
      elemCentroid[2] = 0.0;

      stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elem);
      for (size_t i = 0; i < nodes.size(); ++i) {
        stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[i]);
        const stk::mesh::EntityFieldData<double> nodeCoords = ngpCoordsField(nodeIndex);

        elemCentroid[0] += nodeCoords[0];
        elemCentroid[1] += nodeCoords[1];
        elemCentroid[2] += nodeCoords[2];
      }

      elemCentroid[0] /= nodes.size();
      elemCentroid[1] /= nodes.size();
      elemCentroid[2] /= nodes.size();
    }
  );
};

void device_compute_centroid_entity_access_function(const stk::mesh::Selector& selector,
                                                    stk::mesh::Field<double>& centroidField,
                                                    const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  auto coordsData = coordsField.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
      auto centroidValues = centroidData.entity_values(elem);
      if (not centroidValues.is_field_defined()) {
        return;
      }

      centroidValues(0_comp) = 0.0;
      centroidValues(1_comp) = 0.0;
      centroidValues(2_comp) = 0.0;

      stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elem);
      for (size_t i = 0; i < nodes.size(); ++i) {
        auto coordsValues = coordsData.entity_values(nodes[i]);

        centroidValues(0_comp) += coordsValues(0_comp);
        centroidValues(1_comp) += coordsValues(1_comp);
        centroidValues(2_comp) += coordsValues(2_comp);
      }

      centroidValues(0_comp) /= nodes.size();
      centroidValues(1_comp) /= nodes.size();
      centroidValues(2_comp) /= nodes.size();
    }
  );
};

void device_compute_centroid_bucket_access_function(const stk::mesh::Selector& selector,
                                                    stk::mesh::Field<double>& centroidField,
                                                    const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  auto coordsData = coordsField.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::ELEM_RANK, selector);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto centroidValues = centroidData.bucket_values(bucketId);
      if (not centroidValues.is_field_defined()) {
        return;
      }

      const stk::mesh::EntityIdx numElems = centroidValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numElems),
        [&](stk::mesh::EntityIdx elem) {
          centroidValues(elem, 0_comp) = 0.0;
          centroidValues(elem, 1_comp) = 0.0;
          centroidValues(elem, 2_comp) = 0.0;

          stk::mesh::NgpMesh::ConnectedNodes nodes =
              ngpMesh.get_bucket(stk::topology::ELEM_RANK, bucketId).get_nodes(static_cast<int>(elem));
          for (size_t i = 0; i < nodes.size(); ++i) {
            auto coordsValues = coordsData.entity_values(nodes[i]);

            centroidValues(elem, 0_comp) += coordsValues(0_comp);
            centroidValues(elem, 1_comp) += coordsValues(1_comp);
            centroidValues(elem, 2_comp) += coordsValues(2_comp);
          }

          centroidValues(elem, 0_comp) /= nodes.size();
          centroidValues(elem, 1_comp) /= nodes.size();
          centroidValues(elem, 2_comp) /= nodes.size();
        }
      );
    }
  );
};

auto legacy_device_compute_centroid_entity_component_access = [](const stk::mesh::Selector& selector,
                                                                 stk::mesh::Field<double>& centroidField,
                                                                 const stk::mesh::Field<double>& coordsField)
{
  legacy_device_compute_centroid_entity_component_access_function(selector, centroidField, coordsField);
};


auto legacy_device_compute_centroid_entity_access = [](const stk::mesh::Selector& selector,
                                                       stk::mesh::Field<double>& centroidField,
                                                       const stk::mesh::Field<double>& coordsField)
{
  legacy_device_compute_centroid_entity_access_function(selector, centroidField, coordsField);
};

auto device_compute_centroid_entity_access = [](const stk::mesh::Selector& selector,
                                                stk::mesh::Field<double>& centroidField,
                                                const stk::mesh::Field<double>& coordsField)
{
  device_compute_centroid_entity_access_function(selector, centroidField, coordsField);
};

auto device_compute_centroid_bucket_access = [](const stk::mesh::Selector& selector,
                                                stk::mesh::Field<double>& centroidField,
                                                const stk::mesh::Field<double>& coordsField)
{
  device_compute_centroid_bucket_access_function(selector, centroidField, coordsField);
};

#ifdef STK_ENABLE_GPU
  constexpr int numDeviceIters = 2000;
#else
  constexpr int numDeviceIters = 100;
#endif
//------------------------------------------------------------------------------
TEST_F(LegacyDeviceFieldDataAccess, entity_component_SingleBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test(numDeviceIters, legacy_device_compute_centroid_entity_component_access,
                        device_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(LegacyDeviceFieldDataAccess, entity_component_MultiBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_multiple_block_test(numDeviceIters, legacy_device_compute_centroid_entity_component_access,
                          device_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(LegacyDeviceFieldDataAccess, entity_component_PartialBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_partial_block_test(numDeviceIters, legacy_device_compute_centroid_entity_component_access,
                         device_verify_averaged_centroids_are_center_of_mesh);
}

//------------------------------------------------------------------------------
TEST_F(LegacyDeviceFieldDataAccess, entity_SingleBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test(numDeviceIters, legacy_device_compute_centroid_entity_access,
                        device_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(LegacyDeviceFieldDataAccess, entity_MultiBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_multiple_block_test(numDeviceIters, legacy_device_compute_centroid_entity_access,
                          device_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(LegacyDeviceFieldDataAccess, entity_PartialBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_partial_block_test(numDeviceIters, legacy_device_compute_centroid_entity_access,
                         device_verify_averaged_centroids_are_center_of_mesh);
}

//------------------------------------------------------------------------------
TEST_F(DeviceFieldDataAccess, entity_SingleBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test(numDeviceIters, device_compute_centroid_entity_access,
                        device_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(DeviceFieldDataAccess, entity_MultiBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_multiple_block_test(numDeviceIters, device_compute_centroid_entity_access,
                          device_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(DeviceFieldDataAccess, entity_PartialBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_partial_block_test(numDeviceIters, device_compute_centroid_entity_access,
                         device_verify_averaged_centroids_are_center_of_mesh);
}

//------------------------------------------------------------------------------
TEST_F(DeviceFieldDataAccess, bucket_SingleBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test(numDeviceIters, device_compute_centroid_bucket_access,
                        device_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(DeviceFieldDataAccess, bucket_MultiBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_multiple_block_test(numDeviceIters, device_compute_centroid_bucket_access,
                          device_verify_averaged_centroids_are_center_of_mesh);
}

TEST_F(DeviceFieldDataAccess, bucket_PartialBlock)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_partial_block_test(numDeviceIters, device_compute_centroid_bucket_access,
                         device_verify_averaged_centroids_are_center_of_mesh);
}

}
