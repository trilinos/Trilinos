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
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldDataManager.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_topology/topology.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_io/FillMesh.hpp>
#include <vector>
#include <string>

namespace
{

double initial_value[3] = {-1, 2, -0.3};

class NgpFieldAccessPerformance : public stk::unit_test_util::MeshFixture
{
public:
  using DoubleVecField = stk::mesh::Field<double>;

  NgpFieldAccessPerformance()
    : batchTimer(get_comm())
  { }

  virtual ~NgpFieldAccessPerformance() {
    reset_mesh();
  }

  void setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                                std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager,
                                                unsigned initialBucketCapacity = stk::mesh::get_default_initial_bucket_capacity(),
                                                unsigned maximumBucketCapacity = stk::mesh::get_default_maximum_bucket_capacity())
  {
    stk::mesh::MeshBuilder builder(communicator);
    builder.set_spatial_dimension(m_spatialDim);
    builder.set_entity_rank_names(m_entityRankNames);
    builder.set_aura_option(auraOption);
    builder.set_field_data_manager(std::move(fieldDataManager));
    builder.set_initial_bucket_capacity(initialBucketCapacity);
    builder.set_maximum_bucket_capacity(maximumBucketCapacity);

        if(nullptr == metaData) {
          metaData = builder.create_meta_data();
        }

        if(nullptr == bulkData) {
          bulkData = builder.create(metaData);
          m_auraOption = auraOption;
          m_initialBucketCapacity = initialBucketCapacity;
          m_maximumBucketCapacity = maximumBucketCapacity;
        }
  }

  DoubleVecField * createNodalVectorField(const std::string &field_name)
  {
    const int spatialDimension = static_cast<int>(get_meta().spatial_dimension());
    DoubleVecField &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, field_name);
    stk::mesh::put_field_on_mesh(field, field.mesh_meta_data().universal_part(), spatialDimension, initial_value);
    return &field;
  }

  void createNodalVectorFields()
  {
    dispField  = createNodalVectorField("disp");
    velField   = createNodalVectorField("vel");
    accField   = createNodalVectorField("acc");
    forceField = createNodalVectorField("force");
  }

  void testVectorFieldSum(unsigned NUM_ITERS)
  {
    stk::mesh::NgpField<double> & ngpDispField  = stk::mesh::get_updated_ngp_field<double>(*dispField);
    stk::mesh::NgpField<double> & ngpVelField   = stk::mesh::get_updated_ngp_field<double>(*velField);
    stk::mesh::NgpField<double> & ngpAccField   = stk::mesh::get_updated_ngp_field<double>(*accField);
    stk::mesh::NgpField<double> & ngpForceField = stk::mesh::get_updated_ngp_field<double>(*forceField);
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    const int spatialDimension = static_cast<int>(get_meta().spatial_dimension());

      for (unsigned i = 0; i < NUM_ITERS; ++i) {
        stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, get_meta().locally_owned_part(),
                                       KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex & index)
                                       {
                                         for (int component = 0; component < spatialDimension; ++component) {
                                           ngpForceField(index, component) = alpha * ngpDispField(index, component) +
                                                                             beta * ngpVelField(index, component) +
                                                                             gamma * ngpAccField(index, component);
                                         }
                                       });
      }

    ngpForceField.modify_on_device();
  }

  void testHostVectorFieldSum(unsigned NUM_ITERS)
  {
    stk::mesh::get_updated_ngp_field<double>(*dispField);
    stk::mesh::get_updated_ngp_field<double>(*velField);
    stk::mesh::get_updated_ngp_field<double>(*accField);
    stk::mesh::get_updated_ngp_field<double>(*forceField);
    const int spatialDimension = static_cast<int>(get_meta().spatial_dimension());

      for (unsigned i = 0; i < NUM_ITERS; ++i) {
        const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().locally_owned_part());
        for (const stk::mesh::Bucket * bucket : buckets) {
          for (const stk::mesh::Entity & entity : *bucket) {
            const double * dispFieldData = stk::mesh::field_data(*dispField, entity);
            const double *  velFieldData = stk::mesh::field_data(*velField, entity);
            const double *  accFieldData = stk::mesh::field_data(*accField, entity);
            double * forceFieldData = stk::mesh::field_data(*forceField, entity);
            for (int component = 0; component < spatialDimension; ++component) {
              forceFieldData[component] = alpha * dispFieldData[component] +
                                          beta *   velFieldData[component] +
                                          gamma *  accFieldData[component];
            }
          }
        }
      }
  }

  void testPureHostVectorFieldSum(unsigned NUM_ITERS)
  {
    const int spatialDimension = static_cast<int>(get_meta().spatial_dimension());

      for (unsigned i = 0; i < NUM_ITERS; ++i) {
        const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().locally_owned_part());
        for (const stk::mesh::Bucket * bucket : buckets) {
          for (const stk::mesh::Entity & entity : *bucket) {
            const double * dispFieldData = stk::mesh::field_data(*dispField, entity);
            const double *  velFieldData = stk::mesh::field_data(*velField, entity);
            const double *  accFieldData = stk::mesh::field_data(*accField, entity);
            double * forceFieldData = stk::mesh::field_data(*forceField, entity);
            for (int component = 0; component < spatialDimension; ++component) {
              forceFieldData[component] = alpha * dispFieldData[component] +
                                          beta *   velFieldData[component] +
                                          gamma *  accFieldData[component];
            }
          }
        }
      }
  }

  void checkResult()
  {
    stk::mesh::NgpField<double> & ngpForceField = stk::mesh::get_updated_ngp_field<double>(*forceField);
    ngpForceField.sync_to_host();
    const int spatialDimension = static_cast<int>(get_meta().spatial_dimension());

    std::vector<double> expectedValues(spatialDimension);
    for (int i = 0; i < spatialDimension; ++i) {
      expectedValues[i] = alpha * initial_value[i] + beta * initial_value[i] + gamma * initial_value[i];
    }

    const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().locally_owned_part());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & node : *bucket) {
        const double * forceFieldData = stk::mesh::field_data(*forceField, node);
        for (int component = 0; component < spatialDimension; ++component) {
          EXPECT_DOUBLE_EQ(forceFieldData[component], expectedValues[component]);
        }
      }
    }
  }

  void checkHostResult()
  {
    const int spatialDimension = static_cast<int>(get_meta().spatial_dimension());

    std::vector<double> expectedValues(spatialDimension);
    for (int i = 0; i < spatialDimension; ++i) {
      expectedValues[i] = alpha * initial_value[i] + beta * initial_value[i] + gamma * initial_value[i];
    }

    const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().locally_owned_part());
    for (const stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const double * forceFieldData = stk::mesh::field_data(*forceField, entity);
        for (int component = 0; component < spatialDimension; ++component) {
          EXPECT_DOUBLE_EQ(forceFieldData[component], expectedValues[component]);
        }
      }
    }
  }

  static constexpr double alpha = -1.4;
  static constexpr double beta = 0.3333333;
  static constexpr double gamma = 3.14159;

  stk::unit_test_util::BatchTimer batchTimer;

  DoubleVecField * dispField;
  DoubleVecField * velField;
  DoubleVecField * accField;
  DoubleVecField * forceField;
};

TEST_F(NgpFieldAccessPerformance, pureHost_vectorSum_DefaultFieldDataManager)
{
  if (get_parallel_size() != 1) return;

  unsigned numElemsPerDim = 100;
  const int weKnowThereAreFiveRanks = 5;
  auto fieldDataManager = std::make_unique<stk::mesh::DefaultFieldDataManager>(weKnowThereAreFiveRanks);

  batchTimer.initialize_batch_timer();
  setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
  createNodalVectorFields();
  stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(numElemsPerDim), *bulkData);

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 1000;
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    testPureHostVectorFieldSum(NUM_ITERS);
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);

  checkHostResult();
}

TEST_F(NgpFieldAccessPerformance, host_vectorSum_DefaultFieldDataManager)
{
  if (get_parallel_size() != 1) return;

  unsigned numElemsPerDim = 100;
  const int weKnowThereAreFiveRanks = 5;
  auto fieldDataManager = std::make_unique<stk::mesh::DefaultFieldDataManager>(weKnowThereAreFiveRanks);

  batchTimer.initialize_batch_timer();
  setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
  createNodalVectorFields();
  stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(numElemsPerDim), *bulkData);

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 1000;
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    testHostVectorFieldSum(NUM_ITERS);
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);

  checkHostResult();
}

void fill_mesh(stk::mesh::BulkData& bulk, unsigned numElemsPerDim)
{
  stk::io::StkMeshIoBroker stkIo(MPI_COMM_WORLD);
  stkIo.set_bulk_data(bulk);
  stkIo.add_mesh_database(stk::unit_test_util::get_mesh_spec(numElemsPerDim), stk::io::READ_MESH);
  stkIo.create_input_mesh();
  const bool delayFieldDataAllocation = true;
  stkIo.populate_mesh(delayFieldDataAllocation);
  stkIo.populate_field_data();
}

TEST_F(NgpFieldAccessPerformance, vectorSum_DefaultFieldDataManager)
{
  if (get_parallel_size() != 1) return;

  unsigned numElemsPerDim = 100;
  const int weKnowThereAreFiveRanks = 5;
  auto fieldDataManager = std::make_unique<stk::mesh::DefaultFieldDataManager>(weKnowThereAreFiveRanks);

  batchTimer.initialize_batch_timer();
  setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
  createNodalVectorFields();
  fill_mesh(*bulkData, numElemsPerDim);

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 1000;
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    testVectorFieldSum(NUM_ITERS);
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);

  checkResult();
}

TEST_F(NgpFieldAccessPerformance, vectorSum_ContiguousFieldDataManager)
{
  if (get_parallel_size() != 1) return;

  unsigned numElemsPerDim = 100;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();

  batchTimer.initialize_batch_timer();
  setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
  createNodalVectorFields();
  fill_mesh(*bulkData, numElemsPerDim);

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 2000;
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    testVectorFieldSum(NUM_ITERS);
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);

  checkResult();
}

}

