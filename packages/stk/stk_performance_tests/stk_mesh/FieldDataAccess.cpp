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
#include <stk_util/util/FieldDataAllocator.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

namespace {

class FieldDataAccess : public stk::unit_test_util::MeshFixture
{
public:
  FieldDataAccess()
    : batchTimer(get_comm()),
      m_centroidField(nullptr),
      m_centroidFieldLeft(nullptr),
      m_centroidFieldRight(nullptr)
  { }

protected:
  stk::mesh::Field<double>& declare_vector_field(const std::string& fieldName)
  {
    stk::mesh::Field<double>& newField = get_meta().declare_field<double>(stk::topology::ELEM_RANK, fieldName);
    stk::mesh::put_field_on_mesh(newField, get_meta().universal_part(), 3, nullptr);
    return newField;
  }

  void declare_centroid_field()
  {
    m_centroidField = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "centroid");
    stk::mesh::put_field_on_mesh(*m_centroidField, get_meta().universal_part(), 3, nullptr);
  }

  void declare_centroid_field_left()
  {
    m_centroidFieldLeft = &get_meta().declare_field<double, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "centroid");
    stk::mesh::put_field_on_mesh(*m_centroidFieldLeft, get_meta().universal_part(), 3, nullptr);
  }

  void declare_centroid_field_right()
  {
    m_centroidFieldRight = &get_meta().declare_field<double, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "centroid");
    stk::mesh::put_field_on_mesh(*m_centroidFieldRight, get_meta().universal_part(), 3, nullptr);
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
  void run_single_block_test_left(int numIters, const CENTROID_FUNCTOR& centroidFunctor,
                                  const VERIFIER_FUNCTOR& verifierFunctor)
  {
    const unsigned NUM_RUNS = 5;
    const int ELEMS_PER_DIM = 100;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    declare_centroid_field_left();
    stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(ELEMS_PER_DIM), get_bulk());

    stk::mesh::Selector selector(get_meta().locally_owned_part());
    const stk::mesh::Field<double>& coordsField = *dynamic_cast<const stk::mesh::Field<double>*>(get_meta().coordinate_field());

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        centroidFunctor(selector, *m_centroidFieldLeft, coordsField);
      }
      batchTimer.stop_batch_timer();
    }

    verifierFunctor(ELEMS_PER_DIM, 1, 1, *m_centroidFieldLeft);

    batchTimer.print_batch_timing(numIters);
  }

  template <typename CENTROID_FUNCTOR, typename VERIFIER_FUNCTOR>
  void run_single_block_test_right(int numIters, const CENTROID_FUNCTOR& centroidFunctor,
                                   const VERIFIER_FUNCTOR& verifierFunctor)
  {
    const unsigned NUM_RUNS = 5;
    const int ELEMS_PER_DIM = 100;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    declare_centroid_field_right();
    stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(ELEMS_PER_DIM), get_bulk());

    stk::mesh::Selector selector(get_meta().locally_owned_part());
    const stk::mesh::Field<double>& coordsField = *dynamic_cast<const stk::mesh::Field<double>*>(get_meta().coordinate_field());

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        centroidFunctor(selector, *m_centroidFieldRight, coordsField);
      }
      batchTimer.stop_batch_timer();
    }

    verifierFunctor(ELEMS_PER_DIM, 1, 1, *m_centroidFieldRight);

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

  template <typename FIELD_DATA_FUNCTOR>
  void run_field_data_acquisition_test(int numIters, const FIELD_DATA_FUNCTOR& fieldDataFunctor)
  {
    const unsigned NUM_RUNS = 5;
    const int ELEMS_PER_DIM = 100;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Field<double>& field1 = declare_vector_field("field1");
    stk::mesh::Field<double>& field2 = declare_vector_field("field2");
    stk::mesh::Field<double>& field3 = declare_vector_field("field3");
    stk::mesh::Field<double>& field4 = declare_vector_field("field4");
    stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(ELEMS_PER_DIM), get_bulk());

    stk::mesh::Selector selector(get_meta().locally_owned_part());

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      fieldDataFunctor(numIters, field1, field2, field3, field4);
      batchTimer.stop_batch_timer();
    }

    batchTimer.print_batch_timing(numIters);
  }

  template <typename CENTROID_FUNCTOR, typename VERIFIER_FUNCTOR>
  void run_unified_memory_demo(int numRuns, int numIters, int numExtraFields, int elemsPerDim,
                               const CENTROID_FUNCTOR& centroidFunctor, const VERIFIER_FUNCTOR& verifierFunctor)
  {
    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    declare_centroid_field();

    if (get_bulk().parallel_rank() == 0) {
      std::cout << "Declaring " << numExtraFields << " extra Fields..." << std::endl;
    }
    for (int f = 0; f < numExtraFields; ++f) {
      auto& extraField = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "extraField" + std::to_string(f));
      stk::mesh::put_field_on_mesh(extraField, get_meta().universal_part(), 3, nullptr);
    }

    if (get_bulk().parallel_rank() == 0) {
      std::cout << "Building " << elemsPerDim << "x" << elemsPerDim << "x" << elemsPerDim << " mesh..." << std::endl;
    }
    stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(elemsPerDim), get_bulk());

    stk::mesh::Selector selector(get_meta().locally_owned_part());
    const stk::mesh::Field<double>& coordsField = *dynamic_cast<const stk::mesh::Field<double>*>(get_meta().coordinate_field());

    if (get_bulk().parallel_rank() == 0) {
      std::cout << "Running centroid calculations: " << numRuns << " runs with " << numIters << " iters..." << std::endl;
    }
    for (int run = 0; run < numRuns; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        centroidFunctor(selector, *m_centroidField, coordsField);
      }
      batchTimer.stop_batch_timer();
    }

    if (get_bulk().parallel_rank() == 0) {
      std::cout << "Verifying results..." << std::endl;
    }
    verifierFunctor(elemsPerDim, 1, 1, *m_centroidField);

    batchTimer.print_batch_timing(numIters);

    if (get_bulk().parallel_rank() == 0) {
      std::cout << std::endl;
    }
    stk::parallel_print_time_without_output_and_hwm(get_comm(), batchTimer.get_min_batch_time(), std::cout);
  }

  void run_bucket_allocation_test(int meshSize, int bucketCapacity, int numNodeFields, int numElemFields,
                                  bool printPointers, bool allocateOnDevice)
  {
    batchTimer.initialize_batch_timer();
    stk::FieldDataAllocator<std::byte> fieldDataAllocator;

    const int numNodes = (meshSize+1) * (meshSize+1) * (meshSize+1);
    const int numElems = meshSize * meshSize * meshSize;

    const int numNodeBuckets = (numNodes + bucketCapacity - 1) / bucketCapacity;
    const int numElemBuckets = (numElems + bucketCapacity - 1) / bucketCapacity;

    const int nodeBucketSize = bucketCapacity * numNodeFields * (3 * sizeof(double));  // Vector double Fields
    const int elemBucketSize = bucketCapacity * numElemFields * (3 * sizeof(double));  // Vector double Fields

    std::cout << "Using mesh size of " << meshSize << " to generate " << numNodes << " nodes and "
              << numElems << " elements" << std::endl;
    std::cout << "Bucket capacity is " << bucketCapacity << std::endl;
    std::cout << "Using " << numNodeFields << " Node Fields" << std::endl;
    std::cout << "Using " << numElemFields << " Elem Fields" << std::endl;
    std::cout << "Allocating " << numNodeBuckets << " Node Buckets of " << nodeBucketSize << " bytes each" << std::endl;
    std::cout << "Allocating " << numElemBuckets << " Elem Buckets of " << elemBucketSize << " bytes each" << std::endl;

    const int nodeBucketTotalSize = numNodeBuckets * nodeBucketSize;
    const int elemBucketTotalSize = numElemBuckets * elemBucketSize;
    const double totalAllocationGB = static_cast<double>(nodeBucketTotalSize + elemBucketTotalSize) /
                                     (1024.0 * 1024.0 * 1024.0);

    std::cout << "Allocating " << totalAllocationGB << " GB total RAM..." << std::endl;

    batchTimer.start_batch_timer();
    if (allocateOnDevice) {
      using AllocationType = stk::FieldDataAllocator<std::byte>::DeviceAllocationType;

      std::vector<AllocationType> nodeAllocations;
      nodeAllocations.reserve(numNodeBuckets);
      for (int i = 0; i < numNodeBuckets; ++i) {
        auto nodeAllocation = fieldDataAllocator.device_allocate(nodeBucketSize);
        if (printPointers) {
          printf("###,%p,%i\n", (void*)nodeAllocation.data(), nodeBucketSize);  // Greppable, CSV-formatted output
        }
        nodeAllocations.push_back(nodeAllocation);
      }
      std::vector<AllocationType> elemAllocations;
      elemAllocations.reserve(numElemBuckets);
      for (int i = 0; i < numElemBuckets; ++i) {
        auto elemAllocation = fieldDataAllocator.device_allocate(elemBucketSize);
        if (printPointers) {
          printf("###,%p,%i\n", (void*)elemAllocation.data(), elemBucketSize);  // Greppable, CSV-formatted output
        }
        elemAllocations.push_back(elemAllocation);
      }
      std::cout << "Deallocating storage..." << std::endl;
    }
    else {
      using AllocationType = stk::FieldDataAllocator<std::byte>::HostAllocationType;

      std::vector<AllocationType> nodeAllocations;
      nodeAllocations.reserve(numNodeBuckets);
      for (int i = 0; i < numNodeBuckets; ++i) {
        auto nodeAllocation = fieldDataAllocator.host_allocate(nodeBucketSize);
        if (printPointers) {
          printf("###,%p,%i\n", (void*)nodeAllocation.data(), nodeBucketSize);  // Greppable, CSV-formatted output
        }
        nodeAllocations.push_back(nodeAllocation);
      }
      std::vector<AllocationType> elemAllocations;
      elemAllocations.reserve(numElemBuckets);
      for (int i = 0; i < numElemBuckets; ++i) {
        auto elemAllocation = fieldDataAllocator.host_allocate(elemBucketSize);
        if (printPointers) {
          printf("###,%p,%i\n", (void*)elemAllocation.data(), elemBucketSize);  // Greppable, CSV-formatted output
        }
        elemAllocations.push_back(elemAllocation);
      }
      std::cout << "Deallocating storage..." << std::endl;
    }
    batchTimer.stop_batch_timer();
    batchTimer.print_batch_timing(1);
  }

  void run_memory_allocation_test(int allocationSize, int numAllocations, bool printPointers, bool allocateOnDevice)
  {
    batchTimer.initialize_batch_timer();
    stk::FieldDataAllocator<std::byte> fieldDataAllocator;

    std::cout << "Allocating " << numAllocations << " segments of " << allocationSize << " bytes..." << std::endl;

    batchTimer.start_batch_timer();
    if (allocateOnDevice) {
      using AllocationType = stk::FieldDataAllocator<std::byte>::DeviceAllocationType;

      std::vector<AllocationType> allocations;
      allocations.reserve(numAllocations);
      for (int i = 0; i < numAllocations; ++i) {
        auto allocation = fieldDataAllocator.device_allocate(allocationSize);
        if (printPointers) {
          printf("###,%p,%i\n", (void*)allocation.data(), allocationSize);  // Greppable, CSV-formatted output
        }
        allocations.push_back(allocation);
      }
      std::cout << "Deallocating storage..." << std::endl;
    }
    else {
      using AllocationType = stk::FieldDataAllocator<std::byte>::HostAllocationType;

      std::vector<AllocationType> allocations;
      allocations.reserve(numAllocations);
      for (int i = 0; i < numAllocations; ++i) {
        auto allocation = fieldDataAllocator.host_allocate(allocationSize);
        if (printPointers) {
          printf("###,%p,%i\n", (void*)allocation.data(), allocationSize);  // Greppable, CSV-formatted output
        }
        allocations.push_back(allocation);
      }
      std::cout << "Deallocating storage..." << std::endl;
    }
    batchTimer.stop_batch_timer();
    batchTimer.print_batch_timing(1);
  }

  stk::unit_test_util::BatchTimer batchTimer;
  stk::mesh::Field<double> *m_centroidField;
  stk::mesh::Field<double, stk::mesh::Layout::Left> *m_centroidFieldLeft;
  stk::mesh::Field<double, stk::mesh::Layout::Right> *m_centroidFieldRight;
};

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

auto host_verify_averaged_centroids_are_center_of_mesh_left = [](int elemsPerDim, int numTotalBlocks, int numUsedBlocks,
                                                                 stk::mesh::Field<double, stk::mesh::Layout::Left>& centroidField)
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

auto host_verify_averaged_centroids_are_center_of_mesh_right = [](int elemsPerDim, int numTotalBlocks, int numUsedBlocks,
                                                                 stk::mesh::Field<double, stk::mesh::Layout::Right>& centroidField)
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


#ifdef NDEBUG
  constexpr int numHostIters = 100;
#else
  constexpr int numHostIters = 1;
#endif

#ifdef STK_ENABLE_GPU
  constexpr int numDeviceIters = 2000;
#else
  constexpr int numDeviceIters = 100;
#endif


#ifndef STK_UNIFIED_MEMORY

class LegacyFieldDataAccess : public FieldDataAccess {};

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


class LegacyDeviceFieldDataAccess : public FieldDataAccess {};

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

#endif  // For not STK_UNIFIED_MEMORY


auto compute_centroid_entity_access = [](const stk::mesh::Selector& selector,
                                         stk::mesh::Field<double>& centroidField,
                                         const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite>();
  auto coordsData = coordsField.data();

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

auto compute_centroid_entity_access_left = [](const stk::mesh::Selector& selector,
                                              stk::mesh::Field<double, stk::mesh::Layout::Left>& centroidField,
                                              const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite>();
  auto coordsData = coordsField.data();

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

auto compute_centroid_entity_access_right = [](const stk::mesh::Selector& selector,
                                               stk::mesh::Field<double, stk::mesh::Layout::Right>& centroidField,
                                               const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite>();
  auto coordsData = coordsField.data();

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

auto compute_centroid_entity_access_left_auto = [](const stk::mesh::Selector& selector,
                                                   stk::mesh::Field<double, stk::mesh::Layout::Left>& centroidField,
                                                   const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  stk::mesh::FieldBase& centroidFieldBase = static_cast<stk::mesh::FieldBase&>(centroidField);
  auto centroidData = centroidFieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>();
  auto coordsData = coordsField.data();

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

auto compute_centroid_entity_access_right_auto = [](const stk::mesh::Selector& selector,
                                                    stk::mesh::Field<double, stk::mesh::Layout::Right>& centroidField,
                                                    const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  stk::mesh::FieldBase& centroidFieldBase = static_cast<stk::mesh::FieldBase&>(centroidField);
  auto centroidData = centroidFieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>();
  auto coordsData = coordsField.data();

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

auto compute_centroid_bucket_access = [](const stk::mesh::Selector& selector,
                                         stk::mesh::Field<double>& centroidField,
                                         const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite>();
  auto coordsData = coordsField.data();

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

auto compute_centroid_bucket_access_left = [](const stk::mesh::Selector& selector,
                                              stk::mesh::Field<double, stk::mesh::Layout::Left>& centroidField,
                                              const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite>();
  auto coordsData = coordsField.data();

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

auto compute_centroid_bucket_access_right = [](const stk::mesh::Selector& selector,
                                              stk::mesh::Field<double, stk::mesh::Layout::Right>& centroidField,
                                               const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite>();
  auto coordsData = coordsField.data();

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

auto compute_centroid_bucket_access_left_auto = [](const stk::mesh::Selector& selector,
                                                   stk::mesh::Field<double, stk::mesh::Layout::Left>& centroidField,
                                                   const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  stk::mesh::FieldBase& centroidFieldBase = static_cast<stk::mesh::FieldBase&>(centroidField);
  auto centroidData = centroidFieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>();
  auto coordsData = coordsField.data();

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

auto compute_centroid_bucket_access_right_auto = [](const stk::mesh::Selector& selector,
                                                    stk::mesh::Field<double, stk::mesh::Layout::Right>& centroidField,
                                                    const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  stk::mesh::FieldBase& centroidFieldBase = static_cast<stk::mesh::FieldBase&>(centroidField);
  auto centroidData = centroidFieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>();
  auto coordsData = coordsField.data();

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

auto host_field_data_acquisition = [](int numIters, stk::mesh::Field<double>& field1, stk::mesh::Field<double>& field2,
                                      stk::mesh::Field<double>& field3, stk::mesh::Field<double>& field4)
{
  for (int iter = 0; iter < numIters; ++iter) {
    auto fieldData1 = field1.data();
    auto fieldData2 = field2.data();
    auto fieldData3 = field3.data();
    auto fieldData4 = field4.data();
  }
};

auto host_entity_values_acquisition = [](int numIters, stk::mesh::Field<double>& field1,
                                         stk::mesh::Field<double>& field2, stk::mesh::Field<double>& field3,
                                         stk::mesh::Field<double>& field4)
{
  auto fieldData1 = field1.data();
  auto fieldData2 = field2.data();
  auto fieldData3 = field3.data();
  auto fieldData4 = field4.data();
  const stk::mesh::BulkData& bulk = field1.get_mesh();
  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1);

  for (int iter = 0; iter < numIters; ++iter) {
    auto entityValues1 = fieldData1.entity_values(elem);
    asm volatile("" : : "r,m"(entityValues1) : "memory");  // Prevent from optimizing away by "using" object in volatile do-nothing operation

    auto entityValues2 = fieldData2.entity_values(elem);
    asm volatile("" : : "r,m"(entityValues2) : "memory");

    auto entityValues3 = fieldData3.entity_values(elem);
    asm volatile("" : : "r,m"(entityValues3) : "memory");

    auto entityValues4 = fieldData4.entity_values(elem);
    asm volatile("" : : "r,m"(entityValues4) : "memory");
  }
};

auto host_bucket_values_acquisition = [](int numIters, stk::mesh::Field<double>& field1,
                                         stk::mesh::Field<double>& field2, stk::mesh::Field<double>& field3,
                                         stk::mesh::Field<double>& field4)
{
  auto fieldData1 = field1.data();
  auto fieldData2 = field2.data();
  auto fieldData3 = field3.data();
  auto fieldData4 = field4.data();
  const stk::mesh::BulkData& bulk = field1.get_mesh();
  stk::mesh::Bucket& bucket = bulk.bucket(bulk.get_entity(stk::topology::ELEM_RANK, 1));

  for (int iter = 0; iter < numIters; ++iter) {
    auto bucketValues1 = fieldData1.bucket_values(bucket);
    asm volatile("" : : "r,m"(bucketValues1) : "memory");  // Prevent from optimizing away by "using" object in volatile do-nothing operation

    auto bucketValues2 = fieldData2.bucket_values(bucket);
    asm volatile("" : : "r,m"(bucketValues2) : "memory");

    auto bucketValues3 = fieldData3.bucket_values(bucket);
    asm volatile("" : : "r,m"(bucketValues3) : "memory");

    auto bucketValues4 = fieldData4.bucket_values(bucket);
    asm volatile("" : : "r,m"(bucketValues4) : "memory");
  }
};


//------------------------------------------------------------------------------
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
TEST_F(FieldDataAccess, entity_SingleBlock_layoutLeft)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test_left(numHostIters, compute_centroid_entity_access_left,
                             host_verify_averaged_centroids_are_center_of_mesh_left);
}

TEST_F(FieldDataAccess, entity_SingleBlock_layoutRight)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test_right(numHostIters, compute_centroid_entity_access_right,
                              host_verify_averaged_centroids_are_center_of_mesh_right);
}

//------------------------------------------------------------------------------
TEST_F(FieldDataAccess, entity_SingleBlock_layoutLeftAuto)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test_left(numHostIters, compute_centroid_entity_access_left_auto,
                             host_verify_averaged_centroids_are_center_of_mesh_left);
}

TEST_F(FieldDataAccess, entity_SingleBlock_layoutRightAuto)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test_right(numHostIters, compute_centroid_entity_access_right_auto,
                              host_verify_averaged_centroids_are_center_of_mesh_right);
}


//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
TEST_F(FieldDataAccess, bucket_SingleBlock_layoutLeft)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test_left(numHostIters, compute_centroid_bucket_access_left,
                             host_verify_averaged_centroids_are_center_of_mesh_left);
}

TEST_F(FieldDataAccess, bucket_SingleBlock_layoutRight)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test_right(numHostIters, compute_centroid_bucket_access_right,
                              host_verify_averaged_centroids_are_center_of_mesh_right);
}

//------------------------------------------------------------------------------
TEST_F(FieldDataAccess, bucket_SingleBlock_layoutLeftAuto)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test_left(numHostIters, compute_centroid_bucket_access_left_auto,
                             host_verify_averaged_centroids_are_center_of_mesh_left);
}

TEST_F(FieldDataAccess, bucket_SingleBlock_layoutRightAuto)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_single_block_test_right(numHostIters, compute_centroid_bucket_access_right_auto,
                              host_verify_averaged_centroids_are_center_of_mesh_right);
}


//------------------------------------------------------------------------------
TEST_F(FieldDataAccess, FieldDataAcquisition)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_field_data_acquisition_test(20'000'000, host_field_data_acquisition);
}

//------------------------------------------------------------------------------
TEST_F(FieldDataAccess, EntityValuesAcquisition)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_field_data_acquisition_test(500'000'000, host_entity_values_acquisition);
}

//------------------------------------------------------------------------------
TEST_F(FieldDataAccess, BucketValuesAcquisition)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_field_data_acquisition_test(500'000'000, host_bucket_values_acquisition);
}


//------------------------------------------------------------------------------

class DeviceFieldDataAccess : public FieldDataAccess {};

void device_compute_centroid_entity_access_function(const stk::mesh::Selector& selector,
                                                    stk::mesh::Field<double>& centroidField,
                                                    const stk::mesh::Field<double>& coordsField)
{
  const stk::mesh::BulkData& bulk = centroidField.get_mesh();
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  auto centroidData = centroidField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  auto coordsData = coordsField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

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

  auto centroidData = centroidField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  auto coordsData = coordsField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

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

void device_field_data_acquisition_function(int numIters, stk::mesh::Field<double>& field1,
                                            stk::mesh::Field<double>& field2, stk::mesh::Field<double>& field3,
                                            stk::mesh::Field<double>& field4)
{
  for (int iter = 0; iter < numIters; ++iter) {
    [[maybe_unused]] auto fieldData1 = field1.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    [[maybe_unused]] auto fieldData2 = field2.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    [[maybe_unused]] auto fieldData3 = field3.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    [[maybe_unused]] auto fieldData4 = field4.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  }
}

auto device_field_data_acquisition = [](int numIters, stk::mesh::Field<double>& field1, stk::mesh::Field<double>& field2,
                                        stk::mesh::Field<double>& field3, stk::mesh::Field<double>& field4)
{
  device_field_data_acquisition_function(numIters, field1, field2, field3, field4);
};

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

//------------------------------------------------------------------------------
TEST_F(DeviceFieldDataAccess, FieldDataAcquisition)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_field_data_acquisition_test(20'000'000, device_field_data_acquisition);
}


class UnifiedMemoryDemonstration : public DeviceFieldDataAccess {};

//------------------------------------------------------------------------------
TEST_F(UnifiedMemoryDemonstration, centroidTest)
{
  int numRuns = stk::unit_test_util::get_command_line_option<int>("--num-runs", 1);
  int numIters = stk::unit_test_util::get_command_line_option<int>("--num-iters", 1);
  int numExtraFields = stk::unit_test_util::get_command_line_option<int>("--num-extra-fields", 85);
  int elemsPerDim = stk::unit_test_util::get_command_line_option<int>("--elems-per-dim", 480);

  run_unified_memory_demo(numRuns, numIters, numExtraFields, elemsPerDim,
                          device_compute_centroid_entity_access, device_verify_averaged_centroids_are_center_of_mesh);
}

//------------------------------------------------------------------------------
TEST_F(UnifiedMemoryDemonstration, bucketAllocationTestHost)
{
  int meshSize = stk::unit_test_util::get_command_line_option<int>("--mesh-size", 480);
  int bucketCapacity = stk::unit_test_util::get_command_line_option<int>("--bucket-capacity", 512);
  int numNodeFields = stk::unit_test_util::get_command_line_option<int>("--num-node-fields", 1);
  int numElemFields = stk::unit_test_util::get_command_line_option<int>("--num-elem-fields", 86);
  bool printPointers = stk::unit_test_util::get_command_line_option<bool>("--print-pointers", false);

  const bool allocateOnDevice = false;
  run_bucket_allocation_test(meshSize, bucketCapacity, numNodeFields, numElemFields, printPointers, allocateOnDevice);
}

TEST_F(UnifiedMemoryDemonstration, bucketAllocationTestDevice)
{
  int meshSize = stk::unit_test_util::get_command_line_option<int>("--mesh-size", 480);
  int bucketCapacity = stk::unit_test_util::get_command_line_option<int>("--bucket-capacity", 512);
  int numNodeFields = stk::unit_test_util::get_command_line_option<int>("--num-node-fields", 1);
  int numElemFields = stk::unit_test_util::get_command_line_option<int>("--num-elem-fields", 86);
  bool printPointers = stk::unit_test_util::get_command_line_option<bool>("--print-pointers", false);

  const bool allocateOnDevice = true;
  run_bucket_allocation_test(meshSize, bucketCapacity, numNodeFields, numElemFields, printPointers, allocateOnDevice);
}

//------------------------------------------------------------------------------
TEST_F(UnifiedMemoryDemonstration, memoryAllocationTestHost)
{
  int allocationSize = stk::unit_test_util::get_command_line_option<int>("--allocation-size", 1024);
  int numAllocations = stk::unit_test_util::get_command_line_option<int>("--num-allocations", 100);
  bool printPointers = stk::unit_test_util::get_command_line_option<bool>("--print-pointers", false);

  const bool allocateOnDevice = false;
  run_memory_allocation_test(allocationSize, numAllocations, printPointers, allocateOnDevice);
}

TEST_F(UnifiedMemoryDemonstration, memoryAllocationTestDevice)
{
  int allocationSize = stk::unit_test_util::get_command_line_option<int>("--allocation-size", 1024);
  int numAllocations = stk::unit_test_util::get_command_line_option<int>("--num-allocations", 100);
  bool printPointers = stk::unit_test_util::get_command_line_option<bool>("--print-pointers", false);

  const bool allocateOnDevice = true;
  run_memory_allocation_test(allocationSize, numAllocations, printPointers, allocateOnDevice);
}

}
