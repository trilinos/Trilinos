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

#include <string>
#include <ostream>
#include <gtest/gtest.h>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

class NgpFieldSyncTest : public stk::unit_test_util::MeshFixture
{
public:
  NgpFieldSyncTest() : stk::unit_test_util::MeshFixture()
  {}

  void setup_mesh_with_many_blocks_many_elements(unsigned numBlocks, unsigned numElemPerDim)
  {
    std::string meshDesc = "generated:" + std::to_string(numElemPerDim) + "x"
                                        + std::to_string(numElemPerDim) + "x"
                                        + std::to_string(numElemPerDim);
    stk::performance_tests::setup_multiple_blocks(get_meta(), numBlocks);
    stk::io::fill_mesh(meshDesc, get_bulk());
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemPerDim);
  }

  void setup_fields(unsigned numComponent, unsigned tensorFieldSizePerElem, unsigned vectorFieldSizePerElem)
  {
    const std::vector<int> init(numComponent, 1);
    auto tensorField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "TensorField");
    auto vectorField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "VectorField");
    auto stkIntField = &get_meta().declare_field<int>(stk::topology::ELEMENT_RANK, "intField", 1);

    stk::mesh::put_field_on_mesh(*tensorField, get_meta().universal_part(), tensorFieldSizePerElem, nullptr);
    stk::mesh::put_field_on_mesh(*vectorField, get_meta().universal_part(), vectorFieldSizePerElem, nullptr);
    stk::mesh::put_field_on_mesh(*stkIntField, get_meta().universal_part(), numComponent, init.data());
  }

  stk::mesh::Selector get_non_contiguous_partial_selector(unsigned numBlocksToSync)
  {
    stk::mesh::PartVector parts;
    unsigned partCount = 0;

    stk::mesh::PartVector elemBlockParts;
    stk::mesh::fill_element_block_parts(get_meta(), stk::topology::HEX_8, elemBlockParts);
    const unsigned numBlocks = elemBlockParts.size();

    for(unsigned i = 0; i < numBlocks; i++) {
      if(i % 2 == 0) {
        if(partCount == numBlocksToSync) { break; }
        parts.push_back(elemBlockParts[i]);
        partCount++;
      }
    }

    for(unsigned i = 0; i < numBlocks; i++) {
      if(i % 2 == 1) {
        if(partCount == numBlocksToSync) { break; }
        parts.push_back(elemBlockParts[i]);
        partCount++;
      }
    }

    return stk::mesh::selectUnion(parts);
  }

  stk::mesh::Selector get_contiguous_partial_selector(unsigned numBlocksToSync)
  {
    stk::mesh::PartVector parts;

    stk::mesh::PartVector elemBlockParts;
    stk::mesh::fill_element_block_parts(get_meta(), stk::topology::HEX_8, elemBlockParts);
    const unsigned numBlocks = elemBlockParts.size();
    STK_ThrowRequire(numBlocksToSync <= numBlocks);

    for(unsigned i = 0; i < numBlocksToSync; i++) {
      parts.push_back(elemBlockParts[i]);
    }

    return stk::mesh::selectUnion(parts);
  }
};

class NgpFieldUpdateFixture : public stk::unit_test_util::MeshFixture
{
public:
 NgpFieldUpdateFixture()
   : stk::unit_test_util::MeshFixture(),
      tensorField(nullptr),
      vectorField(nullptr),
      tensorFieldSizePerElem(72),
      vectorFieldSizePerElem(8),
      numElemBlocks(100),
      numElemsPerDim(100),
      numElements(std::pow(numElemsPerDim, 3))
    {}

  virtual void setup_host_mesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    setup_mesh_with_fields("generated:100x100x100");
  }

  std::string generate_stacked_block_mesh_desc(unsigned numBlocks)
  {
    std::string meshDesc;

    for(unsigned i = 0; i < numBlocks; i++) {
      meshDesc += get_nodal_string_for_block(i+1);
      if(i != numBlocks-1) {
        meshDesc += "\n";
      }
    }
    return meshDesc;
  }

  std::string get_nodal_string_for_block(unsigned blockId)
  {
    std::ostringstream blockStr;
    for(unsigned i = 0; i < 2; i++) {
      blockStr << "0," << blockId*2+i-1 << ",HEX_8,";
      for(unsigned j = 0; j < 8; j++) {
        blockStr << j+1 + ((blockId-1)*2 + i) * 8 << ",";
      }
      blockStr << "block_" << blockId;
      if(i != 1) {
       blockStr << "\n";
      }
    }
    return blockStr.str();
  }
  
  void setup_mesh_with_stacked_blocks(unsigned numBlocks)
  {
    double init = 0.0;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    std::string meshDesc = generate_stacked_block_mesh_desc(numBlocks);

    stk::mesh::FieldBase* field = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "FieldA");
    stk::mesh::put_field_on_mesh(*field, get_meta().universal_part(), &init);
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_mesh_with_fields(const std::string &meshSpecification)
  {
    tensorField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "TensorField");
    vectorField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "VectorField");
    stk::mesh::put_field_on_mesh(*tensorField, get_meta().universal_part(), tensorFieldSizePerElem, nullptr);
    stk::mesh::put_field_on_mesh(*vectorField, get_meta().universal_part(), vectorFieldSizePerElem, nullptr);
    stk::performance_tests::setup_multiple_blocks(get_meta(), numElemBlocks);
    stk::io::fill_mesh(meshSpecification, get_bulk());
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemsPerDim);
  }

  void update_fields()
  {
    stk::mesh::NgpField<double>& ngpTensorField = stk::mesh::get_updated_ngp_field<double>(*tensorField);
    ngpTensorField.sync_to_device();

    stk::mesh::NgpField<double>& ngpVectorField = stk::mesh::get_updated_ngp_field<double>(*vectorField);
    ngpVectorField.sync_to_device();
  }

protected:
  stk::mesh::Field<double>* tensorField;
  stk::mesh::Field<double>* vectorField;
  unsigned tensorFieldSizePerElem;
  unsigned vectorFieldSizePerElem;
  unsigned numElemBlocks;
  unsigned numElemsPerDim;
  unsigned numElements;
};

class NgpMeshChangeElementPartMembershipWithFields : public NgpFieldUpdateFixture
{
public:
  NgpMeshChangeElementPartMembershipWithFields()
    : NgpFieldUpdateFixture()
  { }

  void setup_host_mesh() override
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    get_meta().declare_part(newPartName);
    setup_mesh_with_fields("generated:100x100x100");
  }

  void change_element_part_membership(int cycle)
  {
    get_bulk().modification_begin();
    const stk::mesh::Part* part = get_part();
    STK_ThrowRequireMsg(part!=nullptr,"get_part returned nullptr, newPartName="<<newPartName);
    get_bulk().change_entity_parts<stk::mesh::ConstPartVector>(get_element(cycle), {get_part()});
    get_bulk().modification_end();
    stk::mesh::get_updated_ngp_mesh(get_bulk());
  }

private:
  stk::mesh::Entity get_element(int cycle)
  {
    stk::mesh::EntityId elemId = cycle+1;
    return get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
  }

  const stk::mesh::Part* get_part()
  {
    return get_meta().get_part(newPartName);
  }

  std::string newPartName;
};

class NgpMeshCreateEntityWithFields : public NgpFieldUpdateFixture 
{
public:
  NgpMeshCreateEntityWithFields()
    : NgpFieldUpdateFixture()
  { }

  void create_entity(int cycle, int numElemsToCreatePerModCycle)
  {
    get_bulk().modification_begin();
    for(int i=0; i<numElemsToCreatePerModCycle; ++i) {
      const int index = cycle*numElemsToCreatePerModCycle + i;
      get_bulk().declare_element(get_new_entity_id(index));
    }
    get_bulk().modification_end();
    stk::mesh::get_updated_ngp_mesh(get_bulk());
  }

private:
  stk::mesh::EntityId get_new_entity_id(int cycle)
  {
    return numElements + cycle + 1;
  }
};

class NgpMeshGhostingEntityWithFields : public NgpFieldUpdateFixture 
{
public:
  NgpMeshGhostingEntityWithFields()
    : NgpFieldUpdateFixture(),
      ghostingName("myGhosting"),
      ghosting(nullptr)
  { }

  void setup_host_mesh() override
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    setup_mesh_with_fields("generated:100x100x100");
    get_bulk().modification_begin();
    ghosting = &get_bulk().create_ghosting(ghostingName);
    get_bulk().modification_end();
  } 

protected:
  void ghost_element(int cycle)
  {
    get_bulk().modification_begin();
    get_bulk().change_ghosting(*ghosting, element_to_ghost(cycle));
    get_bulk().modification_end();
    stk::mesh::get_updated_ngp_mesh(get_bulk());
  }

private:
  stk::mesh::EntityProcVec element_to_ghost(int cycle)
  {
    stk::mesh::EntityId firstLocalElemId = get_parallel_rank()*numElements/2 + 1;
    stk::mesh::EntityId elemId = firstLocalElemId + cycle;
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    int ghostRank = 1 - get_parallel_rank();
    return {{elem, ghostRank}};
  }

  std::string ghostingName;
  stk::mesh::Ghosting* ghosting;
};

TEST_F( NgpMeshChangeElementPartMembershipWithFields, Timing )
{
  if (get_parallel_size() != 1) return;

  const unsigned NUM_RUNS = 5;
  #ifdef STK_ENABLE_GPU
  const int NUM_ITERS = 50;
  #else
  const int NUM_ITERS = 2500;
  #endif
  const int NUM_FAKE_ITERS = 2500;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());

  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_host_mesh();
    batchTimer.start_batch_timer();

    for (int i = 0; i < NUM_ITERS; i++) {
      change_element_part_membership(i);
      update_fields();
    }
    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_FAKE_ITERS);
}

TEST_F( NgpMeshCreateEntityWithFields, Timing )
{
  if (get_parallel_size() != 1) return;

  const unsigned NUM_RUNS = 5;
  #ifdef STK_ENABLE_GPU
  const int numModCycles = 50;
  #else
  const int numModCycles = 4000;
  #endif
  const int numFakeModCycles = 4000;
  const int numElemsToCreatePerModCycle = 40;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_host_mesh();
    batchTimer.start_batch_timer();


    for (int i=0; i<numModCycles; i++) {
      create_entity(i, numElemsToCreatePerModCycle);
      update_fields();
    }
    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(numFakeModCycles);
};

TEST_F( NgpMeshGhostingEntityWithFields, Timing )
{
  if (get_parallel_size() != 2) return;

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 500;

  stk::unit_test_util::BatchTimer batchTimer(get_comm());

  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_host_mesh();
    batchTimer.start_batch_timer();
    for (int i = 0; i < NUM_ITERS; i++) {
      ghost_element(i);
      update_fields();
    }
    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F(NgpFieldSyncTest, PartialSyncTiming)
{
  if(get_parallel_size() != 1) return;

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 20000000);
  unsigned numComponents = stk::unit_test_util::get_command_line_option("-c", 1);
  unsigned numBlocks = stk::unit_test_util::get_command_line_option("-b", 20);
  unsigned numBlocksToSync = stk::unit_test_util::get_command_line_option("-s", 5);
  unsigned numElemPerDim = stk::unit_test_util::get_command_line_option("-e", 100);
  unsigned tensorFieldSizePerElem = stk::unit_test_util::get_command_line_option("--tensorField", 72);
  unsigned vectorFieldSizePerElem = stk::unit_test_util::get_command_line_option("-vectorField", 8);
  bool justSyncAll = stk::unit_test_util::get_command_line_option("-a", false);
  bool contiguousBlocks = stk::unit_test_util::get_command_line_option("-t", true);
  numBlocksToSync = std::min(numBlocks, numBlocksToSync);
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  setup_fields(numComponents, tensorFieldSizePerElem, vectorFieldSizePerElem);
  setup_mesh_with_many_blocks_many_elements(numBlocks, numElemPerDim);

  stk::mesh::FieldBase* fieldBase = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField");
  stk::mesh::FieldBase* tensorFieldBase = get_meta().get_field(stk::topology::ELEMENT_RANK, "TensorField");
  stk::mesh::FieldBase* vectorFieldBase = get_meta().get_field(stk::topology::ELEMENT_RANK, "VectorField");
  stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(*fieldBase);
  stk::mesh::NgpField<double>& ngpTensorField = stk::mesh::get_updated_ngp_field<double>(*tensorFieldBase);
  stk::mesh::NgpField<double>& ngpVectorField = stk::mesh::get_updated_ngp_field<double>(*vectorFieldBase);
  stk::mesh::Selector selector;
  
  if(contiguousBlocks) {
    selector = get_contiguous_partial_selector(numBlocksToSync);
  } else {
    selector = get_non_contiguous_partial_selector(numBlocksToSync);
  }

  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    for(unsigned i = 0; i < NUM_ITERS; i++) {
      if(justSyncAll) {
        ngpIntField.modify_on_host();
        ngpTensorField.modify_on_host();
        ngpVectorField.modify_on_host();
      } 
      else {
        ngpIntField.modify_on_host(selector);
        ngpTensorField.modify_on_host(selector);
        ngpVectorField.modify_on_host(selector);
      }
      ngpIntField.sync_to_device();
      ngpTensorField.sync_to_device();
      ngpVectorField.sync_to_device();
    }
      batchTimer.stop_batch_timer();
  }

  std::cout << "Blocks: " << numBlocks << " Blocks to sync: "
            << numBlocksToSync << " (" << (numBlocksToSync * 100.0 / numBlocks) << " %)" << std::endl;
  batchTimer.print_batch_timing(NUM_ITERS);
}
