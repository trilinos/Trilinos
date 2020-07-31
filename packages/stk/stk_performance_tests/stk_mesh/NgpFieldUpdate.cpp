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
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_performance_tests/stk_mesh/timer.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

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
    setup_mesh_with_fields("generated:100x100x100", stk::mesh::BulkData::NO_AUTO_AURA);
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

    stk::mesh::FieldBase* field = &get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEMENT_RANK, "FieldA");
    stk::mesh::put_field_on_mesh(*field, get_meta().universal_part(), &init);
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_mesh_with_fields(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    tensorField = &get_meta().declare_field<stk::mesh::Field<double, stk::mesh::Cartesian>>(stk::topology::ELEMENT_RANK, "TensorField");
    vectorField = &get_meta().declare_field<stk::mesh::Field<double, stk::mesh::Cartesian>>(stk::topology::ELEMENT_RANK, "VectorField");
    stk::mesh::put_field_on_mesh(*tensorField, get_meta().universal_part(), tensorFieldSizePerElem, static_cast<double*>(nullptr));
    stk::mesh::put_field_on_mesh(*vectorField, get_meta().universal_part(), vectorFieldSizePerElem, static_cast<double*>(nullptr));
    stk::performance_tests::setup_multiple_blocks(get_meta(), numElemBlocks);
    setup_mesh(meshSpecification, auraOption);
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemsPerDim, numElemBlocks);
  }

  void update_fields()
  {
    stk::mesh::NgpField<double>& ngpTensorField = stk::mesh::get_updated_ngp_field<double>(*tensorField);
    ngpTensorField.sync_to_device();

    stk::mesh::NgpField<double>& ngpVectorField = stk::mesh::get_updated_ngp_field<double>(*vectorField);
    ngpVectorField.sync_to_device();
  }

protected:
  stk::mesh::Field<double, stk::mesh::Cartesian>* tensorField;
  stk::mesh::Field<double, stk::mesh::Cartesian>* vectorField;
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
    get_meta().declare_part(newPartName);
    setup_mesh_with_fields("generated:100x100x100", stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void change_element_part_membership(int cycle)
  {
    get_bulk().modification_begin();
    get_bulk().change_entity_parts<stk::mesh::ConstPartVector>(get_element(cycle), {get_part()});
    get_bulk().modification_end();
    get_bulk().get_updated_ngp_mesh();
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

  void create_entity(int cycle)
  {
    get_bulk().modification_begin();
    get_bulk().declare_element(get_new_entity_id(cycle));
    get_bulk().modification_end();
    get_bulk().get_updated_ngp_mesh();
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
    setup_mesh_with_fields("generated:100x100x100", stk::mesh::BulkData::NO_AUTO_AURA);
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
    get_bulk().get_updated_ngp_mesh();
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

  const int NUM_RUNS = 100;

  stk::performance_tests::Timer timer(get_comm());
  setup_host_mesh();

  for (int i=0; i<NUM_RUNS; i++) {
    change_element_part_membership(i);
    timer.start_timing();
    update_fields();
    timer.update_timing();
  }
  timer.print_timing(NUM_RUNS);
}

TEST_F( NgpMeshCreateEntityWithFields, Timing )
{
  if (get_parallel_size() != 1) return;

  const int NUM_RUNS = 100;

  stk::performance_tests::Timer timer(get_comm());
  setup_host_mesh();

  for (int i=0; i<NUM_RUNS; i++) {
    create_entity(i);
    timer.start_timing();
    update_fields();
    timer.update_timing();
  }
  timer.print_timing(NUM_RUNS);
};

TEST_F( NgpMeshGhostingEntityWithFields, Timing )
{
  if (get_parallel_size() != 2) return;

  const int NUM_RUNS = 100;

  stk::performance_tests::Timer timer(get_comm());
  setup_host_mesh();

  for (int i=0; i<NUM_RUNS; i++) {
    ghost_element(i);
    timer.start_timing();
    update_fields();
    timer.update_timing();
  }
  timer.print_timing(NUM_RUNS);
}