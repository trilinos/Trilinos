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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_NE, etc
#include <stddef.h>                     // for size_t, NULL
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/GetEntities.hpp>
#include <string>                       // for string
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_util/parallel/Parallel.hpp"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_io/WriteMesh.hpp"

#include <stk_mesh/base/Field.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace stk { namespace mesh { class Part; } }
namespace {
using stk::unit_test_util::build_mesh;

void compare_fields(const stk::mesh::BulkData& bulk1,
                    const stk::mesh::BulkData& bulk2,
                    const stk::topology::rank_t rank,
                    const std::string& fieldName)
{
  const stk::mesh::MetaData &meta1 = bulk1.mesh_meta_data();
  const stk::mesh::MetaData &meta2 = bulk2.mesh_meta_data();

  stk::mesh::FieldBase* field1 = meta1.get_field(rank, fieldName);
  stk::mesh::FieldBase* field2 = meta2.get_field(rank, fieldName);

  ASSERT_TRUE(nullptr != field1);
  ASSERT_TRUE(nullptr != field2);

  stk::mesh::EntityVector entityList1;
  stk::mesh::EntityVector entityList2;

  const stk::mesh::Selector s1(*field1);
  const stk::mesh::Selector s2(*field2);

  stk::mesh::get_entities(bulk1, rank, s1, entityList1);
  stk::mesh::get_entities(bulk2, rank, s2, entityList2);

  ASSERT_EQ(entityList1.size(), entityList2.size());

  for(unsigned i=0; i<entityList1.size(); ++i)
  {
    stk::mesh::Entity entity1 = entityList1[i];
    stk::mesh::Entity entity2 = entityList2[i];

    ASSERT_EQ(bulk1.identifier(entity1), bulk2.identifier(entity2));

    double* data1 = static_cast<double*>(stk::mesh::field_data(*field1, entity1) );
    double* data2 = static_cast<double*>(stk::mesh::field_data(*field2, entity2) );

    ASSERT_EQ(*data1, *data2);
  }
}

void verify_field_is_valid(const stk::mesh::MetaData &meta,
                           stk::mesh::Entity entity,
                           double *initialValue,
                           const unsigned fieldLength,
                           const std::string &fieldName)
{
  stk::topology::rank_t rank = meta.mesh_bulk_data().entity_rank(entity);
  stk::mesh::FieldBase* field = meta.get_field(rank, fieldName);
  ASSERT_TRUE(field != nullptr);
  const double* fieldDataForEntity = static_cast<double*>(stk::mesh::field_data(*field, entity) );
  ASSERT_EQ(fieldLength, stk::mesh::field_scalars_per_entity(*field, entity));
  ASSERT_EQ(initialValue[0], fieldDataForEntity[0]);
}

void verify_entity_from_file(stk::mesh::BulkData& input_bulk, stk::mesh::Entity input_entity, stk::mesh::BulkData &bulk, const std::string& bulkEntitySetName)
{
  stk::mesh::Part* bulkEntitySetPart = bulk.mesh_meta_data().get_part(bulkEntitySetName);
  ASSERT_TRUE(bulkEntitySetPart != nullptr);

  stk::topology::rank_t rank = input_bulk.entity_rank(input_entity);
  stk::mesh::EntityVector entities;
  stk::mesh::get_entities(bulk, rank, *bulkEntitySetPart, entities);
  ASSERT_EQ(1u, entities.size());
  EXPECT_EQ(input_bulk.identifier(input_entity), bulk.identifier(entities[0]));
}

void verify_write_and_read_mesh(stk::mesh::BulkData& input_bulk, const std::string &filename, stk::mesh::BulkData &bulkFromFile)
{
  int step = 1;
  double time = 1.0;
  stk::io::write_mesh_with_fields(filename, input_bulk, step, time);

  int numSteps;
  double maxTime;
  stk::io::fill_mesh_save_step_info(filename, bulkFromFile, numSteps, maxTime);
  unlink(filename.c_str());

  ASSERT_EQ(step, numSteps);
  ASSERT_EQ(time, maxTime);
}

void verify_field_exists(stk::mesh::BulkData& bulk,
                         stk::topology::rank_t rank,
                         const std::string& fieldName)
{
  stk::mesh::FieldBase* field = bulk.mesh_meta_data().get_field(rank, fieldName);
  ASSERT_TRUE(field != nullptr);
}

void verify_field_in_file(stk::mesh::BulkData& input_bulk,
                          stk::mesh::Entity input_entity,
                          const std::string& entitySetName,
                          const std::string& fieldName,
                          const std::string &filename)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(input_bulk.parallel());
  stk::mesh::BulkData& bulk = *bulkPtr;
  verify_write_and_read_mesh(input_bulk, filename, bulk);

  stk::topology::rank_t rank = input_bulk.entity_rank(input_entity);
  verify_entity_from_file(input_bulk, input_entity, bulk, entitySetName);
  verify_field_exists(bulk, rank, fieldName);
  compare_fields(input_bulk, bulk, rank, fieldName);
}

void verify_nodesetField_in_file(stk::mesh::BulkData& input_bulk, stk::mesh::Entity node, const std::string& nodesetName, const std::string& fieldName)
{
  std::string filename("nodelist.exo");
  ASSERT_TRUE(input_bulk.entity_rank(node) == stk::topology::NODE_RANK);

  verify_field_in_file(input_bulk, node, nodesetName, fieldName, filename);
}

class MeshWithNodeset : public stk::unit_test_util::MeshFixture
{
};

//BEGINTEST1
TEST_F(MeshWithNodeset, createAndWriteNodesetWithField)
{
  if (stk::parallel_machine_size(get_comm()) == 1)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    std::string nodesetName("nodelist_1");
    stk::mesh::Part& nodesetPart = get_meta().declare_part(nodesetName, stk::topology::NODE_RANK);

    const std::string fieldName = "nodesetField";
    const unsigned fieldLength = 3;
    double initialValue[fieldLength] {0., 0., 0.};
    const int numStates = 1;
    stk::mesh::Field<double> &newField =
        get_meta().declare_field<double>(stk::topology::NODE_RANK, fieldName, numStates);

    stk::mesh::put_field_on_mesh(newField, nodesetPart, fieldLength, initialValue);
    stk::io::set_field_output_type(newField, stk::io::FieldOutputType::VECTOR_3D);

    stk::io::fill_mesh("generated:1x1x1", get_bulk());

    stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);

    get_bulk().modification_begin();
    get_bulk().change_entity_parts(node1, stk::mesh::ConstPartVector{&nodesetPart});
    get_bulk().modification_end();

    stk::io::put_io_part_attribute(nodesetPart);

    verify_field_is_valid(get_meta(), node1, initialValue, fieldLength, fieldName);
    verify_nodesetField_in_file(get_bulk(), node1, nodesetName, fieldName);
  }
}
//ENDTEST1


void verify_sidesetField_in_file(stk::mesh::BulkData& input_bulk, stk::mesh::Entity side, const std::string& sidesetName, const std::string& fieldName)
{
  std::string filename("surface.exo");
  ASSERT_TRUE(input_bulk.entity_rank(side) == input_bulk.mesh_meta_data().side_rank());

  verify_field_in_file(input_bulk, side, sidesetName, fieldName, filename);
}

class MeshWithSideset : public stk::unit_test_util::MeshFixture
{
};

//BEGINTEST2
TEST_F(MeshWithSideset, createAndWriteSidesetWithField)
{
  if (stk::parallel_machine_size(get_comm()) == 1)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    std::string sidesetName("surface_1");
    stk::mesh::Part& sidesetPart = get_meta().declare_part(sidesetName, get_meta().side_rank());

    const std::string fieldName = "sidesetField";
    const unsigned fieldLength = 3;
    double initialValue[fieldLength] {1., 1., 1.};
    const int numStates = 1;
    stk::mesh::Field<double> &newField =
        get_meta().declare_field<double>(get_meta().side_rank(), fieldName, numStates);

    stk::mesh::put_field_on_mesh(newField, sidesetPart, fieldLength, initialValue);
    stk::io::set_field_output_type(newField, stk::io::FieldOutputType::VECTOR_3D);

    stk::io::fill_mesh("generated:1x1x1", get_bulk());

    stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
    unsigned sideOrdinal = 0;

    get_bulk().modification_begin();
    stk::mesh::Entity side = get_bulk().declare_element_side(elem1, sideOrdinal, stk::mesh::PartVector{&sidesetPart});
    get_bulk().modification_end();

    stk::io::put_io_part_attribute(sidesetPart);

    verify_field_is_valid(get_meta(), side, initialValue, fieldLength, fieldName);
    verify_sidesetField_in_file(get_bulk(), side, sidesetName, fieldName);
  }
}
//ENDTEST2

}
