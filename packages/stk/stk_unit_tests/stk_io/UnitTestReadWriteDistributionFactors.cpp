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
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <vector>
#include <string>

namespace {

class DistributionFactor : public stk::unit_test_util::MeshFixture
{
public:
  DistributionFactor()
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }

  void test_df_values(const stk::mesh::BulkData& bulk, const std::vector<std::pair<std::string, std::string>> & dfFieldMapping)
  {
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

    for (const auto & dfFieldInfo : dfFieldMapping) {
      const stk::mesh::Field<double> & dfField = *meta.get_field<double>(stk::topology::FACE_RANK, dfFieldInfo.first);
      const stk::mesh::Part & sidesetPart = *meta.get_part(dfFieldInfo.second);

      const stk::mesh::EntityVector faces = stk::mesh::get_entities(bulk, stk::topology::FACE_RANK, sidesetPart);
      for (stk::mesh::Entity face : faces) {
        double * df = stk::mesh::field_data(dfField, face);
        unsigned fieldLength = stk::mesh::field_scalars_per_entity(dfField, face);
        for (unsigned i = 0; i < fieldLength; ++i) {
          EXPECT_EQ(2.0, df[i]);
        }
      }
    }
  }

  void write_serial_df_mesh(const std::string & meshSpec, const std::string & fileName,
                            const std::vector<std::pair<std::string, std::string>> & dfFieldMapping)
  {
    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
      stk::mesh::MeshBuilder builder(MPI_COMM_SELF);
      builder.set_spatial_dimension(3);
      builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
      auto bulk = builder.create();
      stk::mesh::MetaData & meta = bulk->mesh_meta_data();

      for (const auto & dfFieldInfo : dfFieldMapping) {
        stk::mesh::FieldBase & dfField = meta.declare_field<double>(stk::topology::FACE_RANK, dfFieldInfo.first);
        stk::mesh::Part & sidesetPart = meta.declare_part_with_topology(dfFieldInfo.second, stk::topology::QUAD_4);

        const unsigned numNodes = sidesetPart.topology().num_nodes();
        stk::mesh::put_field_on_mesh(dfField, sidesetPart, numNodes, nullptr);

        stk::io::set_field_role(dfField, Ioss::Field::MESH);
        stk::io::set_distribution_factor_field(sidesetPart, dfField);
      }

      stk::io::fill_mesh(meshSpec, *bulk);

      for (const auto & dfFieldInfo : dfFieldMapping) {
        const stk::mesh::Field<double> & dfField = *meta.get_field<double>(stk::topology::FACE_RANK, dfFieldInfo.first);
        const stk::mesh::Part & sidesetPart = *meta.get_part(dfFieldInfo.second);

        const stk::mesh::EntityVector faces = stk::mesh::get_entities(*bulk, stk::topology::FACE_RANK, sidesetPart);
        for (stk::mesh::Entity face : faces) {
          double * df = stk::mesh::field_data(dfField, face);
          unsigned fieldLength = stk::mesh::field_scalars_per_entity(dfField, face);
          for (unsigned i = 0; i < fieldLength; ++i) {
            df[i] = 2.0;
          }
        }
      }

      stk::io::write_mesh(fileName, *bulk);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
};

std::string get_text_mesh_spec()
{
  const std::string meshSpec = "textmesh:0,1,HEX_8,1,2,5,4,10,11,14,13,block_1\n"
                               "0,2,HEX_8,2,3,6,5,11,12,15,14,block_2\n"
                               "0,3,HEX_8,4,5,8,7,13,14,17,16,block_3\n"
                               "0,4,HEX_8,5,6,9,8,14,15,18,17,block_4\n"
                               "|coordinates: 0,0,0, 1,0,0, 2,0,0, 0,1,0, 1,1,0, 2,1,0, 0,2,0, 1,2,0, 2,2,0, 0,0,1, 1,0,1, 2,0,1, 0,1,1, 1,1,1, 2,1,1, 0,2,1, 1,2,1, 2,2,1"
                               "|sideset:name=internal_surface; data=1,3, 1,2, 2,3, 3,2; split=block";

  return meshSpec;
}

std::string get_text_mesh_spec_2elem()
{
  const std::string meshSpec = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                               "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2\n"
                               "|coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, 0,0,2, 1,0,2, 1,1,2, 0,1,2"
                               "|sideset:name=internal_surface; data=1,6; split=block";

  return meshSpec;
}


TEST_F(DistributionFactor, load_mesh_then_check_DF)
{
  const std::string fileName = "df_mesh.g";
  const std::string meshDesc = get_text_mesh_spec_2elem();

  const std::vector<std::pair<std::string, std::string>> dfFieldMapping = {{"internal_surface_df", "internal_surface"}};
  write_serial_df_mesh(meshDesc, fileName, dfFieldMapping);
  MPI_Barrier(stk::parallel_machine_world());  

  stk::io::StkMeshIoBroker stkIo(stk::parallel_machine_world());
  stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
  stkIo.set_bulk_data(bulkData);
  stkIo.add_mesh_database(fileName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.add_all_mesh_fields_as_input_fields();
  stkIo.populate_bulk_data();

  stk::mesh::BulkData& bulkData = stkIo.bulk_data();

  test_df_values(bulkData, dfFieldMapping);
}

TEST_F(DistributionFactor, load_mesh_then_check_DF_4elem)
{
  const std::string fileName = "df_mesh.g";
  const std::string meshDesc = get_text_mesh_spec();

  const std::vector<std::pair<std::string, std::string>> dfFieldMapping = {{"internal_surface_df", "internal_surface"}};
  write_serial_df_mesh(meshDesc, fileName, dfFieldMapping);
  MPI_Barrier(stk::parallel_machine_world());  

  stk::io::StkMeshIoBroker stkIo(stk::parallel_machine_world());
  stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));
  stkIo.set_bulk_data(bulkData);
  stkIo.add_mesh_database(fileName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.add_all_mesh_fields_as_input_fields();
  stkIo.populate_bulk_data();

  stk::mesh::BulkData& bulkData = stkIo.bulk_data();
  test_df_values(bulkData, dfFieldMapping);
}

TEST_F(DistributionFactor, change_entity_owner_then_check_DF)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) < 3)
  {
    GTEST_SKIP();
  }
  const std::string fileName = "df_mesh.g";
  const std::string meshDesc = get_text_mesh_spec_2elem();

  const std::vector<std::pair<std::string, std::string>> dfFieldMapping = {{"internal_surface_df", "internal_surface"}};
  write_serial_df_mesh(meshDesc, fileName, dfFieldMapping);
  MPI_Barrier(stk::parallel_machine_world());  

  stk::io::StkMeshIoBroker stkIo(stk::parallel_machine_world());
  stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "LINEAR"));
  stkIo.set_bulk_data(bulkData);
  stkIo.add_mesh_database(fileName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.add_all_mesh_fields_as_input_fields();
  stkIo.populate_bulk_data();

  stk::mesh::BulkData& bulkData = stkIo.bulk_data();
  test_df_values(bulkData, dfFieldMapping);

  stk::mesh::EntityProcVec entitiesToChange;
  if (stk::parallel_machine_rank(stk::parallel_machine_world()) == 1)
  {
    stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEMENT_RANK, 2);
    int newOwnerProc = 2;
    entitiesToChange.push_back(std::make_pair(element, newOwnerProc));
    for (int i=9; i <= 12; ++i)
    {
      stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, i);
      entitiesToChange.push_back(std::make_pair(node, newOwnerProc));
    }
  }

  bulkData.change_entity_owner(entitiesToChange);
  test_df_values(bulkData, dfFieldMapping);
}

}
