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

#include <stk_io/FillMesh.hpp>
#include "UnitTestReadWriteEdges.hpp"
#include "UnitTestReadWriteUtils.hpp"
#include <stk_unit_test_utils/BuildMesh.hpp>

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

class StkEdgeIoTestForResultOutput : public StkEdgeIoTest
{
public:

  StkEdgeIoTestForResultOutput() : StkEdgeIoTest() { }

  void setup_mesh_with_edge_field(unsigned numBlocks, unsigned numStates = 1)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Field<double>& edgeField = create_field<double>(get_meta(), stk::topology::EDGE_RANK, edgeFieldName, numStates);
    setup_edge_mesh(numBlocks);

    stk::mesh::EntityVector edges;
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::EDGE_RANK), edges);

    initialize_edge_field(&edgeField, edges);
  }

  void initialize_edge_field(stk::mesh::FieldBase* edgeField, stk::mesh::EntityVector& edges)
  {
    for(auto edge : edges) {
      for(unsigned i = 0; i < edgeField->number_of_states(); i++) {
        stk::mesh::FieldState state = stk::mesh::FieldState(i+stk::mesh::StateNP1);
        stk::mesh::FieldBase* field = edgeField->field_state(state);
        double* data = reinterpret_cast<double*>(stk::mesh::field_data(*field, edge));

        *data = get_bulk().identifier(edge) + i*100;
      }
    }
  }

  void write_mesh(stk::io::DatabasePurpose purpose)
  {
    stk::mesh::FieldBase* edgeField = get_meta().get_field(stk::topology::EDGE_RANK, edgeFieldName);
    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(get_bulk());
    size_t outputFileIndex = stkIo.create_output_mesh(fileName, purpose);

    ASSERT_TRUE(edgeField != nullptr);
    stkIo.add_field(outputFileIndex, *edgeField);
    stkIo.write_output_mesh(outputFileIndex);
    stkIo.begin_output_step(outputFileIndex, 0.0);
    stkIo.write_defined_output_fields(outputFileIndex);
    stkIo.end_output_step(outputFileIndex);
  }

  void test_output_mesh_and_edge_field()
  {
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

    load_output_mesh(*bulk);
    test_output_mesh(*bulk);
    test_edge_field(*bulk);
  }

  void load_output_mesh(stk::mesh::BulkData& bulk) override
  {
    int numSteps;
    double maxTime;
    stk::io::fill_mesh_save_step_info(fileName, bulk, numSteps, maxTime);
    EXPECT_EQ(1, numSteps);
  }

  virtual void test_edge_field(const stk::mesh::BulkData& bulk)
  {
    stk::mesh::FieldBase* field = bulk.mesh_meta_data().get_field(stk::topology::EDGE_RANK, edgeFieldName);
    ASSERT_TRUE(field != nullptr);

    std::vector<const stk::mesh::FieldBase*> fieldVector{field};

    stk::mesh::communicate_field_data(bulk, fieldVector);

    stk::mesh::EntityVector edges;
    stk::mesh::get_entities(bulk, stk::topology::EDGE_RANK, edges);

    for(stk::mesh::Entity& edge : edges) {
      double* data = reinterpret_cast<double*>(stk::mesh::field_data(*field, edge));
      EXPECT_EQ((double)bulk.identifier(edge), *data);
    }
  }

  template <typename T>
  stk::mesh::Field<T> & create_field(stk::mesh::MetaData& meta, stk::topology::rank_t rank, const std::string & name,
                                     unsigned numStates = 1, unsigned numComponent = 1)
  {
    const std::vector<T> init(numComponent, 1);
    stk::mesh::Field<T> & field = meta.declare_field<T>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), numComponent, init.data());
    return field;
  }

protected:
  std::string edgeFieldName = "edgeField";
};

class StkEdgeIoTestForRestart : public StkEdgeIoTestForResultOutput
{
public:
  StkEdgeIoTestForRestart() : StkEdgeIoTestForResultOutput() { }

  void load_output_mesh(stk::mesh::BulkData& bulk) override
  {
    stk::io::StkMeshIoBroker stkIo;
    stk::io::fill_mesh_preexisting(stkIo, fileName, bulk, stk::io::READ_RESTART);
    int numSteps = stkIo.get_num_time_steps();
    EXPECT_EQ(1, numSteps);

    stk::mesh::FieldBase* edgeField = bulk.mesh_meta_data().get_field(stk::topology::EDGE_RANK, edgeFieldName);
    stk::io::MeshField mf(edgeField, edgeField->name());
    mf.set_read_time(0.0);

    stkIo.read_input_field(mf);
    bulk.update_field_data_states(edgeField);
  }

  void test_output_mesh_and_edge_field(unsigned numStates)
  {
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    create_field<double>(meta, stk::topology::EDGE_RANK, edgeFieldName, numStates);
    load_output_mesh(*bulk);
    test_output_mesh(*bulk);
    test_edge_field(*bulk);
  }

  void test_selected_edge_field(const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
  {
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    stk::mesh::FieldBase* field = meta.get_field(stk::topology::EDGE_RANK, edgeFieldName);
    ASSERT_TRUE(field != nullptr);

    stk::mesh::EntityVector edges;
    stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::EDGE_RANK), edges);

    stk::mesh::FieldBase* fieldStateN   = field->field_state(stk::mesh::StateN);
    stk::mesh::FieldBase* fieldStateNM1 = field->field_state(stk::mesh::StateNM1);

    for(auto edge : edges) {
      double* dataN   = reinterpret_cast<double*>(stk::mesh::field_data(*fieldStateN, edge));
      double* dataNM1 = reinterpret_cast<double*>(stk::mesh::field_data(*fieldStateNM1, edge));

      double expectedDataNValue   = bulk.identifier(edge);
      double expectedDataNM1Value = bulk.identifier(edge) + 100;
      EXPECT_EQ(expectedDataNValue, *dataN);
      EXPECT_EQ(expectedDataNM1Value, *dataNM1);
    }
  }

  void test_edge_field(const stk::mesh::BulkData& bulk) override
  {
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    test_selected_edge_field(bulk, meta.locally_owned_part());

    stk::mesh::FieldBase* field = meta.get_field(stk::topology::EDGE_RANK, edgeFieldName);
    std::vector<const stk::mesh::FieldBase*> fieldVector;
    for(unsigned i = 0; i < field->number_of_states(); i++) {
      stk::mesh::FieldState state = (stk::mesh::FieldState)i;
      fieldVector.push_back(field->field_state(state));
    }
    stk::mesh::communicate_field_data(bulk, fieldVector);

    test_selected_edge_field(bulk, meta.universal_part() & !meta.locally_owned_part());
  }
};

TEST_F(StkEdgeIoTestForResultOutput, SerialWriteMeshWithEdgeField)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  io_test_utils::ExpectedValues expectedValues;
  expectedValues.numEdgesPerProc = std::vector<unsigned>{12};
  expectedValues.numLocalEdgesPerProc = std::vector<unsigned>{12};
  expectedValues.numFacesPerProc = std::vector<unsigned>{0};
  expectedValues.numLocalFacesPerProc = std::vector<unsigned>{0};
  expectedValues.numConnectedEdges = 0;
  expectedValues.globalEdgeCount = 12;
  expectedValues.globalElemCount = 1;

  setup_mesh_with_edge_field(1);

  set_expected_values(expectedValues);
  test_edges(get_bulk());
  write_mesh(stk::io::WRITE_RESULTS);

  test_output_mesh_and_edge_field();
}

TEST_F(StkEdgeIoTestForResultOutput, ParallelWriteMeshWithEdgeField)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }

  io_test_utils::ExpectedValues expectedValues;
  expectedValues.numEdgesPerProc = std::vector<unsigned>{12, 12};
  expectedValues.numLocalEdgesPerProc = std::vector<unsigned>{12, 8};
  expectedValues.numFacesPerProc = std::vector<unsigned>{0, 0};
  expectedValues.numLocalFacesPerProc = std::vector<unsigned>{0, 0};
  expectedValues.numConnectedEdges = 0;
  expectedValues.globalEdgeCount = 20;
  expectedValues.globalElemCount = 2;

  setup_mesh_with_edge_field(2);

  set_expected_values(expectedValues);
  test_edges(get_bulk());
  write_mesh(stk::io::WRITE_RESULTS);

  test_output_mesh_and_edge_field();
}

TEST_F(StkEdgeIoTestForRestart, SerialWriteMeshWithEdgeField)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  io_test_utils::ExpectedValues expectedValues;
  expectedValues.numEdgesPerProc = std::vector<unsigned>{12};
  expectedValues.numLocalEdgesPerProc = std::vector<unsigned>{12};
  expectedValues.numFacesPerProc = std::vector<unsigned>{0};
  expectedValues.numLocalFacesPerProc = std::vector<unsigned>{0};
  expectedValues.numConnectedEdges = 0;
  expectedValues.globalEdgeCount = 12;
  expectedValues.globalElemCount = 1;

  unsigned numStates = 3;

  setup_mesh_with_edge_field(1, numStates);

  set_expected_values(expectedValues);
  test_edges(get_bulk());
  write_mesh(stk::io::WRITE_RESTART);

  test_output_mesh_and_edge_field(numStates);
}

TEST_F(StkEdgeIoTestForRestart, ParallelWriteMeshWithEdgeField)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { return; }

  io_test_utils::ExpectedValues expectedValues;
  expectedValues.numEdgesPerProc = std::vector<unsigned>{12, 12, 12};
  expectedValues.numLocalEdgesPerProc = std::vector<unsigned>{12, 8, 8};
  expectedValues.numFacesPerProc = std::vector<unsigned>{0, 0, 0};
  expectedValues.numLocalFacesPerProc = std::vector<unsigned>{0, 0, 0};
  expectedValues.numConnectedEdges = 0;
  expectedValues.globalEdgeCount = 28;
  expectedValues.globalElemCount = 3;

  unsigned numStates = 3;

  setup_mesh_with_edge_field(3, numStates);

  set_expected_values(expectedValues);
  test_edges(get_bulk());
  write_mesh(stk::io::WRITE_RESTART);

  test_output_mesh_and_edge_field(numStates);
}
