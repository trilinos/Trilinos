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
#include "UnitTestReadWriteFaces.hpp"
#include "UnitTestReadWriteUtils.hpp"
#include <stk_unit_test_utils/BuildMesh.hpp>

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

class StkFaceIoTestForResultOutput : public StkFaceIoTest
{
public:

  StkFaceIoTestForResultOutput() : StkFaceIoTest() { }

  void setup_mesh_with_face_field(unsigned numBlocks, unsigned numStates = 1)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Field<double>& faceField = create_field<double>(get_meta(), stk::topology::FACE_RANK, faceFieldName, numStates);
    setup_face_mesh(numBlocks);

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);

    initialize_face_field(&faceField, faces);
  }

  void initialize_face_field(stk::mesh::FieldBase* faceField, stk::mesh::EntityVector& faces)
  {
    for(unsigned i = 0; i < faceField->number_of_states(); i++) {
      stk::mesh::FieldState state = stk::mesh::FieldState(i+stk::mesh::StateNP1);
      stk::mesh::FieldBase* field = faceField->field_state(state);
      auto fieldData = field->data<double, stk::mesh::ReadWrite>();
      for(auto face : faces) {
        auto data = fieldData.entity_values(face);

        data() = get_bulk().identifier(face) + i*100;
      }
    }
  }

  void write_mesh(stk::io::DatabasePurpose purpose)
  {
    stk::mesh::FieldBase* faceField = get_meta().get_field(stk::topology::FACE_RANK, faceFieldName);
    stkIoOutput.set_bulk_data(get_bulk());
    size_t outputFileIndex = stkIoOutput.create_output_mesh(fileName, purpose);

    ASSERT_TRUE(faceField != nullptr);
    stkIoOutput.add_field(outputFileIndex, *faceField);
    stkIoOutput.write_output_mesh(outputFileIndex);
    stkIoOutput.begin_output_step(outputFileIndex, 0.0);
    stkIoOutput.write_defined_output_fields(outputFileIndex);
    stkIoOutput.end_output_step(outputFileIndex);
    stkIoOutput.flush_output();
  }

  void test_output_mesh_and_face_field()
  {
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD);

    load_output_mesh(*bulk);
    test_output_mesh(*bulk);
    test_face_field(*bulk);
  }

  void load_output_mesh(stk::mesh::BulkData& bulk) override
  {
    stk::io::fill_mesh_with_fields(fileName, stkIoInput, bulk);
    int numSteps = stkIoInput.get_num_time_steps();
    EXPECT_EQ(1, numSteps);
  }

  virtual void test_face_field(const stk::mesh::BulkData& bulk)
  {
    stk::mesh::FieldBase* field = bulk.mesh_meta_data().get_field(stk::topology::FACE_RANK, faceFieldName);
    ASSERT_TRUE(field != nullptr);

    std::vector<const stk::mesh::FieldBase*> fieldVector{field};

    stk::mesh::communicate_field_data(bulk, fieldVector);

    stk::mesh::EntityVector faces;
    stk::mesh::get_entities(bulk, stk::topology::FACE_RANK, faces);

    auto fieldData = field->data<double>();
    for(stk::mesh::Entity& face : faces) {
      auto data = fieldData.entity_values(face);
      EXPECT_EQ((double)bulk.identifier(face), data());
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
  std::string faceFieldName = "faceField";
};

class StkFaceIoTestForRestart : public StkFaceIoTestForResultOutput
{
public:
  StkFaceIoTestForRestart() : StkFaceIoTestForResultOutput() { }

  void load_output_mesh(stk::mesh::BulkData& bulk) override
  {
    stk::io::fill_mesh_preexisting(stkIoInput, fileName, bulk, stk::io::READ_RESTART);
    int numSteps = stkIoInput.get_num_time_steps();
    EXPECT_EQ(1, numSteps);

    stk::mesh::FieldBase* faceField = bulk.mesh_meta_data().get_field(stk::topology::FACE_RANK, faceFieldName);
    stk::io::MeshField mf(faceField, faceField->name());
    mf.set_read_time(0.0);

    stkIoInput.read_input_field(mf);
    bulk.update_field_data_states(faceField);
  }

  void test_output_mesh_and_face_field(unsigned numStates)
  {
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    create_field<double>(meta, stk::topology::FACE_RANK, faceFieldName, numStates);
    load_output_mesh(*bulk);
    test_output_mesh(*bulk);
    test_face_field(*bulk);
  }

  void test_selected_face_field(const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
  {
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    stk::mesh::FieldBase* field = meta.get_field(stk::topology::FACE_RANK, faceFieldName);
    ASSERT_TRUE(field != nullptr);

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::EDGE_RANK), faces);

    stk::mesh::FieldBase* fieldStateN   = field->field_state(stk::mesh::StateN);
    stk::mesh::FieldBase* fieldStateNM1 = field->field_state(stk::mesh::StateNM1);
    auto fieldStateNData = fieldStateN->data<double>();
    auto fieldStateNM1Data = fieldStateNM1->data<double>();

    for(auto face : faces) {
      auto dataN   = fieldStateNData.entity_values(face);
      auto dataNM1 = fieldStateNM1Data.entity_values(face);

      double expectedDataNValue   = bulk.identifier(face);
      double expectedDataNM1Value = bulk.identifier(face) + 100;
      EXPECT_EQ(expectedDataNValue, dataN());
      EXPECT_EQ(expectedDataNM1Value, dataNM1());
    }
  }

  void test_face_field(const stk::mesh::BulkData& bulk) override
  {
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    test_selected_face_field(bulk, meta.locally_owned_part());

    stk::mesh::FieldBase* field = meta.get_field(stk::topology::FACE_RANK, faceFieldName);
    std::vector<const stk::mesh::FieldBase*> fieldVector;
    for(unsigned i = 0; i < field->number_of_states(); i++) {
      stk::mesh::FieldState state = (stk::mesh::FieldState)i;
      fieldVector.push_back(field->field_state(state));
    }
    stk::mesh::communicate_field_data(bulk, fieldVector);

    test_selected_face_field(bulk, meta.universal_part() & !meta.locally_owned_part());
  }
};

TEST_F(StkFaceIoTestForResultOutput, SerialWriteMeshWithFaceField)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  set_file_name("SerialWriteMeshWithFaceField.exo");

  io_test_utils::ExpectedValues expectedTestValues;
  expectedTestValues.numEdgesPerProc = std::vector<unsigned>{0};
  expectedTestValues.numLocalEdgesPerProc = std::vector<unsigned>{0};
  expectedTestValues.numFacesPerProc = std::vector<unsigned>{1};
  expectedTestValues.numLocalFacesPerProc = std::vector<unsigned>{1};
  expectedTestValues.numConnectedEdges = 0;
  expectedTestValues.globalEdgeCount = 0;
  expectedTestValues.globalElemCount = 2;

  unsigned numBlocks = 2;
  setup_mesh_with_face_field(numBlocks);

  set_expected_values(expectedTestValues);
  test_faces(get_bulk());
  write_mesh(stk::io::WRITE_RESULTS);

  test_output_mesh_and_face_field();
}

TEST_F(StkFaceIoTestForResultOutput, ParallelWriteMeshWithFaceField)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }

  io_test_utils::ExpectedValues expectedTestValues;
  expectedTestValues.numEdgesPerProc = std::vector<unsigned>{0, 0};
  expectedTestValues.numLocalEdgesPerProc = std::vector<unsigned>{0, 0};
  expectedTestValues.numFacesPerProc = std::vector<unsigned>{1, 1};
  expectedTestValues.numLocalFacesPerProc = std::vector<unsigned>{1, 0};
  expectedTestValues.numConnectedEdges = 0;
  expectedTestValues.globalEdgeCount = 0;
  expectedTestValues.globalElemCount = 2;

  unsigned numBlocks = 2;
  setup_mesh_with_face_field(numBlocks);

  set_expected_values(expectedTestValues);
  test_faces(get_bulk());
  write_mesh(stk::io::WRITE_RESULTS);

  test_output_mesh_and_face_field();
}

TEST_F(StkFaceIoTestForRestart, SerialWriteMeshWithFaceField)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  io_test_utils::ExpectedValues expectedTestValues;
  expectedTestValues.numEdgesPerProc = std::vector<unsigned>{0};
  expectedTestValues.numLocalEdgesPerProc = std::vector<unsigned>{0};
  expectedTestValues.numFacesPerProc = std::vector<unsigned>{1};
  expectedTestValues.numLocalFacesPerProc = std::vector<unsigned>{1};
  expectedTestValues.numConnectedEdges = 0;
  expectedTestValues.globalEdgeCount = 0;
  expectedTestValues.globalElemCount = 2;

  unsigned numStates = 3;
  unsigned numBlocks = 2;

  setup_mesh_with_face_field(numBlocks, numStates);

  set_expected_values(expectedTestValues);
  test_faces(get_bulk());
  write_mesh(stk::io::WRITE_RESTART);

  test_output_mesh_and_face_field(numStates);
}

TEST_F(StkFaceIoTestForRestart, ParallelWriteMeshWithFaceField)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { return; }

  io_test_utils::ExpectedValues expectedTestValues;
  expectedTestValues.numEdgesPerProc = std::vector<unsigned>{0, 0, 0};
  expectedTestValues.numLocalEdgesPerProc = std::vector<unsigned>{0, 0, 0};
  expectedTestValues.numFacesPerProc = std::vector<unsigned>{1, 2, 1};
  expectedTestValues.numLocalFacesPerProc = std::vector<unsigned>{1, 1, 0};
  expectedTestValues.numConnectedEdges = 0;
  expectedTestValues.globalEdgeCount = 0;
  expectedTestValues.globalElemCount = 3;

  unsigned numStates = 3;
  unsigned numBlocks = 3;

  setup_mesh_with_face_field(numBlocks, numStates);

  set_expected_values(expectedTestValues);
  test_faces(get_bulk());
  write_mesh(stk::io::WRITE_RESTART);

  test_output_mesh_and_face_field(numStates);
}
