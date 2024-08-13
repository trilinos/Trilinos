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
#include "gtest/gtest.h"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_balance/internal/Diagnostics.hpp"
#include "stk_balance/internal/DiagnosticsContainer.hpp"
#include "stk_balance/internal/DiagnosticsPrinter.hpp"
#include "stk_balance/internal/privateDeclarations.hpp"
#include "stk_balance/rebalance.hpp"
#include "stk_balance/balance.hpp"
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_unit_test_utils/TextMeshToFile.hpp>

class TestDiagnosticsComputation : public stk::unit_test_util::MeshFixture
{
protected:
  TestDiagnosticsComputation()
  {
  }

  virtual void NGPSetUp() override {
    stk::balance::impl::g_diagnosticsContainer.clear();
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
    testing::internal::CaptureStderr();
  }

  virtual void NGPTearDown() override {
    stk::balance::impl::g_diagnosticsContainer.clear();
    stk::EnvData::instance().m_outputP0 = &std::cout;
    testing::internal::GetCapturedStderr();
  }

  std::string get_output_file_name() {
    return "tempOutputMesh.g";
  }

  void fill_rebalance_settings(stk::balance::BalanceSettings & balanceSettings, const std::string & decompMethod,
                               int numOutputProcessors)
  {
    balanceSettings.set_is_rebalancing(true);
    balanceSettings.set_output_filename(get_output_file_name());
    balanceSettings.set_num_input_processors(get_parallel_size());
    balanceSettings.set_num_output_processors(numOutputProcessors);
    balanceSettings.setDecompMethod(decompMethod);
    balanceSettings.setShouldPrintDiagnostics(true);
  }

  std::vector<const stk::mesh::Field<double>*> create_multi_criteria_fields(int numCriteria) {
    get_meta().enable_late_fields();

    std::vector<const stk::mesh::Field<double>*> multiCriteriaFields;
    for (int i = 0; i < numCriteria; ++i) {
      stk::mesh::Field<double> & field = get_meta().declare_field<double>(stk::topology::ELEM_RANK,
                                                                          "test_criteria_field_" + std::to_string(i));
      stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
      multiCriteriaFields.push_back(&field);
    }

    return multiCriteriaFields;
  }

  void fill_multi_criteria_fields(std::vector<const stk::mesh::Field<double>*> & multiCriteriaFields) {
    for (unsigned i = 0; i < multiCriteriaFields.size(); ++i) {
      const stk::mesh::Field<double> & field = *multiCriteriaFields[i];

      stk::mesh::EntityVector elems;
      stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elems);

      for (const stk::mesh::Entity & elem : elems) {
        double * weight = stk::mesh::field_data(field, elem);
        *weight = get_bulk().identifier(elem) * (i+1);
      }
    }
  }

  void build_mesh(const std::string & meshDesc) {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
    MPI_Barrier(get_comm());
  }

  void build_mesh_file(const std::string & meshDesc, const stk::balance::BalanceSettings & balanceSettings,
                       stk::io::StkMeshIoBroker & ioBroker) {
    const std::string tempFileName = "tempFile.g";

    stk::unit_test_util::TextMeshToFile tMesh(get_comm(), stk::mesh::BulkData::AUTO_AURA);
    tMesh.setup_mesh(meshDesc, tempFileName);
    tMesh.write_mesh();

    MPI_Barrier(get_comm());

    get_meta().set_coordinate_field_name(balanceSettings.getCoordinateFieldName());
    stk::io::fill_mesh_preexisting(ioBroker, tempFileName, get_bulk());
  }

  void rebalanceMesh(stk::io::StkMeshIoBroker & ioBroker, const stk::balance::BalanceSettings & balanceSettings)
  {
    stk::mesh::BulkData & bulk = ioBroker.bulk_data();

    stk::balance::internal::register_internal_fields_and_parts(bulk, balanceSettings);
    stk::balance::set_up_diagnostics(balanceSettings);

    stk::balance::rebalance(ioBroker, balanceSettings);

    stk::balance::DiagnosticsPrinter diagPrinter(get_comm(), balanceSettings.get_num_output_processors());
    diagPrinter.print(sierra::Env::outputP0());
  }

  std::string mesh_desc_four_beams() {
    const std::string meshDesc = "0,1,BEAM_2,1,2\n"
                                 "0,2,BEAM_2,2,3\n"
                                 "0,3,BEAM_2,3,4\n"
                                 "0,4,BEAM_2,3,5";

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 2,0,0, 3,-1,0, 3,1,0
    };

    return stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates);
  }

  std::string mesh_desc_four_shells_in_square() {
    const std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,5,4\n"
                                 "0,2,SHELL_QUAD_4,2,3,6,5\n"
                                 "0,3,SHELL_QUAD_4,4,5,8,7\n"
                                 "0,4,SHELL_QUAD_4,5,6,9,8";

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 2,0,0,
      0,1,0, 1,1,0, 2,1,0,
      0,2,0, 1,2,0, 2,2,0
    };

    return stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates);
  }

  std::string mesh_desc_three_hex_in_row() {
    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                                 "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                                 "0,3,HEX_8,9,10,11,12,13,14,15,16";

    std::vector<double> coordinates = {
      0,0,0, 0,1,0, 0,1,1, 0,0,1,
      1,0,0, 1,1,0, 1,1,1, 1,0,1,
      2,0,0, 2,1,0, 2,1,1, 2,0,1,
      3,0,0, 3,1,0, 3,1,1, 3,0,1
    };

    return stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates);
  }

  std::string mesh_desc_four_hex_in_square() {
    const std::string meshDesc = "0,1,HEX_8,1,2,5,4,7,8,11,10\n"
                                 "0,2,HEX_8,2,3,6,5,8,9,12,11\n"
                                 "0,3,HEX_8,7,8,11,10,13,14,17,16\n"
                                 "0,4,HEX_8,8,9,12,11,14,15,18,17";

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 2,0,0,
      0,1,0, 1,1,0, 2,1,0,
      0,0,1, 1,0,1, 2,0,1,
      0,1,1, 1,1,1, 2,1,1,
      0,0,2, 1,0,2, 2,0,2,
      0,1,2, 1,1,2, 2,1,2
    };

    return stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates);
  }

  std::string mesh_desc_hex_pyramid_tet() {
    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                                 "0,2,PYRAMID_5,5,6,7,8,10\n"
                                 "0,3,TET_4,5,9,8,10\n"
                                 "0,4,TET_4,8,9,12,10\n"
                                 "0,5,TET_4,8,12,10,11\n"
                                 "0,6,TET_4,7,8,10,11";

    std::vector<double> coordinates = {
      0,0,0, 0,1,0, 0,1,1, 0,0,1,
      1,0,0, 1,1,0, 1,1,1, 1,0,1,
      2,0,0, 2,1,0, 2,1,1, 2,0,1
    };

    return stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates);
  }

  template <typename T, typename DT>
  void test_diag_values(const std::vector<DT> & expectedValues) {
    if (get_parallel_rank() == 0) {
      T * diag = stk::balance::get_diagnostic<T>();
      ASSERT_NE(diag, nullptr);

      for (unsigned rank = 0; rank < expectedValues.size(); ++rank) {
        EXPECT_DOUBLE_EQ(diag->get_rank_value(rank), expectedValues[rank]) << "rank=" << rank;
      }
    }
  }

  template <typename T, typename DT>
  void test_diag_multi_values(unsigned column, const std::vector<DT> & expectedValues) {
    if (get_parallel_rank() == 0) {
      T * diag = stk::balance::get_diagnostic<T>();
      ASSERT_NE(diag, nullptr);

      for (unsigned rank = 0; rank < expectedValues.size(); ++rank) {
        EXPECT_DOUBLE_EQ(diag->get_rank_value(column, rank), expectedValues[rank]) << "column=" << column << ", rank=" << rank;
      }
    }
  }
};

TEST_F(TestDiagnosticsComputation, ElementCount_Balance_HexMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_three_hex_in_row());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<unsigned> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {3}; }
  else if (get_parallel_size() == 2) { expectedValues = {2, 1}; }
  else if (get_parallel_size() == 3) { expectedValues = {1, 1, 1}; }
  else if (get_parallel_size() == 4) { expectedValues = {1, 1, 1, 0}; }

  test_diag_values<stk::balance::ElementCountDiagnostic, unsigned>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, ElementCount_Balance_HexMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_three_hex_in_row());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<unsigned> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {3}; }
  else if (get_parallel_size() == 2) { expectedValues = {2, 1}; }
  else if (get_parallel_size() == 3) { expectedValues = {1, 1, 1}; }
  else if (get_parallel_size() == 4) { expectedValues = {1, 1, 0, 1}; }

  test_diag_values<stk::balance::ElementCountDiagnostic, unsigned>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, ElementCount_Balance_HexPyramidTetMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_hex_pyramid_tet());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<unsigned> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {6}; }
  else if (get_parallel_size() == 2) { expectedValues = {3, 3}; }
  else if (get_parallel_size() == 3) { expectedValues = {1, 3, 2}; }
  else if (get_parallel_size() == 4) { expectedValues = {1, 2, 2, 1}; }

  test_diag_values<stk::balance::ElementCountDiagnostic, unsigned>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, ElementCount_Balance_HexPyramidTetMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_hex_pyramid_tet());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<unsigned> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {6}; }
  else if (get_parallel_size() == 2) { expectedValues = {3, 3}; }
  else if (get_parallel_size() == 3) { expectedValues = {1, 1, 4}; }
  else if (get_parallel_size() == 4) { expectedValues = {1, 2, 1, 2}; }

  test_diag_values<stk::balance::ElementCountDiagnostic, unsigned>(expectedValues);
}

class RebalanceNumOutputProcs : public TestDiagnosticsComputation,
                                public testing::WithParamInterface<int> {};

TEST_P(RebalanceNumOutputProcs, ElementCount_Rebalance_HexMesh_GeometricPartitioner)
{
  if (get_parallel_size() != 2) return;

  stk::balance::GraphCreationSettings balanceSettings;
  fill_rebalance_settings(balanceSettings, "rcb", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  build_mesh_file(mesh_desc_three_hex_in_row(), balanceSettings, ioBroker);

  rebalanceMesh(ioBroker, balanceSettings);

  std::vector<unsigned> expectedValues;
  if      (GetParam() == 1) { expectedValues = {3}; }
  else if (GetParam() == 2) { expectedValues = {2, 1}; }
  else if (GetParam() == 3) { expectedValues = {1, 1, 1}; }
  else if (GetParam() == 4) { expectedValues = {1, 1, 1, 0}; }

  test_diag_values<stk::balance::ElementCountDiagnostic, unsigned>(expectedValues);
}

TEST_P(RebalanceNumOutputProcs, ElementCount_Rebalance_HexMesh_GraphPartitioner)
{
  if (get_parallel_size() != 2) return;

  stk::balance::GraphCreationSettings balanceSettings;
  fill_rebalance_settings(balanceSettings, "parmetis", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  build_mesh_file(mesh_desc_three_hex_in_row(), balanceSettings, ioBroker);

  rebalanceMesh(ioBroker, balanceSettings);

  std::vector<unsigned> expectedValues;
  if      (GetParam() == 1) { expectedValues = {3}; }
  else if (GetParam() == 2) { expectedValues = {2, 1}; }
  else if (GetParam() == 3) { expectedValues = {1, 1, 1}; }
  else if (GetParam() == 4) { expectedValues = {1, 1, 0, 1}; }

  test_diag_values<stk::balance::ElementCountDiagnostic, unsigned>(expectedValues);
}

TEST_P(RebalanceNumOutputProcs, ElementCount_Rebalance_HexPyramidTetMesh_GeometricPartitioner)
{
  if (get_parallel_size() != 2) return;

  stk::balance::GraphCreationSettings balanceSettings;
  fill_rebalance_settings(balanceSettings, "rcb", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  build_mesh_file(mesh_desc_hex_pyramid_tet(), balanceSettings, ioBroker);

  rebalanceMesh(ioBroker, balanceSettings);

  std::vector<unsigned> expectedValues;
  if      (GetParam() == 1) { expectedValues = {6}; }
  else if (GetParam() == 2) { expectedValues = {3, 3}; }
  else if (GetParam() == 3) { expectedValues = {1, 4, 1}; }
  else if (GetParam() == 4) { expectedValues = {1, 2, 2, 1}; }

  test_diag_values<stk::balance::ElementCountDiagnostic, unsigned>(expectedValues);
}

TEST_P(RebalanceNumOutputProcs, ElementCount_Rebalance_HexPyramidTetMesh_GraphPartitioner)
{
  if (get_parallel_size() != 2) return;

  stk::balance::GraphCreationSettings balanceSettings;
  fill_rebalance_settings(balanceSettings, "parmetis", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  build_mesh_file(mesh_desc_hex_pyramid_tet(), balanceSettings, ioBroker);

  rebalanceMesh(ioBroker, balanceSettings);

  std::vector<unsigned> expectedValues;
  if      (GetParam() == 1) { expectedValues = {6}; }
  else if (GetParam() == 2) { expectedValues = {3, 3}; }
  else if (GetParam() == 3) { expectedValues = {1, 1, 4}; }
  else if (GetParam() == 4) { expectedValues = {1, 2, 1, 2}; }

  test_diag_values<stk::balance::ElementCountDiagnostic, unsigned>(expectedValues);
}


TEST_F(TestDiagnosticsComputation, TotalElementWeight_Balance_HexMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_three_hex_in_row());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<unsigned> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0}; }
  else if (get_parallel_size() == 2) { expectedValues = {2, 1}; }
  else if (get_parallel_size() == 3) { expectedValues = {1, 1, 1}; }
  else if (get_parallel_size() == 4) { expectedValues = {1, 1, 1, 0}; }

  test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(0, expectedValues);
}

TEST_F(TestDiagnosticsComputation, TotalElementWeight_Balance_HexMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_three_hex_in_row());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<unsigned> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0}; }
  else if (get_parallel_size() == 2) { expectedValues = {2, 1}; }
  else if (get_parallel_size() == 3) { expectedValues = {1, 1, 1}; }
  else if (get_parallel_size() == 4) { expectedValues = {1, 1, 0, 1}; }

  test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(0, expectedValues);
}

TEST_F(TestDiagnosticsComputation, TotalElementWeight_Balance_HexPyramidTetMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_hex_pyramid_tet());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<unsigned> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0}; }
  else if (get_parallel_size() == 2) { expectedValues = {3, 3}; }
  else if (get_parallel_size() == 3) { expectedValues = {1, 3, 2}; }
  else if (get_parallel_size() == 4) { expectedValues = {1, 2, 2, 1}; }

  test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(0, expectedValues);
}

TEST_F(TestDiagnosticsComputation, TotalElementWeight_Balance_HexPyramidTetMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_hex_pyramid_tet());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<unsigned> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0}; }
  else if (get_parallel_size() == 2) { expectedValues = {3, 3}; }
  else if (get_parallel_size() == 3) { expectedValues = {1, 1, 4}; }
  else if (get_parallel_size() == 4) { expectedValues = {1, 2, 1, 2}; }

  test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(0, expectedValues);
}


TEST_P(RebalanceNumOutputProcs, TotalElementWeight_Rebalance_HexMesh_GeometricPartitioner)
{
  if (get_parallel_size() != 2) return;

  stk::balance::GraphCreationSettings balanceSettings;
  fill_rebalance_settings(balanceSettings, "rcb", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  build_mesh_file(mesh_desc_three_hex_in_row(), balanceSettings, ioBroker);

  rebalanceMesh(ioBroker, balanceSettings);

  std::vector<unsigned> expectedValues;
  if      (GetParam() == 1) { expectedValues = {0}; }
  else if (GetParam() == 2) { expectedValues = {2, 1}; }
  else if (GetParam() == 3) { expectedValues = {1, 1, 1}; }
  else if (GetParam() == 4) { expectedValues = {1, 1, 1, 0}; }

  test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(0, expectedValues);
}

TEST_P(RebalanceNumOutputProcs, TotalElementWeight_Rebalance_HexMesh_GraphPartitioner)
{
  if (get_parallel_size() != 2) return;

  stk::balance::GraphCreationSettings balanceSettings;
  fill_rebalance_settings(balanceSettings, "parmetis", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  build_mesh_file(mesh_desc_three_hex_in_row(), balanceSettings, ioBroker);

  rebalanceMesh(ioBroker, balanceSettings);

  std::vector<unsigned> expectedValues;
  if      (GetParam() == 1) { expectedValues = {0}; }
  else if (GetParam() == 2) { expectedValues = {2, 1}; }
  else if (GetParam() == 3) { expectedValues = {1, 1, 1}; }
  else if (GetParam() == 4) { expectedValues = {1, 1, 0, 1}; }

  test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(0, expectedValues);
}

TEST_P(RebalanceNumOutputProcs, TotalElementWeight_Rebalance_HexPyramidTetMesh_GeometricPartitioner)
{
  if (get_parallel_size() != 2) return;

  stk::balance::GraphCreationSettings balanceSettings;
  fill_rebalance_settings(balanceSettings, "rcb", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  build_mesh_file(mesh_desc_hex_pyramid_tet(), balanceSettings, ioBroker);

  rebalanceMesh(ioBroker, balanceSettings);

  std::vector<unsigned> expectedValues;
  if      (GetParam() == 1) { expectedValues = {0}; }
  else if (GetParam() == 2) { expectedValues = {3, 3}; }
  else if (GetParam() == 3) { expectedValues = {1, 4, 1}; }
  else if (GetParam() == 4) { expectedValues = {1, 2, 2, 1}; }

  test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(0, expectedValues);
}

TEST_P(RebalanceNumOutputProcs, TotalElementWeight_Rebalance_HexPyramidTetMesh_GraphPartitioner)
{
  if (get_parallel_size() != 2) return;

  stk::balance::GraphCreationSettings balanceSettings;
  fill_rebalance_settings(balanceSettings, "parmetis", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  build_mesh_file(mesh_desc_hex_pyramid_tet(), balanceSettings, ioBroker);

  rebalanceMesh(ioBroker, balanceSettings);

  std::vector<unsigned> expectedValues;
  if      (GetParam() == 1) { expectedValues = {0}; }
  else if (GetParam() == 2) { expectedValues = {3, 3}; }
  else if (GetParam() == 3) { expectedValues = {1, 1, 4}; }
  else if (GetParam() == 4) { expectedValues = {1, 2, 1, 2}; }

  test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(0, expectedValues);
}


TEST_F(TestDiagnosticsComputation, TotalElementWeight_MultiCriteria_Balance_HexMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_three_hex_in_row());

  const int numCriteria = 2;
  std::vector<const stk::mesh::Field<double>*> multiCriteriaFields = create_multi_criteria_fields(numCriteria);
  fill_multi_criteria_fields(multiCriteriaFields);

  stk::balance::MultipleCriteriaSettings balanceSettings(multiCriteriaFields);
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  for (int column = 0; column < numCriteria; ++column) {
    const unsigned elem1Weight = 1 * (column+1);
    const unsigned elem2Weight = 2 * (column+1);
    const unsigned elem3Weight = 3 * (column+1);
    std::vector<unsigned> expectedValues;
    if      (get_parallel_size() == 1) { expectedValues = {0}; }
    else if (get_parallel_size() == 2) { expectedValues = {elem1Weight+elem2Weight, elem3Weight}; }
    else if (get_parallel_size() == 3) { expectedValues = {elem1Weight,             elem2Weight, elem3Weight}; }
    else if (get_parallel_size() == 4) { expectedValues = {elem1Weight,             elem2Weight, elem3Weight, 0}; }

    test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(column, expectedValues);
  }
}

TEST_F(TestDiagnosticsComputation, TotalElementWeight_MultiCriteria_Balance_HexMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_three_hex_in_row());

  const int numCriteria = 2;
  std::vector<const stk::mesh::Field<double>*> multiCriteriaFields = create_multi_criteria_fields(numCriteria);
  fill_multi_criteria_fields(multiCriteriaFields);

  stk::balance::MultipleCriteriaSettings balanceSettings(multiCriteriaFields);
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  for (int column = 0; column < numCriteria; ++column) {
    const unsigned elem1Weight = 1 * (column+1);
    const unsigned elem2Weight = 2 * (column+1);
    const unsigned elem3Weight = 3 * (column+1);
    std::vector<unsigned> expectedValues;
    if      (get_parallel_size() == 1) { expectedValues = {0}; }
    else if (get_parallel_size() == 2) { expectedValues = {elem1Weight+elem2Weight, elem3Weight}; }
    else if (get_parallel_size() == 3) {
      unsigned numElemsOnLastProc = stk::mesh::count_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part());
      MPI_Bcast(&numElemsOnLastProc, 1, MPI_UNSIGNED, 2, get_comm());
      if (numElemsOnLastProc > 0) {
        expectedValues = {elem3Weight, 0, elem1Weight+elem2Weight};
      }
      else {
        expectedValues = {elem3Weight, elem1Weight+elem2Weight, 0};
      }
    }
    else if (get_parallel_size() == 4) { expectedValues = {elem1Weight, elem2Weight, elem3Weight, 0}; }

    test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(column, expectedValues);
  }
}

TEST_P(RebalanceNumOutputProcs, TotalElementWeight_MultiCriteria_Rebalance_HexMesh_GeometricPartitioner)
{
  if (get_parallel_size() != 2) return;
  if (GetParam() == 4) return;  // RCB with multi-criteria weights with an empty output subdomain hangs

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const int numCriteria = 2;
  std::vector<const stk::mesh::Field<double>*> multiCriteriaFields = create_multi_criteria_fields(numCriteria);

  stk::balance::MultipleCriteriaSettings balanceSettings(multiCriteriaFields);
  fill_rebalance_settings(balanceSettings, "rcb", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  build_mesh_file(mesh_desc_three_hex_in_row(), balanceSettings, ioBroker);
  fill_multi_criteria_fields(multiCriteriaFields);

  rebalanceMesh(ioBroker, balanceSettings);

  //for (int column = 0; column < numCriteria; ++column) {
  for (int column = 0; column < 1; ++column) {
    const unsigned elem1Weight = 1 * (column+1);
    const unsigned elem2Weight = 2 * (column+1);
    const unsigned elem3Weight = 3 * (column+1);
    std::vector<unsigned> expectedValues;
    if      (GetParam() == 1) { expectedValues = {0}; }
    else if (GetParam() == 2) { expectedValues = {elem1Weight+elem2Weight, elem3Weight}; }
    else if (GetParam() == 3) { expectedValues = {elem1Weight,             elem2Weight, elem3Weight}; }
    else if (GetParam() == 4) { expectedValues = {elem1Weight,             elem2Weight, elem3Weight, 0}; }

    test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(column, expectedValues);
  }
}

TEST_P(RebalanceNumOutputProcs, TotalElementWeight_MultiCriteria_Rebalance_HexMesh_GraphPartitioner)
{
  if (get_parallel_size() != 2) return;

  const int numCriteria = 2;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  std::vector<const stk::mesh::Field<double>*> multiCriteriaFields = create_multi_criteria_fields(numCriteria);

  stk::balance::MultipleCriteriaSettings balanceSettings(multiCriteriaFields);
  fill_rebalance_settings(balanceSettings, "parmetis", GetParam());

  stk::io::StkMeshIoBroker ioBroker;
  build_mesh_file(mesh_desc_three_hex_in_row(), balanceSettings, ioBroker);
  fill_multi_criteria_fields(multiCriteriaFields);

  rebalanceMesh(ioBroker, balanceSettings);

  for (int column = 0; column < numCriteria; ++column) {
    const unsigned elem1Weight = 1 * (column+1);
    const unsigned elem2Weight = 2 * (column+1);
    const unsigned elem3Weight = 3 * (column+1);
    std::vector<unsigned> expectedValues;
    if      (GetParam() == 1) { expectedValues = {0}; }
    else if (GetParam() == 2) { expectedValues = {elem1Weight+elem2Weight, elem3Weight}; }
    else if (GetParam() == 3) { expectedValues = {0,                       elem3Weight, elem1Weight+elem2Weight}; }
    else if (GetParam() == 4) { expectedValues = {elem1Weight,             elem2Weight, elem3Weight, 0}; }

    test_diag_multi_values<stk::balance::TotalElementWeightDiagnostic, unsigned>(column, expectedValues);
  }
}


TEST_F(TestDiagnosticsComputation, NodeInterfaceSize_Balance_HexMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_hex_in_square());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0.0/18.0}; }
  else if (get_parallel_size() == 2) { expectedValues = {6.0/12.0, 6.0/12.0}; }
  else if (get_parallel_size() == 3) { expectedValues = {6.0/12.0,  6.0/8.0, 6.0/8.0}; }
  else if (get_parallel_size() == 4) { expectedValues = { 6.0/8.0,  6.0/8.0, 6.0/8.0, 6.0/8.0}; }

  test_diag_values<stk::balance::RelativeNodeInterfaceSizeDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, NodeInterfaceSize_Balance_HexMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_hex_in_square());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0.0/18.0}; }
  else if (get_parallel_size() == 2) { expectedValues = {6.0/12.0, 6.0/12.0}; }
  else if (get_parallel_size() == 3) { expectedValues = { 6.0/8.0, 6.0/8.0, 6.0/12.0}; }
  else if (get_parallel_size() == 4) { expectedValues = { 6.0/8.0, 6.0/8.0,  6.0/8.0, 6.0/8.0}; }

  test_diag_values<stk::balance::RelativeNodeInterfaceSizeDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, NodeInterfaceSize_Balance_HexPyramidTetMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_hex_pyramid_tet());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  stk::io::write_mesh("nodeInterfaceSize_balance_hexPyramidTet_geometric.g", get_bulk());

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0.0/12.0}; }
  else if (get_parallel_size() == 2) { expectedValues = { 4.0/10.0, 4.0/6.0}; }
  else if (get_parallel_size() == 3) { expectedValues = { 4.0/8.0,  6.0/6.0, 4.0/5.0}; }
  else if (get_parallel_size() == 4) { expectedValues = { 4.0/8.0,  6.0/6.0, 4.0/5.0, 4.0/4.0}; }

  test_diag_values<stk::balance::RelativeNodeInterfaceSizeDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, NodeInterfaceSize_Balance_HexPyramidTetMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_hex_pyramid_tet());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  stk::io::write_mesh("nodeInterfaceSize_balance_hexPyramidTet_graph.g", get_bulk());

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0.0/12.0}; }
  else if (get_parallel_size() == 2) { expectedValues = {4.0/10.0, 4.0/6.0}; }
  else if (get_parallel_size() == 3) { expectedValues = {4.0/8.0,  4.0/4.0, 6.0/8.0}; }
  else if (get_parallel_size() == 4) { expectedValues = {4.0/8.0,  6.0/6.0, 4.0/4.0, 4.0/5.0}; }

  test_diag_values<stk::balance::RelativeNodeInterfaceSizeDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, NodeInterfaceSize_Balance_ShellMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_shells_in_square());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0.0/9.0}; }
  else if (get_parallel_size() == 2) { expectedValues = {3.0/6.0, 3.0/6.0}; }
  else if (get_parallel_size() == 3) { expectedValues = {3.0/6.0, 3.0/4.0, 3.0/4.0}; }
  else if (get_parallel_size() == 4) { expectedValues = {3.0/4.0, 3.0/4.0, 3.0/4.0, 3.0/4.0}; }

  test_diag_values<stk::balance::RelativeNodeInterfaceSizeDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, NodeInterfaceSize_Balance_ShellMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_shells_in_square());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0.0/9.0}; }
  else if (get_parallel_size() == 2) { expectedValues = {3.0/6.0, 3.0/6.0}; }
  else if (get_parallel_size() == 3) { expectedValues = {3.0/4.0, 3.0/4.0, 3.0/6.0}; }
  else if (get_parallel_size() == 4) { expectedValues = {3.0/4.0, 3.0/4.0, 3.0/4.0, 3.0/4.0}; }

  test_diag_values<stk::balance::RelativeNodeInterfaceSizeDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, NodeInterfaceSize_Balance_BeamMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_beams());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0.0/5.0}; }
  else if (get_parallel_size() == 2) { expectedValues = {1.0/3.0, 1.0/3.0}; }
  else if (get_parallel_size() == 3) { expectedValues = {1.0/3.0, 1.0/2.0, 1.0/2.0}; }
  else if (get_parallel_size() == 4) { expectedValues = {1.0/2.0, 2.0/2.0, 1.0/2.0, 1.0/2.0}; }

  test_diag_values<stk::balance::RelativeNodeInterfaceSizeDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, NodeInterfaceSize_Balance_BeamMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_beams());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {0.0/5.0}; }
  else if (get_parallel_size() == 2) { expectedValues = {1.0/3.0, 1.0/3.0}; }
  else if (get_parallel_size() == 3) { expectedValues = {1.0/2.0, 2.0/2.0, 1.0/3.0}; }
  else if (get_parallel_size() == 4) { expectedValues = {1.0/2.0, 2.0/2.0, 1.0/2.0, 1.0/2.0}; }

  test_diag_values<stk::balance::RelativeNodeInterfaceSizeDiagnostic, double>(expectedValues);
}


TEST_F(TestDiagnosticsComputation, ConnectivityWeight_Balance_HexMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_hex_in_square());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  const double cornerNode = 8.0;
  const double edgeNode = 12.0;
  const double centerNode = 18.0;
  const double elemWeight = (2*cornerNode + 4*edgeNode + 2*centerNode)/8;
  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {4*elemWeight}; }
  else if (get_parallel_size() == 2) { expectedValues = {2*elemWeight, 2*elemWeight}; }
  else if (get_parallel_size() == 3) { expectedValues = {2*elemWeight,   elemWeight, elemWeight}; }
  else if (get_parallel_size() == 4) { expectedValues = {  elemWeight,   elemWeight, elemWeight, elemWeight}; }

  test_diag_values<stk::balance::ConnectivityWeightDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, ConnectivityWeight_Balance_HexMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_hex_in_square());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  const double cornerNode = 8.0;
  const double edgeNode = 12.0;
  const double centerNode = 18.0;
  const double elemWeight = (2*cornerNode + 4*edgeNode + 2*centerNode)/8;
  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {4*elemWeight}; }
  else if (get_parallel_size() == 2) { expectedValues = {2*elemWeight, 2*elemWeight}; }
  else if (get_parallel_size() == 3) { expectedValues = {  elemWeight,   elemWeight, 2*elemWeight}; }
  else if (get_parallel_size() == 4) { expectedValues = {  elemWeight,   elemWeight,   elemWeight, elemWeight}; }

  test_diag_values<stk::balance::ConnectivityWeightDiagnostic, double>(expectedValues);
}

std::tuple<double, double, double, double, double, double>
get_hex_pyramid_tet_element_connectivity_weights()
{
  const double node1Weight =  8.0;
  const double node2Weight =  8.0;
  const double node3Weight =  8.0;
  const double node4Weight =  8.0;
  const double node5Weight = 10.0;
  const double node6Weight =  9.0;
  const double node7Weight = 10.0;
  const double node8Weight = 12.0;
  const double node9Weight =  5.0;
  const double node10Weight = 8.0;
  const double node11Weight = 5.0;
  const double node12Weight = 5.0;
  const double hexElemsPerNode = 1;
  const double pyrElemsPerNode = 6.0/2.0;
  const double tetElemsPerNode = 6;
  const double elem1Weight = (node1Weight + node2Weight + node3Weight + node4Weight +
                              node5Weight + node6Weight + node7Weight + node8Weight)/8/hexElemsPerNode;
  const double elem2Weight = (node5Weight + node6Weight + node7Weight + node8Weight + node10Weight)/5/pyrElemsPerNode;
  const double elem3Weight = (node5Weight + node9Weight + node8Weight + node10Weight)/4/tetElemsPerNode;
  const double elem4Weight = (node8Weight + node9Weight + node12Weight + node10Weight)/4/tetElemsPerNode;
  const double elem5Weight = (node8Weight + node12Weight + node10Weight + node11Weight)/4/tetElemsPerNode;
  const double elem6Weight = (node7Weight + node8Weight + node10Weight + node11Weight)/4/tetElemsPerNode;

  return std::make_tuple(elem1Weight, elem2Weight, elem3Weight, elem4Weight, elem5Weight, elem6Weight);
}

TEST_F(TestDiagnosticsComputation, ConnectivityWeight_Balance_HexPyramidTetMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_hex_pyramid_tet());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  double e1wt, e2wt, e3wt, e4wt, e5wt, e6wt;
  std::tie(e1wt, e2wt, e3wt, e4wt, e5wt, e6wt) = get_hex_pyramid_tet_element_connectivity_weights();

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {e1wt+e2wt+e3wt+e4wt+e5wt+e6wt}; }
  else if (get_parallel_size() == 2) { expectedValues = {e1wt+e2wt+e3wt, e4wt+e5wt+e6wt}; }
  else if (get_parallel_size() == 3) { expectedValues = {e1wt, e2wt+e3wt+e6wt, e4wt+e5wt}; }
  else if (get_parallel_size() == 4) { expectedValues = {e1wt, e2wt+e3wt, e4wt+e5wt, e6wt}; }

  test_diag_values<stk::balance::ConnectivityWeightDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, ConnectivityWeight_Balance_HexPyramidTetMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_hex_pyramid_tet());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  double e1wt, e2wt, e3wt, e4wt, e5wt, e6wt;
  std::tie(e1wt, e2wt, e3wt, e4wt, e5wt, e6wt) = get_hex_pyramid_tet_element_connectivity_weights();

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {e1wt+e2wt+e3wt+e4wt+e5wt+e6wt}; }
  else if (get_parallel_size() == 2) { expectedValues = {e1wt+e2wt+e3wt, e4wt+e5wt+e6wt}; }
  else if (get_parallel_size() == 3) { expectedValues = {e1wt, e6wt, e2wt+e3wt+e4wt+e5wt}; }
  else if (get_parallel_size() == 4) { expectedValues = {e1wt, e2wt+e3wt, e6wt, e4wt+e5wt}; }

  test_diag_values<stk::balance::ConnectivityWeightDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, ConnectivityWeight_Balance_ShellMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_shells_in_square());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  const double cornerNode = 4.0;
  const double edgeNode = 6.0;
  const double centerNode = 9.0;
  const double quadShellElemsPerNode = 1.0;
  const double elemWeight = (cornerNode + 2*edgeNode + centerNode)/4/quadShellElemsPerNode;
  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {4*elemWeight}; }
  else if (get_parallel_size() == 2) { expectedValues = {2*elemWeight, 2*elemWeight}; }
  else if (get_parallel_size() == 3) { expectedValues = {2*elemWeight,   elemWeight, elemWeight}; }
  else if (get_parallel_size() == 4) { expectedValues = {  elemWeight,   elemWeight, elemWeight, elemWeight}; }

  test_diag_values<stk::balance::ConnectivityWeightDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, ConnectivityWeight_Balance_ShellMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_shells_in_square());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  const double cornerNode = 4.0;
  const double edgeNode = 6.0;
  const double centerNode = 9.0;
  const double quadShellElemsPerNode = 1.0;
  const double elemWeight = (cornerNode + 2*edgeNode + centerNode)/4/quadShellElemsPerNode;
  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {4*elemWeight}; }
  else if (get_parallel_size() == 2) { expectedValues = {2*elemWeight, 2*elemWeight}; }
  else if (get_parallel_size() == 3) { expectedValues = {  elemWeight,   elemWeight, 2*elemWeight}; }
  else if (get_parallel_size() == 4) { expectedValues = {  elemWeight,   elemWeight,   elemWeight, elemWeight}; }

  test_diag_values<stk::balance::ConnectivityWeightDiagnostic, double>(expectedValues);
}

std::tuple<double, double, double, double>
get_beam_element_connectivity_weights()
{
  const double node1Weight = 2.0;
  const double node2Weight = 3.0;
  const double node3Weight = 4.0;
  const double node4Weight = 2.0;
  const double node5Weight = 2.0;
  const double beamElemsPerNode = 1.0;
  const double elem1Weight = (node1Weight + node2Weight)/2/beamElemsPerNode;
  const double elem2Weight = (node2Weight + node3Weight)/2/beamElemsPerNode;
  const double elem3Weight = (node3Weight + node4Weight)/2/beamElemsPerNode;
  const double elem4Weight = (node3Weight + node5Weight)/2/beamElemsPerNode;

  return std::make_tuple(elem1Weight, elem2Weight, elem3Weight, elem4Weight);
}

TEST_F(TestDiagnosticsComputation, ConnectivityWeight_Balance_BeamMesh_GeometricPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_beams());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setDecompMethod("rcb");
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  double e1wt, e2wt, e3wt, e4wt;
  std::tie(e1wt, e2wt, e3wt, e4wt) = get_beam_element_connectivity_weights();

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {e1wt+e2wt+e3wt+e4wt}; }
  else if (get_parallel_size() == 2) { expectedValues = {e1wt+e2wt, e3wt+e4wt}; }
  else if (get_parallel_size() == 3) { expectedValues = {e1wt+e2wt, e3wt, e4wt}; }
  else if (get_parallel_size() == 4) { expectedValues = {e1wt, e2wt, e3wt, e4wt}; }

  test_diag_values<stk::balance::ConnectivityWeightDiagnostic, double>(expectedValues);
}

TEST_F(TestDiagnosticsComputation, ConnectivityWeight_Balance_BeamMesh_GraphPartitioner)
{
  if (get_parallel_size() > 4) return;

  build_mesh(mesh_desc_four_beams());

  stk::balance::GraphCreationSettings balanceSettings;
  balanceSettings.setShouldPrintDiagnostics(true);

  stk::balance::balanceStkMesh(balanceSettings, get_bulk());

  double e1wt, e2wt, e3wt, e4wt;
  std::tie(e1wt, e2wt, e3wt, e4wt) = get_beam_element_connectivity_weights();

  std::vector<double> expectedValues;
  if      (get_parallel_size() == 1) { expectedValues = {e1wt+e2wt+e3wt+e4wt}; }
  else if (get_parallel_size() == 2) { expectedValues = {e1wt+e2wt, e3wt+e4wt}; }
  else if (get_parallel_size() == 3) { expectedValues = {e1wt, e2wt, e3wt+e4wt}; }
  else if (get_parallel_size() == 4) { expectedValues = {e1wt, e2wt, e4wt, e3wt}; }

  test_diag_values<stk::balance::ConnectivityWeightDiagnostic, double>(expectedValues);
}

INSTANTIATE_TEST_SUITE_P(TestDiagnosticsComputation, RebalanceNumOutputProcs, testing::Values(1, 2, 3, 4));
