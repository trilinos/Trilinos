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

#include "gtest/gtest.h"
#include "stk_balance/balanceUtils.hpp"
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace {

using stk::unit_test_util::build_mesh;

class BalanceSettingsTester : public ::testing::Test
{
public:
  stk::balance::BalanceSettings& get_balance_settings()
  {
    return m_settings;
  }

  int numRequiredCommonNodes(const std::vector<stk::topology::topology_t> & topologyListA,
                             const std::vector<stk::topology::topology_t> & topologyListB)
  {
    const int numNodes = m_settings.getNumNodesRequiredForConnection(topologyListA.front(), topologyListB.front());

    for (auto topologyA : topologyListA) {
      for (auto topologyB : topologyListB) {
        int numNodes1 = m_settings.getNumNodesRequiredForConnection(topologyA, topologyB);
        int numNodes2 = m_settings.getNumNodesRequiredForConnection(topologyB, topologyA);
        EXPECT_EQ(numNodes1, numNodes) << "number of nodes for " << topologyA << " and " << topologyB
                                       << " do not match expected value of " << numNodes << " for group";
        EXPECT_EQ(numNodes2, numNodes) << "number of nodes for " << topologyB << " and " << topologyA
                                       << " do not match expected value of " << numNodes << " for group";
      }
    }
    return numNodes;
  }

private:
  stk::balance::StkBalanceSettings m_settings;
};

TEST_F(BalanceSettingsTester, defaultInitialDecompMethod)
{
  stk::balance::BalanceSettings& settings = get_balance_settings();
  EXPECT_EQ(settings.getInitialDecompMethod(), "RIB");
}

TEST_F(BalanceSettingsTester, inputFilename)
{
  stk::balance::BalanceSettings& settings = get_balance_settings();
  settings.set_input_filename("input.g");
  EXPECT_EQ(settings.get_input_filename(), "input.g");
  settings.set_input_filename("cube.g");
  EXPECT_EQ(settings.get_input_filename(), "cube.g");
}

TEST_F(BalanceSettingsTester, outputFilename)
{
  stk::balance::BalanceSettings& settings = get_balance_settings();
  settings.set_output_filename("output.g");
  EXPECT_EQ(settings.get_output_filename(), "output.g");
  settings.set_output_filename("output/cube.g");
  EXPECT_EQ(settings.get_output_filename(), "output/cube.g");
}

TEST_F(BalanceSettingsTester, getGraphVertexWeight_supports_PYRAMID_5)
{
  stk::balance::BalanceSettings& settings = get_balance_settings();
  EXPECT_EQ(1, settings.getGraphVertexWeight(stk::topology::PYRAMID_5));
}

TEST(BalanceSettings, defaultContactSearchStatus)
{
  stk::balance::BalanceSettings balanceSettings;
  EXPECT_FALSE(balanceSettings.includeSearchResultsInGraph());

  stk::balance::StkBalanceSettings stkBalanceSettings;
  EXPECT_TRUE(stkBalanceSettings.includeSearchResultsInGraph());

  stk::balance::BasicGeometricSettings basicGeometricSettings;
  EXPECT_FALSE(basicGeometricSettings.includeSearchResultsInGraph());

  stk::balance::GraphCreationSettings graphCreationSettings;
  EXPECT_TRUE(graphCreationSettings.includeSearchResultsInGraph());

  stk::balance::GraphCreationSettingsWithCustomTolerances  graphCreationSettingsWithCustomTolerances;
  EXPECT_TRUE(graphCreationSettingsWithCustomTolerances.includeSearchResultsInGraph());

  stk::balance::BasicZoltan2Settings basicZoltan2Settings;
  EXPECT_FALSE(basicZoltan2Settings.includeSearchResultsInGraph());

  EXPECT_TRUE(graphCreationSettingsWithCustomTolerances.includeSearchResultsInGraph());

  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, MPI_COMM_WORLD);
  stk::balance::DoubleFieldType& field = bulk->mesh_meta_data().declare_field<double>(stk::topology::NODE_RANK, "dummy", 1);
  stk::balance::FieldVertexWeightSettings fieldVertexWeightSettings(*bulk, field);
  EXPECT_FALSE(fieldVertexWeightSettings.includeSearchResultsInGraph());

  stk::balance::MultipleCriteriaSettings multipleCriteriaSettings({&field});
  EXPECT_FALSE(multipleCriteriaSettings.includeSearchResultsInGraph());

  stk::balance::BasicColoringSettings basicColoringSettings;
  EXPECT_FALSE(basicColoringSettings.includeSearchResultsInGraph());

  stk::balance::BasicColoringByTopologySettings basicColoringByTopologySettings;
  EXPECT_FALSE(basicColoringByTopologySettings.includeSearchResultsInGraph());
}

TEST(BalanceSettings, toggleContactSearchStatus)
{
  stk::balance::BalanceSettings balanceSettings;
  balanceSettings.setIncludeSearchResultsInGraph(true);
  EXPECT_FALSE(balanceSettings.includeSearchResultsInGraph());

  stk::balance::StkBalanceSettings stkBalanceSettings;
  stkBalanceSettings.setIncludeSearchResultsInGraph(false);
  EXPECT_FALSE(stkBalanceSettings.includeSearchResultsInGraph());

  stk::balance::BasicGeometricSettings basicGeometricSettings;
  basicGeometricSettings.setIncludeSearchResultsInGraph(true);
  EXPECT_FALSE(basicGeometricSettings.includeSearchResultsInGraph());

  stk::balance::GraphCreationSettings graphCreationSettings;
  graphCreationSettings.setIncludeSearchResultsInGraph(false);
  EXPECT_FALSE(graphCreationSettings.includeSearchResultsInGraph());

  stk::balance::GraphCreationSettingsWithCustomTolerances  graphCreationSettingsWithCustomTolerances;
  graphCreationSettingsWithCustomTolerances.setIncludeSearchResultsInGraph(false);
  EXPECT_FALSE(graphCreationSettingsWithCustomTolerances.includeSearchResultsInGraph());

  stk::balance::BasicZoltan2Settings basicZoltan2Settings;
  basicZoltan2Settings.setIncludeSearchResultsInGraph(true);
  EXPECT_TRUE(basicZoltan2Settings.includeSearchResultsInGraph());

  stk::balance::GraphCreationSettingsForZoltan2 graphCreationSettingsForZoltan2;
  graphCreationSettingsForZoltan2.setIncludeSearchResultsInGraph(false);
  EXPECT_FALSE(graphCreationSettingsWithCustomTolerances.includeSearchResultsInGraph());

  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, MPI_COMM_WORLD);
  stk::balance::DoubleFieldType& field = bulk->mesh_meta_data().declare_field<double>(stk::topology::NODE_RANK, "dummy", 1);
  stk::balance::FieldVertexWeightSettings fieldVertexWeightSettings(*bulk, field);
  fieldVertexWeightSettings.setIncludeSearchResultsInGraph(true);
  EXPECT_TRUE(fieldVertexWeightSettings.includeSearchResultsInGraph());

  stk::balance::MultipleCriteriaSettings multipleCriteriaSettings({&field});
  multipleCriteriaSettings.setIncludeSearchResultsInGraph(true);
  EXPECT_TRUE(multipleCriteriaSettings.includeSearchResultsInGraph());

  stk::balance::BasicColoringSettings basicColoringSettings;
  basicColoringSettings.setIncludeSearchResultsInGraph(true);
  EXPECT_FALSE(basicColoringSettings.includeSearchResultsInGraph());

  stk::balance::BasicColoringByTopologySettings basicColoringByTopologySettings;
  basicColoringByTopologySettings.setIncludeSearchResultsInGraph(true);
  EXPECT_FALSE(basicColoringByTopologySettings.includeSearchResultsInGraph());
}

std::vector<stk::topology::topology_t> get_0dim_topologies() {
  return {stk::topology::PARTICLE};
}

std::vector<stk::topology::topology_t> get_1dim_topologies() {
  return {stk::topology::LINE_2,
          stk::topology::LINE_2_1D,
          stk::topology::LINE_3_1D,
          stk::topology::BEAM_2,
          stk::topology::BEAM_3,
          stk::topology::SHELL_LINE_2,
          stk::topology::SHELL_LINE_3,
          stk::topology::SPRING_2,
          stk::topology::SPRING_3};
}

std::vector<stk::topology::topology_t> get_2dim_topologies() {
  return {stk::topology::TRI_3_2D,
          stk::topology::TRI_4_2D,
          stk::topology::QUAD_4_2D,
          stk::topology::SHELL_TRI_3,
          stk::topology::SHELL_TRI_4,
          stk::topology::SHELL_QUAD_4};
}

std::vector<stk::topology::topology_t> get_3dim_topologies() {
  return {stk::topology::TET_4,
          stk::topology::PYRAMID_5,
          stk::topology::WEDGE_6,
          stk::topology::HEX_8};
}

std::vector<stk::topology::topology_t> get_2dim_2ndOrder_topologies() {
  return {stk::topology::TRI_6_2D,
          stk::topology::QUAD_8_2D,
          stk::topology::QUAD_9_2D,
          stk::topology::SHELL_TRI_6,
          stk::topology::SHELL_QUAD_8,
          stk::topology::SHELL_QUAD_9};
}

std::vector<stk::topology::topology_t> get_3dim_2ndOrder_topologies() {
  return {stk::topology::TET_8,
          stk::topology::TET_10,
          stk::topology::TET_11,
          stk::topology::PYRAMID_13,
          stk::topology::PYRAMID_14,
          stk::topology::WEDGE_12,
          stk::topology::WEDGE_15,
          stk::topology::WEDGE_18,
          stk::topology::HEX_20,
          stk::topology::HEX_27};
}

std::vector<stk::topology::topology_t> get_super_element_topologies() {
  return {stk::create_superelement_topology(1),
          stk::create_superelement_topology(10),
          stk::create_superelement_topology(100)};
}

constexpr int NoConnection = 1000;

TEST_F(BalanceSettingsTester, numNodesRequiredForConnection)
{
  EXPECT_EQ(numRequiredCommonNodes(get_0dim_topologies(), get_0dim_topologies()),          1);
  EXPECT_EQ(numRequiredCommonNodes(get_0dim_topologies(), get_1dim_topologies()),          1);
  EXPECT_EQ(numRequiredCommonNodes(get_0dim_topologies(), get_2dim_topologies()),          1);
  EXPECT_EQ(numRequiredCommonNodes(get_0dim_topologies(), get_3dim_topologies()),          1);
  EXPECT_EQ(numRequiredCommonNodes(get_0dim_topologies(), get_2dim_2ndOrder_topologies()), 1);
  EXPECT_EQ(numRequiredCommonNodes(get_0dim_topologies(), get_3dim_2ndOrder_topologies()), 1);
  EXPECT_EQ(numRequiredCommonNodes(get_0dim_topologies(), get_super_element_topologies()), NoConnection);

  EXPECT_EQ(numRequiredCommonNodes(get_1dim_topologies(), get_1dim_topologies()),          1);
  EXPECT_EQ(numRequiredCommonNodes(get_1dim_topologies(), get_2dim_topologies()),          1);
  EXPECT_EQ(numRequiredCommonNodes(get_1dim_topologies(), get_3dim_topologies()),          1);
  EXPECT_EQ(numRequiredCommonNodes(get_1dim_topologies(), get_2dim_2ndOrder_topologies()), 1);
  EXPECT_EQ(numRequiredCommonNodes(get_1dim_topologies(), get_3dim_2ndOrder_topologies()), 1);
  EXPECT_EQ(numRequiredCommonNodes(get_1dim_topologies(), get_super_element_topologies()), NoConnection);

  EXPECT_EQ(numRequiredCommonNodes(get_2dim_topologies(), get_2dim_topologies()),          2);
  EXPECT_EQ(numRequiredCommonNodes(get_2dim_topologies(), get_3dim_topologies()),          2);
  EXPECT_EQ(numRequiredCommonNodes(get_2dim_topologies(), get_2dim_2ndOrder_topologies()), 2);
  EXPECT_EQ(numRequiredCommonNodes(get_2dim_topologies(), get_3dim_2ndOrder_topologies()), 2);
  EXPECT_EQ(numRequiredCommonNodes(get_2dim_topologies(), get_super_element_topologies()), NoConnection);

  EXPECT_EQ(numRequiredCommonNodes(get_3dim_topologies(), get_3dim_topologies()),          3);
  EXPECT_EQ(numRequiredCommonNodes(get_3dim_topologies(), get_2dim_2ndOrder_topologies()), 3);
  EXPECT_EQ(numRequiredCommonNodes(get_3dim_topologies(), get_3dim_2ndOrder_topologies()), 3);
  EXPECT_EQ(numRequiredCommonNodes(get_3dim_topologies(), get_super_element_topologies()), NoConnection);

  EXPECT_EQ(numRequiredCommonNodes(get_2dim_2ndOrder_topologies(), get_2dim_2ndOrder_topologies()), 3);
  EXPECT_EQ(numRequiredCommonNodes(get_2dim_2ndOrder_topologies(), get_3dim_2ndOrder_topologies()), 3);
  EXPECT_EQ(numRequiredCommonNodes(get_2dim_2ndOrder_topologies(), get_super_element_topologies()), NoConnection);

  EXPECT_EQ(numRequiredCommonNodes(get_3dim_2ndOrder_topologies(), get_3dim_2ndOrder_topologies()), 4);
  EXPECT_EQ(numRequiredCommonNodes(get_3dim_2ndOrder_topologies(), get_super_element_topologies()), NoConnection);

  EXPECT_EQ(numRequiredCommonNodes(get_super_element_topologies(), get_super_element_topologies()), NoConnection);
}

#define PRINT_AND_RETHROW(CODE_BLOCK, EXCEPT_STREAM)                           \
  try {                                                                        \
    do {                                                                       \
      CODE_BLOCK;                                                              \
    }                                                                          \
    while(0);                                                                  \
  }                                                                            \
  catch(const std::exception& ex) {                                            \
    EXCEPT_STREAM << "std::exception thrown: " << ex.what() << std::endl;      \
    throw;                                                                     \
  }                                                                            \
  catch(...) {                                                                 \
    EXCEPT_STREAM << "Unknown structure thrown" << std::endl;                  \
    throw;                                                                     \
  }

#define EXPECT_NO_THROW_PRINT(CODE_BLOCK) EXPECT_NO_THROW(PRINT_AND_RETHROW(CODE_BLOCK, std::cerr))


TEST(BalanceSettings, getGraphVertexWeight_allValidElementTopologies)
{
  stk::balance::GraphCreationSettings settings;
  for (stk::topology t = stk::topology::BEGIN_ELEMENT_RANK; t < stk::topology::END_ELEMENT_RANK; ++t) {
    EXPECT_NO_THROW_PRINT(settings.getGraphVertexWeight(t));
  }
}

TEST(BalanceSettings, getGraphEdgeWeight_allValidElementTopologies)
{
  stk::balance::GraphCreationSettings settings;
  for (stk::topology t = stk::topology::BEGIN_ELEMENT_RANK; t < stk::topology::END_ELEMENT_RANK; ++t) {
    EXPECT_NO_THROW_PRINT(settings.getGraphEdgeWeight(t, t));
  }
}

TEST(BalanceSettings, getNumNodesRequiredForConnection_allValidElementTopologies)
{
  stk::balance::GraphCreationSettings settings;
  for (stk::topology t = stk::topology::BEGIN_ELEMENT_RANK; t < stk::topology::END_ELEMENT_RANK; ++t) {
    EXPECT_NO_THROW_PRINT(settings.getNumNodesRequiredForConnection(t, t));
  }
}

}
