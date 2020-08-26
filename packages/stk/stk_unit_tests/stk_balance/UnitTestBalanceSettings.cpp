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
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_balance/balanceUtils.hpp"

class BalanceSettingsTester : public stk::unit_test_util::MeshFixture
{
public:
  stk::balance::BalanceSettings& get_balance_settings()
  {
    return m_settings;
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

  stk::balance::UserSpecifiedVertexWeightsSetting userSpecifiedVertexWeightsSetting;
  EXPECT_FALSE(userSpecifiedVertexWeightsSetting.includeSearchResultsInGraph());

  stk::balance::GraphCreationSettingsForZoltan2 graphCreationSettingsForZoltan2;
  EXPECT_TRUE(graphCreationSettingsWithCustomTolerances.includeSearchResultsInGraph());

  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  stk::balance::DoubleFieldType& field = meta.declare_field<stk::balance::DoubleFieldType>(stk::topology::NODE_RANK, "dummy", 1);
  stk::balance::FieldVertexWeightSettings fieldVertexWeightSettings(bulk, field);
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

  stk::balance::UserSpecifiedVertexWeightsSetting userSpecifiedVertexWeightsSetting;
  userSpecifiedVertexWeightsSetting.setIncludeSearchResultsInGraph(true);
  EXPECT_TRUE(userSpecifiedVertexWeightsSetting.includeSearchResultsInGraph());

  stk::balance::GraphCreationSettingsForZoltan2 graphCreationSettingsForZoltan2;
  graphCreationSettingsForZoltan2.setIncludeSearchResultsInGraph(false);
  EXPECT_FALSE(graphCreationSettingsWithCustomTolerances.includeSearchResultsInGraph());

  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  stk::balance::DoubleFieldType& field = meta.declare_field<stk::balance::DoubleFieldType>(stk::topology::NODE_RANK, "dummy", 1);
  stk::balance::FieldVertexWeightSettings fieldVertexWeightSettings(bulk, field);
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

