#include <iosfwd>
#include <fstream>
#include "gtest/gtest.h"

#include "stk_balance/balanceUtils.hpp"
#include "stk_balance/setup/Parser.hpp"

#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_mesh/base/GetEntities.hpp"

class BalanceCommandLine : public stk::unit_test_util::MeshFixture
{
protected:
  BalanceCommandLine()
    : infile("mesh.g"),
      m_execName("stk_balance"),
      m_outputDir(""),
      m_parser(get_comm()),
      m_fullInfileOutputDirOptions(false),
      m_infileOptionName("--infile"),
      m_outputDirOptionName("--output-directory")
  {
    setup_mesh("generated:1x1x4|sideset:x", stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), m_faces);
    ThrowRequire(m_faces.size() > 0);
  }

  void use_full_infile_output_dir_options()
  {
    m_fullInfileOutputDirOptions = true;
  }

  void set_output_directory(const std::string& dir)
  {
    m_outputDir = dir;
  }

  void turn_on_contact_search_before_parsing()
  {
    m_settings.setIncludeSearchResultsInGraph(true);
  }

  void turn_off_contact_search_before_parsing()
  {
    m_settings.setIncludeSearchResultsInGraph(false);
  }

  const stk::balance::BalanceSettings& get_stk_balance_settings(const std::vector<std::string>& options)
  {
    std::vector<const char*> args = assemble_args(options);
    m_parser.parse_command_line_options(args.size(), args.data(), m_settings);
    return m_settings;
  }

  std::vector<const char*> assemble_args(const std::vector<std::string>& options) const
  {
    std::vector<const char*> args = {m_execName.c_str()};

    if (m_fullInfileOutputDirOptions) args.push_back(m_infileOptionName.c_str());
    args.push_back(infile.c_str());

    if (!m_outputDir.empty()) {
      if (m_fullInfileOutputDirOptions) args.push_back(m_outputDirOptionName.c_str());
      args.push_back(m_outputDir.c_str());
    }

    for (const std::string& option : options) {
      args.push_back(option.c_str());
    }

    return args;
  }

  double get_absolute_tolerance_for_unit_mesh(const stk::balance::BalanceSettings& balanceSettings)
  {
    return balanceSettings.getToleranceForFaceSearch(get_bulk(), *get_meta().coordinate_field(),
        get_bulk().begin_nodes(m_faces[0]), get_bulk().num_nodes(m_faces[0]));
  }

  void check_absolute_tolerance_for_face_search(const stk::balance::BalanceSettings& balanceSettings, double tolerance)
  {
    EXPECT_TRUE(balanceSettings.isConstantFaceSearchTolerance());
    EXPECT_DOUBLE_EQ(get_absolute_tolerance_for_unit_mesh(balanceSettings), tolerance);
  }

  void check_relative_tolerance_for_face_search(const stk::balance::BalanceSettings& balanceSettings, double tolerance)
  {
    EXPECT_FALSE(balanceSettings.isConstantFaceSearchTolerance());
    EXPECT_DOUBLE_EQ(get_absolute_tolerance_for_unit_mesh(balanceSettings), tolerance);
  }

  std::string infile;

private:
  std::string m_execName;
  std::string m_outputDir;
  stk::balance::StkBalanceSettings m_settings;
  stk::balance::Parser m_parser;
  stk::mesh::EntityVector m_faces;

  bool m_fullInfileOutputDirOptions;
  std::string m_infileOptionName;
  std::string m_outputDirOptionName;
};

TEST_F(BalanceCommandLine, createBalanceSettings_default)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), "./"+infile);

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_outputDirectory)
{
  const std::string outputDir = "output";
  set_output_directory(outputDir);
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), outputDir+"/"+infile);

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_outputDirectory_fullOptions)
{
  const std::string outputDir = "output";
  set_output_directory(outputDir);
  use_full_infile_output_dir_options();
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), outputDir+"/"+infile);

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  check_relative_tolerance_for_face_search(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_bothDefaults)
{
  EXPECT_THROW(get_stk_balance_settings({"--sd", "--sm"}), std::logic_error);
}


TEST_F(BalanceCommandLine, createBalanceSettings_defaultAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--face-search-abs-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_defaultRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--face-search-rel-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_relative_tolerance_for_face_search(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, createBalanceSettings_bothTolerances)
{
  EXPECT_THROW(get_stk_balance_settings({"--face-search-abs-tol", "--face-search-rel-tol"}), std::logic_error);
}


TEST_F(BalanceCommandLine, createBalanceSettings_faceSearchAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--face-search-abs-tol=0.001"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_faceSearchRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--face-search-rel-tol=0.123"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_relative_tolerance_for_face_search(balanceSettings, 0.123);
}


TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_defaultAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--face-search-abs-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_defaultRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--face-search-rel-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  check_relative_tolerance_for_face_search(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_faceSearchAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--face-search-abs-tol=0.005"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.005);
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_faceSearchRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--face-search-rel-tol=0.123"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  check_relative_tolerance_for_face_search(balanceSettings, 0.123);
}


TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_defaultAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--face-search-abs-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_defaultRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--face-search-rel-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_relative_tolerance_for_face_search(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_faceSearchAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--face-search-abs-tol=0.0005"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0005);
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_faceSearchRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--face-search-rel-tol=0.123"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_relative_tolerance_for_face_search(balanceSettings, 0.123);
}


TEST_F(BalanceCommandLine, contactSearch_respectsInitialDefault_enabled)
{
  turn_on_contact_search_before_parsing();
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
}

TEST_F(BalanceCommandLine, contactSearch_respectsInitialDefault_disabled)
{
  turn_off_contact_search_before_parsing();
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  EXPECT_FALSE(balanceSettings.includeSearchResultsInGraph());
}

TEST_F(BalanceCommandLine, contactSearch_badValue)
{
  EXPECT_THROW(get_stk_balance_settings({"--contact-search=junk"}), std::logic_error);
}

TEST_F(BalanceCommandLine, disableSearch_default_caseInsensitive)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--contact-search=Off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_FALSE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, disableSearch_smDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--contact-search=off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_FALSE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  check_relative_tolerance_for_face_search(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, disableSearch_sdDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--contact-search=off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_FALSE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}


TEST_F(BalanceCommandLine, enableSearch_default_caseInsensitive)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--contact-search=ON"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, enableSearch_whenDisabledByDefault)
{
  turn_off_contact_search_before_parsing();
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--contact-search=ON"});

  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
}

TEST_F(BalanceCommandLine, enableSearch_smDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--contact-search=on"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  check_relative_tolerance_for_face_search(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, enableSearch_sdDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--contact-search=on"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}


TEST_F(BalanceCommandLine, decompMethodRcb)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--decomp-method=rcb"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "rcb");
}

TEST_F(BalanceCommandLine, decompMethodRib)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--decomp-method=rib"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "rib");
}

TEST_F(BalanceCommandLine, decompMethodMultijagged)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--decomp-method=multijagged"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "multijagged");
}

TEST_F(BalanceCommandLine, decompMethodParmetis)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--decomp-method=parmetis"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  check_absolute_tolerance_for_face_search(balanceSettings, 0.0001);
}
