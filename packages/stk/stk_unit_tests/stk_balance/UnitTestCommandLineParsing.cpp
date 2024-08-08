#include <iosfwd>
#include <fstream>
#include "gtest/gtest.h"

#include "stk_balance/balanceUtils.hpp"
#include "stk_balance/setup/Parser.hpp"

#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_balance/setup/DefaultSettings.hpp"

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
    setup_mesh("generated:1x1x4|sideset:x", stk::mesh::BulkData::AUTO_AURA);
    stk::mesh::get_entities(get_bulk(), stk::topology::FACE_RANK, get_meta().universal_part(), m_faces);
    STK_ThrowRequire(m_faces.size() > 0);
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

  void check_vertex_weight_block_multiplier(const stk::balance::BalanceSettings& balanceSettings,
                                            const std::vector<std::pair<std::string, double>> & expectedMultipliers)
  {
    const std::map<std::string, double> & blockMultipliers = balanceSettings.getVertexWeightBlockMultipliers();
    EXPECT_EQ(blockMultipliers.size(), expectedMultipliers.size()) << "Unexpected number of block multipliers";
    for (const auto & expectedMultiplier : expectedMultipliers) {
      const auto foundMultiplier = blockMultipliers.find(expectedMultiplier.first);
      ASSERT_TRUE(foundMultiplier != blockMultipliers.end()) << "No multiplier found for block " << expectedMultiplier.first;
      EXPECT_DOUBLE_EQ(foundMultiplier->second, expectedMultiplier.second)
        << "Multiplier for block " << expectedMultiplier.first << " (" << foundMultiplier->second
        << ") doesn't match expected (" << expectedMultiplier.second << ")";
    }
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

using namespace stk::balance;

TEST_F(BalanceCommandLine, unrecognizedOption)
{
  EXPECT_ANY_THROW(get_stk_balance_settings({"--unrecognizedOption=3"}));
}

TEST_F(BalanceCommandLine, createBalanceSettings_default)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  const int finalNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD); 

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh.1_to_" + std::to_string(finalNumProcs) + ".log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_outputDirectory)
{
  const std::string outputDir = "output";
  set_output_directory(outputDir);
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  const int finalNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD); 

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), outputDir+"/"+infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh.1_to_" + std::to_string(finalNumProcs) + ".log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_outputDirectory_fullOptions)
{
  const std::string outputDir = "output";
  set_output_directory(outputDir);
  use_full_infile_output_dir_options();
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  const int finalNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD); 

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), outputDir+"/"+infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh.1_to_" + std::to_string(finalNumProcs) + ".log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_customLogfile)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--logfile=custom.log"});

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "custom.log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_shortCustomLogfile)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"-l", "custom.log"});

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "custom.log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_coutLogfile)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--logfile=cout"});

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "cout");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_printDiagnostics)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--print-diagnostics"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            true);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_shortPrintDiagnostics)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"-d"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            true);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_rebalanceTo)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--rebalance-to=16"});

  const int initialNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD);

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh." + std::to_string(initialNumProcs) + "_to_16.log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_useNestedDecomp)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--rebalance-to=12", "--use-nested-decomp"});

  const int initialNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD);

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh." + std::to_string(initialNumProcs) + "_to_12.log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             true);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaultsOverrideSpider)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--fix-spiders=on"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  true);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaultsOverrideMechanism)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--fix-mechanisms=off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaultsOverrideSpiders)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--fix-spiders=off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  false);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaultsOverrideMechanisms)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--fix-mechanisms=off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_bothDefaults)
{
  EXPECT_THROW(get_stk_balance_settings({"--sd", "--sm"}), std::logic_error);
}


TEST_F(BalanceCommandLine, createBalanceSettings_defaultAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--face-search-abs-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_absolute_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchAbsTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_defaultRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--face-search-rel-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_bothTolerances)
{
  EXPECT_THROW(get_stk_balance_settings({"--face-search-abs-tol", "--face-search-rel-tol"}), std::logic_error);
}


TEST_F(BalanceCommandLine, createBalanceSettings_faceSearchAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--face-search-abs-tol=0.001"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_absolute_tolerance_for_face_search(balanceSettings,                      0.001);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_faceSearchRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--face-search-rel-tol=0.123"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      0.123);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}


TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_defaultAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--face-search-abs-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_absolute_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchAbsTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_defaultRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--face-search-rel-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_faceSearchAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--face-search-abs-tol=0.005"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_absolute_tolerance_for_face_search(balanceSettings,                      0.005);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_faceSearchRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--face-search-rel-tol=0.123"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      0.123);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}


TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_defaultAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--face-search-abs-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_absolute_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchAbsTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_defaultRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--face-search-rel-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_faceSearchAbsoluteTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--face-search-abs-tol=0.0005"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_absolute_tolerance_for_face_search(balanceSettings,                      0.0005);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_faceSearchRelativeTolerance)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--face-search-rel-tol=0.123"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      0.123);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}


TEST_F(BalanceCommandLine, createBalanceSettings_contactSearchEdgeWeight)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--contact-search-edge-weight=20"});

  const int finalNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD);

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh.1_to_" + std::to_string(finalNumProcs) + ".log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                20);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_contactSearchEdgeWeight)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm",
                                                                                   "--contact-search-edge-weight=20"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                20);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_contactSearchEdgeWeight)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd",
                                                                                   "--contact-search-edge-weight=20"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                20);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}


TEST_F(BalanceCommandLine, createBalanceSettings_contactSearchVertexWeightMultiplier)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--contact-search-vertex-weight-mult=9"});

  const int finalNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD);

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh.1_to_" + std::to_string(finalNumProcs) + ".log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 9);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_contactSearchVertexWeightMultiplier)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm",
                                                                                   "--contact-search-vertex-weight-mult=9"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 9);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_contactSearchVertexWeightMultiplier)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd",
                                                                                   "--contact-search-vertex-weight-mult=9"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 9);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}


TEST_F(BalanceCommandLine, createBalanceSettings_edgeWeightMultiplier)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--edge-weight-mult=3"});

  const int finalNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD);

  EXPECT_EQ(balanceSettings.get_input_filename(), infile);
  EXPECT_EQ(balanceSettings.get_output_filename(), infile);
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh.1_to_" + std::to_string(finalNumProcs) + ".log");

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               3);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_edgeWeightMultiplier)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm",
                                                                                   "--edge-weight-mult=3"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               3);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_edgeWeightMultiplier)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd",
                                                                                   "--edge-weight-mult=3"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               3);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}


TEST_F(BalanceCommandLine, contactSearch_respectsInitialDefault_enabled)
{
  turn_on_contact_search_before_parsing();
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(), DefaultSettings::useContactSearch);
}

TEST_F(BalanceCommandLine, contactSearch_respectsInitialDefault_disabled)
{
  turn_off_contact_search_before_parsing();
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({});

  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(), false);
}

TEST_F(BalanceCommandLine, contactSearch_badValue)
{
  EXPECT_THROW(get_stk_balance_settings({"--contact-search=junk"}), std::logic_error);
}

TEST_F(BalanceCommandLine, disableSearch_default_caseInsensitive)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--contact-search=Off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       false);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, disableSearch_smDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--contact-search=off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       false);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, disableSearch_sdDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--contact-search=off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       false);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}


TEST_F(BalanceCommandLine, enableSearch_default_caseInsensitive)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--contact-search=ON"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, enableSearch_whenDisabledByDefault)
{
  turn_off_contact_search_before_parsing();
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--contact-search=ON"});

  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(), DefaultSettings::useContactSearch);
}

TEST_F(BalanceCommandLine, enableSearch_smDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm", "--contact-search=on"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, enableSearch_sdDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd", "--contact-search=on"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}


TEST_F(BalanceCommandLine, decompMethodRcb)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--decomp-method=rcb"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "rcb");
}

TEST_F(BalanceCommandLine, shortDecompMethodRcb)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"-m", "rcb"});

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

TEST_F(BalanceCommandLine, decompMethodScotch)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--decomp-method=scotch"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "scotch");
}

TEST_F(BalanceCommandLine, decompMethodParmetis)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--decomp-method=parmetis"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, vertexWeightMethodConnectivity)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--vertex-weight-method=connectivity"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             VertexWeightMethod::CONNECTIVITY);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               10.);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, vertexWeightMethodField)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--vertex-weight-method=field"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             VertexWeightMethod::FIELD);
  EXPECT_EQ(balanceSettings.getVertexWeightFieldName(),                          "");  // Not set.  Will throw on access.
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               1.);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, vertexWeightMethodFieldAndFieldName)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--vertex-weight-method=field",
                                                                                   "--vertex-weight-field-name=JunkField"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             VertexWeightMethod::FIELD);
  EXPECT_EQ(balanceSettings.getVertexWeightFieldName(),                          "JunkField");
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               1.);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, vertexWeightFieldName)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--vertex-weight-field-name=JunkField"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             VertexWeightMethod::FIELD);  // Auto-set
  EXPECT_EQ(balanceSettings.getVertexWeightFieldName(),                          "JunkField");
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               1.);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, vertexWeightFieldName_conflictingVertexWeightMethod)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--vertex-weight-method=connectivity",
                                                                                   "--vertex-weight-field-name=JunkField"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             VertexWeightMethod::FIELD);  // Overridden
  EXPECT_EQ(balanceSettings.getVertexWeightFieldName(),                          "JunkField");
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               1.);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {});
}

TEST_F(BalanceCommandLine, userSpecifiedBlockMultiplier_default)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings(
                                                         {"--block-weights=block_1:2.0"});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {{"block_1", 2.0}});
}

TEST_F(BalanceCommandLine, userSpecifiedBlockMultiplier_smDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sm",
                                                         {"--block-weights=block_2:10,block_5:2.5"}});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::smFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::smVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::smFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::smFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::smGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {{"block_2", 10.0}, {"block_5", 2.5}});
}

TEST_F(BalanceCommandLine, userSpecifiedWeights_sdDefaults)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings({"--sd",
                                                         {"--block-weights=block_1:1.5,block_2:5"}});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::sdFixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::sdVertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::sdFaceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::sdFaceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::sdGraphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {{"block_1", 1.5}, {"block_2", 5.0}});
}

TEST_F(BalanceCommandLine, userSpecifiedBlockMultiplier_badFormatting)
{
  const stk::balance::BalanceSettings& balanceSettings = get_stk_balance_settings(
                                                         {"--block-weights= block_1:1.5 ,block_2 : 3 , block_3:1.1 "});

  EXPECT_EQ(balanceSettings.getDecompMethod(),                                   DefaultSettings::decompMethod);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(),                             false);
  EXPECT_EQ(balanceSettings.includeSearchResultsInGraph(),                       DefaultSettings::useContactSearch);
  EXPECT_EQ(balanceSettings.shouldFixSpiders(),                                  DefaultSettings::fixSpiders);
  EXPECT_EQ(balanceSettings.shouldFixMechanisms(),                               false);
  EXPECT_EQ(balanceSettings.shouldPrintDiagnostics(),                            false);
  EXPECT_EQ(balanceSettings.getVertexWeightMethod(),                             (VertexWeightMethod)DefaultSettings::vertexWeightMethod);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(),                DefaultSettings::faceSearchEdgeWeight);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), DefaultSettings::faceSearchVertexMultiplier);
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightMultiplier(),               DefaultSettings::graphEdgeWeightMultiplier);
  check_relative_tolerance_for_face_search(balanceSettings,                      DefaultSettings::faceSearchRelTol);
  check_vertex_weight_block_multiplier(balanceSettings, {{"block_1", 1.5}, {"block_2", 3.0}, {"block_3", 1.1}});
}


