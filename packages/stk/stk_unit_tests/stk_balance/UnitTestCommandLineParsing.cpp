#include <iosfwd>
#include <gtest/gtest.h>
#include <stk_balance/internal/balanceCommandLine.hpp>
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <fstream>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace
{

TEST(BalanceOutputDirectory, checkOutputFileFromInputFileWithoutPath)
{
  const std::string inputFile = "input.e";
  const std::string outputDirectory = "/home/code/results";

  std::string outputFilename = stk::balance::construct_output_file_name(outputDirectory, inputFile);
  EXPECT_EQ("/home/code/results/input.e", outputFilename);
}

TEST(BalanceOutputDirectory, checkOutputFileFromInputFileAndOutputDirectoryWithExtraBackslash)
{
  const std::string inputFile = "input.e";
  const std::string outputDirectory = "/home/code/results/";

  std::string outputFilename = stk::balance::construct_output_file_name(outputDirectory, inputFile);
  EXPECT_EQ("/home/code/results//input.e", outputFilename);
}

TEST(BalanceOutputDirectory, checkOutputFileFromInputFileWithPath)
{
  const std::string inputFile = "/another/directory/input.e";
  const std::string outputDirectory = "/home/code/results";

  std::string outputFilename = stk::balance::construct_output_file_name(outputDirectory, inputFile);
  EXPECT_EQ("/home/code/results/input.e", outputFilename);
}

class BalanceCommandLine : public ::testing::Test
{
protected:
    BalanceCommandLine() { }

    std::vector<const char *> get_argv_c_strings()
    {
        std::vector<const char *> argvCstr(argv.size());
        for(size_t i=0; i<argv.size(); i++)
            argvCstr[i] = argv[i].c_str();
        return argvCstr;
    }

    void set_argv(const std::vector<std::string>& inputArgv)
    {
        argv = inputArgv;
    }

    stk::balance::StkBalanceSettings get_stk_balance_settings(const std::vector<std::string> & options)
    {
        std::vector<std::string> args {"exe", "infile", "foobar"};
        args.insert(args.end(), options.begin(), options.end());
        set_argv(args);
        std::vector<const char *> argvCstr = get_argv_c_strings();
        stk::balance::ParsedOptions parsedOptions = stk::balance::parse_balance_command_line(argvCstr.size(),
                                                                                             argvCstr.data(),
                                                                                             argv[0],
                                                                                             MPI_COMM_WORLD);

        return stk::balance::create_balance_settings(parsedOptions);
    }

    double getAbsoluteToleranceForUnitMesh(const stk::balance::BalanceSettings & balanceSettings)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::io::fill_mesh("generated:1x1x4|sideset:x", bulk);
        stk::mesh::EntityVector faces;
        stk::mesh::get_selected_entities(meta.universal_part(), bulk.buckets(stk::topology::FACE_RANK), faces);
        ThrowRequire(faces.size() > 0);

        return balanceSettings.getToleranceForFaceSearch(bulk, *meta.coordinate_field(),
                                                         bulk.begin_nodes(faces[0]), bulk.num_nodes(faces[0]));
    }

    void checkAbsoluteToleranceForFaceSearch(const stk::balance::BalanceSettings & balanceSettings, double tolerance)
    {
        EXPECT_TRUE(balanceSettings.isConstantFaceSearchTolerance());
        EXPECT_DOUBLE_EQ(getAbsoluteToleranceForUnitMesh(balanceSettings), tolerance);
    }

    void checkRelativeToleranceForFaceSearch(const stk::balance::BalanceSettings & balanceSettings, double tolerance)
    {
        EXPECT_FALSE(balanceSettings.isConstantFaceSearchTolerance());
        EXPECT_DOUBLE_EQ(getAbsoluteToleranceForUnitMesh(balanceSettings), tolerance);
    }

    std::vector<std::string> argv;
};

TEST_F(BalanceCommandLine, createBalanceSettings_default)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sm"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  checkRelativeToleranceForFaceSearch(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sd"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}


TEST_F(BalanceCommandLine, createBalanceSettings_defaultAbsoluteTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--face-search-abs-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_defaultRelativeTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--face-search-rel-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkRelativeToleranceForFaceSearch(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, createBalanceSettings_bothTolerances)
{
  EXPECT_THROW(get_stk_balance_settings({"--face-search-abs-tol", "--face-search-rel-tol"}), std::logic_error);
}


TEST_F(BalanceCommandLine, createBalanceSettings_faceSearchAbsoluteTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--face-search-abs-tol=0.001"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_faceSearchRelativeTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--face-search-rel-tol=0.123"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkRelativeToleranceForFaceSearch(balanceSettings, 0.123);
}


TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_defaultAbsoluteTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sm",
                                                                               "--face-search-abs-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_defaultRelativeTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sm",
                                                                               "--face-search-rel-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  checkRelativeToleranceForFaceSearch(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_faceSearchAbsoluteTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sm",
                                                                               "--face-search-abs-tol=0.005"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.005);
}

TEST_F(BalanceCommandLine, createBalanceSettings_smDefaults_faceSearchRelativeTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sm",
                                                                               "--face-search-rel-tol=0.123"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  checkRelativeToleranceForFaceSearch(balanceSettings, 0.123);
}


TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_defaultAbsoluteTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sd",
                                                                               "--face-search-abs-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_defaultRelativeTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sd",
                                                                               "--face-search-rel-tol"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkRelativeToleranceForFaceSearch(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_faceSearchAbsoluteTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sd",
                                                                               "--face-search-abs-tol=0.0005"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0005);
}

TEST_F(BalanceCommandLine, createBalanceSettings_sdDefaults_faceSearchRelativeTolerance)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sd",
                                                                               "--face-search-rel-tol=0.123"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkRelativeToleranceForFaceSearch(balanceSettings, 0.123);
}


TEST_F(BalanceCommandLine, contactSearch_badValue)
{
  EXPECT_THROW(get_stk_balance_settings({"--contact-search=junk"}), std::runtime_error);
}

TEST_F(BalanceCommandLine, disableSearch_default_caseInsensitive)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--contact-search=OFF"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_FALSE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, disableSearch_smDefaults)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sm", "--contact-search=off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_FALSE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  checkRelativeToleranceForFaceSearch(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, disableSearch_sdDefaults)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sd", "--contact-search=off"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_FALSE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}


TEST_F(BalanceCommandLine, enableSearch_default_caseInsensitive)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--contact-search=ON"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}

TEST_F(BalanceCommandLine, enableSearch_smDefaults)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sm", "--contact-search=on"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 3.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 10.0);
  checkRelativeToleranceForFaceSearch(balanceSettings, 0.15);
}

TEST_F(BalanceCommandLine, enableSearch_sdDefaults)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--sd", "--contact-search=on"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_TRUE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}


TEST_F(BalanceCommandLine, decompMethodRcb)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--decomp-method=rcb"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "rcb");
}

TEST_F(BalanceCommandLine, decompMethodRib)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--decomp-method=rib"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "rib");
}

TEST_F(BalanceCommandLine, decompMethodMultijagged)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--decomp-method=multijagged"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "multijagged");
}

TEST_F(BalanceCommandLine, decompMethodParmetis)
{
  stk::balance::StkBalanceSettings balanceSettings = get_stk_balance_settings({"--decomp-method=parmetis"});

  EXPECT_EQ(balanceSettings.getDecompMethod(), "parmetis");
  EXPECT_TRUE(balanceSettings.includeSearchResultsInGraph());
  EXPECT_FALSE(balanceSettings.shouldFixSpiders());
  EXPECT_DOUBLE_EQ(balanceSettings.getGraphEdgeWeightForSearch(), 15.0);
  EXPECT_DOUBLE_EQ(balanceSettings.getVertexWeightMultiplierForVertexInSearch(), 5.0);
  checkAbsoluteToleranceForFaceSearch(balanceSettings, 0.0001);
}


TEST_F(BalanceCommandLine, defaultContactSearchStatus)
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

TEST_F(BalanceCommandLine, toggleContactSearchStatus)
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

class BalanceCommandLineWithOutputDir : public BalanceCommandLine
{
protected:
    BalanceCommandLineWithOutputDir()
    {
        argv.push_back("exe");
        argv.push_back("infile");
        argv.push_back("outputDir");
    }

    stk::balance::ParsedOptions test_command_line(stk::balance::AppTypeDefaults expectedAppDefaults)
    {
        const std::string quickExample = "";

        MPI_Comm comm(MPI_COMM_WORLD);

        std::vector<const char *> argvCstr = get_argv_c_strings();
        stk::balance::ParsedOptions balanceOptions = stk::balance::parse_balance_command_line(argvCstr.size(),
                                                                                              argvCstr.data(),
                                                                                              argv[0], comm);

        EXPECT_EQ(expectedAppDefaults, balanceOptions.appTypeDefaults);
        EXPECT_EQ(argv[1], balanceOptions.m_inFile);
        EXPECT_EQ(argv[2], balanceOptions.outputDirectory);
        return balanceOptions;
    }

    void add_defaults_flag(const std::string &defaultsName)
    {
        std::string defaultsFlag = stk::dash_it(defaultsName);
        argv.push_back(defaultsFlag.c_str());
    }
};

TEST_F(BalanceCommandLineWithOutputDir, parseSmDefaultOption_smDefaultSet)
{
    add_defaults_flag(stk::balance::CommandLineOptions().smDefaults.name);
    test_command_line(stk::balance::SM_DEFAULTS);
}

}
