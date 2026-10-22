// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ResultsOutput_Parser.hpp>
#include <Akri_Region.hpp>
#include <Akri_Simulation.hpp>
#include <Akri_ResultsOutputOptions.hpp>
#include <Akri_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <string>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

// Unit tests for ResultsOutput_Parser::parse and
// ResultsSurfaceShellOutput_Parser::parse.  Both take a region_node and a
// Region &, writing into the region's ResultsOutputOptions.  Neither the
// Simulation ctor (timers + scalars) nor the Region ctor (StkMeshIoBroker over
// EnvData::parallel_comm() + a ResultsOutputOptions) requires a committed mesh,
// so the whole thing is unit-constructible; assertions read the options back
// via region.get_results_options().  RuntimeDoomedAdHoc increments a monotonic
// process-global counter rather than throwing; error paths capture
// stk::get_doomed_count() before the call and assert it grew.
//
// Ownership note: Region::Region calls my_simulation.add_region(this), which
// emplaces the raw pointer into the Simulation's std::vector<unique_ptr<Region>>.
// The Simulation therefore OWNS every Region and frees it in its own destructor,
// so a Region must be heap-allocated (never a stack object -- that double-frees).
// make_region() news up a Region owned by the fixture's Simulation, which the
// fixture destroys in TearDown.
class ResultsOutputParserTest : public ::testing::Test
{
protected:
  Parser::Node load_node(const std::string & yamlText)
  {
    return Parser::Node(YAML::Load(yamlText));
  }

  // Unique region/simulation name per test so nothing leaks between them.
  std::string unique_name() const
  {
    const ::testing::TestInfo * info = ::testing::UnitTest::GetInstance()->current_test_info();
    return std::string("ResultsOutputParserTest.") + info->name();
  }

  // Create a Region owned by the fixture's Simulation (see ownership note above).
  Region & make_region()
  {
    mSimulation = std::make_unique<Simulation>(unique_name());
    return *(new Region(*mSimulation, unique_name())); // owned by mSimulation
  }

  void TearDown() override { mSimulation.reset(); }

  std::unique_ptr<Simulation> mSimulation;
};

// -----------------------------------------------------------------------------
// ResultsOutput_Parser::parse
// -----------------------------------------------------------------------------

TEST_F(ResultsOutputParserTest, NoOutputBlockLeavesDefaultsAndDoesNotDoom)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
unrelated: 1
)yaml");

  ResultsOutput_Parser::parse(region_node, region);

  // Default filename set by ResultsOutputOptions ctor is untouched.
  EXPECT_EQ("default_filename", region.get_results_options()->get_filename());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, ParsesDatabaseName)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
output:
  database_name: results.e
)yaml");

  ResultsOutput_Parser::parse(region_node, region);

  EXPECT_EQ("results.e", region.get_results_options()->get_filename());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, MissingDatabaseNameDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  // An output block without database_name must doom (not throw).
  Parser::Node region_node = load_node(R"yaml(
output:
  title: My Title
)yaml");

  EXPECT_NO_THROW(ResultsOutput_Parser::parse(region_node, region));
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, ParsesTitle)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
output:
  database_name: results.e
  title: My Simulation Title
)yaml");

  ResultsOutput_Parser::parse(region_node, region);

  EXPECT_EQ("My Simulation Title", region.get_results_options()->get_title());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, ParsesComposeResultsTrueAddsProperty)
{
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
output:
  database_name: results.e
  compose_results: true
)yaml");

  ResultsOutput_Parser::parse(region_node, region);

  EXPECT_TRUE(region.get_results_options()->get_properties().exists("COMPOSE_RESULTS"));
}

TEST_F(ResultsOutputParserTest, ComposeResultsFalseAddsNoProperty)
{
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
output:
  database_name: results.e
  compose_results: false
)yaml");

  ResultsOutput_Parser::parse(region_node, region);

  EXPECT_FALSE(region.get_results_options()->get_properties().exists("COMPOSE_RESULTS"));
}

TEST_F(ResultsOutputParserTest, ParsesNodalVariables)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
output:
  database_name: results.e
  nodal_variables:
    - temperature
    - displacement
)yaml");

  ResultsOutput_Parser::parse(region_node, region);

  const auto & fields = region.get_results_options()->get_nodal_fields();
  ASSERT_EQ(2u, fields.size());
  // add_nodal_field(name, name) stores both internal and output name identically.
  EXPECT_EQ(1u, fields.count(FieldName_OutputName_Pair("temperature", "temperature")));
  EXPECT_EQ(1u, fields.count(FieldName_OutputName_Pair("displacement", "displacement")));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, ParsesElementVariables)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
output:
  database_name: results.e
  element_variables:
    - stress
)yaml");

  ResultsOutput_Parser::parse(region_node, region);

  const auto & fields = region.get_results_options()->get_element_fields();
  ASSERT_EQ(1u, fields.size());
  EXPECT_EQ(1u, fields.count(FieldName_OutputName_Pair("stress", "stress")));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, ParsesOutputFrequencyOntoVolumeScheduler)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
output:
  database_name: results.e
  output_frequency: 5
)yaml");

  ResultsOutput_Parser::parse(region_node, region);

  // Volume output_frequency drives the volume-output scheduler, not the
  // surface-shell one.
  EXPECT_EQ(5, region.get_results_options()->get_output_step_interval());
  EXPECT_EQ(0, region.get_results_options()->get_surface_shell_step_interval());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// ResultsSurfaceShellOutput_Parser::parse
// -----------------------------------------------------------------------------

TEST_F(ResultsOutputParserTest, NoSurfaceShellBlockDoesNotDoom)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
unrelated: 1
)yaml");

  ResultsSurfaceShellOutput_Parser::parse(region_node, region);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, ParsesSurfaceShellDatabaseName)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
surface_shell_output:
  database_name: shells.e
)yaml");

  ResultsSurfaceShellOutput_Parser::parse(region_node, region);

  EXPECT_EQ("shells.e", region.get_results_options()->get_surface_shells_filename());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, MissingSurfaceShellDatabaseNameDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
surface_shell_output:
  output_frequency: 2
)yaml");

  EXPECT_NO_THROW(ResultsSurfaceShellOutput_Parser::parse(region_node, region));
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, SurfaceShellOutputFrequencyDrivesSurfaceSchedulerNotVolume)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  // Regression: surface_shell_output's output_frequency must drive the
  // surface-shell scheduler, NOT the volume-output scheduler.
  Parser::Node region_node = load_node(R"yaml(
surface_shell_output:
  database_name: shells.e
  output_frequency: 7
)yaml");

  ResultsSurfaceShellOutput_Parser::parse(region_node, region);

  EXPECT_EQ(7, region.get_results_options()->get_surface_shell_step_interval());
  EXPECT_EQ(0, region.get_results_options()->get_output_step_interval());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ResultsOutputParserTest, ParsesSurfaceShellSurfaceNames)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Region & region = make_region();

  Parser::Node region_node = load_node(R"yaml(
surface_shell_output:
  database_name: shells.e
  surface_names:
    - surface_1
    - surface_2
)yaml");

  ResultsSurfaceShellOutput_Parser::parse(region_node, region);

  const auto & names = region.get_results_options()->get_surface_shells_surface_names();
  ASSERT_EQ(2u, names.size());
  EXPECT_EQ(1u, names.count("surface_1"));
  EXPECT_EQ(1u, names.count("surface_2"));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

#else

TEST(ResultsOutputParserTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
