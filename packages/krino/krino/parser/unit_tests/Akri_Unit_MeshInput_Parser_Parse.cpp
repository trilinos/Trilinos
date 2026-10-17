// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshInput_Parser.hpp>
#include <Akri_MeshInputOptions.hpp>
#include <Akri_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <string>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

// Unit tests for MeshInput_Parser::parse -- the registry-keyed entry point that
// walks a "finite_element_models" sequence, creates/looks up a MeshInputOptions
// per model via the process-global MeshInputOptions::get_or_create registry,
// delegates generated-mesh handling to parse_generated_mesh (covered separately
// in Akri_Unit_MeshInput_Parser.cpp), and forwards to Phase_Parser::parse.
// Because the registry is keyed by the model "name" from the YAML, each test
// uses a UNIQUE model name (derived from the gtest test name) so state cannot
// leak between tests.  RuntimeDoomedAdHoc increments a monotonic process-global
// counter rather than throwing; error paths capture stk::get_doomed_count()
// before the call and assert it grew.
class MeshInputParseTest : public ::testing::Test
{
protected:
  Parser::Node load_node(const std::string & yamlText)
  {
    return Parser::Node(YAML::Load(yamlText));
  }

  // Unique model name per test so the name-keyed registry stays isolated.
  std::string model_name() const
  {
    const ::testing::TestInfo * info = ::testing::UnitTest::GetInstance()->current_test_info();
    return std::string("MeshInputParseTest.") + info->name();
  }
};

TEST_F(MeshInputParseTest, NoFiniteElementModelsSequenceReturnsWithoutDooming)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node base_node = load_node(R"yaml(
unrelated: 1
)yaml");

  MeshInput_Parser::parse(base_node);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(MeshInputParseTest, ParsesInputMeshFilename)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string model = model_name();

  Parser::Node base_node = load_node(
      "finite_element_models:\n"
      "  - name: " + model + "\n"
      "    mesh: input.g\n");

  MeshInput_Parser::parse(base_node);

  MeshInputOptions * options = MeshInputOptions::get(model);
  ASSERT_NE(nullptr, options);
  EXPECT_EQ("input.g", options->get_filename());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(MeshInputParseTest, MissingModelNameDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  // No "name" -> empty model name dooms (registry key is the empty string).
  Parser::Node base_node = load_node(R"yaml(
finite_element_models:
  - mesh: input.g
)yaml");

  MeshInput_Parser::parse(base_node);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(MeshInputParseTest, NeitherMeshNorGeneratedMeshDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string model = model_name();

  // A named model with no mesh filename and no generated_mesh must doom.
  Parser::Node base_node = load_node(
      "finite_element_models:\n"
      "  - name: " + model + "\n");

  MeshInput_Parser::parse(base_node);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(MeshInputParseTest, ParsesDecompositionMethodUppercased)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string model = model_name();

  Parser::Node base_node = load_node(
      "finite_element_models:\n"
      "  - name: " + model + "\n"
      "    mesh: input.g\n"
      "    decomposition_method: rcb\n");

  MeshInput_Parser::parse(base_node);

  MeshInputOptions * options = MeshInputOptions::get(model);
  ASSERT_NE(nullptr, options);
  EXPECT_EQ("RCB", options->get_decomposition_method());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(MeshInputParseTest, ParsesUseAllSidesForShellsFalse)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string model = model_name();

  Parser::Node base_node = load_node(
      "finite_element_models:\n"
      "  - name: " + model + "\n"
      "    mesh: input.g\n"
      "    use_all_sides_for_shells: false\n");

  MeshInput_Parser::parse(base_node);

  MeshInputOptions * options = MeshInputOptions::get(model);
  ASSERT_NE(nullptr, options);
  EXPECT_FALSE(options->get_use_all_sides_for_shells());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(MeshInputParseTest, UseAllSidesForShellsDefaultsTrue)
{
  const std::string model = model_name();

  Parser::Node base_node = load_node(
      "finite_element_models:\n"
      "  - name: " + model + "\n"
      "    mesh: input.g\n");

  MeshInput_Parser::parse(base_node);

  MeshInputOptions * options = MeshInputOptions::get(model);
  ASSERT_NE(nullptr, options);
  EXPECT_TRUE(options->get_use_all_sides_for_shells());
}

TEST_F(MeshInputParseTest, ParsesGeneratedMeshThroughParse)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string model = model_name();

  // When generated_mesh is present, parse() takes the generated-mesh path and
  // does NOT require a mesh filename.  Verify the delegation populated options.
  Parser::Node base_node = load_node(
      "finite_element_models:\n"
      "  - name: " + model + "\n"
      "    generated_mesh:\n"
      "      mesh_size: 0.5\n"
      "      domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]\n");

  MeshInput_Parser::parse(base_node);

  MeshInputOptions * options = MeshInputOptions::get(model);
  ASSERT_NE(nullptr, options);
  // On the generated-mesh path set_filename is not called, so the filename
  // stays at its ctor default ("%B.g").
  EXPECT_EQ("%B.g", options->get_filename());
  EXPECT_DOUBLE_EQ(0.5, options->get_generated_mesh_size());
  EXPECT_EQ(MeshInputOptions::GENERATED_MESH_FOR_SPECIFIED_DOMAIN, options->get_generated_mesh_domain_type());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(MeshInputParseTest, GeneratedMeshTakesPrecedenceOverMeshFilename)
{
  const std::string model = model_name();

  // With both a generated_mesh and a mesh filename, has_generated_mesh is true
  // so set_filename is NOT called; the "mesh: ignored.g" entry is ignored and
  // the filename stays at its ctor default ("%B.g").
  Parser::Node base_node = load_node(
      "finite_element_models:\n"
      "  - name: " + model + "\n"
      "    mesh: ignored.g\n"
      "    generated_mesh:\n"
      "      mesh_size: 0.5\n"
      "      domain: [0.0, 0.0, 1.0, 1.0]\n");

  MeshInput_Parser::parse(base_node);

  MeshInputOptions * options = MeshInputOptions::get(model);
  ASSERT_NE(nullptr, options);
  EXPECT_EQ("%B.g", options->get_filename());
}

TEST_F(MeshInputParseTest, ParsesMultipleModels)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string modelA = model_name() + "_A";
  const std::string modelB = model_name() + "_B";

  Parser::Node base_node = load_node(
      "finite_element_models:\n"
      "  - name: " + modelA + "\n"
      "    mesh: a.g\n"
      "  - name: " + modelB + "\n"
      "    mesh: b.g\n");

  MeshInput_Parser::parse(base_node);

  MeshInputOptions * optionsA = MeshInputOptions::get(modelA);
  MeshInputOptions * optionsB = MeshInputOptions::get(modelB);
  ASSERT_NE(nullptr, optionsA);
  ASSERT_NE(nullptr, optionsB);
  EXPECT_EQ("a.g", optionsA->get_filename());
  EXPECT_EQ("b.g", optionsB->get_filename());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

#else

TEST(MeshInputParseTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
