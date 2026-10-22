// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshInput_Parser.hpp>
#include <Akri_MeshInputOptions.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_Parser.hpp>

#include <stk_topology/topology.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <string>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

// Unit tests for MeshInput_Parser::parse_generated_mesh.  This is the richest
// self-contained function in the mesh-input parser: it takes a YAML fem_node and
// a caller-owned MeshInputOptions, so tests can drive every branch with a fresh
// MeshInputOptions instance (avoiding the name-keyed static registry used by
// MeshInput_Parser::parse).  RuntimeDoomedAdHoc increments a monotonic
// process-global counter rather than throwing; error paths capture
// stk::get_doomed_count() before the call and assert it grew.
class GeneratedMeshParserTest : public ::testing::Test
{
protected:
  Parser::Node load_node(const std::string & yamlText)
  {
    return Parser::Node(YAML::Load(yamlText));
  }
};

TEST_F(GeneratedMeshParserTest, NoGeneratedMeshMapReturnsFalse)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
mesh: foo.g
)yaml");

  EXPECT_FALSE(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, Parses3DDomainWithDefaults)
{
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
)yaml");

  EXPECT_TRUE(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_DOUBLE_EQ(0.5, options.get_generated_mesh_size());
  EXPECT_EQ(MeshInputOptions::GENERATED_MESH_FOR_SPECIFIED_DOMAIN, options.get_generated_mesh_domain_type());
  ASSERT_EQ(6u, options.get_generated_mesh_domain().size());
  EXPECT_EQ(3, options.get_generated_mesh_spatial_dimension());
  // Default element type for 3D is a tetrahedron.
  EXPECT_EQ(stk::topology::TETRAHEDRON_4, options.get_generated_mesh_element_type());
}

TEST_F(GeneratedMeshParserTest, Parses2DDomainDefaultsToTriangle)
{
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.25
  domain: [0.0, 0.0, 1.0, 1.0]
)yaml");

  EXPECT_TRUE(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_EQ(2, options.get_generated_mesh_spatial_dimension());
  EXPECT_EQ(stk::topology::TRIANGLE_3_2D, options.get_generated_mesh_element_type());
}

TEST_F(GeneratedMeshParserTest, MissingMeshSizeDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
)yaml");

  MeshInput_Parser::parse_generated_mesh(fem_node, options);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, NegativeMeshSizeDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: -1.0
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
)yaml");

  MeshInput_Parser::parse_generated_mesh(fem_node, options);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, NonNumericMeshSizeDoomsWithoutThrowing)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  // A non-numeric mesh_size scalar must be reported as a clean doom rather than
  // aborting the parse with a yaml-cpp conversion exception.
  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: not_a_number
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
)yaml");

  EXPECT_NO_THROW(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, BadDomainSizeDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 1.0]
)yaml");

  MeshInput_Parser::parse_generated_mesh(fem_node, options);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, ParsesInterfaceBoundingBox3D)
{
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  interface_bounding_box_with_dimension: 3
)yaml");

  EXPECT_TRUE(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_EQ(MeshInputOptions::GENERATED_3D_MESH_FOR_INTERFACE_BOUNDING_BOX, options.get_generated_mesh_domain_type());
}

TEST_F(GeneratedMeshParserTest, ParsesInterfaceBoundingBox2D)
{
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  interface_bounding_box_with_dimension: 2
)yaml");

  EXPECT_TRUE(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_EQ(MeshInputOptions::GENERATED_2D_MESH_FOR_INTERFACE_BOUNDING_BOX, options.get_generated_mesh_domain_type());
}

TEST_F(GeneratedMeshParserTest, BadInterfaceBoundingBoxDimensionDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  interface_bounding_box_with_dimension: 4
)yaml");

  MeshInput_Parser::parse_generated_mesh(fem_node, options);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, SpecifyingBothDomainAndInterfaceBoundingBoxDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
  interface_bounding_box_with_dimension: 3
)yaml");

  MeshInput_Parser::parse_generated_mesh(fem_node, options);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, ParsesExplicitElementTypeHex)
{
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
  element_type: HEX
)yaml");

  EXPECT_TRUE(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_EQ(stk::topology::HEXAHEDRON_8, options.get_generated_mesh_element_type());
}

TEST_F(GeneratedMeshParserTest, ElementTypeIsCaseInsensitive)
{
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
  element_type: tetrahedron
)yaml");

  EXPECT_TRUE(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_EQ(stk::topology::TETRAHEDRON_4, options.get_generated_mesh_element_type());
}

TEST_F(GeneratedMeshParserTest, InvalidElementTypeDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
  element_type: DODECAHEDRON
)yaml");

  MeshInput_Parser::parse_generated_mesh(fem_node, options);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, ElementTypeDimensionMismatchDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  // 3D domain but a 2D (triangle) element type -> dimension mismatch.
  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
  element_type: TRIANGLE
)yaml");

  MeshInput_Parser::parse_generated_mesh(fem_node, options);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, ParsesMeshTypeCubic)
{
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
  mesh_type: CUBIC
)yaml");

  EXPECT_TRUE(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_EQ(BoundingBoxMeshStructureType::CUBIC_BOUNDING_BOX_MESH, options.get_generated_mesh_structure_type());
}

TEST_F(GeneratedMeshParserTest, ParsesMeshTypeBccIn3D)
{
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
  mesh_type: BCC
)yaml");

  EXPECT_TRUE(MeshInput_Parser::parse_generated_mesh(fem_node, options));
  EXPECT_EQ(BoundingBoxMeshStructureType::BCC_BOUNDING_BOX_MESH, options.get_generated_mesh_structure_type());
}

TEST_F(GeneratedMeshParserTest, InvalidMeshTypeDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
  mesh_type: NOT_A_MESH_TYPE
)yaml");

  MeshInput_Parser::parse_generated_mesh(fem_node, options);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, BccMeshTypeIn2DDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  MeshInputOptions options("model");

  // BCC only valid in 3D; a 2D (length-4) domain must doom.
  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
  domain: [0.0, 0.0, 1.0, 1.0]
  mesh_type: BCC
)yaml");

  MeshInput_Parser::parse_generated_mesh(fem_node, options);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(GeneratedMeshParserTest, GeneratedMeshWithoutDomainOrInterfaceThrows)
{
  MeshInputOptions options("model");

  // generated_mesh present with a valid mesh_size but neither domain nor
  // interface_bounding_box_with_dimension -> domain type stays NO_GENERATED_MESH.
  // With no explicit element_type, the default-topology branch calls
  // get_generated_mesh_spatial_dimension(), which STK_ThrowRequire's
  // use_generated_mesh(); since the domain type never got set, that requirement
  // fails and throws before the "must specify domain..." doom at the end is
  // reached.  Document the current behavior: this surfaces as a throw.
  Parser::Node fem_node = load_node(R"yaml(
generated_mesh:
  mesh_size: 0.5
)yaml");

  EXPECT_ANY_THROW(MeshInput_Parser::parse_generated_mesh(fem_node, options));
}

#else

TEST(GeneratedMeshParserTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
