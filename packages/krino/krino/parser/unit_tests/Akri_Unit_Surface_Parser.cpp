/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Akri_Surface_Parser.hpp>

#include <Akri_FacetFileSurface.hpp>
#include <Akri_MeshBlockSkinSurface.hpp>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_Parser.hpp>
#include <Akri_Surface.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <memory>
#include <string>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

class SurfaceParserTest : public ::testing::Test
{
protected:
  SurfaceParserTest()
    : meta(3)
  {
  }

  Parser::Node make_node(const std::string & yamlText) const
  {
    return Parser::Node(YAML::Load(yamlText));
  }

  std::shared_ptr<Surface> parse_surface(const std::string & yamlText) const
  {
    Parser::Node node = make_node(yamlText);
    return Surface_Parser::parse(node, meta, lsTimer);
  }

  template <class SurfaceType>
  std::shared_ptr<SurfaceType> parse_and_cast(const std::string & yamlText) const
  {
    std::shared_ptr<Surface> surface = parse_surface(yamlText);
    return std::dynamic_pointer_cast<SurfaceType>(surface);
  }

  stk::mesh::MetaData meta;

  stk::diag::Timer lsTimer{"LS Timer", sierra::Diag::sierraTimer()};
};

// -----------------------------------------------------------------------------
// Scalar syntax
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, ParsesScalarRandom)
{
  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
random
)yaml");

  ASSERT_NE(surface, nullptr);

  std::shared_ptr<Random> random =
      std::dynamic_pointer_cast<Random>(surface);

  EXPECT_NE(random, nullptr);
}

TEST_F(SurfaceParserTest, UnknownScalarReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
sphere
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// Sphere
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, ParsesNestedSphere)
{
  std::shared_ptr<Sphere> sphere = parse_and_cast<Sphere>(R"yaml(
sphere:
  center: [0.0, 0.0, 0.0]
  radius: 1.0
)yaml");

  EXPECT_NE(sphere, nullptr);
}

TEST_F(SurfaceParserTest, ParsesFlatSphere)
{
  std::shared_ptr<Sphere> sphere = parse_and_cast<Sphere>(R"yaml(
sphere:
center: [0.0, 0.0, 0.0]
radius: 1.0
)yaml");

  EXPECT_NE(sphere, nullptr);
}

TEST_F(SurfaceParserTest, SphereMissingCenterReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
sphere:
  radius: 1.0
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, SphereMissingRadiusReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
sphere:
  center: [0.0, 0.0, 0.0]
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, SphereBadCenterLengthReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
sphere:
  center: [0.0, 0.0]
  radius: 1.0
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, SphereNonPositiveRadiusReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
sphere:
  center: [0.0, 0.0, 0.0]
  radius: 0.0
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, SphereBadSignReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
sphere:
  center: [0.0, 0.0, 0.0]
  radius: 1.0
  sign: 2.0
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// Ellipsoid
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, ParsesNestedEllipsoid)
{
  std::shared_ptr<Ellipsoid> ellipsoid = parse_and_cast<Ellipsoid>(R"yaml(
ellipsoid:
  center: [0.0, 0.0, 0.0]
  semiaxes: [1.0, 2.0, 3.0]
  rotation: [0.0, 0.0, 0.0]
)yaml");

  EXPECT_NE(ellipsoid, nullptr);
}

TEST_F(SurfaceParserTest, ParsesFlatEllipsoid)
{
  std::shared_ptr<Ellipsoid> ellipsoid = parse_and_cast<Ellipsoid>(R"yaml(
ellipsoid:
center: [0.0, 0.0, 0.0]
semiaxes: [1.0, 2.0, 3.0]
rotation: [0.0, 0.0, 0.0]
)yaml");

  EXPECT_NE(ellipsoid, nullptr);
}

TEST_F(SurfaceParserTest, EllipsoidBadRotationLengthReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
ellipsoid:
  center: [0.0, 0.0, 0.0]
  semiaxes: [1.0, 2.0, 3.0]
  rotation: [0.0, 0.0]
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// Cuboid
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, ParsesNestedCuboid)
{
  std::shared_ptr<Cuboid> cuboid = parse_and_cast<Cuboid>(R"yaml(
cuboid:
  center: [0.0, 0.0, 0.0]
  dimensions: [10.0, 20.0, 30.0]
)yaml");

  EXPECT_NE(cuboid, nullptr);
}

TEST_F(SurfaceParserTest, ParsesNestedCuboidWithRotation)
{
  std::shared_ptr<Cuboid> cuboid = parse_and_cast<Cuboid>(R"yaml(
cuboid:
  center: [0.0, 0.0, 0.0]
  dimensions: [10.0, 20.0, 30.0]
  rotation: [0.0, 0.0, 0.0]
)yaml");

  EXPECT_NE(cuboid, nullptr);
}

TEST_F(SurfaceParserTest, ParsesFlatCuboid)
{
  std::shared_ptr<Cuboid> cuboid = parse_and_cast<Cuboid>(R"yaml(
cuboid:
center: [0.0, 0.0, 0.0]
dimensions: [10.0, 20.0, 30.0]
rotation: [0.0, 0.0, 0.0]
)yaml");

  EXPECT_NE(cuboid, nullptr);
}

TEST_F(SurfaceParserTest, CuboidMissingDimensionsReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
cuboid:
  center: [0.0, 0.0, 0.0]
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// Plane
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, ParsesNestedPlane)
{
  std::shared_ptr<Plane> plane = parse_and_cast<Plane>(R"yaml(
plane:
  normal: [0.0, 1.0, 0.0]
  offset: -0.1775
)yaml");

  EXPECT_NE(plane, nullptr);
}

TEST_F(SurfaceParserTest, ParsesFlatPlane)
{
  std::shared_ptr<Plane> plane = parse_and_cast<Plane>(R"yaml(
plane:
normal: [0.0, 1.0, 0.0]
offset: -0.1775
)yaml");

  EXPECT_NE(plane, nullptr);
}

TEST_F(SurfaceParserTest, ParsesFlatPlaneWithNameSibling)
{
  std::shared_ptr<Plane> plane = parse_and_cast<Plane>(R"yaml(
name: LS8
plane:
normal: [0.0, 1.0, 0.0]
offset: -0.1775
)yaml");

  EXPECT_NE(plane, nullptr);
}

TEST_F(SurfaceParserTest, PlaneMissingNormalReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
plane:
  offset: 1.0
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, PlaneMissingOffsetReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
plane:
  normal: [0.0, 1.0, 0.0]
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, PlaneBadSignReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
plane:
  normal: [0.0, 1.0, 0.0]
  offset: 0.0
  sign: 3.0
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, ParsesPlaneWithPoint)
{
  std::shared_ptr<Plane> plane = parse_and_cast<Plane>(R"yaml(
plane:
  normal: [0.0, 2.0, 0.0]
  point: [1.0, 0.25, -3.0]
)yaml");

  ASSERT_NE(plane, nullptr);
  const double tol = 1.e-9;
  EXPECT_NEAR(0., plane->point_signed_distance(stk::math::Vector3d(7.0, 0.25, 4.0)), tol);
  EXPECT_NEAR(0.75, plane->point_signed_distance(stk::math::Vector3d(0.0, 1.0, 0.0)), tol);
}

TEST_F(SurfaceParserTest, ParsesPlaneWithPointAndNegativeSign)
{
  std::shared_ptr<Plane> plane = parse_and_cast<Plane>(R"yaml(
plane:
  normal: [0.0, 1.0, 0.0]
  point: [0.0, 0.25, 0.0]
  sign: -1
)yaml");

  ASSERT_NE(plane, nullptr);
  const double tol = 1.e-9;
  EXPECT_NEAR(-0.75, plane->point_signed_distance(stk::math::Vector3d(0.0, 1.0, 0.0)), tol);
}

TEST_F(SurfaceParserTest, PlaneWithOffsetAndPointReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
plane:
  normal: [0.0, 1.0, 0.0]
  offset: 0.0
  point: [0.0, 0.0, 0.0]
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, PlaneBadPointLengthReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
plane:
  normal: [0.0, 1.0, 0.0]
  point: [0.0, 0.0]
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, PlaneZeroNormalReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
plane:
  normal: [0.0, 0.0, 0.0]
  offset: 1.0
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, PlaneNonFiniteValuesReturnNullAndDoom)
{
  for (const std::string & yamlText : {
      std::string("plane:\n  normal: [0.0, .inf, 0.0]\n  offset: 1.0\n"),
      std::string("plane:\n  normal: [0.0, 1.0, 0.0]\n  offset: .nan\n"),
      std::string("plane:\n  normal: [0.0, 1.0, 0.0]\n  point: [0.0, -.inf, 0.0]\n"),
      std::string("plane:\n  normal: [0.0, 1.0, 0.0]\n  offset: 1.0\n  multiplier: .inf\n")})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    EXPECT_EQ(parse_surface(yamlText), nullptr) << yamlText;
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
  }
}

TEST_F(SurfaceParserTest, ParsesPlaneWithFractionOfMeshExtent)
{
  for (const std::string & yamlText : {
      std::string("plane:\n  normal: [0.0, 2.0, 0.0]\n  fraction_of_mesh_extent: 0.25\n"),
      std::string("plane:\nnormal: [0.0, 2.0, 0.0]\nfraction_of_mesh_extent: 0.25\n")})
  {
    std::shared_ptr<Plane> plane = parse_and_cast<Plane>(yamlText);
    ASSERT_NE(plane, nullptr) << yamlText;
    EXPECT_TRUE(plane->needs_offset_from_extent());
    EXPECT_DOUBLE_EQ(0.25, plane->fraction_of_extent());

    // The plane is placed once the extent along its (unit) normal is known: here y in [-1, 3].
    plane->set_offset_from_extent(-1., 3.);
    EXPECT_FALSE(plane->needs_offset_from_extent());
    const double tol = 1.e-9;
    EXPECT_NEAR(0., plane->point_signed_distance(stk::math::Vector3d(5.0, 0.0, -2.0)), tol);
    EXPECT_NEAR(1., plane->point_signed_distance(stk::math::Vector3d(0.0, 1.0, 0.0)), tol);
  }
}

TEST_F(SurfaceParserTest, PlaneWithMoreThanOnePlacementReturnsNullAndDooms)
{
  for (const std::string & placementLines : {
      std::string("  offset: 0.0\n  fraction_of_mesh_extent: 0.5\n"),
      std::string("  point: [0.0, 0.0, 0.0]\n  fraction_of_mesh_extent: 0.5\n"),
      std::string("  offset: 0.0\n  point: [0.0, 0.0, 0.0]\n  fraction_of_mesh_extent: 0.5\n")})
  {
    const std::string yamlText = "plane:\n  normal: [0.0, 1.0, 0.0]\n" + placementLines;
    const unsigned oldDoomedCount = stk::get_doomed_count();
    EXPECT_EQ(parse_surface(yamlText), nullptr) << yamlText;
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
  }
}

TEST_F(SurfaceParserTest, PlaneNonFiniteFractionOfMeshExtentReturnsNullAndDooms)
{
  for (const std::string & value : {std::string(".nan"), std::string(".inf"), std::string("-.inf")})
  {
    const std::string yamlText = "plane:\n  normal: [0.0, 1.0, 0.0]\n  fraction_of_mesh_extent: " + value + "\n";
    const unsigned oldDoomedCount = stk::get_doomed_count();
    EXPECT_EQ(parse_surface(yamlText), nullptr) << yamlText;
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
  }
}

// -----------------------------------------------------------------------------
// Cylinder
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, ParsesNestedCylinder)
{
  std::shared_ptr<Cylinder> cylinder = parse_and_cast<Cylinder>(R"yaml(
cylinder:
  p1: [0.0, 0.0, -1.0]
  p2: [0.0, 0.0,  1.0]
  radius: 0.25
)yaml");

  EXPECT_NE(cylinder, nullptr);
}

TEST_F(SurfaceParserTest, ParsesFlatCylinder)
{
  std::shared_ptr<Cylinder> cylinder = parse_and_cast<Cylinder>(R"yaml(
cylinder:
p1: [0.0, 0.0, -1.0]
p2: [0.0, 0.0,  1.0]
radius: 0.25
)yaml");

  EXPECT_NE(cylinder, nullptr);
}

TEST_F(SurfaceParserTest, CylinderMissingP1ReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
cylinder:
  p2: [0.0, 0.0, 1.0]
  radius: 0.25
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, CylinderNonPositiveRadiusReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
cylinder:
  p1: [0.0, 0.0, -1.0]
  p2: [0.0, 0.0,  1.0]
  radius: -1.0
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// string_function
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, ParsesNestedStringFunction)
{
  std::shared_ptr<LevelSet_String_Function> surface =
      parse_and_cast<LevelSet_String_Function>(R"yaml(
string_function:
  expression: "x*x + y*y + z*z - 1.0"
)yaml");

  EXPECT_NE(surface, nullptr);
}

TEST_F(SurfaceParserTest, ParsesStringFunctionWithBoundingBox)
{
  std::shared_ptr<LevelSet_String_Function> surface =
      parse_and_cast<LevelSet_String_Function>(R"yaml(
string_function:
  expression: "x*x + y*y + z*z - 1.0"
  bounding_box: [-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]
)yaml");

  EXPECT_NE(surface, nullptr);
}

TEST_F(SurfaceParserTest, StringFunctionMissingExpressionReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
string_function:
  bounding_box: [-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, StringFunctionBadBoundingBoxLengthReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
string_function:
  expression: "x"
  bounding_box: [-1.0, -1.0, -1.0]
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// random map syntax
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, ParsesNestedRandomWithSeed)
{
  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
random:
  seed: 123
)yaml");

  ASSERT_NE(surface, nullptr);

  std::shared_ptr<Random> random =
      std::dynamic_pointer_cast<Random>(surface);

  EXPECT_NE(random, nullptr);
}

TEST_F(SurfaceParserTest, ParsesFlatRandomWithSeed)
{
  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
random:
seed: 123
)yaml");

  ASSERT_NE(surface, nullptr);

  std::shared_ptr<Random> random =
      std::dynamic_pointer_cast<Random>(surface);

  EXPECT_NE(random, nullptr);
}

// -----------------------------------------------------------------------------
// Surface marker errors
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, InvalidSurfaceMarkerValueReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
sphere: invalid_marker_value
center: [0.0, 0.0, 0.0]
radius: 1.0
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, UnknownMapReturnsNullWithoutDooming)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
not_a_surface:
  value: 1
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// facets
// -----------------------------------------------------------------------------

TEST_F(SurfaceParserTest, FacetsMissingFilenameReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
facets:
  format: STL
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(SurfaceParserTest, FacetsUnknownFormatReturnsNullAndDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  std::shared_ptr<Surface> surface = parse_surface(R"yaml(
facets:
  filename: dummy.file
  format: BAD_FORMAT
)yaml");

  EXPECT_EQ(surface, nullptr);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// Facet files are read into a connected surface, which reads its file when first prepared, so these
// files need not exist.
TEST_F(SurfaceParserTest, FacetsOfEveryFormatAreConnected)
{
  for (const std::string format : {"STL", "FAC", "PLY", "EXO", "stl"})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    std::shared_ptr<Surface> surface =
        parse_surface("facets:\n  filename: dummy.file\n  format: " + format + "\n");
    EXPECT_NE(nullptr, std::dynamic_pointer_cast<FacetFileSurface>(surface)) << format;
    EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount) << format;
  }
}

// -----------------------------------------------------------------------------
// mesh_block_skin
// -----------------------------------------------------------------------------

// The file is read when the surface is first prepared, so it need not exist to parse.
TEST_F(SurfaceParserTest, MeshBlockSkinWithBlockOrBlocks_isAMeshBlockSkinSurface)
{
  for (const std::string blocks : {"  block: block_1\n", "  blocks: [block_1, block_2]\n"})
  {
    auto surface = parse_and_cast<MeshBlockSkinSurface>(
        "mesh_block_skin:\n  filename: parts.g\n" + blocks + "  sign: -1\n  sharp_feature_angle: 135\n");
    EXPECT_NE(nullptr, surface) << blocks;
  }
}

TEST_F(SurfaceParserTest, MeshBlockSkinWithoutFilenameOrBlocks_returnsNullAndDooms)
{
  for (const std::string yaml : {"mesh_block_skin:\n  block: block_1\n",
           "mesh_block_skin:\n  filename: parts.g\n",
           "mesh_block_skin:\n  filename: parts.g\n  block: block_1\n  blocks: [block_2]\n",
           "mesh_block_skin:\n  filename: parts.g\n  block: block_1\n  sign: 2\n"})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    EXPECT_EQ(nullptr, parse_surface(yaml)) << yaml;
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yaml;
  }
}

TEST_F(SurfaceParserTest, FacetsConnectedWithRemoveCoincidentAndWatertightOptions)
{
  for (const std::string merge : {"true", "false"})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    std::shared_ptr<Surface> surface = parse_surface(R"yaml(
facets:
  filename: dummy.ply
  format: PLY
  remove_coincident_facets: true
  check_watertight: true
  watertight_bounding_box: [0, 0, 0, 1, 1, 1]
  merge_coincident_vertices: )yaml" +
        merge + "\n");
    EXPECT_NE(nullptr, std::dynamic_pointer_cast<FacetFileSurface>(surface)) << merge;
    EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount) << merge;
  }
}

// Sharp features are off unless sharp_feature_angle is given, with or without open boundary edges
// as features.
TEST_F(SurfaceParserTest, FacetsWithSharpFeatureAngle_haveSharpFeatures)
{
  const std::string facets = "facets:\n  filename: dummy.ply\n  format: PLY\n";
  for (const std::string options : {"",
           "  sharp_feature_angle: 135\n",
           "  sharp_feature_angle: 135\n  open_boundary_edges_are_sharp_features: true\n"})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    const auto surface =
        std::dynamic_pointer_cast<FacetFileSurface>(parse_surface(facets + options));
    ASSERT_NE(nullptr, surface) << options;
    EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount) << options;
    const bool hasSharpFeatures = !options.empty();
    EXPECT_EQ(hasSharpFeatures, surface->has_edges()) << options;
    EXPECT_EQ(hasSharpFeatures, surface->has_corners()) << options;
  }
}

// A sharp feature angle out of range, open boundary edges as features without an angle, and sharp
// features without merged vertices (an unconnected surface) are input errors.
TEST_F(SurfaceParserTest, FacetsWithInvalidSharpFeatureOptions_returnsNullAndDooms)
{
  const std::string facets = "facets:\n  filename: dummy.ply\n  format: PLY\n";
  for (const std::string options : {"  sharp_feature_angle: 0\n",
           "  sharp_feature_angle: 180\n",
           "  sharp_feature_angle: -30\n",
           "  open_boundary_edges_are_sharp_features: true\n",
           "  sharp_feature_angle: 135\n  merge_coincident_vertices: false\n"})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    EXPECT_EQ(nullptr, parse_surface(facets + options)) << options;
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << options;
  }
}

// -----------------------------------------------------------------------------
// motion: a surface with a motion is a TransformedSurface placed by it
// -----------------------------------------------------------------------------

const BoundingBox anyBox(stk::math::Vector3d(-10.,-10.,-10.), stk::math::Vector3d(10.,10.,10.));

TEST_F(SurfaceParserTest, SphereWithMotion_placedByItsMotion)
{
  std::shared_ptr<Surface> sphere = parse_surface(R"yaml(
sphere:
  center: [1.0, 0.0, 0.0]
  radius: 0.5
  motion:
    translational_velocity: [0.0, 2.0, 0.0]
    initial_displacement: [0.0, 0.0, 1.0]
)yaml");

  ASSERT_NE(sphere, nullptr);
  EXPECT_EQ(TRANSFORMED_SURFACE, sphere->type());
  sphere->prepare_to_compute(0., anyBox, 0.);
  EXPECT_DOUBLE_EQ(-0.5, sphere->point_signed_distance(stk::math::Vector3d(1., 0., 1.)));
  sphere->prepare_to_compute(1.5, anyBox, 0.);
  EXPECT_DOUBLE_EQ(-0.5, sphere->point_signed_distance(stk::math::Vector3d(1., 3., 1.)));
}

// Flat syntax, and a rotation about a reference point: a quarter turn per unit time about z
// through (1,1,0), after an initial quarter turn the other way, with the reference point moving.
TEST_F(SurfaceParserTest, FlatCuboidWithEveryMotionKey_placedByItsMotion)
{
  std::shared_ptr<Surface> cuboid = parse_surface(R"yaml(
cuboid:
center: [2.0, 1.0, 0.0]
dimensions: [0.2, 0.2, 0.2]
motion:
  reference_point: [1.0, 1.0, 0.0]
  initial_rotation: [0.0, 0.0, -1.5707963267948966]
  initial_displacement: [0.0, 0.0, 3.0]
  rotational_velocity: [0.0, 0.0, 1.5707963267948966]
  translational_velocity: [1.0, 0.0, 0.0]
)yaml");

  ASSERT_NE(cuboid, nullptr);
  EXPECT_EQ(TRANSFORMED_SURFACE, cuboid->type());
  cuboid->prepare_to_compute(0., anyBox, 0.);
  EXPECT_NEAR(-0.1, cuboid->point_signed_distance(stk::math::Vector3d(1., 0., 3.)), 1.e-12);
  // At time 2: turned a quarter turn the first way in all, about the reference point now at (3,1,3).
  cuboid->prepare_to_compute(2., anyBox, 0.);
  EXPECT_NEAR(-0.1, cuboid->point_signed_distance(stk::math::Vector3d(3., 2., 3.)), 1.e-12);
}

TEST_F(SurfaceParserTest, MotionWithTwoComponentVectors_zIsZero)
{
  std::shared_ptr<Surface> sphere = parse_surface(R"yaml(
sphere:
  center: [0.0, 0.0, 0.0]
  radius: 0.5
  motion:
    translational_velocity: [0.9, 0.0]
    initial_displacement: [0.0, 1.0]
)yaml");

  ASSERT_NE(sphere, nullptr);
  sphere->prepare_to_compute(1., anyBox, 0.);
  EXPECT_DOUBLE_EQ(-0.5, sphere->point_signed_distance(stk::math::Vector3d(0.9, 1., 0.)));
}

TEST_F(SurfaceParserTest, SurfaceWithoutMotion_notTransformed)
{
  std::shared_ptr<Surface> sphere = parse_surface(R"yaml(
sphere:
  center: [0.0, 0.0, 0.0]
  radius: 0.5
)yaml");
  ASSERT_NE(sphere, nullptr);
  EXPECT_EQ(SPHERE, sphere->type());
}

TEST_F(SurfaceParserTest, BadMotion_returnsNullAndDooms)
{
  for (const std::string motion : {
           "  motion:\n    translational_velocity: [1.0]\n",
           "  motion:\n    rotational_velocity: [0.0, 1.0]\n",
           "  motion:\n    translation_velocity: [1.0, 0.0, 0.0]\n",
           "  motion: fast\n"})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    std::shared_ptr<Surface> surface = parse_surface("sphere:\n  center: [0.0, 0.0, 0.0]\n  radius: 0.5\n" + motion);
    EXPECT_EQ(surface, nullptr) << motion;
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << motion;
  }
}

#else

TEST(SurfaceParserTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
