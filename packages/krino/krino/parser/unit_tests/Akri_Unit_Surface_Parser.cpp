/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Akri_Surface_Parser.hpp>

#include <Akri_AnalyticSurf.hpp>
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
// facets negative tests only. Positive tests would require actual files.
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

#else

TEST(SurfaceParserTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
