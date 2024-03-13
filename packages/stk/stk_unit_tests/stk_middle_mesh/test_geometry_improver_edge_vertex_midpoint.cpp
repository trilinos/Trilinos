#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/geometry_improver_edge_vertex_midpoints.hpp"
#include "stk_middle_mesh/nonconformal4.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

class GeometryImproverEdgeVertexMidpointsTester : public ::testing::Test
{
  protected:
    void setup(mesh::impl::MeshSpec meshspec1, mesh::impl::MeshSpec meshspec2)
    {
      mesh1 = mesh::impl::create_mesh(meshspec1, [](const utils::Point& pt) { return pt; });
      mesh2 = mesh::impl::create_mesh(meshspec2, [](const utils::Point& pt) { return utils::Point(pt.x, pt.y, 1); });

      double eps = 1e-12;
      NormalProjectionOpts opts;
      opts.classifierTolerances = PointClassifierNormalWrapperTolerances(eps);
      opts.edgeTracerTolerances = middle_mesh::impl::EdgeTracerTolerances(eps);
      opts.geometryImprovers    = {nonconformal4::impl::GeometryImprovers::EdgeVertexMidpoints};
      nonconformal4::impl::Nonconformal4 maker(mesh1, mesh2, opts);

      meshIn = maker.create();
    }

    std::shared_ptr<mesh::Mesh> mesh1;
    std::shared_ptr<mesh::Mesh> mesh2;
    std::shared_ptr<mesh::Mesh> meshIn;
};

mesh::MeshEntityPtr get_closest_vert(std::shared_ptr<mesh::Mesh> mesh, const utils::Point& pt)
{
  double minDist              = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr minVert = nullptr;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      utils::Point ptOrig = vert->get_point_orig(0);
      utils::Point ptXy(ptOrig.x, ptOrig.y, 0);

      double dist = dot(pt - ptXy, pt - ptXy);
      if (dist < minDist)
      {
        minDist = dist;
        minVert = vert;
      }
    }

  return minVert;
}

void test_z_value(std::shared_ptr<mesh::Mesh> mesh, const utils::Point& vertLocation, double expectedZ)
{
  mesh::MeshEntityPtr vert = get_closest_vert(mesh, vertLocation);
  EXPECT_NEAR(vert->get_point_orig(0).z, expectedZ, 1e-13);
}

} // namespace

TEST_F(GeometryImproverEdgeVertexMidpointsTester, TwoToFour)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec meshspec1, meshspec2;
  meshspec1.numelX = 2;
  meshspec2.numelX = 4;
  meshspec1.numelY = 2;
  meshspec2.numelY = 4;
  meshspec1.xmin   = 0;
  meshspec2.xmin   = 0;
  meshspec1.xmax   = 1;
  meshspec2.xmax   = 1;
  meshspec1.ymin   = 0;
  meshspec2.ymin   = 0;
  meshspec1.ymax   = 1;
  meshspec2.ymax   = 1;

  setup(meshspec1, meshspec2);

  double deltaX = 1.0 / 4;
  double deltaY = 1.0 / 4;

  for (int i = 0; i < meshspec2.numelX + 1; ++i)
    for (int j = 0; j < meshspec2.numelY + 1; ++j)
    {
      double x                   = deltaX * i;
      double y                   = deltaY * j;
      mesh::MeshEntityPtr vertIn = get_closest_vert(meshIn, utils::Point(x, y, 0));
      EXPECT_NEAR(vertIn->get_point_orig(0).z, 0, 1e-13);
    }
}

TEST_F(GeometryImproverEdgeVertexMidpointsTester, TwoToFive)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec meshspec1, meshspec2;
  meshspec1.numelX = 2;
  meshspec2.numelX = 5;
  meshspec1.numelY = 2;
  meshspec2.numelY = 5;
  meshspec1.xmin   = 0;
  meshspec2.xmin   = 0;
  meshspec1.xmax   = 1;
  meshspec2.xmax   = 1;
  meshspec1.ymin   = 0;
  meshspec2.ymin   = 0;
  meshspec1.ymax   = 1;
  meshspec2.ymax   = 1;

  setup(meshspec1, meshspec2);

  test_z_value(meshIn, utils::Point(0, 0), 0);
  test_z_value(meshIn, utils::Point(0.2, 0), 0);
  test_z_value(meshIn, utils::Point(0.4, 0), 0);
  test_z_value(meshIn, utils::Point(0.5, 0), 0);
  test_z_value(meshIn, utils::Point(0.6, 0), 0);
  test_z_value(meshIn, utils::Point(0.8, 0), 0);
  test_z_value(meshIn, utils::Point(1, 0), 0);

  test_z_value(meshIn, utils::Point(0, 0.2), 0);
  test_z_value(meshIn, utils::Point(0.2, 0.2), 0);
  test_z_value(meshIn, utils::Point(0.4, 0.2), 0);
  test_z_value(meshIn, utils::Point(0.5, 0.2), 0.5);
  test_z_value(meshIn, utils::Point(0.6, 0.2), 0);
  test_z_value(meshIn, utils::Point(0.8, 0.2), 0);
  test_z_value(meshIn, utils::Point(1, 0.2), 0);

  test_z_value(meshIn, utils::Point(0, 0.4), 0);
  test_z_value(meshIn, utils::Point(0.2, 0.4), 0);
  test_z_value(meshIn, utils::Point(0.4, 0.4), 0);
  test_z_value(meshIn, utils::Point(0.5, 0.4), 0.5);
  test_z_value(meshIn, utils::Point(0.6, 0.4), 0);
  test_z_value(meshIn, utils::Point(0.8, 0.4), 0);
  test_z_value(meshIn, utils::Point(1, 0.4), 0);

  test_z_value(meshIn, utils::Point(0, 0.5), 0);
  test_z_value(meshIn, utils::Point(0.2, 0.5), 0.5);
  test_z_value(meshIn, utils::Point(0.4, 0.5), 0.5);
  test_z_value(meshIn, utils::Point(0.5, 0.5), 0);
  test_z_value(meshIn, utils::Point(0.6, 0.5), 0.5);
  test_z_value(meshIn, utils::Point(0.8, 0.5), 0.5);
  test_z_value(meshIn, utils::Point(1, 0.5), 0);

  test_z_value(meshIn, utils::Point(0, 0.6), 0);
  test_z_value(meshIn, utils::Point(0.2, 0.6), 0);
  test_z_value(meshIn, utils::Point(0.4, 0.6), 0);
  test_z_value(meshIn, utils::Point(0.5, 0.6), 0.5);
  test_z_value(meshIn, utils::Point(0.6, 0.6), 0);
  test_z_value(meshIn, utils::Point(0.8, 0.6), 0);
  test_z_value(meshIn, utils::Point(1, 0.6), 0);

  test_z_value(meshIn, utils::Point(0, 0.8), 0);
  test_z_value(meshIn, utils::Point(0.2, 0.8), 0);
  test_z_value(meshIn, utils::Point(0.4, 0.8), 0);
  test_z_value(meshIn, utils::Point(0.5, 0.8), 0.5);
  test_z_value(meshIn, utils::Point(0.6, 0.8), 0);
  test_z_value(meshIn, utils::Point(0.8, 0.8), 0);
  test_z_value(meshIn, utils::Point(1, 0.8), 0);

  test_z_value(meshIn, utils::Point(0, 1.0), 0);
  test_z_value(meshIn, utils::Point(0.2, 1.0), 0);
  test_z_value(meshIn, utils::Point(0.4, 1.0), 0);
  test_z_value(meshIn, utils::Point(0.5, 1.0), 0);
  test_z_value(meshIn, utils::Point(0.6, 1.0), 0);
  test_z_value(meshIn, utils::Point(0.8, 1.0), 0);
  test_z_value(meshIn, utils::Point(1, 1.0), 0);
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
