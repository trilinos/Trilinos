#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_io.hpp"
#include "stk_middle_mesh/nonconformal4.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

class GeometryImproverEdgeVertexTaylorPatchTester : public ::testing::Test
{
  protected:
    void setup(mesh::impl::MeshSpec meshspec1, mesh::impl::MeshSpec meshspec2)
    {
      auto f = [&](const utils::Point& pt) { return zmapping(pt); };
      mesh1  = mesh::impl::create_mesh(meshspec1, f);
      mesh2  = mesh::impl::create_mesh(meshspec2, f);

      double eps = 1e-12;
      NormalProjectionOpts opts;
      opts.classifierTolerances = PointClassifierNormalWrapperTolerances(eps);
      opts.edgeTracerTolerances = middle_mesh::impl::EdgeTracerTolerances(eps);
      opts.geometryImprovers    = {nonconformal4::impl::GeometryImprovers::RestoreMesh2Verts,
                                   nonconformal4::impl::GeometryImprovers::EdgeVertexTaylorPatchQuadraticZonly};
      nonconformal4::impl::Nonconformal4 maker(mesh1, mesh2, opts);

      meshIn = maker.create();
    }

    utils::Point zmapping(const utils::Point& pt) { return utils::Point(pt.x, pt.y, pt.x + pt.y + 0.5 * pt.x * pt.y); }

    std::shared_ptr<mesh::Mesh> mesh1;
    std::shared_ptr<mesh::Mesh> mesh2;
    std::shared_ptr<mesh::Mesh> meshIn;
};

} // namespace

TEST_F(GeometryImproverEdgeVertexTaylorPatchTester, TwoToFive)
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

  for (auto& vert : meshIn->get_vertices())
    if (vert)
    {
      // the coordinate mapping is quadratic, and so is
      // the taylor patch, so the result should be exact in this case
      auto pt = vert->get_point_orig(0);
      EXPECT_NEAR(pt.z, zmapping(pt).z, 1e-13);
    }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
