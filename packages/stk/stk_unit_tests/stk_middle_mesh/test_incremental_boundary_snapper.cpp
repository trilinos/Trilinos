#include "gtest/gtest.h"

#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/element_operations_2d.hpp"
#include "stk_middle_mesh/incremental_mesh_boundary_snapper.hpp"
#include "util/meshes.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

namespace {

void test_areas_positive(std::shared_ptr<Mesh> mesh)
{
  mesh::impl::ElementOperations2D elemOps;
  for (auto& e : mesh->get_elements())
    if (e)
    {
      EXPECT_GE(elemOps.compute_area(e), 0);
    }
}
} // namespace


TEST(IncrementalMeshBoundarySnapper, QuarterAnnulus)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  double pi   = std::atan(1) * 4;
  int nmeshes = 5;

  for (int i = 0; i < nmeshes; ++i)
  {
    MeshSpec spec, spec2;
    spec.numelX = 5;
    spec.numelY = 5;
    spec.xmin   = 0.5;
    spec.xmax   = 1.5;
    spec.ymin   = 0;
    spec.ymax   = pi / 2;

    spec2.numelX = 5 + i;
    spec2.numelY = 5 + i;
    spec2.xmin   = 0.5;
    spec2.xmax   = 1.5;
    spec2.ymin   = 0;
    spec2.ymax   = pi / 2;

    // remap coordinates to annulus
    auto func = [&](const utils::Point& pt) {
      // interpret x and y as r and theta
      double r     = pt.x;
      double theta = pt.y;

      double x = r * std::cos(theta);
      double y = r * std::sin(theta);
      double z = 0.0;
      utils::Point pt2(x, y, z);
      return pt2;
    };

    std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);
    std::shared_ptr<Mesh> mesh2 = create_mesh(spec2, func);

    // printVertEdges("mesh1_initial", mesh1);
    // printVertEdges("mesh2_initial", mesh2);

    IncrementalBoundarySnapperOpts snapperOpts;
    snapperOpts.boundarySnapNsteps = 2;
    auto snapper                   = make_incremental_boundary_snapper(mesh1, mesh2, MPI_COMM_WORLD, snapperOpts);
    snapper->snap();

    // printVertEdges("mesh1_final", mesh1);
    // printVertEdges("mesh2_final", mesh2);

    test_areas_positive(mesh1);
    test_areas_positive(mesh2);

    mesh::impl::ElementOperations2D elemOps;
    std::array<double, 2> areasLocal, areasGlobal;
    areasLocal[0] = elemOps.compute_area(mesh1);
    areasLocal[1] = elemOps.compute_area(mesh2);

    MPI_Allreduce(areasLocal.data(), areasGlobal.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_FLOAT_EQ(areasGlobal[0], areasGlobal[1]);
  }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
