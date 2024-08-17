#include "gtest/gtest.h"

#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/boundary_fixture.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/create_mesh_quality_improver.hpp"
#include "stk_middle_mesh/incremental_mesh_boundary_snapper.hpp"
#include "stk_middle_mesh/mesh_boundary_snapper.hpp"
#include "stk_middle_mesh/mesh_quality_improver.hpp"
#include "stk_middle_mesh/mesh_quality_statistics.hpp"
#include "stk_middle_mesh/regularized_distortion_metric.hpp"
#include "util/meshes.hpp"

#ifdef STK_BUILT_FOR_SIERRA
#include "stk_middle_mesh_util/create_stk_mesh.hpp"
#endif

#include "stk_middle_mesh/mesh_io.hpp" //TODO: DEBUGGING

namespace stk {
namespace middle_mesh {

using namespace mesh::impl;
using namespace utils::impl;
using stk_interface::StkMeshCreator;

namespace {

// func is a lambda (or any other callable) that maps a given
// point to where it is supposted to be on the surface
// Ideally this would be the projection to the closest
// point on the surface, but sometimes we have to settle
// for the projection along a coordinate axis
template <typename T>
double compute_surface_error(std::shared_ptr<mesh::Mesh> mesh, T func)
{
  double maxError = 0, minError = std::numeric_limits<double>::max(), avgError = 0, errorNorm = 0;

  int nverts = 0;
  for (auto& vert : mesh->get_vertices())
  {
    if (vert)
    {
      auto pt          = vert->get_point_orig(0);
      auto ptOnSurface = func(pt);

      auto error = std::sqrt(dot(pt - ptOnSurface, pt - ptOnSurface));
      maxError   = std::max(error, maxError);
      minError   = std::min(error, minError);
      avgError += error;
      errorNorm += error * error;
      nverts++;
    }
  }

  avgError  = avgError / nverts;
  errorNorm = std::sqrt(errorNorm);

  std::cout << "max_error = " << maxError << ", min_error = " << minError << ", avg error = " << avgError
            << ", error_norm = " << errorNorm << std::endl;

  return maxError;
}
} // namespace

TEST(MeshQualityImprover, AnnulusRotation)
{
  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions
  if (comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  double pi     = std::atan(1) * 4;
  int nmeshes   = 640;
  double dtheta = pi / (16 * nmeshes);

  for (int i = 1; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    std::cout << "dtheta = " << i * dtheta * 180.0 / pi << std::endl;

    // std::cout << "creating mesh1" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh1 = impl::make_annulus_mesh(10, 10, 0.5, 1.5, 0);

    // std::cout << "\ncreating mesh2" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh2 = impl::make_annulus_mesh(13, 13, 0.5, 1.5, i * dtheta);

    auto bsnapper = make_incremental_boundary_snapper(mesh1, mesh2, MPI_COMM_WORLD);
    bsnapper->snap();

    std::cout << "number of invalid points = " << bsnapper->get_mesh2_quality_improver()->count_invalid_points()
              << std::endl;
    EXPECT_EQ(bsnapper->get_mesh2_quality_improver()->count_invalid_points(), 0);
  }
}

TEST(MeshQualityImprover, Ellipsoid)
{
  // project coordinates onto an ellipsoid

  if (comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  int nmeshes   = 30; // 60;
  double zscale = 2;

  for (int i = 0 /* 0 */; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    MeshSpec spec, spec2;
    spec.numelX = 5;
    spec.numelY = 5;
    spec.xmin   = -0.75; //-1;
    spec.xmax   = 0.75;  // 1;
    spec.ymin   = -0.75; //-1;
    spec.ymax   = 0.75;  // 1;

    spec2.numelX = 5 + i;
    spec2.numelY = 5 + i;
    spec2.xmin   = -0.75; //-1;
    spec2.xmax   = 0.75;  // 1;
    spec2.ymin   = -0.75; //-1;
    spec2.ymax   = 0.75;  // 1;

    auto func = [&](const utils::Point& pt) {
      double x = pt.x;
      double y = pt.y;

      double xprime = x * std::sqrt(std::max(1 - y * y / 2, 0.0));
      double yprime = y * std::sqrt(std::max(1 - x * x / 2, 0.0));
      double zprime = zscale * std::sqrt(std::max(1 - x * x - y * y, 0.0));
      utils::Point pt2(xprime, yprime, zprime);
      return pt2;
    };

    // std::cout << "creating mesh1" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);

    // std::cout << "\ncreating mesh2" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);

    // printVertEdges("mesh1_initial", mesh1);
    // printVertEdges("mesh2_initial", mesh2);

    auto bsnapper = make_incremental_boundary_snapper(mesh1, mesh2, MPI_COMM_WORLD);
    bsnapper->snap();
    std::cout << "after boundary snapper, number of invalid verts = "
              << bsnapper->get_mesh2_quality_improver()->count_invalid_points() << std::endl;
    EXPECT_EQ(bsnapper->get_mesh2_quality_improver()->count_invalid_points(), 0);

    // printVertEdges("mesh1_snapped", mesh1);
    // printVertEdges("mesh2_snapped", mesh2);
  }
}

#ifdef STK_BUILT_FOR_SIERRA

TEST(MeshQualityImprover, EllipsoidFromCAD)
{
  // project coordinates onto an ellipsoid

  if (comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  std::string meshPath            = "./";
  std::vector<std::string> fnames = {"half_ellipsoid_2_0.5_quad.g", "half_ellipsoid_2_0.4_quad.g",
                                     "half_ellipsoid_2_0.3_quad.g", "half_ellipsoid_2_0.2_quad.g",
                                     "half_ellipsoid_2_0.1_quad.g"};
  int nmeshes                     = fnames.size();

  for (int i = 0; i < nmeshes; ++i)
  {
    StkMeshCreator creator1(meshPath + fnames[0]);
    std::shared_ptr<mesh::Mesh> mesh1 = creator1.create_mesh_from_part("block_1").mesh;

    StkMeshCreator creator2(meshPath + fnames[i]);
    std::shared_ptr<mesh::Mesh> mesh2 = creator2.create_mesh_from_part("block_1").mesh;

    std::cout << "mesh2 number of elements = " << mesh2->get_elements().size() << std::endl;

    // printVertEdges("mesh1_initial", mesh1);
    // printVertEdges("mesh2_initial", mesh2);

    auto bsnapper = make_incremental_boundary_snapper(mesh1, mesh2, MPI_COMM_WORLD);
    bsnapper->snap();
    std::cout << "after boundary snapper, number of invalid verts = "
              << bsnapper->get_mesh2_quality_improver()->count_invalid_points() << std::endl;
    EXPECT_EQ(bsnapper->get_mesh2_quality_improver()->count_invalid_points(), 0);

    // printVertEdges("mesh1_snapped", mesh1);
    // printVertEdges("mesh2_snapped", mesh2);
  }
}

#endif

TEST(MeshQualityImprover, AnnulusRefining)
{
  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions
  if (comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 17; // 17;

  for (int i = 1 /*1*/; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;

    auto errorFunc = [](const utils::Point& pt) { return utils::Point(pt.x, pt.y, 0); };

    std::shared_ptr<mesh::Mesh> mesh1 = impl::make_annulus_mesh(5, 5, 0.5, 1.5, 0);

    std::shared_ptr<mesh::Mesh> mesh2 = impl::make_annulus_mesh(5 + i, 5 + i, 0.5, 1.5, 0);

    std::cout << "initial mesh1 surface error: ";
    compute_surface_error(mesh1, errorFunc);

    std::cout << "initial mesh2 surface error: ";
    compute_surface_error(mesh2, errorFunc);

    auto bsnapper = make_incremental_boundary_snapper(mesh1, mesh2, MPI_COMM_WORLD);
    bsnapper->snap();

    std::cout << "after snap mesh1 surface error: ";
    compute_surface_error(mesh1, errorFunc);

    std::cout << "after snap mesh2 surface error: ";
    compute_surface_error(mesh2, errorFunc);

    std::cout << "number of invalid points = " << bsnapper->get_mesh2_quality_improver()->count_invalid_points()
              << std::endl;


    std::cout << "final mesh1 surface error: ";
    double maxError1 = compute_surface_error(mesh1, errorFunc);

    std::cout << "final mesh2 surface error: ";
    double maxError2 = compute_surface_error(mesh2, errorFunc);

    EXPECT_EQ(bsnapper->get_mesh2_quality_improver()->count_invalid_points(), 0);
    EXPECT_NEAR(maxError1, 0, 1e-13);
    EXPECT_NEAR(maxError2, 0, 1e-13);
  }
}

TEST(MeshQualityImprover, InclinePlaneRefining)
{
  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions
  if (comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes  = 17; // 17;
  double pi    = std::atan(1) * 4;
  double theta = 30 * pi / 180;

  for (int i = 1 /*1*/; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;

    MeshSpec spec, spec2;
    spec.numelX = 5;
    spec.numelY = 5;
    spec.xmin   = 0;
    spec.xmax   = 2;
    spec.ymin   = 0;
    spec.ymax   = 3;

    spec2.numelX = 5 + i;
    spec2.numelY = 5 + i;
    spec2.xmin   = 0;
    spec2.xmax   = 2;
    spec2.ymin   = 0;
    spec2.ymax   = 3;

    auto func      = [&](const utils::Point& pt) { return utils::Point(pt.x, pt.y, pt.x * std::tan(theta)); };
    auto errorFunc = [&](const utils::Point& pt) {
      utils::Point normal(-std::cos(pi / 2 - theta), 0, std::sin(pi / 2 - theta));
      return pt - dot(pt, normal) * normal;
    };

    // std::cout << "creating mesh1" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);

    // std::cout << "\ncreating mesh2" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);

    // printVertEdges("mesh1_initial", mesh1);
    // printVertEdges("mesh2_initial", mesh2);

    auto bsnapper = make_incremental_boundary_snapper(mesh1, mesh2, MPI_COMM_WORLD);
    bsnapper->snap();

    // printVertEdges("mesh1_snapped", mesh1);
    // printVertEdges("mesh2_snapped", mesh2);
    std::cout << "number of invalid points = " << bsnapper->get_mesh2_quality_improver()->count_invalid_points()
              << std::endl;
    // EXPECT_GE(fixer2->count_invalid_points(), 0);
    // printVertEdges("mesh1_fixed", mesh1);
    // printVertEdges("mesh2_fixed", mesh2);

    std::cout << "final mesh1 surface error: ";
    double maxError1 = compute_surface_error(mesh1, errorFunc);

    std::cout << "final mesh2 surface error: ";
    double maxError2 = compute_surface_error(mesh2, errorFunc);

    EXPECT_EQ(bsnapper->get_mesh2_quality_improver()->count_invalid_points(), 0);
    EXPECT_NEAR(maxError1, 0, 1e-11);
    EXPECT_NEAR(maxError2, 0, 1e-11);
  }
}

TEST(MeshQualityImprover, TorusRotation)
{
  // project annulus onto surface of sphere (which is not really of torus)
  // Rotate the annulus of
  // the mesh2 sector to test different kinds of topology interactions
  if (comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  double pi     = std::atan(1) * 4;
  int nmeshes   = 640;
  double dtheta = pi / (16 * nmeshes);

  double radiusIn  = 0.5;
  double radiusOut = 1.5;

  auto func = [radiusOut](const double x, const double y) {
    // project upward onto surface of sphere
    return std::sqrt(std::max(radiusOut * radiusOut - x * x - y * y, 0.0));
  };

  for (int i = 1; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    std::cout << "dtheta = " << i * dtheta * 180.0 / pi << std::endl;

    // std::cout << "creating mesh1" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh1 = impl::make_annulus_mesh(10, 10, radiusIn, radiusOut, 0, MPI_COMM_WORLD, false, func);

    // std::cout << "\ncreating mesh2" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh2 = impl::make_annulus_mesh(13, 13, radiusIn, radiusOut, i * dtheta, MPI_COMM_WORLD, false, func);

    // printVertEdges("mesh_initial", mesh2);

    auto bsnapper = make_incremental_boundary_snapper(mesh1, mesh2, MPI_COMM_WORLD);
    bsnapper->snap();

    // printVertEdges("mesh_snapped", mesh2);
    std::cout << "number of invalid points = " << bsnapper->get_mesh2_quality_improver()->count_invalid_points()
              << std::endl;
    EXPECT_EQ(bsnapper->get_mesh2_quality_improver()->count_invalid_points(), 0);
  }
}

TEST(MeshQualityImprover, EigthSphereRefining)
{
  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions
  if (comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes   = 60; /*60;*/
  double radius = 1.5;
  for (int i = 0; i < nmeshes; ++i)
  {
    // computes closest point on sphere to given point
    auto errorFunc = [&](const utils::Point& pt) {
      auto cp = pt / std::sqrt(dot(pt, pt));
      return radius * cp;
    };

    std::shared_ptr<mesh::Mesh> mesh1 = create_eigth_sphere(5, 5, 0.5, radius);

    std::shared_ptr<mesh::Mesh> mesh2 = create_eigth_sphere(5 + i, 5 + i, 0.5, radius);

    std::cout << "initial mesh1 surface error: ";
    compute_surface_error(mesh1, errorFunc);

    std::cout << "initial mesh2 surface error: ";
    compute_surface_error(mesh2, errorFunc);

    // printVertEdges("mesh1_initial", mesh1);
    // printVertEdges("mesh2_initial", mesh2);

    auto bsnapper = make_incremental_boundary_snapper(mesh1, mesh2, MPI_COMM_WORLD);
    bsnapper->snap();

    std::cout << "after snap mesh1 surface error: ";
    compute_surface_error(mesh1, errorFunc);

    std::cout << "after mesh2 surface error: ";
    compute_surface_error(mesh2, errorFunc);

    // printVertEdges("mesh1_snapped", mesh1);
    // printVertEdges("mesh2_snapped", mesh2);

    EXPECT_EQ(bsnapper->get_mesh2_quality_improver()->count_invalid_points(), 0);

    // printVertEdges("mesh1_fixed", mesh1);
    // printVertEdges("mesh2_fixed", mesh2);
  }
}

} // namespace middle_mesh
} // namespace stk
