#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/incremental_mesh_boundary_snapper.hpp"
#include "stk_middle_mesh/mesh_io.hpp"
#include "stk_middle_mesh/nonconformal4.hpp"
#include "util/meshes.hpp"
#include "util/nonconformal_interface_helpers.hpp"
#include "stk_middle_mesh/application_interface.hpp"
#include "gtest/gtest.h"
#include <cmath>

#ifdef STK_BUILT_IN_SIERRA
#include "stk_middle_mesh_util/create_stk_mesh.hpp"
#include "stk_middle_mesh_util/exodus_writer.hpp"

namespace stk {
namespace middle_mesh {

using namespace nonconformal4::impl;
using namespace utils::impl;
using namespace mesh::impl;
using stk_interface::StkMeshCreator;

TEST(Interface, EigthSphereNew)
{
  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions

  if (comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 50; // 60;  // Note: this breaks on mesh 51

  for (int i = 0 /* 0 */; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i << " / " << nmeshes << std::endl;
    std::shared_ptr<mesh::Mesh> mesh1 = create_eigth_sphere(5, 5, 0.5, 1.5);

    std::shared_ptr<mesh::Mesh> mesh2 = create_eigth_sphere(5 + i, 5 + i, 0.5, 1.5);

    // print_vert_edges("mesh1_initial", mesh1);
    // print_vert_edges("mesh2_initial", mesh2);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.classifierTolerances.normalInterpolationTolerances.pointClassifierTol = 1e-8;
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::FakeParallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();
    auto middleMesh1 = interface->get_middle_grid_for_mesh1();
    auto middleMesh2 = interface->get_middle_grid_for_mesh2();
    auto mesh1Class                 = interface->get_mesh1_classification();
    auto mesh2Class                 = interface->get_mesh2_classification();    
    auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();

    test_util::test_every_element_classified(middleMesh1, mesh1Class);
    test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    test_util::test_area_per_element(mesh1, mesh1InverseClassification);
  }
}

TEST(Interface, RefiningNew)
{
  if (comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 60;

  for (int i = 0; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    mesh::impl::MeshSpec spec, spec2;
    spec.numelX = 5;
    spec.numelY = 5;
    spec.xmin   = 0;
    spec.xmax   = 1;
    spec.ymin   = 0;
    spec.ymax   = 1;

    spec2.numelX = 5 + i;
    spec2.numelY = 5 + i;
    spec2.xmin   = 0;
    spec2.xmax   = 1;
    spec2.ymin   = 0;
    spec2.ymax   = 1;

    auto func = [&](const utils::Point& pt) { return utils::Point(pt.x, 0, pt.y); };

    // std::cout << "creating mesh1" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);

    // std::cout << "\ncreating mesh2" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);

    // don't boundary snap because mesh2's domain is larger than mesh1's
    mesh::impl::ElementOperations2D elemOps;
    EXPECT_FLOAT_EQ(elemOps.compute_area(mesh1), elemOps.compute_area(mesh2));

    double eps = 1e-12;
    BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = BoundarySnapAndQualityImprovementType::None;

    MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::FakeParallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   snapOpts,
                                                   middleMeshOpts);

    interface->create_middle_grid();
    auto middleMesh1 = interface->get_middle_grid_for_mesh1();
    auto middleMesh2 = interface->get_middle_grid_for_mesh2();
    auto mesh1Class                 = interface->get_mesh1_classification();
    auto mesh2Class                 = interface->get_mesh2_classification();    
    auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();

    test_util::test_every_element_classified(middleMesh1, mesh1Class);
    test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
  }
}

TEST(Interface, AnnulusRotationNew)
{
  if (comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions

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

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::FakeParallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);
    interface->create_middle_grid();
    auto middleMesh1 = interface->get_middle_grid_for_mesh1();
    auto middleMesh2 = interface->get_middle_grid_for_mesh2();
    auto mesh1Class                 = interface->get_mesh1_classification();
    auto mesh2Class                 = interface->get_mesh2_classification();    
    auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();

    test_util::test_every_element_classified(middleMesh1, mesh1Class);
    test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
  }
}

TEST(Interface, EllipsoidNew)
{
  if (comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();  
  // project coordinates onto an ellipsoid

  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes   = 29; // 60;
  double zscale = 2;

  for (int i = 4 /* 0 */; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    mesh::impl::MeshSpec spec, spec2;
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

    // print_vert_edges("mesh1_initial", mesh1);
    // print_vert_edges("mesh2_initial", mesh2);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::FakeParallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();
    auto middleMesh1 = interface->get_middle_grid_for_mesh1();
    auto middleMesh2 = interface->get_middle_grid_for_mesh2();
    auto mesh1Class                 = interface->get_mesh1_classification();
    auto mesh2Class                 = interface->get_mesh2_classification();    
    auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();

    test_util::test_every_element_classified(middleMesh1, mesh1Class);
    test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
  }
}

#ifdef STK_BUILT_IN_SIERRA

TEST(Interface, EllipsoidFromCADNew)
{
  // project coordinates onto an ellipsoid
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  std::string meshPath            = "./";
  std::vector<std::string> fnames = {//"half_ellipsoid_2_0.5_quad.g",
                                     "half_ellipsoid_2_0.4_quad.g", "half_ellipsoid_2_0.3_quad.g",
                                     "half_ellipsoid_2_0.2_quad.g", "half_ellipsoid_2_0.1_quad.g",
                                     "half_ellipsoid_2_0.04_quad.g"};
  int nmeshes                     = fnames.size();

  for (int i = 1; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i << " / " << nmeshes << std::endl;
    stk_interface::StkMeshCreator creator1(meshPath + fnames[0]);
    std::shared_ptr<mesh::Mesh> mesh1 = creator1.create_mesh_from_part("block_1").mesh;

    stk_interface::StkMeshCreator creator2(meshPath + fnames[i]);
    std::shared_ptr<mesh::Mesh> mesh2 = creator2.create_mesh_from_part("block_1").mesh;

    std::cout << "mesh2 number of elements = " << mesh2->get_elements().size() << std::endl;

    // print_vert_edges("mesh1_initial", mesh1);
    // print_vert_edges("mesh2_initial", mesh2);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::FakeParallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();
    auto middleMesh1 = interface->get_middle_grid_for_mesh1();
    auto middleMesh2 = interface->get_middle_grid_for_mesh2();
    auto mesh1Class                 = interface->get_mesh1_classification();
    auto mesh2Class                 = interface->get_mesh2_classification();    
    auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();

    test_util::test_every_element_classified(middleMesh1, mesh1Class);
    test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    test_util::test_area_per_element(mesh1, mesh1InverseClassification);         
  }
}

#endif

TEST(Interface, AnnulusRefiningNew)
{
  if (comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();
      
  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 20; // 20;

  for (int i = 1 /*1 */; i < nmeshes; ++i)
  {
    std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;

    // std::cout << "creating mesh1" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh1 = impl::make_annulus_mesh(5, 5, 0.5, 1.5, 0);

    // std::cout << "\ncreating mesh2" << std::endl;
    std::shared_ptr<mesh::Mesh> mesh2 = impl::make_annulus_mesh(5 + i, 5 + i, 0.5, 1.5, 0);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);

    auto interface = application_interface_factory(ApplicationInterfaceType::FakeParallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();
    auto middleMesh1 = interface->get_middle_grid_for_mesh1();
    auto middleMesh2 = interface->get_middle_grid_for_mesh2();
    auto mesh1Class                 = interface->get_mesh1_classification();
    auto mesh2Class                 = interface->get_mesh2_classification();    
    auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();

    test_util::test_every_element_classified(middleMesh1, mesh1Class);
    test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
  }
}


} // namespace middle_mesh
} // namespace stk

#endif
