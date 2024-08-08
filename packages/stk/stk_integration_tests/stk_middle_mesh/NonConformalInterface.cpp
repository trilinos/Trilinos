#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/incremental_mesh_boundary_snapper.hpp"
#include "stk_middle_mesh/mesh_io.hpp"
#include "stk_middle_mesh/nonconformal4.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "util/meshes.hpp"
#include "util/nonconformal_interface_helpers.hpp"
#include "stk_middle_mesh/application_interface.hpp"
#include "gtest/gtest.h"
#include <cmath>


#ifdef STK_BUILT_FOR_SIERRA
#include "stk_middle_mesh_util/create_stk_mesh.hpp"
#include "stk_middle_mesh_util/exodus_writer.hpp"

namespace stk {
namespace middle_mesh {

using namespace nonconformal4::impl;
using namespace utils::impl;
using namespace mesh::impl;
using stk_interface::StkMeshCreator;

namespace {

class CommSplitter
{
  public:
    CommSplitter()
    {
      split_comm_world();
    }

    ~CommSplitter() { MPI_Comm_free(&m_splitComm); }

    int get_color() const { return m_color; }

    MPI_Comm get_comm() const { return m_splitComm; }

  private:

    void split_comm_world()
    {
      int commWorldSize = utils::impl::comm_size(MPI_COMM_WORLD);
      MPI_Comm newComm;
      int color = 0;
      std::cout << "commWorldSize = " << commWorldSize << std::endl;
      if (commWorldSize == 1)
      {
        MPI_Comm_dup(MPI_COMM_WORLD, &newComm);
      } else
      {
        int nprocs1 = commWorldSize / 2 + (commWorldSize % 2);
        //int nprocs2 = commWorldSize / 2;

        color = utils::impl::comm_rank(MPI_COMM_WORLD) < nprocs1  ? 0 : 1;
        std::cout << "nprocs1 = " << nprocs1 << ", color = " << color << std::endl;

        MPI_Comm_split(MPI_COMM_WORLD, color, 0, &newComm);
      }

      m_color = color;
      m_splitComm = newComm;
    }

      int m_color;
      MPI_Comm m_splitComm;

};


}

TEST(Interface, EigthSphereNewQuad)
{
  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions

  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);
  if (commSize > 4)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 50; // 60;  // Note: this breaks on mesh 51
  CommSplitter commSplitter;

  for (int i = 0 /* 0 */; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i << " / " << nmeshes << std::endl;
    }

    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    // std::cout << "creating mesh1" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)    
    {
      mesh1 = create_eigth_sphere(5, 5, 0.5, 1.5, commSplitter.get_comm());
    }

    if (commSize == 1 || commSplitter.get_color() == 1)
    {
      mesh2 = create_eigth_sphere(5 + i, 5 + i, 0.5, 1.5, commSplitter.get_comm());
    }

    // print_vert_edges("mesh1_initial", mesh1);
    // print_vert_edges("mesh2_initial", mesh2);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.classifierTolerances.normalInterpolationTolerances.pointClassifierTol = 5e-8;
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }
  }
}

TEST(Interface, EigthSphereNewTri)
{
  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions

  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);
  if (commSize > 4)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 50; // 60;  // Note: this breaks on mesh 51
  CommSplitter commSplitter;

  for (int i = 0 /* 0 */; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i << " / " << nmeshes << std::endl;
    }

    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    // std::cout << "creating mesh1" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)    
    {
      mesh1 = create_eigth_sphere(5, 5, 0.5, 1.5, commSplitter.get_comm(), true);
    }

    if (commSize == 1 || commSplitter.get_color() == 1)
    {
      mesh2 = create_eigth_sphere(5 + i, 5 + i, 0.5, 1.5, commSplitter.get_comm(), true);
    }

    // print_vert_edges("mesh1_initial", mesh1);
    // print_vert_edges("mesh2_initial", mesh2);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.classifierTolerances.normalInterpolationTolerances.pointClassifierTol = 5e-8;
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }
  }
}

TEST(Interface, RefiningNewQuad)
{
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);
  if (comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 60;

  CommSplitter commSplitter;

  for (int i = 0; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    }
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
    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    // std::cout << "creating mesh1" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh1 = create_mesh(spec, func, commSplitter.get_comm());

    // std::cout << "\ncreating mesh2" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 1)    
      mesh2 = create_mesh(spec2, func, commSplitter.get_comm());

    mesh::impl::ElementOperations2D elemOps;
    double mesh1Area = stk::get_global_sum(MPI_COMM_WORLD, mesh1 ? elemOps.compute_area(mesh1) : 0.0);
    double mesh2Area = stk::get_global_sum(MPI_COMM_WORLD, mesh2 ? elemOps.compute_area(mesh2) : 0.0);
    EXPECT_FLOAT_EQ(mesh1Area, mesh2Area);

    double eps = 1e-12;
    BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = BoundarySnapAndQualityImprovementType::None;

    MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   snapOpts,
                                                   middleMeshOpts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }
  }
}

TEST(Interface, RefiningNewTri)
{
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);
  if (comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 60;

  CommSplitter commSplitter;

  for (int i = 0; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    }
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
    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    // std::cout << "creating mesh1" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh1 = create_mesh(spec, func, commSplitter.get_comm(), true);

    // std::cout << "\ncreating mesh2" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 1)    
      mesh2 = create_mesh(spec2, func, commSplitter.get_comm(), true);

    mesh::impl::ElementOperations2D elemOps;
    double mesh1Area = stk::get_global_sum(MPI_COMM_WORLD, mesh1 ? elemOps.compute_area(mesh1) : 0.0);
    double mesh2Area = stk::get_global_sum(MPI_COMM_WORLD, mesh2 ? elemOps.compute_area(mesh2) : 0.0);
    EXPECT_FLOAT_EQ(mesh1Area, mesh2Area);

    double eps = 1e-12;
    BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = BoundarySnapAndQualityImprovementType::None;

    MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   snapOpts,
                                                   middleMeshOpts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }
  }
}

TEST(Interface, AnnulusRotationNewQuad)
{
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);  
  if (commSize > 4)
    GTEST_SKIP();

  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions
  CommSplitter commSplitter;

  std::cout << std::setprecision(16) << std::endl;
  double pi     = std::atan(1) * 4;
  int nmeshes   = 640;
  double dtheta = pi / (16 * nmeshes);

  for (int i = 1; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
      std::cout << "dtheta = " << i * dtheta * 180.0 / pi << std::endl;
    }

    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    // std::cout << "creating mesh1" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh1 = impl::make_annulus_mesh(10, 10, 0.5, 1.5, 0, commSplitter.get_comm());

    // std::cout << "\ncreating mesh2" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh2 = impl::make_annulus_mesh(13, 13, 0.5, 1.5, i * dtheta, commSplitter.get_comm());

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);
    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      EXPECT_GE(middleMesh1->get_elements().size(), mesh1->get_elements().size());
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }
  }   
}

TEST(Interface, AnnulusRotationNewTri)
{
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);  
  if (commSize > 4)
    GTEST_SKIP();

  // project coordinates onto sector of a sphere.  Change the size of
  // the mesh2 sector to test different kinds of topology interactions
  CommSplitter commSplitter;

  std::cout << std::setprecision(16) << std::endl;
  double pi     = std::atan(1) * 4;
  int nmeshes   = 640;
  double dtheta = pi / (16 * nmeshes);

  for (int i = 1; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
      std::cout << "dtheta = " << i * dtheta * 180.0 / pi << std::endl;
    }

    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    // std::cout << "creating mesh1" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh1 = impl::make_annulus_mesh(10, 10, 0.5, 1.5, 0, commSplitter.get_comm(), true);

    // std::cout << "\ncreating mesh2" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh2 = impl::make_annulus_mesh(13, 13, 0.5, 1.5, i * dtheta, commSplitter.get_comm(), true);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);
    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      EXPECT_GE(middleMesh1->get_elements().size(), mesh1->get_elements().size());
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }
  }   
}

TEST(Interface, EllipsoidNewQuad)
{
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);   
  if (commSize > 4)
    GTEST_SKIP();  
  // project coordinates onto an ellipsoid

  CommSplitter commSplitter;

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes   = 29; // 60;
  double zscale = 2;

  for (int i = 1 /* 0 */; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    }
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

    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    // std::cout << "creating mesh1" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh1 = create_mesh(spec, func, commSplitter.get_comm());

    // std::cout << "\ncreating mesh2" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh2 = create_mesh(spec2, func, commSplitter.get_comm());

    // print_vert_edges("mesh1_initial", mesh1);
    // print_vert_edges("mesh2_initial", mesh2);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      EXPECT_GE(middleMesh1->get_elements().size(), mesh1->get_elements().size());
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    } 
  }
}


TEST(Interface, EllipsoidNewTri)
{
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);   
  if (commSize > 4)
    GTEST_SKIP();  
  // project coordinates onto an ellipsoid

  CommSplitter commSplitter;

  std::cout << std::setprecision(16) << std::endl;
  int nmeshes   = 29; // 60;
  double zscale = 2;

  for (int i = 1 /* 0 */; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    }
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

    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    // std::cout << "creating mesh1" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh1 = create_mesh(spec, func, commSplitter.get_comm(), true);

    // std::cout << "\ncreating mesh2" << std::endl;
    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh2 = create_mesh(spec2, func, commSplitter.get_comm(), true);

    // print_vert_edges("mesh1_initial", mesh1);
    // print_vert_edges("mesh2_initial", mesh2);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      EXPECT_GE(middleMesh1->get_elements().size(), mesh1->get_elements().size());
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    } 
  }
}

#ifdef STK_BUILT_FOR_SIERRA

TEST(Interface, EllipsoidFromCADNewQuad)
{
  // project coordinates onto an ellipsoid
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);     
  if (commSize > 4)
    GTEST_SKIP();

  CommSplitter commSplitter;


  std::cout << std::setprecision(16) << std::endl;
  std::string meshPath            = "./";
  std::vector<std::string> fnames = {//"half_ellipsoid_2_0.5_quad.g",
                                     "half_ellipsoid_2_0.4_quad.g", "half_ellipsoid_2_0.3_quad.g",
                                     "half_ellipsoid_2_0.2_quad.g", "half_ellipsoid_2_0.1_quad.g",
                                     "half_ellipsoid_2_0.04_quad.g"};
  int nmeshes                     = fnames.size();

  for (int i = 1 /* 1 */; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i << " / " << nmeshes << std::endl;
    }
    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    if (commSize == 1 || commSplitter.get_color() == 0)
    {
      stk_interface::StkMeshCreator creator1(meshPath + fnames[0], "RIB", commSplitter.get_comm());
      mesh1 = creator1.create_mesh_from_part("block_1").mesh;
    }

    if (commSize == 1 || commSplitter.get_color() == 1)
    {
      stk_interface::StkMeshCreator creator2(meshPath + fnames[i], "RIB", commSplitter.get_comm());
      mesh2 = creator2.create_mesh_from_part("block_1").mesh;
    }


    // print_vert_edges("mesh1_initial", mesh1);
    // print_vert_edges("mesh2_initial", mesh2);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      EXPECT_GE(middleMesh1->get_elements().size(), mesh1->get_elements().size());
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }         
  }
}

TEST(Interface, EllipsoidFromCADNewTri)
{
  // project coordinates onto an ellipsoid
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);     
  if (commSize > 4)
    GTEST_SKIP();

  CommSplitter commSplitter;


  std::cout << std::setprecision(16) << std::endl;
  std::string meshPath            = "./";
  std::vector<std::string> fnames = {//"half_ellipsoid_2_0.5_tri.g",
                                     "half_ellipsoid_2_0.4_tri.g", "half_ellipsoid_2_0.3_tri.g",
                                     "half_ellipsoid_2_0.2_tri.g", "half_ellipsoid_2_0.1_tri.g",
                                     "half_ellipsoid_2_0.04_tri.g"};
  int nmeshes                     = fnames.size();

  for (int i = 1 /* 1 */; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i << " / " << nmeshes << std::endl;
    }
    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    if (commSize == 1 || commSplitter.get_color() == 0)
    {
      stk_interface::StkMeshCreator creator1(meshPath + fnames[0], "RIB", commSplitter.get_comm());
      mesh1 = creator1.create_mesh_from_part("block_1").mesh;
    }

    if (commSize == 1 || commSplitter.get_color() == 1)
    {
      stk_interface::StkMeshCreator creator2(meshPath + fnames[i], "RIB", commSplitter.get_comm());
      mesh2 = creator2.create_mesh_from_part("block_1").mesh;
    }


    // print_vert_edges("mesh1_initial", mesh1);
    // print_vert_edges("mesh2_initial", mesh2);

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      EXPECT_GE(middleMesh1->get_elements().size(), mesh1->get_elements().size());
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }         
  }
}

#endif

TEST(Interface, AnnulusRefiningNewQuad)
{ 
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);      
  if (commSize > 4)
    GTEST_SKIP();
      
  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 20; // 20;
  CommSplitter commSplitter;

  for (int i = 1 /*1 */; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    }

    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    if (commSize == 1 || commSplitter.get_color() == 0)
    {
      mesh1 = impl::make_annulus_mesh(5, 5, 0.5, 1.5, 0, commSplitter.get_comm());
    }

    if (commSize == 1 || commSplitter.get_color() == 1)
    {
      mesh2 = impl::make_annulus_mesh(5 + i, 5 + i, 0.5, 1.5, 0, commSplitter.get_comm());
    }

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);

    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      EXPECT_GE(middleMesh1->get_elements().size(), mesh1->get_elements().size());
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }   
  }
}

TEST(Interface, AnnulusRefiningNewTri)
{ 
  int commSize = comm_size(MPI_COMM_WORLD);
  int commRank = comm_rank(MPI_COMM_WORLD);      
  if (commSize > 4)
    GTEST_SKIP();
      
  std::cout << std::setprecision(16) << std::endl;
  int nmeshes = 20; // 20;
  CommSplitter commSplitter;

  for (int i = 1 /*1 */; i < nmeshes; ++i)
  {
    if (commRank == 0)
    {
      std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    }

    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    if (commSize == 1 || commSplitter.get_color() == 0)
    {
      mesh1 = impl::make_annulus_mesh(5, 5, 0.5, 1.5, 0, commSplitter.get_comm(), true);
    }

    if (commSize == 1 || commSplitter.get_color() == 1)
    {
      mesh2 = impl::make_annulus_mesh(5 + i, 5 + i, 0.5, 1.5, 0, commSplitter.get_comm(), true);
    }

    double eps = 1e-12;
    MiddleGridOpts opts;
    opts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    opts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);

    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                   MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                   VolumeSnapOpts(),
                                                   BoundarySnapAndQualityImprovementOpts(),
                                                   opts);

    interface->create_middle_grid();

    if (mesh1)
    {
      auto middleMesh1 = interface->get_middle_grid_for_mesh1();
      auto mesh1Class                 = interface->get_mesh1_classification();
      auto mesh1InverseClassification = interface->compute_mesh1_inverse_classification();
      EXPECT_GE(middleMesh1->get_elements().size(), mesh1->get_elements().size());
      test_util::test_every_element_classified(middleMesh1, mesh1Class);
      test_util::test_area_per_element(mesh1, mesh1InverseClassification);    
    }

    if (mesh2)
    {
      auto middleMesh2 = interface->get_middle_grid_for_mesh2();
      auto mesh2Class                 = interface->get_mesh2_classification();    
      test_util::test_every_element_classified(middleMesh2, mesh2Class);    
    }   
  }
}


} // namespace middle_mesh
} // namespace stk

#endif
