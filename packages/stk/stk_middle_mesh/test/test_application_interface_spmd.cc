#include "application_interface.h"
#include "util/nonconformal_interface_helpers.h"
#include "utils.h"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

using namespace stk::middle_mesh;

class ApplicationInterfaceSPMDTester : public ::testing::Test
{
  protected:
    void runtest(int nprocs, const ParallelSearchOpts& parallelSearchOpts, const VolumeSnapOpts& volumeSnapOpts,
                 const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts, const MiddleGridOpts& middleGridOpts)
    {
      if (nprocs == 0)
        throw std::runtime_error("nprocs cannot be zero");

      if (nprocs > utils::impl::comm_size(MPI_COMM_WORLD))
        throw std::runtime_error("too many processes specified");

      inputMesh1 = nullptr;
      inputMesh2 = nullptr;
      middleGrid = nullptr;

      bool onMeshComm = create_mesh_comm(nprocs);

      if (onMeshComm)
      {
        create_input_meshes();

        auto interface =
            application_interface_spmd_factory(ApplicationInterfaceType::FakeParallel, inputMesh1, inputMesh2,
                                               parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);

        interface->create_middle_grid();
        auto middleGrid1        = interface->get_middle_grid_for_mesh1();
        auto middleGrid2        = interface->get_middle_grid_for_mesh2();
        auto classification1    = interface->get_mesh1_classification();
        auto classification2    = interface->get_mesh2_classification();
        auto invClassification1 = interface->compute_mesh1_inverse_classification();
        auto invClassification2 = interface->compute_mesh2_inverse_classification();

        EXPECT_EQ(classification1->get_mesh(), middleGrid1);
        EXPECT_EQ(classification2->get_mesh(), middleGrid2);
        EXPECT_EQ(invClassification1->get_mesh(), inputMesh1);
        EXPECT_EQ(invClassification2->get_mesh(), inputMesh2);

        test_areas_positive(middleGrid1);
        test_areas_positive(middleGrid2);

        test_every_element_classified(middleGrid1, classification1);
        test_every_element_classified(middleGrid2, classification2);

        test_every_element_classified_inverse(inputMesh1, invClassification1);
        test_every_element_classified_inverse(inputMesh2, invClassification2);

        test_area_per_element(inputMesh1, invClassification1);
        test_area_per_element(inputMesh2, invClassification2);

        MPI_Comm_free(&meshComm);
      }
    }

    bool create_mesh_comm(int nprocs)
    {
      int myrank = utils::impl::comm_rank(MPI_COMM_WORLD);
      int color  = myrank < nprocs ? 0 : MPI_UNDEFINED;
      MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);

      return color == 0;
    }

    void create_input_meshes()
    {
      mesh::impl::MeshSpec spec1, spec2;
      spec1.xmin   = 0;
      spec2.xmin   = 0;
      spec1.xmax   = 1;
      spec2.xmax   = 1;
      spec1.ymin   = 0;
      spec2.ymin   = 0;
      spec1.ymax   = 1;
      spec2.ymax   = 1;
      spec1.numelX = 4;
      spec2.numelX = 5;
      spec1.numelY = 4;
      spec2.numelY = 5;

      auto f     = [](const utils::Point& pt) { return pt; };
      inputMesh1 = create_mesh(spec1, f, meshComm);
      inputMesh2 = create_mesh(spec2, f, meshComm);
    }

    MPI_Comm meshComm = MPI_COMM_NULL;
    std::shared_ptr<mesh::Mesh> inputMesh1;
    std::shared_ptr<mesh::Mesh> inputMesh2;
    std::shared_ptr<mesh::Mesh> middleGrid;
};

} // namespace

TEST_F(ApplicationInterfaceSPMDTester, Defaults)
{
  int commsize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (commsize > 4)
    GTEST_SKIP();

  ParallelSearchOpts parallelSearchOpts;
  VolumeSnapOpts volumeSnapOpts;
  BoundarySnapAndQualityImprovementOpts boundarySnapOpts;
  MiddleGridOpts middleGridOpts;

  for (int nprocs = 1; nprocs <= commsize; ++nprocs)
  {
    runtest(nprocs, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);
  }
}

TEST_F(ApplicationInterfaceSPMDTester, EnableVolumeSnap)
{
  int commsize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (commsize > 4)
    GTEST_SKIP();

  ParallelSearchOpts parallelSearchOpts;
  VolumeSnapOpts volumeSnapOpts;
  volumeSnapOpts.type = VolumeSnapType::Standard;
  BoundarySnapAndQualityImprovementOpts boundarySnapOpts;
  MiddleGridOpts middleGridOpts;

  for (int nprocs = 1; nprocs <= commsize; ++nprocs)
  {
    runtest(nprocs, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);
  }
}

TEST_F(ApplicationInterfaceSPMDTester, BoundarySnapThenQuality)
{
  int commsize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (commsize > 4)
    GTEST_SKIP();

  ParallelSearchOpts parallelSearchOpts;
  VolumeSnapOpts volumeSnapOpts;
  BoundarySnapAndQualityImprovementOpts boundarySnapOpts;
  boundarySnapOpts.type = BoundarySnapAndQualityImprovementType::SnapThenQuality;
  MiddleGridOpts middleGridOpts;

  for (int nprocs = 1; nprocs < commsize; ++nprocs)
  {
    runtest(nprocs, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);
  }
}

TEST(ApplicationInterfaceSPMD, DifferentCommunicatorsError)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec1;
  spec1.xmin   = 0;
  spec1.xmax   = 1;
  spec1.ymin   = 0;
  spec1.ymax   = 1;
  spec1.numelX = 4;
  spec1.numelY = 4;

  MPI_Comm comm1 = MPI_COMM_WORLD;
  MPI_Comm comm2 = utils::impl::comm_dup(comm1);

  auto f          = [](const utils::Point& pt) { return pt; };
  auto inputMesh1 = mesh::impl::create_mesh(spec1, f, comm1);
  auto inputMesh2 = mesh::impl::create_mesh(spec1, f, comm2);

  EXPECT_ANY_THROW(application_interface_spmd_factory(ApplicationInterfaceType::FakeParallel, inputMesh1, inputMesh2));

  MPI_Comm_free(&comm2);
}

TEST(ApplicationInterfaceSPMD, MeshNotCreatedError)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec1;
  spec1.xmin   = 0;
  spec1.xmax   = 1;
  spec1.ymin   = 0;
  spec1.ymax   = 1;
  spec1.numelX = 4;
  spec1.numelY = 4;

  MPI_Comm comm1 = MPI_COMM_WORLD;

  auto f          = [](const utils::Point& pt) { return pt; };
  auto inputMesh1 = mesh::impl::create_mesh(spec1, f, comm1);
  auto inputMesh2 = mesh::impl::create_mesh(spec1, f, comm1);

  auto interface = application_interface_spmd_factory(ApplicationInterfaceType::FakeParallel, inputMesh1, inputMesh2);

  EXPECT_ANY_THROW(interface->get_middle_grid_for_mesh1());
  EXPECT_ANY_THROW(interface->get_middle_grid_for_mesh2());
  EXPECT_ANY_THROW(interface->get_mesh1_classification());
  EXPECT_ANY_THROW(interface->get_mesh2_classification());
  EXPECT_ANY_THROW(interface->compute_mesh1_inverse_classification());
  EXPECT_ANY_THROW(interface->compute_mesh2_inverse_classification());
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
