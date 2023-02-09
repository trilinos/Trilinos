#include "application_interface.h"
#include "util/nonconformal_interface_helpers.h"
#include "utils.h"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

class ApplicationInterfaceMPMDTester : public ::testing::Test
{
  protected:
    void runtest(int nprocs1, int nprocs2, const ParallelSearchOpts& parallelSearchOpts,
                 const VolumeSnapOpts& volumeSnapOpts, const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
                 const MiddleGridOpts& middleGridOpts)
    {
      if (nprocs1 == 0 || nprocs2 == 0)
        throw std::runtime_error("nprocs cannot be zero");

      if (nprocs1 + nprocs2 > utils::impl::comm_size(MPI_COMM_WORLD))
        throw std::runtime_error("too many processes specified");

      inputMesh  = nullptr;
      middleGrid = nullptr;

      int color = create_mesh_comm(nprocs1, nprocs2);
      create_input_meshes(color);
      bool amIMesh1 = color == 0;

      auto interface = application_interface_mpmd_factory(ApplicationInterfaceType::FakeParallel, inputMesh, amIMesh1,
                                                          MPI_COMM_WORLD, parallelSearchOpts, volumeSnapOpts,
                                                          boundarySnapOpts, middleGridOpts);

      auto middleGrid        = interface->create_middle_grid();
      auto classification    = interface->get_mesh_classification();
      auto invClassification = interface->compute_mesh_inverse_classification();

      EXPECT_EQ(classification->get_mesh(), middleGrid);
      EXPECT_EQ(invClassification->get_mesh(), inputMesh);

      test_areas_positive(middleGrid);
      test_every_element_classified(middleGrid, classification);
      test_every_element_classified_inverse(inputMesh, invClassification);
      test_area_per_element(inputMesh, invClassification);

      if (meshComm != MPI_COMM_NULL)
        MPI_Comm_free(&meshComm);
    }

    int create_mesh_comm(int nprocs1, int nprocs2)
    {
      int myrank = utils::impl::comm_rank(MPI_COMM_WORLD);
      int color  = 0;
      if (myrank >= 0 && myrank < nprocs1)
        color = 0;
      else if (myrank >= nprocs1 && myrank <= (nprocs1 + nprocs2))
        color = 1;
      else
        color = MPI_UNDEFINED;

      MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);

      return color;
    }

    void create_input_meshes(int color)
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

      auto f = [](const utils::Point& pt) { return pt; };
      if (color == 0)
        inputMesh = create_mesh(spec1, f, meshComm);
      else if (color == 1)
        inputMesh = create_mesh(spec2, f, meshComm);
    }

    MPI_Comm meshComm = MPI_COMM_NULL;
    std::shared_ptr<mesh::Mesh> inputMesh;
    std::shared_ptr<mesh::Mesh> middleGrid;
};

} // namespace

TEST_F(ApplicationInterfaceMPMDTester, Defaults)
{
  int commsize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (commsize < 2 || commsize > 4)
    GTEST_SKIP();

  ParallelSearchOpts parallelSearchOpts;
  VolumeSnapOpts volumeSnapOpts;
  BoundarySnapAndQualityImprovementOpts boundarySnapOpts;
  MiddleGridOpts middleGridOpts;

  for (int nprocs1 = 1; nprocs1 < commsize; ++nprocs1)
  {
    int nprocs2 = commsize - nprocs1;
    runtest(nprocs1, nprocs2, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);
  }
}

TEST_F(ApplicationInterfaceMPMDTester, EnableVolumeSnap)
{
  int commsize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (commsize < 2 || commsize > 4)
    GTEST_SKIP();

  ParallelSearchOpts parallelSearchOpts;
  VolumeSnapOpts volumeSnapOpts;
  volumeSnapOpts.type = VolumeSnapType::Standard;
  BoundarySnapAndQualityImprovementOpts boundarySnapOpts;
  MiddleGridOpts middleGridOpts;

  for (int nprocs1 = 1; nprocs1 < commsize; ++nprocs1)
  {
    int nprocs2 = commsize - nprocs1;
    runtest(nprocs1, nprocs2, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);
  }
}

TEST_F(ApplicationInterfaceMPMDTester, BoundarySnapThenQuality)
{
  int commsize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (commsize < 2 || commsize > 4)
    GTEST_SKIP();

  ParallelSearchOpts parallelSearchOpts;
  VolumeSnapOpts volumeSnapOpts;
  BoundarySnapAndQualityImprovementOpts boundarySnapOpts;
  boundarySnapOpts.type = BoundarySnapAndQualityImprovementType::SnapThenQuality;
  MiddleGridOpts middleGridOpts;

  for (int nprocs1 = 1; nprocs1 < commsize; ++nprocs1)
  {
    int nprocs2 = commsize - nprocs1;
    runtest(nprocs1, nprocs2, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);
  }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
