#include <stk_middle_mesh/application_interface.hpp>
#include "util/nonconformal_interface_helpers.hpp"
#include <stk_middle_mesh/utils.hpp>
#include "stk_middle_mesh/communication_api_mpmd.hpp"
#include "gtest/gtest.h"
#include <utility>

namespace {

using stk::middle_mesh::utils::Point;

class XiCoordinatesForTest : public stk::middle_mesh::XiCoordinates
{
  public:
    const std::vector<stk::middle_mesh::utils::Point>& get_xi_coords(stk::middle_mesh::mesh::MeshEntityType type) override
    {
      if (type == stk::middle_mesh::mesh::MeshEntityType::Quad)
        return m_quadXiCoords;
      else if (type == stk::middle_mesh::mesh::MeshEntityType::Triangle)
        return m_triangleXiCoords;
      else
        throw std::runtime_error("only triangles and quads supported");
    }

    std::pair<double, double> get_xi_coord_range(stk::middle_mesh::mesh::MeshEntityType type) override { return std::make_pair(0, 2); }

  private:
    std::vector<stk::middle_mesh::utils::Point> m_triangleXiCoords = {Point(0.1, 0.2), Point(1.5, 0.3)};
    std::vector<stk::middle_mesh::utils::Point> m_quadXiCoords = {Point(0.1, 0.2), Point(0.6, 0.7), Point(1.4, 0.4)};
};

class ApplicationInterfaceMPMDTester : public ::testing::Test
{
  protected:
    void runtest(int nprocs1, int nprocs2,
                 const stk::middle_mesh::ParallelSearchOpts& parallelSearchOpts,
                 const stk::middle_mesh::VolumeSnapOpts& volumeSnapOpts,
                 const stk::middle_mesh::BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
                 const stk::middle_mesh::MiddleGridOpts& middleGridOpts)
    {
      if (nprocs1 == 0 || nprocs2 == 0)
        throw std::runtime_error("nprocs cannot be zero");

      if (nprocs1 + nprocs2 > stk::middle_mesh::utils::impl::comm_size(m_unionComm))
        throw std::runtime_error("too many processes specified");

      inputMesh  = nullptr;
      middleGrid = nullptr;

      int color = create_mesh_comm(nprocs1, nprocs2);
      create_input_meshes(color);
      bool amIMesh1 = color == 0;

      auto xiPts = std::make_shared<XiCoordinatesForTest>();
      auto interface = application_interface_mpmd_factory(stk::middle_mesh::ApplicationInterfaceType::FakeParallel, inputMesh, amIMesh1,
                                                          m_unionComm, xiPts, parallelSearchOpts, volumeSnapOpts,
                                                          boundarySnapOpts, middleGridOpts);

      middleGrid        = interface->create_middle_grid();
      auto classification    = interface->get_mesh_classification();
      auto invClassification = interface->compute_mesh_inverse_classification();
      auto remoteInfo        = interface->get_remote_info();
      auto xiPtsOnInputMesh  = interface->get_xi_points_on_input_mesh();

      EXPECT_EQ(classification->get_mesh(), middleGrid);
      EXPECT_EQ(invClassification->get_mesh(), inputMesh);

      stk::middle_mesh::test_util::test_areas_positive(middleGrid);
      stk::middle_mesh::test_util::test_every_element_classified(middleGrid, classification);
      stk::middle_mesh::test_util::test_every_element_classified_inverse(inputMesh, invClassification);
      stk::middle_mesh::test_util::test_area_per_element(inputMesh, invClassification);
      test_remote_info(remoteInfo, amIMesh1);
      test_remote_info(remoteInfo, !amIMesh1);
      test_xi_points(inputMesh, middleGrid, xiPtsOnInputMesh, classification, xiPts);

      if (meshComm != MPI_COMM_NULL)
        MPI_Comm_free(&meshComm);
    }

    int create_mesh_comm(int nprocs1, int nprocs2)
    {
      int myrank = stk::middle_mesh::utils::impl::comm_rank(m_unionComm);
      int color  = 0;
      if (myrank >= 0 && myrank < nprocs1)
        color = 0;
      else if (myrank >= nprocs1 && myrank <= (nprocs1 + nprocs2))
        color = 1;
      else
        color = MPI_UNDEFINED;

      MPI_Comm_split(m_unionComm, color, 0, &meshComm);

      return color;
    }

    void create_input_meshes(int color)
    {
      stk::middle_mesh::mesh::impl::MeshSpec spec1, spec2;
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

      auto f = [](const stk::middle_mesh::utils::Point& pt) { return pt; };
      if (color == 0)
        inputMesh = stk::middle_mesh::mesh::impl::create_mesh(spec1, f, meshComm);
      else if (color == 1)
        inputMesh = stk::middle_mesh::mesh::impl::create_mesh(spec2, f, meshComm);
    }

    void test_remote_info(stk::middle_mesh::mesh::FieldPtr<stk::middle_mesh::mesh::RemoteSharedEntity> remoteInfo, bool amISender)
    {
      stk::middle_mesh::mesh::FieldPtr<stk::middle_mesh::utils::Point> centroidFieldPtr;
      
      if (middleGrid)
        centroidFieldPtr = stk::middle_mesh::mesh::create_field<stk::middle_mesh::utils::Point>(middleGrid,
                                                                             stk::middle_mesh::mesh::impl::FieldShape(0, 0, 1), 1);

      if (amISender && middleGrid)
      {
        auto& centroidField = *centroidFieldPtr;
        for (auto el : middleGrid->get_elements())
          if (el)
          {
            centroidField(el, 0, 0) = stk::middle_mesh::mesh::compute_centroid(el);
          }
      }

      stk::middle_mesh::MiddleMeshFieldCommunicationMPMD<stk::middle_mesh::utils::Point> exchanger(m_unionComm, middleGrid, remoteInfo);
      exchanger.start_exchange(centroidFieldPtr, amISender);
      exchanger.finish_exchange(centroidFieldPtr, amISender);

      if (!amISender && middleGrid)
      {
        auto& centroidField = *centroidFieldPtr;
        for (auto el : middleGrid->get_elements())
          if (el)
          {
            stk::middle_mesh::utils::Point centroidLocal = stk::middle_mesh::mesh::compute_centroid(el);
            stk::middle_mesh::utils::Point centroidRemote = centroidField(el, 0, 0);
            for (int i=0; i < 3; ++i)
              EXPECT_NEAR(centroidLocal[i], centroidRemote[i], 1e-13);
          }
      }
    }

    void test_xi_points(std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh, 
                        std::shared_ptr<stk::middle_mesh::mesh::Mesh> middleGrid,
                        stk::middle_mesh::mesh::FieldPtr<stk::middle_mesh::utils::Point> xiCoordsOnInputMeshPtr,
                        stk::middle_mesh::mesh::FieldPtr<stk::middle_mesh::mesh::MeshEntityPtr> meshInToInputMeshPtr,
                        std::shared_ptr<stk::middle_mesh::XiCoordinates> xiCoords)
    {
      auto& xiCoordsOnInputMesh = *xiCoordsOnInputMeshPtr;
      auto& meshInToInputMesh   = *meshInToInputMeshPtr;
      for (auto& elMiddle : middleGrid->get_elements())
        if (elMiddle)
        {
          stk::middle_mesh::mesh::MeshEntityPtr elInput = meshInToInputMesh(elMiddle, 0, 0);
          auto& xiCoordsEl = xiCoords->get_xi_coords(elMiddle->get_type());
          std::pair<double, double> xiRangeMiddle = xiCoords->get_xi_coord_range(elMiddle->get_type());
          std::pair<double, double> xiRangeInput  = xiCoords->get_xi_coord_range(elInput->get_type());

          for (size_t i=0; i < xiCoordsEl.size(); ++i)
          {
            Point xiCoordsMiddle = stk::middle_mesh::mesh::convert_xi_coords_from_range(xiRangeMiddle.first, xiRangeMiddle.second, xiCoordsEl[i]);
            auto middleMeshCoords = stk::middle_mesh::mesh::compute_coords_from_xi_3d(elMiddle, xiCoordsMiddle);

            Point xiCoordsInput = stk::middle_mesh::mesh::convert_xi_coords_from_range(xiRangeInput.first, xiRangeInput.second,
                                                                            xiCoordsOnInputMesh(elMiddle, i, 0));
            auto inputMeshCoords  = stk::middle_mesh::mesh::compute_coords_from_xi_3d(elInput, xiCoordsInput);

            for (int d=0; d < 3; ++d)
              EXPECT_NEAR(middleMeshCoords[d], inputMeshCoords[d], 1e-12);
          }

        }      
    }    

    MPI_Comm m_unionComm = MPI_COMM_WORLD;
    MPI_Comm meshComm = MPI_COMM_NULL;
    std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh;
    std::shared_ptr<stk::middle_mesh::mesh::Mesh> middleGrid;
};

} // namespace

TEST_F(ApplicationInterfaceMPMDTester, Defaults)
{
  int commsize = stk::middle_mesh::utils::impl::comm_size(MPI_COMM_WORLD);
  if (commsize < 2 || commsize > 4)
    GTEST_SKIP();

  stk::middle_mesh::ParallelSearchOpts parallelSearchOpts;
  stk::middle_mesh::VolumeSnapOpts volumeSnapOpts;
  stk::middle_mesh::BoundarySnapAndQualityImprovementOpts boundarySnapOpts;
  stk::middle_mesh::MiddleGridOpts middleGridOpts;

  for (int nprocs1 = 1; nprocs1 < commsize; ++nprocs1)
  {
    int nprocs2 = commsize - nprocs1;
    runtest(nprocs1, nprocs2, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);
  }
}

TEST_F(ApplicationInterfaceMPMDTester, EnableVolumeSnap)
{
  int commsize = stk::middle_mesh::utils::impl::comm_size(MPI_COMM_WORLD);
  if (commsize < 2 || commsize > 4)
    GTEST_SKIP();

  stk::middle_mesh::ParallelSearchOpts parallelSearchOpts;
  stk::middle_mesh::VolumeSnapOpts volumeSnapOpts;
  volumeSnapOpts.type = stk::middle_mesh::VolumeSnapType::Standard;
  stk::middle_mesh::BoundarySnapAndQualityImprovementOpts boundarySnapOpts;
  stk::middle_mesh::MiddleGridOpts middleGridOpts;

  for (int nprocs1 = 1; nprocs1 < commsize; ++nprocs1)
  {
    int nprocs2 = commsize - nprocs1;
    runtest(nprocs1, nprocs2, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);
  }
}

TEST_F(ApplicationInterfaceMPMDTester, BoundarySnapThenQuality)
{
  int commsize = stk::middle_mesh::utils::impl::comm_size(MPI_COMM_WORLD);
  if (commsize < 2 || commsize > 4)
    GTEST_SKIP();

  stk::middle_mesh::ParallelSearchOpts parallelSearchOpts;
  stk::middle_mesh::VolumeSnapOpts volumeSnapOpts;
  stk::middle_mesh::BoundarySnapAndQualityImprovementOpts boundarySnapOpts;
  boundarySnapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::SnapThenQuality;
  stk::middle_mesh::MiddleGridOpts middleGridOpts;

  for (int nprocs1 = 1; nprocs1 < commsize; ++nprocs1)
  {
    int nprocs2 = commsize - nprocs1;
    runtest(nprocs1, nprocs2, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);
  }
}
