#include <stk_middle_mesh/application_interface.hpp>
#include "util/nonconformal_interface_helpers.hpp"
#include <stk_middle_mesh/utils.hpp>
#include "stk_middle_mesh/communication_api_spmd.hpp"
#include "gtest/gtest.h"

namespace {

using namespace stk::middle_mesh;


class XiCoordinatesForTest : public XiCoordinates
{
  public:
    const std::vector<utils::Point>& get_xi_coords(mesh::MeshEntityType type) override
    {
      if (type == mesh::MeshEntityType::Quad)
        return m_quadXiCoords;
      else if (type == mesh::MeshEntityType::Triangle)
        return m_triangleXiCoords;
      else
        throw std::runtime_error("only triangles and quads supported");
    }

    std::pair<double, double> get_xi_coord_range(mesh::MeshEntityType type) override { return std::make_pair(0, 2); }

  private:
    std::vector<utils::Point> m_triangleXiCoords = {utils::Point(0.1, 0.2), utils::Point(1.5, 0.3)};
    std::vector<utils::Point> m_quadXiCoords = {utils::Point(0.1, 0.2), utils::Point(0.6, 0.7), utils::Point(1.4, 0.4)};
};

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
      middleGrid1 = nullptr;
      middleGrid2 = nullptr;

      bool onMeshComm = create_mesh_comm(nprocs);

      if (onMeshComm)
      {
        create_input_meshes();

        auto xiPts = std::make_shared<XiCoordinatesForTest>();
        auto interface =
            application_interface_spmd_factory(ApplicationInterfaceType::FakeParallel, inputMesh1, inputMesh2, xiPts,
                                               parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts);

        interface->create_middle_grid();
        middleGrid1        = interface->get_middle_grid_for_mesh1();
        middleGrid2        = interface->get_middle_grid_for_mesh2();
        auto classification1    = interface->get_mesh1_classification();
        auto classification2    = interface->get_mesh2_classification();
        auto invClassification1 = interface->compute_mesh1_inverse_classification();
        auto invClassification2 = interface->compute_mesh2_inverse_classification();
        auto remoteInfoOneToTwo = interface->get_remote_info_mesh_one_to_two();
        auto remoteInfoTwoToOne = interface->get_remote_info_mesh_two_to_one();
        auto xiPtsOnMesh1       = interface->get_xi_points_on_mesh1();
        auto xiPtsOnMesh2      = interface->get_xi_points_on_mesh2();


        EXPECT_EQ(classification1->get_mesh(), middleGrid1);
        EXPECT_EQ(classification2->get_mesh(), middleGrid2);
        EXPECT_EQ(invClassification1->get_mesh(), inputMesh1);
        EXPECT_EQ(invClassification2->get_mesh(), inputMesh2);

        test_util::test_areas_positive(middleGrid1);
        test_util::test_areas_positive(middleGrid2);

        test_util::test_every_element_classified(middleGrid1, classification1);
        test_util::test_every_element_classified(middleGrid2, classification2);

        test_util::test_every_element_classified_inverse(inputMesh1, invClassification1);
        test_util::test_every_element_classified_inverse(inputMesh2, invClassification2);

        test_util::test_area_per_element(inputMesh1, invClassification1);
        test_util::test_area_per_element(inputMesh2, invClassification2);

        test_remote_info(remoteInfoOneToTwo, remoteInfoTwoToOne, true);
        test_remote_info(remoteInfoOneToTwo, remoteInfoTwoToOne, false);

        test_xi_points(inputMesh1, middleGrid1, xiPtsOnMesh1, classification1, xiPts);
        test_xi_points(inputMesh2, middleGrid2, xiPtsOnMesh2, classification2, xiPts);

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

    void test_remote_info(mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfoOneToTwo, 
                          mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfoTwoToOne,   
                          bool sendOneToTwo)
    {
      
      auto middleGridSend = sendOneToTwo ? middleGrid1 : middleGrid2;
      auto middleGridRecv = sendOneToTwo ? middleGrid2 : middleGrid1;

      auto centroidFieldSendPtr = mesh::create_field<utils::Point>(middleGridSend,
                                                            mesh::impl::FieldShape(0, 0, 1), 1);

      auto centroidFieldRecvPtr = mesh::create_field<utils::Point>(middleGridRecv,
                                                            mesh::impl::FieldShape(0, 0, 1), 1);

      auto& centroidFieldSend = *centroidFieldSendPtr;
      for (auto el : middleGridSend->get_elements())
        if (el)
        {
          centroidFieldSend(el, 0, 0) = mesh::compute_centroid(el);
        }

      MiddleMeshFieldCommunicationSPMD<utils::Point> exchanger(middleGrid1, middleGrid2,
                                                              remoteInfoOneToTwo, remoteInfoTwoToOne);
      exchanger.start_exchange(centroidFieldSendPtr, centroidFieldRecvPtr);
      exchanger.finish_exchange(centroidFieldRecvPtr);

      auto& centroidFieldRecv = *centroidFieldRecvPtr;
      for (auto el : middleGridRecv->get_elements())
        if (el)
        {
          utils::Point centroidLocal  = mesh::compute_centroid(el);
          utils::Point centroidRemote = centroidFieldRecv(el, 0, 0);
          for (int i=0; i < 3; ++i)
            EXPECT_NEAR(centroidLocal[i], centroidRemote[i], 1e-13);
        }
    }

    void test_xi_points(std::shared_ptr<mesh::Mesh> inputMesh, std::shared_ptr<mesh::Mesh> middleGrid,
                        mesh::FieldPtr<utils::Point> xiCoordsOnInputMeshPtr,
                        mesh::FieldPtr<mesh::MeshEntityPtr> meshInToInputMeshPtr,
                        std::shared_ptr<XiCoordinates> xiCoords)
    {
      auto& xiCoordsOnInputMesh = *xiCoordsOnInputMeshPtr;
      auto& meshInToInputMesh   = *meshInToInputMeshPtr;
      for (auto& elMiddle : middleGrid->get_elements())
        if (elMiddle)
        {
          mesh::MeshEntityPtr elInput = meshInToInputMesh(elMiddle, 0, 0);
          auto& xiCoordsEl = xiCoords->get_xi_coords(elMiddle->get_type());
          std::pair<double, double> xiRangeMiddle = xiCoords->get_xi_coord_range(elMiddle->get_type());
          std::pair<double, double> xiRangeInput  = xiCoords->get_xi_coord_range(elInput->get_type());

          for (size_t i=0; i < xiCoordsEl.size(); ++i)
          {
            utils::Point xiCoordsMiddle = mesh::convert_xi_coords_from_range(xiRangeMiddle.first, xiRangeMiddle.second, xiCoordsEl[i]);
            auto middleMeshCoords = mesh::compute_coords_from_xi_3d(elMiddle, xiCoordsMiddle);

            utils::Point xiCoordsInput = mesh::convert_xi_coords_from_range(xiRangeInput.first, xiRangeInput.second,
                                                                            xiCoordsOnInputMesh(elMiddle, i, 0));
            auto inputMeshCoords  = mesh::compute_coords_from_xi_3d(elInput, xiCoordsInput);

            for (int d=0; d < 3; ++d)
              EXPECT_NEAR(middleMeshCoords[d], inputMeshCoords[d], 1e-12);
          }

        }      
    }

    MPI_Comm meshComm = MPI_COMM_NULL;
    std::shared_ptr<mesh::Mesh> inputMesh1;
    std::shared_ptr<mesh::Mesh> inputMesh2;
    std::shared_ptr<mesh::Mesh> middleGrid1;
    std::shared_ptr<mesh::Mesh> middleGrid2;
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
