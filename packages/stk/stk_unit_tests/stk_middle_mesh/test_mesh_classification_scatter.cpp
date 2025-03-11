#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_classification_scatter.hpp"
#include "stk_middle_mesh/mesh_entity.hpp"
#include "stk_middle_mesh/mesh_gather_to_root.hpp"
#include "stk_middle_mesh/mesh_scatter_from_root.hpp"
#include "gtest/gtest.h"
#include "stk_middle_mesh/variable_size_field.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

class MeshClassificationScatterTester : public ::testing::Test
{
  protected:
    void setup(MPI_Comm inputCommArg, MPI_Comm meshComm, int rootRankOnInputCommArg,
               const mesh::impl::MeshSpec& meshspec, int /*meshCommSize*/)
    {
      inputGridSerial               = nullptr;
      inputGridParallel             = nullptr;
      middleGridSerial              = nullptr;
      middleGridParallel            = nullptr;
      inputGridSerialElementOrigins = nullptr;

      inputComm           = inputCommArg;
      rootRankOnInputComm = rootRankOnInputCommArg;

      auto f = [](const utils::Point& pt) { return pt; };

      if (meshComm != MPI_COMM_NULL)
      {
        inputGridParallel = mesh::impl::create_mesh(meshspec, f, meshComm);
      }

      mesh::impl::MeshGatherToRoot gather(inputComm, rootRankOnInputComm, inputGridParallel);
      inputGridSerial               = gather.gather();
      inputGridSerialElementOrigins = gather.get_element_origins();

      mesh::FieldPtr<int> middleGridDestRanks = nullptr;
      if (utils::impl::comm_rank(inputComm) == rootRankOnInputComm)
      {
        middleGridSerial = mesh::impl::create_mesh(meshspec, f, MPI_COMM_SELF);
        middleGridDestRanks =
            get_element_destinations(middleGridSerial, inputGridSerial, inputGridSerialElementOrigins);
      } else
      {
        middleGridSerial    = nullptr;
        middleGridDestRanks = nullptr;
      }

      mesh::impl::MeshScatterFromRoot scatterMiddleGrid(inputComm, middleGridSerial, meshComm, middleGridDestRanks);
      middleGridParallel = scatterMiddleGrid.scatter();
      middleGridSerialEntityDestinations = scatterMiddleGrid.get_entity_destinations();
    }

    void runtest(double tol)
    {
      mesh::FieldPtr<mesh::MeshEntityPtr> middleGridClassificationSerial;
      if (inputGridSerial)
        middleGridClassificationSerial = create_classification(inputGridSerial, middleGridSerial);

      mesh::impl::MeshClassificationScatter scatter(inputComm, middleGridSerial, middleGridSerialEntityDestinations,
                                                    middleGridParallel,
                                                    inputGridParallel,
                                                    inputGridSerialElementOrigins, middleGridClassificationSerial);
      auto middleGridClassificationParallel = scatter.scatter();

      if (middleGridParallel)
        check_parallel_classification(middleGridClassificationParallel, tol);
    }

    mesh::FieldPtr<int> get_element_destinations(std::shared_ptr<mesh::Mesh> middleGridSerialArg,
                                                 std::shared_ptr<mesh::Mesh> inputGridSerialArg,
                                                 mesh::FieldPtr<mesh::RemoteSharedEntity> inputGridSerialElementOriginsArg)
    {
      auto classification       = create_classification(inputGridSerialArg, middleGridSerialArg);
      auto& classificationField = *classification;

      auto destRanks       = mesh::create_field<int>(middleGridSerialArg, mesh::impl::FieldShape(0, 0, 1), 1);
      auto& elementOrigins = *inputGridSerialElementOriginsArg;
      for (auto& middleGridEl : middleGridSerialArg->get_elements())
        if (middleGridEl)
        {
          mesh::MeshEntityPtr inputGridEl  = classificationField(middleGridEl, 0, 0);
          (*destRanks)(middleGridEl, 0, 0) = elementOrigins(inputGridEl, 0, 0).remoteRank;
        }

      return destRanks;
    }

    mesh::FieldPtr<mesh::MeshEntityPtr> create_classification(std::shared_ptr<mesh::Mesh> inputGridSerialArg,
                                                              std::shared_ptr<mesh::Mesh> middleGridSerialArg)
    {
      auto classificationPtr =
          mesh::create_field<mesh::MeshEntityPtr>(middleGridSerialArg, mesh::impl::FieldShape(0, 0, 1), 1);
      auto& classificationField = *classificationPtr;

      for (auto& middleGridEl : middleGridSerialArg->get_elements())
      {
        if (middleGridEl)
        {
          utils::Point middleGridElCentroid  = mesh::compute_centroid(middleGridEl);
          double minDist                     = std::numeric_limits<double>::max();
          mesh::MeshEntityPtr minInputGridEl = nullptr;
          for (auto& inputGridEl : inputGridSerialArg->get_elements())
            if (inputGridEl)
            {
              utils::Point inputGridElCentroid = mesh::compute_centroid(inputGridEl);
              utils::Point disp                = middleGridElCentroid - inputGridElCentroid;
              double dist                      = std::sqrt(dot(disp, disp));
              if (dist < minDist)
              {
                minDist        = dist;
                minInputGridEl = inputGridEl;
              }
            }

          EXPECT_NEAR(minDist, 0, 1e-13);
          classificationField(middleGridEl, 0, 0) = minInputGridEl;
        }
      }

      return classificationPtr;
    }

    void check_parallel_classification(mesh::FieldPtr<mesh::MeshEntityPtr> middleGridClassificationParallel, double tol)
    {
      auto& classificationField = *middleGridClassificationParallel;
      for (auto& middleGridEl : middleGridParallel->get_elements())
        if (middleGridEl)
        {
          mesh::MeshEntityPtr inputGridEl   = classificationField(middleGridEl, 0, 0);
          utils::Point middleGridElCentroid = mesh::compute_centroid(middleGridEl);
          utils::Point inputGridElCentroid  = mesh::compute_centroid(inputGridEl);

          utils::Point disp = middleGridElCentroid - inputGridElCentroid;
          double dist       = std::sqrt(dot(disp, disp));

          EXPECT_NEAR(dist, 0, tol);
        }
    }

    MPI_Comm inputComm;
    int rootRankOnInputComm;
    std::shared_ptr<mesh::Mesh> middleGridSerial;
    std::shared_ptr<mesh::Mesh> middleGridParallel;
    std::shared_ptr<mesh::Mesh> inputGridSerial;
    std::shared_ptr<mesh::Mesh> inputGridParallel;
    mesh::FieldPtr<mesh::RemoteSharedEntity> inputGridSerialElementOrigins;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> middleGridSerialEntityDestinations;
};

} // namespace

TEST_F(MeshClassificationScatterTester, BothCommsSame)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};
  int maxRoot               = std::min(2, utils::impl::comm_size(MPI_COMM_WORLD));

  for (int rootRank = 0; rootRank < maxRoot; ++rootRank)
  {
    setup(MPI_COMM_WORLD, MPI_COMM_WORLD, rootRank, spec, utils::impl::comm_size(MPI_COMM_WORLD));
    runtest(1e-13);
  }
}

TEST_F(MeshClassificationScatterTester, DifferentComms)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  for (int rootRank = 0; rootRank < 2; ++rootRank)
  {
    int color = utils::impl::comm_rank(MPI_COMM_WORLD) == rootRank ? MPI_UNDEFINED : 0;
    MPI_Comm meshComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);
    setup(MPI_COMM_WORLD, meshComm, rootRank, spec, utils::impl::comm_size(MPI_COMM_WORLD) - 1);

    runtest(1e-13);

    // if (mesh_comm != MPI_COMM_NULL)
    //   MPI_Comm_free(&mesh_comm);
  }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
