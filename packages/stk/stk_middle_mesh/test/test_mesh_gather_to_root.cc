#include "create_mesh.h"
#include "mesh_gather_to_root.h"
#include "util/mesh_comparer.h"
#include "util/meshes.h"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

class MeshGatherToRootTester : public ::testing::Test
{
  protected:
    template <typename T>
    void setup(MPI_Comm inputCommArg, MPI_Comm meshCommArg, int rootRankOnInputCommArg,
               const mesh::impl::MeshSpec& spec, T func)
    {
      inputComm           = inputCommArg;
      meshComm            = meshCommArg;
      rootRankOnInputComm = rootRankOnInputCommArg;
      meshspec            = spec;
      if (meshComm != MPI_COMM_NULL)
        meshParallel = mesh::impl::create_mesh(spec, func, meshComm); // TODO: mesh_comm is not used?

      if (utils::impl::comm_rank(inputComm) == rootRankOnInputComm)
        meshSerial = mesh::impl::create_mesh(meshspec, func, MPI_COMM_SELF);
    }

    void setup(MPI_Comm inputCommArg, MPI_Comm meshCommArg, int rootRankOnInputCommArg,
               const mesh::impl::MeshSpec& spec)
    {
      auto f = [](const utils::Point& pt) { return pt; };
      setup(inputCommArg, meshCommArg, rootRankOnInputCommArg, spec, f);
    }

    void runtest(double tol)
    {
      mesh::impl::MeshGatherToRoot gatherer(inputComm, rootRankOnInputComm, meshParallel);
      std::shared_ptr<mesh::Mesh> meshGathered = gatherer.gather();

      if (utils::impl::comm_rank(inputComm) == rootRankOnInputComm)
      {
        compare_meshes(meshGathered, meshSerial, tol);
      }

      test_element_origins(meshGathered, gatherer.get_element_origins(), tol);
    }

    struct CentroidData
    {
        int localIdOnOriginProcess;
        utils::Point centroid;
    };

    void test_element_origins(std::shared_ptr<mesh::Mesh> meshGathered,
                              mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins, double tol)
    {
      utils::impl::ParallelExchange<CentroidData> exchanger(inputComm, 75);
      if (meshGathered)
      {
        // have root process send (lid, centroid) of elements to the destination
        // described by element_origins
        auto& elementOriginsField = *elementOrigins;
        for (auto& el : meshGathered->get_elements())
          if (el)
          {
            mesh::RemoteSharedEntity remote = elementOriginsField(el, 0, 0);
            utils::Point centroid           = mesh::compute_centroid_3d(el);

            auto& buf = exchanger.get_send_buffer(remote.remoteRank);
            buf.push_back({remote.remoteId, centroid});
          }
      }

      // have receiving processes get element by lid, compute centroid, compare
      if (meshComm != MPI_COMM_NULL)
      {
        exchanger.set_recv_buffer_size(rootRankOnInputComm, count_valid(meshParallel->get_elements()));
      }

      exchanger.start_sends();
      exchanger.start_recvs();
      exchanger.complete_recvs();
      exchanger.complete_sends();

      if (meshComm != MPI_COMM_NULL)
      {
        auto& recvBuf = exchanger.get_recv_buffer(rootRankOnInputComm);

        for (auto& data : recvBuf)
        {
          mesh::MeshEntityPtr el = meshParallel->get_elements()[data.localIdOnOriginProcess];
          utils::Point centroid  = mesh::compute_centroid_3d(el);

          utils::Point disp = centroid - data.centroid;
          double dist       = std::sqrt(dot(disp, disp));
          EXPECT_NEAR(dist, 0, tol);
        }
      }
    }

    MPI_Comm inputComm;
    MPI_Comm meshComm;
    int rootRankOnInputComm;
    mesh::impl::MeshSpec meshspec;
    std::shared_ptr<mesh::Mesh> meshParallel;
    std::shared_ptr<mesh::Mesh> meshSerial;
};

} // namespace

TEST_F(MeshGatherToRootTester, BothCommsSameRootProc0)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  int rootRank              = 0;
  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  setup(MPI_COMM_WORLD, MPI_COMM_WORLD, rootRank, spec);
  runtest(1e-13);
}

TEST_F(MeshGatherToRootTester, BothCommsSameRootProc1)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  int rootRank              = 1;
  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  setup(MPI_COMM_WORLD, MPI_COMM_WORLD, rootRank, spec);
  runtest(1e-13);
}

TEST_F(MeshGatherToRootTester, SplitCommRootProc0)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  int rootRank              = 0;
  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  int color = utils::impl::comm_rank(MPI_COMM_WORLD) == utils::impl::comm_size(MPI_COMM_WORLD) - 1 ? MPI_UNDEFINED : 0;
  MPI_Comm meshComm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);
  setup(MPI_COMM_WORLD, meshComm, rootRank, spec);
  runtest(1e-13);
}

TEST_F(MeshGatherToRootTester, SplitCommRootProc1)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  int rootRank              = 1;
  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  int color = utils::impl::comm_rank(MPI_COMM_WORLD) == utils::impl::comm_size(MPI_COMM_WORLD) - 1 ? MPI_UNDEFINED : 0;
  MPI_Comm meshComm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);
  setup(MPI_COMM_WORLD, meshComm, rootRank, spec);
  runtest(1e-13);
}

TEST_F(MeshGatherToRootTester, SplitCommRootProc0ReversedRanks)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  int rootRank              = 0;
  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  int color = utils::impl::comm_rank(MPI_COMM_WORLD) == utils::impl::comm_size(MPI_COMM_WORLD) - 1 ? MPI_UNDEFINED : 0;
  int key   = utils::impl::comm_size(MPI_COMM_WORLD) - utils::impl::comm_rank(MPI_COMM_WORLD);
  MPI_Comm meshComm;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &meshComm);
  setup(MPI_COMM_WORLD, meshComm, rootRank, spec);
  runtest(1e-13);
}

TEST_F(MeshGatherToRootTester, SplitCommNotIncludingRootProc0)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  int rootRank              = 0;
  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  int color = utils::impl::comm_rank(MPI_COMM_WORLD) == rootRank ? MPI_UNDEFINED : 0;
  MPI_Comm meshComm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);
  setup(MPI_COMM_WORLD, meshComm, rootRank, spec);
  runtest(1e-13);
}

TEST_F(MeshGatherToRootTester, SplitCommNotIncludingRootProc1)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  int rootRank              = 1;
  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  int color = utils::impl::comm_rank(MPI_COMM_WORLD) == rootRank ? MPI_UNDEFINED : 0;
  MPI_Comm meshComm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);
  setup(MPI_COMM_WORLD, meshComm, rootRank, spec);
  runtest(1e-13);
}

TEST_F(MeshGatherToRootTester, BothCommsSameRootProc0Periodic)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  int rootRank = 0;

  int nelemR     = 4;
  int nelemTheta = 4;
  double rIn     = 0.5;
  double rOut    = 1.5;
  mesh::impl::MeshSpec spec;
  spec.numelX    = nelemR;
  spec.numelY    = nelemTheta;
  spec.xmin      = rIn;
  spec.xmax      = rOut;
  spec.ymin      = 0;
  spec.ymax      = 2 * PI;
  spec.yPeriodic = true;

  auto func = [&](const utils::Point& pt) {
    // interpret x and y as r and theta
    double r     = pt.x;
    double theta = pt.y;

    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    double z = 0;
    utils::Point pt2(x, y, z);
    return pt2;
  };

  setup(MPI_COMM_WORLD, MPI_COMM_WORLD, rootRank, spec, func);
  runtest(1e-13);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
