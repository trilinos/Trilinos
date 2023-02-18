#include "create_mesh.h"
#include "mesh_scatter_from_root.h"
#include "util/mesh_comparer.h"
#include "util/meshes.h"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

void trim_mesh(std::shared_ptr<mesh::Mesh> meshSerial, MPI_Comm meshComm, mesh::FieldPtr<int> destRanksOnMeshComm);

void check_remotes_unique(std::shared_ptr<mesh::Mesh> mesh);

class MeshScatterFromRootTester : public ::testing::Test
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

      bool makeTriangles = false;
      bool reverseXEdges = true;
      meshSerial         = mesh::impl::create_mesh(meshspec, func, MPI_COMM_SELF, makeTriangles, reverseXEdges);
    }

    void setup(MPI_Comm inputCommArg, MPI_Comm meshCommArg, int rootRankOnInputCommArg,
               const mesh::impl::MeshSpec& spec)
    {
      auto f = [](const utils::Point& pt) { return pt; };
      setup(inputCommArg, meshCommArg, rootRankOnInputCommArg, spec, f);
    }

    void runtest(mesh::FieldPtr<int> destRanksOnInputComm, double tol)
    {
      auto inputMesh = utils::impl::comm_rank(inputComm) == rootRankOnInputComm ? meshSerial : nullptr;
      mesh::impl::MeshScatterFromRoot scatterer(inputComm, inputMesh, meshComm, destRanksOnInputComm);
      std::shared_ptr<mesh::Mesh> meshScattered = scatterer.scatter();

      if (meshComm != MPI_COMM_NULL)
      {
        trim_mesh(meshSerial, inputComm, destRanksOnInputComm);

        mesh::check_topology(meshScattered);
        compare_meshes(meshScattered, meshSerial, tol);

        if (!meshspec.xPeriodic && !meshspec.yPeriodic)
          check_remotes_unique(meshScattered);
      }
    }

    MPI_Comm inputComm;
    MPI_Comm meshComm;
    int rootRankOnInputComm;
    mesh::impl::MeshSpec meshspec;
    std::shared_ptr<mesh::Mesh> meshSerial;
};

void check_remotes_unique(std::shared_ptr<mesh::Mesh> mesh)
{
  std::set<int> remoteRanks;
  for (int dim = 0; dim < 2; ++dim)
    for (auto& entity : mesh->get_mesh_entities(dim))
      if (entity && entity->count_remote_shared_entities() > 0)
      {
        remoteRanks.clear();
        for (int i = 0; i < entity->count_remote_shared_entities(); ++i)
        {
          int remoteRank = entity->get_remote_shared_entity(i).remoteRank;
          EXPECT_EQ(remoteRanks.count(remoteRank), 0u);
          remoteRanks.insert(remoteRank);
        }
      }
}

void trim_mesh(std::shared_ptr<mesh::Mesh> meshSerial, MPI_Comm inputComm, mesh::FieldPtr<int> destRanksOnInputComm)
{
  // delete any element that has destination rank difference than my rank on mesh_comm
  int myrank           = utils::impl::comm_rank(inputComm);
  auto& destRanksField = *destRanksOnInputComm;
  for (auto& el : meshSerial->get_elements())
    if (el && destRanksField(el, 0, 0) != myrank)
      meshSerial->delete_face(el);
}

// creates a destination rank field.  All ranks greater than or
// equal to remove_rank are incremented by 1, so that remove_rank
// does not appear in the field values
mesh::FieldPtr<int> get_dest_rank_field(std::shared_ptr<mesh::Mesh> mesh, const std::vector<double>& xBoundaries,
                                        int removeRank = -1)
{
  assert(xBoundaries.size() >= 1);
  auto fieldPtr = mesh::create_field<int>(mesh, mesh::impl::FieldShape(0, 0, 1), 1);
  auto& field   = *fieldPtr;

  for (auto& el : mesh->get_elements())
    if (el)
    {
      utils::Point centroid = mesh::compute_centroid(el);
      int destRank          = -1;
      for (size_t i = 0; i < xBoundaries.size() - 1; ++i)
        if (centroid.x >= xBoundaries[i] && centroid.x < xBoundaries[i + 1])
        {
          destRank = i;
          break;
        }

      if (destRank == -1)
        throw std::runtime_error("centroid not within range");

      if (destRank >= removeRank)
        destRank++;

      field(el, 0, 0) = destRank;
    }

  return fieldPtr;
}

std::vector<double> linspace(double start, double end, int nvalues)
{
  assert(nvalues >= 2);
  double deltaX = (end - start) / (nvalues - 1);
  std::vector<double> vals;
  for (int i = 0; i < nvalues; ++i)
    vals.push_back(start + deltaX * i);

  return vals;
}

} // namespace

TEST_F(MeshScatterFromRootTester, BothCommsSame)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  for (int rootRank = 0; rootRank < 2; ++rootRank)
  {
    for (int i = 1; i < utils::impl::comm_size(MPI_COMM_WORLD); ++i)
    {
      setup(MPI_COMM_WORLD, MPI_COMM_WORLD, rootRank, spec);

      auto destRankField = get_dest_rank_field(meshSerial, linspace(0, 1, i + 1));
      runtest(destRankField, 1e-13);
    }
  }
}

TEST_F(MeshScatterFromRootTester, DifferentComms)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec = {.numelX = 8, .numelY = 8, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  for (int rootRank = 0; rootRank < 2; ++rootRank)
  {
    for (int i = 1; i < utils::impl::comm_size(MPI_COMM_WORLD); ++i)
    {
      int color = utils::impl::comm_rank(MPI_COMM_WORLD) == rootRank ? MPI_UNDEFINED : 0;
      MPI_Comm meshComm;
      MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);
      setup(MPI_COMM_WORLD, meshComm, rootRank, spec);

      auto destRankField = get_dest_rank_field(meshSerial, linspace(0, 1, i + 1), rootRank);
      runtest(destRankField, 1e-13);

      if (meshComm != MPI_COMM_NULL)
        MPI_Comm_free(&meshComm);
    }
  }
}

TEST_F(MeshScatterFromRootTester, BothCommsSameRootProc0Periodic)
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

  for (int i = 1; i < utils::impl::comm_size(MPI_COMM_WORLD); ++i)
  {
    setup(MPI_COMM_WORLD, MPI_COMM_WORLD, rootRank, spec, func);
    auto destRankField = get_dest_rank_field(meshSerial, linspace(-rOut, rOut, i + 1));
    runtest(destRankField, 1e-13);
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
