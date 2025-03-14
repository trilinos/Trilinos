#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/mesh_entity.hpp"
#include "stk_middle_mesh/mesh_scatter_from_root.hpp"
#include "stk_middle_mesh/parallel_exchange.hpp"
#include "util/mesh_comparer.hpp"
#include "util/meshes.hpp"
#include "gtest/gtest.h"
#include "stk_middle_mesh/variable_size_field.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

void trim_mesh(std::shared_ptr<mesh::Mesh> meshSerial, MPI_Comm meshComm, mesh::FieldPtr<int> destRanksOnMeshComm);


void check_entity_destinations(MPI_Comm unionComm, int rootRankOnUnionComm,
                               std::shared_ptr<mesh::Mesh> inputMesh,
                               std::shared_ptr<mesh::Mesh> meshScattered,
                               mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> entityDestinations);

class MeshScatterFromRootTester : public ::testing::Test
{
  protected:
    template <typename T>
    void setup(MPI_Comm unionCommArg, MPI_Comm meshCommArg, int rootRankOnUnionCommArg,
               const mesh::impl::MeshSpec& spec, T func)
    {
      unionComm           = unionCommArg;
      meshComm            = meshCommArg;
      rootRankOnUnionComm = rootRankOnUnionCommArg;
      meshspec            = spec;

      bool makeTriangles = false;
      bool reverseXEdges = true;
      meshSerial         = mesh::impl::create_mesh(meshspec, func, MPI_COMM_SELF, makeTriangles, reverseXEdges);
    }

    void setup(MPI_Comm unionCommArg, MPI_Comm meshCommArg, int rootRankOnUnionCommArg,
               const mesh::impl::MeshSpec& spec)
    {
      auto f = [](const utils::Point& pt) { return pt; };
      setup(unionCommArg, meshCommArg, rootRankOnUnionCommArg, spec, f);
    }

    void runtest(mesh::FieldPtr<int> destRanksOnInputComm, double tol)
    {
      auto inputMesh = utils::impl::comm_rank(unionComm) == rootRankOnUnionComm ? meshSerial : nullptr;
      mesh::impl::MeshScatterFromRoot scatterer(unionComm, inputMesh, meshComm, destRanksOnInputComm);
      std::shared_ptr<mesh::Mesh> meshScattered = scatterer.scatter();
      mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> entityDestinations = scatterer.get_entity_destinations();

      check_entity_destinations(unionComm, rootRankOnUnionComm,
                                inputMesh, meshScattered, entityDestinations);

      if (meshComm != MPI_COMM_NULL)
      {
        trim_mesh(meshSerial, unionComm, destRanksOnInputComm);

        mesh::check_topology(meshScattered);
        compare_meshes(meshScattered, meshSerial, tol);
      }
    }

    MPI_Comm unionComm;
    MPI_Comm meshComm;
    int rootRankOnUnionComm;
    mesh::impl::MeshSpec meshspec;
    std::shared_ptr<mesh::Mesh> meshSerial;
};


void trim_mesh(std::shared_ptr<mesh::Mesh> meshSerial, MPI_Comm unionComm, mesh::FieldPtr<int> destRanksOnInputComm)
{
  // delete any element that has destination rank difference than my rank on mesh_comm
  int myrank           = utils::impl::comm_rank(unionComm);
  auto& destRanksField = *destRanksOnInputComm;
  for (auto& el : meshSerial->get_elements())
    if (el && destRanksField(el, 0, 0) != myrank)
      meshSerial->delete_face(el);
}

// creates a destination rank field.  All ranks greater than or
// equal to remove_rank are incremented by 1, so that remove_rank
// does not appear in the field values
mesh::FieldPtr<int> get_dest_rank_field(std::shared_ptr<mesh::Mesh> mesh, const std::vector<double>& xBoundaries,
                                        int removeRank = 999999)
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

struct EntityCentroid
{
  int localId;
  utils::Point pt;
};

void check_entity_destinations(MPI_Comm unionComm, int rootRankOnUnionComm,
                               std::shared_ptr<mesh::Mesh> inputMesh,
                               std::shared_ptr<mesh::Mesh> meshScattered,
                               mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> entityDestinations)
{

  for (int dim=0; dim < 3; ++dim)
  {
    utils::impl::ParallelExchange<EntityCentroid> exchanger(unionComm, 2024 + dim);

    if (inputMesh)
    {
      for (auto entity : inputMesh->get_mesh_entities(dim))
      {
        utils::Point centroid = mesh::compute_centroid(entity);
        for (int i=0; i < entityDestinations->get_num_comp(entity, 0); ++i)
        {
          mesh::RemoteSharedEntity remote = (*entityDestinations)(entity, 0, i);
          exchanger.get_send_buffer(remote.remoteRank).push_back(EntityCentroid{remote.remoteId, centroid});
        }
      }
    }

    if (meshScattered)
    {
      exchanger.set_recv_buffer_size(rootRankOnUnionComm, mesh::count_valid(meshScattered->get_mesh_entities(dim)));
    }

    exchanger.start_sends();
    exchanger.start_recvs();
    exchanger.complete_recvs();

    if (meshScattered)
    {
      for (auto& centroid : exchanger.get_recv_buffer(rootRankOnUnionComm))
      {
        mesh::MeshEntityPtr entity = meshScattered->get_mesh_entities(dim)[centroid.localId];
        utils::Point localCentroid = mesh::compute_centroid(entity);

        for (int i=0; i < 3; ++i)
          EXPECT_NEAR(centroid.pt[i], localCentroid[i], 1e-13);
      }
    }

    exchanger.complete_sends();
  }
}

} // namespace

TEST_F(MeshScatterFromRootTester, BothCommsSame)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  int commSize = utils::impl::comm_size(MPI_COMM_WORLD);
  mesh::impl::MeshSpec spec = {.numelX = 4, .numelY = 4, .xmin = 0, .xmax = 1, .ymin = 0, .ymax = 1};

  for (int rootRank = 0; rootRank < commSize; ++rootRank)
  {
    setup(MPI_COMM_WORLD, MPI_COMM_WORLD, rootRank, spec);

    auto destRankField = get_dest_rank_field(meshSerial, linspace(0, 1, commSize + 1));
    runtest(destRankField, 1e-13);
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
      MPI_Comm meshCommTest;
      MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshCommTest);
      setup(MPI_COMM_WORLD, meshCommTest, rootRank, spec);

      auto destRankField = get_dest_rank_field(meshSerial, linspace(0, 1, i + 1), rootRank);
      runtest(destRankField, 1e-13);

      if (meshCommTest != MPI_COMM_NULL)
        MPI_Comm_free(&meshCommTest);
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
    auto destRankField = get_dest_rank_field(meshSerial, linspace(-rOut, rOut, i + 2));
    runtest(destRankField, 1e-13);
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
