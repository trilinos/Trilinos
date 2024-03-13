#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/gathered_mesh_coordinate_update.hpp"
#include "stk_middle_mesh/mesh_gather_to_root.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

class GatheredMeshCoordinateUpdateTester : public ::testing::Test
{
  protected:
    void setup(MPI_Comm unionComm, int rootRankOnUnionComm, MPI_Comm meshComm, const mesh::impl::MeshSpec& spec)
    {
      parallelMesh   = nullptr;
      serialMesh     = nullptr;
      elementOrigins = nullptr;

      // create a parallel mesh on mesh_comm
      if (meshComm != MPI_COMM_NULL)
      {
        auto f       = [](const utils::Point& pt) { return pt; };
        parallelMesh = mesh::impl::create_mesh(spec, f, meshComm);
      }

      // gather it to root on union_comm
      mesh::impl::MeshGatherToRoot gatherer(unionComm, rootRankOnUnionComm, parallelMesh);
      serialMesh     = gatherer.gather();
      elementOrigins = gatherer.get_element_origins();
    }

    std::shared_ptr<mesh::Mesh> parallelMesh;
    std::shared_ptr<mesh::Mesh> serialMesh;
    mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins;
};

mesh::FieldPtr<utils::Point> copy_coordinate_field(std::shared_ptr<mesh::Mesh> mesh)
{
  auto coordsOrig = mesh::create_field<utils::Point>(mesh, mesh::impl::FieldShape(1, 0, 0), 1);
  for (auto& v : mesh->get_vertices())
    (*coordsOrig)(v, 0, 0) = v->get_point_orig(0);

  return coordsOrig;
}

void scale_coordinate_field(std::shared_ptr<mesh::Mesh> mesh, const std::array<double, 3>& factors)
{
  for (auto& v : mesh->get_vertices())
  {
    utils::Point pt = v->get_point_orig(0);
    pt.x *= factors[0];
    pt.y *= factors[1];
    pt.z *= factors[2];
    v->set_point_orig(0, pt);
  }
}

void check_coordinate_field_scaled(std::shared_ptr<mesh::Mesh> mesh, const std::array<double, 3>& factors,
                                   mesh::FieldPtr<utils::Point> coordsOrig)
{
  for (auto& v : mesh->get_vertices())
  {
    utils::Point ptNew  = v->get_point_orig(0);
    utils::Point ptOrig = (*coordsOrig)(v, 0, 0);

    EXPECT_NEAR(ptNew.x, ptOrig.x * factors[0], 1e-13);
    EXPECT_NEAR(ptNew.y, ptOrig.y * factors[1], 1e-13);
    EXPECT_NEAR(ptNew.z, ptOrig.z * factors[2], 1e-13);
  }
}

} // namespace

TEST_F(GatheredMeshCoordinateUpdateTester, RootOnMeshComm)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  MPI_Comm unionComm = MPI_COMM_WORLD;
  mesh::impl::MeshSpec spec;
  spec.xmin                     = 0;
  spec.xmax                     = 1;
  spec.ymin                     = 0;
  spec.ymax                     = 1;
  spec.numelX                   = 4;
  spec.numelY                   = 4;
  std::array<double, 3> factors = {2, 3, 4};

  for (int rootRank = 0; rootRank < std::min(2, utils::impl::comm_size(unionComm)); ++rootRank)
    for (int meshCommSize = 1; meshCommSize <= utils::impl::comm_size(unionComm); ++meshCommSize)
    {
      MPI_Comm meshComm;
      int color = utils::impl::comm_rank(unionComm) < meshCommSize ? 0 : MPI_UNDEFINED;
      MPI_Comm_split(unionComm, color, 0, &meshComm);

      setup(unionComm, rootRank, meshComm, spec);

      mesh::FieldPtr<utils::Point> coordsOrig;
      if (color == 0)
      {
        coordsOrig = copy_coordinate_field(parallelMesh);
      }

      if (serialMesh)
        scale_coordinate_field(serialMesh, factors);

      mesh::impl::GatheredMeshCoordinateUpdate updator(unionComm, serialMesh, elementOrigins, parallelMesh);
      updator.update();

      if (color == 0)
        check_coordinate_field_scaled(parallelMesh, factors, coordsOrig);

      if (color == 0)
        MPI_Comm_free(&meshComm);
    }
}

TEST_F(GatheredMeshCoordinateUpdateTester, RootNotOnMeshComm)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) < 2 || utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  MPI_Comm unionComm = MPI_COMM_WORLD;
  mesh::impl::MeshSpec spec;
  spec.xmin                     = 0;
  spec.xmax                     = 1;
  spec.ymin                     = 0;
  spec.ymax                     = 1;
  spec.numelX                   = 4;
  spec.numelY                   = 4;
  std::array<double, 3> factors = {2, 3, 4};

  for (int rootRank = 0; rootRank < std::min(2, utils::impl::comm_size(unionComm)); ++rootRank)
    for (int meshCommSize = 1; meshCommSize <= utils::impl::comm_size(unionComm) - 1; ++meshCommSize)
    {
      MPI_Comm meshComm;

      int color  = 0;
      int myrank = utils::impl::comm_rank(unionComm);
      if (rootRank == 0)
      {
        color = (myrank >= 1 && myrank <= meshCommSize) ? 0 : MPI_UNDEFINED;
      } else if (rootRank == 1)
      {
        color = (myrank <= meshCommSize + 1 && myrank != rootRank) ? 0 : MPI_UNDEFINED;
      } else
        throw std::runtime_error("root rank can only be 0 or 1");

      MPI_Comm_split(unionComm, color, 0, &meshComm);

      setup(unionComm, rootRank, meshComm, spec);

      mesh::FieldPtr<utils::Point> coordsOrig;
      if (color == 0)
      {
        coordsOrig = copy_coordinate_field(parallelMesh);
      }

      if (serialMesh)
        scale_coordinate_field(serialMesh, factors);

      mesh::impl::GatheredMeshCoordinateUpdate updator(unionComm, serialMesh, elementOrigins, parallelMesh);
      updator.update();

      if (color == 0)
        check_coordinate_field_scaled(parallelMesh, factors, coordsOrig);

      if (color == 0)
        MPI_Comm_free(&meshComm);
    }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
