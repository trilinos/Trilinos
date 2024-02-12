#include "gtest/gtest.h"
#include "stk_middle_mesh/predicates/average_normal_field.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_scatter.hpp"

using namespace stk::middle_mesh;

TEST(AveragedNormalField, PlaneQuad)
{
  mesh::impl::MeshSpec spec;
  spec.xmin = 0; spec.ymin = 0;
  spec.xmax = 1; spec.ymax = 1;
  spec.numelX = 4, spec.numelY = 4;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh = mesh::impl::create_mesh(spec, f);

  predicates::impl::AveragedNormalField averagedNormalField(mesh);
  auto& normalField = *(averagedNormalField.get_field());

  double edgeLength = (spec.xmax - spec.xmin)/spec.numelX;
  for (auto& vert : mesh->get_vertices())
  {
    if (vert)
    {
      utils::Point normal = normalField(vert, 0, 0);
      
      EXPECT_NEAR(normal.x, 0.0, 1e-13);
      EXPECT_NEAR(normal.y, 0.0, 1e-13);
      EXPECT_NEAR(normal.z, edgeLength, 1e-13);
    }
  }

}

TEST(AveragedNormalField, PlaneTriangle)
{
  mesh::impl::MeshSpec spec;
  spec.xmin = 0; spec.ymin = 0;
  spec.xmax = 1; spec.ymax = 1;
  spec.numelX = 4, spec.numelY = 4;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh = mesh::impl::create_mesh(spec, f, MPI_COMM_WORLD, true);

  predicates::impl::AveragedNormalField averagedNormalField(mesh);
  auto& normalField = *(averagedNormalField.get_field());

  double deltaX = (spec.xmax - spec.xmin)/spec.numelX;
  double xyEdgeLength = deltaX;
  double diagonalEdgeLength = std::sqrt(2*xyEdgeLength*xyEdgeLength);  
  for (auto& vert : mesh->get_vertices())
  {
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      bool isOnXLowerBoundary = std::abs(pt.x - spec.xmin) < 1e-13;
      bool isOnXUpperBoundary = std::abs(pt.x - spec.xmax) < 1e-13;
      bool isOnYLowerBoundary = std::abs(pt.y - spec.ymin) < 1e-13;
      bool isOnYUpperBoundary = std::abs(pt.y - spec.ymax) < 1e-13;

      utils::Point normal = normalField(vert, 0, 0);

      double expectedEdgeLength = 0;
      if ((isOnXLowerBoundary && isOnYLowerBoundary) ||
          (isOnXUpperBoundary && isOnYUpperBoundary))
      {
        expectedEdgeLength = (2*xyEdgeLength + diagonalEdgeLength)/3;
      } else if ( (isOnXLowerBoundary && isOnYUpperBoundary) ||
                  (isOnXUpperBoundary && isOnYLowerBoundary))
      {
        expectedEdgeLength = 2*xyEdgeLength/2;
      } else if (isOnXLowerBoundary || isOnXUpperBoundary || isOnYLowerBoundary || isOnYUpperBoundary)
      {
        expectedEdgeLength = (3*xyEdgeLength + diagonalEdgeLength)/4;
      } else
      {
        expectedEdgeLength = (4*xyEdgeLength + 2*diagonalEdgeLength)/6;
      }
      
      EXPECT_NEAR(normal.x, 0.0, 1e-13);
      EXPECT_NEAR(normal.y, 0.0, 1e-13);
      EXPECT_NEAR(normal.z, expectedEdgeLength, 1e-13);
    }
  }
}

TEST(AveragedNormalField, ScatteredMesh)
{
  mesh::impl::MeshSpec spec;
  spec.xmin = 0; spec.ymin = 0;
  spec.xmax = 1; spec.ymax = 1;
  spec.numelX = 4, spec.numelY = 4;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh = mesh::impl::create_mesh(spec, f);

  int commsize = utils::impl::comm_size(MPI_COMM_WORLD);
  int commrank = utils::impl::comm_rank(MPI_COMM_WORLD);
  auto meshScatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(MPI_COMM_WORLD, mesh);
  for (auto& el : mesh->get_elements())
    if (el)
    {
      int destRank1 = commrank;
      int destRank2 = (commrank + 1) % commsize;
      meshScatterSpec->add_destination(el, destRank1);
      meshScatterSpec->add_destination(el, destRank2);
    }

  mesh::impl::MeshScatter scatterer(meshScatterSpec, mesh, MPI_COMM_WORLD);
  auto meshScattered = scatterer.scatter();

  predicates::impl::AveragedNormalField averagedNormalField(meshScattered);
  auto& normalField = *(averagedNormalField.get_field());

  double edgeLength = (spec.xmax - spec.xmin)/spec.numelX;
  for (auto& vert : meshScattered->get_vertices())
  {
    if (vert)
    {
      utils::Point normal = normalField(vert, 0, 0);
      
      EXPECT_NEAR(normal.x, 0.0, 1e-13);
      EXPECT_NEAR(normal.y, 0.0, 1e-13);
      EXPECT_NEAR(normal.z, edgeLength, 1e-13);
    }
  }

}