#include "util/mesh_comparer.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

MeshComparer::MeshComparer(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2, double tol)
  : m_mesh1(mesh1)
  , m_mesh2(mesh2)
  , m_tol(tol)
{}

void MeshComparer::compare()
{
  compare_entity_counts();
  compare_centroids();
}

void MeshComparer::compare_entity_counts()
{
  for (int dim = 0; dim < 3; ++dim)
    EXPECT_EQ(mesh::count_valid(m_mesh1->get_mesh_entities(dim)), mesh::count_valid(m_mesh2->get_mesh_entities(dim)));
}

void MeshComparer::compare_centroids()
{
  for (int dim = 0; dim <= 2; ++dim)
    for (auto& mesh1Entity : m_mesh1->get_mesh_entities(dim))
      if (mesh1Entity)
      {
        utils::Point mesh1Centroid = mesh::compute_centroid(mesh1Entity); // TODO: generalize to edges, verts
        auto mesh2Entity           = find_closest_mesh2_entity(dim, mesh1Centroid);
        utils::Point mesh2Centroid = mesh::compute_centroid(mesh2Entity);

        EXPECT_NEAR(compute_dist(mesh1Centroid, mesh2Centroid), 0, m_tol);
        EXPECT_EQ(mesh1Entity->get_type(), mesh2Entity->get_type());
        EXPECT_EQ(mesh1Entity->count_down(), mesh2Entity->count_down());
        EXPECT_EQ(mesh1Entity->count_up(), mesh2Entity->count_up());

        if (dim > 0)
          compare_downward(mesh1Entity, mesh2Entity);

        if (dim < 2)
          compare_upward(mesh1Entity, mesh2Entity);
      }
}

mesh::MeshEntityPtr MeshComparer::find_closest_mesh2_entity(int dim, const utils::Point& centroid)
{
  double minDist                = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr minEntity = nullptr;
  for (auto& mesh2Entity : m_mesh2->get_mesh_entities(dim))
    if (mesh2Entity)
    {
      auto mesh2Centroid = mesh::compute_centroid(mesh2Entity);
      double dist        = compute_dist(centroid, mesh2Centroid);
      if (dist < minDist)
      {
        minDist   = dist;
        minEntity = mesh2Entity;
      }
    }

  return minEntity;
}

void MeshComparer::compare_downward(mesh::MeshEntityPtr mesh1Entity, mesh::MeshEntityPtr mesh2Entity)
{
  EXPECT_EQ(mesh1Entity->count_down(), mesh2Entity->count_down());
  for (int i = 0; i < mesh1Entity->count_down(); ++i)
  {
    mesh::MeshEntityPtr mesh1DownEntity = mesh1Entity->get_down(i);
    mesh::MeshEntityPtr mesh2DownEntity = mesh2Entity->get_down(i);

    double dist = compute_dist(mesh::compute_centroid(mesh1DownEntity), mesh::compute_centroid(mesh2DownEntity));
    EXPECT_NEAR(dist, 0, m_tol);
    EXPECT_EQ(mesh1Entity->get_down_orientation(i), mesh2Entity->get_down_orientation(i));
  }
}

void MeshComparer::compare_upward(mesh::MeshEntityPtr mesh1Entity, mesh::MeshEntityPtr mesh2Entity)
{
  EXPECT_EQ(mesh1Entity->count_up(), mesh2Entity->count_up());
  for (int i = 0; i < mesh1Entity->count_up(); ++i)
  {
    mesh::MeshEntityPtr mesh1UpEntity = mesh1Entity->get_up(i);
    utils::Point mesh1UpCentroid      = mesh::compute_centroid(mesh1UpEntity);

    mesh::MeshEntityPtr mesh2UpEntity = find_upward_entity(mesh2Entity, mesh1UpCentroid);

    double dist = compute_dist(mesh1UpCentroid, mesh::compute_centroid(mesh2UpEntity));
    EXPECT_NEAR(dist, 0, m_tol);
  }
}

mesh::MeshEntityPtr MeshComparer::find_upward_entity(mesh::MeshEntityPtr entity, const utils::Point& centroid)
{
  double minDist                = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr minEntity = nullptr;
  for (int i = 0; i < entity->count_up(); ++i)
  {
    double dist = compute_dist(mesh::compute_centroid(entity->get_up(i)), centroid);
    if (dist < minDist)
    {
      minDist   = dist;
      minEntity = entity->get_up(i);
    }
  }

  return minEntity;
}

double MeshComparer::compute_dist(const utils::Point& pt1, const utils::Point& pt2)
{
  auto disp = pt1 - pt2;
  return std::sqrt(dot(disp, disp));
}

void compare_meshes(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2, double tol)
{
  MeshComparer comparer(mesh1, mesh2, tol);
  comparer.compare();
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
