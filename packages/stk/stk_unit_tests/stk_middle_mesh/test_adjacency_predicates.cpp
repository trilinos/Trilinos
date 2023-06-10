#include "gtest/gtest.h"

#include "stk_middle_mesh/predicates/adjacency_predicates.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(AdjacencyPredicates, any_edge_intersect)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(2, 0);
  auto v3   = mesh->create_vertex(2, 2);
  auto v4   = mesh->create_vertex(0, 2);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);
  mesh::MeshEntityPtr edges1[mesh::MAX_DOWN];
  mesh::get_downward(el, 1, edges1);

  auto v5 = mesh->create_vertex(1, 1);
  auto v6 = mesh->create_vertex(1.5, 0.5);
  auto v7 = mesh->create_vertex(3, 1);
  auto v8 = mesh->create_vertex(1, 3);
  auto e1 = mesh->create_edge(v5, v6);
  auto e2 = mesh->create_edge(v6, v7);
  mesh->create_edge(v7, v8);
  auto e4 = mesh->create_edge(v5, v8);
  mesh->create_quad_from_verts(v5, v6, v7, v8);

  predicates::impl::AdjacencyPredicates preds(mesh);

  EXPECT_EQ(preds.any_edge_intersect(el, e2), edges1[1]);
  EXPECT_EQ(preds.any_edge_intersect(el, e4), edges1[2]);
  EXPECT_EQ(preds.any_edge_intersect(el, e1), nullptr);
}

TEST(AdjacencyPredicates, any_edges_intersect)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(2, 0);
  auto v3   = mesh->create_vertex(2, 2);
  auto v4   = mesh->create_vertex(0, 2);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);

  auto v5  = mesh->create_vertex(5, 0);
  auto v6  = mesh->create_vertex(7, 0);
  auto v7  = mesh->create_vertex(7, 2);
  auto v8  = mesh->create_vertex(5, 2);
  auto el1 = mesh->create_quad_from_verts(v5, v6, v7, v8);

  // this it tricky: any_edges_intersect *any* edge of el2 that intersect with
  // el1, no guarantee which one it is
  auto v9  = mesh->create_vertex(1, 1);
  auto v10 = mesh->create_vertex(1.5, 1);
  auto v11 = mesh->create_vertex(1.5, 3);
  auto v12 = mesh->create_vertex(1, 3);
  auto el2 = mesh->create_quad_from_verts(v9, v10, v11, v12);
  mesh::MeshEntityPtr edges2[mesh::MAX_DOWN];
  mesh::get_downward(el2, 1, edges2);
  predicates::impl::AdjacencyPredicates preds(mesh);

  EXPECT_EQ(preds.any_edges_intersect(el, el1), nullptr);
  auto edgeIntersect = preds.any_edges_intersect(el, el2);
  EXPECT_TRUE(edgeIntersect == edges2[1] || edgeIntersect == edges2[3]);
}

TEST(AdjacencyPredicates, any_vertices_contained)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(2, 0);
  auto v3   = mesh->create_vertex(2, 2);
  auto v4   = mesh->create_vertex(0, 2);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);

  auto v5  = mesh->create_vertex(5, 0);
  auto v6  = mesh->create_vertex(7, 0);
  auto v7  = mesh->create_vertex(7, 2);
  auto v8  = mesh->create_vertex(5, 2);
  auto el1 = mesh->create_quad_from_verts(v5, v6, v7, v8);

  // this it tricky: any_edges_intersect *any* edge of el2 that intersect with
  // el1, no guarantee which one it is
  auto v9  = mesh->create_vertex(1, 1);
  auto v10 = mesh->create_vertex(3, 1);
  auto v11 = mesh->create_vertex(3, 3);
  auto v12 = mesh->create_vertex(1, 3);
  auto el2 = mesh->create_quad_from_verts(v9, v10, v11, v12);
  mesh::MeshEntityPtr edges2[mesh::MAX_DOWN];
  mesh::get_downward(el2, 1, edges2);
  predicates::impl::AdjacencyPredicates preds(mesh);

  EXPECT_EQ(preds.any_vertices_contained(el, el1), nullptr);
  EXPECT_EQ(preds.any_vertices_contained(el, el2), v9);
}

TEST(AdjacencyPredicates, is_point_contained)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(2, 0);
  auto v3   = mesh->create_vertex(2, 2);
  auto v4   = mesh->create_vertex(0, 2);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);
  predicates::impl::AdjacencyPredicates preds(mesh);

  auto v5 = utils::Point(1, 1);
  auto v6 = utils::Point(0, 1);
  auto v7 = utils::Point(0, 0);
  auto v8 = utils::Point(1, 3);

  EXPECT_TRUE(preds.is_point_contained(el, v5));
  EXPECT_TRUE(preds.is_point_contained(el, v6));
  EXPECT_TRUE(preds.is_point_contained(el, v7));
  EXPECT_FALSE(preds.is_point_contained(el, v8));
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
