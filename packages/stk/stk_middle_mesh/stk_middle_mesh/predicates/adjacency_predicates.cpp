#include "stk_middle_mesh/predicates/adjacency_predicates.hpp"
#include "stk_middle_mesh/mesh_entity.hpp"
#include "stk_middle_mesh/predicates/edge_intersection_primitive.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

mesh::MeshEntityPtr AdjacencyPredicates::any_edge_intersect(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr edge2)
{
  assert(get_type_dimension(el1->get_type()) == 2);
  assert(get_type_dimension(edge2->get_type()) == 1);

  for (int i = 0; i < el1->count_down(); ++i)
    if (do_edges_intersect(el1->get_down(i), edge2))
      return el1->get_down(i);

  return nullptr;
}

mesh::MeshEntityPtr AdjacencyPredicates::any_edges_intersect(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2)
{
  assert(get_type_dimension(el1->get_type()) == 2);
  assert(get_type_dimension(el2->get_type()) == 2);

  for (int i = 0; i < el1->count_down(); ++i)
    for (int j = 0; j < el2->count_down(); ++j)
      if (do_edges_intersect(el1->get_down(i), el2->get_down(j)))
        return el2->get_down(j);

  return nullptr;
}

mesh::MeshEntityPtr AdjacencyPredicates::any_vertices_contained(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2)
{
  assert(get_type_dimension(el1->get_type()) == 2);
  assert(get_type_dimension(el2->get_type()) == 2);

  mesh::MeshEntityPtr verts2[mesh::MAX_DOWN];
  int nverts = get_downward(el2, 0, verts2);

  for (int i = 0; i < nverts; ++i)
  {
    auto r = m_classifier.classify(el1, el2, verts2[i]->get_point_orig(0));
    if (r.type != PointClassification::Exterior)
      return verts2[i];
  }

  return nullptr;
}

bool AdjacencyPredicates::is_point_contained(mesh::MeshEntityPtr el1, const utils::Point& pt)
{
  return m_classifier.classify_reverse(el1, pt, true).type != PointClassification::Exterior;
}

bool AdjacencyPredicates::is_point_contained(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2, const utils::Point& pt)
{
  return m_classifier.classify(el1, el2, pt).type != PointClassification::Exterior;
}

bool AdjacencyPredicates::do_edges_intersect(mesh::MeshEntityPtr edge1, mesh::MeshEntityPtr edge2)
{
  utils::Point edge1Pt1 = edge1->get_down(0)->get_point_orig(0);
  utils::Point edge1Pt2 = edge1->get_down(1)->get_point_orig(0);

  utils::Point edge2Pt1  = edge2->get_down(0)->get_point_orig(0);
  utils::Point edge2Pt2  = edge2->get_down(1)->get_point_orig(0);
  utils::Point pt1Normal = (*m_normalField)(edge2->get_down(0), 0, 0);
  utils::Point pt2Normal = (*m_normalField)(edge2->get_down(1), 0, 0);

  auto result = compute_edge_intersection(edge1Pt1, edge1Pt2, edge2Pt1, edge2Pt2, pt1Normal, pt2Normal, 1e-12, 0.25);

  bool firstIntersectionInRange = result.intersection1_found() && in_range(result.get_alpha1(), 0, 1, m_edgeTol) &&
                                  in_range(result.get_beta1(), 0, 1, m_edgeTol);

  bool secondIntersectionInRange = result.intersection2_found() && in_range(result.get_alpha2(), 0, 1, m_edgeTol) &&
                                   in_range(result.get_beta2(), 0, 1, m_edgeTol);

  return firstIntersectionInRange || secondIntersectionInRange;
}

bool AdjacencyPredicates::in_range(double val, double valMin, double valMax, double eps)
{
  return (val >= (valMin - eps)) && (val <= valMax + eps);
}

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
