#ifndef ADJACENCY_PREDICATES_H
#define ADJACENCY_PREDICATES_H

#include <stk_middle_mesh/field.hpp>
#include "intersection.hpp"
#include "point_classifier_normal_wrapper.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

class AdjacencyPredicates
{
  public:
    // For all API functions, el2/edge2 must be part of mesh2
    AdjacencyPredicates(std::shared_ptr<mesh::Mesh> mesh2)
      : m_classifier(mesh2)
      , m_normalField(m_classifier.get_normal_field())
    {}

    // checks if edge2 intersects with any edge of el1
    // Returns the first edge of el1 that intersects, otherwise nullptr
    mesh::MeshEntityPtr any_edge_intersect(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr edge2);

    // checks if any edge of el1 intersect any edge of el2
    // returns the edge of el2 that intersects el1 if such an edge exists,
    // nullptr otherwise
    mesh::MeshEntityPtr any_edges_intersect(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2);

    // Checks the edges of el2 to see if they intersect the edges of el1
    // If any edges intersect, returns the vertex inside el1, otherwise
    // nullptr
    // mesh::MeshEntityPtr getCrossingVertex(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2);

    // if any of el2's vertices are contained in el1, returns one of the
    // vertices, otherwise nullptr
    mesh::MeshEntityPtr any_vertices_contained(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2);

    // given a point on mesh1, returns true if it projects onto the given element
    // on mesh2
    bool is_point_contained(mesh::MeshEntityPtr el2, const utils::Point& pt);

    // given a point on mesh2 and the element it is part of (el2), returns true
    // if the point projects onto el1
    bool is_point_contained(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2, const utils::Point& pt);

  private:
    bool do_edges_intersect(mesh::MeshEntityPtr edge1, mesh::MeshEntityPtr edge2);

    bool in_range(double val, double valMin, double valMax, double eps);

    PointClassifierNormalWrapper m_classifier;
    mesh::FieldPtr<utils::Point> m_normalField;
    double m_edgeTol = 1e-6;
};

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
#endif
