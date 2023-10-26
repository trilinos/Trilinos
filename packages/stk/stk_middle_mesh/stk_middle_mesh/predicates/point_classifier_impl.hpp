#ifndef POINT_CLASSIFIER_IMPL_H
#define POINT_CLASSIFIER_IMPL_H

#include <stk_middle_mesh/element_operations_2d.hpp>
#include "intersection_common.hpp"
#include <stk_middle_mesh/mesh.hpp>
#include <stk_middle_mesh/predicates/element_xi_classifier.hpp>
#include "triangle_coord_utils.hpp"
#include <array>
#include <ostream>

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

// used for returning data on a possible intersection with an edge.
// If there is no intersection, record.id = -1.  If the intersection
// is on a vertex, record.type is PointClassification::Vert and the
// record.id identifies which vert, and the edge_xi value is undefined.
// If the intersection is on the interior of an edge, edge_xi gives
// xi coordinate of the intersection on the edge.
struct PossibleEdgeIntersectionForTriangle
{
    PointRecordForTriangle record;
    double edgeXi;
};

// Returns information on the possible overlap between an edge formed by
// two points, both outside the domain.
// has_intersection tells if there is any intersection.  If there is, the
// other fields identify which edges and the xi coordinates on those edges
// of the intersection
struct CornerRecordForTriangle
{
    explicit CornerRecordForTriangle(const bool hasIntersection_       = false,
                                     const PointRecordForTriangle& r1_ = PointRecordForTriangle(),
                                     const PointRecordForTriangle& r2_ = PointRecordForTriangle(),
                                     const double xi1_ = -1, const double xi2_ = -1)
      : hasIntersection(hasIntersection_)
      , record1(r1_)
      , record2(r2_)
      , xi1(xi1_)
      , xi2(xi2_)
    {}
    bool hasIntersection;
    PointRecordForTriangle record1;
    PointRecordForTriangle record2;
    double xi1;
    double xi2;
};

inline std::ostream& operator<<(std::ostream& os, const CornerRecordForTriangle& record)
{
  os << "Corner Record: has intersection = ";
  if (record.hasIntersection)
    os << "true";
  else
    os << "false";

  if (record.hasIntersection)
  {
    os << ", r1 = " << record.record1 << ", record = " << record.record2;
    os << ", xi1 = " << record.xi1 << ", xi2 = " << record.xi2;
  }

  return os;
}

class PointClassifierForTriangle
{
  public:
    explicit PointClassifierForTriangle(const double eps = 1e-13, const double smallSlopeTol = 1e-13)
      : m_eps(eps)
      , m_smallSlopeTol(smallSlopeTol)
      , m_xiClassifier(eps)
    {}

    void set_projection(utils::impl::Projection* proj);

    // returns the point classification.  If classification is Vert or Edge,
    // the id is the local index of the entity of face
    PointRecordForTriangle classify(mesh::MeshEntityPtr face, const utils::Point& pt);

    // For a point classified on an edge, computes the xi coordinate of the
    // point on that edge
    double get_edge_xi(const PointRecordForTriangle& record);

    // given two points, one inside el and one outside the element, computes
    // intersection with boundary of el.  The PointRecord describes whether
    // the intersection is on a vert or an edge, and the second part of the
    // pair is the xi coordinate of the intersection, if it is on an edge
    PossibleEdgeIntersectionForTriangle get_edge_xi(const PointRecordForTriangle& record1,
                                                    const PointRecordForTriangle& record2);

    CornerRecordForTriangle get_edge_xi_corner(const PointRecordForTriangle& record1,
                                               const PointRecordForTriangle& record2);

    // computes possible intersection on a line defined by a point classified
    // on either a vert or edge and an exterior point.
    // PointRecord::id = -1 if no intersection.
    // If there is an intersection, it must be classified on either a vertex
    // or edge.  If it is an edge, the second element of the std::pair
    // gives the xi coordinate of the intersection on the edge.
    // exclude_edge_verts: if true (default), if one of the points is
    // classified on a vert, both that vert and the verts on the other
    // end of the adjacent edges are excluded from possible intersections.
    // If false, only that vert is excluded, and not the edge-adjacent
    // verts.
    // If the point is classified on an edge, exclude_edge_verts excludes
    // the verts on the ends of the edge
    PossibleEdgeIntersectionForTriangle get_edge_intersection_xi(const PointRecordForTriangle& record1,
                                                                 const PointRecordForTriangle& record2,
                                                                 const bool excludeEdgeVerts = true);

    utils::Point compute_xyz_coords(const PointRecordForTriangle& record);

    // computes the length of the orthogonal projection line through record 1
    // onto the edge
    double compute_orthogonal_dist(const PointRecordForTriangle& record1, const int id);

    // computes the xi coordinate on the specified edge of the point computed
    // by doing an orthogonal projection of record1 onto the edge
    double get_edge_xi_orthogonal(const PointRecordForTriangle& record1, const int id);

    double get_eps() const { return m_eps; }

    void set_eps(const double val) { m_eps = val; }

  private:
    // handles a line between two exterior point that is not degenerate
    // (horizontal or vertical)
    void compute_intersection_general_corner(mesh::MeshEntityPtr el, const PointRecordForTriangle& record1,
                                             const PointRecordForTriangle& record2, CornerRecordForTriangle& recordOut);

    // handles a line that can only possibly intersect with 2 edges (ie.
    // a horizontal or vertical line)
    void compute_intersection_degenerate_line(mesh::MeshEntityPtr el, const PointRecordForTriangle& record1,
                                              const PointRecordForTriangle& record2, const int id1, const int id2,
                                              CornerRecordForTriangle& recordOut);

    void set_corner_record(mesh::MeshEntityPtr el, int id1, int id2, double xi1, double xi2,
                           CornerRecordForTriangle& recordOut);

    // computes the xi coordinates for the given point, which lies on the
    // specified edge
    double get_edge_xi(mesh::MeshEntityType type, const int id, const utils::Point& ptXi,
                       mesh::EntityOrientation orient);

    // helper function for public getEdgeXi, computes intersection
    // of line with edge
    std::pair<int, double> get_edge_xi_triangle(mesh::MeshEntityPtr el, const PointRecordForTriangle& record1,
                                                const PointRecordForTriangle& record2);

    // helper function for public getEdgeXi, computes intersection of line with
    // edge
    std::pair<int, double> get_edge_xi_quad(mesh::MeshEntityPtr el, const PointRecordForTriangle& record1,
                                            const PointRecordForTriangle& record2);

    // takes an (edge_id, edge_xi) pair and determines if it is on a vertex or edge interior
    // Returns a PointRecord and a xi value.  The xi value is only meanful if the intersection
    // is on the edge interior
    PossibleEdgeIntersectionForTriangle
    classify_intersection_on_vert_if_possible(mesh::MeshEntityPtr el, const std::pair<int, double>& edgeIntersection);

    // computes a point on line given by alpha * (record2.pt_xi - record1.pt_xi) + record1.pt_xi
    utils::Point compute_point_on_line(const PointRecordForTriangle& record1, const PointRecordForTriangle& record2,
                                       const double alpha);

    // uses the xi coordinates in record1 and record2 to define a line, and computes
    // the intersection of that line with the specified edge of the element.  Any point
    // on the line can be expressed as p = alpha * (record2.pt_xi - record1.pt_xi) + record1.pt_xi.
    // This function returns the alpha value at the intersection with the specified edge
    double compute_edge_intersection(mesh::MeshEntityType type, const PointRecordForTriangle& record1,
                                     const PointRecordForTriangle& record2, const int id);

    // compute xi coordinate of line that passes through specified edge
    double compute_edge_intersection_quad(const PointRecordForTriangle& record1, const PointRecordForTriangle& record2,
                                          const int id);

    double compute_edge_intersection_triangle(const PointRecordForTriangle& record1,
                                              const PointRecordForTriangle& record2, const int id);

    // for case where both points are close to an edge
    std::pair<int, double> compute_xi_closest(mesh::MeshEntityPtr el, const PointRecordForTriangle& record1,
                                              const PointRecordForTriangle& record2);

    double clamp(const double val, const double vL, const double vR) { return std::min(std::max(val, vL), vR); }

    // excluded_verts and excluded_edges are masks: true if the ith entity
    // should be excluded
    PointRecordForTriangle get_edge_intersection_xi(mesh::MeshEntityPtr el1, const PointRecordForTriangle& recordI,
                                                    const PointRecordForTriangle& recordE,
                                                    const bool excludedVerts[mesh::MAX_DOWN],
                                                    const bool excludedEdges[mesh::MAX_DOWN]);

    double m_eps;
    double m_smallSlopeTol = 1e-8;
    ElementXiClassifier m_xiClassifier;
    TriangleCoordUtils m_triangleCoordUtils;
    mesh::impl::ElementOperations2D m_elemOps;
};

PointRecordForTriangle classify_triangle(const utils::Point& ptXi);

PointRecordForTriangle classify_quad(const utils::Point& ptXi);

} // namespace impl
} // namespace predicates
} // namespace middle_mesh
} // namespace stk
#endif
