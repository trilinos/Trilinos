#ifndef POINT_CLASSIFIER_H
#define POINT_CLASSIFIER_H

#include <stk_middle_mesh/element_operations_2d.hpp>
#include "intersection_common.hpp"
#include <stk_middle_mesh/mesh.hpp>
#include "point_classifier_impl.hpp"
#include "quad_to_triangles.hpp"

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
struct PossibleEdgeIntersection
{
    PointRecord record;
    double edgeXi = -1;
};

struct CornerRecord
{
    explicit CornerRecord(const bool hasIntersection_ = false, const PointRecord& r1_ = PointRecord(),
                          const PointRecord& r2_ = PointRecord(), const double xi1_ = -1, const double xi2_ = -1)
      : hasIntersection(hasIntersection_)
      , record1(r1_)
      , record2(r2_)
      , xi1(xi1_)
      , xi2(xi2_)
    {}
    bool hasIntersection;
    PointRecord record1;
    PointRecord record2;
    double xi1;
    double xi2;
};

class PointClassifier
{
  public:
    virtual ~PointClassifier(){};

    explicit PointClassifier(const double eps = 1e-13, const double smallSlopeTol = 1e-13);

    void set_projection(utils::impl::Projection* proj);

    // returns the point classification.  If classification is Vert or Edge,
    // the id is the local index of the entity of face
    PointRecord classify(mesh::MeshEntityPtr face, const utils::Point& pt);

    // For a point classified on an edge, computes the xi coordinate of the
    // point on that edge
    double get_edge_xi(const PointRecord& record);

    // given two points, one inside el and one outside the element, computes
    // intersection with boundary of el.  The PointRecord describes whether
    // the intersection is on a vert or an edge.  The edge_xi is gives
    // the xi coordinate of the intersection, if it is on an edge
    PossibleEdgeIntersection get_edge_xi(const PointRecord& record1, const PointRecord& record2);

    // Returns information on the possible overlap between an edge formed by
    // two points, both outside the domain.
    CornerRecord get_edge_xi_corner(const PointRecord& record1, const PointRecord& record2);

    // computes possible intersection on a line defined by a point classified
    // on either a vert or edge and an exterior point
    // PointRecord::id = -1 if no intersection
    // If there is an intersection, it must be classified on either a vertex
    // or edge.  If it is an edge, the xi coordinate of the intersection is also returned
    PossibleEdgeIntersection get_edge_intersection_xi(const PointRecord& record1, const PointRecord& record2);

    // for a point classified on either a vertex, edge, or interior of an element,
    // compute the XYZ coordinates of the intersection
    utils::Point compute_xyz_coords(const PointRecord& record1);

    // computes the length of the orthogonal projection line through record 1
    // onto the edge
    double compute_orthogonal_dist(const PointRecord& record1, const int id);

    // computes the xi coordinate on the specified edge of the point computed
    // by doing an orthogonal projection of record1 onto the edge
    double get_edge_xi_orthogonal(const PointRecord& record1, const int id);

    double get_eps() const { return m_classifier.get_eps(); }

    void set_eps(const double val) { m_classifier.set_eps(val); }

  protected:
    // given the xi coordinate on an edge of m_el1 or m_el2, computes
    // the xi coordinate on the edge of r.el (in the coordinate system
    // of that edge, not the reference element)
    // el is the quad element, and id is the edge id of the quad that
    // tri_xi lies on
    double convert_edge_xi(double triXi, mesh::MeshEntityPtr el, const int id);

    // given the result of a call to PointClassifier::getEdgeXi for a given
    // triangle, tries to classify points on the exterior edges if
    // at all possible (ie. if the intersection is on the vertex of the
    // interior edge, classify the intersection on a non-interior edge
    // std::pair<int, double> normalizeEdgeIntersect(std::pair<int, double> p, const int tri_id);

    //-------------------------------------------------------------------------
    // Helper functions for getEdgeXiCorner

    PointRecord convert_point_record(const PointRecordForTriangle& record);

    CornerRecord convert_tri_corners_to_quad_corner(mesh::MeshEntityPtr quadEl, const CornerRecordForTriangle& c1,
                                                    const CornerRecordForTriangle& c2);

    PossibleEdgeIntersection get_exterior_quad_point_record(mesh::MeshEntityPtr quadEl,
                                                            const CornerRecordForTriangle& corner);

    PossibleEdgeIntersection get_exterior_quad_point_record(mesh::MeshEntityPtr quadEl,
                                                            const PointRecordForTriangle& record, double xi);

    //-------------------------------------------------------------------------
    // helper functions for getEdgeIntersectXi

    // excludes some vertex intersection matches
    // p is the result of a call to impl::PointClassifier::getEdgeIntersectXi
    // with include_edge_verts = false
    // r: is the PointRecord of the point
    // r_inner: is the PointRecordForTriangle for the triangle element that
    // was passed into getEdgeIntersectXi
    // std::pair<impl::PointRecord, double> excludeVertsAdjacentInQuad(std::pair<impl::PointRecord, double> p, const
    // PointRecord& r, const PointRecordForTriangle r_inner);

    // creates the two PointRecordForTriangles for a given vert of tri_el
    std::pair<PointRecordForTriangle, PointRecordForTriangle>
    make_triangle_records_for_quad_vert_intersection(mesh::MeshEntityPtr triEl, const int vertId);

    // creates the two PointRecordForTriangles for the point on the given edge of
    // tri_el
    std::pair<PointRecordForTriangle, PointRecordForTriangle>
    make_triangle_records_for_quad_edge_intersection(mesh::MeshEntityPtr triEl, const int edgeId, const double xi);

    // converts the output of impl::PointClassifier::getEdgeIntersectionXi
    // to the output of this class
    PossibleEdgeIntersection
    make_triangle_records_for_intersection_that_must_exist(mesh::MeshEntityPtr quadEl, mesh::MeshEntityPtr triEl,
                                                           PossibleEdgeIntersectionForTriangle p);

    PointClassifierForTriangle m_classifier;
    QuadToTriangles m_quadToTriangles;
    mesh::impl::ElementOperations2D m_elemOps;
    utils::impl::Projection m_defaultProjection;
};

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
#endif
