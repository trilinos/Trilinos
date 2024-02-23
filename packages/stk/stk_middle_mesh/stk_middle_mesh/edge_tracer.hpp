#ifndef EDGE_TRACER_H
#define EDGE_TRACER_H

// #include "predicates/point_classifier_base.hpp"
#include "edge_tracer_opts.hpp"
#include "field.hpp"
#include "mesh.hpp"
#include "predicates/intersection_common.hpp"
#include "projection.hpp"
#include <iterator>
#include <limits>
#include <vector>

namespace stk {
namespace middle_mesh {

namespace impl::testing {
class EdgeTracerTesterBase;
}

namespace nonconformal4 {
namespace impl {

struct EdgeIntersection
{
    predicates::impl::PointRecord record;
    double alpha; // distance along the input edge, in range [0, 1]
    double beta;  // distance along the edge found to intersect the input edge, in range [0, 1]
                  // value is only defined if record.type == PointClassification::Edge, except
                  // for the endpoint of the line, when this value is undefined
};

inline std::ostream& operator<<(std::ostream& os, const EdgeIntersection& intersection)
{
  os << intersection.record << ", alpha = " << intersection.alpha << ", beta = " << intersection.beta;
  return os;
}

class EdgeTracer
{
  public:
    // the edge_error_tol is used to determine which of the two possible edge intersection solutions is a solution to
    // the original problem.  Emperically, it can be fairly large because the incorrect solution has a significant
    // error.
    EdgeTracer(std::shared_ptr<mesh::Mesh> mesh, mesh::FieldPtr<utils::Point> normalField,
               ::stk::middle_mesh::impl::EdgeTracerTolerances tolerances)
      : m_mesh(mesh)
      , m_normalField(normalField)
      , m_tolerances(tolerances)
    {}

    // convenience method for tracing an edge.  edge must be part of the
    // mesh passed into the constructor.
    void trace_edge(mesh::MeshEntityPtr edge, const predicates::impl::PointRecord& recordStart,
                    const predicates::impl::PointRecord& recordEnd, std::vector<EdgeIntersection>& intersections);

    // overwrites intersections vector with all the intersections of the
    // given edge with edges of the other mesh.
    // Note that this does not include the first or last vertex of the
    // edge, even if they happen to be on an edge of mesh2.
    // pt_start: xyz coordinates of first vert of the edge
    // normal_start: normal vector at first vert of the edge
    // record_start: the PointRecord describing the projection of pt_start onto
    //               the other mesh (not the mesh passed into the constructor)
    // pt_end: xyz coordinates of the second vert of the edge
    // normal_end: normal vector at second vert of the edge
    // record_end: the PointRecord describing the projection of pt_end onto
    //             the other mesh (not the mesh passed into the constructor)
    void trace_edge(const utils::Point& ptStart, const utils::Point& normalStart,
                    const predicates::impl::PointRecord& recordStart, const utils::Point& ptEnd,
                    const utils::Point& normalEnd, const predicates::impl::PointRecord& recordEnd,
                    std::vector<EdgeIntersection>& intersections);

  private:
    // Given any of:
    //  1) the starting point of the edge
    //  2) the starting point of the edge and the first intersection of the edge with another
    //     mesh entity
    //  3) the two most recent intersections of the edge with other mesh entities
    // computes the next intersection of the edge when traversing the edge from start to end.
    // When the next intersection is the endpoint of the line, the return value is undefined.
    // The at_end_of_line variable determines when the end of the line is reached
    // prev_intersection: the previous intersection point, or a default constructed PointRecord if
    //                    this is the first intersection computed for the edge
    // current_intersection: the most recent intersection point (or the starting point of the edge if
    //                       this is the first intersection computed for the edge)
    // end_vert: the PointRecord for the endpoint of the edge
    // start_pt: the original coordinates of the beginning of the edge
    // end_pt: the original coordinates of the end of the edge
    // TODO: the use of PointRecord is unfortunate because it is a much
    //      larger struct than needed
    EdgeIntersection compute_next_intersection(const predicates::impl::PointRecord& prevIntersection,
                                               const predicates::impl::PointRecord& currentIntersection,
                                               const predicates::impl::PointRecord& endVert,
                                               const utils::Point& ptStart, const utils::Point& ptEnd,
                                               const utils::Point& normalStart, const utils::Point& normalEnd,
                                               bool& atEndOfLine);

    std::vector<mesh::MeshEntityPtr> get_included_elements(mesh::MeshEntityPtr entity,
                                                           const predicates::impl::PointRecord& prevIntersection,
                                                           const predicates::impl::PointRecord& currentIntersection);

    std::vector<mesh::MeshEntityPtr> get_included_edges(const std::vector<mesh::MeshEntityPtr>& includedElements,
                                                        const predicates::impl::PointRecord& currentIntersection);

    std::vector<mesh::MeshEntityPtr> get_elements(mesh::MeshEntityPtr entity);

    std::vector<mesh::MeshEntityPtr> get_excluded_elements(const predicates::impl::PointRecord& prevIntersection,
                                                           const predicates::impl::PointRecord& currentIntersection);

    std::vector<mesh::MeshEntityPtr> get_common_elements(mesh::MeshEntityPtr entity1, mesh::MeshEntityPtr entity2);

    std::vector<mesh::MeshEntityPtr> get_common_entities(std::vector<mesh::MeshEntityPtr>& entities1,
                                                         std::vector<mesh::MeshEntityPtr>& entities2);

    std::vector<mesh::MeshEntityPtr> get_included_entities(std::vector<mesh::MeshEntityPtr>& allEntities,
                                                           const std::vector<mesh::MeshEntityPtr>& excludedEntities);

    bool at_end_of_line(const std::vector<mesh::MeshEntityPtr>& includedElements,
                        const predicates::impl::PointRecord& endVert);

    bool contains_entity(const std::vector<mesh::MeshEntityPtr>& elements, mesh::MeshEntityPtr entity);

    std::vector<mesh::MeshEntityPtr> get_edges(const std::vector<mesh::MeshEntityPtr>& elements);

    std::vector<mesh::MeshEntityPtr> get_excluded_edges(const predicates::impl::PointRecord& record);

    int compute_best_intersection(const utils::Point& startPt, const utils::Point& endPt,
                                  const utils::Point& startNormal, const utils::Point& endNormal,
                                  const std::vector<mesh::MeshEntityPtr>& edges, double& alpha, double& beta);

    int choose_best_intersection(const std::vector<double>& alphas, const std::vector<double>& betas,
                                 const std::vector<mesh::MeshEntityPtr>& edges,
                                 const std::vector<int>& intersectionIdxToEdgeIdx, const utils::Point& edgeStartPt,
                                 const utils::Point& edgeEndPt);

    // int chooseBestIntersectionWhenZeroInRange(const std::vector<double>& alphas,
    //                                           const std::vector<double>& betas);

    // returns a value describing how far outside the range [val_min, val_max] the
    // given value is.  The value is always positive
    double compute_deviation(double val, double valMin, double valMax);

    // when there are multiple possible intersections (as defined by intersections_in_range),
    // choose the one that gives a line with the smallest angle between the edge
    // (edge_start_point to edge_end_point) and the edge that would be formed
    // by the intersection (edge_start_point to the point defined by beta)
    int choose_intersection_with_minimum_angle(const std::vector<double>& alphas, const std::vector<double>& betas,
                                               const std::vector<int>& intersectionsInRange,
                                               const std::vector<mesh::MeshEntityPtr>& edges,
                                               const std::vector<int>& intersectionIdxToEdgeIdx,
                                               const utils::Point& edgeStartPt, const utils::Point& edgeEndPt);

    // int chooseBestIntersectionWhenMultipleInRange(const std::vector<double>& alphas,
    //                                               const std::vector<double>& betas,
    //                                               const std::vector<int>& intersections_in_range);

    EdgeIntersection create_edge_intersection(mesh::MeshEntityPtr edge,
                                              const std::vector<mesh::MeshEntityPtr>& includedElements, double alpha,
                                              double beta);

    mesh::MeshEntityPtr get_intersected_entity(mesh::MeshEntityPtr edge, double beta);

    double clamp(double val, double valMin, double valMax);

    mesh::MeshEntityPtr get_parent_element(mesh::MeshEntityPtr edge,
                                           const std::vector<mesh::MeshEntityPtr>& includedElements);

    bool is_entity_on_mesh(mesh::MeshEntityPtr entity);

    void redistribute_points_near_endpoints(std::vector<EdgeIntersection>& intersections, double eps);

    void redistribute_points_near_endpoint(std::vector<EdgeIntersection>& intersections, bool isLeftEndpoint,
                                           double eps);

    std::shared_ptr<mesh::Mesh> m_mesh;
    mesh::FieldPtr<utils::Point> m_normalField;
    // EdgeIntersectionCalculator m_edge_intersector;
    ::stk::middle_mesh::impl::EdgeTracerTolerances m_tolerances;

    friend class stk::middle_mesh::impl::testing::EdgeTracerTesterBase;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
