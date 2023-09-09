#ifndef POINT_CLASSIFIER_NORMAL_WRAPPER
#define POINT_CLASSIFIER_NORMAL_WRAPPER

#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/mesh_entity.hpp"
#include "point_classifier_normal_interpolation.hpp"
#include "intersection_common.hpp"
#include "quad_to_triangles.hpp"
#include "triangle_coord_utils.hpp"
#include "average_normal_field.hpp"

namespace stk {
namespace middle_mesh {

namespace impl::testing {
class PointClassifierNormalWrapperTester;
}

namespace predicates {
namespace impl {

using stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances;

class PointClassifierNormalWrapper
{
  public:
    explicit PointClassifierNormalWrapper(
        std::shared_ptr<mesh::Mesh> mesh,
        const PointClassifierNormalWrapperTolerances& tolerances = PointClassifierNormalWrapperTolerances())
      : m_normalClassifier(tolerances.normalInterpolationTolerances)
      , m_mesh(mesh)
      , m_normalField(AveragedNormalField(mesh).get_field())
      , m_tolerances(tolerances)
    {}

    QuadToTriangles& get_quad_to_triangles() { return m_quadToTriangles; }

    // src_face must be a face on the mesh that was passed to the constructor
    PointRecord classify(mesh::MeshEntityPtr destFace, mesh::MeshEntityPtr srcFace, const utils::Point& pt);

    // project point pt on to face dest_face, which must be part of the mesh
    // that was passed into the constructor
    PointRecord classify_reverse(mesh::MeshEntityPtr destFace, const utils::Point& pt, bool logicalResultOnly = false);

    // for a point classified on either a vertex, edge, or interior of an element,
    // compute the XYZ coordinates of the intersection.
    // Points classifies on the exterior are allowed if allowExterior=true
    utils::Point compute_xyz_coords(const PointRecord& record1, bool allowExterior=false);

    // compute the parametric coordinates of the point on record1.el.  Point
    // must be classifier on a vertex, edge, or interior of the element.  Can
    // be on the exterior if allowExterior=true
    utils::Point compute_xi_coords(const PointRecord& record1, bool allowExterior=false);

    // For a point classified on an edge, computes the xi coordinate of the
    // point on that edge
    double get_edge_xi(const PointRecord& record);

    // computes the length of the orthogonal projection line through record 1
    // onto the edge
    double compute_orthogonal_dist(const PointRecord& record1, const int id);

    mesh::FieldPtr<utils::Point> get_normal_field() { return m_normalField; }

    PointRecord create_vert_record(mesh::MeshEntityPtr el, int vertId);

    PointRecord create_edge_record(mesh::MeshEntityPtr el, int edgeId, double edgeXiOnReferenceEl);

    PointRecord classify_onto(const PointRecord& record, mesh::MeshEntityPtr el);

  private:
    std::array<utils::Point, 4> get_normal_vectors(mesh::MeshEntityPtr face);

    utils::Point compute_normal_vector(mesh::MeshEntityPtr face, const utils::Point& pt);

    utils::Point compute_normal_vector_quad(mesh::MeshEntityPtr face, const utils::Point& pt);

    utils::Point compute_normal_vector_triangle(mesh::MeshEntityPtr face, const utils::Point& pt);

    utils::Point compute_triangle_xi_coords(std::array<mesh::MeshEntityPtr, 3>& verts, const utils::Point& pt);

    utils::Point compute_triangle_xi_coords(std::array<utils::Point, 3>& verts, const utils::Point& pt);

    bool is_entity_on_mesh(mesh::MeshEntityPtr entity);

    PointClassifierNormalInterpolation m_normalClassifier;
    std::shared_ptr<mesh::Mesh> m_mesh;
    mesh::FieldPtr<utils::Point> m_normalField;
    TriangleCoordUtils m_triangleCoordUtils;
    PointClassifierNormalWrapperTolerances m_tolerances;
    QuadToTriangles m_quadToTriangles;

    friend middle_mesh::impl::testing::PointClassifierNormalWrapperTester;
};

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
#endif
