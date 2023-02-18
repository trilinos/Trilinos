#ifndef POINT_CLASSIFIER_NORMAL_WRAPPER
#define POINT_CLASSIFIER_NORMAL_WRAPPER

#include "field.h"
#include "mesh.h"
#include "point_classifier_normal_interpolation.h"
#include "predicates/intersection_common.h"
#include "predicates/quad_to_triangles.h"
#include "triangle_coord_utils.h"

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
      , m_tolerances(tolerances)
    {
      compute_averaged_normal_vectors(mesh);
    }

    // src_face must be a face on the mesh that was passed to the constructor
    PointRecord classify(mesh::MeshEntityPtr destFace, mesh::MeshEntityPtr srcFace, const utils::Point& pt);

    // project point pt on to face dest_face, which must be part of the mesh
    // that was passed into the constructor
    PointRecord classify_reverse(mesh::MeshEntityPtr destFace, const utils::Point& pt, bool logicalResultOnly = false);

    // for a point classified on either a vertex, edge, or interior of an element,
    // compute the XYZ coordinates of the intersection
    utils::Point compute_xyz_coords(const PointRecord& record1);

    // For a point classified on an edge, computes the xi coordinate of the
    // point on that edge
    double get_edge_xi(const PointRecord& record);

    // computes the length of the orthogonal projection line through record 1
    // onto the edge
    double compute_orthogonal_dist(const PointRecord& record1, const int id);

    mesh::FieldPtr<utils::Point> get_normal_field() { return m_normalField; }

  private:
    std::array<utils::Point, 4> get_normal_vectors(mesh::MeshEntityPtr face);

    utils::Point compute_normal_vector(mesh::MeshEntityPtr face, const utils::Point& pt);

    utils::Point compute_normal_vector_quad(mesh::MeshEntityPtr face, const utils::Point& pt);

    utils::Point compute_normal_vector_triangle(mesh::MeshEntityPtr face, const utils::Point& pt);

    utils::Point compute_triangle_xi_coords(std::array<mesh::MeshEntityPtr, 3>& verts, const utils::Point& pt);

    utils::Point compute_triangle_xi_coords(std::array<utils::Point, 3>& verts, const utils::Point& pt);

    void compute_averaged_normal_vectors(std::shared_ptr<mesh::Mesh> mesh);

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
