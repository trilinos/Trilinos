#ifndef MESH_PROJECTION_CALCULATOR
#define MESH_PROJECTION_CALCULATOR

#include "edge_tracer.hpp"
#include "mesh_relational_data.hpp"
#include "predicates/point_classifier_normal_wrapper.hpp"
#include "bounding_box_search_opts.hpp"
#include <set> //TODO: DEBUGGING
#include <unordered_set>

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

class MeshProjectionCalculator
{
    std::set<int> m_vertIds = {1067, 2900};

  public:
    MeshProjectionCalculator(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                             std::shared_ptr<MeshRelationalData> relationalData,
                             std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> pointClassifier,
                             const middle_mesh::impl::EdgeTracerTolerances& edgeTracerTolerances,
                             [[maybe_unused]] const search::BoundingBoxSearchOpts& searchOpts = search::BoundingBoxSearchOpts())
      : m_mesh1(mesh1)
      , m_mesh2(mesh2)
      , m_pointClassifier(pointClassifier)
      , m_edgeTracer(m_mesh2, m_pointClassifier->get_normal_field(), edgeTracerTolerances)
      , m_relationalData(relationalData)
    {}

    void project();

  private:
    template <typename T>
    using SetType = std::set<T, mesh::MeshEntityCompare>;

    using Bool = int_least8_t;

    void create_mesh1_fake_verts();

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh2_to_mesh1_element_maps();

    void project_mesh2_vertices_onto_mesh1();

    void choose_unique_vert_projection(
        mesh::MeshEntityPtr vert2, const std::vector<predicates::impl::PointRecord>& verts2AllClassificationsPtr);

    void record_vert2_classification(mesh::MeshEntityPtr vert2, const predicates::impl::PointRecord& record,
                                     const utils::Point& pt);

    void record_edge2_vertex_association(mesh::VariableSizeField<VertOnEdge>& edges2ToFakeVertsIn,
                                         mesh::MeshEntityPtr vert2, FakeVert fv);

    int get_closest_projection(const std::vector<predicates::impl::PointRecord>& verts2AllClassificationsPtr,
                               mesh::MeshEntityPtr vert2, utils::Point& closestPt);

    void intersect_mesh2_edges_with_mesh1();

    void record_mesh2_edge_on_mesh1_vert(mesh::MeshEntityPtr edge2, const EdgeIntersection& intersection);

    void record_mesh2_edge_intersects_mesh1_edge(mesh::MeshEntityPtr edge2, const EdgeIntersection& intersection);

    void sort_edge2_fake_verts(mesh::MeshEntityPtr edge2);

    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifier;
    EdgeTracer m_edgeTracer;
    FakeVertGenerator m_fakeVertGen;
    std::shared_ptr<MeshRelationalData> m_relationalData;
    search::BoundingBoxSearchOpts m_searchOpts;

    const bool m_output = false;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
