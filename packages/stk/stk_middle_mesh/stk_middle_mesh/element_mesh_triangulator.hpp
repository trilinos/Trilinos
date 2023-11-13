#ifndef ELEMENT_MESH_TRIANGULATOR_H
#define ELEMENT_MESH_TRIANGULATOR_H

#include "element_mesh_extractor.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

class ElementMeshTriangulator
{
  public:
    ElementMeshTriangulator(std::shared_ptr<MeshRelationalData> relationalData)
      : m_relationalData(relationalData)
    {}

    void triangulate(ElementMeshData elementMeshData);

  private:
    void perturb_boundary_nodes(double fac);

    void triangulate_element_mesh();

    void fix_triangulation(int numConstraintEdges);

    bool should_delete_edges(const predicates::impl::PointRecord& class1, const predicates::impl::PointRecord& class2,
                             mesh::MeshEntityPtr v1ElementMeshIn, mesh::MeshEntityPtr v2ElementMeshIn);

    bool should_delete_edge_between_vert_and_edge(const predicates::impl::PointRecord& classV,
                                                  const predicates::impl::PointRecord& classE,
                                                  mesh::MeshEntityPtr elementMeshVertOnMesh1Edge);

    void fix_triangulation_with_boundary_ordering();

    // TODO: DEBUGGING
    void check_for_cap_elements();

    int get_edge_split_index(mesh::MeshEntityPtr edge1, mesh::MeshEntityPtr elementMeshVert);

    std::shared_ptr<MeshRelationalData> m_relationalData;
    ElementMeshData m_elementMeshData;
    static constexpr bool M_OUTPUT = false;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
