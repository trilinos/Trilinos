#ifndef ELEMENT_MESH_EXTRACTOR_H
#define ELEMENT_MESH_EXTRACTOR_H

#include "field.hpp"
#include "mesh.hpp"
#include "mesh_relational_data.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

namespace {
struct EdgeVertexIdPair
{
    EdgeVertexIdPair(int id1_, int id2_)
      : id1(std::min(id1_, id2_))
      , id2(std::max(id1_, id2_))
    {}

    int id1;
    int id2;
};

inline std::ostream& operator<<(std::ostream& os, const EdgeVertexIdPair& ids)
{
  os << ids.id1 << ", " << ids.id2;
  return os;
}

inline bool operator<(const EdgeVertexIdPair& lhs, const EdgeVertexIdPair& rhs)
{
  if (lhs.id1 != rhs.id1)
    return lhs.id1 < rhs.id1;
  else
    return lhs.id2 < rhs.id2;
}

inline bool operator==(const EdgeVertexIdPair& lhs, const EdgeVertexIdPair& rhs)
{
  return lhs.id1 == rhs.id1 && lhs.id2 == rhs.id2;
}

inline bool operator!=(const EdgeVertexIdPair& lhs, const EdgeVertexIdPair& rhs)
{
  return !(lhs == rhs);
}

} // namespace

struct ElementMeshData
{
    mesh::MeshEntityPtr el1 = nullptr;
    std::shared_ptr<mesh::Mesh> elementMeshIn;
    mesh::FieldPtr<predicates::impl::PointRecord> elementMeshInVertClass;
    mesh::FieldPtr<mesh::MeshEntityPtr> elementMeshEntitiesToMeshInEntities;
};

// given mesh_in (created from the intersection of 2 meshes), it creates a mesh with
// the entities that are on a given element of mesh1
class ElementMeshExtractor
{
  public:
    ElementMeshExtractor(std::shared_ptr<mesh::Mesh> meshIn, std::shared_ptr<MeshRelationalData> relationalData)
      : m_meshIn(meshIn)
      , m_relationalData(relationalData)
    {}

    ElementMeshData extract_element_mesh(mesh::MeshEntityPtr el1);

    void write_elements_back_to_middle_grid(ElementMeshData& elementMeshData, mesh::MeshEntityPtr el1);

  private:
    void get_el1_edges_in(mesh::MeshEntityPtr el1, std::vector<EdgeVertexIdPair>& edgesIn);

    void sort_mesh_in_verts(mesh::MeshEntityPtr el1);

    int get_vert_in_mesh1_el_index(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr vertIn);

    predicates::impl::PointRecord convert_classification_to_element(mesh::MeshEntityPtr el1,
                                                                    const predicates::impl::PointRecord& record);

    std::shared_ptr<mesh::Mesh> m_meshIn;
    std::shared_ptr<MeshRelationalData> m_relationalData;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
