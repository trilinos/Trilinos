#ifndef MIDDLE_GRID_CONSTRAINT_GENERATOR_H
#define MIDDLE_GRID_CONSTRAINT_GENERATOR_H

#include "mesh_relational_data.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"

#include <set> //TODO: DEBUGGING

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

class MiddleGridConstraintGenerator
{
    std::set<int> m_vertIds = {2486, 2810}; // TODO: DEBUGGING

  public:
    MiddleGridConstraintGenerator(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                                  std::shared_ptr<mesh::Mesh> meshIn,
                                  std::shared_ptr<MeshRelationalData> relationalData,
                                  std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> pointClassifier)
      : m_mesh1(mesh1)
      , m_mesh2(mesh2)
      , m_meshIn(meshIn)
      , m_relationalData(relationalData)
      , m_pointClassifier(pointClassifier)
    {}

    void generate();

  private:
    using ExchangerKnown = stk::DataExchangeKnownPatternNonBlockingBuffer<int>;

    void create_mesh1_vertices();

    void set_mesh1_vert_shared_entities(ExchangerKnown& exchanger);

    void create_mesh2_interior_vertices();

    void create_mesh1_edges();

    void set_mesh1_edge_shared_entities(ExchangerKnown& exchanger);

    void split_edges();

    void pack_edge_split_shared_info(ExchangerKnown& exchanger,
          mesh::VariableSizeFieldPtr<int> sharedEntityInfoPtr, mesh::MeshEntityPtr edge1);

    void set_edge_split_shared_entities(ExchangerKnown& exchanger, 
           mesh::VariableSizeFieldPtr<int> sharedEntityInfoPtr);

    void sort_edge_splits(mesh::MeshEntityPtr edge1);

    mesh::MeshEntityPtr get_mesh_in_edge_from_mesh1_edge(mesh::MeshEntityPtr edge1);

    void create_internal_edges();

    void sort_verts_on_edge2(mesh::MeshEntityPtr edge2);

    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    std::shared_ptr<mesh::Mesh> m_meshIn;
    std::shared_ptr<MeshRelationalData> m_relationalData;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifier;   
    const bool m_output = false;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
