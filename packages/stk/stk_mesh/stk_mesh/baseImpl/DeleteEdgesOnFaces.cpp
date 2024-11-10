#include "stk_mesh/baseImpl/DeleteEdgesOnFaces.hpp"

#include <vector>
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace mesh {
namespace impl {

EdgesOnFacesDeleter::EdgesOnFacesDeleter(BulkData& bulk, Part* edgePart) :
  m_bulk(bulk),
  m_edgePart(edgePart)
{
  STK_ThrowRequire(edgePart->topology().rank() == stk::topology::EDGE_RANK);
}

void EdgesOnFacesDeleter::delete_surface_edges()
{
  stk::mesh::Selector edgesSelector = *m_edgePart;
  std::vector<stk::mesh::Entity> edges;
  stk::mesh::get_entities(m_bulk, stk::topology::EDGE_RANK, edgesSelector, edges);

  std::vector<stk::mesh::Entity> nodes, faces;
  std::vector<int> ordinals;
  for (stk::mesh::Entity edge : edges)
  {
    nodes.assign(m_bulk.begin_nodes(edge), m_bulk.end_nodes(edge));
    ordinals.assign(m_bulk.begin_node_ordinals(edge), m_bulk.end_node_ordinals(edge));
    delete_edge_to_node_relations(edge, nodes, ordinals);

    faces.assign(m_bulk.begin_faces(edge), m_bulk.end_faces(edge));
    delete_face_to_edge_relations(faces, edge);

    bool didDelete = m_bulk.destroy_entity(edge);
    STK_ThrowRequireMsg(didDelete, "failed to destroy entity");
  }
}

void EdgesOnFacesDeleter::delete_edge_to_node_relations(Entity edge, const std::vector<Entity>& nodes, const std::vector<int>& ordinals)
{
  for (size_t i=0; i < nodes.size(); ++i)
  {
    bool didDelete = m_bulk.destroy_relation(edge, nodes[i], ordinals[i]);
    STK_ThrowRequireMsg(didDelete, "failed to destroy relation");
  }
}

void EdgesOnFacesDeleter::delete_face_to_edge_relations(const std::vector<Entity>& faces, Entity edge)
{
  for (stk::mesh::Entity face : faces)
  {
    int idx = get_edge_idx(face, edge);
    bool didDelete = m_bulk.destroy_relation(face, edge, idx);
    STK_ThrowRequireMsg(didDelete, "failed to destroy relation");
  }
}

int EdgesOnFacesDeleter::get_edge_idx(Entity face, Entity edge)
{
  int idx = -1;
  stk::mesh::ConnectivityOrdinal const* downEdgeOrdinal = m_bulk.begin_edge_ordinals(face);
  for (stk::mesh::Entity const * downEdge = m_bulk.begin_edges(face); downEdge != m_bulk.end_edges(face); ++downEdge)
  {
    if (*downEdge == edge)
    {
      idx = *downEdgeOrdinal;
      break;
    }

    downEdgeOrdinal++;
  }
  STK_ThrowRequireMsg(idx >= 0, "unable to find edge");

  return idx;
}

}
}
}
