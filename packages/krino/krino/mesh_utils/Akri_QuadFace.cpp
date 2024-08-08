#include "Akri_QuadFace.hpp"

#include <type_traits>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include "Akri_MeshHelpers.hpp"

namespace krino {

static QuadFace quad_face_from_quad_node_offsets(stk::mesh::Entity::entity_value_type quadNodeOffset0,
    stk::mesh::Entity::entity_value_type quadNodeOffset1,
    stk::mesh::Entity::entity_value_type quadNodeOffset2,
    stk::mesh::Entity::entity_value_type quadNodeOffset3)
{
  static_assert(std::is_same<stk::mesh::Entity::entity_value_type, uint32_t>::value, "stk::mesh::Entity must be 32 bit.");
  const std::bitset<128> quadFaceValue = (std::bitset<128>(quadNodeOffset3) << 96) | (std::bitset<128>(quadNodeOffset2) << 64) | (std::bitset<128>(quadNodeOffset1) << 32) | std::bitset<128>(quadNodeOffset0);
  return QuadFace(quadFaceValue);
}

QuadFace quad_face_from_unordered_nodes(const stk::mesh::BulkData & mesh, const std::array<stk::mesh::Entity,4> & faceNodes)
{
  return quad_face_from_unordered_nodes(mesh, faceNodes[0], faceNodes[1], faceNodes[2], faceNodes[3]);
}

QuadFace quad_face_from_unordered_nodes(const stk::mesh::BulkData & mesh, stk::mesh::Entity faceNode0, stk::mesh::Entity faceNode1, stk::mesh::Entity faceNode2, stk::mesh::Entity faceNode3)
{
  std::array<stk::mesh::Entity,4> faceNodes{faceNode0, faceNode1, faceNode2, faceNode3};
  std::sort(faceNodes.begin(), faceNodes.end(), stk::mesh::EntityLess(mesh));
  return quad_face_from_quad_node_offsets(faceNodes[0].local_offset(), faceNodes[1].local_offset(), faceNodes[2].local_offset(), faceNodes[3].local_offset());
}

std::array<stk::mesh::Entity,4> get_quad_face_nodes_sorted_by_id(const QuadFace quadFace)
{
  static_assert(std::is_same<stk::mesh::Entity::entity_value_type, uint32_t>::value, "stk::mesh::Entity must be 32 bit.");
  static constexpr std::bitset<128> all32bit(0xFFFFFFFF);
  return std::array<stk::mesh::Entity, 4>{stk::mesh::Entity((quadFace.value() & all32bit).to_ulong()),
    stk::mesh::Entity(((quadFace.value() >> 32) & all32bit).to_ulong()),
    stk::mesh::Entity(((quadFace.value() >> 64) & all32bit).to_ulong()),
    stk::mesh::Entity((quadFace.value() >> 96).to_ulong())};
}

void append_entity_quad_faces(const stk::mesh::BulkData & mesh, const stk::topology entityTopology, const stk::mesh::Entity entity, std::vector<QuadFace> & entityQuadFaces)
{
  const unsigned numFaces = entityTopology.num_faces();

  const stk::mesh::Entity * entityNodes = mesh.begin_nodes(entity);
  std::array<stk::mesh::Entity,4> faceNodes;

  for (unsigned iFace = 0; iFace < numFaces; ++iFace)
  {
    entityTopology.side_nodes(entityNodes, iFace, faceNodes.data());
    entityQuadFaces.push_back(quad_face_from_unordered_nodes(mesh, faceNodes));
  }
}

void append_entity_quad_faces(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, std::vector<QuadFace> & entityQuadFaces)
{
  append_entity_quad_faces(mesh, mesh.bucket(entity).topology(), entity, entityQuadFaces);
}

void fill_entity_quad_faces(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, std::vector<QuadFace> & entityQuadFaces)
{
  const stk::topology entityTopology = mesh.bucket(entity).topology();
  const unsigned numFaces = entityTopology.num_edges();

  entityQuadFaces.clear();
  entityQuadFaces.reserve(numFaces);

  append_entity_quad_faces(mesh, entityTopology, entity, entityQuadFaces);
}

}
