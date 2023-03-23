#include "stk_middle_mesh/predicates/intersection_common.hpp"
#include "stk_middle_mesh/mesh_entity.hpp"

#include <stdexcept>

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

std::string get_enum_string(IntersectionType type)
{
  switch (type)
  {
    case IntersectionType::NONE: {
      return "none";
    }
    case IntersectionType::POINT: {
      return "point";
    }
    case IntersectionType::OVERLAP: {
      return "overlap";
    }
    default:
      throw std::invalid_argument("invalid IntersectionType Enum");
  }
}

template <typename T>
mesh::MeshEntityPtr get_entity_t(const T& record)
{
  assert(record.type != PointClassification::Exterior);
  if (record.type == PointClassification::Vert)
  {
    std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
    get_downward(record.el, 0, verts.data());
    return verts[record.id];
  } else if (record.type == PointClassification::Edge)
  {
    return record.el->get_down(record.id);
  } else
    return record.el;
}

mesh::MeshEntityPtr get_entity(const PointRecord& record)
{
  return get_entity_t(record);
}

mesh::MeshEntityPtr get_entity(const PointRecordForTriangle& record)
{
  return get_entity_t(record);
}



int get_entity_id(mesh::MeshEntityPtr el, mesh::MeshEntityPtr entity)
{
  assert(entity->get_type() == stk::middle_mesh::mesh::MeshEntityType::Vertex ||
         entity->get_type() == stk::middle_mesh::mesh::MeshEntityType::Edge);
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> entities;
  int nentities = get_downward(el, get_type_dimension(entity->get_type()), entities.data());
  for (int i = 0; i < nentities; ++i)
    if (entities[i] == entity)
      return i;

  throw std::runtime_error("could not find entity");
}

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
