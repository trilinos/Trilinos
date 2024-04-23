#include "mesh_entity.hpp"
#include "mesh.hpp"

#include <stdexcept>

namespace stk {
namespace middle_mesh {
namespace mesh {

std::string enum_to_string(MeshEntityType type)
{
  switch (type)
  {
    case MeshEntityType::Vertex:
      return "vertex";
    case MeshEntityType::Edge:
      return "edge";
    case MeshEntityType::Triangle:
      return "triangle";
    case MeshEntityType::Quad:
      return "quad";
    default:
      throw std::invalid_argument("unrecognized MeshEntityType");
  }
}

int MeshEntity::count_down() const
{
  return m_down.size();
}

MeshEntityPtr MeshEntity::get_down(const int i)
{
#ifdef NDEBUG
  return m_down[i];
#else
  return m_down.at(i);
#endif
}

EntityOrientation MeshEntity::get_down_orientation(const int i)
{
#ifdef NDEBUG
  return m_orientation[i];
#else
  return m_orientation.at(i);
#endif
}

void MeshEntity::set_down_orientation(int i, EntityOrientation orient)
{
#ifdef NDEBUG
  m_orientation[i] = orient;
#else
  m_orientation.at(i) = orient;
#endif
}


void MeshEntity::replace_down(const int i, MeshEntityPtr e, EntityOrientation orient)
{
  if (!e)
    throw std::invalid_argument("entity must not be null");

  m_down[i] = e;
  if (m_type == MeshEntityType::Edge)
  {
    m_orientation[i] = orient;
  }
}

void MeshEntity::delete_down(const int i)
{
#ifdef NDEBUG
  m_down[i] = nullptr;
#else
  m_down.at(i) = nullptr;
#endif
}

void MeshEntity::delete_down(MeshEntityPtr e)
{
  for (unsigned int i = 0; i < m_down.size(); ++i)
    if (m_down[i] == e)
    {
      m_down[i] = nullptr;
      return;
    }

  throw std::invalid_argument("cannot delete downward adjacency that does not exist");
}

int MeshEntity::count_up() const
{
  return m_up.size();
}

MeshEntityPtr MeshEntity::get_up(const int i) const
{
#ifdef NDEBUG
  return m_up[i];
#else
  return m_up.at(i);
#endif
}

void MeshEntity::set_up(MeshEntityPtr e)
{
  if (!e)
    throw std::invalid_argument("entity must not be null");

  assert(get_type_dimension(e->get_type()) == get_type_dimension(m_type) + 1);
  if (std::find(m_up.begin(), m_up.end(), e) == m_up.end())
    m_up.push_back(e);
}

void MeshEntity::replace_up(const int i, MeshEntityPtr e)
{
#ifdef NDEBUG
  m_up[i] = e;
#else
  if (!e)
    throw std::invalid_argument("entity must not be null");

  if (e != m_up.at(i) && std::find(m_up.begin(), m_up.end(), e) != m_up.end())
    throw std::invalid_argument("cannot replace with duplicate entity");

  m_up.at(i)         = e;
#endif
}

void MeshEntity::replace_up(MeshEntityPtr eOld, MeshEntityPtr eNew)
{
  if (!eNew)
    throw std::invalid_argument("entity must not be null");

  for (unsigned int i = 0; i < m_up.size(); ++i)
    if (m_up[i] == eOld)
      replace_up(i, eNew);
}

const RemoteSharedEntity& MeshEntity::get_remote_shared_entity(int i)
{
  assert(i >= 0 && i < count_remote_shared_entities());
  return m_remoteEntities[i];
}

void MeshEntity::delete_remote_shared_entity(int i)
{
  assert(i >= 0 && i < count_remote_shared_entities());
  auto it = m_remoteEntities.begin() + i;
  m_remoteEntities.erase(it);
}

void MeshEntity::delete_entity(const int i, std::vector<MeshEntityPtr>& entities)
{
  auto it = entities.begin() + i;
  entities.erase(it);
}

void MeshEntity::delete_entity(MeshEntityPtr e, std::vector<MeshEntityPtr>& entities)
{
  for (unsigned int i = 0; i < entities.size(); ++i)
    if (entities[i] == e)
    {
      delete_entity(i, entities);
      return;
    }

  throw std::invalid_argument("could not find MeshEntity to delete");
}

void MeshEntity::set_point_orig(const int i, const utils::Point& pt)
{
#ifdef NDEBUG
  m_pointsOrig[i] = pt;
#else
  m_pointsOrig.at(i) = pt;
#endif
}

std::ostream& operator<<(std::ostream& os, const MeshEntityPtr& e)
{
  MeshEntityPtr verts[MAX_DOWN];
  int nverts = 1;
  verts[0]   = e;
  if (e->get_type() != MeshEntityType::Vertex)
    nverts = get_downward(e, 0, verts);

  os << "mesh " << enum_to_string(e->get_type()) << " with vertices";
  if (nverts > 1)
    os << "\n";

  for (int i = 0; i < nverts; ++i)
  {
    // os << "  projected point: " << verts[i]->getNode(0).pt << ", original pt " << verts[i]->get_point_orig(0);
    os << ", coords = " << verts[i]->get_point_orig(0);

    if (i != nverts - 1)
      os << std::endl;
  }

  return os;
}

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
