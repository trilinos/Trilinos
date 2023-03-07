#ifndef MESH_ENTITY_H
#define MESH_ENTITY_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include <array>

#include "projection.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {

class MeshEntity;
using MeshEntityPtr = MeshEntity*;

const int MAX_DOWN = 4;

enum class MeshEntityType
{
  Vertex = 0,
  Edge,
  Triangle,
  Quad
};

constexpr int get_type_dimension(MeshEntityType type)
{
  std::array<int, 4> dims = {0, 1, 2, 2};
  return dims[int(type)];
}

std::string enum_to_string(MeshEntityType type);

inline std::ostream& operator<<(std::ostream& os, MeshEntityType type)
{
  os << enum_to_string(type);

  return os;
}

enum class EntityOrientation
{
  Standard = 0,
  Reversed
};

inline std::ostream& operator<<(std::ostream& os, EntityOrientation orient)
{
  if (orient == EntityOrientation::Standard)
    os << "Standard";
  else
    os << "Reversed";

  return os;
}

struct RemoteSharedEntity
{
    int remoteRank;
    int remoteId;
};

class MeshEntity
{
  public:
    MeshEntity(MeshEntityType type, int id, const std::vector<utils::Point>& pts = std::vector<utils::Point>(),
               const std::vector<MeshEntityPtr>& down       = std::vector<MeshEntityPtr>(),
               const std::vector<EntityOrientation>& orient = std::vector<EntityOrientation>())
      : m_type(type)
      , m_id(id)
      , m_down(down)
      , m_orientation(orient)
      , m_pointsOrig(pts)
    {
      switch (type)
      {
        case MeshEntityType::Vertex: {
          assert(down.size() == 0);
          break;
        }
        case MeshEntityType::Edge: {
          assert(down.size() == 2);
          break;
        }
        case MeshEntityType::Triangle: {
          assert(down.size() == 3);
          break;
        }
        case MeshEntityType::Quad: {
          assert(down.size() == 4);
          break;
        }
        default:
          throw std::invalid_argument("unrecognized MeshEntityType");
      }

      if (orient.size() == 0)
        for (unsigned int i = 0; i < m_down.size(); ++i)
          m_orientation.push_back(EntityOrientation::Standard);

      for (auto& downI : down)
        if (!downI)
          throw std::invalid_argument("entity cannot be null");
    }

    MeshEntityType get_type() const { return m_type; }

    int get_id() const { return m_id; }

    int count_down() const;

    MeshEntityPtr get_down(const int i);

    EntityOrientation get_down_orientation(const int i);

    void replace_down(const int i, MeshEntityPtr e);

    void delete_down(const int i);

    void delete_down(MeshEntityPtr e);

    int count_up() const;

    MeshEntityPtr get_up(const int i) const;

    void set_up(MeshEntityPtr e);

    void replace_up(const int i, MeshEntityPtr e);

    void replace_up(MeshEntityPtr eOld, MeshEntityPtr eNew);

    void delete_up(const int i) { delete_entity(i, m_up); }

    void delete_up(MeshEntityPtr e) { delete_entity(e, m_up); }

    void add_remote_shared_entity(const RemoteSharedEntity& remote) { m_remoteEntities.push_back(remote); }

    int count_remote_shared_entities() const { return m_remoteEntities.size(); }

    const RemoteSharedEntity& get_remote_shared_entity(int i);

    void delete_remote_shared_entity(int i);

    // gets the non-projected coordinates
    utils::Point get_point_orig(const int i)
    {
      assert(i >= 0 && i < int(m_pointsOrig.size()));
      return m_pointsOrig[i]; 
    }

    void set_point_orig(const int i, const utils::Point& pt);

    int count_points() { return m_pointsOrig.size(); }

    friend class Mesh;

  private:
    void set_id(const int id) { m_id = id; }

    void delete_entity(const int i, std::vector<MeshEntityPtr>& entities);

    void delete_entity(MeshEntityPtr e, std::vector<MeshEntityPtr>& entities);

    MeshEntityType m_type;
    int m_id;
    std::vector<MeshEntityPtr> m_down;
    std::vector<EntityOrientation> m_orientation;
    std::vector<MeshEntityPtr> m_up;
    std::vector<RemoteSharedEntity> m_remoteEntities;
    std::vector<utils::Point> m_pointsOrig;
};

inline bool is_less(const MeshEntityPtr& lhs, const MeshEntityPtr& rhs)
{
  int lhsDim = get_type_dimension(lhs->get_type());
  int rhsDim = get_type_dimension(rhs->get_type());
  if (lhsDim != rhsDim)
    return lhsDim < rhsDim;
  else
    return lhs->get_id() < rhs->get_id();
}

class MeshEntityCompare
{
  public:
    bool operator()(const MeshEntityPtr& lhs, const MeshEntityPtr& rhs) const { return is_less(lhs, rhs); }
};

std::ostream& operator<<(std::ostream& os, const MeshEntityPtr& pt);

} // namespace mesh

} // namespace middle_mesh
} // namespace stk
#endif
