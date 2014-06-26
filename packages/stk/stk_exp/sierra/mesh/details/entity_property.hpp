#ifndef STK_SIERRA_MESH_DETAILS_ENTITY_PROPERTY_HPP
#define STK_SIERRA_MESH_DETAILS_ENTITY_PROPERTY_HPP

#include <sierra/mesh/details/entity_rank.hpp>
#include <sierra/mesh/details/entity_id.hpp>

namespace sierra {
namespace mesh {
namespace details {

struct entity_property {

  entity_property(entity_rank r = entity_rank(), entity_id id = entity_id(), int proc=0)
   : m_rank(r), m_id(id), m_proc(proc) {}

  entity_rank m_rank;
  entity_id   m_id;
  int         m_proc;

  bool operator==(entity_property rhs) const
  { return m_rank == rhs.m_rank && m_id == rhs.m_id; }

  bool operator!=(entity_property rhs) const
  { return m_rank != rhs.m_rank || m_id != rhs.m_id; }

  bool operator<(entity_property rhs) const
  {
    if (m_rank != rhs.m_rank) return m_rank < rhs.m_rank;
    return m_id < rhs.m_id;
  }

  bool operator<=(entity_property rhs) const
  {
    return (*this < rhs) || (*this == rhs);
  }

  bool operator>(entity_property rhs) const
  {
    if (m_rank != rhs.m_rank) return m_rank > rhs.m_rank;
    return m_id > rhs.m_id;
  }

  bool operator>=(entity_property rhs) const
  { return (*this > rhs) || (*this == rhs);}
};

} // details
} // mesh
} // sierra

#endif //STK_SIERRA_MESH_DETAILS_ENTITY_PROPERTY_HPP

