#ifndef SIERRA_SIERRA_MESH_INTRUSIVE_RELATION_HPP
#define SIERRA_SIERRA_MESH_INTRUSIVE_RELATION_HPP

#include <sierra/mesh/details/entity_key.hpp>
#include <sierra/mesh/details/relation_position.hpp>

namespace sierra {
namespace mesh {
namespace details {

struct entity_keyrelation {
  entity_keyrelation(
                     const relation_position & prop = relation_position(),
                     entity_key to = entity_key()
                   )
    : m_position(prop), m_to(to)
  {}

  bool operator==(entity_keyrelation rhs) const
  { return m_position == rhs.m_position &&
           m_to == rhs.m_to;
  }

  bool operator!=(entity_keyrelation rhs) const
  { return !(*this == rhs) ; }

  bool operator<(entity_keyrelation rhs) const
  { return m_position < rhs.m_position ||
           (! (rhs.m_position < m_position) && m_to < rhs.m_to);  }

  bool operator<=(entity_keyrelation rhs) const
  { return *this < rhs || *this == rhs; }

  relation_position  m_position;
  entity_key         m_to;
};

struct entity_descriptorrelation {
  entity_descriptorrelation(
                           const relation_position & prop = relation_position(),
                           entity_descriptor to = entity_descriptor()
                          )
    : m_position(prop), m_to(to)
  {}

  relation_position  m_position;
  entity_descriptor  m_to;
};

} // details
} // mesh
} // sierra

#endif //SIERRA_SIERRA_MESH_INTRUSIVE_RELATION_HPP
