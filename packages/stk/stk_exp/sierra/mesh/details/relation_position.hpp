#ifndef SIERRA_SIERRA_MESH_RELATION_POSITION_HPP
#define SIERRA_SIERRA_MESH_RELATION_POSITION_HPP

namespace sierra {
namespace mesh {
namespace details {

struct relation_position {
  size_t m_rank;
  size_t m_id;

  relation_position( size_t r=static_cast<size_t>(-1), size_t i=static_cast<size_t>(-1))
    : m_rank(r), m_id(i)
  {}

  bool operator==(relation_position rhs) const
  { return m_rank == rhs.m_rank &&
           m_id == rhs.m_id;
  }

  bool operator!=(relation_position rhs) const
  { return !(*this == rhs) ; }

  bool operator<(relation_position rhs) const
  { return m_rank < rhs.m_rank ||
           (! (rhs.m_rank < m_rank) && m_id < rhs.m_id);  }

  bool operator<=(relation_position rhs) const
  { return *this < rhs || *this == rhs; }

};

} // details
} // mesh
} // sierra

#endif //SIERRA_SIERRA_MESH_RELATION_POSITION_HPP
