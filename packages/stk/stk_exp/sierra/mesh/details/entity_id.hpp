#ifndef STK_SIERRA_MESH_ENTITY_ID_HPP
#define STK_SIERRA_MESH_ENTITY_ID_HPP

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace sierra {
namespace mesh {
namespace details {

class entity_id {
  public:
    typedef int   id_type;
    static const     id_type INVALID_ID = -1;

    entity_id()
      : m_id( INVALID_ID )
    {}

    explicit
    entity_id(id_type id)
      : m_id(id)
    {}

    operator id_type() const { return m_id; }

    bool operator == ( const entity_id rhs ) const
    { return m_id == rhs.m_id; }

    bool operator != ( const entity_id rhs ) const
    { return m_id != rhs.m_id; }

    bool operator < ( const entity_id rhs ) const
    { return m_id < rhs.m_id; }

    bool operator > ( const entity_id rhs ) const
    { return m_id > rhs.m_id; }

    bool operator <= ( const entity_id rhs ) const
    { return m_id <= rhs.m_id; }

    bool operator >= ( const entity_id rhs ) const
    { return m_id >= rhs.m_id; }

    bool is_valid() const
    { return m_id != INVALID_ID; }

    // allow serialization
  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & m_id;
    }

  private:
    id_type m_id;
};

inline
bool is_valid( entity_id id) { return id.is_valid(); }

} // details
} // mesh
} // sierra

#endif // SIERRA_SIERRA_MESH_ENTITY_RANK_HPP
