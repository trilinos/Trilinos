#ifndef SIERRA_SIERRA_MESH_ENTITY_RANK_HPP
#define SIERRA_SIERRA_MESH_ENTITY_RANK_HPP

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace sierra {
namespace mesh {
namespace details {

class entity_rank {
  public:
    typedef size_t   rank_type;
    // section 4.7.2 of the C++ standard implies that
    // the maximum value of an unsigned type is equal to
    // static_cast<UnsignedType>(-1)
    // --the static_cast is necessary to suppress compiler warnings
    static const     rank_type INVALID_RANK = static_cast<rank_type>(-1);

    entity_rank()
      : m_rank( INVALID_RANK )
    {}

    explicit
    entity_rank(rank_type rank)
      : m_rank(rank)
    {}

    operator rank_type() const { return m_rank; }

    bool operator == ( const entity_rank rhs ) const
    { return m_rank == rhs.m_rank; }

    bool operator != ( const entity_rank rhs ) const
    { return m_rank != rhs.m_rank; }

    bool operator < ( const entity_rank rhs ) const
    { return m_rank < rhs.m_rank; }

    bool operator > ( const entity_rank rhs ) const
    { return m_rank > rhs.m_rank; }

    bool operator <= ( const entity_rank rhs ) const
    { return m_rank <= rhs.m_rank; }

    bool operator >= ( const entity_rank rhs ) const
    { return m_rank >= rhs.m_rank; }

    bool is_valid() const
    { return m_rank != INVALID_RANK; }

    // allow serialization
  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & m_rank;
    }

  private:
    rank_type m_rank;
};

inline
bool is_valid( entity_rank rank) { return rank.is_valid(); }

} // details
} // mesh
} // sierra

#endif // SIERRA_SIERRA_MESH_ENTITY_RANK_HPP
