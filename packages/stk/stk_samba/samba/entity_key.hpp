#ifndef SAMBA_SAMBA_ENTITY_KEY_HPP
#define SAMBA_SAMBA_ENTITY_KEY_HPP

#include <samba/utility.hpp>

#include <samba/entity_topology.hpp>
#include <samba/process_id.hpp>
#include <samba/entity_local_id.hpp>

#include <boost/static_assert.hpp>

namespace samba {

/**
 * entity_keys are opaque handles to entities that will refer to
 * the same entity for the duration of its life; entity_keys are
 * never invalidated except if the entity is destroyed or has its
 * owner changed (change-owner is considered to be a destroy/recreate).
 *
 * Certain bits of information are embedded in the entity_key to
 * allow for high-performing queries for commonly-accessed entity
 * properties. entity_key aggregates entity properties that are
 * fixed.
 *
 * entity_key is an aggregate POD type and does not support arithmetic
 * operations.
 */
struct entity_key
{
  struct tag
  {
    friend inline std::ostream& operator<<(std::ostream& out, tag)
    { return out << "entity_key"; }
  };

  typedef uint64_t value_type;
  static const int bits = 64;

  static entity_key invalid()
  {
    static const entity_key i = create( entity_topology::invalid()
                                       ,process_id::invalid()
                                       ,entity_local_id::invalid()
                                      );
    return i;
  }

  static const int topology_bits = entity_topology::num_bits;
  static const int process_bits  = process_id::num_bits;
  static const int local_id_bits = 32;

  BOOST_STATIC_ASSERT( topology_bits + process_bits + local_id_bits == sizeof(value_type)*8 );

  // pack the entity_topology, process_id, and entity_local_id
  // into the 64 bit entity_key.
  // entity_topology : 8  -> r
  // process_id      : 24 -> s
  // entity_local_id : 32 -> o
  static entity_key create(entity_topology t, process_id p, entity_local_id l)
  {
    uint64_t tt = t();
    tt <<= process_bits + local_id_bits;
    uint64_t pp = p();
    pp <<= local_id_bits;
    uint64_t ll = l();

    entity_key key = {{tt|pp|ll}};
    return key;
  }

  static entity_key create(value_type v)
  {
    entity_key i = {{v}};
    return i;
  }

  //get the topology
  entity_topology topology() const
  {
    static const uint64_t mask = 0xFF00000000000000ULL;
    entity_topology t = { static_cast<entity_topology::value_type>((m_value & mask) >> (process_bits + local_id_bits)) };
    return t;
  };

  //get the process ordinal
  process_id process() const
  {
    static const uint64_t mask = 0x00FFFFFF00000000ULL;
    process_id p = { static_cast<process_id::value_type>((m_value & mask) >> local_id_bits) };
    return p;
  }

  //get the ordinal
  entity_local_id local_id() const
  {
    static const uint64_t mask = 0x00000000FFFFFFFFULL;
    entity_local_id l = { static_cast<entity_local_id::value_type>(m_value & mask) };
    return l;
  }

  value_type   operator()() const { return m_value; }

  SAMBA_PRIMITIVE_COMPARISONS(entity_key,value_type)

  //serialize the entity_key
  template <typename Archive>
  void serialize(Archive &ar, unsigned version)
  { ar & m_value; }

  friend inline std::ostream& operator<<(std::ostream& out, entity_key d)
  {
    // {entity_key:topology,process,local_id}
    out << "{" << tag() << ":"
        << d.topology() << ","
        << d.process() << ","
        << d.local_id() << "}";
    return out;
  }

  struct big_endian_view
  {
    entity_topology::value_type topology : 8;
    process_id::value_type      process  : 24;
    entity_local_id::value_type local_id : 32;
  };

  struct little_endian_view
  {
    entity_local_id::value_type local_id : 32;
    process_id::value_type      process  : 24;
    entity_topology::value_type topology : 8;
  };

  union {
    value_type m_value;
    little_endian_view m_litte_endian; // for debugging
    big_endian_view m_big_endian;     //for debugging
  };
};

} // namespace samba

SAMBA_IS_PRIMITIVE(samba::entity_key)

#endif //SAMBA_SAMBA_ENTITY_KEY_HPP
