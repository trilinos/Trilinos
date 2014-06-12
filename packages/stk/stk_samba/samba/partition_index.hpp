#ifndef SAMBA_SAMBA_PARTITION_INDEX_HPP
#define SAMBA_SAMBA_PARTITION_INDEX_HPP

#include <samba/partition_id.hpp>
#include <samba/partition_offset.hpp>
#include <samba/entity_rank.hpp>

#include <samba/utility/primitive_is_valid.hpp>

#include <boost/static_assert.hpp>

namespace samba {

/**
 * partition_indexs are opaque handles that refer to entities. These can
 * become stale/invalid if the mesh is modified.
 *
 * Certain bits of information are embedded in the partition_index to
 * allow for high-performing queries for commonly-accessed entity
 * properties. partition_index aggregates entity properties that
 * are not fixed.
 *
 * partition_index is an aggregate POD type and it does not support arithmetic
 * operations.
 */
struct partition_index
{
  struct tag
  {
    friend inline std::ostream& operator<<(std::ostream& out, tag)
    { return out << "partition_index"; }
  };

  typedef uint64_t value_type;
  static const int bits = 64;

  static const partition_index invalid()
  {
    static const partition_index i = create( entity_rank::invalid()
                                            ,partition_id::invalid()
                                            ,partition_offset::invalid()
                                           );
    return i;
  }

  static const int rank_bits      = entity_rank::num_bits;
  static const int partition_bits = partition_id::num_bits;
  static const int offset_bits    = 32;

  BOOST_STATIC_ASSERT( rank_bits + partition_bits + offset_bits == sizeof(value_type)*8 );

  // pack the entity_rank, partition_id, and partition_offset
  // into the 64 bit partition_index.
  // entity_rank    : 8  -> r
  // partition_id : 24 -> s
  // partition_offset  : 32 -> o
  static partition_index create(entity_rank r, partition_id p, partition_offset o)
  {
    uint64_t rr = r();
    rr <<= partition_bits + offset_bits;
    uint64_t pp = p();
    pp <<= offset_bits;
    uint64_t oo = o();

    partition_index descriptor = {{rr|pp|oo}};
    return descriptor;
  }

  static partition_index create(value_type v)
  {
    partition_index i = {{v}};
    return i;
  }

  //get the rank
  entity_rank rank() const
  {
    static const uint64_t mask = 0xFF00000000000000ULL;
    entity_rank t = { static_cast<entity_rank::value_type>((m_value & mask) >> (partition_bits + offset_bits)) };
    return t;
  };

  //get the partition ordinal
  partition_id partition() const
  {
    static const uint64_t mask = 0x00FFFFFF00000000ULL;
    partition_id s = { static_cast<partition_id::value_type>((m_value & mask) >> offset_bits) };
    return s;
  }

  //get the offset
  partition_offset offset() const
  {
    static const uint64_t mask = 0x00000000FFFFFFFFULL;
    partition_offset o = { static_cast<partition_offset::value_type>(m_value & mask) };
    return o;
  }

  value_type   operator()() const { return m_value; }

  SAMBA_PRIMITIVE_COMPARISONS(partition_index,value_type)

  //serialize the partition_index
  template <typename Archive>
  void serialize(Archive &ar, unsigned version)
  { ar & m_value; }

  friend inline std::ostream& operator<<(std::ostream& out, partition_index d)
  {
    // {partition_index:rank,partition,offset}
    out << "{" << tag() << ":"
        << d.rank() << ","
        << d.partition() << ","
        << d.offset() << "}";
    return out;
  }

  struct big_endian_view
  {
    entity_rank::value_type          rank      : 8;
    partition_id::value_type partition : 24;
    partition_offset::value_type        offset    : 32;
  };

  struct little_endian_view
  {
    partition_offset::value_type        offset    : 32;
    partition_id::value_type partition : 24;
    entity_rank::value_type          rank      : 8;
  };

  union {
    value_type m_value;
    little_endian_view m_litte_endian; // for debugging
    big_endian_view m_big_endian;     //for debugging
  };
};

} //namespace samba

SAMBA_IS_PRIMITIVE(samba::partition_index)

#endif //SAMBA_SAMBA_PARTITION_INDEX_HPP
