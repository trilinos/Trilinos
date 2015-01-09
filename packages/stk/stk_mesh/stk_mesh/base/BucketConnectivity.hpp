// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_MESH_BUCKET_CONNECTIVITY_HPP
#define STK_MESH_BUCKET_CONNECTIVITY_HPP

#include "stk_mesh/base/Types.hpp"      // for ConnectivityOrdinal, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include "stk_util/environment/ReportHandler.hpp"
#include <stk_util/util/TrackingAllocator.hpp>

namespace stk {
namespace mesh {

namespace impl {

struct LowerConnectivityCompare
{
  bool operator()(Entity first_entity, ConnectivityOrdinal first_ordinal, Entity second_entity, ConnectivityOrdinal second_ordinal) const
  {
    // only compare ordinals
    return first_ordinal < second_ordinal;
  }
};

template <typename BulkData>
struct LowerConnectivitityRankSensitiveCompare
{
  LowerConnectivitityRankSensitiveCompare(const BulkData &bulk_data) : m_mesh(bulk_data) { }

  const BulkData &m_mesh;

  bool operator()(Entity first_entity, ConnectivityOrdinal first_ordinal,
                  Entity second_entity, ConnectivityOrdinal second_ordinal) const;
};

struct HigherConnectivityCompare
{
  bool operator()(Entity first_entity, ConnectivityOrdinal first_ordinal, Entity second_entity, ConnectivityOrdinal second_ordinal) const
  {
    // Needs to match LessRelation in BulkData.hpp
    return std::make_pair(first_ordinal,  first_entity.is_local_offset_valid() ?  first_entity.local_offset()  : Entity::MaxEntity) <
           std::make_pair(second_ordinal, second_entity.is_local_offset_valid() ? second_entity.local_offset() : Entity::MaxEntity);
  }
};

template <typename BulkData>
struct HigherConnectivityRankSensitiveCompare
{
  HigherConnectivityRankSensitiveCompare(const BulkData &bulk_data) : m_mesh(bulk_data) { }

  const BulkData &m_mesh;

  bool operator()(Entity first_entity, ConnectivityOrdinal first_ordinal, Entity second_entity,
                  ConnectivityOrdinal second_ordinal) const;
};

template <typename Connectivity>
inline void check_bucket_ordinal(unsigned bucket_ordinal, Connectivity const* connectivity)
{
  ThrowAssertMsg(bucket_ordinal < connectivity->size(),
                 "bucket_ordinal " << bucket_ordinal << " is out of range, bucket size is " << connectivity->size());
}

template<EntityRank TargetRank, ConnectivityType >
class BucketConnectivity;

template<EntityRank TargetRank >
class BucketConnectivity<TargetRank, FIXED_CONNECTIVITY>
{
 public:
  typedef BucketConnectivity<TargetRank, FIXED_CONNECTIVITY> SelfType;
  typedef BucketConnectivity<TargetRank, DYNAMIC_CONNECTIVITY> OtherType;

  static const EntityRank target_rank = TargetRank;
  static const ConnectivityType connectivity_type = FIXED_CONNECTIVITY;

  typedef TrackedVectorMetaFunc<Entity, BucketRelationTag>::type              EntityVector;
  typedef TrackedVectorMetaFunc<ConnectivityOrdinal, BucketRelationTag>::type ConnectivityOrdinalVector;
  typedef TrackedVectorMetaFunc<Permutation, BucketRelationTag>::type         PermutationVector;

  BucketConnectivity() //default constructed BucketConnectivity implies connectivity is not used
    : m_num_connectivity(0u)
    , m_targets()
    , m_ordinals()
    , m_permutations()
  {}

  BucketConnectivity(unsigned arg_num_connectivity)
    : m_num_connectivity(0)
    , m_targets()
    , m_ordinals()
    , m_permutations()
  {
    set_num_connectivity(arg_num_connectivity);
  }

  void set_num_connectivity(unsigned arg_num_connectivity)
  {
    ThrowAssertMsg(m_num_connectivity == 0, "Cannot reset num_connectivity");
    ThrowAssertMsg(arg_num_connectivity != 0, "Cannot set num connectivity to 0 for fixed connectivity");

    m_num_connectivity = arg_num_connectivity;

    // Ordinal is the same for all fixed, just counts up
    m_ordinals.resize(m_num_connectivity);
    for (ConnectivityOrdinal ord=static_cast<ConnectivityOrdinal>(0); ord < static_cast<ConnectivityOrdinal>(m_num_connectivity); ++ord) {
      m_ordinals[ord] = ord;
    }
  }

  // Entity iterator

  Entity const* begin(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[bucket_ordinal * m_num_connectivity]; }

  Entity      * begin(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[bucket_ordinal * m_num_connectivity]; }

  Entity const* end(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[(bucket_ordinal + 1) * m_num_connectivity]; }

  Entity      * end(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[(bucket_ordinal + 1) * m_num_connectivity]; }

  // Ordinal iterator

  ConnectivityOrdinal const* begin_ordinals(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_ordinals[0]; }

  ConnectivityOrdinal      * begin_ordinals(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_ordinals[0]; }

  ConnectivityOrdinal const* end_ordinals(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_ordinals[0] + m_num_connectivity; }

  ConnectivityOrdinal      * end_ordinals(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_ordinals[0] + m_num_connectivity; }

  // Permutation iterator

  Permutation const* begin_permutations(unsigned bucket_ordinal) const
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[bucket_ordinal * m_num_connectivity];
  }

  Permutation      * begin_permutations(unsigned bucket_ordinal)
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[bucket_ordinal * m_num_connectivity];
  }

  Permutation const* end_permutations(unsigned bucket_ordinal) const
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[(bucket_ordinal + 1) * m_num_connectivity];
  }

  Permutation      * end_permutations(unsigned bucket_ordinal)
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[(bucket_ordinal + 1) * m_num_connectivity];
  }

  // Queries

  unsigned num_connectivity(unsigned /*bucket_ordinal*/) const
  { return m_num_connectivity; }

  // return number of entities
  unsigned size() const
  { return m_targets.size() / m_num_connectivity; }

  // Modification API

  bool add_connectivity(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal, Permutation permutation = INVALID_PERMUTATION)
  {
    ThrowAssertMsg(ordinal < m_num_connectivity,
                   "Ordinal " <<  ordinal << " exceeds topological limit: " << m_num_connectivity);
    impl::check_bucket_ordinal(bucket_ordinal, this);
#ifndef NDEBUG
    // TODO - Add topology invariant; target entity should match a sub topology
#endif

    unsigned index = m_num_connectivity*bucket_ordinal + ordinal;

    if (m_targets[index] == to) {
      ThrowAssert(!has_permutation() || m_permutations[index] == permutation);
      // Already exists
      return false;
    }

    m_targets[index] = to;

    if (has_permutation()) {
      m_permutations[index] = permutation;
    }

    invariant_check_helper(bucket_ordinal);

    return true;
  }

  bool remove_connectivity(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal)
  {
    ThrowAssertMsg(ordinal < m_num_connectivity,
                   "Ordinal " <<  ordinal << " exceeds topological limit: " << m_num_connectivity);
    impl::check_bucket_ordinal(bucket_ordinal, this);

    unsigned index = m_num_connectivity*bucket_ordinal + ordinal;
    if (m_targets[index] != to) {
      return false;
    }

    // Clear
    m_targets[index] = Entity::InvalidEntity;
    if (has_permutation()) {
      m_permutations[index] = INVALID_PERMUTATION;
    }

    invariant_check_helper(bucket_ordinal);

    return true;
  }

  void begin_modification()
  {}

  template <typename BulkData> // hack to get around dependency
  void end_modification(BulkData* mesh = NULL);

  void add_entity()
  {
    const unsigned new_conn_size = m_targets.size() + m_num_connectivity;
    Entity invalid = {Entity::InvalidEntity};
    m_targets.resize(new_conn_size, invalid); // Not a perf issue: vectors are smart when resizing
    if (has_permutation()) {
      m_permutations.resize(new_conn_size, INVALID_PERMUTATION);
    }

    invariant_check_helper();
  }

  // Always removes last entity
  void remove_entity()
  {
    ThrowAssertMsg(size() > 0, "Cannot remove, connectivity is already empty");

    const unsigned new_conn_size = m_targets.size() - m_num_connectivity;
    m_targets.resize(new_conn_size);
    if (has_permutation()) {
      m_permutations.resize(new_conn_size);
    }

    invariant_check_helper();
  }

  void copy_entity(unsigned from_ordinal, SelfType& to, unsigned to_ordinal=-1u)
  {
    ThrowAssertMsg(size() > 0, "Cannot move, connectivity is empty");

    if (to_ordinal == -1u) { // check if we should just append
      to_ordinal = to.size();
      to.add_entity(); // make room for new entity
    }

    // Copy connectivity to other BucketConnectivity
    copy_connectivity(from_ordinal, to, to_ordinal);

    invariant_check_helper();
    to.invariant_check_helper();
  }

  void copy_to_fixed(unsigned from_ordinal, SelfType& to)
  { ThrowAssert(false); }

  void copy_to_fixed(unsigned from_ordinal, OtherType& to)
  { ThrowAssert(false); }

  bool has_permutation() const
  {
    const static bool rv = TargetRank != stk::topology::NODE_RANK;
    return rv;
  }

  void debug_dump(std::ostream& out) const
  {
    out << "For fixed connectivity to rank: " << TargetRank << "\n";
    out << "  size is: " << size() << "\n";
    for (int i = 0, ie = size(); i < ie; ++i) {
      out << "    At ordinal " << i << "\n";
      debug_dump(out, i);
    }
    out << std::endl;
  }

  void debug_dump(std::ostream& out, unsigned ordinal) const
  {
    for (int j = m_num_connectivity * ordinal, je = m_num_connectivity*(ordinal+1); j < je; ++j) {
      out << "        target:" << m_targets[j].local_offset() << "\n";
    }
  }

private:

  void copy_connectivity(unsigned from_ordinal, SelfType& to, unsigned to_ordinal)
  {
    unsigned to_offset   = to_ordinal * m_num_connectivity;
    unsigned from_offset = from_ordinal * m_num_connectivity;

    std::copy(m_targets.begin() + from_offset,
              m_targets.begin() + from_offset + m_num_connectivity,
              to.m_targets.begin() + to_offset);

    if (has_permutation()) {
      std::copy(m_permutations.begin() + from_offset,
                m_permutations.begin() + from_offset + m_num_connectivity,
                to.m_permutations.begin() + to_offset);
    }
  }

  void invariant_check_helper(unsigned bucket_ordinal) const
  {

  #ifdef STK_INVARIANCE_CHECK
    const Entity* keys_begin = begin(bucket_ordinal);
    const Entity* keys_end   = end(bucket_ordinal);
    const ConnectivityOrdinal* ordinals_begin = begin_ordinals(bucket_ordinal);
    const ConnectivityOrdinal* ordinals_end   = end_ordinals(bucket_ordinal);
    const Permutation* permutations_begin = begin_permutations(bucket_ordinal);
    const Permutation* permutations_end   = end_permutations(bucket_ordinal);

    ThrowAssertMsg(keys_end - keys_begin == m_num_connectivity,
                   "Expected data to be of size " << m_num_connectivity << ", " << bucket_ordinal << " has keys " << keys_end - keys_begin);

    ThrowAssertMsg(keys_end - keys_begin == ordinals_end - ordinals_begin,
                   "Num keys, " << keys_end - keys_begin << ", does not match num ordinals, " << ordinals_end - ordinals_begin);
    if (has_permutation()) {
      ThrowAssertMsg(keys_end - keys_begin == permutations_end - permutations_begin,
                     "Num keys, " << keys_end - keys_begin << ", does not match num permutations, " << permutations_end - permutations_begin);
    }
    else {
      ThrowAssertMsg(permutations_end - permutations_begin == 0,
                     "Expected 0 permutations for node connectivity, found: " << permutations_end - permutations_begin);
    }

    const Entity*               kitr = keys_begin;
    const ConnectivityOrdinal*  oitr = ordinals_begin;
    const Permutation*          pitr = permutations_begin;
    for ( ; kitr != keys_end; ++kitr, ++oitr) {
      if (*kitr != Entity()) {
        ThrowAssertMsg(*oitr == kitr - keys_begin,
                       "For bucket_ordinal " << bucket_ordinal << ", connectivity to entity " << kitr->local_offset() <<
                       ", found out-of-order connectivity at index " << kitr - keys_begin << ", its ordinal is " << *oitr);
        // TODO
        //entity_rank to_rank  = topology_rank(kitr->topology(), m_spatial_dimension);
        //ThrowAssertMsg(to_rank() == TargetRank,
        //                 (debug_message() << "Found connectivity to wrong rank " << to_rank << ", expected " << entity_rank::create(TargetRank)));
      }
      else {
        if (has_permutation()) {
          ThrowAssertMsg(*pitr == INVALID_PERMUTATION, "If key is invalid, then permutation should be invalid");
        }
      }

      if (has_permutation()) {
        ++pitr;
      }
      // TODO - Anything else we can check here?
    }
#endif

  }

  void invariant_check_helper() const
  {

  #ifdef STK_INVARIANCE_CHECK
    ThrowAssertMsg(static_cast<unsigned>(m_ordinals.size()) == m_num_connectivity,
                   "Total size of ordinals " << m_ordinals.size() << " does not match num_connectivity " << m_num_connectivity);

    if (has_permutation()) {
      ThrowAssertMsg(m_permutations.size() == m_targets.size(),
                     "Total size of permutations " << m_permutations.size() << " does not match size of keys " << m_targets.size());
    }
    else {
      ThrowAssertMsg(m_permutations.empty(), "Permutations should be empty for nodal connectivity");
    }
  #endif

  }

  // Call this at the end of modification
  template <typename BulkData>
  void invariant_check_helper(BulkData* mesh = NULL) const
  {

  #ifdef STK_INVARIANCE_CHECK
    invariant_check_helper();

    if (mesh != NULL) {
      for (size_t i = 0, e = m_targets.size(); i < e; ++i) {
        Entity entity = m_targets[i];
        if (mesh->is_valid(entity)) {
          //        // TODO
          //        const EntityKey key_converted_from_partition_index =
          //            mesh->convert<EntityKey>(m_target_by_partition_index[i]);
          //        ThrowAssertMsg(key == key_converted_from_partition_index,
          //            (debug_message() << "Key converted from partition index " << key_converted_from_partition_index
          //                << " does not match expected key " << key));
        }
      }
    }
  #endif

  }

  // Illegal
  BucketConnectivity(const SelfType&);
  SelfType& operator=(const SelfType&);

  // MEMBERS

  unsigned m_num_connectivity;

  // connectivity data
  EntityVector              m_targets;
  ConnectivityOrdinalVector m_ordinals; // shared for all entities
  PermutationVector         m_permutations;

  // friend OtherType; // 1337! Will have to wait for c++11
  friend class BucketConnectivity<TargetRank, DYNAMIC_CONNECTIVITY>;
};

// Want a way for all dynamic connectivity instantiations to share the same id space
struct Counter
{
  static int counter;
};

// Uncomment to enable profiling
//#define STK_MESH_ANALYZE_DYN_CONN


// Profiling data for an individual dynamic connectivity object
struct DynConnData
{
  // from-rank for the associated dynamic connectivity
  EntityRank m_from_rank;

  // to-rank for the associated dynamic connectivity
  EntityRank m_to_rank;

  // the maximum capacity ever achieved by connectivity vectors
  size_t m_max_capacity;

  // at the point at which maximum capacity (member above) was achieved, how much memory
  // was lost due to "abandoned space".
  // "abandoned space" - When an entity overflows its current chunks, it gets additional
  // chunks but must be copied to the end. The space left behind is abandoned and will
  // not be reused until the next compress (resize_and_order_by_index).
  size_t m_abandoned_space;

  // at the point at which maximum capacity (member above) was achieved, how much memory
  // was lost due to unused chunk capacity. If chunk size is > 1, it's possible that an
  // entity is not using all the space available in it's chunk. For example, if chunk size
  // is 8 and an entity has 5 connectivities, then unused chunk capacity is 3 for that
  // entity. This member stores the sum over all entities.
  size_t m_unused_chunk_capacity;

  // The number of times this dynamic connectivity had to be grown
  size_t m_num_growths;

  // The number of times any entity overflowed it's chunk allocation and had to be
  // copied to the end
  size_t m_num_entity_relocations;

  // at the point at which maximum capacity (member above) was achieved, what is the
  // total amount of wasted memory
  size_t m_total_unused_memory;

  // at the point at which maximum capacity (member above) was achieved, what is the
  // amount of memory that is wasted due to vector capacity growth over-provisioning.
  size_t m_unused_capacity;

  // at the point at which maximum capacity (member above) was achieved, what is the
  // number of connectivity being stored.
  size_t m_total_num_conn;

  DynConnData(EntityRank from_rank, EntityRank to_rank) :
    m_from_rank(from_rank),
    m_to_rank(to_rank),
    m_max_capacity(0),
    m_abandoned_space(0),
    m_unused_chunk_capacity(0),
    m_num_growths(0),
    m_num_entity_relocations(0),
    m_total_unused_memory(0),
    m_unused_capacity(0),
    m_total_num_conn(0)
  {}
};

struct DynConnMetrics
{
  static std::vector<DynConnData> m_data;
};

template<EntityRank TargetRank >
class BucketConnectivity<TargetRank, DYNAMIC_CONNECTIVITY>
{
  enum connectivity_direction { Lower=0,Higher=1,Adjacent=2 };

public:
  typedef BucketConnectivity<TargetRank, DYNAMIC_CONNECTIVITY> SelfType;
  typedef BucketConnectivity<TargetRank, FIXED_CONNECTIVITY> OtherType;

  static const EntityRank target_rank = TargetRank;
  static const ConnectivityType connectivity_type = DYNAMIC_CONNECTIVITY;

  typedef typename DynamicConnectivityTagSelector<TargetRank>::type TagType;

  typedef typename TrackedVectorMetaFunc<Entity,              TagType>::type              EntityVector;
  typedef typename TrackedVectorMetaFunc<ConnectivityOrdinal, TagType>::type ConnectivityOrdinalVector;
  typedef typename TrackedVectorMetaFunc<Permutation,         TagType>::type         PermutationVector;
  typedef typename TrackedVectorMetaFunc<uint32_t,            TagType>::type              UInt32Vector;
  typedef typename TrackedVectorMetaFunc<uint16_t,            TagType>::type              UInt16Vector;

  static const unsigned chunk_size = 1u;

  BucketConnectivity(EntityRank from_rank, BulkData *bulk_data)
    : m_from_rank(from_rank)
    , m_direction( (m_from_rank > TargetRank) ? Lower : ((m_from_rank == TargetRank) ? Adjacent : Higher))
    , m_active(false)
    , m_needs_shrink_to_fit(false)
    , m_num_inactive(0)
    , m_indices()
    , m_num_connectivities()
    , m_total_connectivities(0)
    , m_targets()
    , m_ordinals()
    , m_permutations()
    , m_bulk_data(bulk_data)
    , m_id(Counter::counter++)
    , m_rank_sensitive_higher_connectivity_cmp(*m_bulk_data)
    , m_rank_sensitive_lower_connectivity_cmp(*m_bulk_data)
    , m_last_capacity(0)
  {
#ifdef STK_MESH_ANALYZE_DYN_CONN
    DynConnMetrics::m_data.push_back(DynConnData(from_rank, target_rank));
    m_data_idx = DynConnMetrics::m_data.size() - 1;
#endif
  }

  // Entity iterator

  Entity const* begin(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[m_active ? m_indices[bucket_ordinal] : 0]; }

  Entity      * begin(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[m_active ? m_indices[bucket_ordinal] : 0]; }

  Entity const* end(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return begin(bucket_ordinal) + num_connectivity(bucket_ordinal); }

  Entity      * end(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return begin(bucket_ordinal) + num_connectivity(bucket_ordinal); }

  // Ordinal iterator

  ConnectivityOrdinal const* begin_ordinals(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_ordinals[m_active ? m_indices[bucket_ordinal] : 0]; }

  ConnectivityOrdinal      * begin_ordinals(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_ordinals[m_active ? m_indices[bucket_ordinal] : 0]; }

  ConnectivityOrdinal const* end_ordinals(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return begin_ordinals(bucket_ordinal) + num_connectivity(bucket_ordinal); }

  ConnectivityOrdinal      * end_ordinals(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return begin_ordinals(bucket_ordinal) + num_connectivity(bucket_ordinal); }

  // Permutation iterator

  Permutation const* begin_permutations(unsigned bucket_ordinal) const
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[m_active ? m_indices[bucket_ordinal] : 0];
  }

  Permutation      * begin_permutations(unsigned bucket_ordinal)
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[m_active ? m_indices[bucket_ordinal] : 0];
  }

  Permutation const* end_permutations(unsigned bucket_ordinal) const
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return begin_permutations(bucket_ordinal) + num_connectivity(bucket_ordinal);
  }

  Permutation      * end_permutations(unsigned bucket_ordinal)
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return begin_permutations(bucket_ordinal) + num_connectivity(bucket_ordinal);
  }

  // Queries

  unsigned num_connectivity(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return m_active ? m_num_connectivities[bucket_ordinal] : 0; }

  // return number of entities
  unsigned size() const
  { return m_active ? m_indices.size() : m_num_inactive; }

  // Modification API

  bool add_connectivity(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal, Permutation permutation = INVALID_PERMUTATION)
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);

    m_needs_shrink_to_fit = true;

    if (!m_active) {
      activate();
    }

    if (target_rank <= stk::topology::ELEMENT_RANK) {
      switch(m_direction)
      {
      case Lower: return add_helper(bucket_ordinal, to, ordinal, permutation, LowerConnectivityCompare());
      case Higher: return add_helper(bucket_ordinal, to, ordinal, permutation, HigherConnectivityCompare());
      case Adjacent: return add_helper(bucket_ordinal, to, ordinal, permutation, LowerConnectivityCompare()); // same comparing as lower
      default:
        ThrowAssertMsg(false, "What type of connectivity are you trying to add? " << m_direction);
        return false;
      }
    }
    else {
      switch(m_direction)
      {
      case Lower: return add_helper(bucket_ordinal, to, ordinal, permutation, m_rank_sensitive_lower_connectivity_cmp);
      case Higher: return add_helper(bucket_ordinal, to, ordinal, permutation, m_rank_sensitive_higher_connectivity_cmp);
      case Adjacent: return add_helper(bucket_ordinal, to, ordinal, permutation, m_rank_sensitive_lower_connectivity_cmp);
      default:
        ThrowAssertMsg(false, "What type of connectivity are you trying to add? " << m_direction);
        return false;
      }
    }
  }

  bool remove_connectivity(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal)
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);

    if (!m_active) return false;

    uint32_t found_idx = ~0u;
    const uint32_t end_i = m_indices[bucket_ordinal]+m_num_connectivities[bucket_ordinal];
    for (uint32_t i = m_indices[bucket_ordinal]; i < end_i; ++i)
    {
      //remove connectivity
      if ( m_targets[i] == to && m_ordinals[i] == ordinal ) {
        found_idx = i;
        --m_num_connectivities[bucket_ordinal];
        --m_total_connectivities;
        break;
      }
    }

    //slide memory down
    if (found_idx != ~0u) {
      m_needs_shrink_to_fit = true;
      for (uint32_t i = found_idx; i < end_i - 1; ++i) {
        m_targets[i] = m_targets[i+1];
        m_ordinals[i]        = m_ordinals[i+1];
        if (has_permutation()) {
          m_permutations[i]  = m_permutations[i+1];
        }
      }
    }

    invariant_check_helper(bucket_ordinal);

    return found_idx != ~0u;
  }

  void begin_modification()
  {}

  template <typename BulkData>
  void end_modification(BulkData* mesh = NULL);

  void add_entity()
  {
    if (m_active) {
      m_indices.push_back(m_targets.size());
      m_num_connectivities.push_back(0);
      m_needs_shrink_to_fit = true;
    }
    else {
      ++m_num_inactive;
    }

    invariant_check_helper();
  }

  void remove_entity()
  {
    ThrowAssertMsg(size() > 0, "Cannot remove, connectivity is already empty");

    if (m_active) {
      m_indices.pop_back();
      m_total_connectivities -= m_num_connectivities.back();
      m_num_connectivities.pop_back();
      m_needs_shrink_to_fit = true;
    }
    else {
      --m_num_inactive;
    }

    invariant_check_helper();
  }

  void copy_entity(unsigned from_ordinal, SelfType& to, unsigned to_ordinal=-1u)
  {
    ThrowAssert(m_from_rank == to.m_from_rank);
    impl::check_bucket_ordinal(from_ordinal, this);
    invariant_check_helper(from_ordinal);

    if (to_ordinal == -1u) {
      to_ordinal = to.size();
      to.add_entity();
    }
    impl::check_bucket_ordinal(to_ordinal, &to);

    // Manage activation state
    if (!m_active) {
      if (to.m_active) {
        to.m_total_connectivities -= to.m_num_connectivities[to_ordinal];
        to.m_num_connectivities[to_ordinal] = 0;
        to.m_needs_shrink_to_fit = true;
      }
      return;
    }
    if (m_active && !to.m_active) {
      to.activate();
    }

    // Copy data
    if (&to == this) {
      // easy
      // note this implements swap semantics instead of copy, but this is necessary
      // to avoid aliasing in certain situations when sorting Partitions
      std::swap(m_indices[to_ordinal], m_indices[from_ordinal]);
      std::swap(m_num_connectivities[to_ordinal], m_num_connectivities[from_ordinal]);

      m_needs_shrink_to_fit = true;
      to.m_needs_shrink_to_fit = true;
    }
    else {
      // much harder
      const unsigned from_num    = m_num_connectivities[from_ordinal];
      const unsigned to_num   = to.m_num_connectivities[to_ordinal];
      const int delta_num     = from_num - to_num;
      if (delta_num > 0) {
        // If adding additional connectivity, need to reserve space
        to.add_connectivity_helper(to_ordinal, delta_num);
      }
      else {
        to.m_num_connectivities[to_ordinal] = from_num;
        to.m_total_connectivities += delta_num;
      }
      if (delta_num != 0) {
        to.m_needs_shrink_to_fit = true;
      }
      copy_connectivity(from_ordinal, to, to_ordinal);
    }

    invariant_check_helper();
    to.invariant_check_helper();
    to.invariant_check_helper(to_ordinal);
  }

  void copy_to_fixed(unsigned from_ordinal, OtherType& to)
  {
    const unsigned num_conn_to_move = m_active ? m_num_connectivities[from_ordinal] : 0;

    ThrowAssert(OtherType::connectivity_type == FIXED_CONNECTIVITY);
    ThrowAssertMsg(size() > 0, "Cannot move, connectivity is empty");
    ThrowAssertMsg(num_conn_to_move <= to.num_connectivity(666 /*any unsigned, doesn't matter*/), "Incompatible");

    const unsigned to_offset = to.m_targets.size();
    to.add_entity(); // make room for new entity

    const unsigned from_offset = m_active ? m_indices[from_ordinal] : 0;

#ifndef NDEBUG
    // Check the ordinals are compatible with fixed connectivity
    ConnectivityOrdinal const* ordinals = &m_ordinals[0] + from_offset;
    for (unsigned i = 0; i < num_conn_to_move; ++i) {
      ThrowAssert(ordinals[i] == i);
    }
#endif

    std::copy(m_targets.begin() + from_offset,
              m_targets.begin() + from_offset + num_conn_to_move,
              to.m_targets.begin() + to_offset);

    if (has_permutation()) {
      std::copy(m_permutations.begin() + from_offset,
                m_permutations.begin() + from_offset + num_conn_to_move,
                to.m_permutations.begin() + to_offset);
    }

    invariant_check_helper();
    to.invariant_check_helper();
  }

  void copy_to_fixed(unsigned from_ordinal, SelfType& to)
  { ThrowAssert(false); }

  bool has_permutation() const
  { return TargetRank != stk::topology::NODE_RANK && m_from_rank != stk::topology::NODE_RANK; }

  void debug_dump(std::ostream& out) const
  {
    out << "For dynamic connectivity to rank: " << TargetRank << ", with id: " << m_id << "\n";
    if (m_active) {
      out << "  size is: " << m_indices.size() << "\n";
      for (int i = 0, ie = m_indices.size(); i < ie; ++i) {
        out << "    At ordinal " << i << "\n";
        debug_dump(out, i, false);
      }
      out << std::endl;
    }
    else {
      out << "  size is: " << m_num_inactive << ", but inactive" << std::endl;
    }
  }

  void debug_dump(std::ostream& out, unsigned ordinal, bool add_context=true) const
  {
    if (m_active) {
      int idx = m_indices[ordinal];
      int num = m_num_connectivities[ordinal];
      if (add_context) {
        out << "For dynamic connectivity to rank: " << TargetRank << ", with id: " << m_id << "\n";
      }
      out << "      Index is: " << idx << ", Num is: " << num << "\n";
      for (int j = idx, je = idx + num; j < je; ++j) {
        out << "        (target:" << m_targets[j].local_offset() << ", ordinal:" << m_ordinals[j] << ")\n";
      }
    }
    else {
      out << "      Index is: 0, Num is: 0\n";
    }
  }

private:

  DynConnData& profile_data() const
  {
    return DynConnMetrics::m_data[m_data_idx];
  }

  void copy_connectivity(unsigned from_ordinal, SelfType& to, unsigned to_ordinal)
  {
    unsigned num_conn    = m_num_connectivities[from_ordinal];
    unsigned to_offset   = to.m_indices[to_ordinal];
    unsigned from_offset = m_indices[from_ordinal];
    ThrowAssert(to.m_num_connectivities[to_ordinal] == num_conn);

    std::copy(m_targets.begin() + from_offset,
              m_targets.begin() + from_offset + num_conn,
              to.m_targets.begin() + to_offset);

    std::copy(m_ordinals.begin() + from_offset,
              m_ordinals.begin() + from_offset + num_conn,
              to.m_ordinals.begin() + to_offset);

    if (has_permutation()) {
      std::copy(m_permutations.begin() + from_offset,
                m_permutations.begin() + from_offset + num_conn,
                to.m_permutations.begin() + to_offset);
    }
  }

  static unsigned num_chunks(unsigned num)
  { return (num + chunk_size -1)/chunk_size; }

  void activate()
  {
    ThrowAssert(!m_active);

    m_indices.resize(m_num_inactive, 0);
    m_num_connectivities.resize(m_num_inactive, 0);

    m_active = true;
    m_num_inactive = 0;
  }

  template <typename Vector>
  void resize_and_order_by_index_helper(Vector & data, unsigned capacity, bool update_index = false)
  {
    Vector temp;
    temp.reserve(capacity);

    uint32_t current_index=0;
    for(size_t i=0, e=m_indices.size(); i<e; ++i)
    {
      const uint32_t entity_data_size = num_chunks(m_num_connectivities[i]) * chunk_size;

      const uint32_t begin_offset = m_indices[i];
      const uint32_t end_offset   = begin_offset + entity_data_size;

      if (begin_offset != end_offset) {
        temp.insert(temp.end(), data.begin() + begin_offset, data.begin() + end_offset);
      }

      if (update_index) {
        m_indices[i] = current_index;
      }

      current_index += entity_data_size;
    }

    temp.swap(data);
    ThrowAssert(data.capacity() == capacity); // no growths took place
    m_last_capacity = capacity;
  }

  void resize_and_order_by_index(unsigned capacity = 0u)
  {
#ifdef STK_MESH_ANALYZE_DYN_CONN
    if (capacity != 0u) {
      ++profile_data().m_num_growths;
    }

    if (m_targets.capacity() > profile_data().m_max_capacity) {
      profile_data().m_max_capacity = m_targets.capacity();
      profile_data().m_total_unused_memory = m_targets.capacity() - m_total_connectivities;
      profile_data().m_unused_capacity = m_targets.capacity() - m_targets.size();
      profile_data().m_total_num_conn = m_total_connectivities;
      const size_t total_unused_active = m_targets.size() - m_total_connectivities;
      size_t total_unused_chunk_capacity = 0;
      for (size_t i = 0, e = m_num_connectivities.size(); i < e; ++i) {
        total_unused_chunk_capacity += num_chunks(m_num_connectivities[i]) * chunk_size - m_num_connectivities[i];
      }
      profile_data().m_abandoned_space = total_unused_active - total_unused_chunk_capacity;
      profile_data().m_unused_chunk_capacity = total_unused_chunk_capacity;
    }
#endif

    //compute needed capacity
    if (capacity == 0u) {
      for( size_t i=0, e=m_indices.size(); i<e; ++i) {
        capacity += num_chunks(m_num_connectivities[i]) * chunk_size;
      }
    }

    //move permutation
    if (has_permutation()) {
      resize_and_order_by_index_helper(m_permutations, capacity);
    }

    //move ordinal
    resize_and_order_by_index_helper(m_ordinals, capacity);

    //move target_by_key
    resize_and_order_by_index_helper(m_targets, capacity, true /*update index*/);
  }

  // The default capacity_ratio enables the most of the speed of the simple doubling
  // policy, while preventing the O(N^2) memory growth that would occur for
  // pathologically ordered add_connectivity calls under that policy.
  unsigned compute_new_connectivity_capacity(unsigned minimum_size, unsigned capacity_ratio = 8)
  {
    const unsigned old_capacity = m_targets.capacity();
    static const unsigned careful_threshold = 2048;

    if (old_capacity < careful_threshold)
    {
      unsigned new_capacity = old_capacity > 0 ? 2 * old_capacity : 8*chunk_size;
      while (new_capacity < minimum_size) {
        new_capacity *= 2;
      }
      return new_capacity;
    }

    // The old capacity is at or above the threshold for being careful about growing
    // the connectivity representation(and memory footprint).  Only grow the capacity
    // if compressing the representation will not yield sufficient unused capacity.
    if (capacity_ratio * m_total_connectivities > old_capacity)
    {
      return 2 * old_capacity;
    }
    else
    {
      return old_capacity;
    }
  }

  void add_connectivity_helper(unsigned bucket_ordinal, unsigned num_to_add=1)
  {
    const unsigned chunks_needed_by_entity = num_chunks(m_num_connectivities[bucket_ordinal]+num_to_add);
    const unsigned chunks_used_by_entity   = num_chunks(m_num_connectivities[bucket_ordinal]);

    if (chunks_needed_by_entity == chunks_used_by_entity)
    {
      m_total_connectivities               += num_to_add;
      m_num_connectivities[bucket_ordinal] += num_to_add;
      return;
    }

    const unsigned chunks_available = num_chunks(m_targets.capacity() - m_targets.size());

    if (chunks_available < chunks_needed_by_entity)
    {
      const unsigned new_capacity = compute_new_connectivity_capacity(m_targets.size() + chunks_needed_by_entity * chunk_size);
      resize_and_order_by_index(new_capacity);
    }

    const bool last_entity_by_index = (chunks_used_by_entity > 0) &&
      (m_indices[bucket_ordinal] + chunks_used_by_entity*chunk_size == m_targets.size());
    Entity invalid = {Entity::InvalidEntity};

    //copy to end
    if (!last_entity_by_index)
    {
#ifdef STK_MESH_ANALYZE_DYN_CONN
      if (chunks_used_by_entity > 0) {
        ++profile_data().m_num_entity_relocations;
      }
#endif

      uint32_t new_index = static_cast<uint32_t>(m_targets.size());

      m_targets.insert(m_targets.end(), chunks_needed_by_entity*chunk_size, invalid);
      std::copy(begin(bucket_ordinal), end(bucket_ordinal), m_targets.begin() + new_index);

      m_ordinals.insert(m_ordinals.end(), chunks_needed_by_entity*chunk_size, INVALID_CONNECTIVITY_ORDINAL);
      std::copy(begin_ordinals(bucket_ordinal), end_ordinals(bucket_ordinal), m_ordinals.begin() + new_index);

      if (has_permutation()) {
        m_permutations.insert(m_permutations.end(), chunks_needed_by_entity*chunk_size, INVALID_PERMUTATION);
        std::copy(begin_permutations(bucket_ordinal), end_permutations(bucket_ordinal), m_permutations.begin() + new_index);
      }

      m_indices[bucket_ordinal] = new_index;
    }
    //add new chunk to end
    else {
      const unsigned extra_chunks_needed = chunks_needed_by_entity - chunks_used_by_entity;
      m_targets.insert(m_targets.end(), extra_chunks_needed*chunk_size, invalid);
      m_ordinals.insert(m_ordinals.end(), extra_chunks_needed*chunk_size, INVALID_CONNECTIVITY_ORDINAL);
      if (has_permutation()) {
        m_permutations.insert(m_permutations.end(), extra_chunks_needed*chunk_size, INVALID_PERMUTATION);
      }
    }

    m_total_connectivities               += num_to_add;
    m_num_connectivities[bucket_ordinal] += num_to_add;
  }

  template <typename ConnectivityComparator>
  bool add_helper(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal, Permutation permutation,
                  const ConnectivityComparator &compare)
  {
#ifndef NDEBUG
    // TODO - If downward conn, check to's rank and topology
#endif
    bool rv = true;

    add_connectivity_helper(bucket_ordinal);

    const uint32_t begin_index = m_indices[bucket_ordinal] + m_num_connectivities[bucket_ordinal] - 1;

    if (m_num_connectivities[bucket_ordinal] == 1) {
      m_targets[begin_index] = to;
      m_ordinals[begin_index] = ordinal;
      if (has_permutation()) {
        m_permutations[begin_index] = permutation;
      }
      return true;
    }

    for (uint32_t i = begin_index, e = m_indices[bucket_ordinal]; i > e; --i)
    {
      //slide up
      if ( compare(to, ordinal, m_targets[i-1], m_ordinals[i-1u]) ) {
        m_targets[i] = m_targets[i-1u];
        m_ordinals[i] = m_ordinals[i-1u];
        if (has_permutation()) {
          m_permutations[i] = m_permutations[i-1u];
        }
        //insert if on last iteration
        if ((i-1)==e) {
          m_targets[i-1u] = to;
          m_ordinals[i-1u] = ordinal;
          if (has_permutation()) {
            m_permutations[i-1u] = permutation;
          }
        }
      }
      //insert
      else if ( compare(m_targets[i-1], m_ordinals[i-1u], to, ordinal) ) {
        m_targets[i] = to;
        m_ordinals[i] = ordinal;
        if (has_permutation()) {
          m_permutations[i] = permutation;
        }
        break;
      }
      //duplicate -- insert new and remove the original
      else
      {
        m_targets[i] = to;
        m_ordinals[i] = ordinal;
        if (has_permutation()) {
          ThrowAssert(permutation == m_permutations[i-1u]);
          m_permutations[i] = permutation;
        }
        remove_connectivity(bucket_ordinal, to, ordinal);
        rv = false;
        break;
      }
    }

    invariant_check_helper(bucket_ordinal);

    return rv;
  }

  void invariant_check_helper(unsigned bucket_ordinal) const
  {
  #ifdef STK_INVARIANCE_CHECK
    const Entity* keys_begin = begin(bucket_ordinal);
    const Entity* keys_end   = end(bucket_ordinal);
    const ConnectivityOrdinal* ordinals_begin = begin_ordinals(bucket_ordinal);
    const ConnectivityOrdinal* ordinals_end   = end_ordinals(bucket_ordinal);
    const Permutation* permutations_begin = begin_permutations(bucket_ordinal);
    const Permutation* permutations_end   = end_permutations(bucket_ordinal);

    ThrowAssertMsg(m_last_capacity == m_targets.capacity(), "Expected " << m_last_capacity << " found " << m_targets.capacity());
    ThrowAssertMsg((keys_end - keys_begin) == num_connectivity(bucket_ordinal),
                   "Expected data to be of size " << num_connectivity(bucket_ordinal) << ", " << bucket_ordinal << " has keys " << keys_end - keys_begin);

    ThrowAssertMsg(keys_end - keys_begin == ordinals_end - ordinals_begin,
                   "Num keys, " << keys_end - keys_begin << ", does not match num ordinals, " << ordinals_end - ordinals_begin);
    if (has_permutation()) {
      ThrowAssertMsg(keys_end - keys_begin == permutations_end - permutations_begin,
                     "Num keys, " << keys_end - keys_begin << ", does not match num permutations, " << permutations_end - permutations_begin);
    }
    else {
      ThrowAssertMsg(permutations_end - permutations_begin == 0,
                     "Expected 0 permutations for node connectivity, found: " << permutations_end - permutations_begin);
    }

    const Entity*               kitr = keys_begin;
    const ConnectivityOrdinal*  oitr = ordinals_begin;
    for ( ; kitr != keys_end; ++kitr, ++ oitr) {
      ThrowAssertMsg(*kitr != Entity(),
                     "Should not have invalid connectivity for dynamic connectivity");
      // TODO
      //entity_rank to_rank  = topology_rank(kitr->topology(), m_spatial_dimension);
      //ThrowAssertMsg(to_rank() == ToEntityRank::value,
      //               (debug_message() << "Found connectivity to wrong rank " << to_rank << ", expected " << entity_rank::create(ToEntityRank::value)));
      if (kitr + 1 != keys_end) {
        if (m_direction == Higher) { // back rel
          if (target_rank <= stk::topology::ELEMENT_RANK) {
            ThrowAssertMsg(HigherConnectivityCompare()(*kitr, *oitr, *(kitr + 1), *(oitr + 1)),
                           "Connectivity out of order for bucket_ordinal " << bucket_ordinal << "; data at " << kitr - keys_begin <<
                           "\nis (ordinal:" << *oitr << ", local_offset:" << kitr->local_offset() << ")" <<
                           ",\ndata at next slot is (ordinal:" << *(oitr + 1) << ", local_offset:" << (kitr + 1)->local_offset() << ")" <<
                           "\nConnectivity id is: " << m_id);
          }
          else {
            ThrowAssertMsg(m_rank_sensitive_higher_connectivity_cmp(*kitr, *oitr, *(kitr + 1), *(oitr + 1)),
                           "Connectivity out of order for bucket_ordinal " << bucket_ordinal << "; data at " << kitr - keys_begin <<
                           "\nis (" << *oitr << ", " << kitr->local_offset() << ")" <<
                           ",\ndata at next slot is (" << *(oitr + 1) << ", " << (kitr + 1)->local_offset() << ")" <<
                           "\nConnectivity id is: " << m_id);
          }
        }
        else {
          if (target_rank <= stk::topology::ELEMENT_RANK) {
            ThrowAssertMsg(LowerConnectivityCompare()(*kitr, *oitr, *(kitr + 1), *(oitr + 1)),
                           "Connectivity out of order for bucket_ordinal " << bucket_ordinal << "; data at " << kitr - keys_begin <<
                           "\nis (" << *oitr << ", " << kitr->local_offset() << ")" <<
                           ",\ndata at next slot is (" << *(oitr + 1) << ", " << (kitr + 1)->local_offset() << ")" <<
                           "\nConnectivity id is: " << m_id);
          }
          else {
            ThrowAssertMsg(m_rank_sensitive_lower_connectivity_cmp(*kitr, *oitr, *(kitr + 1), *(oitr + 1)),
                           "Connectivity out of order for bucket_ordinal " << bucket_ordinal << "; data at " << kitr - keys_begin <<
                           "\nis (" << *oitr << ", " << kitr->local_offset() << ")" <<
                           ",\ndata at next slot is (" << *(oitr + 1) << ", " << (kitr + 1)->local_offset() << ")" <<
                           "\nConnectivity id is: " << m_id);
          }
        }
      }
      // TODO - Anything else we can check here?
    }

    invariant_check_helper();
  #endif

  }

  void invariant_check_helper() const
  {

  #ifdef STK_INVARIANCE_CHECK
    if (!m_active) {
      ThrowAssertMsg(m_num_connectivities.size() == 0, "Expect empty data if inactive");
    }

    ThrowAssertMsg(m_last_capacity == m_targets.capacity(), "Expected " << m_last_capacity << " found " << m_targets.capacity());

    unsigned connectivities_sum = 0;
    size_t lim = m_num_connectivities.size();
    for (size_t i = 0; i < lim; ++i)
    {
      connectivities_sum += m_num_connectivities[i];
    }
    ThrowAssertMsg(m_total_connectivities == connectivities_sum,
                   "Expected m_total_connectivities == " << connectivities_sum << ", found " << m_total_connectivities);

    ThrowAssertMsg(m_num_connectivities.size() == m_indices.size(),
                   "Expected m_num_connectivities to be of size " << m_indices.size() << ", found " << m_num_connectivities.size());

    ThrowAssertMsg(m_targets.size() == m_ordinals.size(),
                   "Total size of keys " << m_targets.size() << " does not match size of ordinals " << m_ordinals.size());

    if (has_permutation()) {
      ThrowAssertMsg(m_permutations.size() == m_targets.size(),
                     "Total size of permutationss " << m_permutations.size() << " does not match size of keys " << m_targets.size());
    }
    else {
      ThrowAssertMsg(m_permutations.empty(), "Permutations should be empty for nodal connectivity");
    }

    for (size_t o = 1, e = m_indices.size(); o < e; ++o) {
      const size_t curr_index = m_indices[o];
      ThrowAssertMsg(curr_index <= m_targets.size(),
                     "Index is wrong, " << curr_index << " is beyond max " << m_targets.size());
      if (!m_needs_shrink_to_fit) {
        const size_t index_diff     = curr_index - m_indices[o-1];
        const size_t prior_num_conn = num_chunks(num_connectivity(o-1)) * chunk_size;
        ThrowAssertMsg(prior_num_conn == index_diff,
                       "For offset " << o << ", num_connectivity/index mismatch, index_diff is " << index_diff << ", num conn is " << prior_num_conn);
      }
    }
  #endif

  }

  // Call after modification end
  template <typename BulkData>
  void invariant_check_helper(BulkData* mesh = NULL) const
  {

  #ifdef STK_INVARIANCE_CHECK
    invariant_check_helper();

    ThrowAssert(!m_active || !m_needs_shrink_to_fit);

    for (size_t o = 0, e = m_indices.size(); o < e; ++o) {
      invariant_check_helper(o);
      if (o > 0) {
        const size_t curr_index     = m_indices[o];
        const size_t index_diff     = curr_index - m_indices[o-1];
        const size_t prior_num_conn = num_connectivity(o-1);
        ThrowAssertMsg(prior_num_conn == index_diff,
                       "For offset " << o << ", num_connectivity/index mismatch, index_diff is " << index_diff << ", num conn is " << prior_num_conn);
      }
    }

    // Check that connectivity is in-sync
    ThrowAssertMsg(m_targets.size() == m_ordinals.size(),
                   "Total size of partition indices " << m_targets.size() << " does not match size of ordinals " << m_ordinals.size());
  #endif

  }

  // Illegal
  BucketConnectivity(const SelfType&);
  SelfType& operator=(const SelfType&);

  // MEMBERS

  EntityRank m_from_rank;
  connectivity_direction m_direction;

  bool m_active; // In many cases, uses will not make use of dynamic connectivity, so don't even waste the memory unless it looks like they want it
  bool m_needs_shrink_to_fit; // True if this object potentially has partially full vectors or out-of-order entities
  unsigned  m_num_inactive;

  // meta data
  UInt32Vector m_indices;  // Common index into vectors below that stores where connectivity starts for a partition_offset (entity).
  UInt16Vector m_num_connectivities;
  unsigned     m_total_connectivities;

  // connectivity data
  EntityVector              m_targets;
  ConnectivityOrdinalVector m_ordinals;
  PermutationVector         m_permutations;

  BulkData * m_bulk_data;
  int        m_id;

  impl::HigherConnectivityRankSensitiveCompare<BulkData>  m_rank_sensitive_higher_connectivity_cmp;
  impl::LowerConnectivitityRankSensitiveCompare<BulkData>  m_rank_sensitive_lower_connectivity_cmp;

  size_t m_last_capacity;

  size_t m_data_idx;
};

}

//
// BucketConnectivity
//

template <EntityRank TargetRank>
template <typename BulkData> // hack to get around dependency
inline
void impl::BucketConnectivity<TargetRank, FIXED_CONNECTIVITY>::end_modification(BulkData* mesh)
{
  //TODO: If bucket is blocked, no longer need to shrink to fit!

  if (m_targets.size() < m_targets.capacity()) {

    {
      EntityVector temp(m_targets.begin(), m_targets.end());
      m_targets.swap(temp);
    }

    {
      PermutationVector temp(m_permutations.begin(), m_permutations.end());
      m_permutations.swap(temp);
    }
  }

  invariant_check_helper(mesh);
}

template <EntityRank TargetRank>
template <typename BulkData>
inline
void impl::BucketConnectivity<TargetRank, DYNAMIC_CONNECTIVITY>::end_modification(BulkData* mesh)
{
  if (m_active && m_needs_shrink_to_fit) {
    resize_and_order_by_index();

    {
      UInt32Vector temp(m_indices.begin(), m_indices.end());
      m_indices.swap(temp);
    }

    {
      UInt16Vector temp(m_num_connectivities.begin(), m_num_connectivities.end());
      m_num_connectivities.swap(temp);
    }

    m_needs_shrink_to_fit = false;
  }

  invariant_check_helper(mesh);
}

template <typename BulkData>
inline
bool impl::LowerConnectivitityRankSensitiveCompare<BulkData>::operator()(Entity first_entity, ConnectivityOrdinal first_ordinal,
                                                                         Entity second_entity, ConnectivityOrdinal second_ordinal) const
{
  const EntityRank first_rank = m_mesh.entity_rank(first_entity);
  const EntityRank second_rank = m_mesh.entity_rank(second_entity);

  return (first_rank < second_rank)
         || ((first_rank == second_rank) && (first_ordinal < second_ordinal));
}

template <typename BulkData>
inline
bool impl::HigherConnectivityRankSensitiveCompare<BulkData>::operator()(Entity first_entity, ConnectivityOrdinal first_ordinal, Entity second_entity, ConnectivityOrdinal second_ordinal) const
{
  const EntityRank first_rank = m_mesh.entity_rank(first_entity);
  const EntityRank second_rank = m_mesh.entity_rank(second_entity);

  if (first_rank < second_rank) {
    return true;
  }
  if (first_rank > second_rank) {
    return false;
  }
  // Needs to match LessRelation in BulkData.hpp
  return std::make_pair(first_ordinal,  first_entity.is_local_offset_valid() ?  first_entity.local_offset()  : Entity::MaxEntity) <
         std::make_pair(second_ordinal, second_entity.is_local_offset_valid() ? second_entity.local_offset() : Entity::MaxEntity);
}





}} //namespace stk::mesh::impl

#endif
