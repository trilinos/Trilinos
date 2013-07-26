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

  typedef std::vector<Entity,              tracking_allocator<Entity, BucketRelationTag> >              EntityVector;
  typedef std::vector<ConnectivityOrdinal, tracking_allocator<ConnectivityOrdinal, BucketRelationTag> > ConnectivityOrdinalVector;
  typedef std::vector<Permutation,         tracking_allocator<Permutation, BucketRelationTag> >         PermutationVector;

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

  bool add_connectivity(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal, Permutation permutation)
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

  void swap(unsigned my_ordinal, SelfType& to, unsigned to_ordinal)
  {
    ThrowAssertMsg(m_num_connectivity == to.m_num_connectivity,
                   "Incompatible connectivitiess, src connectivity expects " << m_num_connectivity <<
                   " per entity, target connectivity expects " << to.m_num_connectivity);
    impl::check_bucket_ordinal(my_ordinal, this);
    impl::check_bucket_ordinal(to_ordinal, &to);

    const unsigned first_index  = m_num_connectivity * my_ordinal;
    const unsigned second_index = m_num_connectivity * to_ordinal;

    std::swap_ranges( m_targets.begin()+first_index
                     ,m_targets.begin()+first_index+m_num_connectivity
                     ,to.m_targets.begin()+second_index
                    );

    if (has_permutation()) {
      std::swap_ranges( m_permutations.begin()+first_index
                       ,m_permutations.begin()+first_index+m_num_connectivity
                       ,to.m_permutations.begin()+second_index
                      );
    }

    invariant_check_helper(my_ordinal);
    to.invariant_check_helper(to_ordinal);
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

  void move_entity(SelfType& to)
  {
    ThrowAssertMsg(size() > 0, "Cannot move, connectivity is empty");

    unsigned to_offset = to.m_targets.size();
    to.add_entity(); // make room for new entity

    unsigned from_offset = m_targets.size() - m_num_connectivity;

    std::copy(m_targets.begin()+from_offset, m_targets.end(), to.m_targets.begin() + to_offset);

    if (has_permutation()) {
      std::copy(m_permutations.begin()+from_offset, m_permutations.end(), to.m_permutations.begin() + to_offset);
    }

    remove_entity();

    invariant_check_helper();
    to.invariant_check_helper();
  }

  void move_to_fixed(SelfType& to)
  { ThrowAssert(false); }

  void move_to_fixed(OtherType& to)
  { ThrowAssert(false); }

  bool has_permutation() const
  {
    const static bool rv = TargetRank != stk::topology::NODE_RANK;
    return rv;
  }

private:

  void invariant_check_helper(unsigned bucket_ordinal) const
  {
#ifndef NDEBUG
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
#ifndef NDEBUG
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
#ifndef NDEBUG
    invariant_check_helper();

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

template<EntityRank TargetRank >
class BucketConnectivity<TargetRank, DYNAMIC_CONNECTIVITY>
{
  enum connectivity_direction { Lower=0,Higher=1,Adjacent=2 };

public:
  typedef BucketConnectivity<TargetRank, DYNAMIC_CONNECTIVITY> SelfType;
  typedef BucketConnectivity<TargetRank, FIXED_CONNECTIVITY> OtherType;

  static const EntityRank target_rank = TargetRank;
  static const ConnectivityType connectivity_type = DYNAMIC_CONNECTIVITY;

  typedef std::vector<Entity,              tracking_allocator<Entity, DynamicBucketRelationTag> >              EntityVector;
  typedef std::vector<ConnectivityOrdinal, tracking_allocator<ConnectivityOrdinal, DynamicBucketRelationTag> > ConnectivityOrdinalVector;
  typedef std::vector<Permutation,         tracking_allocator<Permutation, DynamicBucketRelationTag> >         PermutationVector;
  typedef std::vector<uint32_t,            tracking_allocator<uint32_t, DynamicBucketRelationTag> >            UInt32Vector;
  typedef std::vector<uint16_t,            tracking_allocator<uint16_t, DynamicBucketRelationTag> >            UInt16Vector;

  static const unsigned chunk_size = 1u;

  BucketConnectivity() //default constructed BucketConnectivity implies connectivity is not used
    : m_from_rank(InvalidEntityRank)
    , m_direction(Lower)
    , m_active(false)
    , m_num_inactive(0)
    , m_indices()
    , m_num_connectivities()
    , m_targets()
    , m_ordinals()
    , m_permutations()
    , m_bulk_data(0)
  {}

  BucketConnectivity(EntityRank from_rank, BulkData *bulk_data)
    : m_from_rank(from_rank)
    , m_direction( (m_from_rank > TargetRank) ? Lower : ((m_from_rank == TargetRank) ? Adjacent : Higher))
    , m_active(false)
    , m_num_inactive(0)
    , m_indices()
    , m_num_connectivities()
    , m_targets()
    , m_ordinals()
    , m_permutations()
    , m_bulk_data(bulk_data)
    , m_rank_sensitive_higher_connectivity_cmp(*m_bulk_data)
    , m_rank_sensitive_lower_connectivity_cmp(* m_bulk_data)
  {}

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

  bool add_connectivity(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal, Permutation permutation)
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);

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
        break;
      }
    }

    //slide memory down
    if (found_idx != ~0u) {
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

  void swap(unsigned my_ordinal, SelfType& to, unsigned to_ordinal)
  {
    impl::check_bucket_ordinal(my_ordinal, this);
    impl::check_bucket_ordinal(to_ordinal, &to);

    // Manage activation state
    if (!m_active && !to.m_active) {
      return;
    }
    else if (!m_active && to.m_active) {
      activate();
    }
    else if (m_active && !to.m_active) {
      to.activate();
    }

    // Swap data
    if (&to == this) {
      // easy
      std::swap( m_indices[my_ordinal], to.m_indices[to_ordinal] );
      std::swap( m_num_connectivities[my_ordinal], to.m_num_connectivities[to_ordinal]);
    }
    else if (m_num_connectivities[my_ordinal] == to.m_num_connectivities[to_ordinal]) {
      // harder, but still relatively easy
      const unsigned first_index  = m_indices[my_ordinal];
      const unsigned second_index = to.m_indices[to_ordinal];

      std::swap_ranges( m_targets.begin() + first_index
                        ,m_targets.begin() + first_index + m_num_connectivities[my_ordinal]
                        ,to.m_targets.begin() + second_index
                        );

      std::swap_ranges( m_ordinals.begin() + first_index
                        ,m_ordinals.begin() + first_index + m_num_connectivities[my_ordinal]
                        ,to.m_ordinals.begin() + second_index
                        );

      if (has_permutation()) {
        std::swap_ranges( m_permutations.begin() + first_index
                          ,m_permutations.begin() + first_index + to.m_num_connectivities[to_ordinal]
                          ,to.m_permutations.begin() + second_index
                          );
      }
    }
    else {
      // very hard, clearly not the intended use of this data structure, would row-storage be better?
      const unsigned first_num    = m_num_connectivities[my_ordinal];
      const unsigned second_num   = to.m_num_connectivities[to_ordinal];

      if (first_num < second_num) {
        diff_length_connectivity_swap(my_ordinal, *this, to_ordinal, to);
      }
      else {
        diff_length_connectivity_swap(to_ordinal, to, my_ordinal, *this);
      }

    }

    invariant_check_helper();
    to.invariant_check_helper();
  }

  void begin_modification()
  {}

  template <typename BulkData>
  void end_modification(BulkData* mesh = NULL);

  void add_entity()
  {
    if (m_active) {
      m_indices.push_back(0); // it's OK for entities with no connectivity to have wrong indices
      m_num_connectivities.push_back(0);
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
      m_num_connectivities.pop_back();
    }
    else {
      --m_num_inactive;
    }

    invariant_check_helper();
  }

  void move_entity(SelfType & to)
  {
    ThrowAssertMsg(size() > 0, "Cannot move, connectivity is empty");

    const unsigned to_num_connectivity_before_move = to.m_targets.size();
    const unsigned my_start_index                  = m_active ? m_indices.back() : 0;
    const unsigned num_connectivity_to_move        = m_active ? m_num_connectivities.back() : 0;

    if (!m_active && !to.m_active) {
      remove_entity();
      to.add_entity();
      return;
    }
    if (m_active && !to.m_active) {
      to.activate();
    }

    //move over index and num_connectivity
    to.m_indices.push_back(to_num_connectivity_before_move);

    to.m_num_connectivities.push_back(num_connectivity_to_move);

    //move over target, ordinal, and permutation
    to.m_targets.insert( to.m_targets.end()
                                 ,m_targets.begin() + my_start_index
                                 ,m_targets.begin() + my_start_index + num_connectivity_to_move
                               );

    to.m_ordinals.insert( to.m_ordinals.end()
                          ,m_ordinals.begin() + my_start_index
                          ,m_ordinals.begin() + my_start_index + num_connectivity_to_move
                        );


    if (has_permutation()) {
      to.m_permutations.insert( to.m_permutations.end()
                               ,m_permutations.begin() + my_start_index
                               ,m_permutations.begin() + my_start_index + num_connectivity_to_move
                               );
    }

    remove_entity();

    invariant_check_helper();
    to.invariant_check_helper();
  }

  void move_to_fixed(OtherType& to)
  {
    const unsigned num_conn_to_move = m_active ? m_num_connectivities.back() : 0;

    ThrowAssert(OtherType::connectivity_type == FIXED_CONNECTIVITY);
    ThrowAssertMsg(size() > 0, "Cannot move, connectivity is empty");
    ThrowAssertMsg(num_conn_to_move <= to.num_connectivity(666 /*any unsigned, doesn't matter*/), "Incompatible");

    const unsigned to_offset = to.m_targets.size();
    to.add_entity(); // make room for new entity

    const unsigned from_offset = m_active ? m_indices.back() : 0;

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

    remove_entity();

    invariant_check_helper();
    to.invariant_check_helper();
  }

  void move_to_fixed(SelfType& to)
  { ThrowAssert(false); }

  bool has_permutation() const
  { return TargetRank != stk::topology::NODE_RANK && m_from_rank != stk::topology::NODE_RANK; }

private:

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

  void diff_length_connectivity_swap(unsigned ordinal_for_conn_getting_more, SelfType& conn_getting_more,
                                     unsigned ordinal_for_conn_getting_less, SelfType& conn_getting_less)
  {
    const unsigned getting_more_index  = conn_getting_more.m_indices[ordinal_for_conn_getting_more];
    const unsigned getting_less_index  = conn_getting_less.m_indices[ordinal_for_conn_getting_less];
    const unsigned getting_more_num    = conn_getting_more.m_num_connectivities[ordinal_for_conn_getting_more];
    const unsigned getting_less_num    = conn_getting_less.m_num_connectivities[ordinal_for_conn_getting_less];

    ThrowAssert(getting_more_num < getting_less_num);

    std::copy( conn_getting_less.m_targets.begin() + getting_less_index,
               conn_getting_less.m_targets.begin() + getting_less_index + getting_less_num,
               std::back_inserter(conn_getting_more.m_targets) );
    std::copy( conn_getting_more.m_targets.begin() + getting_more_index,
               conn_getting_more.m_targets.begin() + getting_more_index + getting_more_num,
               conn_getting_less.m_targets.begin() + getting_less_index); // no back inserting needed, just overwrite

    std::copy( conn_getting_less.m_ordinals.begin() + getting_less_index,
               conn_getting_less.m_ordinals.begin() + getting_less_index + getting_less_num,
               std::back_inserter(conn_getting_more.m_ordinals) );
    std::copy( conn_getting_more.m_ordinals.begin() + getting_more_index,
               conn_getting_more.m_ordinals.begin() + getting_more_index + getting_more_num,
               conn_getting_less.m_ordinals.begin() + getting_less_index );

    if (has_permutation()) {
      std::copy( conn_getting_less.m_permutations.begin() + getting_less_index,
                 conn_getting_less.m_permutations.begin() + getting_less_index + getting_less_num,
                 std::back_inserter(conn_getting_more.m_permutations) );
      std::copy( conn_getting_more.m_permutations.begin() + getting_more_index,
                 conn_getting_more.m_permutations.begin() + getting_more_index + getting_more_num,
                 conn_getting_less.m_permutations.begin() + getting_less_index );
    }

    conn_getting_more.m_indices[ordinal_for_conn_getting_more] = conn_getting_more.m_targets.size() - getting_less_num;
    // don't need to change the index of the conn that got less

    std::swap( conn_getting_more.m_num_connectivities[ordinal_for_conn_getting_more],
               conn_getting_less.m_num_connectivities[ordinal_for_conn_getting_less]);
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
  }

  void resize_and_order_by_index(unsigned capacity = 0u)
  {
    //compute needed capacity
    if (capacity == 0u) {
      for( size_t i=0, e=m_indices.size(); i<e; ++i) {
        capacity += num_chunks(m_num_connectivities[i]);
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

  void add_connectivity_helper(unsigned bucket_ordinal)
  {
    const unsigned chunks_needed_by_entity = num_chunks(m_num_connectivities[bucket_ordinal]+1);
    const unsigned chunks_used_by_entity   = num_chunks(m_num_connectivities[bucket_ordinal]);

    if (chunks_needed_by_entity == chunks_used_by_entity)
    {
      ++m_num_connectivities[bucket_ordinal];
      return;
    }

    const unsigned chunks_available = num_chunks(m_targets.capacity() - m_targets.size());

    if (chunks_available < chunks_needed_by_entity)
    {
      // Important line, defines how capacity grows. We do doublings for now.
      const unsigned new_capacity = m_targets.capacity() > 0 ? 2*m_targets.capacity() : 8*chunk_size;
      resize_and_order_by_index(new_capacity);
    }

    const bool last_entity_by_index = (chunks_used_by_entity > 0) &&
      (m_indices[bucket_ordinal] + chunks_used_by_entity*chunk_size == m_targets.size());
    Entity invalid = {Entity::InvalidEntity};

    //copy to end
    if (!last_entity_by_index)
    {
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
      m_targets.insert(m_targets.end(), chunk_size, invalid);
      m_ordinals.insert(m_ordinals.end(), chunk_size, INVALID_CONNECTIVITY_ORDINAL);
      if (has_permutation()) {
        m_permutations.insert(m_permutations.end(), chunk_size, INVALID_PERMUTATION);
      }
    }

    ++m_num_connectivities[bucket_ordinal];
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
#ifndef NDEBUG
    const Entity* keys_begin = begin(bucket_ordinal);
    const Entity* keys_end   = end(bucket_ordinal);
    const ConnectivityOrdinal* ordinals_begin = begin_ordinals(bucket_ordinal);
    const ConnectivityOrdinal* ordinals_end   = end_ordinals(bucket_ordinal);
    const Permutation* permutations_begin = begin_permutations(bucket_ordinal);
    const Permutation* permutations_end   = end_permutations(bucket_ordinal);

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
                           "Connectivity out of order; data at " << kitr - keys_begin <<
                           "\nis (" << *oitr << ", " << kitr->local_offset() << ") rank=?" <<
                           ",\ndata at next slot is (" << *(oitr + 1) << ", " << (kitr + 1)->local_offset() << ") rank=?");
          }
          else {
          ThrowAssertMsg(m_rank_sensitive_higher_connectivity_cmp(*kitr, *oitr, *(kitr + 1), *(oitr + 1)),
                         "Connectivity out of order; data at " << kitr - keys_begin <<
                         "\nis (" << *oitr << ", " << kitr->local_offset() << ") rank=?" <<
                         ",\ndata at next slot is (" << *(oitr + 1) << ", " << (kitr + 1)->local_offset() << ") rank=?");
          }
        }
        else {
          if (target_rank <= stk::topology::ELEMENT_RANK) {
            ThrowAssertMsg(LowerConnectivityCompare()(*kitr, *oitr, *(kitr + 1), *(oitr + 1)),
                           "Connectivity out of order; data at " << kitr - keys_begin <<
                           "\nis (" << *oitr << ", " << kitr->local_offset() << ") rank=?" <<
                           ",\ndata at next slot is (" << *(oitr + 1) << ", " << (kitr + 1)->local_offset() << ") rank=?");
          }
          else {
            ThrowAssertMsg(m_rank_sensitive_lower_connectivity_cmp(*kitr, *oitr, *(kitr + 1), *(oitr + 1)),
                           "Connectivity out of order; data at " << kitr - keys_begin <<
                           "\nis (" << *oitr << ", " << kitr->local_offset() << "),\ndata at next slot is (" << *(oitr + 1) << ", " << (kitr + 1)->local_offset() << ")");
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
#ifndef NDEBUG
    if (!m_active) {
      ThrowAssertMsg(m_num_connectivities.size() == 0, "Expect empty data if inactive");
    }

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
    }
#endif
  }

  // Call after modification end
  template <typename BulkData>
  void invariant_check_helper(BulkData* mesh = NULL) const
  {
#ifndef NDEBUG
    invariant_check_helper();

    for (size_t o = 1, e = m_indices.size(); o < e; ++o) {
      const size_t curr_index     = m_indices[o];
      const size_t index_diff     = curr_index - m_indices[o-1];
      const size_t prior_num_conn = num_connectivity(o-1);
      ThrowAssertMsg(prior_num_conn == index_diff,
                     "For offset " << o << ", num_connectivity/index mismatch, index_diff is " << index_diff << ", num conn is " << prior_num_conn);
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
  unsigned  m_num_inactive;

  // meta data
  UInt32Vector m_indices;  // Common index into vectors below that stores where connectivity starts for a partition_offset (entity).
  UInt16Vector m_num_connectivities;

  // connectivity data
  EntityVector              m_targets;
  ConnectivityOrdinalVector m_ordinals;
  PermutationVector         m_permutations;

  BulkData                       * m_bulk_data;
  impl::HigherConnectivityRankSensitiveCompare<BulkData>  m_rank_sensitive_higher_connectivity_cmp;
  impl::LowerConnectivitityRankSensitiveCompare<BulkData>  m_rank_sensitive_lower_connectivity_cmp;
};

}}} //namespace stk::mesh::impl
