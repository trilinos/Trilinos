#ifndef SAMBA_SAMBA_MESH_FIXED_PARTITION_CONNECTIVITY_TCC
#define SAMBA_SAMBA_MESH_FIXED_PARTITION_CONNECTIVITY_TCC

#include <samba/types.hpp>
#include <samba/rank_index.hpp>

namespace samba { namespace detail {

template< typename ToEntityRank >
class partition_connectivity<ToEntityRank, connectivity_kind::fixed_type>
{
private:
  template <typename IndexType>
  struct index_helper
  {
    static index_helper<IndexType> create(partition_offset arg_value)
    { index_helper tmp = {arg_value()}; return tmp; }

    uint32_t operator ()() const { return m_value; }

    uint32_t m_value;
  };

  typedef typename entity_rank_to_index_type<ToEntityRank>::type rank_index;

  typedef index_helper<entity_key> entity_key_helper;
  typedef index_helper<partition_index> partition_index_helper;
  typedef index_helper<rank_index> rank_index_helper;
  typedef index_helper<connectivity_ordinal> connectivity_ordinal_helper;
  typedef index_helper<connectivity_orientation> connectivity_orientation_helper;

public:

  typedef partition_connectivity<ToEntityRank, connectivity_kind::fixed_type> self_type;

  partition_connectivity() //default constructed partition_connectivity implies connectivity is not used
    : m_spatial_dimension(spatial_dimension::invalid())
    , m_from_topology(entity_topology::invalid())
    , m_num_connectivity(0u)
    , m_ordinal()
    , m_orientation()
    , m_target_by_key()
    , m_target_by_partition_index()
    , m_target_by_rank_index()
  {}

  partition_connectivity(entity_topology from_topology, spatial_dimension arg_spatial_dimension)
    : m_spatial_dimension(arg_spatial_dimension)
    , m_from_topology(from_topology)
    , m_num_connectivity( ToEntityRank::value == entity_rank::node_type::value ? num_nodes(from_topology)
        : (ToEntityRank::value == entity_rank::edge_type::value ? num_edges(from_topology)
          : (ToEntityRank::value == entity_rank::face_type::value ? num_faces(from_topology)
            : (ToEntityRank::value == entity_rank::element_type::value ? num_sides(from_topology) : 0u)
            )
          )
        )
    , m_ordinal(m_num_connectivity)
    , m_orientation()
    , m_target_by_key()
    , m_target_by_partition_index()
    , m_target_by_rank_index()
  {
    connectivity_ordinal ord;
    for ( ord=0; ord < m_num_connectivity; ++ord) m_ordinal[ord()] = ord;
  }


  // range
  template< typename Index >
  std::pair<const Index*, const Index *> range(partition_offset offset) const
  { return std::make_pair( begin<Index>(offset), end<Index>(offset)); }

  //begin
  template< typename Index >
  const Index * begin(partition_offset offset) const
  { return begin(index_helper<Index>::create(offset)); }

  //end
  template< typename Index >
  const Index * end(partition_offset offset) const
  { return begin<Index>(offset) + m_num_connectivity; }

  size_t num_connectivity(partition_offset /*offset*/) const
  { return m_num_connectivity; }

  size_t size() const
  { return m_target_by_key.size() / m_num_connectivity; }

  // Modification API
  void add_connectivity(partition_offset from, entity_key to, connectivity_ordinal ordinal, connectivity_orientation orientation)
  {
#ifndef NDEBUG
    BOOST_ASSERT_MSG(ordinal() < m_num_connectivity,
                     (debug_message() << "Ordinal " << ordinal << " exceeds topological limit: " << m_num_connectivity));
    if ( topology_rank(m_from_topology,m_spatial_dimension) > topology_rank(to.topology(),m_spatial_dimension) )
    {
      entity_topology expected_to_topology = connectivity_topology(m_from_topology, topology_rank(to.topology(),m_spatial_dimension), ordinal);
      BOOST_ASSERT_MSG( expected_to_topology == entity_topology::invalid() || expected_to_topology == to.topology(),
                        (debug_message() << "Expected connectivity to entity with topology " << expected_to_topology << ", found " << to.topology()));
    }
#endif

    size_t index = m_num_connectivity*from() + ordinal();

    m_target_by_key[index] = to;
    if (ToEntityRank::value != entity_rank::node_type::value) {
      m_orientation[index] = orientation;
    }

    invariant_check_helper(from);
  }

  void remove_connectivity(partition_offset from, entity_key to, connectivity_ordinal ordinal)
  {
    BOOST_ASSERT_MSG(ordinal() < m_num_connectivity,
                     (debug_message() << "Ordinal " << ordinal << " exceeds topological limit: " << m_num_connectivity));

    size_t index = m_num_connectivity*from() + ordinal();
    m_target_by_key[index] = entity_key::invalid();
    if (ToEntityRank::value != entity_rank::node_type::value) {
      m_orientation[index] = connectivity_orientation::invalid();
    }

    invariant_check_helper(from);
  }

  void swap(partition_offset first, partition_offset second)
  {
    size_t first_index = m_num_connectivity*first();
    size_t second_index = m_num_connectivity*second();

    std::swap_ranges( m_target_by_key.begin()+first_index
                     ,m_target_by_key.begin()+first_index+m_num_connectivity
                     ,m_target_by_key.begin()+second_index
                    );

    if (ToEntityRank::value != entity_rank::node_type::value) {
      std::swap_ranges( m_orientation.begin()+first_index
                       ,m_orientation.begin()+first_index+m_num_connectivity
                       ,m_orientation.begin()+second_index
                      );
    }

    invariant_check_helper(first);
    invariant_check_helper(second);
  }

  void begin_modification()
  {
    std::vector<partition_index> ptemp;
    m_target_by_partition_index.swap(ptemp);

    std::vector<rank_index> rtemp;
    m_target_by_rank_index.swap(rtemp);
  }

  void end_modification(mesh_impl * mesh)
  {
    {
      std::vector<entity_key> temp_target_by_key(m_target_by_key.begin(),m_target_by_key.end());
      m_target_by_key.swap(temp_target_by_key);
    }

    {
      std::vector<connectivity_orientation> temp_orientation(m_orientation.begin(),m_orientation.end());
      m_orientation.swap(temp_orientation);
    }

    // Populate partition_index and rank_index structures
    m_target_by_partition_index.resize(m_target_by_key.size(), partition_index::invalid());
    m_target_by_rank_index.resize(m_target_by_key.size(),      rank_index::invalid());
    for (size_t i = 0, end = m_target_by_key.size(); i < end; ++i) {
      const entity_key target_key = m_target_by_key[i];
      if (target_key != entity_key::invalid()) {
        m_target_by_partition_index[i] = mesh->convert(target_key);
        m_target_by_rank_index[i]      = mesh->convert<rank_index>(target_key);
      }
    }

    invariant_check_helper(mesh);
  }

  void add_entities(size_t how_many)
  {
    const size_t new_size = m_target_by_key.size() + how_many*m_num_connectivity;
    m_target_by_key.resize(new_size, entity_key::invalid());
    if (ToEntityRank::value != entity_rank::node_type::value) {
      m_orientation.resize(new_size, connectivity_orientation::invalid());
    }

    invariant_check_helper();
  }

  void remove_entities(size_t how_many)
  {
    const size_t new_size = m_target_by_key.size() - how_many*m_num_connectivity;
    m_target_by_key.resize(new_size);
    if (ToEntityRank::value != entity_rank::node_type::value) {
      m_orientation.resize(new_size);
    }

    invariant_check_helper();
  }

  void move_entities(self_type& to, size_t how_many)
  {
    size_t to_offset = to.m_target_by_key.size();
    to.add_entities(how_many);

    size_t from_offset = m_target_by_key.size() - how_many*m_num_connectivity;

    std::copy(m_target_by_key.begin()+from_offset, m_target_by_key.end(), to.m_target_by_key.begin() + to_offset);

    if (ToEntityRank::value != entity_rank::node_type::value) {
      std::copy(m_orientation.begin()+from_offset, m_orientation.end(), to.m_orientation.begin() + to_offset);
    }

    remove_entities(how_many);

    invariant_check_helper();
    to.invariant_check_helper();
  }

private:

  //begin entity_key
  const entity_key * begin(entity_key_helper offset) const
  { return &*m_target_by_key.begin() + m_num_connectivity*offset(); }

  //begin partition_index
  const partition_index * begin(partition_index_helper offset) const
  { return &*m_target_by_partition_index.begin() + m_num_connectivity*offset(); }

  //begin rank_index
  const rank_index * begin(rank_index_helper offset) const
  { return &*m_target_by_rank_index.begin() + m_num_connectivity*offset(); }

  //begin connectivity_ordinal
  const connectivity_ordinal * begin(connectivity_ordinal_helper offset) const
  { return &*m_ordinal.begin(); }

  //begin connectivity_orientation
  const connectivity_orientation * begin(connectivity_orientation_helper offset) const
  {
    if (ToEntityRank::value != entity_rank::node_type::value) {
      return &*m_orientation.begin() + m_num_connectivity*offset();
    }
    else {
      return NULL;
    }
  }

  void invariant_check_helper(partition_offset offset) const
  {
#ifndef NDEBUG
    const entity_key* keys_begin = begin<entity_key>(offset);
    const entity_key* keys_end   = end<entity_key>(offset);
    const connectivity_ordinal* ordinals_begin = begin<connectivity_ordinal>(offset);
    const connectivity_ordinal* ordinals_end   = end<connectivity_ordinal>(offset);
    const connectivity_orientation* orientations_begin = begin<connectivity_orientation>(offset);
    const connectivity_orientation* orientations_end   = end<connectivity_orientation>(offset);

    BOOST_ASSERT_MSG(static_cast<size_t>(keys_end - keys_begin) == m_num_connectivity,
                       (debug_message() << "Expected data to be of size " << m_num_connectivity << ", " << offset << " has keys " << keys_end - keys_begin));

    BOOST_ASSERT_MSG(keys_end - keys_begin == ordinals_end - ordinals_begin,
                     (debug_message() << "Num keys, " << keys_end - keys_begin << ", does not match num ordinals, " << ordinals_end - ordinals_begin));
    if (ToEntityRank::value != entity_rank::node_type::value) {
      BOOST_ASSERT_MSG(keys_end - keys_begin == orientations_end - orientations_begin,
                       (debug_message() << "Num keys, " << keys_end - keys_begin << ", does not match num orientations, " << orientations_end - orientations_begin));
    }
    else {
      BOOST_ASSERT_MSG(orientations_end - orientations_begin == 0,
                       (debug_message() << "Expected 0 orientations for node connectivity, found: " << orientations_end - orientations_begin));
    }

    const entity_key*               kitr = keys_begin;
    const connectivity_ordinal*     oitr = ordinals_begin;
    const connectivity_orientation* ritr = orientations_begin;
    for ( ; kitr != keys_end; ++kitr, ++ oitr) {
      if (*kitr != entity_key::invalid()) {
        BOOST_ASSERT_MSG(*oitr == kitr - keys_begin,
                         (debug_message() << "Found out-of-order connectivity at index " << kitr - keys_begin << ", its ordinal is " << *oitr));
        entity_rank to_rank  = topology_rank(kitr->topology(), m_spatial_dimension);
        BOOST_ASSERT_MSG(to_rank() == ToEntityRank::value,
                         (debug_message() << "Found connectivity to wrong rank " << to_rank << ", expected " << entity_rank::create(ToEntityRank::value)));
      }
      else {
        if (ToEntityRank::value != entity_rank::node_type::value) {
          BOOST_ASSERT_MSG(*ritr == connectivity_orientation::invalid(), "If key is invalid, then orientation should be invalid");
        }
      }

      if (ToEntityRank::value != entity_rank::node_type::value) {
        ++ritr;
      }
      // TODO - Anything else we can check here?
    }
#endif
  }

  void invariant_check_helper(mesh_impl * mesh = NULL) const
  {
#ifndef NDEBUG
    const bool compressed = mesh != NULL;

    BOOST_ASSERT_MSG(m_ordinal.size() == m_num_connectivity,
                     (debug_message() << "Total size of ordinals " << m_ordinal.size() << " does not match chunk size " << m_num_connectivity));

    if (ToEntityRank::value == entity_rank::node_type::value) {
      BOOST_ASSERT_MSG(m_orientation.empty(), "Orientations should be empty for nodal connectivity");
    }
    else {
      BOOST_ASSERT_MSG(m_orientation.size() == m_target_by_key.size(),
                       (debug_message() << "Total size of orientations " << m_orientation.size() << " does not match size of keys " << m_target_by_key.size()));
    }

    if (compressed) {
      BOOST_ASSERT_MSG(m_target_by_partition_index.size() == m_target_by_key.size(),
                       (debug_message() << "Total size of partition indices " << m_target_by_partition_index.size() << " does not match size of keys " << m_target_by_key.size()));
      BOOST_ASSERT_MSG(m_target_by_rank_index.size() == m_target_by_key.size(),
                       (debug_message() << "Total size of rank indices " << m_target_by_rank_index.size() << " does not match size of keys " << m_target_by_key.size()));

      for (size_t i = 0, e = m_target_by_key.size(); i < e; ++i) {
        const entity_key key = m_target_by_key[i];
        if (key != entity_key::invalid()) {
          const entity_key key_converted_from_partition_index = mesh->convert<entity_key>(m_target_by_partition_index[i]);
          const entity_key key_converted_from_rank_index      = mesh->convert<entity_key>(m_target_by_partition_index[i]);
          BOOST_ASSERT_MSG(key == key_converted_from_partition_index,
                           (debug_message() << "Key converted from partition index " << key_converted_from_partition_index <<
                            " does not match expected key " << key));
          BOOST_ASSERT_MSG(key == key_converted_from_rank_index,
                           (debug_message() << "Key converted from rank index " << key_converted_from_partition_index <<
                            " does not match expected key " << key));
        }
      }
    }
#endif
  }

  // MEMBERS

  spatial_dimension m_spatial_dimension;
  entity_topology m_from_topology;  // Partition topology
  unsigned m_num_connectivity;

  // connectivity data
  std::vector<connectivity_ordinal> m_ordinal;
  std::vector<connectivity_orientation> m_orientation;
  std::vector<entity_key> m_target_by_key;

  // Only updated at end_modification;
  std::vector<partition_index> m_target_by_partition_index;
  std::vector<rank_index> m_target_by_rank_index;
};


//begin connectivity_orientation (entity_rank::node_type)
template <>
template <>
inline
const connectivity_orientation *
partition_connectivity<entity_rank::node_type, connectivity_kind::fixed_type>::begin<connectivity_orientation>(partition_offset /*offset*/) const
{ return NULL; }


//end connectivity_orientation (entity_rank::node_type)
template <>
template <>
inline
const connectivity_orientation *
partition_connectivity<entity_rank::node_type, connectivity_kind::fixed_type>::end<connectivity_orientation>(partition_offset /*offset*/) const
{ return NULL; }

}} //namespace samba::detail

#endif // SAMBA_SAMBA_MESH_FIXED_PARTITION_CONNECTIVITY_TCC
