#ifndef SAMBA_SAMBA_MESH_DYNAMIC_PARTITION_CONNECTIVITY_HPP
#define SAMBA_SAMBA_MESH_DYNAMIC_PARTITION_CONNECTIVITY_HPP

#include <samba/rank_index.hpp>

namespace samba { namespace detail {

template< typename ToEntityRank >
class partition_connectivity<ToEntityRank, connectivity_kind::dynamic_type>
{
  enum connectivity_type { Lower=0,Higher=1,Adjacent=2 };

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
  typedef partition_connectivity<ToEntityRank, connectivity_kind::dynamic_type> self_type;

  static const size_t chunk_size = 1u;

  partition_connectivity() //default constructed partition_connectivity implies connectivity is not used
    : m_spatial_dimension(spatial_dimension::invalid())
    , m_from_topology(entity_topology::invalid())
    , m_from_rank(entity_rank::invalid())
    , m_type(Lower)
    , m_index()
    , m_num_connectivity()
    , m_ordinal()
    , m_orientation()
    , m_target_by_key()
    , m_target_by_partition_index()
    , m_target_by_rank_index()
  {}

  partition_connectivity(entity_topology from_topology, spatial_dimension arg_spatial_dimension)
    : m_spatial_dimension(arg_spatial_dimension)
    , m_from_topology(from_topology)
    , m_from_rank(topology_rank(m_from_topology,m_spatial_dimension))
    , m_type( (m_from_rank > ToEntityRank::value) ? Lower : ((m_from_rank < ToEntityRank::value) ? Higher : Adjacent))
    , m_index()
    , m_num_connectivity()
    , m_ordinal()
    , m_orientation()
    , m_target_by_key()
    , m_target_by_partition_index()
    , m_target_by_rank_index()
  {}

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
  { return begin<Index>(offset) + m_num_connectivity[offset()]; }

  size_t num_connectivity(partition_offset offset) const
  { return m_num_connectivity[offset()]; }

  size_t size() const
  { return m_index.size(); }

  // Modification API
  void add_connectivity(partition_offset from, entity_key to, connectivity_ordinal ordinal, connectivity_orientation orientation)
  {
    BOOST_ASSERT_MSG( topology_rank(to.topology(),m_spatial_dimension) == ToEntityRank::value,
        (debug_message() << "Target rank must be " << entity_rank::create(ToEntityRank::value) << ", " << to.topology()
        << " has " << topology_rank(to.topology(),m_spatial_dimension)));

    switch(m_type)
    {
    case Lower: add_lower_helper(from,to,ordinal,orientation); break;
    case Higher: add_higher_helper(from,to,ordinal,orientation); break;
    case Adjacent: add_adjacent_helper(from,to,ordinal,orientation); break;
    default: BOOST_ASSERT_MSG(false, "What type of connectivity are you trying to add?"); break;
    }

    invariant_check_helper(from);
  }

  void remove_connectivity(partition_offset from, entity_key to, connectivity_ordinal ordinal)
  {
    bool found = false;
    for (uint32_t i = m_index[from()], e = m_index[from()]+m_num_connectivity[from()]; i < e; ++i)
    {
      //remove connectivity
      if ( m_target_by_key[i] == to && m_ordinal[i] == ordinal ) {
        found = true;
        --m_num_connectivity[from()];
      }
      //slide memory down
      if (found && i < e-1u) {
        m_target_by_key[i] = m_target_by_key[i+1u];
        m_ordinal[i] = m_ordinal[i+1u];
        if (has_orientation()) {
          m_orientation[i] = m_orientation[i+1u];
        }
      }
    }

    BOOST_ASSERT_MSG(found,
        (debug_message() << "Unable to remove connectivity to " << to << " at ordinal " << ordinal));

    invariant_check_helper(from);
  }

  void swap(partition_offset first, partition_offset second)
  {
    std::swap( m_index[first()], m_index[second()] );
    std::swap( m_num_connectivity[first()], m_num_connectivity[second()]);
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
    order_by_index();

    {
      std::vector<uint32_t> temp_index(m_index);
      m_index.swap(temp_index);
    }

    {
      std::vector<uint16_t> temp_num_connectivity(m_num_connectivity);
      m_num_connectivity.swap(temp_num_connectivity);
    }

    // Populate partition_index and rank_index structures
    m_target_by_partition_index.resize(m_target_by_key.size(), partition_index::invalid());
    m_target_by_rank_index.resize(m_target_by_key.size(), rank_index::invalid());
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
    const size_t new_size = how_many + m_index.size();

    m_index.resize(new_size, 0); // it's OK for entities with no connectivity to have wrong indices
    m_num_connectivity.resize(new_size, 0);

    invariant_check_helper();
  }

  void remove_entities(size_t how_many)
  {
    const size_t new_size = m_index.size() - how_many;

    m_index.resize(new_size);
    m_num_connectivity.resize(new_size);

    invariant_check_helper();
  }

  void move_entities(self_type & to, size_t how_many)
  {
    order_by_index();
    const size_t num_entities_after_move = size() - how_many;

    //move over index and num_connectivity
    const size_t to_size_before_move = to.size();
    const size_t to_num_connectivity_before_move = to.m_target_by_key.size();
    const size_t num_connectivity_after_move = m_index[num_entities_after_move];

    to.m_index.insert( to.m_index.end()
                      ,m_index.begin() + num_entities_after_move
                      ,m_index.end());

    to.m_num_connectivity.insert( to.m_num_connectivity.end()
                                 ,m_num_connectivity.begin() + num_entities_after_move
                                 ,m_num_connectivity.end());

    //adjust the 'to' index by the difference in the connectivity between target partition before
    //move and source partition after move
    const int64_t index_diff = to_num_connectivity_before_move - num_connectivity_after_move;
    for (size_t i=to_size_before_move, e=to.size(); i<e; ++i)
    { to.m_index[i] += index_diff; }

    //move over target, ordinal, and orientation

    BOOST_ASSERT_MSG( ((m_target_by_key.size() - num_connectivity_after_move) % chunk_size) == 0, "Error, moving partial chunk");
    BOOST_ASSERT_MSG( ((to.m_target_by_key.size()) % chunk_size) == 0, "Error, moving partial chunk");

    to.m_target_by_key.insert( to.m_target_by_key.end()
                              ,m_target_by_key.begin() + num_connectivity_after_move
                              ,m_target_by_key.end()
                             );

    to.m_ordinal.insert( to.m_ordinal.end()
                        ,m_ordinal.begin() + num_connectivity_after_move
                        ,m_ordinal.end()
                       );


    if (has_orientation()) {
      to.m_orientation.insert( to.m_orientation.end()
                              ,m_orientation.begin() + num_connectivity_after_move
                              ,m_orientation.end()
                             );
    }

    remove_entities(how_many);

    invariant_check_helper();
    to.invariant_check_helper();
  }

private:

  static size_t num_chunks(size_t num)
  { return (num + chunk_size -1u)/chunk_size; }

  bool has_orientation() const
  {
    const static bool rv = ToEntityRank::value != entity_rank::node_type::value && m_from_rank != entity_rank::node();
    return rv;
  }

  //begin entity_key
  const entity_key * begin(entity_key_helper offset) const
  { return &*m_target_by_key.begin() + m_index[offset()]; }

  //begin partition_index
  const partition_index * begin(partition_index_helper offset) const
  { return &*m_target_by_partition_index.begin() + m_index[offset()]; }

  //begin rank_index
  const rank_index * begin(rank_index_helper offset) const
  { return &*m_target_by_rank_index.begin() + m_index[offset()]; }

  //begin connectivity_ordinal
  const connectivity_ordinal * begin(connectivity_ordinal_helper offset) const
  { return &*m_ordinal.begin() + m_index[offset()]; }

  //begin connectivity_orientation
  const connectivity_orientation * begin(connectivity_orientation_helper offset) const
  {
    if (has_orientation()) {
      return &*m_orientation.begin() + m_index[offset()];
    }
    else {
      return NULL;
    }
  }

  void order_by_index(size_t capacity = 0u)
  {
    //compute needed capacity
    if (capacity == 0u) {
      for( size_t i=0, e=m_index.size(); i<e; ++i) {
        capacity += num_chunks(m_num_connectivity[i]);
      }
    }

    //move orientation
    if (has_orientation()) {
      std::vector<connectivity_orientation> temp_orientation;
      temp_orientation.reserve(capacity);

      uint32_t current_index=0;
      for( size_t i=0, e=m_index.size(); i<e; ++i)
      {
        const size_t begin_offset = m_index[i];
        const size_t end_offset = begin_offset + num_chunks(m_num_connectivity[i])*chunk_size;

        if (begin_offset != end_offset) {
          temp_orientation.insert(temp_orientation.end(), m_orientation.begin()+begin_offset, m_orientation.begin() + end_offset);
        }
        current_index += static_cast<uint32_t>(chunk_size*num_chunks(m_num_connectivity[i]));
      }

      temp_orientation.swap(m_orientation);
    }

    //move ordinal
    {
      std::vector<connectivity_ordinal> temp_ordinal;
      temp_ordinal.reserve(capacity);

      uint32_t current_index=0;
      for( size_t i=0, e=m_index.size(); i<e; ++i)
      {
        const size_t begin_offset = m_index[i];
        const size_t end_offset = begin_offset + num_chunks(m_num_connectivity[i])*chunk_size;

        if (begin_offset != end_offset) {
          temp_ordinal.insert(temp_ordinal.end(), m_ordinal.begin()+begin_offset, m_ordinal.begin() + end_offset);
        }
        current_index += static_cast<uint32_t>(chunk_size*num_chunks(m_num_connectivity[i]));
      }
      temp_ordinal.swap(m_ordinal);
    }

    //move target_by_key
    {
      std::vector<entity_key> temp_target_by_key;
      temp_target_by_key.reserve(capacity);

      uint32_t current_index=0;
      for( size_t i=0, e=m_index.size(); i<e; ++i)
      {
        const size_t begin_offset = m_index[i];
        const size_t end_offset = begin_offset + num_chunks(m_num_connectivity[i])*chunk_size;

        if (begin_offset != end_offset) {
          temp_target_by_key.insert(temp_target_by_key.end(), m_target_by_key.begin()+begin_offset, m_target_by_key.begin() + end_offset);
        }
        m_index[i] = current_index; //update index
        current_index += static_cast<uint32_t>(chunk_size*num_chunks(m_num_connectivity[i]));
      }
      temp_target_by_key.swap(m_target_by_key);
    }
  }

  void add_connectivity_helper(partition_offset from)
  {
    const size_t chunks_needed = num_chunks(m_num_connectivity[from()]+1);
    const size_t chunks_used = num_chunks(m_num_connectivity[from()]);

    if (chunks_needed == chunks_used)
    {
      ++m_num_connectivity[from()];
      return;
    }

    const size_t chunks_available = num_chunks(m_target_by_key.capacity() - m_target_by_key.size());

    if (chunks_available < chunks_needed)
    {
      const size_t new_capacity = m_target_by_key.capacity() > 0u ? 2u*m_target_by_key.capacity() : 8u*chunk_size;
      order_by_index(new_capacity);
    }

    const bool last_entity_by_index = (chunks_used > 0u) && (m_index[from()] + chunks_used*chunk_size == m_target_by_key.size());

    //copy to end
    if (!last_entity_by_index)
    {
      uint32_t new_index = static_cast<uint32_t>(m_target_by_key.size());

      m_target_by_key.insert(m_target_by_key.end(),chunks_needed*chunk_size,entity_key::invalid());
      std::copy(begin<entity_key>(from), end<entity_key>(from), m_target_by_key.begin() + new_index);

      m_ordinal.insert(m_ordinal.end(),chunks_needed*chunk_size,connectivity_ordinal::invalid());
      std::copy(begin<connectivity_ordinal>(from), end<connectivity_ordinal>(from), m_ordinal.begin() + new_index);

      if (has_orientation()) {
        m_orientation.insert(m_orientation.end(),chunks_needed*chunk_size,connectivity_orientation::invalid());
        std::copy(begin<connectivity_orientation>(from), end<connectivity_orientation>(from), m_orientation.begin() + new_index);
      }

      m_index[from()] = new_index;
    }
    //add new chunk to end
    else {
      m_target_by_key.insert(m_target_by_key.end(),chunk_size,entity_key::invalid());
      m_ordinal.insert(m_ordinal.end(),chunk_size,connectivity_ordinal::invalid());
      if (has_orientation()) {
        m_orientation.insert(m_orientation.end(),chunk_size,connectivity_orientation::invalid());
      }
    }

    ++m_num_connectivity[from()];
  }

  void add_lower_helper(partition_offset from, entity_key to, connectivity_ordinal ordinal, connectivity_orientation orientation)
  {

#ifndef NDEBUG
    unsigned num_ordinals = 0;
    switch(ToEntityRank::value)
    {
    case entity_rank::node_type::value: num_ordinals = num_nodes(m_from_topology); break;
    case entity_rank::edge_type::value: num_ordinals = num_edges(m_from_topology); break;
    case entity_rank::face_type::value: num_ordinals = num_faces(m_from_topology); break;
    case entity_rank::element_type::value: num_ordinals = num_sides(m_from_topology); break;
    }

    BOOST_ASSERT_MSG(ordinal < num_ordinals,
        (debug_message() << ordinal << " exceeds topological limit: " << num_ordinals));

    entity_topology downward_topology = connectivity_topology(m_from_topology, entity_rank::create(ToEntityRank::value), ordinal);
    BOOST_ASSERT_MSG( downward_topology == entity_topology::invalid() || downward_topology == to.topology(),
        (debug_message() << "Expected connectivity to entity with " << downward_topology << ", found " << to.topology()));

#endif

    add_connectivity_helper(from);

    const uint32_t begin_index = m_index[from()]+m_num_connectivity[from()] -1u;

    if (m_num_connectivity[from()] == 1u) {
      m_target_by_key[begin_index] = to;
      m_ordinal[begin_index] = ordinal;
      if (has_orientation()) {
        m_orientation[begin_index] = orientation;
      }
      return;
    }

    for (uint32_t i = begin_index, e = m_index[from()]; i > e; --i)
    {
      //slide up
      if (ordinal < m_ordinal[i-1u] ) {
        m_target_by_key[i] = m_target_by_key[i-1u];
        m_ordinal[i] = m_ordinal[i-1u];
        if (has_orientation()) {
          m_orientation[i] = m_orientation[i];
        }
        //insert if on last iteration
        if ((i-1)==e) {
          m_target_by_key[i-1u] = to;
          m_ordinal[i-1u] = ordinal;
          if (has_orientation()) {
            m_orientation[i-1u] = orientation;
          }
        }
      }
      //insert
      else if (ordinal > m_ordinal[i-1u]) {
        m_target_by_key[i] = to;
        m_ordinal[i] = ordinal;
        if (has_orientation()) {
          m_orientation[i] = orientation;
        }
        break;
      }
      //duplicate -- insert new and remove the original
      else
      {
        m_target_by_key[i] = to;
        m_ordinal[i] = ordinal;
        if (has_orientation()) {
          m_orientation[i] = orientation;
        }
        remove_connectivity(from,to,ordinal);
        break;
      }
    }
  }

  void add_higher_helper(partition_offset from, entity_key to, connectivity_ordinal ordinal, connectivity_orientation orientation)
  {

#ifndef NDEBUG
    unsigned num_ordinals = 0;
    switch(m_from_rank())
    {
    case entity_rank::node_type::value: num_ordinals = num_nodes(to.topology()); break;
    case entity_rank::edge_type::value: num_ordinals = num_edges(to.topology()); break;
    case entity_rank::face_type::value: num_ordinals = num_faces(to.topology()); break;
    case entity_rank::element_type::value: num_ordinals = num_sides(to.topology()); break;
    }

    BOOST_ASSERT_MSG(ordinal < num_ordinals,
        (debug_message() << ordinal << " exceeds topological limit: " << num_ordinals));

    entity_topology downward_topology = connectivity_topology(to.topology(), m_from_rank, ordinal);
    BOOST_ASSERT_MSG( downward_topology == entity_topology::invalid() || downward_topology== m_from_topology,
        (debug_message() << "Expected connectivity to entity with " << downward_topology << ", found " << m_from_topology));

#endif

    add_connectivity_helper(from);

    const uint32_t begin_index = m_index[from()]+m_num_connectivity[from()] -1u;

    if (m_num_connectivity[from()] == 1u) {
      m_target_by_key[begin_index] = to;
      m_ordinal[begin_index] = ordinal;
      if (has_orientation()) {
        m_orientation[begin_index] = orientation;
      }
      return;
    }

    for (uint32_t i = begin_index, e = m_index[from()]; i > e; --i)
    {
      //slide up
      if (std::make_pair(to,ordinal) < std::make_pair(m_target_by_key[i-1u],m_ordinal[i-1u]) ) {
        m_target_by_key[i] = m_target_by_key[i-1u];
        m_ordinal[i] = m_ordinal[i-1u];
        if (has_orientation()) {
          m_orientation[i] = m_orientation[i];
        }
        //insert if on last iteration
        if ((i-1u)==e) {
          m_target_by_key[i-1u] = to;
          m_ordinal[i-1u] = ordinal;
          if (has_orientation()) {
            m_orientation[i-1u] = orientation;
          }
        }
      }
      //insert
      else if (std::make_pair(to,ordinal) > std::make_pair(m_target_by_key[i-1u],m_ordinal[i-1u]) ) {
        m_target_by_key[i] = to;
        m_ordinal[i] = ordinal;
        if (has_orientation()) {
          m_orientation[i] = orientation;
        }
        break;
      }
      //duplicate -- insert new and remove the original
      else
      {
        m_target_by_key[i] = to;
        m_ordinal[i] = ordinal;
        if (has_orientation()) {
          m_orientation[i] = orientation;
        }
        remove_connectivity(from,to,ordinal);
        break;
      }
    }
  }

  void add_adjacent_helper(partition_offset from, entity_key to, connectivity_ordinal ordinal, connectivity_orientation orientation)
  {

#ifndef NDEBUG
    entity_topology from_side = side_topology(m_from_topology, ordinal());

    BOOST_ASSERT_MSG(orientation < num_sides(m_from_topology) || num_sides(m_from_topology)==0,
        (debug_message() << "Ordinal " << ordinal << " exceeds topological limit: " << num_sides(m_from_topology) ));

    BOOST_ASSERT_MSG(orientation < num_sides(to.topology()) || num_sides(to.topology())==0,
        (debug_message() << "Ordinal " << orientation.ordinal() << " exceeds topological limit: " << num_sides(to.topology())));

    entity_topology to_side = side_topology(to.topology(), orientation());

    BOOST_ASSERT_MSG( from_side == to_side,
        (debug_message() << "Adjacent entities must have matching side topologies" << from_side << ", " << to_side));
#endif

    add_connectivity_helper(from);

    const uint32_t begin_index = m_index[from()]+m_num_connectivity[from()] -1u;

    if (m_num_connectivity[from()] == 1u) {
      m_target_by_key[begin_index] = to;
      m_ordinal[begin_index] = ordinal;
      if (has_orientation()) {
        m_orientation[begin_index] = orientation;
      }
      return;
    }

    for (uint32_t i = begin_index, e = m_index[from()]; i > e; --i)
    {
      //slide up
      if (ordinal < m_ordinal[i-1u] ) {
        m_target_by_key[i] = m_target_by_key[i-1u];
        m_ordinal[i] = m_ordinal[i-1u];
        if (has_orientation()) {
          m_orientation[i] = m_orientation[i];
        }
        //insert if on last iteration
        if ((i-1u)==e) {
          m_target_by_key[i-1u] = to;
          m_ordinal[i-1u] = ordinal;
          if (has_orientation()) {
            m_orientation[i-1u] = orientation;
          }
        }
      }
      //insert
      else if (ordinal > m_ordinal[i-1u]) {
        m_target_by_key[i] = to;
        m_ordinal[i] = ordinal;
        if (has_orientation()) {
          m_orientation[i] = orientation;
        }
        break;
      }
      //duplicate -- insert new and remove the original
      else
      {
        m_target_by_key[i] = to;
        m_ordinal[i] = ordinal;
        if (has_orientation()) {
          m_orientation[i] = orientation;
        }
        remove_connectivity(from,to,ordinal);
        break;
      }
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

    BOOST_ASSERT_MSG(static_cast<size_t>(keys_end - keys_begin) == num_connectivity(offset),
                     (debug_message() << "Expected data to be of size " << num_connectivity(offset) << ", " << offset << " has keys " << keys_end - keys_begin));

    BOOST_ASSERT_MSG(keys_end - keys_begin == ordinals_end - ordinals_begin,
                     (debug_message() << "Num keys, " << keys_end - keys_begin << ", does not match num ordinals, " << ordinals_end - ordinals_begin));
    if (has_orientation()) {
      BOOST_ASSERT_MSG(keys_end - keys_begin == orientations_end - orientations_begin,
                       (debug_message() << "Num keys, " << keys_end - keys_begin << ", does not match num orientations, " << orientations_end - orientations_begin));
    }
    else {
      BOOST_ASSERT_MSG(orientations_end - orientations_begin == 0,
                       (debug_message() << "Expected 0 orientations for node connectivity, found: " << orientations_end - orientations_begin));
    }

    const entity_key*               kitr = keys_begin;
    const connectivity_ordinal*     oitr = ordinals_begin;
    for ( ; kitr != keys_end; ++kitr, ++ oitr) {
      BOOST_ASSERT_MSG(*kitr != entity_key::invalid(),
                       "Should not have invalid connectivity for dynamic connectivity");
      entity_rank to_rank  = topology_rank(kitr->topology(), m_spatial_dimension);
      BOOST_ASSERT_MSG(to_rank() == ToEntityRank::value,
                       (debug_message() << "Found connectivity to wrong rank " << to_rank << ", expected " << entity_rank::create(ToEntityRank::value)));
      if (kitr + 1 != keys_end) {
        if (m_type == Higher) { // back rel
          BOOST_ASSERT_MSG(std::make_pair(*kitr, *oitr) < std::make_pair(*(kitr + 1), *(oitr + 1)),
                           (debug_message() << "Connectivity out of order; data at " << kitr - keys_begin <<
                            "\nis (" << *oitr << ", " << *kitr << "),\ndata at next slot is (" << *(oitr + 1) << ", " << *(kitr + 1) << ")"));
        }
        else {
          BOOST_ASSERT_MSG(*oitr < *(oitr + 1),
                           (debug_message() << "Connectivity out of order; data at " << kitr - keys_begin <<
                            "\nis (" << *oitr << ", " << *kitr << "),\ndata at next slot is (" << *(oitr + 1) << ", " << *(kitr + 1) << ")"));
        }
      }
      // TODO - Anything else we can check here?
    }

    invariant_check_helper();
#endif
  }

  void invariant_check_helper(mesh_impl * mesh = NULL) const
  {
#ifndef NDEBUG
    const bool compressed = mesh != NULL;

      // Check sizes
    BOOST_ASSERT_MSG(m_num_connectivity.size() == m_index.size(),
                     (debug_message() << "Expected m_num_connectivity to be of size " << m_index.size() << ", found " << m_num_connectivity.size()));

    BOOST_ASSERT_MSG(m_target_by_key.size() == m_ordinal.size(),
                     (debug_message() << "Total size of keys " << m_target_by_key.size() << " does not match size of ordinals " << m_ordinal.size()));

    if (has_orientation()) {
      BOOST_ASSERT_MSG(m_orientation.size() == m_target_by_key.size(),
                       (debug_message() << "Total size of orientations " << m_orientation.size() << " does not match size of keys " << m_target_by_key.size()));
    }
    else {
      BOOST_ASSERT_MSG(m_orientation.empty(), "Orientations should be empty for nodal connectivity");
    }

    for (size_t o = 1, e = m_index.size(); o < e; ++o) {
      partition_offset current_offset = partition_offset::create(o);
      partition_offset prior_offset   = partition_offset::create(o - 1);
      const size_t curr_index         = m_index[o];
      BOOST_ASSERT_MSG(curr_index <= m_target_by_key.size(),
                       (debug_message() << "Index is wrong, " << curr_index << " is beyond max " << m_target_by_key.size()));
      if (compressed) {
        const size_t index_diff     = curr_index - m_index[o-1];
        const size_t prior_num_conn = num_connectivity(prior_offset);
        BOOST_ASSERT_MSG(prior_num_conn == index_diff,
                         (debug_message() << "For offset " << current_offset << ", num_connectivity/index mismatch, index_diff is " << index_diff << ", num conn is " << prior_num_conn));
      }
    }

    // Check that connectivity is in-sync
    if (compressed) {
      BOOST_ASSERT_MSG(m_target_by_partition_index.size() == m_ordinal.size(),
                       (debug_message() << "Total size of partition indices " << m_target_by_partition_index.size() << " does not match size of ordinals " << m_ordinal.size()));
      BOOST_ASSERT_MSG(m_target_by_rank_index.size() == m_ordinal.size(),
                       (debug_message() << "Total size of rank indices " << m_target_by_rank_index.size() << " does not match size of ordinals " << m_ordinal.size()));

      for (size_t i = 0, e = m_target_by_key.size(); i < e; ++i) {
        const entity_key key = m_target_by_key[i];
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
#endif
  }

  // MEMBERS

  spatial_dimension m_spatial_dimension;
  entity_topology m_from_topology;
  entity_rank m_from_rank;
  connectivity_type m_type;

  // meta data
  std::vector<uint32_t> m_index;  // Common index into vectors below that stores where connectivity starts for a partition_offset (entity).
  std::vector<uint16_t> m_num_connectivity;

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
partition_connectivity<entity_rank::node_type, connectivity_kind::dynamic_type>::begin<connectivity_orientation>(partition_offset /*offset*/) const
{
  return NULL;
}

//end connectivity_orientation (entity_rank::node_type)
template <>
template <>
inline
const connectivity_orientation *
partition_connectivity<entity_rank::node_type, connectivity_kind::dynamic_type>::end<connectivity_orientation>(partition_offset /*offset*/) const
{
  return NULL;
}

template <>
template <>
inline
const connectivity_orientation *
partition_connectivity<entity_rank::edge_type, connectivity_kind::dynamic_type>::end<connectivity_orientation>(partition_offset offset) const
{
  if (has_orientation()) {
    return begin<connectivity_orientation>(offset) + m_num_connectivity[offset()];
  }
  else {
    return NULL;
  }
}

template <>
template <>
inline
const connectivity_orientation *
partition_connectivity<entity_rank::face_type, connectivity_kind::dynamic_type>::end<connectivity_orientation>(partition_offset offset) const
{
  if (has_orientation()) {
    return begin<connectivity_orientation>(offset) + m_num_connectivity[offset()];
  }
  else {
    return NULL;
  }
}

template <>
template <>
inline
const connectivity_orientation *
partition_connectivity<entity_rank::element_type, connectivity_kind::dynamic_type>::end<connectivity_orientation>(partition_offset offset) const
{
  if (has_orientation()) {
    return begin<connectivity_orientation>(offset) + m_num_connectivity[offset()];
  }
  else {
    return NULL;
  }
}


}} //namespace samba::detail

#endif // SAMBA_SAMBA_MESH_DYNAMIC_PARTITION_CONNECTIVITY_HPP
