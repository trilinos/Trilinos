#ifndef SAMBA_SAMBA_MESH_MESH_IMPL_TCC
#define SAMBA_SAMBA_MESH_MESH_IMPL_TCC

#include <samba/utility/stream_vector.hpp>

#include <samba/rank_index.hpp>

#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/copy.hpp>

#include <iostream>

namespace samba { namespace detail {

struct order_partition_ids_by_parts
{
  order_partition_ids_by_parts(const std::vector<mesh_impl::partition_handle>& partition_impls)
    : m_partition_impls(partition_impls)
  {}

  bool operator()(partition_id lhs, partition_id rhs)
  {
    if (lhs == partition_id::invalid()) {
      return false;
    }
    else if (rhs == partition_id::invalid()) {
      return true;
    }
    return m_partition_impls[lhs()]->parts() < m_partition_impls[rhs()]->parts();
  }

  std::vector<mesh_impl::partition_handle> const & m_partition_impls;
};


//*************************************************************************
//constructor
//*************************************************************************
inline
mesh_impl::mesh_impl( samba::connectivity_map const& arg_connectivity_map
                     ,process_id arg_process
                    )
  : m_process(arg_process)
  , m_is_in_modification_cycle(true)
  , m_connectivity_map(arg_connectivity_map)
  , m_set_specs()
  , m_inducable_parts()
  , m_parts_to_partition()
  , m_partitions()
  , m_local_ids()
  , m_deleted_keys()
  , m_non_local_key_to_descriptor_map()
  , m_key_to_descriptor_map()
  , m_descriptor_to_key_map()
  , m_signals()
  , m_node_index_to_partition_index_map()
  , m_edge_index_to_partition_index_map()
  , m_face_index_to_partition_index_map()
  , m_element_index_to_partition_index_map()
{}

//*************************************************************************
//signals
//*************************************************************************

inline mesh_signals & mesh_impl::signals()
{ return m_signals; }

inline mesh_signals const& mesh_impl::signals() const
{ return m_signals; }


//*************************************************************************
//modification
//*************************************************************************

inline void mesh_impl::add_set_expression_by_name(const std::string& name, set_expression expr)
{ m_expr_db[name] = expr; }

inline bool mesh_impl::begin_modification()
{
  if (!m_is_in_modification_cycle) {
    m_is_in_modification_cycle = true;

    // notify partitions that we have entered a modification cycle
    for (size_t p = 0, end_partition = m_partitions.size(); p < end_partition; ++p) {
      m_partitions[p]->begin_modification_impl();
    }

    m_signals.begin_modification_signal(); // signal

    // Clear rank -> partition index maps. Rank indices are not reliable during modification.
    {
      std::vector<partition_index> empty_vec;
      m_node_index_to_partition_index_map.swap(empty_vec);
    }

    {
      std::vector<partition_index> empty_vec;
      m_edge_index_to_partition_index_map.swap(empty_vec);
    }

    {
      std::vector<partition_index> empty_vec;
      m_face_index_to_partition_index_map.swap(empty_vec);
    }

    {
      std::vector<partition_index> empty_vec;
      m_element_index_to_partition_index_map.swap(empty_vec);
    }

    return true;
  }
  else {
    return false;
  }
}

inline void mesh_impl::remove_relations_to_deleted_entities_helper()
{
  for (std::vector<partition_handle>::const_iterator itr = m_partitions.begin(), end = m_partitions.end();
       itr != end;
       ++itr) {
    partition_handle curr_partition = *itr;

    // Iterate backwards wherever possible to increase liklihood that this algorithm
    // will work regardless of partition-connectivity implementation
    for (size_t o = curr_partition->size(), end = 0; o > end; --o) {
      partition_offset offset = partition_offset::create(o-1);
      for (entity_rank rank = entity_rank::node(); rank <= entity_rank::element(); ++rank) {
        if (curr_partition->connectivity_kind(rank) != connectivity_kind::invalid()) {
          entity_key_iterator begin = curr_partition->begin_connectivity<entity_key>(rank, offset);
          entity_key_iterator end   = curr_partition->end_connectivity<entity_key>(rank, offset);
          ordinal_iterator oitr     = curr_partition->end_connectivity<connectivity_ordinal>(rank, offset);

          for (entity_key_iterator itr = end; itr != begin; --itr, --oitr) {
            entity_key target            = *(itr - 1);
            connectivity_ordinal ordinal = *(oitr - 1);

            // We know an entity has been deleted if the key maps to an invalid descriptor
            if (target != entity_key::invalid()) {
              if (has_been_deleted_helper(target)) {
                curr_partition->remove_connectivity(offset, target, ordinal);
              }
            }
            else {
              BOOST_ASSERT_MSG(curr_partition->connectivity_kind(rank) == connectivity_kind::fixed(),
                               "Should only have invalid connectivity for fixed kind");
            }
          }
        }
      }
    }
  }
}

inline bool mesh_impl::has_been_deleted_helper(entity_key key) const
{
  if (key.process() == m_process) {
    return m_key_to_descriptor_map[key.topology()()][key.local_id()()] == partition_index::invalid();
  }
  else {
    return m_non_local_key_to_descriptor_map.find(key) == m_non_local_key_to_descriptor_map.end();
  }
}

struct bitset_compare
{
  template <size_t N>
  bool operator()(std::bitset<N> const& lhs, std::bitset<N> const& rhs) const
  {
    return lhs.to_ulong() < rhs.to_ulong();
  }
};

struct order_keys_by_partition
{
  order_keys_by_partition(mesh_impl const* mesh_arg) : m_mesh(mesh_arg) {}

  bool operator()(entity_key lhs, entity_key rhs) const
  {
    partition_index lhs_pindex = m_mesh->convert(lhs);
    partition_index rhs_pindex = m_mesh->convert(rhs);
    return lhs_pindex.partition() < rhs_pindex.partition();
  }

  mesh_impl const* m_mesh;
};

// alternate induced part idea: Have a separate contain representing unique collections of inducable parts, every entity
// key will map to an index in this vector, representing what parts it needs to be added to at the end.
// Not having an entry in this map means you should be removed from all induced parts and added to none.

inline void mesh_impl::apply_induced_part_membership_helper()
{
  // Use bitset to represent part collections to improve performance when there
  // are a small number of inducable parts
  BOOST_ASSERT_MSG(entity_topology::invalid()() + m_inducable_parts.size() < MAX_INDUCED_PARTS,
                   "There are too many induced parts to use the fast algorithm");

  // Set up bitset maps
  std::vector< std::vector<induced_part_set> > key_to_induced_part_map(m_key_to_descriptor_map.size());
  for (size_t i = 0, e = key_to_induced_part_map.size(); i < e; ++i) {
    entity_topology topo = entity_topology::create(i);
    if (topology_rank(topo, spatial_dimension()) < entity_rank::element()) {
      key_to_induced_part_map[i].resize(m_key_to_descriptor_map[i].size());
    }
  }

  // Since topologies are inducable, all partitions have entities in inducable parts
  //   for each partition, iterate over all downward connectivity
  //     for each related entity, add/modify the bitset map above, adding the new induced parts

  for (std::vector<partition_handle>::const_iterator itr = m_partitions.begin(), end = m_partitions.end();
       itr != end;
       ++itr) {
    partition_handle p_handle = *itr;
    if (p_handle->rank() != entity_rank::node() && !p_handle->empty()) {
      induce_parts_through_downward_connectivity_helper(p_handle, key_to_induced_part_map);
    }
  }

  // Invert the new_parts_map, giving us collections of entities being added to the same induced parts.
  // We must store keys because the handles must stay valid while things are being moved around
  typedef std::vector<entity_key> entity_key_vector;
  typedef std::vector<partition_index> partition_index_vector;
  typedef std::map<induced_part_set, entity_key_vector, bitset_compare> parts_to_entities_map;
  parts_to_entities_map new_parts_to_entities_map;
  for (size_t topology_id = 0; topology_id < key_to_induced_part_map.size(); ++topology_id) {
    entity_topology topology = entity_topology::create(topology_id);
    for (size_t local_id_int = 0; local_id_int < key_to_induced_part_map[topology_id].size(); ++local_id_int) {
      induced_part_set& induced_parts = key_to_induced_part_map[topology_id][local_id_int];
      entity_local_id local_id = entity_local_id::create(local_id_int);
      entity_key entity = entity_key::create(topology, m_process, local_id);
      if (!has_been_deleted_helper(entity)) {
        parts_to_entities_map::iterator fitr = new_parts_to_entities_map.find(induced_parts);
        if (fitr != new_parts_to_entities_map.end()) {
          fitr->second.push_back(entity);
        }
        else {
          entity_key_vector temp(1, entity);
          new_parts_to_entities_map[induced_parts] = temp;
        }
      }
    }
  }

  // For each collection of partition-indices being moved into a common destination partition, group into
  // common source partitions and do the move in bulk
  entity_part_vector new_parts; // declared up here to avoid lots of little allocations
  entity_part_vector parts_to_add; // declared up here to avoid lots of little allocations
  partition_index_vector indices_to_move; // declared up here to avoid lots of little allocations
  for (parts_to_entities_map::iterator mitr = new_parts_to_entities_map.begin(), end = new_parts_to_entities_map.end();
       mitr != end;
       ++mitr) {
    induced_part_set const& set_of_parts_to_add = mitr->first;
    entity_key_vector& entities = mitr->second;
    std::sort(entities.begin(), entities.end(), order_keys_by_partition(this));

    parts_to_add.clear();
    parts_to_add.reserve(set_of_parts_to_add.count());
    for (size_t b = 0, be = MAX_INDUCED_PARTS; b < be; ++b) {
      if (set_of_parts_to_add.test(b)) {
        if (b < entity_topology::invalid()()) {
          parts_to_add.push_back(entity_part(entity_topology::create(b), true /*induced*/));
        }
        else {
          parts_to_add.push_back(entity_part(m_inducable_parts[b - entity_topology::invalid()()], true /*induced*/));
        }
      }
    }

    partition_id common_partition = convert_helper<partition_index>(entities.front()).partition();
    size_t common_partition_span = 0;
    for (size_t i = 0, e = entities.size(); i <= e; ++i) {
      partition_id curr_partition = i == e ? partition_id::invalid() : convert_helper<partition_index>(entities[i]).partition();
      if (curr_partition != common_partition) {
        entity_part_vector const& old_parts = m_partitions[common_partition()]->parts();
        size_t induced_idx = find_induced_parts_helper(old_parts);
        new_parts.clear();
        std::copy(old_parts.begin(), old_parts.begin() + induced_idx, std::back_inserter(new_parts));
        std::copy(parts_to_add.begin(), parts_to_add.end(), std::back_inserter(new_parts));
        partition_id new_partition = get_partition_non_const(new_parts);

        if (new_partition != common_partition) {
          indices_to_move.clear();
          indices_to_move.reserve(i - common_partition_span);
          for (size_t k = common_partition_span; k < i; ++k) {
            indices_to_move.push_back(convert_helper<partition_index>(entities[k]));
            BOOST_ASSERT(indices_to_move.back().partition() == common_partition);
          }

          move_entities_helper(common_partition, new_partition, indices_to_move.begin(), indices_to_move.end(), true /*during end mod*/);
        }

        common_partition = curr_partition;
        common_partition_span = i;
      }
    }
  }
}

inline void mesh_impl::compress_data_structure_helper(std::vector<partition_id>& new_partition_order)
{
  //
  // partitions - remove empty paritiions, sort based on parts
  //

  // create new descriptor vector, setting empty partitions to invalid
  new_partition_order.reserve(m_partitions.size());
  for (size_t p = 0, partitions_end = m_partitions.size(); p < partitions_end; ++p) {
    if (m_partitions[p]->empty()) {
      new_partition_order.push_back(partition_id::invalid());
    }
    else {
      partition_id descriptor = m_partitions[p]->partition();
      BOOST_ASSERT_MSG(descriptor() == p, "Invalid partition, out of order.");
      new_partition_order.push_back(descriptor);
    }
  }

  // remove invalid (empty) partitions, sort remaining partitions by their parts;
  // the result is that new_partition_order holds the old descriptors in the updated ordering.
  std::sort(new_partition_order.begin(), new_partition_order.end(),
            order_partition_ids_by_parts(m_partitions));
  new_partition_order.erase(
    std::find(new_partition_order.begin(), new_partition_order.end(), partition_id::invalid()),
    new_partition_order.end());

  // update m_partitions based on new-order. old descriptors become invalid here
  std::vector<partition_handle> new_partitions;
  new_partitions.reserve(new_partition_order.size());
  for (size_t p = 0, partitions_end = new_partition_order.size(); p < partitions_end; ++p) {
    partition_id old_descriptor = new_partition_order[p];
    new_partitions.push_back(m_partitions[old_descriptor()]);
    new_partitions.back()->update_partition_id(partition_id::create(p));
  }
  m_partitions.swap(new_partitions);

  // update m_parts_to_partition based on new partition_ids
  part_vector_to_partition_map new_parts_to_partition;
  for (std::vector<partition_handle>::const_iterator itr = m_partitions.begin(), end = m_partitions.end();
       itr != end;
       ++itr) {
    new_parts_to_partition[(*itr)->parts()] = (*itr)->partition();
  }
  m_parts_to_partition.swap(new_parts_to_partition);
  BOOST_ASSERT_MSG(m_parts_to_partition.size() == m_partitions.size(), "Size mismatch");

  //
  // key/descriptor maps need to be updated
  //

  // update key-map. We can do this in-bulk since the only thing that has changed is that
  // partition descriptors, not the offsets of the entities within the paritions.
  std::vector<std::vector<entity_key> > new_descriptor_to_key_map(new_partition_order.size());
  for (size_t p = 0, partitions_end = new_partition_order.size(); p < partitions_end; ++p) {
    partition_id old_descriptor = new_partition_order[p];
    new_descriptor_to_key_map[p].swap(m_descriptor_to_key_map[old_descriptor()]);
  }
  m_descriptor_to_key_map.swap(new_descriptor_to_key_map);
  BOOST_ASSERT_MSG(m_descriptor_to_key_map.size() == m_partitions.size(), "Size mismatch");

  // update descriptors. This is potentially expensive since a naive implementation would
  // iterate over every entity in the system. To reduce the expense, we take advantage of the
  // fact that, if a partition remains at the same index in the m_partitions vector, then the
  // partition_indexs for entities in that partition will remain the same.
  // This must be done after updating the key-map because we use the key-map in the inner loop below.
  for (size_t p = 0, end_partition = m_partitions.size(); p < end_partition; ++p) {
    if (new_partition_order[p]() != p) { // partition moved
      for(size_t offset = 0, offset_end = m_partitions[p]->size(); offset < offset_end; ++offset) {
        entity_key key = m_descriptor_to_key_map[p][offset];
        m_key_to_descriptor_map[key.topology()()][key.local_id()()] =
          partition_index::create(topology_rank(key.topology(), spatial_dimension()),
                                  partition_id::create(p),
                                  partition_offset::create(offset));
      }
    }
  }
}

inline void mesh_impl::populate_rank_index_helper()
{
  //
  // Populate rank index maps. Needs to happen before end-modification on individual partitions
  // so they can use these maps for conversions to set up their rank-index connectivity.
  //

  // Sizing run
  uint32_t num_nodes = 0, num_edges = 0, num_faces = 0, num_elements = 0;
  for (size_t p = 0, end_partition = m_partitions.size(); p < end_partition; ++p) {
    partition_handle partition = m_partitions[p];
    switch (partition->rank()()) {
    case entity_rank::node_type::value:
      num_nodes += partition->size();
      break;
    case entity_rank::edge_type::value:
      num_edges += partition->size();
      break;
    case entity_rank::face_type::value:
      num_faces += partition->size();
      break;
    case entity_rank::element_type::value:
      num_elements += partition->size();
      break;
    default:
      BOOST_ASSERT_MSG(false, "Should never make it here");
    }
  }

  m_node_index_to_partition_index_map.reserve(num_nodes);
  m_edge_index_to_partition_index_map.reserve(num_edges);
  m_face_index_to_partition_index_map.reserve(num_faces);
  m_element_index_to_partition_index_map.reserve(num_elements);

  uint32_t node_count = 0, edge_count = 0, face_count = 0, element_count = 0;
  for (size_t p = 0, end_partition = m_partitions.size(); p < end_partition; ++p) {
    partition_handle partition = m_partitions[p];

    std::vector<partition_index>* rank_index_to_partition_index_map = NULL;
    uint32_t* count = NULL;
    switch (partition->rank()()) {
    case entity_rank::node_type::value:
      rank_index_to_partition_index_map = &m_node_index_to_partition_index_map;
      count = &node_count;
      break;
    case entity_rank::edge_type::value:
      rank_index_to_partition_index_map = &m_edge_index_to_partition_index_map;
      count = &edge_count;
      break;
    case entity_rank::face_type::value:
      rank_index_to_partition_index_map = &m_face_index_to_partition_index_map;
      count = &face_count;
      break;
    case entity_rank::element_type::value:
      rank_index_to_partition_index_map = &m_element_index_to_partition_index_map;
      count = &element_count;
      break;
    default:
      BOOST_ASSERT_MSG(false, "Should never make it here");
    }

    for (partition_offset offset = {0}; offset < partition->size(); ++offset) {
      rank_index_to_partition_index_map->push_back(partition->descriptor(offset));
    }

    partition->set_starting_rank_count(*count);
    *count += partition->size();
  }
}

inline bool mesh_impl::end_modification()
{
  if (m_is_in_modification_cycle) {

    remove_relations_to_deleted_entities_helper();

    apply_induced_part_membership_helper();

    std::vector<partition_id> new_partition_order;
    compress_data_structure_helper(new_partition_order);

    populate_rank_index_helper();

    // compress partitions

    for (size_t p = 0, end_partition = m_partitions.size(); p < end_partition; ++p) {
      m_partitions[p]->end_modification_impl(); // compress
    }

    // recycle id's of deleted keys

    for (interval_set<entity_key>::const_iterator itr = m_deleted_keys.begin(), end = m_deleted_keys.end();
         itr != end;
         ++itr) {
      entity_key deleted_key = *itr;
      if (deleted_key.process() != m_process) {
        m_non_local_key_to_descriptor_map.erase(deleted_key);
      }
      else {
        m_local_ids[deleted_key.topology()()] -= deleted_key.local_id();
      }
    }

    //update max_local_ids
    for (entity_topology t={0}, e=entity_topology::create(m_local_ids.size()); t<e; ++t) {
      max_local_id_changed_helper(t, max_local_id(t));
    }

    // ==================================================
    // wrap-up
    // ==================================================

    m_is_in_modification_cycle = false;

    // signal fields
    m_signals.reorder_partitions_signal(new_partition_order);
    m_signals.end_modification_signal();

    // TODO - Really could use some invariant-checking now that we have
    // set up all these complex data structures!

    return true;
  }
  else {
    return false;
  }
}

inline entity_block_key mesh_impl::add_entity_block(std::string const& name,
                                                    entity_rank rank,
                                                    bool inducable)
{
  BOOST_ASSERT_MSG(m_is_in_modification_cycle, "Not in modification cycle");

  entity_block_key s = find_entity_block(name);
  if (s != entity_block_key::invalid())
  {
    BOOST_ASSERT_MSG(m_set_specs[s()].rank == rank,
                     (debug_message() << "Cannot change rank of entity_block to " << rank << ", was " << m_set_specs[s()].rank));
    BOOST_ASSERT_MSG(m_set_specs[s()].inducable == inducable,
                     (debug_message() << "Cannot change inducability of entity_block to " << inducable));
    return s;
  }
  //create set
  s = entity_block_key::create(m_set_specs.size());
  int inducable_index = inducable ? static_cast<int>(m_inducable_parts.size()) : -1;
  entity_block_spec spec = {name, rank, inducable, inducable_index};
  m_set_specs.push_back(spec);
  if (inducable) {
    m_inducable_parts.push_back(s);
  }
  m_signals.add_entity_block_signal(s);
  return s;
}


template <typename EntitySetIterator>
inline
interval<entity_key> mesh_impl::add_entities( entity_topology topology
                                             ,size_t how_many
                                             ,EntitySetIterator first
                                             ,EntitySetIterator last
                                             ,bool is_owned
                                            )
{
  BOOST_ASSERT_MSG(m_is_in_modification_cycle, "Not in modification cycle");
  validate_helper(topology);
  validate_helper(first, last);

  entity_rank rank = topology_rank(topology,spatial_dimension());

  // Compute part-vector for new entities
  entity_part_vector parts;
  {
    parts.reserve(4 + last-first);

    parts.push_back(rank);
    parts.push_back(topology);
    parts.push_back(entity_state::universe());
    if (is_owned) parts.push_back(entity_state::owned());
    parts.insert(parts.end(),first,last);

    std::sort(parts.begin(),parts.end());
  }

  // Find/create partition_id for this part vector
  partition_id sub = get_partition_non_const(parts);

  // Create keys for the new entities
  entity_key_interval keys = create_keys_helper(topology, how_many);

  partition_offset first_offset = partition_offset::create(m_descriptor_to_key_map[sub()].size());
  partition_index first_descriptor = partition_index::create(rank, sub, first_offset);
  partition_index last_descriptor  = partition_index::create(rank, sub, first_offset + how_many);
  interval<partition_index> descriptors(first_descriptor, last_descriptor);

  m_descriptor_to_key_map[sub()].insert( m_descriptor_to_key_map[sub()].end()
                                        ,keys.begin()
                                        ,keys.end()
                                       );
  for (size_t i = 0; i < how_many; ++i) {
    m_key_to_descriptor_map[topology()][keys[i].local_id()()] = descriptors[i];
  }

  m_partitions[sub()]->add_entities_impl(how_many);

  m_signals.add_entity_keys_signal(keys);
  m_signals.add_entities_signal(sub,how_many);

  return keys;
}

template <typename EntityKeyIterator>
inline
void mesh_impl::remove_entities( EntityKeyIterator first
                                ,EntityKeyIterator last
                               )
{
  BOOST_ASSERT_MSG(m_is_in_modification_cycle, "Not in modification cycle");
  validate_helper(first, last);

  if( first == last ) return;

  //copy to a vector
  std::vector<entity_key> removed_keys(last-first);
  std::copy(first,last,removed_keys.begin());
  entity_key_iterator kbegin = &removed_keys[0];
  entity_key_iterator kend   = kbegin + (last - first);

  // convert to partition_indexs
  std::vector<partition_index> descriptors;
  descriptors.reserve(removed_keys.size());

  for(entity_key_iterator itr = kbegin; itr != kend; ++itr) {
    descriptors.push_back(convert(*itr));
    m_deleted_keys += *itr;
  }

  boost::sort(descriptors);

  // Iterate over all removed descriptors. For each batch of descriptors
  // in the same partition_id, remove them from the partition_id all at once.

  std::vector<partition_index>::iterator partition_first = descriptors.begin();
  std::vector<partition_index>::iterator partition_last  = descriptors.begin();
  std::vector<partition_index>::iterator end             = descriptors.end();

  while (partition_last != end) {
    partition_first = partition_last;
    partition_id current_partition = partition_first->partition();

    // iterate until we run out of descriptors in this partition
    for (; partition_last != end && partition_last->partition() == current_partition; ++partition_last);

    remove_entities_helper(current_partition, partition_first, partition_last);
  }

  for(entity_key_iterator itr = kbegin; itr != kend; ++itr) {
    m_key_to_descriptor_map[itr->topology()()][itr->local_id()()] = partition_index::invalid();
  }

  m_signals.remove_entity_keys_signal(kbegin, kend);
}

template <typename EntityKeyIterator, typename EntitySetIterator, typename RemoveEntitySetIterator>
inline
void mesh_impl::move_entities( EntityKeyIterator first_key
                              ,EntityKeyIterator last_key
                              ,EntitySetIterator first_add_set
                              ,EntitySetIterator last_add_set
                              ,RemoveEntitySetIterator first_remove_set
                              ,RemoveEntitySetIterator last_remove_set
                             )
{
  BOOST_ASSERT_MSG(m_is_in_modification_cycle, "Not in modification cycle");
  validate_helper(first_key, last_key);
  validate_helper(first_add_set, last_add_set);
  validate_helper(first_remove_set, last_remove_set);

  // nothing to do
  if(  (first_key == last_key) ||
     ((first_add_set == last_add_set) && (first_remove_set==last_remove_set))
    )
  {
    return;
  }

  // convert to partition_indexs
  std::vector<partition_index> descriptors;
  descriptors.reserve(last_key-first_key);
  for(; first_key != last_key; ++first_key)
    descriptors.push_back(convert(*first_key));

  // sort descriptors, this will group-together descriptors with the same partition
  boost::sort(descriptors);

  // maintain two iterators and walk.

  std::vector<partition_index>::iterator       last_entity_in_partition = descriptors.begin();
  const std::vector<partition_index>::iterator end_of_entities          = descriptors.end();

  while (last_entity_in_partition != end_of_entities) {
    const std::vector<partition_index>::iterator first_entity_in_partition = last_entity_in_partition;
    partition_id current_partition = first_entity_in_partition->partition();

    //get the parts for the new partition
    entity_part_vector all_parts;

    boost::copy(m_partitions[current_partition()]->parts(), std::back_inserter(all_parts));

    all_parts.insert(all_parts.end(),first_add_set,last_add_set);

    entity_part_vector remove_parts;
    remove_parts.reserve(last_remove_set - first_remove_set);
    remove_parts.insert(remove_parts.end(), first_remove_set, last_remove_set);

    boost::sort(all_parts);

    entity_part_vector::iterator unique_end = std::unique(all_parts.begin(), all_parts.end());

    boost::sort(remove_parts);

    entity_part_vector parts;
    std::set_difference( all_parts.begin(), unique_end
                        ,remove_parts.begin(), remove_parts.end()
                        ,std::back_inserter(parts)
                       );

    partition_id new_partition = get_partition_non_const(parts);

    //advance to end of current partition
    for (; last_entity_in_partition != end_of_entities && last_entity_in_partition->partition() == current_partition;
         ++last_entity_in_partition);

    if (current_partition != new_partition) {
      move_entities_helper(current_partition, new_partition, first_entity_in_partition, last_entity_in_partition);
    }
  }
}

inline void mesh_impl::add_connectivity( entity_key from
                                        ,entity_key to
                                        ,connectivity_ordinal ordinal
                                        ,connectivity_orientation orientation
                                        )
{
  BOOST_ASSERT_MSG(m_is_in_modification_cycle, "Not in modification cycle");
  validate_helper(from);
  validate_helper(to);

#ifndef NDEBUG
  entity_rank from_rank = topology_rank(from.topology(), spatial_dimension());
  entity_rank to_rank   = topology_rank(to.topology(),   spatial_dimension());
#endif
  BOOST_ASSERT_MSG(m_connectivity_map(from_rank, to_rank) != connectivity_kind::invalid(),
                   (debug_message() << "Relations from " << from_rank << " to " << to_rank << " not supported"));

  partition_index from_descriptor = convert(from);

  m_partitions[from_descriptor.partition()()]->add_connectivity( from_descriptor.offset(), to, ordinal, orientation);
}

inline void mesh_impl::remove_connectivity( entity_key from
                                           ,entity_key to
                                           ,connectivity_ordinal ordinal
                                           )
{
  BOOST_ASSERT_MSG(m_is_in_modification_cycle, "Not in modification cycle");
  validate_helper(from);
  validate_helper(to);

#ifndef NDEBUG
  entity_rank from_rank = topology_rank(from.topology(), spatial_dimension());
  entity_rank to_rank   = topology_rank(to.topology(),   spatial_dimension());
#endif
  BOOST_ASSERT_MSG(m_connectivity_map(from_rank, to_rank) != connectivity_kind::invalid(),
                   (debug_message() << "Relations from " << from_rank << " to " << to_rank << " not supported"));

  partition_index from_descriptor = convert(from);

  m_partitions[from_descriptor.partition()()]->remove_connectivity( from_descriptor.offset(), to, ordinal);
}


//*****************************************************************************
//internal helpers
//*****************************************************************************

inline void mesh_impl::add_topology_helper(entity_topology topology)
{
  validate_helper(topology);

  //is this a new topology ?
  if(m_local_ids.size() <= topology) {
    m_local_ids.resize(topology()+1);
    m_key_to_descriptor_map.resize(topology()+1);
    // TODO Should this event signal???
  }
}

inline void mesh_impl::max_local_id_changed_helper( entity_topology topology
                                                    ,entity_local_id max_id
                                                    )
{
  validate_helper(topology);

  std::vector<partition_index> & descriptors = m_key_to_descriptor_map[topology()];
  if(descriptors.size() != max_id()) {
    descriptors.resize(max_id(), partition_index::invalid());
  }
}

inline
interval<entity_key> mesh_impl::create_keys_helper(entity_topology topology, size_t how_many)
{
  validate_helper(topology);

  add_topology_helper(topology);

  interval<entity_local_id> ids = get_local_ids_helper(topology,how_many);

  max_local_id_changed_helper(topology,ids.upper());

  interval<entity_key> keys =
    interval<entity_key>( entity_key::create(topology, m_process, ids.lower() )
                         ,entity_key::create(topology, m_process ,ids.upper() )
                        );
  return keys;
}

inline
interval<entity_local_id> mesh_impl::get_local_ids_helper(entity_topology topology, size_t how_many)
{
  validate_helper(topology);

  // The local-id interval for the new entities
  interval<entity_local_id> ids;

  // The local-id interval-set that the new ids must be added to
  interval_set<entity_local_id> & id_set = m_local_ids[topology()];

  if ( id_set.empty() || ( id_set.interval_begin()->lower()() >= how_many ) ) {
    // The new ids can fit into the beginning of the interval-set
    ids = interval<entity_local_id>(0, how_many);
  }
  else {
    interval_set<entity_local_id>::const_interval_iterator curr = id_set.interval_begin();
    interval_set<entity_local_id>::const_interval_iterator next = id_set.interval_begin();
    interval_set<entity_local_id>::const_interval_iterator end  = id_set.interval_end();

    // Scan the interval-set, looking for an open interval large enough to hold
    // the local-id interval for the new entities. This implementation allows the
    // local-ids to become fragmented, but we believe it's worth it to allow
    // for created batches of entities to be contiguous.
    for(++next;next != end && ((next->lower()() - curr->upper()()) < how_many); ++curr,++next);

    uint32_t value = curr->upper()();
    ids = interval<entity_local_id>(value,value+how_many);
  }

  id_set += ids;

  return ids;
}

// Note: the first and second entry in the part vector should always be a
// an entity_rank and an entity_topology respectively
// TODO throw error if part.size()<2 || part[0] != entity_rank || part[1] != entity_topology
inline partition_id mesh_impl::get_partition_non_const( entity_part_vector const& parts )
{
  validate_helper(parts.begin(), parts.end());

  part_vector_to_partition_map::const_iterator itr = m_parts_to_partition.find(parts);
  if (itr != m_parts_to_partition.end()) {
    return itr->second;
  }
  return create_partition(parts);
}

inline partition_id mesh_impl::create_partition( entity_part_vector const& parts )
{
  validate_helper(parts.begin(), parts.end());

  partition_id sub = partition_id::create(m_partitions.size());

  m_descriptor_to_key_map.resize(sub()+1);

  m_partitions.push_back( partition_handle(new partition_impl(this,sub,parts)) );

  m_parts_to_partition[parts] = sub;

  m_signals.add_partition_signal(sub);

  return sub;
}

template <typename Iterator>
inline void mesh_impl::move_entities_helper( partition_id from
                                            ,partition_id to
                                            ,Iterator first
                                            ,Iterator last
                                            ,bool during_end_modification
                                           )
{
  validate_helper(from);
  validate_helper(to);
  validate_helper(first, last);

  if (from == to) return; // nothing to do!

  const size_t how_many = last-first;

  if (how_many == 0) return;

  entity_rank rank = first->rank();

#ifndef NDEBUG
  // If this call did not come from end_modification, make sure user is not trying
  // to put entities into a part of differing rank
  if (!during_end_modification) {
    entity_part_vector const& parts = m_partitions[to()]->parts();
    for (entity_part_vector::const_iterator itr = parts.begin(), end = parts.end(); itr != end; ++itr) {
      if (itr->which() == entity_part::Topology) {
        BOOST_ASSERT_MSG(topology_rank(itr->topology(), spatial_dimension()) == rank,
                         (debug_message() << "Cannot add entity of rank " << rank << " to topology " << itr->topology()));
      }
      else if (itr->which() == entity_part::Set) {
        entity_block_key eset = itr->set();
        entity_block_spec eset_spec = m_set_specs[eset()];
        if (eset_spec.rank != entity_rank::invalid()) {
          BOOST_ASSERT_MSG(eset_spec.rank == rank,
                           (debug_message() << "Cannot add entity of rank " << rank <<
                            " to set " << eset_spec.name << " which has rank " << eset_spec.rank));
        }
      }
    }
  }

  // All entities being moved should be of the same rank
  for (Iterator itr = first; itr != last; ++itr) {
    BOOST_ASSERT_MSG(itr->rank() == rank,
                     (debug_message() << "entity " << *itr << " does not have rank matching " << rank));
  }
#endif

  std::vector<entity_key> & from_map = m_descriptor_to_key_map[from()];
  std::vector<entity_key> & to_map = m_descriptor_to_key_map[to()];

  std::vector<partition_offset> move_offsets;
  move_offsets.reserve(last-first);

  for(;first!=last; ++first) {
    move_offsets.push_back(first->offset());
  }

  boost::sort(move_offsets,std::greater<partition_offset>());

  partition_offset last_offset = partition_offset::create( from_map.size()-1);

  partition_handle from_impl = m_partitions[from()];

  //swap entities to end
  for(std::vector<partition_offset>::const_iterator itr=move_offsets.begin(), end_itr=move_offsets.end();
      itr != end_itr;
      ++itr, --last_offset
     )
  {
    if(*itr != last_offset) {
      entity_key & a = from_map[(*itr)()];
      entity_key & b = from_map[last_offset()];
      m_key_to_descriptor_map[a.topology()()][a.local_id()()] = partition_index::create(rank,from,last_offset);
      m_key_to_descriptor_map[b.topology()()][b.local_id()()] = partition_index::create(rank,from,*itr);
      std::swap(a, b);
      from_impl->swap(*itr, last_offset);
      m_signals.swap_offset_signal(from, *itr, last_offset);
    }
  }
  //since last_offset wraps an unsigned it is possible for it to be max unsigned
  //if moving the entity partition.  increment he to avoid integer overflow
  ++last_offset;

  partition_offset to_offset = partition_offset::create(to_map.size());

  //update d2k map
  to_map.insert(to_map.end(), from_map.begin()+last_offset(), from_map.end());

  //update k2d map
  for (partition_offset i=to_offset,e=to_offset+how_many; i<e; ++i)
  {
    entity_key k = to_map[i()];
    m_key_to_descriptor_map[k.topology()()][k.local_id()()] = partition_index::create(rank,to,i);
  }

  from_map.erase(from_map.begin()+last_offset(), from_map.end());

  from_impl->move_entities_impl(*m_partitions[to()], how_many);

  m_signals.move_entities_signal(from,to,how_many);
}

template <typename Iterator>
inline void mesh_impl::remove_entities_helper( partition_id from
                                              ,Iterator first
                                              ,Iterator last
                                             )
{
  validate_helper(from);
  validate_helper(first, last);

  const size_t how_many = last-first;

  if (how_many == 0) return;

  entity_rank rank = first->rank();

  std::vector<entity_key> & from_map = m_descriptor_to_key_map[from()];

  std::vector<partition_offset> move_offsets;
  move_offsets.reserve(last-first);

  for(;first!=last; ++first) {
    move_offsets.push_back(first->offset());
  }

  boost::sort(move_offsets,std::greater<partition_offset>());

  partition_offset last_offset = partition_offset::create( from_map.size()-1);

  partition_handle from_impl = m_partitions[from()];

  //swap entities to end
  for(std::vector<partition_offset>::const_iterator itr=move_offsets.begin(), end_itr=move_offsets.end();
      itr != end_itr;
      ++itr, --last_offset
     )
  {
    if(*itr != last_offset) {
      entity_key & a = from_map[(*itr)()];
      entity_key & b = from_map[last_offset()];
      m_key_to_descriptor_map[a.topology()()][a.local_id()()] = partition_index::create(rank,from,last_offset);
      m_key_to_descriptor_map[b.topology()()][b.local_id()()] = partition_index::create(rank,from,*itr);
      std::swap(a, b);
      from_impl->swap(*itr, last_offset);
      m_signals.swap_offset_signal(from, *itr, last_offset);
    }
  }
  //since last_offset wraps an unsigned it is possible for it to be max unsigned
  //if moving the entity partition.  increment he to avoid integer overflow
  ++last_offset;

  from_map.erase(from_map.begin()+last_offset(), from_map.end());

  from_impl->remove_entities_impl(how_many);

  m_signals.remove_entities_signal(from,how_many);
}

inline
size_t mesh_impl::find_induced_parts_helper(entity_part_vector const& orig_parts) const
{
  validate_helper(orig_parts.begin(), orig_parts.end());

  entity_part_vector::const_iterator
    begin             = orig_parts.begin(),
    end               = orig_parts.end(),
    first_induced_idx = begin;
  for ( ; first_induced_idx != end && !first_induced_idx->induced(); ++first_induced_idx);

  return first_induced_idx - begin;
}

inline
void mesh_impl::add_induced_parts_helper(entity_part_vector const& parent_parts,
                                         const size_t first_induced_idx,
                                         induced_part_set& induced_parts)
{
  validate_helper(parent_parts.begin(), parent_parts.end());

  for (size_t i = 0; i < first_induced_idx; ++i) {
    entity_part part = parent_parts[i];
    if (part.which() == entity_part::Topology) {
      induced_parts.set(part.topology()());
    }
    else if (part.which() == entity_part::Set && m_set_specs[part.set()()].inducable) {
      induced_parts.set(entity_topology::invalid()() + m_set_specs[part.set()()].inducable_index);
    }
  }
}

inline
void mesh_impl::induce_parts_through_downward_connectivity_helper(partition_handle                              p_handle,
                                                                  std::vector< std::vector<induced_part_set> >& key_to_induced_part_map)
{
  entity_rank partition_rank = p_handle->rank();
  BOOST_ASSERT(partition_rank != entity_rank::node());

  entity_part_vector const& parent_parts = p_handle->parts();
  const size_t first_induced_idx = find_induced_parts_helper(parent_parts);
  entity_rank r = partition_rank;
  do {
    --r;
    if (connectivity_map()(partition_rank, r) != connectivity_map::invalid()) {
      for (partition_offset o = {0}; o < p_handle->size(); ++o) {
        induce_parts_helper(parent_parts,
                            first_induced_idx,
                            p_handle->begin_connectivity<entity_key>(r, o),
                            p_handle->end_connectivity<entity_key>(r, o),
                            key_to_induced_part_map);
      }
    }
  } while (r != entity_rank::node());
}

template <typename RelationIterator>
inline
void mesh_impl::induce_parts_helper(entity_part_vector const&                     parent_parts,
                                    const size_t                                  first_induced_idx,
                                    RelationIterator                              begin,
                                    RelationIterator                              end,
                                    std::vector< std::vector<induced_part_set> >& key_to_induced_part_map)
{
  for (RelationIterator rel_itr = begin; rel_itr != end; ++rel_itr) {
    entity_key to = *rel_itr;
    if (to != entity_key::invalid()) {
      induced_part_set& induced_parts = key_to_induced_part_map[to.topology()()][to.local_id()()];
      add_induced_parts_helper(parent_parts, first_induced_idx, induced_parts);
    }
  }
}

// get_rank_index_map specializations

template <>
inline
const std::vector<partition_index>& mesh_impl::get_rank_index_map<entity_rank::node_type>() const
{ return m_node_index_to_partition_index_map; }

template <>
inline
const std::vector<partition_index>& mesh_impl::get_rank_index_map<entity_rank::edge_type>() const
{ return m_edge_index_to_partition_index_map; }

template <>
inline
const std::vector<partition_index>& mesh_impl::get_rank_index_map<entity_rank::face_type>() const
{ return m_face_index_to_partition_index_map; }

template <>
inline
const std::vector<partition_index>& mesh_impl::get_rank_index_map<entity_rank::element_type>() const
{ return m_element_index_to_partition_index_map; }

} } //namespace samba::detail

#endif //SAMBA_SAMBA_MESH_MESH_IMPL_TCC
