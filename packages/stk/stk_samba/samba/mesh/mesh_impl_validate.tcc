#ifndef SAMBA_SAMBA_MESH_MESH_IMPL_VALIDATE_TCC
#define SAMBA_SAMBA_MESH_MESH_IMPL_VALIDATE_TCC

namespace samba { namespace detail {

inline
void mesh_impl::validate_helper(partition_id p) const
{
  BOOST_ASSERT_MSG(p != partition_id::invalid(), "invalid");
  BOOST_ASSERT_MSG(p() < m_partitions.size(),
                   (debug_message() << "partition_id " << p() << " out-of-bounds, max is: " << m_partitions.size() - 1));
}

inline
void mesh_impl::validate_helper(partition_index d) const
{
  BOOST_ASSERT_MSG(d != partition_index::invalid(), "invalid");
  BOOST_ASSERT_MSG(d.partition()() < m_partitions.size(),
                   (debug_message() << "invalid partition_index, has partition " << d.partition() <<
                    " that is out-of-bounds, max is: " << m_partitions.size() - 1));
  BOOST_ASSERT_MSG(d.offset()() < m_partitions[d.partition()()]->size(),
                   (debug_message() << "invalid partition_index, has offset " << d.offset() <<
                    " that is out-of-bounds, max is: " << m_partitions[d.partition()()]->size() - 1));
  BOOST_ASSERT_MSG(d.partition()() < m_descriptor_to_key_map.size(),
                   (debug_message() << "m_descriptor_to_key_map is wrong, missing partition " << d.partition()));
  BOOST_ASSERT_MSG(d.offset()() < m_descriptor_to_key_map[d.partition()()].size(),
                   (debug_message() << "m_descriptor_to_key_map is wrong, missing offset " << d.offset()));
  BOOST_ASSERT_MSG(d == convert_helper<partition_index>(convert_helper<entity_key>(d)),
                   "Problem with desc<->key mapping");
}

inline
void mesh_impl::validate_helper(entity_key k) const
{
  BOOST_ASSERT_MSG(k != entity_key::invalid(), "invalid");
  BOOST_ASSERT_MSG(k.topology()() < m_key_to_descriptor_map.size(),
                   (debug_message() << "m_key_to_descriptor_map is wrong, missing topology " << k.topology()));
  BOOST_ASSERT_MSG(k.local_id()() < m_key_to_descriptor_map[k.topology()()].size(),
                   (debug_message() << "m_key_to_descriptor_map is wrong, missing local_id " << k.local_id()));
//  BOOST_ASSERT_MSG(k.process() == process_id::invalid(),
//                   (debug_message() << "Bad process id " << k.process()));
  BOOST_ASSERT_MSG(!contains(m_deleted_keys, k),
                   (debug_message() << "Trying to perform operations with deleted key: " << k));
  BOOST_ASSERT_MSG(k == convert_helper<entity_key>(convert_helper<partition_index>(k)),
                   "Problem with desc<->key mapping");
}

inline
void mesh_impl::validate_helper(entity_topology t) const
{
  BOOST_ASSERT_MSG(t != entity_topology::invalid(), "invalid");
}

inline
void mesh_impl::validate_helper(entity_block_key k) const
{
  BOOST_ASSERT_MSG(k != entity_block_key::invalid(), "invalid");
  BOOST_ASSERT_MSG(k() < m_set_specs.size(),
                   (debug_message() << "invalid entity_block_key, has value " << k() << ", max is: " << m_set_specs.size() - 1));
}

inline
void mesh_impl::validate_helper(entity_rank r) const
{
  BOOST_ASSERT_MSG(r != entity_rank::invalid(), "invalid");
}

inline
void mesh_impl::validate_helper(entity_state s) const
{
  BOOST_ASSERT_MSG(s != entity_state::invalid(), "invalid");
}

inline
void mesh_impl::validate_helper(entity_part p) const
{
  BOOST_ASSERT_MSG(p.which() != entity_part::Invalid, "invalid");
#ifndef NDEBUG
  switch (p.which()) {
  case entity_part::Rank:
    validate_helper(p.rank());
    break;
  case entity_part::Topology:
    validate_helper(p.topology());
    break;
  case entity_part::State:
    validate_helper(p.state());
    break;
  case entity_part::Set:
    validate_helper(p.set());
    break;
  default:
    BOOST_ASSERT_MSG(false, (debug_message() << "Unhandled entity_part type: " << p.which()));
  }
#endif
}

inline
void mesh_impl::validate_helper(node_index i) const
{
  BOOST_ASSERT_MSG(i != node_index::invalid(), "invalid");
  BOOST_ASSERT_MSG(i < m_node_index_to_partition_index_map.size(),
                   (debug_message() << i << " is beyond max value " << m_node_index_to_partition_index_map.size()));
  BOOST_ASSERT_MSG(i == convert_helper<node_index>(convert_helper<partition_index>(i)),
                   "Problem with desc<->index mapping");
}

inline
void mesh_impl::validate_helper(edge_index i) const
{
  BOOST_ASSERT_MSG(i != edge_index::invalid(), "invalid");
  BOOST_ASSERT_MSG(i < m_edge_index_to_partition_index_map.size(),
                   (debug_message() << i << " is beyond max value " << m_edge_index_to_partition_index_map.size()));
  BOOST_ASSERT_MSG(i == convert_helper<edge_index>(convert_helper<partition_index>(i)),
                   "Problem with desc<->index mapping");
}

inline
void mesh_impl::validate_helper(face_index i) const
{
  BOOST_ASSERT_MSG(i != face_index::invalid(), "invalid");
  BOOST_ASSERT_MSG(i < m_face_index_to_partition_index_map.size(),
                   (debug_message() << i << " is beyond max value " << m_face_index_to_partition_index_map.size()));
  BOOST_ASSERT_MSG(i == convert_helper<face_index>(convert_helper<partition_index>(i)),
                   "Problem with desc<->index mapping");
}

inline
void mesh_impl::validate_helper(element_index i) const
{
  BOOST_ASSERT_MSG(i != element_index::invalid(), "invalid");
  BOOST_ASSERT_MSG(i < m_element_index_to_partition_index_map.size(),
                   (debug_message() << i << " is beyond max value " << m_element_index_to_partition_index_map.size()));
  BOOST_ASSERT_MSG(i == convert_helper<element_index>(convert_helper<partition_index>(i)),
                   "Problem with desc<->index mapping");
}

template <typename IteratorT>
inline
void mesh_impl::validate_helper(IteratorT first, IteratorT last) const
{
#ifndef NDEBUG
  for (IteratorT itr = first; itr != last; ++itr) {
    validate_helper(*itr);
  }
#endif
}


}} //namespace samba::detail

#endif //SAMBA_SAMBA_MESH_MESH_IMPL_VALIDATE_TCC

