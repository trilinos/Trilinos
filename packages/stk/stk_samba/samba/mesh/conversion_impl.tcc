namespace samba { namespace detail {

template <typename MeshImpl>
inline partition_index conversion_impl<MeshImpl>::convert(entity_key k) const
{
  mesh_impl().validate_helper(k);
  return mesh_impl().convert_helper<partition_index>(k);
}

template <typename MeshImpl>
inline entity_key conversion_impl<MeshImpl>::convert(partition_index d) const
{
  mesh_impl().validate_helper(d);
  return mesh_impl().convert_helper<entity_key>(d);
}

template <typename MeshImpl>
template <typename ToIndex, typename FromIndex>
inline
ToIndex conversion_impl<MeshImpl>::convert(FromIndex from) const
{
  mesh_impl().validate_helper(from);
  return mesh_impl().convert_helper<ToIndex>(from);
}

//*************************************************************************
// Conversion specializations
//*************************************************************************

template <>
template <>
inline
partition_index conversion_impl<mesh_impl>::convert_helper<partition_index>(entity_key k) const
{
  if (k.process() == mesh_impl().m_process) {
    BOOST_ASSERT_MSG(mesh_impl().m_key_to_descriptor_map.size() > k.topology()(),
                     (debug_message() << "Key " << k << " has an undefined topology"));
    BOOST_ASSERT_MSG(mesh_impl().m_key_to_descriptor_map[k.topology()()].size() > k.local_id()(),
                     (debug_message() << "Key " << k << " has out of bounds local id"));

    return mesh_impl().m_key_to_descriptor_map[k.topology()()][k.local_id()()];
  }
  else {
    return mesh_impl().m_non_local_key_to_descriptor_map.at(k);
  }
}

template <>
template <>
inline
entity_key conversion_impl<mesh_impl>::convert_helper<entity_key>(partition_index d) const
{
  BOOST_ASSERT_MSG(d.partition() < mesh_impl().m_descriptor_to_key_map.size(),
                   (debug_message() << "Key " << d << " has out-of-bounds partition " << d.partition()));
  BOOST_ASSERT_MSG(d.offset() < mesh_impl().m_descriptor_to_key_map[d.partition()()].size(),
                   (debug_message() << "Key " << d << " has out-of-bounds offset " << d.offset()));

  return mesh_impl().m_descriptor_to_key_map[d.partition()()][d.offset()()];
}

// rank indexes to partion_index conversions...

template <>
template <>
inline
partition_index conversion_impl<mesh_impl>::convert_helper<partition_index>(node_index from) const
{
  BOOST_ASSERT_MSG(from < mesh_impl().m_node_index_to_partition_index_map.size(),
                   (debug_message() << "Rank index " << from << " is out-of-bounds, max is " << mesh_impl().m_node_index_to_partition_index_map.size()));

  return mesh_impl().m_node_index_to_partition_index_map[from()];
}

template <>
template <>
inline
partition_index conversion_impl<mesh_impl>::convert_helper<partition_index>(edge_index from) const
{
  BOOST_ASSERT_MSG(from < mesh_impl().m_edge_index_to_partition_index_map.size(),
                   (debug_message() << "Rank index " << from << " is out-of-bounds, max is " << mesh_impl().m_edge_index_to_partition_index_map.size()));

  return mesh_impl().m_edge_index_to_partition_index_map[from()];
}

template <>
template <>
inline
partition_index conversion_impl<mesh_impl>::convert_helper<partition_index>(face_index from) const
{
  BOOST_ASSERT_MSG(from < mesh_impl().m_face_index_to_partition_index_map.size(),
                   (debug_message() << "Rank index " << from << " is out-of-bounds, max is " << mesh_impl().m_face_index_to_partition_index_map.size()));

  return mesh_impl().m_face_index_to_partition_index_map[from()];
}

template <>
template <>
inline
partition_index conversion_impl<mesh_impl>::convert_helper<partition_index>(element_index from) const
{
  BOOST_ASSERT_MSG(from < mesh_impl().m_element_index_to_partition_index_map.size(),
                   (debug_message() << "Rank index " << from << " is out-of-bounds, max is " << mesh_impl().m_element_index_to_partition_index_map.size()));

  return mesh_impl().m_element_index_to_partition_index_map[from()];
}

// rank indexes to entity_key conversions...

template <>
template <>
inline entity_key conversion_impl<mesh_impl>::convert_helper<entity_key>(node_index from) const
{ return convert_helper<entity_key>(convert_helper<partition_index>(from)); }

template <>
template <>
inline entity_key conversion_impl<mesh_impl>::convert_helper<entity_key>(edge_index from) const
{ return convert_helper<entity_key>(convert_helper<partition_index>(from)); }

template <>
template <>
inline entity_key conversion_impl<mesh_impl>::convert_helper<entity_key>(face_index from) const
{ return convert_helper<entity_key>(convert_helper<partition_index>(from)); }

template <>
template <>
inline entity_key conversion_impl<mesh_impl>::convert_helper<entity_key>(element_index from) const
{ return convert_helper<entity_key>(convert_helper<partition_index>(from)); }


// partition_index to rank index conversions...

template <>
template <>
inline
node_index conversion_impl<mesh_impl>::convert_helper<node_index>(partition_index from) const
{
  mesh_impl::partition_handle partition = mesh_impl().m_partitions[from.partition()()];

  BOOST_ASSERT_MSG(partition->rank() == entity_rank::node(),
                   (debug_message() << "Invalid conversion; trying to get node_index for entity of rank " << partition->rank()));

  return node_index::create(partition->starting_rank_count() + from.offset()());
}

template <>
template <>
inline
edge_index conversion_impl<mesh_impl>::convert_helper<edge_index>(partition_index from) const
{
  mesh_impl::partition_handle partition = mesh_impl().m_partitions[from.partition()()];

  BOOST_ASSERT_MSG(partition->rank() == entity_rank::edge(),
                   (debug_message() << "Invalid conversion; trying to get edge_index for entity of rank " << partition->rank()));

  return edge_index::create(partition->starting_rank_count() + from.offset()());
}

template <>
template <>
inline
face_index conversion_impl<mesh_impl>::convert_helper<face_index>(partition_index from) const
{
  mesh_impl::partition_handle partition = mesh_impl().m_partitions[from.partition()()];

  BOOST_ASSERT_MSG(partition->rank() == entity_rank::face(),
                   (debug_message() << "Invalid conversion; trying to get face_index for entity of rank " << partition->rank()));

  return face_index::create(partition->starting_rank_count() + from.offset()());
}

template <>
template <>
inline
element_index conversion_impl<mesh_impl>::convert_helper<element_index>(partition_index from) const
{
  mesh_impl::partition_handle partition = mesh_impl().m_partitions[from.partition()()];

  BOOST_ASSERT_MSG(partition->rank() == entity_rank::element(),
                   (debug_message() << "Invalid conversion; trying to get element_index for entity of rank " << partition->rank()));

  return element_index::create(partition->starting_rank_count() + from.offset()());
}

// entity_key to rank index conversions...

template <>
template <>
inline node_index conversion_impl<mesh_impl>::convert_helper<node_index>(entity_key from) const
{ return convert_helper<node_index>(convert_helper<partition_index>(from)); }

template <>
template <>
inline edge_index conversion_impl<mesh_impl>::convert_helper<edge_index>(entity_key from) const
{ return convert_helper<edge_index>(convert_helper<partition_index>(from)); }

template <>
template <>
inline face_index conversion_impl<mesh_impl>::convert_helper<face_index>(entity_key from) const
{ return convert_helper<face_index>(convert_helper<partition_index>(from)); }

template <>
template <>
inline element_index conversion_impl<mesh_impl>::convert_helper<element_index>(entity_key from) const
{ return convert_helper<element_index>(convert_helper<partition_index>(from)); }

} } // namespace samba::detail
