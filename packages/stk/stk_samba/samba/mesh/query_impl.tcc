namespace samba { namespace detail {

template <typename MeshImpl>
inline
entity_block_key query_impl<MeshImpl>::find_entity_block(std::string const& name) const
{
  typename std::vector<typename MeshImpl::entity_block_spec>::const_iterator itr = std::find_if(mesh_impl().m_set_specs.begin(), mesh_impl().m_set_specs.end(), typename MeshImpl::entity_block_match_name(name));
  entity_block_key set = entity_block_key::invalid();
  if (itr != mesh_impl().m_set_specs.end())
    set = entity_block_key::create(itr - mesh_impl().m_set_specs.begin());
  return set;
}

template <typename MeshImpl>
inline
entity_rank query_impl<MeshImpl>::get_rank(entity_block_key key) const
{
  BOOST_ASSERT_MSG(key() < mesh_impl().m_set_specs.size(),
                   (debug_message() << "entity_block_key " << key() << " out-of-bounds, max is: " << mesh_impl().m_set_specs.size() - 1));
  return mesh_impl().m_set_specs[key()].rank;
}

template <typename MeshImpl>
inline
const std::string& query_impl<MeshImpl>::get_name(entity_block_key key) const
{
  BOOST_ASSERT_MSG(key() < mesh_impl().m_set_specs.size(),
                   (debug_message() << "entity_block_key " << key() << " out-of-bounds, max is: " << mesh_impl().m_set_specs.size() - 1));
  return mesh_impl().m_set_specs[key()].name;
}

template <typename MeshImpl>
inline
size_t query_impl<MeshImpl>::num_entity_blocks() const
{ return mesh_impl().m_set_specs.size(); }

template <typename MeshImpl>
inline
size_t query_impl<MeshImpl>::num_partitions() const
{ return mesh_impl().m_partitions.size(); }

template <typename MeshImpl>
inline
process_id query_impl<MeshImpl>::process() const
{ return mesh_impl().m_process; }

template <typename MeshImpl>
inline
entity_local_id query_impl<MeshImpl>::max_local_id(entity_topology topology) const
{
  mesh_impl().validate_helper(topology);
  if (mesh_impl().m_local_ids.size() > topology() )
    return mesh_impl().m_local_ids[topology()].upper();
  return entity_local_id::create(0);
}

template <typename MeshImpl>
inline
partition_id query_impl<MeshImpl>::get_partition( entity_part_vector const& parts ) const
{
  mesh_impl().validate_helper(parts.begin(), parts.end());

  typename MeshImpl::part_vector_to_partition_map::const_iterator itr = mesh_impl().m_parts_to_partition.find(parts);
  if (itr != mesh_impl().m_parts_to_partition.end()) {
    return itr->second;
  }
  return partition_id::invalid();
}

template <typename MeshImpl>
inline
samba::spatial_dimension query_impl<MeshImpl>::spatial_dimension() const
{ return mesh_impl().m_connectivity_map.spatial_dimension(); }

template <typename MeshImpl>
inline
const samba::connectivity_map& query_impl<MeshImpl>::connectivity_map() const
{ return mesh_impl().m_connectivity_map; }

template <typename MeshImpl>
inline
set_expression query_impl<MeshImpl>::get_set_expression_by_name(const std::string& name) const
{
  typename MeshImpl::expr_db_type::const_iterator fitr = mesh_impl().m_expr_db.find(name);
  BOOST_ASSERT_MSG(fitr != mesh_impl().m_expr_db.end(),
                   (debug_message() << "No set-expression registered for name: " << name));
  return fitr->second;
}

template <typename MeshImpl>
inline
std::ostream& query_impl<MeshImpl>::streamit(std::ostream &os, size_t verbosity)
{
  os << "m_process = " << mesh_impl().m_process << std::endl;
  os << "spatial_dimension = " << mesh_impl().spatial_dimension() << std::endl;
  os << "m_is_in_modification_cycle = " <<  mesh_impl().m_is_in_modification_cycle << std::endl;

  os << "m_set_specs = ";
  size_t num_set_specs = mesh_impl().m_set_specs.size();
  for (size_t i = 0; i < num_set_specs; ++i)
  {
    os << "{entity_block_spec " << mesh_impl().m_set_specs[i].name << " " << mesh_impl().m_set_specs[i].rank << " "
       << (mesh_impl().m_set_specs[i].inducable ? "INDUCIBLE" : "NOT INDUCIBLE") << "}\n";
  }

  // Step through partitions and output their contents
  size_t num_partitions = mesh_impl().m_partitions.size();
  for (size_t i = 0; i < num_partitions; ++i)
  {
    typename MeshImpl::partition_handle p_handle = mesh_impl().m_partitions[i];
    os << "{" << p_handle->partition() << " " << p_handle->topology()
       << " rank=" << p_handle->rank()
       << " size=" << p_handle->size()
       << " parts=(";
    entity_part_vector const& parts = p_handle->parts();
    for (size_t i = 0; i < parts.size(); ++i) {
      os << parts[i] << " ";
    }
    os << ") {\n";

    size_t partition_size = p_handle->size();
    for (partition_offset offset = {0}; offset < partition_size; ++offset)
    {
      os << "  " <<  p_handle->key(offset) << "\n";

      for (entity_rank rank = entity_rank::node(); rank <= entity_rank::element(); ++rank) {
        if (p_handle->connectivity_kind(rank) != connectivity_kind::invalid()) {
          os << "    connectivity for rank: " << rank << "\n";
          entity_key_iterator begin = p_handle->template begin_connectivity<entity_key>(rank, offset);
          entity_key_iterator end   = p_handle->template end_connectivity<entity_key>(rank, offset);
          for (entity_key_iterator itr = begin; itr != end; ++itr) {
            os << "      " << *itr << "\n";
          }
        }
      }
    }
    os << "}}\n";
  }

  os << "m_key_to_descriptor_map = {\n";
  for (size_t top = 0; top < mesh_impl().m_key_to_descriptor_map.size(); ++top) {
    os << "  " << entity_topology::create(top) << "\n";
    for (size_t loc = 0; loc < mesh_impl().m_key_to_descriptor_map[top].size(); ++loc) {
      os << "    " << loc << ": " << mesh_impl().m_key_to_descriptor_map[top][loc] << "\n";
    }
  }
  os << "}\n";

  os << "m_descriptor_to_key_map = {\n";
  for (size_t par = 0; par < mesh_impl().m_descriptor_to_key_map.size(); ++par) {
    os << "  " << partition_id::create(par) << "\n";
    for (size_t off = 0; off < mesh_impl().m_descriptor_to_key_map[par].size(); ++off) {
      os << "    " << off << ": " << mesh_impl().m_descriptor_to_key_map[par][off] << "\n";
    }
  }
  os << "}\n";

  os << std::endl;

  return os;
}

template <typename MeshImpl>
inline
bool query_impl<MeshImpl>::modifiable() const
{ return mesh_impl().m_is_in_modification_cycle; }

template <typename MeshImpl>
inline convertable_to_partition_proxy query_impl<MeshImpl>::operator[](partition_id p) const
{
  mesh_impl().validate_helper(p);
  return mesh_impl().m_partitions[p()].get();
}

template <typename MeshImpl>
inline convertable_to_entity_proxy query_impl<MeshImpl>::operator[](entity_key k) const
{
  mesh_impl().validate_helper(k);
  return (*this)[mesh_impl().convert(k)];
}

template <typename MeshImpl>
inline convertable_to_entity_proxy query_impl<MeshImpl>::operator[](partition_index d) const
{
  mesh_impl().validate_helper(d);
  return convertable_to_entity_proxy(mesh_impl().m_partitions[d.partition()()].get(),d.offset());
}

} } // namespace samba::detail
