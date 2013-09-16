namespace samba {

template <typename MeshHandle>
inline
samba::connectivity_map const& query_interface<MeshHandle>::connectivity_map() const
{ return mesh_impl().connectivity_map(); }

template <typename MeshHandle>
inline
partition_proxy query_interface<MeshHandle>::operator[](partition_id so) const
{ return mesh_impl()[so]; }

template <typename MeshHandle>
inline
entity_proxy query_interface<MeshHandle>::operator[](entity_key key) const
{ return mesh_impl()[key]; }

template <typename MeshHandle>
inline
entity_proxy query_interface<MeshHandle>::operator[](partition_index descriptor) const
{ return mesh_impl()[descriptor]; }

template <typename MeshHandle>
inline
entity_block_key query_interface<MeshHandle>::find_entity_block(std::string const& name) const
{ return mesh_impl().find_entity_block(name); }

template <typename MeshHandle>
inline
void query_interface<MeshHandle>::get_entity_blocks(entity_rank rank, std::vector<entity_block_key> &sets) const
{
  sets.clear();
  for (size_t block = 0, end = num_entity_blocks(); block < end; ++block)
  {
    entity_block_key block_key = entity_block_key::create(block);
    if (get_rank(block_key) == rank)
    {
      sets.push_back(block_key);
    }
  }
}

template <typename MeshHandle>
template <typename Index>
inline
void query_interface<MeshHandle>::get_entities(set_expression const & selector_expr, std::vector<Index> &entities_out) const
{
  entities_out.clear();
  for (partition_id pd = {0}; pd() < num_partitions(); ++pd)
  {
    partition_proxy ppx = (*this)[pd];
    if (contains(selector_expr, ppx.parts() ))
    {
      for (partition_proxy::const_iterator iter = ppx.begin(), end = ppx.end() ;
           iter != end;
           ++iter)
      {
        entities_out.push_back(iter->get<Index>());
      }
    }
  }
}

template <typename MeshHandle>
inline
void query_interface<MeshHandle>::get_partitions(set_expression const & selector_expr, std::vector<partition_id> &partitions_out) const
{
  partitions_out.clear();
  for (partition_id pd = {0}; pd() < num_partitions(); ++pd)
  {
    partition_proxy ppx = (*this)[pd];
    if (contains(selector_expr, ppx.parts() ))
    {
      partitions_out.push_back(pd);
    }
  }
}

template <typename MeshHandle>
inline
size_t query_interface<MeshHandle>::num_entity_blocks() const
{ return mesh_impl().num_entity_blocks(); }

template <typename MeshHandle>
inline
entity_rank query_interface<MeshHandle>::get_rank(entity_block_key key) const
{ return mesh_impl().get_rank(key); }

template <typename MeshHandle>
inline
const std::string& query_interface<MeshHandle>::get_name(entity_block_key key) const
{ return mesh_impl().get_name(key); }

template <typename MeshHandle>
inline
size_t query_interface<MeshHandle>::num_partitions() const
{ return mesh_impl().num_partitions(); }

template <typename MeshHandle>
inline
process_id query_interface<MeshHandle>::process() const
{ return mesh_impl().process(); }

template <typename MeshHandle>
inline
samba::spatial_dimension query_interface<MeshHandle>::spatial_dimension() const
{ return mesh_impl().spatial_dimension(); }

template <typename MeshHandle>
inline
bool query_interface<MeshHandle>::modifiable() const
{ return mesh_impl().modifiable(); }

template <typename MeshHandle>
inline
partition_id query_interface<MeshHandle>::get_partition( entity_part_vector const& parts ) const
{ return mesh_impl().get_partition(parts); }

template <typename MeshHandle>
inline
partition_index query_interface<MeshHandle>::convert(entity_key key) const
{ return mesh_impl().convert(key); }

template <typename MeshHandle>
inline
entity_key query_interface<MeshHandle>::convert(partition_index descriptor) const
{ return mesh_impl().convert(descriptor); }

template <typename MeshHandle>
template <typename ToIndex, typename FromIndex>
inline
ToIndex query_interface<MeshHandle>::convert(FromIndex from) const
{ return mesh_impl().template convert<ToIndex>(from); }

template <typename MeshHandle>
inline
size_t query_interface<MeshHandle>::num_entities(set_expression const & expr) const
{
  size_t count = 0;
  for(partition_id i={0}; i<num_partitions(); ++i)
  {
    const partition_proxy p = (*this)[i];
    if (contains(expr, p))
    {
      count += p.size();
    }
  }
  return count;
}

template <typename MeshHandle>
inline
size_t query_interface<MeshHandle>::num_elements() const
{ return num_entities(entity_rank::element()); }

template <typename MeshHandle>
inline
size_t query_interface<MeshHandle>::num_faces() const
{ return num_entities(entity_rank::face()); }

template <typename MeshHandle>
inline
size_t query_interface<MeshHandle>::num_edges() const
{ return num_entities(entity_rank::edge()); }

template <typename MeshHandle>
inline
size_t query_interface<MeshHandle>::num_nodes() const
{ return num_entities(entity_rank::node()); }

template <typename MeshHandle>
inline
set_expression query_interface<MeshHandle>::get_set_expression_by_name(const std::string& name) const
{ return mesh_impl().get_set_expression_by_name(name); }

template <typename MeshHandle>
inline
detail::mesh_impl * query_interface<MeshHandle>::detail_get_raw_mesh_impl()
{ return static_cast<MeshHandle*>(this)->m_mesh_impl; }

} // namespace samba
