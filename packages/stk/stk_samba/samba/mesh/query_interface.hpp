#ifndef SAMBA_SAMBA_MESH_QUERY_INTERFACE_HPP
#define SAMBA_SAMBA_MESH_QUERY_INTERFACE_HPP

namespace samba {

/**
 * The interface class for this package. Clients will interact with this class
 * to perform mesh operations.
 */
template <typename MeshHandle>
class query_interface
{
 public:

  samba::connectivity_map const& connectivity_map() const;

  partition_proxy operator[](partition_id so) const;

  entity_proxy operator[](entity_key key) const;

  entity_proxy operator[](partition_index descriptor) const;

  entity_block_key find_entity_block(std::string const& name) const;

  void get_entity_blocks(entity_rank rank, std::vector<entity_block_key> &sets) const;

  template <typename Index>
  void get_entities(set_expression const & selector_expr, std::vector<Index> &entities_out) const;

  void get_partitions(set_expression const & selector_expr, std::vector<partition_id> &partitions_out) const;

  size_t num_entity_blocks() const;

  entity_rank get_rank(entity_block_key key) const;

  const std::string& get_name(entity_block_key key) const;

  size_t num_partitions() const;

  process_id process() const;

  samba::spatial_dimension spatial_dimension() const;

  bool modifiable() const;

  partition_id get_partition( entity_part_vector const& parts ) const;

  partition_index convert(entity_key key) const;

  entity_key convert(partition_index descriptor) const;

  template <typename ToIndex, typename FromIndex>
  ToIndex convert(FromIndex from) const;

  size_t num_entities(set_expression const & expr) const;

  size_t num_elements() const;

  size_t num_faces() const;

  size_t num_edges() const;

  size_t num_nodes() const;

  set_expression get_set_expression_by_name(const std::string& name) const;

  //need to break shared_ptr cycle when the mesh holds a field
  detail::mesh_impl * detail_get_raw_mesh_impl();

 protected:
  ~query_interface() {}

 private:
  const detail::mesh_impl& mesh_impl() const
  { return *(static_cast<const MeshHandle*>(this)->m_mesh_impl); }
};

} //namespace samba

#endif // SAMBA_SAMBA_MESH_QUERY_INTERFACE_HPP
