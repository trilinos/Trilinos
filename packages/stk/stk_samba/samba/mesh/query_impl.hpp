#ifndef SAMBA_SAMBA_MESH_QUERY_IMPL_HPP
#define SAMBA_SAMBA_MESH_QUERY_IMPL_HPP

namespace samba { namespace detail {

/**
 * An experiment to see what it would look like if we used the "curiously
 * recurring template pattern" to group/implement categories of
 * methods for mesh_impl. Based on current thinking, these categories
 * will be: queries, conversion, and modification. Setting-up the one
 * for modification is a TODO item.
 */
template <typename MeshImpl>
class query_impl
{
 public:
  entity_block_key find_entity_block(std::string const& name) const;

  entity_rank get_rank(entity_block_key key) const;

  const std::string& get_name(entity_block_key key) const;

  size_t num_entity_blocks() const;

  size_t num_partitions() const;

  process_id process() const;

  entity_local_id max_local_id(entity_topology topology) const;

  partition_id get_partition( entity_part_vector const& parts ) const;

  samba::spatial_dimension spatial_dimension() const;

  const samba::connectivity_map& connectivity_map() const;

  set_expression get_set_expression_by_name(const std::string& name) const;

  std::ostream &streamit(std::ostream &os, size_t verbosity = 10);

  bool modifiable() const;

  // TODO - Remove convertable_to_X classes

  convertable_to_partition_proxy operator[](partition_id so) const;

  convertable_to_entity_proxy operator[](entity_key key) const;

  convertable_to_entity_proxy operator[](partition_index descriptor) const;

 protected:
  ~query_impl() {}

 private:
  const MeshImpl& mesh_impl() const
  { return *static_cast<const MeshImpl*>(this); }
};

} } // namespace samba::detail

#endif // SAMBA_SAMBA_MESH_QUERY_IMPL_HPP
