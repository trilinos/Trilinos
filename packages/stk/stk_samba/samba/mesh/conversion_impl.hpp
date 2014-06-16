#ifndef SAMBA_SAMBA_MESH_CONVERSION_IMPL_HPP
#define SAMBA_SAMBA_MESH_CONVERSION_IMPL_HPP

namespace samba { namespace detail {

/**
 * An experiment to see what it would look like if we used the "curiously
 * recurring template pattern" to group/implement categories of
 * methods for mesh_impl. Based on current thinking, these categories
 * will be: queries, conversion, and modification. Setting-up the one
 * for modification is a TODO item.
 */
template <typename MeshImpl>
class conversion_impl
{
 public:
  partition_index convert(entity_key key) const;

  entity_key convert(partition_index descriptor) const;

  template <typename ToIndex, typename FromIndex>
  ToIndex convert(FromIndex from) const;

 protected:
  ~conversion_impl() {}

  template <typename ToIndex, typename FromIndex>
  ToIndex convert_helper(FromIndex from) const;

 private:
  const MeshImpl& mesh_impl() const
  { return *static_cast<const MeshImpl*>(this); }
};

} } // namespace samba::detail

#endif // SAMBA_SAMBA_MESH_CONVERSION_IMPL_HPP
