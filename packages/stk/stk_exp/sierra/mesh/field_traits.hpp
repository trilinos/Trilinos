#ifndef SIERRA_MESH_GENERIC_FIELD_TRAITS_HPP
#define SIERRA_MESH_GENERIC_FIELD_TRAITS_HPP

namespace sierra {
namespace mesh {

// General version, all associated types are expected to be defined
// in the Field class.
template <class Field>
struct field_traits
{
  typedef typename Field::value_type          data_type;
};

} // namespace mesh
} // namespace sierra

#endif // SIERRA_MESH_GENERIC_FIELD_TRAITS_HPP
