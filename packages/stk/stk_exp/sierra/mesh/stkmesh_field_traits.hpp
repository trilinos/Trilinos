#ifndef SIERRA_MESH_STKMESH_FIELD_TRAITS_HPP
#define SIERRA_MESH_STKMESH_FIELD_TRAITS_HPP

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/DataTraits.hpp>
#include <stk_mesh/base/FieldTraits.hpp>

#include <sierra/mesh/field_traits.hpp>

namespace sierra {
namespace mesh {

template <typename Scalar, typename Tag1, typename Tag2, typename Tag3, typename Tag4, typename Tag5, typename Tag6, typename Tag7>
struct field_traits<stk::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
{
  typedef typename stk::mesh::FieldTraits<stk::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >::data_type data_type;
};

} // namespace mesh
} // namespace sierra

#endif // SIERRA_MESH_STKMESH_FIELD_TRAITS_HPP
