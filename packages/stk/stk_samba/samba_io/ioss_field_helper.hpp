#ifndef SAMBA_SAMBA_IO_IOSS_FIELD_HELPER_HPP
#define SAMBA_SAMBA_IO_IOSS_FIELD_HELPER_HPP

#include <Ioss_SubSystem.h>
#include <samba/mesh.hpp>

namespace samba {
namespace io {
namespace detail {

/**
 * Templated struct to enable use of specialization to define functions that
 * compute arguments for Ioss::Field constructor.
 */
template<typename DataT, typename DimFunctorT>
struct IossFieldHelper
{
  Ioss::Field::BasicType basicType() const { return Ioss::Field::INVALID; }
  std::string storage(mesh mesh_arg) const { return "UNKNOWN";}
};

template<>
struct IossFieldHelper<double, scalar_functor>
{
  Ioss::Field::BasicType basicType() const { return Ioss::Field::REAL; }
  std::string storage(mesh mesh_arg) const { return "scalar";}
};


#if 0
//
// Need to get the right strings for Ioss_Field.
//
template<>
struct IossFieldHelper<int64_t, scalar_functor>
{
  Ioss::Field::BasicType basicType() const {
//     std::cout << "IossFieldHelper<int64_t, scalar_functor> " << std::endl;
    return Ioss::Field::INT64;
  }
  std::string storage(mesh mesh_arg) const { return "INT64";}
};
#endif


template<>
struct IossFieldHelper<double, spatial_dimension_functor>
{
  Ioss::Field::BasicType basicType() const {  return Ioss::Field::REAL; }

  std::string storage(mesh mesh_arg) const
  {
    switch (mesh_arg.spatial_dimension()())
    {
    case 2:
      return "vector_2d";
    case 3:
      return "vector_3d";
    default:
      return "INVALID";
    }
    return "UNKNOWN";
  }
};

} // namespace detail
} // namespace io
} // namespace samba
#endif
