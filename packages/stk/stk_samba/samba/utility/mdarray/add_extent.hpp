#ifndef SAMBA_SAMBA_UTILITY_MDARRAY_ADD_EXTENT_HPP
#define SAMBA_SAMBA_UTILITY_MDARRAY_ADD_EXTENT_HPP

namespace samba { namespace detail {

template <typename T, size_t Extent>
struct add_extent
{ typedef T type[Extent]; };


}} //namespace samba::detail

#endif //SAMBA_SAMBA_UTILITY_MDARRAY_ADD_EXTENT_HPP

