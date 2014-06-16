#ifndef SAMBA_SAMBA_UTILITY_IS_PRIMITIVE_HPP
#define SAMBA_SAMBA_UTILITY_IS_PRIMITIVE_HPP

#include <boost/mpl/bool.hpp>

namespace samba {

template <typename Descriptor>
struct is_primitive
  : public boost::mpl::false_
{};

} //namespace samba

#endif //SAMBA_SAMBA_UTILITY_IS_PRIMITIVE_HPP

