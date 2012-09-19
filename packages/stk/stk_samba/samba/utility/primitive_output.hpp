#ifndef SAMBA_SAMBA_UTILITY_PRIMITIVE_OUTPUT_HPP
#define SAMBA_SAMBA_UTILITY_PRIMITIVE_OUTPUT_HPP

#include <samba/utility/is_primitive.hpp>
#include <samba/utility/tag_of.hpp>
#include <samba/utility/value_of.hpp>
#include <samba/utility/primitive_is_valid.hpp>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <ostream>

namespace samba {

template <typename T>
inline typename boost::enable_if<is_primitive<T>, std::ostream&>::type
operator<<(std::ostream& out, T const& t)
{
  typedef typename tag_of<T>::type tag;
  typedef typename boost::integral_promotion<
    typename value_of<T>::type
  >::type value_type;
  out << "{" << tag() << ":";
  if (is_valid(t))
    out << static_cast<value_type>(t()) << "}";
  else
    out << "invalid}";
  return out;
}


} // namespace samba

#endif //SAMBA_SAMBA_UTILITY_PRIMITIVE_OUTPUT_HPP

