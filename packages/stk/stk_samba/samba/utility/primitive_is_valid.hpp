#ifndef SAMBA_SAMBA_UTILITY_PRIMITIVE_IS_VALID_HPP
#define SAMBA_SAMBA_UTILITY_PRIMITIVE_IS_VALID_HPP

#include <samba/utility/is_primitive.hpp>
#include <samba/utility/tag_of.hpp>
#include <samba/utility/value_of.hpp>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace samba {

template <typename T>
inline typename boost::enable_if< is_primitive<T>,bool>::type
is_valid(T const& t)
{ return t() != T::invalid()(); }

} // namespace samba

#endif //SAMBA_SAMBA_UTILITY_PRIMITIVE_IS_VALID_HPP


