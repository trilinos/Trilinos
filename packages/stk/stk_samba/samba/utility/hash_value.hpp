#ifndef SAMBA_SAMBA_UTILITY_HASH_VALUE_HPP
#define SAMBA_SAMBA_UTILITY_HASH_VALUE_HPP

#include <samba/utility/is_primitive.hpp>
#include <samba/utility/tag_of.hpp>
#include <samba/utility/tag_is_hashable.hpp>
#include <samba/utility/value_of.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>

#include <boost/functional/hash.hpp>

namespace samba {

// overload the hash_value function for samba primitives
template <typename T>
inline typename
boost::enable_if<is_primitive<T> ,size_t>::type
hash_value(T t)
{
  typedef typename T::tag tag;
  size_t seed = tag_hash<tag>::value;
  boost::hash_combine(seed,t());
  return seed;
}

} // namespace samba

#endif  //SAMBA_SAMBA_UTILITY_HASH_VALUE_HPP

