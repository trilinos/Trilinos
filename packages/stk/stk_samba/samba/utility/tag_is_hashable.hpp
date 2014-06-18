#ifndef SAMBA_SAMBA_UTILITY_TAG_IS_HASHABLE_HPP
#define SAMBA_SAMBA_UTILITY_TAG_IS_HASHABLE_HPP

#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>

#include <boost/utility/enable_if.hpp>

namespace samba {

template <typename T>
struct tag_is_hashable
{
  typedef char   yes;
  typedef size_t no;

  template <typename U> struct check;

  template <typename C> static yes is_comparable(check< typename C::value_type>*);
  template <typename C> static no  is_comparable(...);

  typedef bool value_type;
  static bool const value = (sizeof(is_comparable<T>(0)) == sizeof(yes));

  typedef boost::mpl::bool_<value> type;
};

template <typename T, typename Enable = void>
struct tag_hash
  : boost::mpl::integral_c<size_t,0>
{};

template <typename T>
struct tag_hash< T, typename boost::enable_if< tag_is_hashable<T> >::type>
  : boost::mpl::integral_c<size_t,T::value>
{};

} //namespace samba

#endif //SAMBA_SAMBA_UTILITY_TAG_IS_HASHABLE_HPP
