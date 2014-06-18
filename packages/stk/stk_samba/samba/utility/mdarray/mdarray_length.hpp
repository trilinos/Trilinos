#ifndef SAMBA_SAMBA_UTILITY_MDARRAY_MDARRAY_LENGTH_HPP
#define SAMBA_SAMBA_UTILITY_MDARRAY_MDARRAY_LENGTH_HPP

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace samba { namespace detail {

template <typename Array, typename Enable = void>
struct mdarray_length;

template< typename Array >
struct mdarray_length<Array, typename boost::enable_if_c< (boost::rank<Array>::value > 1u) >::type>
{
  typedef Array array;
  typedef typename boost::remove_extent<array>::type sub_array;

  typedef size_t value_type;
  static const size_t value = boost::extent<array>::value * mdarray_length<sub_array>::value;
};

template< typename Array >
struct mdarray_length<Array, typename boost::enable_if_c< (boost::rank<Array>::value == 1u) >::type>
{
  typedef Array array;
  typedef typename boost::remove_extent<array>::type sub_array;

  typedef size_t value_type;
  static const size_t value = boost::extent<array,0u>::value;
};

template< typename Array >
struct mdarray_length<Array, typename boost::enable_if_c< (boost::rank<Array>::value == 0u) >::type>
{
  typedef Array array;
  typedef void sub_array;

  typedef size_t value_type;
  static const size_t value = 1u;
};

}} //namespace samba::detail


#endif //SAMBA_SAMBA_UTILITY_MDARRAY_MDARRAY_LENGTH_HPP

