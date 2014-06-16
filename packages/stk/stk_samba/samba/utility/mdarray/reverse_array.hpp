#ifndef SAMBA_SAMBA_UTILITY_MDARRAY_REVERSE_ARRAY_HPP
#define SAMBA_SAMBA_UTILITY_MDARRAY_REVERSE_ARRAY_HPP

#include <samba/utility/mdarray/add_extent.hpp>

#include <boost/type_traits.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/utility/enable_if.hpp>


namespace samba { namespace detail {

template <typename Array, typename ReversedArray, typename Enable = void>
struct reverse_array_helper;

template <typename Array, typename ReversedArray>
struct reverse_array_helper< Array
                            ,ReversedArray
                            ,typename boost::enable_if< boost::is_array<Array> >::type
                           >
{
  BOOST_MPL_ASSERT_MSG( boost::extent<Array>::value > 0u
                       ,CANNOT_REVERSE_ARRAY_WITH_EXTENT_OF_ZERO
                       ,(Array)
                      );

  static const size_t extent = boost::extent<Array>::value;

  typedef typename boost::remove_extent<Array>::type array;
  typedef typename add_extent<ReversedArray,extent>::type reversed_array;

  typedef typename reverse_array_helper<array,reversed_array>::type type;

};

template <typename Array, typename ReversedArray>
struct reverse_array_helper< Array
                            ,ReversedArray
                            ,typename boost::disable_if< boost::is_array<Array> >::type
                           >
{
  typedef ReversedArray type;
};


template <typename Array>
struct reverse_array
{
  typedef Array array;
  typedef typename boost::remove_all_extents<array>::type value_type;
  typedef typename detail::reverse_array_helper<array,value_type>::type type;
};

}} //namespace samba::detail

#endif //SAMBA_SAMBA_UTILITY_MDARRAY_REVERSE_ARRAY_HPP

