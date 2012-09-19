#ifndef SAMBA_SAMBA_UTILITY_MDARRAY_DIMENSION_HELPER_HPP
#define SAMBA_SAMBA_UTILITY_MDARRAY_DIMENSION_HELPER_HPP

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/assert.hpp>

namespace samba { namespace detail {

template <typename Array>
inline size_t dimension_helper(size_t dim)
{
  BOOST_ASSERT( (dim < boost::rank<Array>::value) );

  switch(dim)
  {
    case 0: return boost::extent<Array,0>::value;
    case 1: return boost::extent<Array,1>::value;
    case 2: return boost::extent<Array,2>::value;
    case 3: return boost::extent<Array,3>::value;
    case 4: return boost::extent<Array,4>::value;
    case 5: return boost::extent<Array,5>::value;
    case 6: return boost::extent<Array,6>::value;
    case 7: return boost::extent<Array,7>::value;
    case 8: return boost::extent<Array,8>::value;
    case 9: return boost::extent<Array,9>::value;
    default: return static_cast<size_t>(-1);
  }
  return static_cast<size_t>(-1);
}

}} //namespace samba::detail


#endif //SAMBA_SAMBA_UTILITY_MDARRAY_DIMENSION_HELPER_HPP


