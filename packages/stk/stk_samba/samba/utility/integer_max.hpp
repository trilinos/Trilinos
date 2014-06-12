#ifndef SAMBA_SAMBA_UTILITY_INTEGER_MAX_HPP
#define SAMBA_SAMBA_UTILITY_INTEGER_MAX_HPP

#include <boost/integer_traits.hpp>


namespace samba {

template <typename Integer>
struct integer_max
{
  typedef integer_max<Integer> type;
  typedef Integer value_type;
  static const value_type value = boost::integer_traits<Integer>::const_max;
  value_type operator()() const { return value; }
};



} // namespace samba

#endif //SAMBA_SAMBA_UTILITY_INTEGER_MAX_HPP
