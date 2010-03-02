/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_UTIL_CI_STRING_H
#define STK_UTIL_UTIL_CI_STRING_H

#include <stk_util/util/ci_traits.hpp>
#include <iosfwd>

// #include <boost/unordered_set.hpp> // boost 1.36.0 unordered_set not in the tr1

// namespace std {
// namespace tr1 {
//   using ::boost::unordered_set;
// }
// }

typedef std::basic_string<char,ignorecase_traits> ci_string;

std::ostream &operator<<(std::ostream &os, const ci_string &s);

std::istream &operator>>(std::istream &is, ci_string &s);

std::string operator+(const std::string &s1, const ci_string &s2);

ci_string operator+(const ci_string &s1, const std::string &s2);

ci_string operator+(const char *cs1, const ci_string &cs2);

// namespace boost {

// template <>
// struct hash<ci_string>
// {
//   std::size_t operator()(const ci_string &s) const;
// };
  
// } // namespace boost

#endif // STK_UTIL_UTIL_CI_STRING_H
