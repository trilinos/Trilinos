#include <iostream>

#include <stk_util/util/ci_string.hpp>

std::ostream &
operator<<(
  std::ostream &        os,
  const ci_string &      s)
{
  return os << std::string(s.data(),s.length());
}
 

std::istream &
operator>>(
  std::istream &        is,
  ci_string &            s)
{
  std::string t;
  is >> t;
  s = ci_string(t.begin(), t.end());
  
  return is;
}

std::string
operator+(
  const std::string &   s1,
  const ci_string &      s2) 
{
  std::string s(s1);
  
  return s.append(s2.begin(), s2.end());
}


ci_string
operator+(
  const char *          s1,
  const ci_string &      s2) 
{
  ci_string s(s1);
  
  return s.append(s2.begin(), s2.end());
}


ci_string
operator+(
  const ci_string &      s1,
  const std::string &   s2) 
{
  ci_string s(s1);
  
  return s.append(s2.begin(), s2.end());
}


// namespace boost {

// std::size_t
// hash<ci_string>::operator()(
//   const ci_string &      s) const
// {
//   std::size_t seed = 0;
    
//   for(ci_string::const_iterator first = s.begin(); first != s.end(); ++first) {
//     std::tr1::hash<char> hasher;
//     seed ^= hasher(std::tolower(*first)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
//   }
    
//   return seed;
// }

// } // namespace boost
