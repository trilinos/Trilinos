#ifndef AMESOS2_UTIL_IS_SAME_HPP
#define AMESOS2_UTIL_IS_SAME_HPP

namespace Amesos {

namespace Util {

/*==================== Simple Template Metaprogramming ==================== */

template <class T, T val>
struct integral_constant
{
  typedef integral_constant<T, val>  type;
  typedef T                          value_type;
  static const T value;
};

/* Some compilers support initializing static const members alongside the
 * definition, but others do not, so we go we the safe method of external
 * initialization.
 */
template <class T, T val>
const T integral_constant<T,val>::value = val;

typedef integral_constant<bool, true>  true_type;
typedef integral_constant<bool, false> false_type;

template <typename, typename>
struct is_same : public false_type 
{};

template <typename T>
struct is_same<T,T> : public true_type
{};

} // end namespace Util

} // end namespace Amesos

#endif  // AMESOS2_UTIL_IS_SAME_HPP
