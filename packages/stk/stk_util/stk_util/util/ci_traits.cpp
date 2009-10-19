#include <stk_util/util/ci_traits.hpp>

int
ignorecase_traits::compare(
  const char *          s1,
  const char *          s2,
  std::size_t           n)
{
  for (std::size_t i = 0; i < n; ++i)
    if (!eq(s1[i], s2[i]))
      return lt(s1[i], s2[i]) ? -1 : 1;

  return 0;
}
  
const char *
ignorecase_traits::find(
  const char *          s,
  std::size_t           n,
  const char &          c)
{
  for (std::size_t i = 0; i < n; ++i)
    if (eq(s[i], c))
      return &(s[i]);

  return 0;
}
