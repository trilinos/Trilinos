/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

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
