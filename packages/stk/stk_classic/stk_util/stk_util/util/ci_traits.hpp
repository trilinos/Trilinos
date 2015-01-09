/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_UTIL_CI_TRAITS_H
#define STK_UTIL_UTIL_CI_TRAITS_H

#include <cctype>
#include <string>

/**
 * @brief Class <b>ignorecase_traits</b> is a character traits class that ignores case during compares.
 *
 * Replace functions of the standard char_traits<char> so that strings behave in a case-insensitive
 * way.
 *
 */
struct ignorecase_traits : public std::char_traits<char>
{
  /**
   * @brief Member function <b>eq</b> return true is c1 and c2 are equal.
   *
   * @param c1			a <b>char</b> const reference to character to compare. 
   *
   * @param c2			a <b>char</b> const reference to character to compare. 
   *
   * @return			a <b>bool</b> of true if c1 and c2 are equal
   */
  static bool eq(const char &c1, const char &c2) {
    return std::toupper(c1) == std::toupper(c2);
  }
  
  /**
   * @brief Member function <b>lt</b> return true is c1 less than c2.
   *
   * @param c1			a <b>char</b> const ...
   *
   * @param c2			a <b>char</b> const ...
   *
   * @return			a <b>bool</b> ...
   */
  static bool lt(const char &c1, const char &c2) {
    return std::toupper(c1) < std::toupper(c2);
  }
  
  /**
   * @brief Member function <b>compare</b> compares up to n characters of s1 and s2 and returns -1 if
   * s1 is less then s2, 0 if they are equal, and 1 if s1 is greater than s2.
   *
   * @param s1			a <b>char</b> const pointer to string to compare.
   *
   * @param s2			a <b>char</b> const pointer to string to compare.
   *
   * @param n			a <b>std::size_t</b> maxiumum number of character to compare. 
   *
   * @return                    an <b>int</b> value of -1 if s1 is less then s2, 0 if they are
   *                            equal, and 1 if s1 is greater than s2.
   */
  static int compare(const char *s1, const char *s2, std::size_t n);
  
  /**
   * @brief Member function <b>find</b> returns char pointer to first occurrence of character
   * <b>c</b> in first <b>n</b> characters of string <b>s</b> or 0 if not found.
   *
   * @param s			a <b>char</b> const pointer to string to search in.
   *
   * @param n			a <b>std::size_t</b> value of the maximum number of characters to
   *                            compare. 
   *
   * @param c			a <b>char</b> const reference to the character to search.
   *
   * @return			a <b>char</b> pointer to first occurrence or 0 is not found.
   */
  static const char *find(const char *s, std::size_t n, const char &c);
};

#endif // STK_UTIL_UTIL_CI_TRAITS_H
