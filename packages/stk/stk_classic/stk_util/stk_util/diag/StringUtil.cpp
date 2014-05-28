/**   ------------------------------------------------------------
 *    Copyright 2002 - 2008 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cstring>
#include <vector>
#include <typeinfo>

#include <float.h>
#include <time.h>

#include <stk_util/diag/StringUtil.hpp>

// Set for unit testing

#define SIERRA_STRING_UNIT_TEST 0

//----------------------------------------------------------------------

namespace sierra {

std::istream &
getline(
  std::istream &		is,
  sierra::String &		s,
  char				eol)
{
  std::string std_string;

  getline(is, std_string, eol);
  s = std_string.c_str();
  return is;
}


int
case_strcmp(
  const char *		c1,
  const char *		c2)
{
  for ( ; ; c1++, c2++) {
    if ( std::tolower(*c1) != std::tolower(*c2) )
      return ( std::tolower(*c1) - std::tolower(*c2) ) ;
    if (*c1 == '\0')
      return 0 ;
  }
}


int
case_strncmp(
  const char *		c1,
  const char *		c2,
  size_t		n)
{

  if (n == 0)
    return 0;

  do {
    if (*c1 != *c2++)
      return std::tolower(*c1) - std::tolower(*(c2 - 1));
    if (*c1++ == 0)
      break;
  } while (--n != 0);
  return 0;
}


const char *
case_strstr(
  const char *		s,
  const char *		find)
{
  if (!*find)
    return s;

  const char *cp = s;
  while (*cp) {
    const char *t1 = cp;
    const char *t2 = find;

    for ( ; std::tolower(*t1) && std::tolower(*t2) && !(std::tolower(*t1) - std::tolower(*t2)); ++t1, ++t2)
      ;

    if (!*t2)
      return cp;

    cp++;
  }

  return NULL;
}


std::string
title(
  const std::string &	s)
{
  std::string t(s);

  bool all_upper = true;
  bool all_lower = true;

  bool next_upper = true;
  for (std::string::iterator c = t.begin(); c != t.end(); ++c) {
    all_upper &= (*c == std::toupper(*c));
    all_lower &= (*c == std::tolower(*c));
    if (next_upper)
      *c = std::toupper(*c);
    else
      *c = std::tolower(*c);
    next_upper = !isalpha(*c);
  }
  
  if (all_upper || all_lower) 
    return t;
  else
    return s;
}


template <typename T>
std::string
to_string(
  const T &		t)
{
  std::ostringstream os;
  os << t;
  return os.str();
}

template std::string to_string<double>(const double &);
template std::string to_string<float>(const float &);
template std::string to_string<int>(const int &);
template std::string to_string<unsigned>(const unsigned &);
template std::string to_string<long>(const long &);
template std::string to_string<unsigned long>(const unsigned long &);

std::string
to_string(
  const double &	r,
  int			precision)
{
  std::ostringstream os;
  os << std::setprecision(precision) << r;
  return std::string(os.str());
}


std::string
to_string(
  const float &		r,
  int			precision)
{
  std::ostringstream os;
  os << std::setprecision(precision) << r;
  return std::string(os.str());
}


std::string
format_time(
  double	t,
  const char *	format)
{
  time_t	time = (time_t) t;
  char s[128];

  ::strftime(s, sizeof(s), format, ::localtime(&time));

  return std::string(s);
}


std::ostream &
object_phrase::print(
  std::ostream &	os) const
{
  if (m_n == 0)
    os << m_plural << " no " << m_noun << "s";
  else if (m_n == 1)
    os << m_singular << " 1 " << m_noun;
  else
    os << m_plural << " " << m_n << " " << m_noun << "s";

  return os;
}


object_phrase::operator std::string() const
{
  std::ostringstream strout;
  strout << *this;

  return strout.str();
}


namespace {

std::string::const_iterator
find_next_char(
  std::string::const_iterator	p,
  std::string::const_iterator	end,
  char				c)
{
  while (p != end && *p != c)
    p++;
  return p;
}

std::string::const_iterator
find_next_not_char(
  std::string::const_iterator	p,
  std::string::const_iterator	end,
  char				c)
{
  while (p != end && *p == c)
    p++;
  return p;
}

inline std::string::const_iterator find_next_space(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_char(p, end, ' ');
}

inline std::string::const_iterator find_next_endl(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_char(p, end, '\n');
}

inline std::string::const_iterator find_next_nonspace(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_not_char(p, end, ' ');
}

} // namespace <null>

std::string
word_wrap(
  const std::string &	s,
  unsigned int		line_length,
  const std::string &	prefix,
  const std::string &	prefix_first_line)
{
  std::string t;
  const std::string *u = &prefix_first_line;

  std::string::const_iterator p0, p1, p2, p3;
  p0 = p1 = p2 = s.begin();

  while (p2 != s.end() ) {

    // skip preceeding whitespace
    p1 = find_next_nonspace(p0, s.end());
    p3 = find_next_endl(p0, s.end());
    p2 = p1 = find_next_space(p1, s.end());
    do {
      p1 = find_next_nonspace(p1, s.end());
      p1 = find_next_space(p1, s.end());
      if (p3 < p1) {
	p2 = p3;
	break;
      }
      if ((unsigned int) (p1 - p0) > (line_length - u->size()))
	break;
      p2 = p1;
    } while (p2 != s.end());

    t.append(*u).append(p0, p2).append("\n");

    if (p2 == p3)
      u = &prefix_first_line;
    else
      u = &prefix;

    p0 = p2 + 1;
  }

  return t;
}


template <class T>
T convert_cast(const String &s)
{
  /* %TRACE% */  /* %TRACE% */
  std::istringstream is(s.c_str());
  T t = 0;

  is >> t;

  if (!is) {
    std::ostringstream msg;
    msg << "Unable to convert \"" << s << "\" to type " << typeid(T).name();
    throw std::runtime_error(msg.str().c_str());
  }
  

  return t;
}

template double convert_cast<double>(const String &);
template float convert_cast<float>(const String &);
template int convert_cast<int>(const String &);
template unsigned convert_cast<unsigned>(const String &);
template long convert_cast<long>(const String &);
template unsigned long convert_cast<unsigned long>(const String &);


// void
// get_function_spec_parts(
//   const std::string &		spec,
//   std::string &			namespace_name,
//   std::string &			class_name,
//   std::string &			function_name,
//   std::vector<std::string> &	arglist)
// {
//   namespace_name.erase(namespace_name.begin(), namespace_name.end());
//   class_name.erase(class_name.begin(), class_name.end());
//   function_name.erase(function_name.begin(), function_name.end());
//   arglist.erase(arglist.begin(), arglist.end());

//   std::string::const_iterator it_paren = find_next_open_paren(spec.begin(), spec.end());
//   std::string::const_iterator it_func_name = find_prev_double_colon(spec.begin(), it_paren);
//   function_name = std::string(it_func_name, it_paren);
//   if (it_func_name != spec.begin()) {
//     it_func_name -= 2;
//     std::string::const_iterator it_class_name = find_prev_double_colon(spec.begin(), it_func_name);
//     class_name = std::string(it_class_name, it_func_name);
//     if (it_class_name != spec.begin()) {
//       it_class_name -= 2;
//       namespace_name = std::string(spec.begin(), it_class_name);
//     }
//   }
// }

} // namespace sierra
