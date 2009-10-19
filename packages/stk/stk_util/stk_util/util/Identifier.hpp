/*--------------------------------------------------------------------*/
/*    Copyright 2008 - 2008 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_UTIL_UTIL_IDENTIFIER_H
#define STK_UTIL_UTIL_IDENTIFIER_H

#include <string>
#include <stk_util/util/ci_string.hpp>

#include <boost/unordered_set.hpp> // boost 1.36.0 unordered_set not in the tr1

namespace std {
namespace tr1 {
  using ::boost::unordered_set;
}
}

namespace stk {

class IdentifierA
{
public:
  static int compare(const char *s1, const size_t s1_length, const char *s2, const size_t s2_length);
  
  IdentifierA()
  {}

  IdentifierA(const IdentifierA &identifier)
    : m_string(identifier.m_string)
  {}

  explicit IdentifierA(const std::string &string) 
    : m_string(string)
  {}
  
  explicit IdentifierA(const char *string) 
    : m_string(string)
  {}
  
  IdentifierA &operator=(const IdentifierA &identifier) {
    if (&identifier != this)
      m_string = identifier.m_string;
    return *this;
  }
  
  IdentifierA &operator=(const std::string &string) {
    m_string = string;
    return *this;
  }
  
  IdentifierA &operator+=(const IdentifierA &identifier) {
    m_string.append(identifier.m_string);
    return *this;
  }

  IdentifierA &operator+=(const std::string &string) {
    m_string.append(string);
    return *this;
  }

  const char *c_str() const {
    return m_string.c_str();
  }
  
  const size_t length() const {
    return m_string.length();
  }
  
  operator const std::string &() const {
    return m_string;
  }

private:
  std::string           m_string;
};


IdentifierA operator+(const IdentifierA &identifier1, const IdentifierA &identifier2);

IdentifierA operator+(const IdentifierA &identifier1, const std::string &string2);

std::string operator+(const std::string &string1, const IdentifierA &identifier2);

std::ostream &operator<<(std::ostream &os, const IdentifierA &identifier);

std::istream &operator>>(std::istream &is, IdentifierA &identifier);

bool operator<(const std::string &s1, const IdentifierA &s2);

bool operator<(const IdentifierA &s1, const std::string &s2);

bool operator<(const IdentifierA &s1, const IdentifierA &s2);

bool operator<(const IdentifierA &s1, const char *s2);

bool operator==(const std::string &s1, const IdentifierA &s2);

bool operator==(const IdentifierA &s1, const std::string &s2);

bool operator==(const IdentifierA &s1, const IdentifierA &s2);

bool operator==(const IdentifierA &s1, const char *s2);

bool operator<=(const std::string &s1, const IdentifierA &s2);

bool operator<=(const IdentifierA &s1, const std::string &s2);

bool operator<=(const IdentifierA &s1, const IdentifierA &s2);

bool operator<=(const IdentifierA &s1, const char *s2);

bool operator>(const std::string &s1, const IdentifierA &s2);

bool operator>(const IdentifierA &s1, const std::string &s2);

bool operator>(const IdentifierA &s1, const IdentifierA &s2);

bool operator>(const IdentifierA &s1, const char *s2);

bool operator>=(const std::string &s1, const IdentifierA &s2);

bool operator>=(const IdentifierA &s1, const std::string &s2);

bool operator>=(const IdentifierA &s1, const IdentifierA &s2);

bool operator>=(const IdentifierA &s1, const char *s2);

bool operator!=(const std::string &s1, const IdentifierA &s2);

bool operator!=(const IdentifierA &s1, const std::string &s2);

bool operator!=(const IdentifierA &s1, const IdentifierA &s2);

bool operator!=(const IdentifierA &s1, const char *s2);




class IdentifierB : public std::string
{
public:
  static int compare(const char *s1, const size_t s1_length, const char *s2, const size_t s2_length);
  
  IdentifierB()
  {}

  IdentifierB(const IdentifierB &identifier)
    : std::string(static_cast<std::string>(identifier))
  {}

  explicit IdentifierB(const std::string &string) 
    : std::string(string)
  {}
  
  explicit IdentifierB(const char *string) 
    : std::string(string)
  {}
  
  IdentifierB &operator=(const IdentifierB &identifier) {
    std::string::operator=(identifier);
    return *this;
  }

  IdentifierB &operator=(const std::string &string) {
    std::string::operator=(string);
    return *this;
  }

  IdentifierB &operator=(const char *string) {
    std::string::operator=(string);
    return *this;
  }
};

IdentifierB operator+(const IdentifierB &identifier1, const IdentifierB &identifier2);

IdentifierB operator+(const IdentifierB &identifier1, const std::string &string2);

IdentifierB operator+(const IdentifierB &identifier1, const char *string2);

std::string operator+(const std::string &string1, const IdentifierB &identifier2);

std::ostream &operator<<(std::ostream &os, const IdentifierB &identifier);

std::istream &operator>>(std::istream &is, IdentifierB &identifier);

bool operator<(const std::string &s1, const IdentifierB &s2);

bool operator<(const IdentifierB &s1, const std::string &s2);

bool operator<(const IdentifierB &s1, const IdentifierB &s2);

bool operator<(const IdentifierB &s1, const char *s2);

bool operator==(const std::string &s1, const IdentifierB &s2);

bool operator==(const IdentifierB &s1, const std::string &s2);

bool operator==(const IdentifierB &s1, const IdentifierB &s2);

bool operator==(const IdentifierB &s1, const char *s2);

bool operator<=(const std::string &s1, const IdentifierB &s2);

bool operator<=(const IdentifierB &s1, const std::string &s2);

bool operator<=(const IdentifierB &s1, const IdentifierB &s2);

bool operator<=(const IdentifierB &s1, const char *s2);

bool operator>(const std::string &s1, const IdentifierB &s2);

bool operator>(const IdentifierB &s1, const std::string &s2);

bool operator>(const IdentifierB &s1, const IdentifierB &s2);

bool operator>(const IdentifierB &s1, const char *s2);

bool operator>=(const std::string &s1, const IdentifierB &s2);

bool operator>=(const IdentifierB &s1, const std::string &s2);

bool operator>=(const IdentifierB &s1, const IdentifierB &s2);

bool operator>=(const IdentifierB &s1, const char *s2);

bool operator!=(const std::string &s1, const IdentifierB &s2);

bool operator!=(const IdentifierB &s1, const std::string &s2);

bool operator!=(const IdentifierB &s1, const IdentifierB &s2);

bool operator!=(const IdentifierB &s1, const char *s2);

} // namespace stk

namespace boost {

template <>
struct hash<stk::IdentifierA>
{
  std::size_t operator()(const stk::IdentifierA &s) const;
};
  
template <>
struct hash<stk::IdentifierB>
{
  std::size_t operator()(const stk::IdentifierB &s) const;
};
  
} // namespace boost

#endif // STK_UTIL_UTIL_IDENTIFIER_H
