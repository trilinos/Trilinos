#include <stk_util/util/Identifier.hpp>
#include <cstring>

namespace stk {

namespace {

int
compare(
  const char *  s1,
  const size_t  s1_length,
  const char *  s2,
  const size_t  s2_length)
{
  const size_t length = std::min(s1_length, s2_length);
  int result = ignorecase_traits::compare(s1, s2, length);
  if (!result)
    result = s1_length - s2_length;
  return result;
}

} // namespace <empty>


int
IdentifierA::compare(
  const char *  s1,
  const size_t  s1_length,
  const char *  s2,
  const size_t  s2_length)
{
  const size_t length = std::min(s1_length, s2_length);
  int result = ignorecase_traits::compare(s1, s2, length);
  if (!result)
    result = s1_length - s2_length;
  return result;
}


IdentifierA operator+(const IdentifierA &identifier1, const IdentifierA &identifier2) {
  IdentifierA identifier(identifier1);
    
  return identifier += identifier2;
}

IdentifierA operator+(const IdentifierA &identifier1, const std::string &string2) {
  IdentifierA identifier(identifier1);
    
  return identifier += IdentifierA(string2);
}

std::string operator+(const std::string &string1, const IdentifierA &identifier2) {
  std::string string(string1);
    
  return string += identifier2;
}

std::ostream &operator<<(std::ostream &os, const IdentifierA &identifier) {
  return os << identifier.c_str();
}

std::istream &operator>>(std::istream &is, IdentifierA &identifier) {
  std::string s;

  is >> s;
  identifier = s;
  
  return is;
}

bool operator<(const std::string &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) < 0;
}

bool operator<(const IdentifierA &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) < 0;
}

bool operator<(const IdentifierA &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) < 0;
}

bool operator<(const IdentifierA &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) < 0;
}

bool operator==(const std::string &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) == 0;
}

bool operator==(const IdentifierA &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) == 0;
}

bool operator==(const IdentifierA &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) == 0;
}

bool operator==(const IdentifierA &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) == 0;
}

bool operator<=(const std::string &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) <= 0;
}

bool operator<=(const IdentifierA &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) <= 0;
}

bool operator<=(const IdentifierA &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) <= 0;
}

bool operator<=(const IdentifierA &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) <= 0;
}

bool operator>(const std::string &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) > 0;
}

bool operator>(const IdentifierA &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) > 0;
}

bool operator>(const IdentifierA &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) > 0;
}

bool operator>(const IdentifierA &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) > 0;
}

bool operator>=(const std::string &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) >= 0;
}

bool operator>=(const IdentifierA &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) >= 0;
}

bool operator>=(const IdentifierA &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) >= 0;
}

bool operator>=(const IdentifierA &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) >= 0;
}

bool operator!=(const std::string &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) != 0;
}

bool operator!=(const IdentifierA &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) != 0;
}

bool operator!=(const IdentifierA &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) != 0;
}

bool operator!=(const IdentifierA &s1, const IdentifierA &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) != 0;
}



int
IdentifierB::compare(
  const char *  s1,
  const size_t  s1_length,
  const char *  s2,
  const size_t  s2_length)
{
  const size_t length = std::min(s1_length, s2_length);
  int result = ignorecase_traits::compare(s1, s2, length);
  if (!result)
    result = s1_length - s2_length;
  return result;
}


std::ostream &operator<<(std::ostream &os, const IdentifierB &identifier) {
  return os << identifier.c_str();
}

std::istream &operator>>(std::istream &is, IdentifierB &identifier) {
  std::string s;

  is >> s;
  identifier = s;
  
  return is;
}

IdentifierB operator+(const IdentifierB &identifier1, const IdentifierB &identifier2) {
  std::string identifier(identifier1);
    
  return IdentifierB(identifier += identifier2);
}

IdentifierB operator+(const IdentifierB &identifier1, const std::string &string2) {
  std::string identifier(identifier1);
    
  return IdentifierB(identifier += string2);
}

IdentifierB operator+(const IdentifierB &identifier1, const char *string2) {
  IdentifierB identifier(identifier1);
    
  return IdentifierB(identifier += string2);
}

std::string operator+(const std::string &string1, const IdentifierB &identifier2) {
  std::string string(string1);
    
  return string += identifier2;
}

bool operator<(const IdentifierB &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) < 0;
}

bool operator<(const IdentifierB &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) < 0;
}

bool operator<(const IdentifierB &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) < 0;
}

bool operator==(const std::string &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) == 0;
}

bool operator==(const IdentifierB &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) == 0;
}

bool operator==(const IdentifierB &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) == 0;
}

bool operator==(const IdentifierB &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) == 0;
}

bool operator<=(const std::string &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) <= 0;
}

bool operator<=(const IdentifierB &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) <= 0;
}

bool operator<=(const IdentifierB &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) <= 0;
}

bool operator<=(const IdentifierB &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) <= 0;
}

bool operator>(const std::string &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) > 0;
}

bool operator>(const IdentifierB &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) > 0;
}

bool operator>(const IdentifierB &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) > 0;
}

bool operator>(const IdentifierB &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) > 0;
}

bool operator>=(const std::string &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) >= 0;
}

bool operator>=(const IdentifierB &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) >= 0;
}

bool operator>=(const IdentifierB &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) >= 0;
}

bool operator>=(const IdentifierB &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) >= 0;
}

bool operator!=(const std::string &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) != 0;
}

bool operator!=(const IdentifierB &s1, const std::string &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) != 0;
}

bool operator!=(const IdentifierB &s1, const char *s2) {
  return compare(s1.c_str(), s1.length(), s2, std::strlen(s2)) != 0;
}

bool operator!=(const IdentifierB &s1, const IdentifierB &s2) {
  return compare(s1.c_str(), s1.length(), s2.c_str(), s2.length()) != 0;
}

} // namespace stk
