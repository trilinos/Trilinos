#ifndef _Reporter_cpp_
#define _Reporter_cpp_

#include "Reporter.hpp"

namespace ngp_testing {

NGP_TEST_INLINE void bounded_strcpy(const char* src, const int maxChar, char* dest) {
  int idx = 0;
  if(src) {
    while(src[idx] != '\0' && idx < maxChar) {
      dest[idx] = src[idx];
      ++idx;
    }
  }
  dest[idx] = '\0';
}

NGP_TEST_INLINE
TruncatedString::TruncatedString(const char* src) {
  const int srcLen = get_string_length(src);
  if(should_truncate(srcLen)) {
    src += get_truncation_offset(srcLen);
  }
  bounded_strcpy(src, maxNumChar, string);
  if(should_truncate(srcLen)) {
    prepend_truncation_indicator(string);
  }
}

NGP_TEST_INLINE
int TruncatedString::get_string_length(const char* str) const {
  int len = 0;
  if(str) {
    while (str[len] != '\0') ++len;
  }
  return len;
}

NGP_TEST_INLINE
bool TruncatedString::should_truncate(const int strLen) const {
  return strLen > maxNumChar;
}

NGP_TEST_INLINE
int TruncatedString::get_truncation_offset(const int strLen) const {
  return strLen - maxNumChar + 1;
}

NGP_TEST_INLINE
void TruncatedString::prepend_truncation_indicator(char* dest) {
  dest[0] = '.';
  dest[1] = '.';
  dest[2] = '.';
}

}

#endif

