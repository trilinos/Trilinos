// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "smart_assert.h" // for SMART_ASSERT
#include "stringx.h"
#include <cctype>  // for tolower, isspace
#include <cstring> // for strspn, strcspn
#include <string>  // for string, operator==
#include <vector>  // for vector

bool abbreviation(const std::string &s, const std::string &master, unsigned min_length)
{
  SMART_ASSERT(min_length > 0);

  if (s.size() > master.size()) {
    return false;
  }

  if (s.size() < min_length) {
    return false;
  }

  for (unsigned i = 0; i < s.size(); ++i) {
    if (s[i] != master[i]) {
      return false;
    }
  }
  return true;
}

bool no_case_equals(const std::string &s1, const std::string &s2)
{
  if (s1.size() != s2.size()) {
    return false;
  }

  for (unsigned i = 0; i < s1.size(); ++i) {
    if (tolower(s1[i]) != tolower(s2[i])) {
      return false;
    }
  }
  return true;
}

std::string &chop_whitespace(std::string &s)
{
  if (!s.empty()) {
    int i = (int)s.size() - 1;
    for (; i >= 0; --i) {
      if (isspace(static_cast<int>(s[i])) == 0) {
        break;
      }

      s.resize(i + 1);
    }
  }
  return s;
}

std::string extract_token(std::string &s, const char *delimiters)
{
  if (!s.empty()) {
    SMART_ASSERT(delimiters != nullptr && !std::string(delimiters).empty());

    // Move through initial delimiters.
    auto p = s.find_first_not_of(delimiters);

    if (p >= s.size()) {
      // no tokens
      s = "";
      return "";
    }
    if (s[p] == '"') {
      // Special case of a quoted variable name which likely contains
      // whitespace, but it should work for any quoted variable
      // name. Some of this is a bit redundant but it makes this section
      // completely self contained so that this code is only executed
      // when a quoted variable name is found and it doesn't affect any
      // action outside of this block of code.

      // Find the closing quote
      auto cq = s.find_first_of('\"', p + 1);

      // No closing quote found. Error out.
      SMART_ASSERT(cq < s.size());

      std::string tok = s.substr(p + 1, cq - (p + 1));

      auto r = s.find_first_not_of(delimiters, cq + 1);

      if (r >= s.size()) {
        s = "";
      }
      else {
        s.erase(0, r);
      }
      return tok;
    }

    // move to end of first token
    auto q = s.find_first_of(delimiters, p);

    if (q >= s.size()) {
      // no more delimiters
      std::string tok = s.substr(p);
      s               = "";
      return tok;
    }

    SMART_ASSERT(q > p);
    std::string tok = s.substr(p, q - p);

    // move to start of the second token
    auto r = s.find_first_not_of(delimiters, q);

    if (r >= s.size()) {
      // no second token
      s = "";
    }
    else {
      s.erase(0, r);
    }

    return tok;
  }

  return "";
}

int count_tokens(const std::string &s, const char *delimiters)
{
  if (!s.empty()) {
    const char *str_ptr = s.c_str();

    // Move through initial delimiters.
    const char *p = &str_ptr[strspn(str_ptr, delimiters)];

    int num_toks = 0;
    while (p[0] != '\0') {
      ++num_toks;
      p = &p[strcspn(p, delimiters)]; // Move through token.
      p = &p[strspn(p, delimiters)];  // Move through delimiters.
    }

    return num_toks;
  }

  return 0;
}

int max_string_length(const std::vector<std::string> &names)
{
  if (names.empty()) {
    return 0;
  }
  auto len = names[0].size();
  for (size_t i = 1; i < names.size(); i++) {
    if (names[i].size() > len) {
      len = names[i].size();
    }
  }
  return len;
}

void to_lower(std::string &s)
{
  for (auto &elem : s) {
    elem = (char)tolower(elem);
  }
}

char first_character(const std::string &s)
{
  for (const auto &elem : s) {
    if (isspace(static_cast<int>(elem)) == 0) {
      return elem;
    }
  }
  return 0;
}

int find_string(const std::vector<std::string> &lst, const std::string &s, bool nocase)
{
  if (nocase) {
    for (unsigned i = 0; i < lst.size(); ++i) {
      if (no_case_equals(lst[i], s)) {
        return i;
      }
    }
  }
  else {
    for (unsigned i = 0; i < lst.size(); ++i) {
      if (lst[i] == s) {
        return i;
      }
    }
  }
  return -1;
}
