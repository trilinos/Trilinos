// Copyright(C) 1999-2020, 2022, 2023, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#pragma once

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <cctype>                                   // for toupper
#include <cstddef>                                  // for size_t
#include <algorithm>                                 // for remove, etc
#include <iterator>                                  // for insert_iterator
#include <set>                                       // for set
#include <string>                                    // for basic_string, etc
#include <utility>                                   // for pair
#include <vector>                                    // for vector
#include <sstream>                       // for ostringstream
#include <iostream>
#include <stdexcept>
#include <numeric>

#if defined(_WIN32) && !defined(__MINGW32__)
#include <string.h>
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
namespace Iotm {
  namespace text_mesh {

    template <class EXCEPTION> void handle_error(const std::ostringstream &message)
    {
      throw EXCEPTION((message).str());
    }

    inline void default_error_handler(const std::ostringstream &message)
    {
      handle_error<std::logic_error>(message);
    }

    template <class ForwardIt, class T>
    ForwardIt bound_search(ForwardIt first, ForwardIt last, const T &value)
    {
      first = std::lower_bound(first, last, value);
      if (!(first == last) && !(value < *first))
        return first;

      return last;
    }

    template <class ForwardIt, class T, class Compare>
    ForwardIt bound_search(ForwardIt first, ForwardIt last, const T &value, Compare comp)
    {
      first = std::lower_bound(first, last, value, comp);
      if (!(first == last) && !(comp(value, *first)))
        return first;

      return last;
    }

    inline std::string strip_whitespace(const std::string &inpt)
    {
      auto start_it = inpt.begin();
      auto end_it   = inpt.rbegin();
      while (std::isspace(*start_it))
        ++start_it;
      while (std::isspace(*end_it))
        ++end_it;
      return std::string(start_it, end_it.base());
    }

    inline std::vector<std::string> get_tokens(const std::string &str,
                                               const std::string &separators)
    {
      std::vector<std::string> tokens;
      auto                     first = std::begin(str);
      while (first != std::end(str)) {
        const auto second =
            std::find_first_of(first, std::end(str), std::begin(separators), std::end(separators));
        if (first != second) {
          std::string token = strip_whitespace(std::string(first, second));
          tokens.emplace_back(token);
        }
        if (second == std::end(str)) {
          break;
        }
        first = std::next(second);
      }
      return tokens;
    }

    inline void convert_to_uppercase(std::string &str)
    {
      std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    }

    inline void convert_to_lowercase(std::string &str)
    {
      std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    }

    inline bool is_positive_number(const std::string &str)
    {
      for (char const &c : str) {
        if (std::isdigit(c) == 0)
          return false;
      }
      return true;
    }

    template <typename T> std::set<T> transform_to_set(const std::vector<T> &dataAsVector)
    {
      std::set<T> dataAsSet;

      for (const T &data : dataAsVector) {
        dataAsSet.insert(data);
      }

      return dataAsSet;
    }

    inline std::pair<unsigned, bool> get_id_from_part_name(const std::string &name,
                                                           const std::string &prefix)
    {
      const unsigned prefixLength = prefix.length();

      if (name.length() < prefixLength + 1)
        return std::make_pair(0, false);

      const std::string namePrefix = name.substr(0, prefixLength);
      const std::string nameSuffix = name.substr(prefixLength);

      if (strcasecmp(namePrefix.c_str(), prefix.c_str()) != 0)
        return std::make_pair(0, false);

      unsigned           id;
      std::istringstream nameSuffixStream(nameSuffix);
      nameSuffixStream >> id;
      if (nameSuffixStream.fail()) {
        return std::make_pair(0, false);
      }
      return std::make_pair(id, true);
    }

  } // namespace text_mesh
} // namespace Iotm
