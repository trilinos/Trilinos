//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSKERNELS_TESTSTRINGUTILS_HPP
#define KOKKOSKERNELS_TESTSTRINGUTILS_HPP

#include <string>
#include <cstring>

namespace Test {

inline int string_compare_no_case(const char* str1, const char* str2) {
  std::string str1_s(str1);
  std::string str2_s(str2);
  for (size_t i = 0; i < str1_s.size(); i++) str1_s[i] = std::tolower(str1_s[i]);
  for (size_t i = 0; i < str2_s.size(); i++) str2_s[i] = std::tolower(str2_s[i]);
  return std::strcmp(str1_s.c_str(), str2_s.c_str());
}

inline int string_compare_no_case(const std::string& str1, const std::string& str2) {
  return string_compare_no_case(str1.c_str(), str2.c_str());
}

}  // namespace Test
#endif
