// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include <array>
#include <cstdlib>
#include <vector>

// We have been relying on the assumption that calling `.data()` on an empty vector
// will return `nullptr`.  However, according to cppreference (based on the standard):
//
// `If size() is ​0​, data() may or may not return a null pointer.`
//
// We don't have any systems on which we have found that (yet?), but this is proactive
// in removing our use of `.data()` on potentially empty vectors...
template <typename T> constexpr T *Data(std::vector<T> &vec)
{
  if (vec.empty()) {
    return nullptr;
  }
  return vec.data();
}

template <typename T> constexpr const T *Data(const std::vector<T> &vec)
{
  if (vec.empty()) {
    return nullptr;
  }
  return vec.data();
}

template <typename T, size_t N> constexpr T *Data(std::array<T, N> &arr)
{
  return N == 0 ? nullptr : arr.data();
}

template <typename T, size_t N> constexpr const T *Data(const std::array<T, N> &arr)
{
  return N == 0 ? nullptr : arr.data();
}
