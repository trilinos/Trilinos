// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_VECTOR_HPP
#define TEUCHOS_VECTOR_HPP

#include <vector>
#include <Teuchos_Assert.hpp>

namespace Teuchos {

/* just some wrappers over std::vector to let us
   do all indexing with int */

template <typename T>
inline int size(std::vector<T> const& v) {
  return int(v.size());
}

template <typename T>
inline typename std::vector<T>::reference at(std::vector<T>& v, int i) {
  TEUCHOS_DEBUG_ASSERT(0 <= i);
  TEUCHOS_DEBUG_ASSERT(i < int(v.size()));
  return v[std::size_t(i)];
}

template <typename T>
inline typename std::vector<T>::const_reference at(std::vector<T> const& v, int i) {
  TEUCHOS_DEBUG_ASSERT(0 <= i);
  TEUCHOS_DEBUG_ASSERT(i < int(v.size()));
  return v[std::size_t(i)];
}

template <typename T>
inline void resize(std::vector<T>& v, int n) {
  TEUCHOS_DEBUG_ASSERT(0 <= n);
  v.resize(std::size_t(n));
}

template <typename T>
inline void reserve(std::vector<T>& v, int n) {
  TEUCHOS_DEBUG_ASSERT(0 <= n);
  v.reserve(std::size_t(n));
}

template <typename T>
inline std::vector<T> make_vector(int n, T const& init_val = T()) {
  return std::vector<T>(std::size_t(n), init_val);
}

template <typename T>
void add_back(std::vector<T>& vector, T& value) {
  using std::swap;
  vector.push_back(T());
  swap(vector.back(), value);
}

}

#endif

