// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_SET_HPP
#define TEUCHOS_SET_HPP

#include <set>
#include <iosfwd>
#include <Teuchos_Assert.hpp>

namespace Teuchos {

template <typename T>
void unite_with(std::set<T>& a, std::set<T> const& b) {
  for (typename std::set<T>::const_iterator it = b.begin(); it != b.end(); ++it) {
    const T& x = *it;
    a.insert(x);
  }
}

template <typename T>
void unite(std::set<T>& result, std::set<T> const& a, std::set<T> const& b) {
  result = a;
  unite_with(result, b);
}

template <typename T>
void subtract_from(std::set<T>& a, std::set<T> const& b) {
  for (typename std::set<T>::const_iterator it = b.begin(); it != b.end(); ++it) {
    const T& x = *it;
    typename std::set<T>::iterator it2 = a.find(x);
    if (it2 == a.end()) continue;
    a.erase(it2);
  }
}

template <typename T>
bool intersects(std::set<T> const& a, std::set<T> const& b) {
  for (typename std::set<T>::const_iterator it = b.begin(); it != b.end(); ++it) {
    const T& x = *it;
    typename std::set<T>::const_iterator it2 = a.find(x);
    if (it2 != a.end()) return true;
  }
  return false;
}

}

#endif
