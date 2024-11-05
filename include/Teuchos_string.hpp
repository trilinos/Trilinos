// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_STRING_HPP
#define TEUCHOS_STRING_HPP

#include <string>
#include <Teuchos_Assert.hpp>

namespace Teuchos {

/* some wrappers over std::string to let us
   do all indexing with int */

inline int size(std::string const& s) {
  return int(s.size());
}

inline std::string::reference at(std::string& s, int i) {
  TEUCHOS_DEBUG_ASSERT(0 <= i);
  TEUCHOS_DEBUG_ASSERT(i < int(s.size()));
  return s[std::size_t(i)];
}

inline std::string::const_reference at(std::string const& s, int i) {
  TEUCHOS_DEBUG_ASSERT(0 <= i);
  TEUCHOS_DEBUG_ASSERT(i < int(s.size()));
  return s[std::size_t(i)];
}

inline void resize(std::string& s, int n) {
  TEUCHOS_DEBUG_ASSERT(0 <= n);
  s.resize(std::size_t(n));
}

}

#endif


