// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_STACK_HPP
#define TEUCHOS_STACK_HPP

#include <stack>

namespace Teuchos {

template <typename T>
int size(std::stack<T> const& s) { return static_cast<int>(s.size()); }

}

#endif
