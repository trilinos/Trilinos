// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_TABLE_HPP
#define TEUCHOS_TABLE_HPP

#include "Teuchos_TableDecl.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_vector.hpp"

namespace Teuchos {

/* pretty simple 2D array */
template <typename T>
Table<T>::Table() {
}

template <typename T>
Table<T>::Table(int ncols_init, int nrows_reserve):ncols(ncols_init) {
  TEUCHOS_ASSERT(0 <= ncols_init);
  reserve(data, ncols * nrows_reserve);
}

template <typename T>
void swap(Table<T>& a, Table<T>& b) {
  using std::swap;
  swap(a.data, b.data);
  swap(a.ncols, b.ncols);
}

template <typename T>
int get_nrows(Table<T> const& t) {
  TEUCHOS_DEBUG_ASSERT(t.ncols > 0);
  TEUCHOS_DEBUG_ASSERT(Teuchos::size(t.data) % t.ncols == 0);
  return Teuchos::size(t.data) / t.ncols;
}

template <typename T>
int get_ncols(Table<T> const& t) { return t.ncols; }

template <typename T>
void resize(Table<T>& t, int new_nrows, int new_ncols) {
  TEUCHOS_ASSERT(new_ncols == t.ncols); // pretty specialized right now
  Teuchos::resize(t.data, new_nrows * t.ncols);
}

template <typename T>
typename Table<T>::Ref at(Table<T>& t, int row, int col) {
  TEUCHOS_DEBUG_ASSERT(0 <= col);
  TEUCHOS_DEBUG_ASSERT(col < t.ncols);
  TEUCHOS_DEBUG_ASSERT(0 <= row);
  TEUCHOS_DEBUG_ASSERT(row < get_nrows(t));
  return Teuchos::at(t.data, row * t.ncols + col);
}

template <typename T>
typename Table<T>::ConstRef at(Table<T> const& t, int row, int col) {
  TEUCHOS_DEBUG_ASSERT(0 <= col);
  TEUCHOS_DEBUG_ASSERT(col < t.ncols);
  TEUCHOS_DEBUG_ASSERT(0 <= row);
  TEUCHOS_DEBUG_ASSERT(row < get_nrows(t));
  return Teuchos::at(t.data, row * t.ncols + col);
}

}

#endif
