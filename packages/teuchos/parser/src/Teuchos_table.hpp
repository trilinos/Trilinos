#ifndef TEUCHOS_TABLE_HPP
#define TEUCHOS_TABLE_HPP

#include <Teuchos_vector.hpp>

namespace Teuchos {

/* pretty simple 2D array */
template <typename T>
struct Table {
  std::vector<T> data;
  int ncols;
  using Ref = typename std::vector<T>::reference;
  using ConstRef = typename std::vector<T>::const_reference;
  Table() = default;
  Table(int ncols_init, int nrows_reserve):ncols(ncols_init) {
    assert(0 <= ncols_init);
    reserve(data, ncols * nrows_reserve);
  }
};

template <typename T>
int get_nrows(Table<T> const& t) {
  assert(t.ncols > 0);
  assert(size(t.data) % t.ncols == 0);
  return size(t.data) / t.ncols;
}

template <typename T>
int get_ncols(Table<T> const& t) { return t.ncols; }

template <typename T>
void resize(Table<T>& t, int new_nrows, int new_ncols) {
  assert(new_ncols == t.ncols); // pretty specialized right now
  Teuchos::resize(t.data, new_nrows * t.ncols);
}

template <typename T>
typename Table<T>::Ref at(Table<T>& t, int row, int col) {
  assert(0 <= col);
  assert(col < t.ncols);
  assert(0 <= row);
  assert(row < get_nrows(t));
  return Teuchos::at(t.data, row * t.ncols + col);
}

template <typename T>
typename Table<T>::ConstRef at(Table<T> const& t, int row, int col) {
  assert(0 <= col);
  assert(col < t.ncols);
  assert(0 <= row);
  assert(row < get_nrows(t));
  return Teuchos::at(t.data, row * t.ncols + col);
}

}

#endif

