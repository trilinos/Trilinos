#ifndef TEUCHOS_TABLE_DECL_HPP
#define TEUCHOS_TABLE_DECL_HPP

#include <vector>

namespace Teuchos {

template <typename T>
struct Table {
  std::vector<T> data;
  int ncols;
  typedef typename std::vector<T>::reference Ref;
  typedef typename std::vector<T>::const_reference ConstRef;
  Table();
  Table(int ncols_init, int nrows_reserve);
};

}

#endif
