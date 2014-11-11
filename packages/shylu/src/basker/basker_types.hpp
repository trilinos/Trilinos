#ifndef BASKER_TYPES_HPP
#define BASKER_TYPES_HPP


template <class Int, class Entry>
struct basker_matrix
{
  Int nrow, ncol, nnz;
  Int *col_ptr, *row_idx;
  Entry *val;
  
};

#endif
