#pragma once
#ifndef __CRS_ROW_VIEW_HPP__
#define __CRS_ROW_VIEW_HPP__

namespace Example { 

  using namespace std;

  template<typename ValueType, 
           typename OrdinalType>
  class CrsRowView {
  public:
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;

  private:
    ordinal_type  _nnz;     // # of nonzeros in row
    ordinal_type *_aj;      // column index compressed format in row
    value_type   *_ax;      // values 

  public:
    ordinal_type NumNonZerosInRow() const    { return _nnz;   }
    ordinal_type ColIndex(const int j) const { return _aj[j]; }
    value_type   Value(const int j) const    { return _ax[j]; }

    value_type at(int col) const {
      auto j = lower_bound(_aj, _aj+_nnz, col);
      return (_aj[j] == col ? _ax[j] : value_type(0));
    }
    
    CrsRowView(const ordinal_type nnz,
               const ordinal_type *aj,
               const value_type   *ax) 
      : _nnz(nnz),
        _aj(aj),
        _ax(ax) 
    { }

  };
}

#endif
