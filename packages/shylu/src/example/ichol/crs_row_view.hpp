#pragma once
#ifndef __CRS_ROW_VIEW_HPP__
#define __CRS_ROW_VIEW_HPP__

#include <Kokkos_Core.hpp>

/// \file crs_row_view.hpp
/// \brief A view to a row extracted from CrsMatrixView.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  /// \class CrsRowView
  template<typename CrsMatBaseType>
  class CrsRowView : public Disp {
  public:
    typedef typename CrsMatBaseType::ordinal_type            ordinal_type;
    typedef typename CrsMatBaseType::value_type              value_type;
    typedef typename CrsMatBaseType::ordinal_type_array_view ordinal_type_array_view;
    typedef typename CrsMatBaseType::value_type_array_view   value_type_array_view;
    
  private:
    ordinal_type _offn;    // offset in columns
    ordinal_type _n;       // column size  
    ordinal_type _nnz;     // # of nonzeros in row

    // this assumes a contiguous memory buffer
    ordinal_type_array_view _aj;      // column index compressed format in row
    value_type_array_view   _ax;      // values 

  public:
    KOKKOS_INLINE_FUNCTION
    ordinal_type NumCols() const { return _n; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumNonZeros() const { return _nnz; } 

    KOKKOS_INLINE_FUNCTION
    ordinal_type Col(const ordinal_type j) const { return _aj[j] - _offn; }

    KOKKOS_INLINE_FUNCTION
    value_type& Value(const ordinal_type j) { return _ax[j]; }

    KOKKOS_INLINE_FUNCTION
    value_type Value(const ordinal_type j) const { return _ax[j]; }
    
    KOKKOS_INLINE_FUNCTION    
    ordinal_type Index(const ordinal_type col) const {
      if (_aj.dimension_0() == 0)
        return -1;
      
      // Qeustion : parallel scan ? 
      ordinal_type base_col = _offn + col;
      ordinal_type *ptr = _aj.ptr_on_device();
      ordinal_type j = lower_bound(ptr, ptr+_nnz, base_col) - ptr;
      
      return (_aj[j] == base_col ? j : -1);
    }

    KOKKOS_INLINE_FUNCTION
    value_type ValueAtColumn(const ordinal_type col) const {
      ordinal_type j = Index(col);
      return (j < 0 ? value_type(0) : _ax[j]);
    }

    CrsRowView()
      : _offn(0),
        _n(0),
        _nnz(0),
        _aj(),
        _ax() 
    { }

    CrsRowView(const ordinal_type            offn,
               const ordinal_type            n,
               const ordinal_type            nnz,
               const ordinal_type_array_view aj,
               const value_type_array_view   ax) 
      : _offn(offn),
        _n(n),
        _nnz(nnz),
        _aj(aj),
        _ax(ax) 
    { }

    ostream& showMe(ostream &os) const {                                                
      os << "  offset = " << _offn
         << ", nnz = " << _nnz 
         << endl; 
      for (ordinal_type j=0;j<_nnz;++j) {
        const value_type val = _ax[j];
        os << "(" << Col(j) << ", "
           << val << ")"
           << endl;
      }                                                                                                   
      return os;
    }
  };
}

#endif
