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
  template<typename CrsMatType>
  class CrsRowView : public Disp {
  public:
    typedef CrsMatType::value_type   value_type;
    typedef CrsMatType::ordinal_type ordinal_type;

  private:
    ordinal_type  _offn;    // offset in columns
    ordinal_type  _n;       // column size  
    ordinal_type  _nnz;     // # of nonzeros in row
    ordinal_type *_aj;      // column index compressed format in row
    value_type   *_ax;      // values 

  public:
    KOKKOS_INLINE_FUNCTION
    ordinal_type NumCols() const                   { return _n;             }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumNonZeros() const               { return _nnz;           } 

    KOKKOS_INLINE_FUNCTION
    ordinal_type Col(const ordinal_type j) const   { return _aj[j] - _offn; }

    KOKKOS_INLINE_FUNCTION
    value_type&  Value(const ordinal_type j)       { return _ax[j];         }

    KOKKOS_INLINE_FUNCTION
    value_type   Value(const ordinal_type j) const { return _ax[j];         }

    KOKKOS_INLINE_FUNCTION    
    ordinal_type Index(const ordinal_type col) const {
      // empty; return 0
      if (_aj == NULL) return -1;
      
      // search 
      ordinal_type base_col = (_offn + col);
      ordinal_type j = (lower_bound(_aj, _aj+_nnz, base_col) - _aj);
      
      return (_aj[j] == base_col ? j : -1);
    }

    KOKKOS_INLINE_FUNCTION
    value_type get(const ordinal_type col) const {
      ordinal_type j = Index(col);
      return (j < 0 ? value_type(0) : _ax[j]);
    }

    KOKKOS_INLINE_FUNCTION
    CrsRowView(ordinal_type offn,
               ordinal_type n,
               ordinal_type nnz,
               ordinal_type *aj,
               value_type   *ax) 
      : _offn(offn),
        _n(n),
        _nnz(nnz),
        _aj(aj),
        _ax(ax) 
    { }

    KOKKOS_INLINE_FUNCTION    
    ostream& showMe(ostream &os) const {                                                
      const int w = 6;
      os << " offset = " << setw(w) << _offn << endl; 
      for (ordinal_type j=0;j<_nnz;++j) {
        os << setw(w) << Col(j) << "  "
           << _ax[j] << endl;
      }                                                                                                   
      return os;
    }
  };
}

#endif
