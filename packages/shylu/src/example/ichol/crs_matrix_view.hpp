#pragma once
#ifndef __CRS_MATRIX_VIEW_HPP__
#define __CRS_MATRIX_VIEW_HPP__

/// \file crs_matrix_view.hpp
/// \brief CRS matrix view object creates 2D view to setup a computing region.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <Kokkos_Core.hpp>

#include "util.hpp"

namespace Example { 

  using namespace std;

  // forward declaration
  template <typename CrsMatBaseType>
  class CrsRowView;

  template<typename CrsMatBaseType>
  class CrsMatrixView : public Disp {
  public:
    typedef typename CrsMatBaseType::value_type    value_type;
    typedef typename CrsMatBaseType::ordinal_type  ordinal_type;
    typedef typename CrsMatBaseType::space_type    space_type;
    typedef typename CrsMatBaseType::memory_traits memory_traits;

    typedef CrsRowView<CrsMatBaseType> crs_row_view_type;

    template<typename T> using range_type = pair<T,T>;

  private:
    ordinal_type  _offm;    // offset in rows
    ordinal_type  _offn;    // offset in cols
    ordinal_type  _m;       // # of rows
    ordinal_type  _n;       // # of cols

    CrsMatBaseType *_base;   // pointer to the base object
    
  public:
    //
    KOKKOS_INLINE_FUNCTION
    void setView(const CrsMatBaseType *base,
                 const ordinal_type offm, const ordinal_type m,
                 const ordinal_type offn, const ordinal_type n) {
      _base = base;

      _offm = offm; _m = m;
      _offn = offn; _n = n;
    }

    KOKKOS_INLINE_FUNCTION
    CrsMatBaseType* BaseObject() const { return _base; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type  OffsetRows() const { return _offm; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type  OffsetCols() const { return _offn; }

    KOKKOS_INLINE_FUNCTION    
    ordinal_type  NumRows() const { return _m; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type  NumCols() const { return _n; }

    KOKKOS_INLINE_FUNCTION
    crs_row_view_type extractRow(const ordinal_type i) const { 
      typedef typename crs_row_view_type::value_type_array   value_type_array;
      typedef typename crs_row_view_type::ordinal_type_array ordinal_type_array;

      // grep a row from base
      ordinal_type ii = _offm + i;  
      ordinal_type_array cols = _base->ColsInRow(ii);
      value_type_array   vals = _base->ValuesInRow(ii);

      ordinal_type *ptr_begin = cols.ptr_on_device();
      ordinal_type *ptr_end   = cols.ptr_on_device() + cols.dimension_0();

      ordinal_type view_begin = lower_bound(ptr_begin, ptr_end, _offn     ) - ptr_begin; 
      ordinal_type view_end   = upper_bound(ptr_begin, ptr_end, _offn+_n-1) - ptr_begin; 

      range_type<ordinal_type> range(view_begin, view_end);

      cols = Kokkos::subview<ordinal_type_array>(cols, range);
      vals = Kokkos::subview<value_type_array>(vals, range);

      return crs_row_view_type(_offn, _n, view_end - view_begin, cols, vals);
    }

    CrsMatrixView()
      : _base(NULL),
        _offm(0),
        _offn(0),
        _m(0),
        _n(0)
    { } 

    CrsMatrixView(const CrsMatrixView &b)
      : _base(b._base),
        _offm(b._offm),
        _offn(b._offn),
        _m(b._m),
        _n(b._n)
    { } 
    
    CrsMatrixView(CrsMatBaseType &b)
      : _base(&b),
        _offm(0),
        _offn(0),
        _m(b.NumRows()),
        _n(b.NumCols())
    { } 

    CrsMatrixView(CrsMatBaseType &b,
                  const ordinal_type offm, const ordinal_type m,
                  const ordinal_type offn, const ordinal_type n) 
      : _base(&b),
        _offm(offm),
        _offn(offn),
        _m(m),
        _n(n) 
    { } 

    ostream& showMe(ostream &os) const {

      if (_base != NULL) {
        os << endl
           << " -- CrsMatrixView::" << _base->Label() << " --" << endl
           << "    Offset in Rows = " << _offm << endl
           << "    Offset in Cols = " << _offn << endl
           << "    # of Rows      = " << _m << endl
           << "    # of Cols      = " << _n << endl;
        
        const int w = 6;
        for (ordinal_type i=0;i<_m;++i) {
          os << endl;
          crs_row_view_type row = extractRow(i);
          row.showMe(os);
        }
      } else {
        os << endl
           << " Base object is null " << endl;
      }
      
      return os;
    }

  };
}

#endif
