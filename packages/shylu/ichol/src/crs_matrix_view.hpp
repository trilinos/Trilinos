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

  template<typename CrsMatBaseType>  
  class CrsRowView;

  template<typename CrsMatBaseType>
  class CrsMatrixView : public Disp {
  public:
    typedef typename CrsMatBaseType::space_type    space_type;
    typedef typename CrsMatBaseType::memory_traits memory_traits;
    
    typedef typename CrsMatBaseType::value_type    value_type;
    typedef typename CrsMatBaseType::ordinal_type  ordinal_type;
    typedef typename CrsMatBaseType::size_type     size_type;

    typedef CrsMatBaseType             mat_base_type;
    typedef CrsRowView<mat_base_type>  row_view_type;

    // be careful this use rcp and atomic operation
    // - use setView to create a view if _rows is not necessary
    // - copy constructor and assignment operator will do soft copy of the object
    typedef Kokkos::View<row_view_type*,space_type,memory_traits> row_view_type_array;
    
  private:
    CrsMatBaseType *_base;   // pointer to the base object

    ordinal_type  _offm;     // offset in rows
    ordinal_type  _offn;     // offset in cols
    ordinal_type  _m;        // # of rows
    ordinal_type  _n;        // # of cols

    row_view_type_array _rows;
    
  public:
    KOKKOS_INLINE_FUNCTION
    void fillRowViewArray(const bool flag = true) {
      if (flag) {
        if (_rows.dimension_0() < _m)
          _rows = row_view_type_array(_base->Label() + "::View::RowViewArray", _m);
        
        for (ordinal_type i=0;i<_m;++i)
          _rows[i].setView(*this, i);
      } else {
        _rows = row_view_type_array();
      }
    }

    KOKKOS_INLINE_FUNCTION
    row_view_type& RowView(const ordinal_type i) const { return _rows[i]; }

    KOKKOS_INLINE_FUNCTION
    void setView(CrsMatBaseType *base,
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

    CrsMatrixView()
      : _base(NULL),
        _offm(0),
        _offn(0),
        _m(0),
        _n(0),
        _rows()
    { } 

    CrsMatrixView(const CrsMatrixView &b)
      : _base(b._base),
        _offm(b._offm),
        _offn(b._offn),
        _m(b._m),
        _n(b._n),
        _rows(b._rows)
    { } 

    CrsMatrixView(CrsMatBaseType *b)
      : _base(b),
        _offm(0),
        _offn(0),
        _m(b->NumRows()),
        _n(b->NumCols())
    { } 

    CrsMatrixView(CrsMatBaseType *b,
                  const ordinal_type offm, const ordinal_type m,
                  const ordinal_type offn, const ordinal_type n) 
      : _base(b),
        _offm(offm),
        _offn(offn),
        _m(m),
        _n(n) 
    { } 

    ostream& showMe(ostream &os) const {
      const int w = 4;
      if (_base != NULL) 
        os << _base->Label() << "::View, "
           << " Offs ( " << setw(w) << _offm << ", " << setw(w) << _offn << " ); "
           << " Dims ( " << setw(w) << _m    << ", " << setw(w) << _n    << " ); ";
      else 
        os << "-- Base object is null --";
      
      return os;
    }

  };
}

#endif
