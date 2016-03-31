#ifndef __TACHO_CRS_MATRIX_VIEW_HPP__
#define __TACHO_CRS_MATRIX_VIEW_HPP__

/// \file Tacho_CrsMatrixView.hpp
/// \brief CRS matrix view object creates 2D view to setup a computing region.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho { 


  template<typename CrsMatBaseType>
  class CrsRowView;

  template<typename CrsMatBaseType>
  class CrsMatrixView {
  public:
    typedef typename CrsMatBaseType::space_type    space_type;
    
    typedef typename CrsMatBaseType::value_type    value_type;
    typedef typename CrsMatBaseType::ordinal_type  ordinal_type;
    typedef typename CrsMatBaseType::size_type     size_type;

    typedef CrsMatBaseType             mat_base_type;
    typedef CrsRowView<mat_base_type>  row_view_type;

    // range type
    template<typename T> using range_type = Kokkos::pair<T,T>;

    // be careful this use rcp and atomic operation
    // - use setView to create a view if _rows is not necessary
    // - copy constructor and assignment operator will do soft copy of the object
    typedef Kokkos::View<row_view_type*,space_type,Kokkos::MemoryUnmanaged> row_view_type_array;
    
  private:
    CrsMatBaseType _base;   // shallow copy of the base matrix

    ordinal_type   _offm;     // offset in rows
    ordinal_type   _offn;     // offset in cols

    ordinal_type   _m;        // # of rows
    ordinal_type   _n;        // # of cols

    size_type _nnz;
    ordinal_type _lc, _rc;

    row_view_type_array _rows;

  public:

    /// Interface functions
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (o)
    KOKKOS_INLINE_FUNCTION
    void setView(const CrsMatBaseType &base,
                 const ordinal_type offm, const ordinal_type m,
                 const ordinal_type offn, const ordinal_type n) {
      _base = base;

      _offm = offm; _m = m;
      _offn = offn; _n = n;
    }

    KOKKOS_INLINE_FUNCTION
    void setRowViewArray(const row_view_type_array rows) {
      // use user provided array
      _rows = rows;

      // fill the array
      _nnz = 0; 
      const bool is_view_to_base = ( _offm == 0 && _m == _base.NumRows() &&
                                     _offn == 0 && _n == _base.NumCols() );
      for (ordinal_type i=0;i<_m;++i) {
        auto &row = _rows(i);

        if (is_view_to_base)
          row.setView(_base, i);
        else
          row.setView(*this, i);

        // statistics
        const size_type tmp = row.NumNonZeros();
        if (tmp) {
          _lc = _nnz ? Util::min(_lc, row.Col(0))     : row.Col(0);
          _rc = _nnz ? Util::max(_rc, row.Col(tmp-1)) : row.Col(tmp-1);
          _nnz += tmp;
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    row_view_type& RowView(const ordinal_type i) const { return _rows[i]; }

    KOKKOS_INLINE_FUNCTION
    const CrsMatBaseType & BaseObject() const { return _base; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type  OffsetRows() const { return _offm; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type  OffsetCols() const { return _offn; }

    KOKKOS_INLINE_FUNCTION    
    ordinal_type  NumRows() const { return _m; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type  NumCols() const { return _n; }

    KOKKOS_INLINE_FUNCTION
    size_type NumNonZeros() const { 
      return _nnz;
    }

    KOKKOS_INLINE_FUNCTION
    void ColRange(ordinal_type &lc, 
                  ordinal_type &rc) {
      lc = _lc; rc = _rc;
    }

    KOKKOS_INLINE_FUNCTION
    bool isNull() const { 
      return (_m == 0 || _n == 0);
    }

    KOKKOS_INLINE_FUNCTION
    bool isEmpty() const { 
      return (_nnz == 0);
    }

    /// ------------------------------------------------------------------

    /// Constructors
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (o)
    KOKKOS_INLINE_FUNCTION
    CrsMatrixView()
      : _base(),
        _offm(0),
        _offn(0),
        _m(0),
        _n(0),
        _nnz(0),
        _lc(-1),
        _rc(-1),
        _rows()
    { } 

    template<typename MT>
    KOKKOS_INLINE_FUNCTION
    CrsMatrixView(const CrsMatrixView<MT> &b)
      : _base(b._base),
        _offm(b._offm),
        _offn(b._offn),
        _m(b._m),
        _n(b._n),
        _nnz(b._nnz),
        _lc(b._lc),
        _rc(b._rc),
        _rows(b._rows)
    { } 

    KOKKOS_INLINE_FUNCTION
    CrsMatrixView(const CrsMatBaseType &b)
      : _base(b),
        _offm(0),
        _offn(0),
        _m(b.NumRows()),
        _n(b.NumCols()),
        _nnz(0),
        _lc(-1),
        _rc(-1),
        _rows()
    { } 

    KOKKOS_INLINE_FUNCTION
    CrsMatrixView(const CrsMatBaseType &b,
                  const ordinal_type offm, const ordinal_type m,
                  const ordinal_type offn, const ordinal_type n) 
      : _base(b),
        _offm(offm),
        _offn(offn),
        _m(m),
        _n(n),
        _nnz(0),
        _lc(-1),
        _rc(-1),
        _rows()
    { } 

    /// Destructor
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (o)
    KOKKOS_INLINE_FUNCTION
    ~CrsMatrixView() = default;

    /// Print out
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (x),
    /// - Callable in KokkosFunctors (x)
    std::ostream& showMe(std::ostream &os) const {
      std::streamsize prec = os.precision();
      os.precision(3);
      os << std::scientific;
      const int w = 4;
      if (_base.isValueArrayNull()) {
        os << "-- Base object is null --;";
      } else {
        const size_type nnz = (_rows.dimension_0() == _m ? _nnz : -1);
        os << _base.Label() << "::View, "
           << " Offs ( " << std::setw(w) << _offm << ", " << std::setw(w) << _offn << " ); "
           << " Dims ( " << std::setw(w) << _m    << ", " << std::setw(w) << _n    << " ); "
           << " NumNonZeros = " << std::setw(w) << nnz << ";"
           << " Density = " << std::setw(w) << (nnz == -1 ? double(0) : double(nnz)/_m/(_rc-_lc+1)) << ";";
      }
      os.unsetf(std::ios::scientific);
      os.precision(prec);
      return os;
    }

  };
}


#endif
