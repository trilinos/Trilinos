#ifndef __TACHO_CRS_ROW_VIEW_HPP__
#define __TACHO_CRS_ROW_VIEW_HPP__

/// \file Tacho_CrsRowView.hpp
/// \brief A view to a row extracted from CrsMatrixView.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 

  /// \class CrsRowView
  /// \brief Extract row view partitioned by matrix view
  template<typename CrsMatBaseType>
  class CrsRowView {
  public:
    typedef typename CrsMatBaseType::space_type   space_type;

    typedef typename CrsMatBaseType::value_type   value_type;
    typedef typename CrsMatBaseType::ordinal_type ordinal_type;

    // range type
    template<typename T> using range_type = Kokkos::pair<T,T>;

    // 1D subview of the base matrix 
    typedef typename Kokkos::View<value_type*,  space_type,Kokkos::MemoryUnmanaged> value_type_array;
    typedef typename Kokkos::View<ordinal_type*,space_type,Kokkos::MemoryUnmanaged> ordinal_type_array; 
    
  private:
    // partitioned row info
    ordinal_type _offn, _n;    

    // this assumes a contiguous memory buffer
    ordinal_type_array _aj;
    value_type_array _ax;

  public:

    /// Interface functions
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (o)

    KOKKOS_INLINE_FUNCTION
    void setView(const CrsMatBaseType &A, 
                 const ordinal_type i) {
      _offn = 0;
      _n    = A.NumCols();
      _aj   = A.ColsInRow(i);
      _ax   = A.ValuesInRow(i);
    }

    KOKKOS_INLINE_FUNCTION
    void setView(const CrsMatrixView<CrsMatBaseType> &A, 
                 const ordinal_type i) {
      _offn = A.OffsetCols();
      _n    = A.NumCols();

      const ordinal_type ii = A.OffsetRows() + i;

      const auto cols = A.BaseObject().ColsInRow(ii);
      const auto vals = A.BaseObject().ValuesInRow(ii);
      
      ordinal_type begin = 0, end = cols.extent(0);

      begin = Util::getLowerBound(cols, begin, end, _offn);
      end   = Util::getLowerBound(cols, begin, end, _offn+_n);
      
      const auto range = range_type<ordinal_type>(begin, end);

      _aj = Kokkos::subview(cols, range);
      _ax = Kokkos::subview(vals, range);
    }

    KOKKOS_INLINE_FUNCTION
    ordinal_type OffsetCols() const { return _offn; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumCols() const { return _n; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumNonZeros() const { return _aj.extent(0); } 

    KOKKOS_INLINE_FUNCTION
    ordinal_type Col(const ordinal_type j) const { return _aj[j] - _offn; }

    KOKKOS_INLINE_FUNCTION
    value_type& Value(const ordinal_type j) { return _ax[j]; }

    KOKKOS_INLINE_FUNCTION
    value_type Value(const ordinal_type j) const { return _ax[j]; }
    
    KOKKOS_INLINE_FUNCTION
    ordinal_type Index(const ordinal_type col,
                       const ordinal_type prev = 0) const {
      const ordinal_type loc = _offn + col;

      // linear search from previous value
      ordinal_type j = prev;
      for (;j<_aj.extent(0) && _aj[j]<loc; ++j); 

      // if found, return index for the location, 
      // otherwise return -1 (not found), -2 (end of array)
      return (j < _aj.extent(0) ? (_aj[j] == loc ? j : -1) : -2);
    }

    // this operation is expensive and should be avoided.
    KOKKOS_INLINE_FUNCTION
    value_type ValueAtColumn(const ordinal_type col) const {
      const ordinal_type j = Index(col);
      return (j < 0 ? value_type(0) : _ax[j]);
    }

    /// ------------------------------------------------------------------

    /// Constructors
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors
    ///   - Default and copy constructors are allowed in KokkosFunctors.
    ///   - Creating internal workspace is not allowed in KokkosFunctors.
    KOKKOS_INLINE_FUNCTION
    CrsRowView()
      : _offn(0),
        _n(0),
        _aj(),
        _ax() 
    { }

    KOKKOS_INLINE_FUNCTION
    CrsRowView(const CrsRowView &b)
      : _offn(b._offn),
        _n(b._n),
        _aj(b._aj),
        _ax(b._ax) 
    { }

    KOKKOS_INLINE_FUNCTION
    CrsRowView(const CrsMatrixView<CrsMatBaseType> &A, 
               const ordinal_type i) {
      this->setView(A, i);
    }

    KOKKOS_INLINE_FUNCTION
    CrsRowView(const CrsMatBaseType &A, 
               const ordinal_type i) {
      this->setView(A, i);
    }

    /// Destructor
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (o)
    KOKKOS_INLINE_FUNCTION
    ~CrsRowView() = default;

    /// Print out
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (x),
    /// - Callable in KokkosFunctors (x)
    std::ostream& showMe(std::ostream &os) const {                                                
      const ordinal_type nnz = this->NumNonZeros();
      const ordinal_type offset = this->OffsetCols();
      os << "  offset = " << offset
         << ", nnz = " << nnz
         << std::endl; 
      const int w = 4;
      for (ordinal_type j=0;j<nnz;++j) {
        const auto val = _ax[j];
        os << "(" << std::setw(w) << Col(j) << ", "
           << std::setw(w) 
           << std::showpos << val << std::noshowpos << ")"
           << std::endl;
      }
      return os;
    }
  };
}

#endif
