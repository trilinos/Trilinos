#ifndef __TACHO_CRS_MATRIX_BASE_HPP__
#define __TACHO_CRS_MATRIX_BASE_HPP__

/// \file Tacho_CrsMatrixBase.hpp
/// \brief CRS matrix base object interfaces to user provided input matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho { 

  /// \class CrsMatrixBase
  /// \breif CRS matrix base object using Kokkos view and subview
  template<typename ValueType,
           typename OrdinalType, 
           typename SizeType = OrdinalType,
           typename SpaceType = void>
  class CrsMatrixBase {
  public:
    typedef ValueType    value_type;
    typedef OrdinalType  ordinal_type;
    typedef SpaceType    space_type;
    typedef SizeType     size_type;

    // 1D view, layout does not matter; no template parameters for that
    typedef Kokkos::View<value_type*,  space_type> value_type_array;
    typedef Kokkos::View<ordinal_type*,space_type> ordinal_type_array;
    typedef Kokkos::View<size_type*,   space_type> size_type_array;

    // range type
    template<typename T> using range_type = Kokkos::pair<T,T>;

    // external interface
    typedef Coo<ordinal_type,value_type> ijv_type;
    
    template<typename, typename, typename, typename>
    friend class CrsMatrixBase;
    
    //friend class CrsMatrixTools;

  private:
    char               _label[Util::LabelSize];   //!< object label

    ordinal_type       _m;          //!< # of rows
    ordinal_type       _n;          //!< # of cols

    size_type          _nnz;        //!< # of nonzeros

    //size_type_array    _ap;   //!< pointers to column index and values

    // if I have more time,
    size_type_array    _ap_begin, _ap_end;     //!< pointers to column index and values

    ordinal_type_array _aj;         //!< column index compressed format
    value_type_array   _ax;         //!< values

  protected:

    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (x)
    KOKKOS_INLINE_FUNCTION
    void createInternalArrays(const ordinal_type m, 
                              const ordinal_type n,
                              const size_type nnz) {
      _m = m;
      _n = n;
      _nnz = nnz;

      if (static_cast<ordinal_type>(_ap_begin.extent(0)) < m) {
        _ap_begin = size_type_array("CrsMatrixBase::RowPtrBeginArray", m);
      } else {
        // otherwise initialize it
        Kokkos::Impl::ViewFill<size_type_array>(_ap_begin, size_type());
      }

      if (static_cast<ordinal_type>(_ap_end.extent(0)) < m) {
        _ap_end = size_type_array("CrsMatrixBase::RowPtrEndArray", m);
      } else {
        // otherwise initialize it
        Kokkos::Impl::ViewFill<size_type_array>(_ap_end, size_type());
      }

      if (static_cast<size_type>(_aj.extent(0)) < nnz) {
        _aj = ordinal_type_array("CrsMatrixBase::ColsArray", nnz);
      } else {
        // otherwise initialize it
        Kokkos::Impl::ViewFill<ordinal_type_array>(_aj, ordinal_type());
      }

      if (static_cast<size_type>(_ax.extent(0)) < nnz) {
        _ax = value_type_array("CrsMatrixBase::ValuesArray", nnz);
      } else {
        // otherwise initialize it
        Kokkos::Impl::ViewFill<value_type_array>(_ax, value_type());
      }
    }
    
  public:

    /// Interface functions
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (o)
    /// \brief Constructor to attach external arrays to the matrix.
    KOKKOS_INLINE_FUNCTION    
    void setExternalMatrix(const ordinal_type m, 
                           const ordinal_type n, 
                           const ordinal_type nnz,
                           const size_type_array &ap_begin,
                           const size_type_array &ap_end,
                           const ordinal_type_array &aj,
                           const value_type_array &ax) {
      _m = m;
      _n = n;
      _nnz = nnz;
      _ap_begin = ap_begin;
      _ap_end = ap_end;
      _aj = aj;
      _ax = ax;
    }

    KOKKOS_INLINE_FUNCTION
    bool isValueArrayNull() const {
      return !_ax.extent(0);
    }

    KOKKOS_INLINE_FUNCTION
    void setLabel(const char *label) { 
      strncpy(_label, label, Util::min(strlen(label)+1, Util::LabelSize));
    }

    KOKKOS_INLINE_FUNCTION
    void setNumNonZeros() { 
      // the nnz here is the size of value and column arrays
      if (_m) _nnz = _ap_end(_m-1);
    }

    KOKKOS_INLINE_FUNCTION
    const char* Label() const { return _label; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumRows() const { return _m; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumCols() const { return _n; }

    KOKKOS_INLINE_FUNCTION
    size_type NumNonZeros() const { return _nnz; }

    KOKKOS_INLINE_FUNCTION
    size_type_array RowPtrBegin() { return _ap_begin; }

    KOKKOS_INLINE_FUNCTION
    size_type_array RowPtrEnd() { return _ap_end; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type_array Cols() const { return _aj; }

    KOKKOS_INLINE_FUNCTION
    value_type_array Values() const { return _ax; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumNonZerosInRow(const ordinal_type i) const { return (_ap_end[i] - _ap_begin[i]); } 

    KOKKOS_INLINE_FUNCTION
    ordinal_type_array ColsInRow(const ordinal_type i) const { 
      return Kokkos::subview(_aj, range_type<size_type>(_ap_begin[i], _ap_end[i]));
    }
    
    KOKKOS_INLINE_FUNCTION
    value_type_array ValuesInRow(const ordinal_type i) const { 
      return Kokkos::subview(_ax, range_type<size_type>(_ap_begin[i], _ap_end[i]));
    }

    KOKKOS_INLINE_FUNCTION
    size_type& RowPtrBegin(const ordinal_type i) { return _ap_begin[i]; }

    KOKKOS_INLINE_FUNCTION
    size_type RowPtrBegin(const ordinal_type i) const { return _ap_begin[i]; }

    KOKKOS_INLINE_FUNCTION
    size_type& RowPtrEnd(const ordinal_type i) { return _ap_end[i]; }
    
    KOKKOS_INLINE_FUNCTION
    size_type RowPtrEnd(const ordinal_type i) const { return _ap_end[i]; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type& Col(const ordinal_type k) { return _aj[k]; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type Col(const ordinal_type k) const { return _aj[k]; }

    KOKKOS_INLINE_FUNCTION
    value_type& Value(const ordinal_type k) { return _ax[k]; }

    KOKKOS_INLINE_FUNCTION
    value_type Value(const ordinal_type k) const { return _ax[k]; }

    /// ------------------------------------------------------------------

    /// Constructors
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors
    ///   - Default and copy constructors are allowed in KokkosFunctors.
    ///   - Creating internal workspace is not allowed in KokkosFunctors.

    /// \brief Default constructor.
    KOKKOS_INLINE_FUNCTION
    CrsMatrixBase() 
      : _m(0),
        _n(0),
        _nnz(0),
        _ap_begin(),
        _ap_end(),
        _aj(),
        _ax()
    { 
      setLabel("CrsMatrixBase");
    }

    /// \brief Constructor with label
    KOKKOS_INLINE_FUNCTION
    CrsMatrixBase(const char *label) 
      : _m(0),
        _n(0),
        _nnz(0),
        _ap_begin(),
        _ap_end(),
        _aj(),
        _ax()
    { 
      setLabel(label); 
    }

    /// \brief Constructor with label
    KOKKOS_INLINE_FUNCTION
    CrsMatrixBase(const char *label,
                  const ordinal_type m,
                  const ordinal_type n,
                  const size_type nnz) 
      : _m(0),
        _n(0),
        _nnz(0),
        _ap_begin(),
        _ap_end(),
        _aj(),
        _ax()
    { 
      setLabel(label); 
      createInternalArrays(m, n, nnz);
    }

    /// \brief Copy constructor (shallow copy), for deep-copy use a method copy
    template<typename VT,
             typename OT,
             typename ST,
             typename SpT>
    KOKKOS_INLINE_FUNCTION    
    CrsMatrixBase(const CrsMatrixBase<VT,OT,ST,SpT> &b) 
      : _m(b._m),
        _n(b._n),
        _nnz(b._nnz),
        _ap_begin(b._ap_begin), 
        _ap_end(b._ap_end), 
        _aj(b._aj),
        _ax(b._ax) 
    { 
      setLabel(b._label); 
    }
    
    /// \brief Constructor to attach external arrays to the matrix.
    KOKKOS_INLINE_FUNCTION    
    CrsMatrixBase(const char *label,
                  const ordinal_type m, 
                  const ordinal_type n, 
                  const ordinal_type nnz,
                  const size_type_array &ap_begin,
                  const size_type_array &ap_end,
                  const ordinal_type_array &aj,
                  const value_type_array &ax) 
      : _m(m),
        _n(n),
        _nnz(nnz),
        _ap_begin(ap_begin), 
        _ap_end(ap_end), 
        _aj(aj),
        _ax(ax) 
    { 
      setLabel(label);
    }

    /// ------------------------------------------------------------------

    /// Destructor
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (x)
    KOKKOS_INLINE_FUNCTION
    ~CrsMatrixBase() = default;

    /// Create and mirror
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (x)
    KOKKOS_INLINE_FUNCTION
    void
    clear() {
      _m = 0; _n = 0; _nnz = 0;

      _ap_begin = size_type_array();
      _ap_end   = size_type_array();
      _aj       = ordinal_type_array();
      _ax       = value_type_array();
    }

    KOKKOS_INLINE_FUNCTION
    void
    create(const ordinal_type m,
           const ordinal_type n,
           const size_type nnz) {
      createInternalArrays(m, n, nnz);
    }

    template<typename SpT>
    KOKKOS_INLINE_FUNCTION
    void
    createConfTo(const CrsMatrixBase<value_type,ordinal_type,size_type,SpT> &b) {
      createInternalArrays(b._m, b._n, b._nnz);
    }
    
    /// \brief deep copy of matrix b
    /// Callable: Device (o), KokkosFunctors (x), Blocking (o)
    template<typename SpT>
    KOKKOS_INLINE_FUNCTION
    void
    mirror(const CrsMatrixBase<value_type,ordinal_type,size_type,SpT> &b) {
      if (Kokkos::Impl::is_same<SpT,space_type>::value) {
        // when the space is same, everything is shallow copy
        // setLabel(b._label);
        _m        = b._m;
        _n        = b._n;
        _nnz      = b._nnz;
        _ap_begin = b._ap_begin;
        _ap_end   = b._ap_end;
        _aj       = b._aj;
        _ax       = b._ax;
      } else {
        // when the space is different, perform deep copy
        createInternalArrays(b._m, b._n, b._nnz);

        const auto ap_begin_range = range_type<ordinal_type>(0, Util::min(_ap_begin.extent(0), b._ap_begin.extent(0)));
        const auto ap_end_range   = range_type<ordinal_type>(0, Util::min(_ap_end.extent(0),   b._ap_end.extent(0)));
        const auto aj_range = range_type<size_type>   (0, Util::min(_aj.extent(0), b._aj.extent(0)));
        const auto ax_range = range_type<size_type>   (0, Util::min(_ax.extent(0), b._ax.extent(0)));

        space_type::execution_space::fence();
        Kokkos::deep_copy(Kokkos::subview(_ap_begin, ap_begin_range), Kokkos::subview(b._ap_begin, ap_begin_range));
        Kokkos::deep_copy(Kokkos::subview(_ap_end,   ap_end_range),   Kokkos::subview(b._ap_end,   ap_end_range));
        Kokkos::deep_copy(Kokkos::subview(_aj, aj_range), Kokkos::subview(b._aj, aj_range));
        Kokkos::deep_copy(Kokkos::subview(_ax, ax_range), Kokkos::subview(b._ax, ax_range));
        space_type::execution_space::fence();
      }
    }
    /// ------------------------------------------------------------------

    /// Print out
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (x),
    /// - Callable in KokkosFunctors (x)

    /// \brief print out to stream
    std::ostream& showMe(std::ostream &os) const {
      std::streamsize prec = os.precision();
      os.precision(16);
      os << std::scientific;

      os << " -- " << _label << " -- " << std::endl
         << "    # of Rows          = " << _m << std::endl
         << "    # of Cols          = " << _n << std::endl
         << "    # of NonZeros      = " << _nnz << std::endl
         << "    Ap End             = " << _ap_end(_m-1) << std::endl
         << std::endl
         << "    RowPtrBeginArray length = " << _ap_begin.extent(0) << std::endl
         << "    RowPtrEndArray length   = " << _ap_end.extent(0) << std::endl
         << "    ColArray length         = " << _aj.extent(0) << std::endl 
         << "    ValueArray length       = " << _ax.extent(0) << std::endl
         << std::endl;
      
      const int w = 10;
      if (_ap_begin.size() && _ap_end.size() && _aj.size() && _ax.size()) {
        os << std::setw(w) <<  "Row" << "  " 
           << std::setw(w) <<  "Col" << "  " 
           << std::setw(w) <<  "Val" << std::endl;
        for (ordinal_type i=0;i<_m;++i) {
          const size_type jbegin = _ap_begin[i], jend = _ap_end[i];
          for (size_type j=jbegin;j<jend;++j) {
            value_type val = _ax[j];
            os << std::setw(w) <<      i << "  " 
               << std::setw(w) << _aj[j] << "  " 
               << std::setw(w) 
               << std::showpos <<    val << std::noshowpos
               << std::endl;
          }
        }
      }

      os.unsetf(std::ios::scientific);
      os.precision(prec);

      return os;
    }

  };

}

#endif

    // /// \brief add the matrix b into this non-zero entires
    // template<typename VT,
    //          typename OT,
    //          typename ST,
    //          typename SpT,
    //          typename MT>
    // int 
    // add(const CrsMatrixBase<VT,OT,ST,SpT,MT> &b) { 

    //   const ordinal_type m = min(b._m, _m);
    //   for (ordinal_type i=0;i<m;++i) {
    //     const size_type jaend = _ap[i+1];
    //     const size_type jbend = b._ap[i+1];

    //     size_type ja = _ap[i];
    //     size_type jb = b._ap[i];
        
    //     for ( ;jb<jbend;++jb) {
    //       for ( ;(_aj[ja]<b._aj[jb] && ja<jaend);++ja);
    //       _ax[ja] += (_aj[ja] == b._aj[jb])*b._ax[jb];
    //     }
    //   }

    //   return 0;
    // }

    // int symmetrize(const int uplo, 
    //                const bool conjugate = false) {
    //   vector<ijv_type> mm;
    //   mm.reserve(_nnz*2);

    //   for (ordinal_type i=0;i<_m;++i) {
    //     const size_type jbegin = _ap[i];
    //     const size_type jend   = _ap[i+1];
    //     for (size_type jj=jbegin;jj<jend;++jj) {
    //       const ordinal_type j = _aj[jj];
    //       const value_type val = (conjugate ? conj(_ax[j]) : _ax[j]);
    //       if        (uplo == Uplo::Lower && i > j) {
    //         mm.push_back(ijv_type(i, j, val));
    //         mm.push_back(ijv_type(j, i, val));
    //       } else if (uplo == Uplo::Upper && i < j) {
    //         mm.push_back(ijv_type(i, j, val));
    //         mm.push_back(ijv_type(j, i, val));
    //       } else if (i == j) {
    //         mm.push_back(ijv_type(i, i, val));
    //       }
    //     }
    //   }
    //   sort(mm.begin(), mm.end(), less<ijv_type>());

    //   createInternalArrays(_m, _n, mm.size());
      
    //   ijv2crs(mm);
      
    //   return 0;
    // }

    // int hermitianize(int uplo) {
    //   return symmetrize(uplo, true);
    // }
