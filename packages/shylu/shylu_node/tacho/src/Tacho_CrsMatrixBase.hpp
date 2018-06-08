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
             typename SpaceType>
    class CrsMatrixBase {
    public:
      typedef ValueType    value_type;
      typedef SpaceType    space_type;

      typedef Kokkos::View<value_type*,  space_type> value_type_array;
      typedef Kokkos::View<ordinal_type*,space_type> ordinal_type_array;
      typedef Kokkos::View<size_type*,   space_type> size_type_array;
    
      template<typename, typename>
      friend class CrsMatrixBase;
    
    private:
      ordinal_type       _m;          //!< # of rows
      ordinal_type       _n;          //!< # of cols
      size_type          _nnz;        //!< # of nonzeros

      size_type_array    _ap;         //!< pointers to column index and values
      ordinal_type_array _aj;         //!< column index compressed format
      value_type_array   _ax;         //!< values

    protected:

      void createInternal(const ordinal_type m, 
                          const ordinal_type n,
                          const size_type nnz) {
        _m = m;
        _n = n;
        _nnz = nnz;

        if (static_cast<ordinal_type>(_ap.extent(0)) < (m+1)) 
          _ap = size_type_array("CrsMatrixBase::RowPtrArray", m+1);
        else 
          Kokkos::Impl::ViewFill<size_type_array>(_ap, size_type());

        if (static_cast<size_type>(_aj.extent(0)) < nnz) 
          _aj = ordinal_type_array("CrsMatrixBase::ColsArray", nnz);
        else 
          Kokkos::Impl::ViewFill<ordinal_type_array>(_aj, ordinal_type());

        if (static_cast<size_type>(_ax.extent(0)) < nnz) 
          _ax = value_type_array("CrsMatrixBase::ValuesArray", nnz);
        else 
          Kokkos::Impl::ViewFill<value_type_array>(_ax, value_type());
      }
    
    public:

      /// Interface functions
      /// ------------------------------------------------------------------
      void setExternalMatrix(const ordinal_type m, 
                             const ordinal_type n, 
                             const size_type nnz,
                             const size_type_array &ap,
                             const ordinal_type_array &aj,
                             const value_type_array &ax) {
        // TODO:: check each array's space
        _m = m;
        _n = n;
        _nnz = nnz;
        _ap = ap;
        _aj = aj;
        _ax = ax;
      }

      void setExternalMatrix(const ordinal_type m, 
                             const ordinal_type n, 
                             const ordinal_type nnz,
                             const size_type* ap,
                             const ordinal_type* aj,
                             const value_type* ax) {
        _m = m;
        _n = n;
        _nnz = nnz;
        _ap = size_type_array(ap, _m+1);
        _aj = ordinal_type_array(aj, _nnz);
        _ax = value_type_array(ax, _nnz);
      }

      KOKKOS_INLINE_FUNCTION
      void setNumNonZeros() { 
        if (_m) _nnz = _ap(_m);
      }

      KOKKOS_INLINE_FUNCTION
      size_type_array& RowPtr() { return _ap; }

      KOKKOS_INLINE_FUNCTION
      ordinal_type_array& Cols() { return _aj; }

      KOKKOS_INLINE_FUNCTION
      value_type_array& Values() { return _ax; }
    
      KOKKOS_INLINE_FUNCTION
      ordinal_type NumRows() const { return _m; }

      KOKKOS_INLINE_FUNCTION
      ordinal_type NumCols() const { return _n; }

      KOKKOS_INLINE_FUNCTION
      size_type NumNonZeros() const { return _nnz; }

      KOKKOS_INLINE_FUNCTION
      size_type& RowPtrBegin(const ordinal_type i) { return _ap[i]; }

      KOKKOS_INLINE_FUNCTION
      size_type RowPtrBegin(const ordinal_type i) const { return _ap[i]; }

      KOKKOS_INLINE_FUNCTION
      size_type& RowPtrEnd(const ordinal_type i) { return _ap[i+1]; }
    
      KOKKOS_INLINE_FUNCTION
      size_type RowPtrEnd(const ordinal_type i) const { return _ap[i+1]; }

      KOKKOS_INLINE_FUNCTION
      ordinal_type& Col(const ordinal_type k) { return _aj[k]; }

      KOKKOS_INLINE_FUNCTION
      ordinal_type Col(const ordinal_type k) const { return _aj[k]; }

      KOKKOS_INLINE_FUNCTION
      value_type& Value(const ordinal_type k) { return _ax[k]; }

      KOKKOS_INLINE_FUNCTION
      value_type Value(const ordinal_type k) const { return _ax[k]; }


      /// Constructors
      /// ------------------------------------------------------------------

      /// \brief Default constructor.
      CrsMatrixBase() 
        : _m(0),
          _n(0),
          _nnz(0),
          _ap(),
          _aj(),
          _ax()
      { 
      }

      /// \brief Constructor with label
      CrsMatrixBase(const ordinal_type m,
                    const ordinal_type n,
                    const size_type nnz) 
        : _m(0),
          _n(0),
          _nnz(0),
          _ap(),
          _aj(),
          _ax()
      { 
        createInternal(m, n, nnz);
      }

      /// \brief Copy constructor (shallow copy), for deep-copy use a method copy
      CrsMatrixBase(const CrsMatrixBase &b) 
        : _m(b._m),
          _n(b._n),
          _nnz(b._nnz),
          _ap(b._ap), 
          _aj(b._aj),
          _ax(b._ax) 
      { 
      }
    
      KOKKOS_INLINE_FUNCTION
      ~CrsMatrixBase() = default;


      /// Create 
      /// ------------------------------------------------------------------

      void
      clear() {
        _m = 0; _n = 0; _nnz = 0;

        _ap = size_type_array();
        _aj = ordinal_type_array();
        _ax = value_type_array();
      }

      void
      create(const ordinal_type m,
             const ordinal_type n,
             const size_type nnz) {
        createInternal(m, n, nnz);
      }

      template<typename SpT>
      void
      createConfTo(const CrsMatrixBase<value_type,SpT> &b) {
        createInternal(b._m, b._n, b._nnz);
      }

    
      /// Create 
      /// ------------------------------------------------------------------

      /// \brief deep copy of matrix b
      template<typename SpT>
      inline
      void
      createMirror(const CrsMatrixBase<value_type,SpT> &b) {
        _m        = b._m;
        _n        = b._n;
        _nnz      = b._nnz;
      
        _ap = Kokkos::create_mirror_view(typename space_type::memory_space(), b._ap);
        _aj = Kokkos::create_mirror_view(typename space_type::memory_space(), b._aj);
        _ax = Kokkos::create_mirror_view(typename space_type::memory_space(), b._ax);
      }

      /// \brief deep copy of matrix b
      template<typename SpT>
      inline
      void
      copy(const CrsMatrixBase<value_type,SpT> &b) {
        Kokkos::deep_copy(_ap, b._ap);
        Kokkos::deep_copy(_aj, b._aj);
        Kokkos::deep_copy(_ax, b._ax);
      }

      /// \brief print out to stream
      inline
      std::ostream& showMe(std::ostream &os, const bool detail = false) const {
        std::streamsize prec = os.precision();
        os.precision(16);
        os << std::scientific;

        os << " -- CrsMatrixBase -- " << std::endl
           << "    # of Rows          = " << _m << std::endl
           << "    # of Cols          = " << _n << std::endl
           << "    # of NonZeros      = " << _nnz << std::endl
           << std::endl
           << "    RowPtrArray length = " << _ap.extent(0) << std::endl
           << "    ColArray length    = " << _aj.extent(0) << std::endl 
           << "    ValueArray length  = " << _ax.extent(0) << std::endl
           << std::endl;
      
        if (detail) {
          const int w = 10;
          if ( (_ap.size() >  _m  ) && 
               (_aj.size() >= _nnz) &&
               (_ax.size() >= _nnz) ) {
            os << std::setw(w) <<  "Row" << "  " 
               << std::setw(w) <<  "Col" << "  " 
               << std::setw(w) <<  "Val" << std::endl;
            for (ordinal_type i=0;i<_m;++i) {
              const size_type jbegin = _ap[i], jend = _ap[i+1];
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
        }

        os.unsetf(std::ios::scientific);
        os.precision(prec);

        return os;
      }

    };

} // namespace tacho
#endif
