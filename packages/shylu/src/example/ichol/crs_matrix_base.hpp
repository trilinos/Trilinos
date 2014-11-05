#pragma once
#ifndef __CRS_MATRIX_BASE_HPP__
#define __CRS_MATRIX_BASE_HPP__

/// \file crs_matrix_base.hpp
/// \brief CRS matrix base object interfaces to user provided input matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <Kokkos_Core.hpp>

#include "util.hpp"
#include "coo.hpp"

namespace Example { 

  using namespace std;
  
  /// \class CrsMatrixBase
  /// \breif CRS matrix base object using Kokkos view and subview
  template<typename ValueType,
           typename OrdinalType, 
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  class CrsMatrixBase : public Disp {
  public:
    typedef ValueType    value_type;
    typedef OrdinalType  ordinal_type;
    typedef SpaceType    space_type;
    typedef SizeType     size_type;
    typedef MemoryTraits memory_traits;

    // 1D view, layout does not matter; no template parameters for that
    typedef Kokkos::View<size_type*,   space_type,memory_traits> size_type_array;
    typedef Kokkos::View<ordinal_type*,space_type,memory_traits> ordinal_type_array;
    typedef Kokkos::View<value_type*,  space_type,memory_traits> value_type_array;

    typedef Kokkos::View<size_type*,   space_type,Kokkos::MemoryUnmanaged> size_type_array_view;
    typedef Kokkos::View<ordinal_type*,space_type,Kokkos::MemoryUnmanaged> ordinal_type_array_view;
    typedef Kokkos::View<value_type*,  space_type,Kokkos::MemoryUnmanaged> value_type_array_view;

    // range type
    template<typename T> using range_type = pair<T,T>;

    // external interface
    typedef Coo<CrsMatrixBase> ijv_type;
    
    friend class CrsMatrixHelper;

  private:
    string             _label;   //!< object label

    ordinal_type       _m;       //!< # of rows
    ordinal_type       _n;       //!< # of cols
    size_type          _nnz;     //!< # of nonzeros
    size_type_array    _ap;      //!< pointers to column index and values
    ordinal_type_array _aj;      //!< column index compressed format
    value_type_array   _ax;      //!< values

  protected:
    void createInternalArrays(const ordinal_type m, 
                              const ordinal_type n,
                              const size_type nnz) {
      _m = m;
      _n = n;
      _nnz = nnz;

      if (_ap.dimension_0() < m+1)
        _ap = size_type_array(_label+"::RowPtrArray", m+1);
      
      if (_aj.dimension_0() < nnz)
        _aj = ordinal_type_array(_label+"::ColsArray", nnz);
      
      if (_ax.dimension_0() < nnz)
        _ax = value_type_array(_label+"::ValuesArray", nnz);
    }

  public:

    KOKKOS_INLINE_FUNCTION
    void setLabel(string label) { _label = label; }

    KOKKOS_INLINE_FUNCTION
    string Label() const { return _label; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumRows() const { return _m; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumCols() const { return _n; }

    KOKKOS_INLINE_FUNCTION
    size_type NumNonZeros() const { return _nnz; }

    KOKKOS_INLINE_FUNCTION
    size_type
    RowPtr(const ordinal_type i) const { 
      return _ap[i];
    }
    
    KOKKOS_INLINE_FUNCTION
    ordinal_type_array_view
    ColsInRow(const ordinal_type i) const {
      return Kokkos::subview<ordinal_type_array_view>(_aj, range_type<size_type>(_ap[i],_ap[i+1]));
    }
    
    KOKKOS_INLINE_FUNCTION
    value_type_array_view   
    ValuesInRow(const ordinal_type i) const {
      return Kokkos::subview<value_type_array_view>(_ax, range_type<size_type>(_ap[i],_ap[i+1]));
    }

    KOKKOS_INLINE_FUNCTION
    ordinal_type  
    NumNonZerosInRow(const ordinal_type i) const { 
      return (_ap[i+1] - _ap[i]);
    } 

    KOKKOS_INLINE_FUNCTION
    value_type& Value(const ordinal_type k) { return _ax[k]; }

    KOKKOS_INLINE_FUNCTION
    value_type Value(const ordinal_type k) const { return _ax[k]; }

    /// \brief Default constructor.
    CrsMatrixBase() 
      : _label("CrsMatrixBase"),
        _m(0),
        _n(0),
        _nnz(0),
        _ap(),
        _aj(),
        _ax()
    { }

    /// \brief Constructor with label
    CrsMatrixBase(const string label) 
      : _label(label),
        _m(0),
        _n(0),
        _nnz(0),
        _ap(),
        _aj(),
        _ax()
    { }

    /// \brief Copy constructor (shallow copy), for deep-copy use a method copy
    template<typename VT,
             typename OT,
             typename ST,
             typename SpT,
             typename MT>
    CrsMatrixBase(const CrsMatrixBase<VT,OT,ST,SpT,MT> &b) 
      : _label(b._label),
        _m(b._m),
        _n(b._n),
        _nnz(b._nnz),
        _ap(b._ap), 
        _aj(b._aj),
        _ax(b._ax) 
    { }

    /// \brief Constructor to allocate internal data structures.
    CrsMatrixBase(const string label,
                  const ordinal_type m, 
                  const ordinal_type n, 
                  const ordinal_type nnz) 
      : _label(label), 
        _m(m),
        _n(n),
        _nnz(nnz),
        _ap(_label+"::RowPtrArray", m+1),
        _aj(_label+"::ColsArray", nnz),
        _ax(_label+"::ValuesArray", nnz)
    { }

    /// \brief Constructor to attach external arrays to the matrix.
    CrsMatrixBase(const string label,
                  const ordinal_type m, 
                  const ordinal_type n, 
                  const ordinal_type nnz,
                  const size_type_array &ap,
                  const ordinal_type_array &aj,
                  const value_type_array &ax) 
      : _label(label), 
        _m(m),
        _n(n),
        _nnz(nnz),
        _ap(ap), 
        _aj(aj),
        _ax(ax) 
    { }
    
    /// \brief deep copy of matrix b
    template<typename VT,
             typename OT,
             typename ST,
             typename SpT,
             typename MT>
    int copy(const CrsMatrixBase<VT,OT,ST,SpT,MT> &b) {
      createInternalArrays(b._m, b._n, b._nnz);

      Kokkos::deep_copy(_ap, b._ap);
      Kokkos::deep_copy(_aj, b._aj);
      Kokkos::deep_copy(_ax, b._ax);

      return 0;
    }

    /// \brief deep copy of lower/upper triangular of matrix b
    template<typename VT,
             typename OT,
             typename ST,
             typename SpT,
             typename MT>
    int 
    copy(const int uplo, 
         const CrsMatrixBase<VT,OT,ST,SpT,MT> &b) { 

      switch (uplo) {
      case Uplo::Lower: {
        // counting nnz for lower triangular
        _nnz = 0;
        for (ordinal_type i=0;i<b._m;++i) {
          ordinal_type jbegin = b._ap[i];
          ordinal_type jend   = b._ap[i+1];
          for (ordinal_type j=jbegin;j<jend;++j) 
            _nnz += (i >= b._aj[j]);
        }

        createInternalArrays(b._m, b._n, _nnz);
        
        // fill the matrix
        _nnz = 0;
        for (ordinal_type i=0;i<_m;++i) {
          ordinal_type jbegin = b._ap[i];
          ordinal_type jend   = b._ap[i+1];
          _ap[i] = _nnz;
          for (ordinal_type j=jbegin;j<jend;++j) {
            if (i >= b._aj[j]) {
              _aj[_nnz] = b._aj[j];
              _ax[_nnz] = b._ax[j]; 
              ++_nnz;
            }
          }
        }
        _ap[_m] = _nnz;
        break;
      }
      case Uplo::Upper:
        // not yet implemented
        break;
      }

      return 0;
    }

    /// \brief deep copy of matrix b with given permutation vectors
    template<typename VT,
             typename OT,
             typename ST,
             typename SpT,
             typename MT>
    int
    copy(const typename CrsMatrixBase<VT,OT,ST,SpT,MT>::ordinal_type_array &p,
         const typename CrsMatrixBase<VT,OT,ST,SpT,MT>::ordinal_type_array &ip,
         const CrsMatrixBase<VT,OT,ST,SpT,MT> &b) {

      createInternalArrays(b._m, b._n, b._nnz);

      // Question:: do I need to use Kokkos::vector ? 
      //            in other words, where do we permute matrix in factoriztion ?
      //            permuting a matrix is a kernel ? 
      vector<ijv_type> tmp;

      // any chance to use parallel_for ?
      _nnz = 0;
      for (ordinal_type i=0;i<_m;++i) {
        ordinal_type ii = ip[i];

        ordinal_type jbegin = b._ap[ii];
        ordinal_type jend   = b._ap[ii+1];

        _ap[i] = _nnz;
        for (ordinal_type j=jbegin;j<jend;++j) {
          ordinal_type jj = p[b._aj[j]];
          ijv_type aij(i, jj, b._ax[j]);
          tmp.push_back(aij);
        }
        // does std algorithm universally work on Kokkos any memory space ?
        sort(tmp.begin(), tmp.end(), less<ijv_type>());
        for (auto it=tmp.begin();it<tmp.end();++it) {
          ijv_type aij = (*it);

          _aj[_nnz] = aij.Col();
          _ax[_nnz] = aij.Val();
          ++_nnz;
        }
        tmp.clear();
      }
      _ap[_m] = _nnz;

      return 0;
    }

    ostream& showMe(ostream &os) const {
      os << " -- " << _label << " -- " << endl
         << "    # of Rows          = " << _m << endl
         << "    # of Cols          = " << _n << endl
         << "    # of NonZeros      = " << _nnz << endl
         << endl
         << "    RowPtrArray length = " << _ap.dimension_0() << endl
         << "    ColsArray   length = " << _aj.dimension_0() << endl 
         << "    ValuesArray length = " << _ax.dimension_0() << endl
         << endl;
      
      const int w = 15;
      if (_ap.size() && _aj.size() && _ax.size()) {
        os << setw(w) <<  "Row" << "  " 
           << setw(w) <<  "Col" << "  " 
           << setw(w) <<  "Val" << endl;
        for (ordinal_type i=0;i<_m;++i) {
          size_type jbegin = _ap[i], jend = _ap[i+1];
          for (size_type j=jbegin;j<jend;++j) {
            value_type val = _ax[j];
            os << setw(w) <<      i << "  " 
               << setw(w) << _aj[j] << "  " 
               << setw(w) <<    val << endl;
          }
        }
      }

      return os;
    }

    int convertGraph(size_type &nnz,
                     size_type_array rptr,
                     ordinal_type_array cidx) const {
      ordinal_type ii = 0;
      size_type jj = 0;
      
      for (ordinal_type i=0;i<_m;++i) {
        size_type jbegin = _ap[i], jend = _ap[i+1];
        rptr[ii++] = jj;
        for (size_type j=jbegin;j<jend;++j) 
          if (i != _aj[j]) 
            cidx[jj++] = _aj[j];
      }
      rptr[ii] = nnz = jj;
      
      return 0;
    }

    int importMatrixMarket(ifstream &file);
  };

}

#include "crs_matrix_base_import.hpp"

#endif
