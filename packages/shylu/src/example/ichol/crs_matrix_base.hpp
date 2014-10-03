#pragma once
#ifndef __CRS_MATRIX_BASE_HPP__
#define __CRS_MATRIX_BASE_HPP__

#include "util.hpp"
#include "coo.hpp"

namespace Example { 

  using namespace std;
  
  template<typename ValueType,
           typename OrdinalType, 
           typename SizeType>
  class CrsMatrixBase {
  public:
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef Coo<value_type,ordinal_type> ijv_type;
    
  private:
    ordinal_type  _m;       // # of rows
    ordinal_type  _n;       // # of cols
    size_type     _nnz;     // # of nonzeros
    size_type    *_ap;      // pointer to column index and values
    ordinal_type *_aj;      // column index compressed format
    value_type   *_ax;      // values

    bool _is_initialized;
    bool _is_symmetry;

    int finalize() {
      if (_is_initialized) {
        _m    = 0;
        _n    = 0;
        _nnz  = 0;
          
        _ap = NULL;
        _aj = NULL;
        _ax = NULL;
          
        _is_initialized = false;
        _is_symmetry = false;
      }
      return 0;
    }
    int init(const ordinal_type m,
             const ordinal_type n,
             const size_type    nnz) {
      if (_is_initialized) 
        finalize();

      _m    = m;
      _n    = n;
      _nnz  = nnz;
        
      _ap = new size_type[m+1]();
      _aj = new ordinal_type[nnz]();
      _ax = new value_type[nnz]();
        
      _is_initialized = true;
      _is_symmetry = false;
      
      return 0;
    }
      
  public:
    ordinal_type  NumRows() const { return _m; }
    ordinal_type  NumCols() const { return _n; }

    size_type*    RowPtr(const ordinal_type i=0)   const { return &_ap[i]; }
    ordinal_type* ColIndex(const ordinal_type i=0) const { return &_aj[_ap[i]]; }
    value_type*   Value(const ordinal_type i=0)    const { return &_ax[_ap[i]]; }

    size_type     NumNonZeros() const { return _nnz;  }
    ordinal_type  NumNonZerosInRow(const ordinal_type i) const { return (_ap[i+1] - _ap[i]); } 

    CrsMatrixBase() 
      : _m(0),
        _n(0),
        _nnz(0),
        _ap(NULL),
        _aj(NULL),
        _ax(NULL),
        _is_initialized(false),
        _is_symmetry(false) { }

    CrsMatrixBase(CrsMatrixBase &b) {
      init(b._m, b._n, b._nnz);

      copy(b._ap, b._ap+_m+1, _ap);
      copy(b._aj, b._aj+_nnz, _aj); 
      copy(b._ax, b._ax+_nnz, _ax); 

      _is_symmetry = b._is_symmetry;
    }

    CrsMatrixBase(CrsMatrixBase &b, const int uplo) {
      _is_symmetry = b._is_symmetry;
      switch (uplo) {
      case Uplo::Lower: {
        init(b._m, b._n, b._nnz);

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

      // maybe resize
    }

    virtual~CrsMatrixBase() {
      finalize();
    }

    int importMatrixMarket(ifstream &file);
    int showMe(ostream &os) const;
    int convertGraph(size_type &nnz,
                     size_type *rptr,
                     ordinal_type *cidx) const;
  };
}

#include "crs_matrix_base_impl.hpp"

#endif
