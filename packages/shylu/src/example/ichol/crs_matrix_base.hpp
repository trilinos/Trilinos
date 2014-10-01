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
    ordinal_type  _base;    // base index 0 or 1
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
        _base = 0;
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
    int init(const ordinal_type base,
             const ordinal_type m,
             const ordinal_type n,
             const size_type    nnz) {
      if (_is_initialized) 
        finalize();

      _base = base;
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
    ordinal_type  BaseVal()     const { return _base; }
    ordinal_type  NumRows()     const { return _m;    }
    ordinal_type  NumCols()     const { return _n;    }
    size_type     NumNonZeros() const { return _nnz;  }
    ordinal_type  NumNonZerosInRow(const int i) { return (_ap[i+1] - _ap[i]); } 
    size_type*    RowPtr(const int i=0)   const { return &_ap[i]; }
    ordinal_type* ColIndex(const int i=0) const { return &_aj[_ap[i]]; }
    value_type*   Value(const int i=0)    const { return &_ax[_ap[i]]; }

    CrsMatrixBase() 
      : _base(0),
        _m(0),
        _nnz(0),
        _ap(NULL),
        _aj(NULL),
        _ax(NULL),
        _is_initialized(false),
        _is_symmetry(false) { }

    CrsMatrixBase(const CrsMatrixBase &b) {
      init(b._base, b._m, b._n, b._nnz);

      copy(b._ap, b._ap+_m+1, _ap);
      copy(b._aj, b._aj+_nnz, _aj); 
      copy(b._ax, b._ax+_nnz, _ax); 

      _is_symmetry = b._is_symmetry;
    }

    virtual~CrsMatrixBase() {
      finalize();
    }

    int importMatrixMarket(ifstream &file);
    int showMe(ostream &os, ordinal_type *perm = NULL) const;
    int convertGraph(size_type &nnz,
                     size_type *rptr,
                     ordinal_type *cidx) const;
  };
}

#include "crs_matrix_base_impl.hpp"

#endif
