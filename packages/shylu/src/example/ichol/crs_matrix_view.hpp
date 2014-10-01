#pragma once
#ifndef __CRS_MATRIX_VIEW_HPP__
#define __CRS_MATRIX_VIEW_HPP__

namespace Example { 

  using namespace std;

  template<typename CrsMatrixBase>
  class CrsMatrixView {
  public:
    typedef CrsMatrixBase::value_type   value_type;
    typedef CrsMatrixBase::ordinal_type ordinal_type;
    typedef CrsMatrixBase::size_type    size_type;

  private:
    typedef valute_type vt;
    typedef ordinal_type ot;
    typedef size_type st;

    ordinal_type  _offm;    // offset in rows
    ordinal_type  _offn;    // offset in cols
    ordinal_type  _m;       // # of rows
    ordinal_type  _n;       // # of cols

    CrsMatrixBase &_base;   // pointer to the base object
    
  public:
    ordinal_type  BaseVal()     const { return _base.BaseVal(); }
    ordinal_type  NumRows()     const { return _m; }
    ordinal_type  NumNonZerosInRow(const int i) { return (_ap[i+1] - _ap[i]); } 
    size_type*    RowPtr(const int i=0)   const { return &_ap[i]; }
    ordinal_type* ColIndex(const int i=0) const { return &_aj[_ap[i]]; }
    value_type*   Value(const int i=0)    const { return &_ax[_ap[i]]; }
    
    CrsRowView<vt,ot,st> ExtractRow(const int i) { 
      return CrsRowView<vt,ot,st>(NumNonZerosInRow(i), 
                                  RowPtr(i), 
                                  ColIndex(i), 
                                  Value(i));
    }

    CrsMatrixBase(const CrsMatrixBase &b) 
      : _base(b),
        _offm(0),
        _offn(0),
        _m(b.NumRows()),
        _n(b.NumCols()) 
    { } 

    CrsMatrixBase(const CrsMatrixBase &b,
                  const ordinal_type offm, const ordinal_type m,
                  const ordinal_type offn, const ordinal_type n) 
      : _base(b),
        _offm(offm),
        _offn(offn),
        _m(m),
        _n(n) 
    { } 
    virtual~CrsMatrixBase() {
    }

  };
}

#endif
