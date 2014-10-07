#pragma once
#ifndef __CRS_MATRIX_VIEW_HPP__
#define __CRS_MATRIX_VIEW_HPP__

namespace Example { 

  using namespace std;

  // forward declaration
  template <typename ValueType, 
            typename OrdinalType>
  class CrsRowView;

  template<typename CrsMatBaseType>
  class CrsMatrixView : public Disp {
  public:
    typedef typename CrsMatBaseType::value_type   value_type;
    typedef typename CrsMatBaseType::ordinal_type ordinal_type;
    typedef typename CrsMatBaseType::size_type    size_type;

  private:
    ordinal_type  _offm;    // offset in rows
    ordinal_type  _offn;    // offset in cols
    ordinal_type  _m;       // # of rows
    ordinal_type  _n;       // # of cols

    CrsMatBaseType *_base;   // pointer to the base object
    
  public:
    CrsMatBaseType* BaseObject() const { return _base; }
    
    void setView(CrsMatBaseType *base,
                 ordinal_type offm, ordinal_type m,
                 ordinal_type offn, ordinal_type n) {
      _base = base;
      _offm = offm; _m = m;
      _offn = offn; _n = n;
    }
    
    ordinal_type  OffsetRows() const { return _offm; }
    ordinal_type  OffsetCols() const { return _offn; }
    
    ordinal_type  NumRows() const { return _m; }
    ordinal_type  NumCols() const { return _n; }
    
    CrsRowView<value_type,ordinal_type> extractRow(const int i) const { 
      ordinal_type ii = _offm + i;  // i at base
      
      // grep a row in base
      ordinal_type *ci_base = _base->ColIndex(ii); 
      value_type   *val_base = _base->Value(ii); 
      ordinal_type nnz_base_row = _base->NumNonZerosInRow(ii);

      // count
      ordinal_type *ci_view = NULL;
      value_type   *val_view = NULL;
      ordinal_type nnz_view_row = 0;

      ordinal_type nn = (_offn + _n);
      for (ordinal_type k=0;k<nnz_base_row;++k) {
        if (ci_base[k] >= nn)
          break;

        if (ci_base[k] >= _offn)
          ++nnz_view_row;

        if (nnz_view_row == 1) {
          ci_view = &ci_base[k];
          val_view = &val_base[k];
        }
      }

      // global view
      return CrsRowView<value_type,ordinal_type>(_offn, _n, nnz_view_row, ci_view, val_view);
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
      streamsize prec = os.precision();
      os.precision(15);
      os << scientific;
      
      os << endl
         << " -- CrsMatrixView -- " << endl
         << "    Offset in Rows = " << _offm << endl
         << "    Offset in Cols = " << _offn << endl
         << "    # of Rows      = " << _m << endl
         << "    # of Cols      = " << _n << endl;

      const int w = 6;
      for (ordinal_type i=0;i<_m;++i) {
        os << endl;
        CrsRowView<value_type,ordinal_type> row = extractRow(i);
        row.showMe(os);
      }

      os.unsetf(ios::scientific);
      os.precision(prec);
      
      return os;
    }

  };
}

#endif
