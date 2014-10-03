#pragma once
#ifndef __CRS_MATRIX_BASE_IMPL_HPP__
#define __CRS_MATRIX_BASE_IMPL_HPP__

namespace Example { 

  using namespace std;

  template<typename ValueType,
           typename OrdinalType, 
           typename SizeType>
  inline int 
  CrsMatrixBase<ValueType,OrdinalType,SizeType>
  ::importMatrixMarket(ifstream &file) {
    // skip initial title comments
    {
      ordinal_type m, n;
      size_type nnz;
          
      while (file.good()) {
        char c = file.peek();
        if (c == '%' || c == '\n') {
          file.ignore(256, '\n');
          continue;
        }
        break;
      }
          
      // read matrix specification
      file >> m >> n >> nnz;
          
      // construct workspace and set variables
      init(m, n, nnz);
    }

    // read the coordinate format (matrix-market)
    vector<ijv_type> mm(_nnz);
    {
      // matrix market use one base index
      const ordinal_type mm_base = 1; 

      for (size_type i=0;i<_nnz;++i) {
        ijv_type aij;
        file >> aij._i >> aij._j >> aij._val;
            
        // one base to zero base
        aij._i -= mm_base;
        aij._j -= mm_base;
            
        mm.push_back(aij);
      }
      sort(mm.begin(), mm.end(), less<ijv_type>());
    }

    // change mm to crs
    {
      ordinal_type ii = 0;
      size_type jj = 0;

      ijv_type prev = mm[0];
      _ap[ii++] = 0;
      _aj[jj] = prev._j;
      _ax[jj] = prev._val;
      ++jj;

      for (typename vector<ijv_type>::iterator it=(mm.begin()+1);it<mm.end();++it) {
        ijv_type aij = (*it);
        
        // row index
        if (aij._i != prev._i) {
          _ap[ii++] = jj; 
        }
            
        if (aij == prev) {
          --jj;
          _aj[jj]  = aij._j;
          _ax[jj] += aij._val;
        } else {
          _aj[jj] = aij._j;
          _ax[jj] = aij._val;
        }
        ++jj;
          
        prev = aij;
      }
          
      // add the last index to terminate the storage
      _ap[ii++] = jj;
      _nnz = jj;
    }
      
    return 0;
  }

  template<typename ValueType,
           typename OrdinalType, 
           typename SizeType>
  inline int 
  CrsMatrixBase<ValueType,OrdinalType,SizeType>
  ::showMe(ostream &os) const {
    streamsize prec = os.precision();
    os.precision(15);
    os << scientific;

    if (!_is_initialized) 
      os << " -- CrsMatrixBase is not initialized -- " << endl;

    os << " -- CrsMatrixBase -- " << endl
       << "    # of Rows      = " << _m << endl
       << "    # of Cols      = " << _n << endl
       << "    # of NonZeros  = " << _nnz << endl;
    const int w = 6;
    for (ordinal_type i=0;i<_m;++i) {
      size_type jbegin = _ap[i], jend = _ap[i+1];
      
      os << endl;
      for (size_type j=jbegin;j<jend;++j) 
        os << setw(w) << i << "  " 
           << setw(w) << _aj[j] << "  " 
           << _ax[j] << endl;
    }

    os.unsetf(ios::scientific);
    os.precision(prec);

    return 0;
  }

  template<typename ValueType,
           typename OrdinalType, 
           typename SizeType>
  inline int 
  CrsMatrixBase<ValueType,OrdinalType,SizeType>
  ::convertGraph(size_type &nnz,
                 size_type *rptr,
                 ordinal_type *cidx) const {
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

}


#endif
