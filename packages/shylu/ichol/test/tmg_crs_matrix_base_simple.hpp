#pragma once
#ifndef __TMG_CRS_MATRIX_BASE_SIMPLE__
#define __TMG_CRS_MATRIX_BASE_SIMPLE__

#include "crs_matrix_base.hpp"

/// \file tmg_crs_matrix_base.hpp
/// \brief Simple test matrix generation.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {
  
  using namespace std;

  template<typename CrsMatrixBaseType>
  class Tmg_CrsMatrixBase_Simple {
  public:
    typedef typename CrsMatrixBaseType::ordinal_type ordinal_type;
    typedef typename CrsMatrixBaseType::size_type size_type;
    typedef typename CrsMatrixBaseType::value_type value_type;

  private:
    ordinal_type _m, _n;
    
  public:
    Tmg_CrsMatrixBase_Simple(const ordinal_type m,
                             const ordinal_type n)
      : _m(m), _n(n) { }

    size_type NumNonZeros() const { return _m*_n; }

    value_type Value(const ordinal_type i, 
                     const ordinal_type j) {
      return value_type(i*_n + j);
    }

    int fill(CrsMatrixBaseType &A) {
      int r_val = 0;

      const size_type nnz = A.NumNonZeros();
      __ASSERT_TRUE__(nnz == NumNonZeros());
      
      typename CrsMatrixBaseType::size_type_array_ptr ptr = A.RowPtr();
      __ASSERT_TRUE__(ptr != NULL);
      
      for (ordinal_type i=0;i<(_m+1);++i)
        ptr[i] = _n*i;
      
      ordinal_type cnt = 0;
      for (ordinal_type i=0;i<_m;++i) {
        typename CrsMatrixBaseType::ordinal_type_array_ptr col = A.ColsInRow(i);
        typename CrsMatrixBaseType::value_type_array_ptr val = A.ValuesInRow(i);

        for (ordinal_type j=0;j<_n;++j,++cnt) {
          col[j] = j;
          val[j] = Value(i, j);
        }
      }
      return r_val;
    }

    int check(CrsMatrixBaseType &A) {
      int r_val = 0;
      const auto epsilon = sqrt(NumericTraits<value_type>::epsilon());
      for (ordinal_type i=0;i<_m;++i) {
        const ordinal_type j_nnz = A.NumNonZerosInRow(i);
        
        typename CrsMatrixBaseType::ordinal_type_array_ptr col = A.ColsInRow(i);
        typename CrsMatrixBaseType::value_type_array_ptr val = A.ValuesInRow(i);

        for (ordinal_type j=0;j<j_nnz;++j) {
          const auto tmp = abs(val[j] - Value(i, col[j]));
          __ASSERT_TRUE__(tmp < epsilon);          
        }
      }
      return r_val;
    }
  };

}

#endif
