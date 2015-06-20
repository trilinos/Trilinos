#pragma once
#ifndef __TMG_DENSE_MATRIX_BASE_SIMPLE__
#define __TMG_DENSE_MATRIX_BASE_SIMPLE__

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

/// \file tmg_dense_matrix_base.hpp
/// \brief Simple test matrix generation.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {
  
  using namespace std;

  template<typename DenseMatrixBaseType>
  class Tmg_DenseMatrixBase_Simple {
  public:
    typedef typename DenseMatrixBaseType::ordinal_type ordinal_type;
    typedef typename DenseMatrixBaseType::size_type size_type;
    typedef typename DenseMatrixBaseType::value_type value_type;

  private:
    ordinal_type _m, _n;
    
  public:
    Tmg_DenseMatrixBase_Simple(const ordinal_type m,
                               const ordinal_type n)
      : _m(m), _n(n) { }

    value_type Value(const ordinal_type i, 
                     const ordinal_type j) {
      return value_type(i*_n + j);
    }

    int fill(DenseMatrixBaseType &A) {
      int r_val = 0;
      
      for (ordinal_type i=0;i<_m;++i) 
        for (ordinal_type j=0;j<_n;++j)
          A.Value(i, j) = Value(i, j);
      
      return r_val;
    }

    int check(DenseMatrixBaseType &A) {
      int r_val = 0;
      const auto epsilon = sqrt(NumericTraits<value_type>::epsilon());
      for (ordinal_type i=0;i<_m;++i) 
        for (ordinal_type j=0;j<_n;++j) {
          const auto tmp = abs(A.Value(i, j) - Value(i, j));
          __ASSERT_TRUE__(tmp < epsilon);          
        }
      return r_val;
    }

    int check(DenseMatrixView<DenseMatrixBaseType> &A) {
      int r_val = 0;
      const auto epsilon = sqrt(NumericTraits<value_type>::epsilon());
      for (ordinal_type i=0;i<A.NumRows();++i) 
        for (ordinal_type j=0;j<A.NumCols();++j) {
          const ordinal_type row_at_i = A.OffsetRows() + i;
          const ordinal_type col_at_j = A.OffsetCols() + j;

          const auto tmp = abs(A.Value(i, j) - Value(row_at_i, col_at_j));
          __ASSERT_TRUE__(tmp < epsilon);          
        }
      return r_val;
    }
  };

}

#endif
