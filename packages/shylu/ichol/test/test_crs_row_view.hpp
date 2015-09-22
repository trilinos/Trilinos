#pragma once
#ifndef __TEST_CRS_ROW_VIEW_HPP__
#define __TEST_CRS_ROW_VIEW_HPP__

/// \file test_crs_matrix_view.hpp
/// \brief Test CrsMatrixBase class
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "tmg_crs_matrix_base_simple.hpp"

namespace Example {
  
  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testCrsRowView(const OrdinalType mbase, 
                     const OrdinalType nbase,
                     const OrdinalType offm,
                     const OrdinalType offn,
                     const OrdinalType m,
                     const OrdinalType n) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    typedef CrsRowView<CrsMatrixBaseType> CrsRowViewType;
    
    typedef Tmg_CrsMatrixBase_Simple<CrsMatrixBaseType> TmgType;

    __DOT_LINE__;
    cout << "testCrsRowView:: mbase = " << mbase << ", nbase = " << nbase 
         << ", offm = " << offm << ", offn = " << offn 
         << ", m = " << m << ", n = " << n << endl;
    __DOT_LINE__;

    TmgType tmg(mbase, nbase);
    
    int r_val = 0, r_val_prev = 0;
    string eval;

    { // Constructor with allocation
      r_val_prev = r_val;
      cout << "testCrsRowView::Begin - " << r_val << endl;

      const ordinal_type nnz_base = tmg.NumNonZeros();
      CrsMatrixBaseType AA("A, Allocated", mbase, nbase, nnz_base);

      r_val += tmg.fill(AA);
      r_val += tmg.check(AA);

      const auto epsilon = sqrt(NumericTraits<value_type>::epsilon());
      CrsMatrixViewType A(&AA, offm, m,
                          /**/ offn, n);
      A.fillRowViewArray();
      cout << A << endl;

      for (ordinal_type i=0;i<A.NumRows();++i) {
        CrsRowViewType a = A.RowView(i);
        cout << a << endl;
        for (ordinal_type j=0;j<a.NumNonZeros();++j) {
          const ordinal_type row_at_i = i + offm;
          const ordinal_type col_at_j = a.Col(j) + offn;
          const auto tmp = abs(tmg.Value(row_at_i, col_at_j) - a.Value(j)); 
          __ASSERT_TRUE__(tmp < epsilon);
        }
      }
      cout << "testCrsRowView::End - " << r_val << endl;      

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testCrsRowView::Eval - " << eval << endl;
    }
    
    __DOT_LINE__;
    
    return r_val;
  }
}

#endif
