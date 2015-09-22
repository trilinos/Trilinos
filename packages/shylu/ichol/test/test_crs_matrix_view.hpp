#pragma once
#ifndef __TEST_CRS_MATRIX_VIEW_HPP__
#define __TEST_CRS_MATRIX_VIEW_HPP__

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
  int testCrsMatrixView(const OrdinalType mbase, 
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
    
    typedef Tmg_CrsMatrixBase_Simple<CrsMatrixBaseType> TmgType;

    __DOT_LINE__;
    cout << "testCrsMatrixView:: mbase = " << mbase << ", nbase = " << nbase 
         << ", offm = " << offm << ", offn = " << offn 
         << ", m = " << m << ", n = " << n << endl;
    __DOT_LINE__;

    TmgType tmg(mbase, nbase);
    
    int r_val = 0, r_val_prev = 0;
    string eval;

    { // Constructor 
      r_val_prev = r_val;
      cout << "testCrsMatrixView::Constructor::Begin - " << r_val << endl;

      const ordinal_type nnz_base = tmg.NumNonZeros();
      CrsMatrixBaseType AA("A, Allocated", mbase, nbase, nnz_base);

      r_val += tmg.fill(AA);

      {
        CrsMatrixViewType A(&AA);
        
        __ASSERT_TRUE__(A.BaseObject() == &AA);
        __ASSERT_TRUE__(A.OffsetRows() == 0);
        __ASSERT_TRUE__(A.OffsetCols() == 0);
        __ASSERT_TRUE__(A.NumRows() == mbase);
        __ASSERT_TRUE__(A.NumCols() == nbase);

        cout << A << endl;
      }
      {
        CrsMatrixViewType A(&AA, offm, m,
                            /**/ offn, n);
        
        __ASSERT_TRUE__(A.BaseObject() == &AA);
        __ASSERT_TRUE__(A.OffsetRows() == offm);
        __ASSERT_TRUE__(A.OffsetCols() == offn);
        __ASSERT_TRUE__(A.NumRows() == m);
        __ASSERT_TRUE__(A.NumCols() == n);

        cout << A << endl;

        CrsMatrixViewType B(A);
        
        __ASSERT_TRUE__(A.BaseObject() == B.BaseObject());
        __ASSERT_TRUE__(A.OffsetRows() == B.OffsetRows());
        __ASSERT_TRUE__(A.OffsetCols() == B.OffsetCols());
        __ASSERT_TRUE__(A.NumRows() == B.NumRows());
        __ASSERT_TRUE__(A.NumCols() == B.NumCols());

        cout << B << endl;
      }

      cout << "testCrsMatrixView::Constructor::End - " << r_val << endl;      

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testCrsMatrixView::Constructor::Eval - " << eval << endl;
    }
    
    __DOT_LINE__;
    
    return r_val;
  }
}

#endif
