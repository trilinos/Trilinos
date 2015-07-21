#pragma once
#ifndef __TEST_DENSE_MATRIX_VIEW_HPP__
#define __TEST_DENSE_MATRIX_VIEW_HPP__

/// \file test_dense_matrix_view.hpp
/// \brief Test DenseMatrixBase class
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "tmg_dense_matrix_base_simple.hpp"

namespace Example {
  
  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testDenseMatrixView(const OrdinalType mbase, 
                          const OrdinalType nbase,
                          const OrdinalType offm,
                          const OrdinalType offn,
                          const OrdinalType m,
                          const OrdinalType n) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;
    typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;
    
    typedef Tmg_DenseMatrixBase_Simple<DenseMatrixBaseType> TmgType;

    __DOT_LINE__;
    cout << "testDenseMatrixView:: mbase = " << mbase << ", nbase = " << nbase 
         << ", offm = " << offm << ", offn = " << offn 
         << ", m = " << m << ", n = " << n << endl;
    __DOT_LINE__;
    
    TmgType tmg(mbase, nbase);
    
    int r_val = 0, r_val_prev = 0;
    string eval;

    { // Constructor 
      r_val_prev = r_val;
      cout << "testDenseMatrixView::Constructor::Begin - " << r_val << endl;

      DenseMatrixBaseType AA("A, Allocated", mbase, nbase);

      r_val += tmg.fill(AA);

      {
        DenseMatrixViewType A(&AA);
        
        __ASSERT_TRUE__(A.BaseObject() == &AA);
        __ASSERT_TRUE__(A.OffsetRows() == 0);
        __ASSERT_TRUE__(A.OffsetCols() == 0);
        __ASSERT_TRUE__(A.NumRows() == mbase);
        __ASSERT_TRUE__(A.NumCols() == nbase);

        r_val += tmg.check(A);
        A.showMeDetail(cout);
      }
      {
        DenseMatrixViewType A(&AA, offm, m,
                              /**/ offn, n);
        
        __ASSERT_TRUE__(A.BaseObject() == &AA);
        __ASSERT_TRUE__(A.OffsetRows() == offm);
        __ASSERT_TRUE__(A.OffsetCols() == offn);
        __ASSERT_TRUE__(A.NumRows() == m);
        __ASSERT_TRUE__(A.NumCols() == n);
        
        r_val += tmg.check(A);
        A.showMeDetail(cout);

        DenseMatrixViewType B(A);
        
        __ASSERT_TRUE__(A.BaseObject() == B.BaseObject());
        __ASSERT_TRUE__(A.OffsetRows() == B.OffsetRows());
        __ASSERT_TRUE__(A.OffsetCols() == B.OffsetCols());
        __ASSERT_TRUE__(A.NumRows() == B.NumRows());
        __ASSERT_TRUE__(A.NumCols() == B.NumCols());

        r_val += tmg.check(B);
        B.showMeDetail(cout);
      }
      
      cout << "testDenseMatrixView::Constructor::End - " << r_val << endl;      
      
      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testDenseMatrixView::Constructor::Eval - " << eval << endl;
    }
    
    __DOT_LINE__;
    
    return r_val;
  }
}

#endif
