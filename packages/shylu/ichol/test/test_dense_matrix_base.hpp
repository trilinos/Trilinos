#pragma once
#ifndef __TEST_DENSE_MATRIX_BASE_HPP__
#define __TEST_DENSE_MATRIX_BASE_HPP__

/// \file test_dense_matrix_base.hpp
/// \brief Test DenseMatrixBase class
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "dense_matrix_base.hpp"
#include "tmg_dense_matrix_base_simple.hpp"

namespace Example {
  
  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testDenseMatrixBase(const OrdinalType m, 
                        const OrdinalType n) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;
    typedef Tmg_DenseMatrixBase_Simple<DenseMatrixBaseType> TmgType;

    __DOT_LINE__;
    cout << "testDenseMatrixBase:: m = " << m << ", n = " << n << endl;
    __DOT_LINE__;

    TmgType tmg(m, n);
    
    int r_val = 0, r_val_prev = 0;
    string eval;

    { // Default constructor
      r_val_prev = r_val;
      cout << "testDenseMatrixBase::DefaultConstructor::Begin - " << r_val << endl;

      DenseMatrixBaseType A;
      
      __ASSERT_TRUE__(A.NumRows() == 0);
      __ASSERT_TRUE__(A.NumCols() == 0);
      __ASSERT_TRUE__(A.ColStride() == 0);
      __ASSERT_TRUE__(A.RowStride() == 0);

      cout << "testDenseMatrixBase::DefaultConstructor::End - " << r_val << endl;

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testDenseMatrixBase::DefaultConstructor::Eval - " << eval << endl;
    }

    { // Constructor with allocation
      r_val_prev = r_val;
      cout << "testDenseMatrixBase::ConstructorInternalAllocation::Begin - " << r_val << endl;

      DenseMatrixBaseType A("A, Allocated", m, n);

      r_val += tmg.fill(A);

      __ASSERT_TRUE__(A.NumRows() == m);
      __ASSERT_TRUE__(A.NumCols() == n);
      __ASSERT_TRUE__(A.ColStride() == m);
      __ASSERT_TRUE__(A.RowStride() == 1);
      
      r_val += tmg.check(A); 

      cout << A << endl;
      cout << "testDenseMatrixBase::ConstructorInternalAllocation::End - " << r_val << endl;      

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testDenseMatrixBase::ConstructorInternalAllocation::Eval - " << eval << endl;
    }

    { // Constructor with external buffers
      r_val_prev = r_val;
      cout << "testDenseMatrixBase::ConstructorExternalAllocation::ColumMajor::Begin - " << r_val << endl;
      
      typename DenseMatrixBaseType::value_type_array   ax("External::ValueArray", m*n);

      DenseMatrixBaseType A("A, Wrapped", 
                            m, n, -1, -1,
                            ax);

      r_val += tmg.fill(A);
    
      __ASSERT_TRUE__(A.NumRows() == m);
      __ASSERT_TRUE__(A.NumCols() == n);
      __ASSERT_TRUE__(A.ColStride() == m);
      __ASSERT_TRUE__(A.RowStride() == 1);

      r_val += tmg.check(A);

      cout << "testDenseMatrixBase::ConstructorExternalAllocation::ColumnMajor::End - " << r_val << endl;

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testDenseMatrixBase::ConstructorExternalAllocation::ColumnMajor::Eval - " << eval << endl;
    }
    { // Constructor with external buffers
      r_val_prev = r_val;
      cout << "testDenseMatrixBase::ConstructorExternalAllocation::RowMajor::Begin - " << r_val << endl;
      
      typename DenseMatrixBaseType::value_type_array ax("External::ValueArray", m*n);

      DenseMatrixBaseType A("A, Wrapped", 
                            m, n, 1, n,
                            ax);

      r_val += tmg.fill(A);
    
      __ASSERT_TRUE__(A.NumRows() == m);
      __ASSERT_TRUE__(A.NumCols() == n);
      __ASSERT_TRUE__(A.ColStride() == 1);
      __ASSERT_TRUE__(A.RowStride() == n);

      r_val += tmg.check(A);

      cout << "testDenseMatrixBase::ConstructorExternalAllocation::RowMajor::End - " << r_val << endl;

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testDenseMatrixBase::ConstructorExternalAllocation::RowMajor::Eval - " << eval << endl;
    }
    { // Constructor with external buffers
      r_val_prev = r_val;
      cout << "testDenseMatrixBase::ConstructorExternalAllocation::IrregularStride::Begin - " << r_val << endl;
      
      typename DenseMatrixBaseType::value_type_array ax("External::ValueArray", 2*m*n);

      DenseMatrixBaseType A("A, Wrapped", 
                            m, n, 2*m, 2,
                            ax);

      r_val += tmg.fill(A);
    
      __ASSERT_TRUE__(A.NumRows() == m);
      __ASSERT_TRUE__(A.NumCols() == n);
      __ASSERT_TRUE__(A.ColStride() == 2*m);
      __ASSERT_TRUE__(A.RowStride() == 2);

      r_val += tmg.check(A);

      cout << "testDenseMatrixBase::ConstructorExternalAllocation::IrregularStrideEnd - " << r_val << endl;

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testDenseMatrixBase::ConstructorExternalAllocation::IrregularStride::Eval - " << eval << endl;
    }

    { // Copy constructor and copy operations
      r_val_prev = r_val;
      cout << "testDenseMatrixBase::Copy::Begin - " << r_val << endl;

      DenseMatrixBaseType A("A, Source", m, n);

      r_val += tmg.fill(A);

      // ---------------------------------------------------------------------------
      cout << "testDenseMatrixBase::Copy::Constructor::Begin - " << r_val << endl;
      DenseMatrixBaseType B(A);
      B.setLabel("B, Shallow-copy of A");

      r_val += tmg.check(B);
      cout << "testDenseMatrixBase::Copy::Constructor::End - " << r_val << endl;
      // ---------------------------------------------------------------------------
      cout << "testDenseMatrixBase::Copy::DeepCopy::Begin - " << r_val << endl;
      DenseMatrixBaseType C("C, Deep-copy of A");
      C.copy(A);

      r_val += tmg.check(C);
      cout << "testDenseMatrixBase::Copy::DeepCopy::End - " << r_val << endl;
      // ---------------------------------------------------------------------------
      cout << "testDenseMatrixBase::Copy::DeepCopyPermutation::Begin - " << r_val << endl;
      typename DenseMatrixBaseType::ordinal_type_array p ("External::PermVector", n);
      typename DenseMatrixBaseType::ordinal_type_array ip("External::InvPermVector", m);    

      for (ordinal_type i=0;i<m;++i)
        ip[i] = (m - i - 1);

      for (ordinal_type j=0;j<n;++j)
        p[j] = (n - j - 1);

      DenseMatrixBaseType E("E, permuted A");
      E.copy(p, A);
      cout << E << endl;

      DenseMatrixBaseType F("F, permuted E");
      F.copy(ip, E);
      r_val += tmg.check(F);
      cout << F << endl;

      cout << "testDenseMatrixBase::Copy::DeepCopyPermutation::End - " << r_val << endl;
      // ---------------------------------------------------------------------------
      cout << "testDenseMatrixBase::Copy::End - " << r_val << endl;

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testDenseMatrixBase::Copy::Eval - " << eval << endl;
    }

    __DOT_LINE__;
    
    return r_val;
  }
}

#endif
