#pragma once
#ifndef __TEST_CRS_MATRIX_BASE_HPP__
#define __TEST_CRS_MATRIX_BASE_HPP__

/// \file test_crs_matrix_base.hpp
/// \brief Test CrsMatrixBase class
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "crs_matrix_base.hpp"
#include "tmg_crs_matrix_base_simple.hpp"

namespace Example {
  
  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testCrsMatrixBase(const OrdinalType m, 
                        const OrdinalType n) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef Tmg_CrsMatrixBase_Simple<CrsMatrixBaseType> TmgType;

    __DOT_LINE__;
    cout << "testCrsMatrixBase:: m = " << m << ", n = " << n << endl;
    __DOT_LINE__;

    TmgType tmg(m, n);
    
    int r_val = 0, r_val_prev = 0;
    string eval;

    { // Default constructor
      r_val_prev = r_val;
      cout << "testCrsMatrixBase::DefaultConstructor::Begin - " << r_val << endl;

      CrsMatrixBaseType A;
      
      __ASSERT_TRUE__(A.NumRows() == 0);
      __ASSERT_TRUE__(A.NumCols() == 0);
      __ASSERT_TRUE__(A.NumNonZeros() == 0);

      //__ASSERT_TRUE__(A.RowPtr() == NULL);
      cout << "testCrsMatrixBase::DefaultConstructor::End - " << r_val << endl;

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testCrsMatrixBase::DefaultConstructor::Eval - " << eval << endl;
    }

    { // Constructor with allocation
      r_val_prev = r_val;
      cout << "testCrsMatrixBase::ConstructorInternalAllocation::Begin - " << r_val << endl;

      const ordinal_type nnz = tmg.NumNonZeros();
      CrsMatrixBaseType A("A, Allocated", m, n, nnz);

      r_val += tmg.fill(A);

      __ASSERT_TRUE__(A.NumRows() == m);
      __ASSERT_TRUE__(A.NumCols() == n);
      __ASSERT_TRUE__(A.NumNonZeros() == nnz);
      
      r_val += tmg.check(A); 

      cout << A << endl;
      cout << "testCrsMatrixBase::ConstructorInternalAllocation::End - " << r_val << endl;      

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testCrsMatrixBase::ConstructorInternalAllocation::Eval - " << eval << endl;
    }
    
    { // Constructor with external buffers
      r_val_prev = r_val;
      cout << "testCrsMatrixBase::ConstructorExternalAllocation::Begin - " << r_val << endl;
      
      const ordinal_type nnz = tmg.NumNonZeros();

      typename CrsMatrixBaseType::size_type_array    ap("External::RowPtrArray", m+1);
      typename CrsMatrixBaseType::ordinal_type_array aj("External::ColIdxArray", nnz);
      typename CrsMatrixBaseType::value_type_array   ax("External::ValueArray",  nnz);

      CrsMatrixBaseType A("A, Wrapped", 
                          m, n, nnz,
                          ap, aj, ax);

      r_val += tmg.fill(A);
    
      __ASSERT_TRUE__(A.NumRows() == m);
      __ASSERT_TRUE__(A.NumCols() == n);
      __ASSERT_TRUE__(A.NumNonZeros() == nnz);

      r_val += tmg.check(A);

      cout << "testCrsMatrixBase::ConstructorExternalAllocation::End - " << r_val << endl;

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testCrsMatrixBase::ConstructorExternalAllocation::Eval - " << eval << endl;
    }

    { // Copy constructor and copy operations
      r_val_prev = r_val;
      cout << "testCrsMatrixBase::Copy::Begin - " << r_val << endl;

      const ordinal_type nnz = tmg.NumNonZeros();
      CrsMatrixBaseType A("A, Source", m, n, nnz);

      r_val += tmg.fill(A);

      // ---------------------------------------------------------------------------
      cout << "testCrsMatrixBase::Copy::Constructor::Begin - " << r_val << endl;
      CrsMatrixBaseType B(A);
      B.setLabel("B, Shallow-copy of A");

      r_val += tmg.check(B);
      cout << "testCrsMatrixBase::Copy::Constructor::End - " << r_val << endl;
      // ---------------------------------------------------------------------------
      cout << "testCrsMatrixBase::Copy::DeepCopy::Begin - " << r_val << endl;
      CrsMatrixBaseType C("C, Deep-copy of A");
      C.copy(A);

      r_val += tmg.check(C);
      cout << "testCrsMatrixBase::Copy::DeepCopy::End - " << r_val << endl;
      // ---------------------------------------------------------------------------
      cout << "testCrsMatrixBase::Copy::DeepCopyLower::Begin - " << r_val << endl;
      
      CrsMatrixBaseType Dl("D, deep-copy of A lower triangular");
      Dl.copy(Uplo::Lower, A);

      r_val += tmg.check(Dl);
      cout << "testCrsMatrixBase::Copy::DeepCopyLower::End - " << r_val << endl;
      // ---------------------------------------------------------------------------
      cout << "testCrsMatrixBase::Copy::DeepCopyUpper::Begin - " << r_val << endl;

      CrsMatrixBaseType Du("D, deep-copy of A upper triangular");
      Du.copy(Uplo::Upper, A);

      r_val += tmg.check(Du);
      cout << "testCrsMatrixBase::Copy::DeepCopyUpper::End - " << r_val << endl;
      // ---------------------------------------------------------------------------
      cout << "testCrsMatrixBase::Copy::DeepCopyPermutation::Begin - " << r_val << endl;
      typename CrsMatrixBaseType::ordinal_type_array p ("External::PermVector", n);
      typename CrsMatrixBaseType::ordinal_type_array ip("External::InvPermVector", m);    

      for (ordinal_type i=0;i<m;++i)
        ip[i] = (m - i - 1);

      for (ordinal_type j=0;j<n;++j)
        p[j] = (n - j - 1);

      CrsMatrixBaseType E("E, permuted A");
      E.copy(p, ip, A);
      cout << E << endl;

      CrsMatrixBaseType F("F, permuted E");
      F.copy(ip, p, E);
      r_val += tmg.check(F);
      cout << F << endl;

      cout << "testCrsMatrixBase::Copy::DeepCopyPermutation::End - " << r_val << endl;
      // ---------------------------------------------------------------------------
      cout << "testCrsMatrixBase::Copy::End - " << r_val << endl;

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testCrsMatrixBase::Copy::Eval - " << eval << endl;
    }

    __DOT_LINE__;
    
    return r_val;
  }
}

#endif
