#pragma once
#ifndef __ICHOL_RIGHT_BY_BLOCKS_SERIAL_HPP__
#define __ICHOL_RIGHT_BY_BLOCKS_SERIAL_HPP__

/// \file ichol_right_by_blocks_serial.hpp
/// \brief Sparse incomplete Cholesky factorization by blocks.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// This naively generates tasks without any merging of task blocks.

namespace Example { 

  using namespace std;
  
  template<typename CrsTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int genScalarTask_UpperRightByBlocks(typename CrsTaskViewType::policy_type &policy,
                                       const CrsTaskViewType A) {
    typedef typename CrsTaskViewType::value_type        value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;

    row_view_type a(A, 0); 
    value_type &aa = a.Value(0);

    int r_val = 0;
    IChol<Uplo::Upper,AlgoIChol::RightUnblockedOpt1>
      ::TaskFunctor<value_type>(aa).apply(r_val);

    return 0;
  }

  template<typename CrsTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int genTrsmTasks_UpperRightByBlocks(typename CrsTaskViewType::policy_type &policy,
                                     const CrsTaskViewType A,
                                     const CrsTaskViewType B) {
    typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
    typedef typename CrsTaskViewType::value_type        value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;

    row_view_type a(A,0), b(B,0); 
    value_type &aa = a.Value(0);

    int r_val = 0;
    for (ordinal_type j=0;j<b.NumNonZeros();++j) {
      value_type &bb = b.Value(j);

      Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,AlgoTrsm::ForRightBlocked>
        ::TaskFunctor<double,value_type>(Diag::NonUnit, 1.0, aa, bb).apply(r_val);
    }

    return 0;
  }

  template<typename CrsTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int genHerkTasks_UpperRightByBlocks(typename CrsTaskViewType::policy_type &policy,
                                      const CrsTaskViewType A,
                                      const CrsTaskViewType C) {
    typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
    typedef typename CrsTaskViewType::value_type        value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;

    // case that X.transpose, A.no_transpose, Y.no_transpose

    row_view_type a(A,0), c; 

    int r_val = 0;

    const ordinal_type nnz = a.NumNonZeros();

    // update herk on diagonals
    for (ordinal_type i=0;i<nnz;++i) {
      const ordinal_type row_at_i = a.Col(i);
      value_type &val_at_i = a.Value(i);

      c.setView(C, row_at_i);

      ordinal_type idx = c.Index(row_at_i, 0);
      if (idx >= 0) {
        value_type &cc = c.Value(idx);
        Herk<Uplo::Upper,Trans::ConjTranspose,AlgoHerk::ForRightBlocked>
          ::TaskFunctor<double,value_type>(-1.0, val_at_i, 1.0, cc).apply(r_val);
      }
    }

    // update gemm on off-diagonals 
    for (ordinal_type i=0;i<nnz;++i) {
      const ordinal_type row_at_i = a.Col(i);
      value_type &val_at_i = a.Value(i);

      c.setView(C, row_at_i);

      ordinal_type idx = 0;
      for (ordinal_type j=0;j<nnz && (idx > -2);++j) {
        const ordinal_type col_at_j = a.Col(j);
        value_type &val_at_j = a.Value(j);

        if (row_at_i != col_at_j) {
          idx = c.Index(col_at_j, idx);
          if (idx >= 0) {
            value_type &cc = c.Value(idx);
            Gemm<Trans::ConjTranspose,Trans::NoTranspose,AlgoGemm::ForRightBlocked>
              ::TaskFunctor<double,value_type>(-1.0, val_at_i, val_at_j, 1.0, cc).apply(r_val);
          }
        }
      }
    }
    
    return 0;
  }

}

#endif
