#pragma once
#ifndef __TRI_SOLVE_U_CT_BY_BLOCKS_VAR1_HPP__
#define __TRI_SOLVE_U_CT_BY_BLOCKS_VAR1_HPP__

/// \file ichol_right_by_blocks_var1.hpp
/// \brief Sparse incomplete Cholesky factorization by blocks.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// This naively generates tasks without any merging of task blocks.

namespace Example { 

  using namespace std;

  template<typename ParallelForType,
           typename CrsTaskViewTypeA,
           typename DenseTaskViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int genTrsmTasks_UpperByBlocks(typename CrsTaskViewType::policy_type &policy,
                                 const CrsTaskViewType &A,
                                 const DenseTaskViewType &B) {
    typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
    typedef typename CrsTaskViewType::value_type        crs_value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;

    typedef typename DenseTaskViewType::value_type      dense_value_type;

    typedef typename CrsTaskViewType::future_type       future_type;
    typedef typename CrsTaskViewType::task_factory_type task_factory_type;

    row_view_type a(A,0);
    crs_value_type &aa = a.Value(0);
    
    for (ordinal_type j=0;j<B.NumCols();++j) {
      dense_value_type &bb = b.Value(0, j);

      future_type f = task_factory_type
        ::create(policy, 
                 Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,AlgoTrsm::ForTriSolveBlocked>
                 ::TaskFunctor<ParallelForType,double,
                 crs_value_type,dense_value_type>(Diag::NonUnit, 1.0, aa, bb));
      
      // trsm dependence
      task_factory_type::addDependence(policy, f, aa.Future());

      // self
      task_factory_type::addDependence(policy, f, bb.Future());

      // place task signature on b
      bb.setFuture(f);
      
      // spawn a task
      task_factory_type::spawn(policy, f);              
    }

    return 0;
  }

  template<typename ParallelForType,
           typename CrsTaskViewType,
           typename DenseTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int genGemmTasks_UpperByBlocks(typename CrsTaskViewType::policy_type &policy,
                                 const CrsTaskViewType &A,
                                 const DenseTaskViewType &B,
                                 const DenseTaskViewType &C) {
    typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
    typedef typename CrsTaskViewType::value_type        crs_value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;

    typedef typename DenseTaskViewType::value_type      dense_value_type;

    typedef typename CrsTaskViewType::future_type       future_type;
    typedef typename CrsTaskViewType::task_factory_type task_factory_type;

    // A ct, B nt, C

    row_view_type a(A,0); 

    const ordinal_type nnz = a.NumNonZeros();

    for (ordinal_type i=0;i<nnz;++i) {
      const ordinal_type row_at_i = a.Col(i);
      crs_value_type &val_at_i = a.Value(i);

      for (ordinal_type j=0;j<C.NumCols();++j) {
        const ordinal_type col_at_j = j;
        dense_value_type &val_at_j = B.Value(row_at_i, col_at_j);

        dense_value_type &cc = C.Value(0, col_at_j);

        future_type f = task_factory_type
          ::create(policy, 
                   Gemm<Trans::ConjTranspose,Trans::NoTranspose,AlgoGemm::ForTriSolveBlocked>
                     ::TaskFunctor<ParallelForType,double,
                     crs_value_type,dense_value_type,dense_value_type>(-1.0, val_at_i, val_at_j, 1.0, cc));
          
          // dependence
          task_factory_type::addDependence(policy, f, val_at_i.Future());
          task_factory_type::addDependence(policy, f, val_at_j.Future());
          
          // self
          task_factory_type::addDependence(policy, f, cc.Future());
          
          // place task signature on y
          cc.setFuture(f);
          
          // spawn a task
          task_factory_type::spawn(policy, f);
        }
      }
    }
    
    return 0;
  }
  



}

#endif
