#pragma once
#ifndef __ICHOL_RIGHT_BY_BLOCKS_VAR1_HPP__
#define __ICHOL_RIGHT_BY_BLOCKS_VAR1_HPP__

/// \file ichol_right_by_blocks_var1.hpp
/// \brief Sparse incomplete Cholesky factorization by blocks.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// This naively generates tasks without any merging of task blocks.

namespace Example { 

  using namespace std;
  
  template<typename ParallelForType,
           typename CrsTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int genScalarTask_UpperRightByBlocks(typename CrsTaskViewType::policy_type &policy,
                                       const CrsTaskViewType &A) {
    typedef typename CrsTaskViewType::value_type        value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;

    typedef typename CrsTaskViewType::future_type       future_type;
    typedef typename CrsTaskViewType::task_factory_type task_factory_type;

    row_view_type a(A, 0); 
    value_type &aa = a.Value(0);

    // construct a task
    future_type f = task_factory_type::create(policy,
                                              IChol<Uplo::Upper,AlgoIChol::RightUnblockedOpt1>
                                              ::TaskFunctor<ParallelForType,value_type>(aa));

    // manage dependence
    task_factory_type::addDependence(policy, f, aa.Future());
    aa.setFuture(f);

    // spawn a task
    task_factory_type::spawn(policy, f);

    return 0;
  }

  template<typename ParallelForType,
           typename CrsTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int genTrsmTasks_UpperRightByBlocks(typename CrsTaskViewType::policy_type &policy,
                                      const CrsTaskViewType &A,
                                      const CrsTaskViewType &B) {
    typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
    typedef typename CrsTaskViewType::value_type        value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;

    typedef typename CrsTaskViewType::future_type       future_type;
    typedef typename CrsTaskViewType::task_factory_type task_factory_type;

    row_view_type a(A,0), b(B,0); 
    value_type &aa = a.Value(0);

    const ordinal_type nnz = b.NumNonZeros();
    for (ordinal_type j=0;j<nnz;++j) {
      value_type &bb = b.Value(j);

      future_type f = task_factory_type
        ::create(policy, 
                 Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,AlgoTrsm::ForRightBlocked>
                 ::TaskFunctor<ParallelForType,double,value_type>(Diag::NonUnit, 1.0, aa, bb));
      
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
           typename CrsTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int genHerkTasks_UpperRightByBlocks(typename CrsTaskViewType::policy_type &policy,
                                      const CrsTaskViewType &A,
                                      const CrsTaskViewType &C) {
    typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
    typedef typename CrsTaskViewType::value_type        value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;

    typedef typename CrsTaskViewType::future_type       future_type;
    typedef typename CrsTaskViewType::task_factory_type task_factory_type;

    // case that X.transpose, A.no_transpose, Y.no_transpose

    row_view_type a(A,0), c; 

    const ordinal_type nnz = a.NumNonZeros();

    // update herk
    for (ordinal_type i=0;i<nnz;++i) {
      const ordinal_type row_at_i = a.Col(i);
      value_type &val_at_i = a.Value(i);

      c.setView(C, row_at_i);

      ordinal_type idx = 0;
      for (ordinal_type j=i;j<nnz && (idx > -2);++j) {
        const ordinal_type col_at_j = a.Col(j);
        value_type &val_at_j = a.Value(j);

        if (row_at_i == col_at_j) {
          idx = c.Index(row_at_i, idx);
          if (idx >= 0) {
            value_type &cc = c.Value(idx);
            future_type f = task_factory_type
              ::create(policy, 
                       Herk<Uplo::Upper,Trans::ConjTranspose,AlgoHerk::ForRightBlocked>
                       ::TaskFunctor<ParallelForType,double,value_type>(-1.0, val_at_i, 1.0, cc));
            
            // dependence
            task_factory_type::addDependence(policy, f, val_at_i.Future());              
            
            // self
            task_factory_type::addDependence(policy, f, cc.Future());
            
            // place task signature on y
            cc.setFuture(f);
            
            // spawn a task
            task_factory_type::spawn(policy, f);
          }
        } else {
          idx = c.Index(col_at_j, idx);
          if (idx >= 0) {
            value_type &cc = c.Value(idx);
            future_type f = task_factory_type
              ::create(policy, 
                       Gemm<Trans::ConjTranspose,Trans::NoTranspose,AlgoGemm::ForRightBlocked>
                       ::TaskFunctor<ParallelForType,double,value_type>(-1.0, val_at_i, val_at_j, 1.0, cc));
            
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
    }
    
    return 0;
  }

}

#endif
