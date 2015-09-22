#pragma once
#ifndef __ICHOL_LEFT_BY_BLOCKS_LOWER_VAR1_HPP__
#define __ICHOL_LEFT_BY_BLOCKS_LOWER_VAR1_HPP__

/// \file ichol_left_by_blocks_lower_var1.hpp
/// \brief Sparse incomplete Cholesky factorization by blocks.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// This naively generates tasks without any merging of task blocks.

namespace Example { 

  using namespace std;
  
  template<>
  template<typename CrsTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int 
  ICholLeftByBlocks<Uplo::Lower>
  ::genScalarTask(typename CrsTaskViewType::policy_type &policy,
                  const CrsTaskViewType A) {
    typedef typename CrsTaskViewType::value_type        value_type;
    typedef typename CrsTaskViewType::future_type       future_type;
    typedef typename CrsTaskViewType::task_factory_type task_factory_type;
    
    // extract the block matrix
    value_type &aa = A.extractRow(0).Value(0);

    // construct a task
    future_type f = task_factory_type
      ::create(policy,
               IChol<Uplo::Lower,Algo::LeftUnblocked>
               ::TaskFunctor<value_type>(aa));
    
    // manage dependence
    task_factory_type::addDependence(policy, aa.Future(), f);
    aa.setFuture(f);

    // spawn a task
    task_factory_type::spawn(policy, f);

    return 0;
  }

  template<>
  template<typename CrsTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int 
  ICholLeftByBlocks<Uplo::Lower>
  ::genGemmTasks(typename CrsTaskViewType::policy_type &policy,
                 const CrsTaskViewType A,
                 const CrsTaskViewType X,
                 const CrsTaskViewType Y) {
    typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
    typedef typename CrsTaskViewType::value_type        value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;
    typedef typename CrsTaskViewType::future_type       future_type;
    typedef typename CrsTaskViewType::task_factory_type task_factory_type;
    
    // case that X.transpose, A.no_transpose, Y.no_transpose

    row_view_type x = X.extractRow(0);

    for (ordinal_type i=0;i<Y.NumRows();++i) {
      row_view_type y = Y.extractRow(i);
      row_view_type a = A.extractRow(i);

      ordinal_type idx_y = y.Index(0);
      
      if (idx_y >= 0 && a.NumNonZeros()) {
        value_type &yy = y.Value(idx_y);

        for (ordinal_type j=0;j<x.NumNonZeros();++j) {
          ordinal_type idx_a = a.Index(x.Col(j));

          if (idx_a >= 0) {
            value_type &aa = a.Value(idx_a);
            value_type &xx = x.Value(j);

            future_type f = task_factory_type
              ::create(policy, 
                       Gemm<Trans::NoTranspose,Trans::Transpose>
                       ::TaskFunctor<double,value_type>(-1.0, aa, xx, 1.0, yy));
            
            // gemm dependence
            task_factory_type::addDependence(policy, aa.Future(), f);              
            
            // self
            task_factory_type::addDependence(policy, yy.Future(), f);
            
            // place task signature on y
            yy.setFuture(f);

            // spawn a task
            task_factory_type::spawn(policy, f);
          }
        }
      }
    }

    return 0;
  }


  template<>
  template<typename CrsTaskViewType>
  KOKKOS_INLINE_FUNCTION
  int 
  ICholLeftByBlocks<Uplo::Lower>
  ::genTrsmTasks(typename CrsTaskViewType::policy_type &policy,
                 const CrsTaskViewType A,
                 const CrsTaskViewType B) {
    typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
    typedef typename CrsTaskViewType::value_type        value_type;
    typedef typename CrsTaskViewType::row_view_type     row_view_type;
    typedef typename CrsTaskViewType::future_type       future_type;
    typedef typename CrsTaskViewType::task_factory_type task_factory_type;

    value_type &aa = A.extractRow(0).Value(0);

    for (ordinal_type i=0;i<B.NumRows();++i) {
      row_view_type b = B.extractRow(i);
      ordinal_type idx_b = b.Index(0);

      if (idx_b >= 0) {
        value_type &bb = b.Value(idx_b);

        future_type f = task_factory_type
          ::create(policy, 
                   Trsm<Side::Right,Uplo::Lower,Trans::Transpose>
                   ::TaskFunctor<double,value_type>(Diag::NonUnit, 1.0, aa, bb));
        
        // trsm dependence
        task_factory_type::addDependence(policy, aa.Future(), f);

        // self
        task_factory_type::addDependence(policy, bb.Future(), f);

        // place task signature on b
        bb.setFuture(f);

        // spawn a task
        task_factory_type::spawn(policy, f);              
      }
    }

    return 0;
  }

}

#endif
