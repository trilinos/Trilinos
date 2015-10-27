#pragma once
#ifndef __GEMM_NT_NT_DENSE_BY_BLOCKS_HPP__
#define __GEMM_NT_NT_DENSE_BY_BLOCKS_HPP__

/// \file gemm_nt_nt_dense_by_blocks.hpp
/// \brief Dense matrix-matrix multiplication 
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  using namespace std;

  // Gemm-By-Blocks
  // ==============
  template<>
  template<typename ScalarType,
           typename DenseTaskViewTypeA,
           typename DenseTaskViewTypeB,
           typename DenseTaskViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Gemm<Trans::NoTranspose,Trans::NoTranspose,
       AlgoGemm::DenseMatrixByBlocks>
  ::invoke(typename DenseTaskViewTypeA::policy_type &policy,
           const typename DenseTaskViewTypeA::policy_type::member_type &member,
           const ScalarType alpha,
           DenseTaskViewTypeA &A,
           DenseTaskViewTypeB &B,
           const ScalarType beta,
           DenseTaskViewTypeC &C) {
    typedef typename DenseTaskViewTypeA::ordinal_type ordinal_type;
    typedef typename DenseTaskViewTypeA::value_type   value_type;

    typedef typename CrsTaskViewType::future_type       future_type;
    typedef typename CrsTaskViewType::task_factory_type task_factory_type;

    if (member.team_rank() == 0) {
      for (ordinal_type p=0;p<A.NumCols();++p) {
        const ScalarType beta_select = (p > 0 ? 1.0 : beta);
        for (ordinal_type k2=0;k2<C.NumCols();++k2) {
          value_type &bb = B.Value( p, k2);
          for (ordinal_type k1=0;k1<C.NumRows();++k1) {
            value_type &aa = A.Value(k1,  p);
            value_type &cc = C.Value(k1, k2);

            future_type f = task_factory_type::create(policy,
                                                      Gemm<Trans::NoTranspose,Trans::NoTranspose,
                                                      AlgoGemm::ExternalBlas>
                                                      ::TaskFunctor<double,value_type,value_type,value_type>
                                                      (alpha, aa, bb, beta_select, cc));

            // dependence
            task_factory_type::addDependence(policy, f, aa.Future());
            task_factory_type::addDependence(policy, f, bb.Future());

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
