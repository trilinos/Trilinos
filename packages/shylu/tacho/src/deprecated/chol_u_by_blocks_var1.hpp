#pragma once
#ifndef __CHOL_U_BY_BLOCKS_VAR1_HPP__
#define __CHOL_U_BY_BLOCKS_VAR1_HPP__

/// \file chol_u_by_blocks_var1.hpp
/// \brief Sparse incomplete Cholesky factorization by blocks.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// This naively generates tasks without any merging of task blocks.

namespace Tacho { 

  using namespace std;

  template<template<int> class ControlType,
           typename ParallelForType,
           typename CrsTaskViewType>
  class CholUpperByBlocksVar1 {
  public:
    KOKKOS_INLINE_FUNCTION
    static int genScalarTask(typename CrsTaskViewType::policy_type &policy,
                             CrsTaskViewType &A) {
      typedef typename CrsTaskViewType::value_type        value_type;
      typedef typename CrsTaskViewType::row_view_type     row_view_type;
      
      typedef typename CrsTaskViewType::future_type       future_type;
      typedef typename CrsTaskViewType::task_factory_type task_factory_type;
      
      row_view_type a(A, 0); 
      value_type &aa = a.Value(0);
      
      // construct a task
      future_type f = task_factory_type::create(policy,
                                                Chol<Uplo::Upper,Control<AlgoChol::ByBlocksVar1>::Chol>
                                                ::TaskFunctor<ParallelForType,value_type>(aa));
      
      // manage dependence
      task_factory_type::addDependence(policy, f, aa.Future());
      aa.setFuture(f);
      
      // spawn a task
      task_factory_type::spawn(policy, f);
      
      return 0;
    }
    
    KOKKOS_INLINE_FUNCTION
    static int genTrsmTasks(typename CrsTaskViewType::policy_type &policy,
                            CrsTaskViewType &A,
                            CrsTaskViewType &B) {
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
                   Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Control<AlgoChol::ByBlocksVar1>::Trsm>
                   ::TaskFunctor<ParallelForType,double,
                   value_type,value_type>(Diag::NonUnit, 1.0, aa, bb));
        
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
    
    KOKKOS_INLINE_FUNCTION
    static int genHerkTasks(typename CrsTaskViewType::policy_type &policy,
                            CrsTaskViewType &A,
                            CrsTaskViewType &C) {
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
        value_type &aa = a.Value(i);
        
        c.setView(C, row_at_i);
        
        ordinal_type idx = 0;
        for (ordinal_type j=i;j<nnz && (idx > -2);++j) {
          const ordinal_type col_at_j = a.Col(j);
          value_type &bb = a.Value(j);
          
          if (row_at_i == col_at_j) {
            idx = c.Index(row_at_i, idx);
            if (idx >= 0) {
              value_type &cc = c.Value(idx);
              future_type f = task_factory_type
                ::create(policy, 
                         Herk<Uplo::Upper,Trans::ConjTranspose,Control<AlgoChol::ByBlocksVar1>::Herk>
                         ::TaskFunctor<ParallelForType,double,
                         value_type,value_type>(-1.0, aa, 1.0, cc));
            
              // dependence
              task_factory_type::addDependence(policy, f, aa.Future());              
            
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
                         Gemm<Trans::ConjTranspose,Trans::NoTranspose,Control<AlgoChol::ByBlocksVar1>::Gemm>
                         ::TaskFunctor<ParallelForType,double,
                         value_type,value_type,value_type>(-1.0, aa, bb, 1.0, cc));
            
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
  
  };

  template<template<int> class ControlType,
           typename ParallelForType,
           typename CrsTaskViewType>
  class FunctionChol<Uplo::Upper,AlgoChol::ByBlocksVar1,
                     ControlType,
                     ParallelForType,
                     CrsTaskViewType> {
  private:
    int _r_val;

  public:
    KOKKOS_INLINE_FUNCTION
    int getReturnValue() const {
      return _r_val;
    }

    KOKKOS_INLINE_FUNCTION
    FunctionChol(typename CrsTaskViewType::policy_type &policy,
                 const typename CrsTaskViewType::policy_type::member_type &member,
                 CrsTaskViewType &A) {
      // this task generation should be done by a root
      // ---------------------------------------------
      if (member.team_rank() == 0) {
        CrsTaskViewType ATL, ATR,      A00, A01, A02,
          /**/          ABL, ABR,      A10, A11, A12,
          /**/                         A20, A21, A22;

        Part_2x2(A,  ATL, ATR,
                 /**/ABL, ABR,
                 0, 0, Partition::TopLeft);

        while (ATL.NumRows() < A.NumRows()) {
          Part_2x2_to_3x3(ATL, ATR, /**/  A00, A01, A02,
                          /*******/ /**/  A10, A11, A12,
                          ABL, ABR, /**/  A20, A21, A22,
                          1, 1, Partition::BottomRight);
          // -----------------------------------------------------

          // A11 = chol(A11)
          CholUpperByBlocksVar1<ControlType,ParallelForType,CrsTaskViewType>
            ::genScalarTask(policy, A11);

          // A12 = inv(triu(A11)') * A12
          CholUpperByBlocksVar1<ControlType,ParallelForType,CrsTaskViewType>
            ::genTrsmTasks(policy, A11, A12);

          // A22 = A22 - A12' * A12
          CholUpperByBlocksVar1<ControlType,ParallelForType,CrsTaskViewType>
            ::genHerkTasks(policy, A12, A22);

          // -----------------------------------------------------
          Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                           A10, A11, A12, /**/ /******/
                           A20, A21, A22, /**/ ABL, ABR,
                           Partition::TopLeft);
        }
      }

      _r_val = 0;
    }
  };
}

#endif
