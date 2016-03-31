#ifndef __TACHO_CHOL_UPPER_BY_BLOCKS_HPP__
#define __TACHO_CHOL_UPPER_BY_BLOCKS_HPP__

/// \file Tacho_Chol_Upper_ByBlocks.hpp
/// \brief Cholesky factorization by-blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 
  
  // Chol-By-Blocks
  // =================================================================================
  template<int ArgVariant, template<int,int> class ControlType>
  class Chol<Uplo::Upper,
             AlgoChol::ByBlocks,ArgVariant,ControlType> {
  public:

    // function interface
    // ==================
    template<typename PolicyType,
             typename MemberType,
             typename CrsTaskViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      CrsTaskViewTypeA &A) {

      typedef typename CrsTaskViewTypeA::ordinal_type      ordinal_type;
      typedef typename CrsTaskViewTypeA::value_type        value_type;
      typedef typename CrsTaskViewTypeA::row_view_type     row_view_type;

      typedef typename CrsTaskViewTypeA::future_type       future_type;

      TaskFactory factory;

      if (member.team_rank() == 0) {
        CrsTaskViewTypeA ATL, ATR,      A00, A01, A02,
          /**/           ABL, ABR,      A10, A11, A12,
          /**/                          A20, A21, A22;

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
          {
            row_view_type a(A11, 0); 
            value_type &aa = a.Value(0);

            if (aa.isEmpty()) {
              continue;
            } else {
              // construct a task
              future_type f = factory.create<future_type>
                (policy,
                 Chol<Uplo::Upper,
                 CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Chol)>
                 ::createTaskFunctor(policy, aa), 1);
              
              // manage dependence
              factory.depend(policy, f, aa.Future());
              aa.setFuture(f);
              
              // spawn a task
              factory.spawn(policy, f);
            } 
          }

          // A12 = inv(triu(A11)') * A12
          {
            row_view_type a(A11, 0), b(A12, 0); 
            value_type &aa = a.Value(0);
            
            const ordinal_type nnz = b.NumNonZeros();
            for (ordinal_type j=0;j<nnz;++j) {
              value_type &bb = b.Value(j);

              if (bb.isEmpty()) {
                continue;
              } else {
                future_type f = factory.create<future_type>
                  (policy, 
                   Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
                   CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Trsm)>
                   ::createTaskFunctor(policy, Diag::NonUnit, 1.0, aa, bb), 2);
                
                // trsm dependence
                factory.depend(policy, f, aa.Future());
                factory.depend(policy, f, bb.Future());
                
                // place task signature on b
                bb.setFuture(f);
                
                // spawn a task
                factory.spawn(policy, f);              
              } 
            }
          }

          // A22 = A22 - A12' * A12
          {
            // case that X.transpose, A.no_transpose, Y.no_transpose
            
            row_view_type a(A12,0);
            
            const ordinal_type nnz = a.NumNonZeros();
            
            // update herk
            for (ordinal_type i=0;i<nnz;++i) {
              const ordinal_type row_at_i = a.Col(i);
              value_type &aa = a.Value(i);
              
              if (aa.isEmpty()) {
                continue;
              } else {
                row_view_type c(A22, row_at_i);
                
                ordinal_type idx = 0;
                for (ordinal_type j=i;j<nnz && (idx > -2);++j) {
                  const ordinal_type col_at_j = a.Col(j);
                  value_type &bb = a.Value(j);
                  
                  if (bb.isEmpty()) {
                    continue;
                  } else {
                    if (row_at_i == col_at_j) {
                      idx = c.Index(row_at_i, idx);
                      if (idx >= 0) {
                        value_type &cc = c.Value(idx);

                        if (cc.isEmpty()) {
                          continue;
                        } else {
                          future_type f = factory.create<future_type>
                            (policy, 
                             Herk<Uplo::Upper,Trans::ConjTranspose,
                             CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Herk)>
                             ::createTaskFunctor(-1.0, aa, 1.0, cc), 2);
                          
                          // dependence
                          factory.depend(policy, f, aa.Future());              
                          factory.depend(policy, f, cc.Future());
                          
                          // place task signature on y
                          cc.setFuture(f);
                          
                          // spawn a task
                          factory.spawn(policy, f);
                        } 
                      }
                    }
                  } 
                } else {
                  idx = c.Index(col_at_j, idx);
                  if (idx >= 0) {
                    value_type &cc = c.Value(idx);
                    
                    if (cc.isEmpty()) {
                      continue;
                    } else {
                      future_type f = factory.create<future_type>
                        (policy, 
                         Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                         CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Gemm)>
                         ::createTaskFunctor(-1.0, aa, bb, 1.0, cc), 3);
                      
                      // dependence
                      factory.depend(policy, f, aa.Future());
                      factory.depend(policy, f, bb.Future());
                      factory.depend(policy, f, cc.Future());
                      
                      // place task signature on y
                      cc.setFuture(f);
                      
                      // spawn a task
                      factory.spawn(policy, f);
                    } 
                  }
                }
              } 
            }
          }

          // -----------------------------------------------------
          Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                           A10, A11, A12, /**/ /******/
                           A20, A21, A22, /**/ ABL, ABR,
                           Partition::TopLeft);
        }
      }
      
      return 0;
    }
    
    // task-data parallel interface
    // ============================
    template<typename PolicyType,
             typename ExecViewTypeA>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;
      
    private:
      ExecViewTypeA _A;
      
      policy_type _policy;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ExecViewTypeA A)
        : _A(A),
          _policy(policy)
      { } 
      
      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "CholByBlocks"; }
      
      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Chol::invoke(_policy, _policy.member_single(), 
                             _A);
        _A.setFuture(typename ExecViewTypeA::future_type());
      }
      
      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        const int ierr = Chol::invoke(_policy, member,
                                      _A);

        // return for only team leader
        if (member.team_rank() == 0) {
          _A.setFuture(typename ExecViewTypeA::future_type());
          r_val = ierr;
        }
      }

    };

    template<typename PolicyType,
             typename ExecViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static
    TaskFunctor<PolicyType,ExecViewTypeA>
    createTaskFunctor(const PolicyType &policy,
                      const ExecViewTypeA &A) {
      return TaskFunctor<PolicyType,ExecViewTypeA>
        (policy, A);
    }

  };
}

#endif
