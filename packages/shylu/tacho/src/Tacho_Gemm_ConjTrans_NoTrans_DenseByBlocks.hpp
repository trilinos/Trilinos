#ifndef __TACHO_GEMM_CONJTRANS_NOTRANS_DENSE_BY_BLOCKS_HPP__
#define __TACHO_GEMM_CONJTRANS_NOTRANS_DENSE_BY_BLOCKS_HPP__

/// \file Tacho_Gemm_ConjTrans_NoTrans_DenseByBlocks.hpp
/// \brief Dense Gemm By Blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  // Gemm-By-Blocks
  // ==============
  template<int ArgVariant, template<int,int> class ControlType>
  class Gemm<Trans::ConjTranspose,Trans::NoTranspose,
             AlgoGemm::DenseByBlocks,ArgVariant,ControlType> {
  public:
    template<typename PolicyType,
             typename MemberType,
             typename ScalarType,
             typename DenseTaskViewTypeA,
             typename DenseTaskViewTypeB,
             typename DenseTaskViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const ScalarType alpha,
                      DenseTaskViewTypeA &A,
                      DenseTaskViewTypeB &B,
                      const ScalarType beta,
                      DenseTaskViewTypeC &C) {
      // static_assert( Kokkos::Impl::is_same<
      //                typename DenseMatrixTypeA::space_type,
      //                typename DenseMatrixTypeB::space_type
      //                >::value &&
      //                Kokkos::Impl::is_same<
      //                typename DenseMatrixTypeB::space_type,
      //                typename DenseMatrixTypeC::space_type
      //                >::value,
      //                "Space type of input matrices does not match" );
      
      // typedef ScalarType scalar_type;
      typedef typename DenseTaskViewTypeA::ordinal_type ordinal_type;
      typedef typename DenseTaskViewTypeA::value_type   value_type;
      typedef typename DenseTaskViewTypeA::future_type  future_type;

      TaskFactory factory;

      const int max_task_dependence = 3;
      if (member.team_rank() == 0) {
        for (ordinal_type p=0;p<A.NumRows();++p) {
          const ScalarType beta_select = (p > 0 ? ScalarType(1.0) : beta);
          for (ordinal_type k2=0;k2<C.NumCols();++k2) {
            value_type &bb = B.Value(p, k2);
            for (ordinal_type k1=0;k1<C.NumRows();++k1) {
              value_type &aa = A.Value(p,  k1);
              value_type &cc = C.Value(k1, k2);

              future_type f = factory.create<future_type>
                (policy,
                 Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                 CtrlDetail(ControlType,AlgoGemm::DenseByBlocks,ArgVariant,Gemm)>
                 ::createTaskFunctor(policy, alpha, aa, bb, beta_select, cc), 
                 max_task_dependence);

              // dependence
              factory.depend(policy, f, aa.Future());
              factory.depend(policy, f, bb.Future());

              // self
              factory.depend(policy, f, cc.Future());
              
              // place task signature on y
              cc.setFuture(f);
              
              // spawn a task
              factory.spawn(policy, f);
            }
          }
        }
      }
      
      return 0;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeC>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      ScalarType _alpha, _beta;
      ExecViewTypeA _A;
      ExecViewTypeB _B;
      ExecViewTypeC _C;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ScalarType alpha,
                  const ExecViewTypeA &A,
                  const ExecViewTypeB &B,
                  const ScalarType beta,
                  const ExecViewTypeC &C)
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _B(B),
          _C(C),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Dense::GemmByBlocks"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Gemm::invoke(_policy, _policy.member_single(),
                             _alpha, _A, _B, _beta, _C);
        _C.setFuture(typename ExecViewTypeC::future_type());
      }

      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        const int ierr = Gemm::invoke(_policy, member,
                                      _alpha, _A, _B, _beta, _C);

        // return for only team leader
        if (member.team_rank() == 0) { 
          _C.setFuture(typename ExecViewTypeC::future_type());
          r_val = ierr; 
        }
      }

    };

    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static
    TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB,ExecViewTypeC>
    createTaskFunctor(const PolicyType &policy,
                      const ScalarType alpha,
                      const ExecViewTypeA &A,
                      const ExecViewTypeB &B,
                      const ScalarType beta,
                      const ExecViewTypeC &C) {
      return TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB,ExecViewTypeC>
        (policy, alpha, A, B, beta, C);
    }

  };

}

#endif
