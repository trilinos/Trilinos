#ifndef __TACHO_GEMM_NOTRANS_NOTRANS_DENSE_BY_BLOCKS_HPP__
#define __TACHO_GEMM_NOTRANS_NOTRANS_DENSE_BY_BLOCKS_HPP__

/// \file Tacho_Gemm_NoTrans_NoTrans_DenseByBlocks.hpp
/// \brief Dense Gemm By Blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  // Gemm-By-Blocks
  // ==============
  template<int ArgVariant, template<int,int> class ControlType>
  class Gemm<Trans::NoTranspose,Trans::NoTranspose,
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
      
#ifdef TACHO_EXECUTE_TASKS_SERIAL
#else
      typedef typename DenseTaskViewTypeA::value_type::future_type future_type;
      TaskFactory factory;
#endif

      if (member.team_rank() == 0) {
        const auto pend = A.NumCols();
        for (auto p=0;p<pend;++p) {
          const auto beta_select = (p > 0 ? ScalarType(1.0) : beta);
          const auto k2end = C.NumCols();
          for (auto k2=0;k2<k2end;++k2) {
            auto &bb = B.Value(p, k2);
            const auto k1end = C.NumRows();
            for (auto k1=0;k1<k1end;++k1) {
              auto &aa = A.Value(k1, p );
              auto &cc = C.Value(k1, k2);

#ifdef TACHO_EXECUTE_TASKS_SERIAL
              Gemm<Trans::NoTranspose,Trans::NoTranspose,
                CtrlDetail(ControlType,AlgoGemm::DenseByBlocks,ArgVariant,Gemm)>
                ::invoke(policy, member, alpha, aa, bb, beta_select, cc);
#else
              future_type f = factory.create<future_type>
                (policy,
                 Gemm<Trans::NoTranspose,Trans::NoTranspose,
                 CtrlDetail(ControlType,AlgoGemm::DenseByBlocks,ArgVariant,Gemm)>
                 ::createTaskFunctor(policy, alpha, aa, bb, beta_select, cc), 3);

              // dependence
              factory.depend(policy, f, aa.Future());
              factory.depend(policy, f, bb.Future());

              // self
              factory.depend(policy, f, cc.Future());
              
              // place task signature on y
              cc.setFuture(f);
              
              // spawn a task
              factory.spawn(policy, f);
#endif

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
