#ifndef __HERK_UPPER_CONJTRANS_DENSE_BY_BLOCKS_HPP__
#define __HERK_UPPER_CONJTRANS_DENSE_BY_BLOCKS_HPP__

/// \file Tacho_Herk_Upper_ConjTrans_DenseByBlocks.hpp
/// \brief Hermitian rank-k By Blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  // Herk-By-Blocks
  // ==============
  template<int ArgVariant, template<int,int> class ControlType>
  class Herk<Uplo::Upper,Trans::ConjTranspose,
             AlgoHerk::DenseByBlocks,ArgVariant,ControlType> {
  public:
    template<typename PolicyType,
             typename MemberType,
             typename ScalarType,
             typename DenseTaskViewTypeA,
             typename DenseTaskViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const ScalarType alpha,
                      DenseTaskViewTypeA &A,
                      const ScalarType beta,
                      DenseTaskViewTypeC &C) {
#ifdef TACHO_EXECUTE_TASKS_SERIAL
#else      
      typedef typename DenseTaskViewTypeA::value_type::future_type future_type;
      TaskFactory factory;
#endif

      if (member.team_rank() == 0) {
        const auto pend = A.NumRows();
        for (auto p=0;p<pend;++p) {
          const auto beta_select = (p > 0 ? ScalarType(1.0) : beta);
          const auto k2end = C.NumCols();
          for (auto k2=0;k2<k2end;++k2) {
            for (auto k1=0;k1<(k2+1);++k1) {
              auto &aa = A.Value(p,  k1);
              auto &cc = C.Value(k1, k2);
              if (k1 == k2) {
#ifdef TACHO_EXECUTE_TASKS_SERIAL
                Herk<Uplo::Upper,Trans::ConjTranspose,
                  CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Herk)>
                  ::invoke(policy, member, alpha, aa, beta_select, cc);
#else
                future_type f = factory.create<future_type>
                  (policy,
                   Herk<Uplo::Upper,Trans::ConjTranspose,
                   CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Herk)>
                   ::createTaskFunctor(policy, alpha, aa, beta_select, cc), 2);
                
                // dependence
                factory.depend(policy, f, aa.Future());
                factory.depend(policy, f, cc.Future());

                // place task signature on y
                cc.setFuture(f);

                // spawn a task
                factory.spawn(policy, f);
#endif
              } else {
                auto &bb = A.Value(p, k2);
#ifdef TACHO_EXECUTE_TASKS_SERIAL
                Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                  CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Gemm)>
                  ::invoke(policy, member, alpha, aa, bb, beta_select, cc);
#else
                future_type f = factory.create<future_type>
                  (policy,
                   Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                   CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Gemm)>
                   ::createTaskFunctor(policy, alpha, aa, bb, beta_select, cc), 3);
                
                // dependence
                factory.depend(policy, f, aa.Future());
                factory.depend(policy, f, bb.Future());
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
      }
      
      return 0;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeC>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      ScalarType _alpha, _beta;
      ExecViewTypeA _A;
      ExecViewTypeC _C;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ScalarType alpha,
                  const ExecViewTypeA &A,
                  const ScalarType beta,
                  const ExecViewTypeC &C)
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _C(C),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Herk"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Herk::invoke(_policy, _policy.member_single(),
                             _alpha, _A, _beta, _C);
        _C.setFuture(typename ExecViewTypeC::future_type());
      }
      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        const int ierr = Herk::invoke(_policy, member,
                                      _alpha, _A, _beta, _C);

        // return for only team leader
        if (!member.team_rank()) {
          _C.setFuture(typename ExecViewTypeC::future_type());
          r_val = ierr;
        }
      }

    };

    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static
    TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeC>
    createTaskFunctor(const PolicyType &policy,
                      const ScalarType alpha,
                      const ExecViewTypeA &A,
                      const ScalarType beta,
                      const ExecViewTypeC &C) {
      return TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeC>
        (policy, alpha, A, beta, C);
    }


  };

}

#endif
