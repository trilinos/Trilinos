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
      
      typedef typename DenseTaskViewTypeA::ordinal_type ordinal_type;
      typedef typename DenseTaskViewTypeA::value_type   value_type;
      typedef typename value_type::future_type          future_type;

      TaskFactory factory;

      if (member.team_rank() == 0) {
        for (ordinal_type p=0;p<A.NumRows();++p) {
          const ScalarType beta_select = (p > 0 ? ScalarType(1.0) : beta);
          for (ordinal_type k2=0;k2<C.NumCols();++k2) {
            for (ordinal_type k1=0;k1<(k2+1);++k1) {
              value_type &aa = A.Value(p,  k1);
              value_type &cc = C.Value(k1, k2);
              if (k1 == k2) {
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
              } else {
                value_type &bb = A.Value(p, k2);
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
