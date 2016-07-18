#ifndef __TRSM_LEFT_UPPER_CONJTRANS_DENSE_BY_BLOCKS_HPP__
#define __TRSM_LEFT_UPPER_CONJTRANS_DENSE_BY_BLOCKS_HPP__

/// \file Tacho_Trsm_Left_Upper_ConjTrans_DenseByBlocks.hpp
/// \brief TRSM-By-Blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<int ArgVariant, template<int,int> class ControlType>
  class Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
             AlgoTrsm::DenseByBlocks,ArgVariant,ControlType> {
  public:
    // data-parallel interface with nested task generation
    // ===================================================
    template<typename PolicyType,
             typename MemberType,
             typename ScalarType,
             typename DenseTaskViewTypeA,
             typename DenseTaskViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const int diagA,
                      const ScalarType alpha,
                      DenseTaskViewTypeA &A,
                      DenseTaskViewTypeB &B) {

#ifdef TACHO_EXECUTE_TASKS_SERIAL
#else
      typedef typename DenseTaskViewTypeA::value_type::future_type future_type;
      TaskFactory factory;
#endif
      
      if (member.team_rank() == 0) {
        DenseTaskViewTypeA ATL, ATR,      A00, A01, A02,
          /**/             ABL, ABR,      A10, A11, A12,
          /**/                            A20, A21, A22;
        
        DenseTaskViewTypeA BT,            B0,
          /**/             BB,            B1,
          /**/                            B2;

        Part_2x2(A,  ATL, ATR,
                 /**/ABL, ABR,
                 0, 0, Partition::TopLeft);

        Part_2x1(B,  BT,
                 /**/BB, 
                 0, Partition::Top);     

        while (ATL.NumRows() < A.NumRows()) {
          const auto alpha_select = (ATL.NumRows() == 0 ? ScalarType(1.0) : alpha);
          
          Part_2x2_to_3x3(ATL, ATR, /**/ A00, A01, A02,
                          /*******/ /**/ A10, A11, A12,
                          ABL, ABR, /**/ A20, A21, A22,
                          1, 1, Partition::BottomRight);
          
          Part_2x1_to_3x1(BT,  /**/ B0,
                          /**/ /**/ B1,
                          BB,  /**/ B2,
                          1, Partition::Bottom);  
          
          //------------------------------------------------------------
          auto &aa = A11.Value(0, 0);
          const auto ncols = B1.NumCols();
          for (auto j=0;j<ncols;++j) {
            auto &bb = B1.Value(0, j);

#ifdef TACHO_EXECUTE_TASKS_SERIAL
            Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
              CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Trsm)>
              ::invoke(policy, member, diagA, alpha_select, aa, bb);
#else
            future_type f = factory.create<future_type>
              (policy,
               Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
               CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Trsm)>
               ::createTaskFunctor(policy, diagA, alpha_select, aa, bb), 2);

            // trsm dependence
            factory.depend(policy, f, aa.Future());
            factory.depend(policy, f, bb.Future());

            // place task signature on b
            bb.setFuture(f);

            // spawn a task
            factory.spawn(policy, f);
#endif
          }

          Gemm<Trans::ConjTranspose,Trans::NoTranspose,
            AlgoGemm::DenseByBlocks,Variant::One>
            ::invoke(policy, member,
                     -1.0, A12, B1, alpha_select, B2);

          //------------------------------------------------------------
          Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                           A10, A11, A12, /**/ /******/
                           A20, A21, A22, /**/ ABL, ABR,
                           Partition::TopLeft);

          Merge_3x1_to_2x1(B0, /**/ BT,
                           B1, /**/ /**/
                           B2, /**/ BB,
                           Partition::Top);
        }
      }
      return 0;
    }



    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;
      
    private:
      int _diagA;
      ScalarType _alpha;

      ExecViewTypeA _A;
      ExecViewTypeB _B;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const int diagA,
                  const ScalarType alpha,
                  const ExecViewTypeA &A,
                  const ExecViewTypeB &B)
        : _diagA(diagA),
          _alpha(alpha),
          _A(A),
          _B(B),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Trsm"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Trsm::invoke(_policy, _policy.member_single(),
                             _diagA, _alpha, _A, _B);
        _B.setFuture(typename ExecViewTypeB::future_type());
      }
      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        const int ierr = Trsm::invoke(_policy, member,
                                      _diagA, _alpha, _A, _B);

        // return for only team leader
        if (!member.team_rank()) {
          _B.setFuture(typename ExecViewTypeB::future_type());
          r_val = ierr;
        }
      }

    };

    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static
    TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB>
    createTaskFunctor(const PolicyType &policy,
                      const int diagA,
                      const ScalarType alpha,
                      const ExecViewTypeA &A,
                      const ExecViewTypeB &B) {
      return TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB>
        (policy, diagA, alpha, A, B);
    }

  };

}

#endif
