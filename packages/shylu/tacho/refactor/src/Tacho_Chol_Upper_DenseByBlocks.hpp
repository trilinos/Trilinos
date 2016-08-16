#ifndef __TACHO_CHOL_UPPER_DENSE_BY_BLOCKS_HPP__
#define __TACHO_CHOL_UPPER_DENSE_BY_BLOCKS_HPP__

/// \file Tacho_Chol_Upper_DenseByBlocks.hpp
/// \brief Dense Chol By Blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  // Chol-By-Blocks
  // ==============
  template<int ArgVariant, template<int,int> class ControlType>
  class Chol<Uplo::Upper,
             AlgoChol::DenseByBlocks,ArgVariant,ControlType> {
  public:
    template<typename PolicyType,
             typename MemberType,
             typename DenseTaskViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      MemberType &member,
                      DenseTaskViewTypeA &A) {
#ifdef TACHO_EXECUTE_TASKS_SERIAL
#else
      typedef typename DenseTaskViewTypeA::value_type::future_type future_type;
#endif

      if (member.team_rank() == 0) {
        DenseTaskViewTypeA ATL, ATR,      A00, A01, A02,
          /**/             ABL, ABR,      A10, A11, A12,
          /**/                            A20, A21, A22;

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
            auto &aa = A11.Value(0, 0);

#ifdef TACHO_EXECUTE_TASKS_SERIAL
            Chol<Uplo::Upper,
              CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Chol)>
              ::invoke(policy, member, aa);
#else
            const auto task_type     = Kokkos::TaskTeam;
            const auto task_priority = Kokkos::TaskHighPriority;

            const future_type f =
              policy.task_spawn(Chol<Uplo::Upper,
                                CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Chol)>
                                ::createTaskFunctor(policy, aa), 
                                aa.Future(),
                                task_type, task_priority);
            TACHO_TEST_FOR_ABORT(f.is_null(),
                                 ">> Tacho::DenseCholByBlocks(Upper) returns a null future (out of memory)");
            aa.setFuture(f);
#endif
          }

          // A12 = inv(triu(A11)') * A12
          {
            auto &aa = A11.Value(0, 0);
            
            const auto ncols = A12.NumCols();
            for (auto j=0;j<ncols;++j) {
              auto &bb = A12.Value(0, j);

#ifdef TACHO_EXECUTE_TASKS_SERIAL
              Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
                CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Trsm)>
                ::invoke(policy, member, Diag::NonUnit, 1.0, aa, bb);
#else
              const auto task_type     = Kokkos::TaskTeam;
              const auto task_priority = Kokkos::TaskRegularPriority;

              const future_type dep[] = { aa.Future(), bb.Future() };
              const future_type f =
                policy.task_spawn(Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
                                  CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Trsm)>
                                  ::createTaskFunctor(policy, Diag::NonUnit, 1.0, aa, bb), 
                                  policy.when_all(2,dep),
                                  task_type, task_priority);
              TACHO_TEST_FOR_ABORT(f.is_null(),
                                   ">> Tacho::DenseCholByBlocks(Upper) returns a null future (out of memory)");
              bb.setFuture(f);
#endif
            }
          }

          // A22 = A22 - A12' * A12
          {
            const auto ncols = A22.NumCols();
            for (auto j=0;j<ncols;++j) {
              {
                auto &aa = A12.Value(0, j);
                auto &cc = A22.Value(j, j);
#ifdef TACHO_EXECUTE_TASKS_SERIAL
                Herk<Uplo::Upper,Trans::ConjTranspose,
                  CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Herk)>
                  ::invoke(policy, member, -1.0, aa, 1.0, cc);
#else
                const auto task_type     = Kokkos::TaskTeam;
                const auto task_priority = Kokkos::TaskHighPriority;
                
                const future_type dep[] = { aa.Future(), cc.Future() };
                const future_type f =
                  policy.task_spawn(Herk<Uplo::Upper,Trans::ConjTranspose,
                                    CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Herk)>
                                    ::createTaskFunctor(policy, -1.0, aa, 1.0, cc), 
                                    policy.when_all(2,dep),
                                    task_type, task_priority);
                TACHO_TEST_FOR_ABORT(f.is_null(),
                                     ">> Tacho::DenseCholByBlocks(Upper) returns a null future (out of memory)");
                cc.setFuture(f);
#endif
              }
              for (auto i=0;i<j;++i) {
                auto &aa = A12.Value(0, i);
                auto &bb = A12.Value(0, j);
                auto &cc = A22.Value(i, j);

#ifdef TACHO_EXECUTE_TASKS_SERIAL
                Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                  CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Gemm)>
                  ::invoke(policy, member, -1.0, aa, bb, 1.0, cc);
#else
                const auto task_type     = Kokkos::TaskTeam;
                const auto task_priority = Kokkos::TaskHighPriority;
                
                const future_type dep[] = { aa.Future(), bb.Future(), cc.Future() };
                const future_type f =
                  policy.task_spawn(Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                                    CtrlDetail(ControlType,AlgoChol::DenseByBlocks,ArgVariant,Gemm)>
                                    ::createTaskFunctor(policy, -1.0, aa, bb, 1.0, cc), 
                                    policy.when_all(3,dep),
                                    task_type, task_priority);
                TACHO_TEST_FOR_ABORT(f.is_null(),
                                     ">> Tacho::DenseCholByBlocks(Upper) returns a null future (out of memory)");
                cc.setFuture(f);
#endif
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
    // ===================\=========
    template<typename PolicyType,
             typename ExecViewTypeA>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      ExecViewTypeA _A;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor() {}

      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ExecViewTypeA &A)
        : _A(A),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Dense::CholByBlocks"; }

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
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
