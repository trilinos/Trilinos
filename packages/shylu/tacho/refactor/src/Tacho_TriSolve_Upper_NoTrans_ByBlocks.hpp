#ifndef __TACHO_TRI_SOLVE_UPPER_NOTRANS_BY_BLOCKS_HPP__
#define __TACHO_TRI_SOLVE_UPPER_NOTRANS_BY_BLOCKS_HPP__

/// \file Tacho_TriSolve_Upper_NoTrans_ByBlocks.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Tacho {

  template<int ArgVariant,
           template<int,int> class ControlType>
  class TriSolve<Uplo::Upper,Trans::NoTranspose,
                 AlgoTriSolve::ByBlocks,ArgVariant,ControlType>{
  public:

    // data-parallel interface
    // =======================
    template<typename PolicyType,
             typename MemberType,
             typename CrsTaskViewTypeA,
             typename DenseTaskViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      MemberType &member,
                      const int diagA,
                      CrsTaskViewTypeA &A,
                      DenseTaskViewTypeB &B,
                      unsigned int &part) {
#ifdef TACHO_EXECUTE_TASKS_SERIAL
#else
      typedef typename CrsTaskViewTypeA::future_type future_type;
#endif
      typedef typename CrsTaskViewTypeA::row_view_type row_view_type;

      // ---------------------------------------------
      if (member.team_rank() == 0) {
        const unsigned int ntasks_window  = TaskWindow::TriSolveByBlocks;
        /**/  unsigned int ntasks_spawned = 0;

        CrsTaskViewTypeA ATL, ATR,      A00, A01, A02,
          /**/           ABL, ABR,      A10, A11, A12,
          /**/                          A20, A21, A22;

        DenseTaskViewTypeB BT,      B0,
          /**/             BB,      B1,
          /**/                      B2;

        Part_2x2(A,  ATL, ATR,
                 /**/ABL, ABR,
                 part, part, Partition::BottomRight);

        Part_2x1(B,  BT,
                 /**/BB,
                 part, Partition::Bottom);

        while (ABR.NumRows() < A.NumRows()) {
          Part_2x2_to_3x3(ATL, ATR, /**/  A00, A01, A02,
                          /*******/ /**/  A10, A11, A12,
                          ABL, ABR, /**/  A20, A21, A22,
                          1, 1, Partition::TopLeft);

          Part_2x1_to_3x1(BT,  /**/  B0,
                          /**/ /**/  B1,
                          BB,  /**/  B2,
                          1, Partition::Top);

          // -----------------------------------------------------

          // B1 = B1 - A12*B2;
          {
            row_view_type a(A12,0);
            const auto nnz = a.NumNonZeros();

            for (auto i=0;i<nnz;++i) {
              const auto row_at_i = a.Col(i);
              auto &aa = a.Value(i);
              
              if (!aa.isNull()) {
                for (auto j=0;j<B1.NumCols();++j) {
                  auto &bb = B2.Value(row_at_i, j);
                  auto &cc = B1.Value(0, j);
                  
#ifdef TACHO_EXECUTE_TASKS_SERIAL
                  Gemm<Trans::NoTranspose,Trans::NoTranspose,
                    CtrlDetail(ControlType,AlgoTriSolve::ByBlocks,ArgVariant,Gemm)>
                    ::invoke(policy, member, -1.0, aa, bb, 1.0, cc);
#else
                  switch (ArgVariant) {
                  case Variant::Three: 
                    {
                      Gemm<Trans::NoTranspose,Trans::NoTranspose,
                        CtrlDetail(ControlType,AlgoTriSolve::ByBlocks,ArgVariant,Gemm)>
                        ::invoke(policy, member, -1.0, aa, bb, 1.0, cc);
                      break;
                    }
                  case Variant::One:
                  case Variant::Two: 
                    {
                      const auto task_type     = Kokkos::TaskTeam;
                      const auto task_priority = Kokkos::TaskRegularPriority;

                      const future_type dep[] = { aa.Future(), bb.Future(), cc.Future() };
                      const auto f = 
                        policy.task_spawn(Gemm<Trans::NoTranspose,Trans::NoTranspose,
                                          CtrlDetail(ControlType,AlgoTriSolve::ByBlocks,ArgVariant,Gemm)>
                                          ::createTaskFunctor(policy, -1.0, aa, bb, 1.0, cc),
                                          policy.when_all(3,dep),
                                          task_type, task_priority);
                      TACHO_TEST_FOR_ABORT(f.is_null(),
                                           ">> Tacho::DenseGemmByBlocks(NoTrans,NoTrans) returns a null future (out of memory)");
                      cc.setFuture(f);
                      break;
                    }
                  }
                  ++ntasks_spawned;
#endif
                }
              }
            }
          }

          // B1 = inv(triu(A11))*B1
          {
            row_view_type a(A11,0);
            auto &aa = a.Value(0);

            if (!aa.isNull()) {
              for (auto j=0;j<B1.NumCols();++j) {
                auto &bb = B1.Value(0, j);
                
#ifdef TACHO_EXECUTE_TASKS_SERIAL
                Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
                  CtrlDetail(ControlType,AlgoTriSolve::ByBlocks,ArgVariant,Trsm)>
                  ::invoke(policy, member, diagA, 1.0, aa, bb);
#else
                switch (ArgVariant) {
                case Variant::Three:
                  {
                    Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
                      CtrlDetail(ControlType,AlgoTriSolve::ByBlocks,ArgVariant,Trsm)>
                      ::invoke(policy, member, diagA, 1.0, aa, bb);
                    break;
                  }
                case Variant::One:
                case Variant::Two:
                  {
                    const auto task_type     = Kokkos::TaskTeam;
                    const auto task_priority = Kokkos::TaskHighPriority;

                    const future_type dep[] = { aa.Future(), bb.Future() };
                    const auto f = 
                      policy.task_spawn(Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
                                        CtrlDetail(ControlType,AlgoTriSolve::ByBlocks,ArgVariant,Trsm)>
                                        ::createTaskFunctor(policy, diagA, 1.0, aa, bb), 
                                        policy.when_all(2,dep),
                                        task_type, task_priority);
                    TACHO_TEST_FOR_ABORT(f.is_null(),
                                         ">> Tacho::TriSolveByBlocks(Upper,NoTrans) returns a null future (out of memory)");
                    bb.setFuture(f);
                    break;
                  }
                }
                ++ntasks_spawned;
#endif
              }
            }
          }

          // -----------------------------------------------------
          Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                           A10, A11, A12, /**/ /******/
                           A20, A21, A22, /**/ ABL, ABR,
                           Partition::BottomRight);

          Merge_3x1_to_2x1(B0, /**/   BT,
                           B1, /**/  /**/
                           B2, /**/   BB,
                           Partition::Bottom);

          if (ntasks_spawned > ntasks_window)
            break;
        }

        part = BB.NumRows();
      }
      return 0;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      int _diagA;

      ExecViewTypeA _A;
      ExecViewTypeB _B;

      PolicyType _policy;

      unsigned int _part;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor() {}
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const int diagA,
                  const ExecViewTypeA &A,
                  const ExecViewTypeB &B)
        : _diagA(diagA),
          _A(A),
          _B(B),
          _policy(policy),
          _part(0)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "TriSolve::Upper::NoTrans::ByBlocks"; }

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {

        if (member.team_rank() == 0) {
          //_policy.clear_dependence(this);

          const int ierr = TriSolve::invoke(_policy, member,
                                            _diagA, _A, _B, _part);

          if (_part < _B.NumRows()) {
            _policy.respawn(this, Kokkos::TaskLowPriority);
          } else {
            _part = 0;
            _B.setFuture(typename ExecViewTypeB::future_type());
          }
          r_val = ierr;
        }
      }

    };

    template<typename PolicyType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static
    TaskFunctor<PolicyType,ExecViewTypeA,ExecViewTypeB>
    createTaskFunctor(const PolicyType &policy,
                      const int diagA,
                      const ExecViewTypeA &A,
                      const ExecViewTypeB &B) {
      return TaskFunctor<PolicyType,ExecViewTypeA,ExecViewTypeB>
        (policy, diagA, A, B);
    }

  };
}
#endif
